#ifndef RS_HPP
#define RS_HPP

#include <string.h>
#include <stdint.h>
#include <vector>
#include <algorithm>

#include "poly.hpp"
#include "gf.hpp"

#if !defined RS_DEBUG && !defined __CC_ARM && !defined RS_NO_ASSERT
#include <assert.h>
#else
#define assert(dummy)
#endif

#define RS_DEBUG false

namespace RS
{

#define MSG_CNT 3   ///< number of message-length polynomials
#define POLY_CNT 14 ///< number of auxiliary polynomials

    template <const uint8_t msg_length, ///< message length excluding ECC
            const uint8_t ecc_length  ///< ECC length
    >
    class ReedSolomon
    {
    public:
        ReedSolomon()
        {
            /* Use 16-bit sizes to avoid overflow when ecc_length is large */
            const uint16_t enc_len  = (uint16_t)msg_length + (uint16_t)ecc_length;
            const uint16_t poly_len = (uint16_t)ecc_length * 2u;
            uint8_t** memptr = &memory;
            uint16_t offset = 0;

            /* Initialize polynomial slots (layout depends on template params) */
            polynoms[0].Init(ID_MSG_IN, offset, enc_len, memptr);
            offset += enc_len;

            polynoms[1].Init(ID_MSG_OUT, offset, enc_len, memptr);
            offset += enc_len;

            for(uint8_t i = ID_GENERATOR; i < ID_MSG_E; i++)
            {
                polynoms[i].Init(i, offset, poly_len, memptr);
                offset += poly_len;
            }

            polynoms[5].Init(ID_MSG_E, offset, enc_len, memptr);
            offset += enc_len;

            for(uint8_t i = ID_TPOLY3; i < ID_ERR_EVAL+2; i++)
            {
                polynoms[i].Init(i, offset, poly_len, memptr);
                offset += poly_len;
            }
        }

        ~ReedSolomon()
        {
            /* Dummy destructor — clear memory pointer (stack memory will be reclaimed automatically) */
            memory = NULL;
        }

        /// Encode a message block (systematic) - produce ECC only into dst
        void EncodeBlock(const void* src, void* dst)
        {
            assert(msg_length + ecc_length < 256);

            static uint8_t generator_cache[ecc_length+1] = {0};
            static bool    generator_cached = false;

            /* allocate polynomial workspace on stack */
            uint8_t stack_memory[MSG_CNT * msg_length + POLY_CNT * ecc_length * 2];
            this->memory = stack_memory;

            const uint8_t* src_ptr = (const uint8_t*) src;
            uint8_t* dst_ptr = (uint8_t*) dst;

            Poly *msg_in  = &polynoms[ID_MSG_IN];
            Poly *msg_out = &polynoms[ID_MSG_OUT];
            Poly *gen     = &polynoms[ID_GENERATOR];

            msg_in->Reset();
            msg_out->Reset();

            if(generator_cached)
            {
                gen->Set(generator_cache, gen->length);
            }
            else
            {
                GeneratorPoly();
                memcpy(generator_cache, gen->ptr(), gen->length);
                generator_cached = true;
            }

            msg_in->Set(src_ptr, msg_length);
            msg_out->Set(src_ptr, msg_length);
            msg_out->length = (uint16_t)(msg_in->length + ecc_length);

            uint8_t coef = 0;
            for(uint16_t i = 0; i < msg_length; i++)
            {
                coef = msg_out->at((uint16_t)i);
                if(coef != 0)
                {
                    for(uint32_t j = 1; j < gen->length; j++)
                    {
                        msg_out->at((uint16_t)(i+j)) ^= gf::mul(gen->at((uint16_t)j), coef);
                    }
                }
            }

            memcpy(dst_ptr, msg_out->ptr()+msg_length, ecc_length * sizeof(uint8_t));
        }

        /// Encode message full [msg|ecc]
        void Encode(const void* src, void* dst)
        {
            uint8_t* dst_ptr = (uint8_t*) dst;
            memcpy(dst_ptr, src, msg_length * sizeof(uint8_t));
            EncodeBlock(src, dst_ptr + msg_length);
        }

        /// Decode systematic block (existing routine)
        int DecodeBlock(const void* src, const void* ecc, void* dst, uint8_t* erase_pos = NULL, size_t erase_count = 0)
        {
            assert(msg_length + ecc_length < 256);

            const uint8_t *src_ptr = (const uint8_t*) src;
            const uint8_t *ecc_ptr = (const uint8_t*) ecc;
            uint8_t *dst_ptr = (uint8_t*) dst;

            const uint16_t src_len = (uint16_t)msg_length + (uint16_t)ecc_length;
            const uint16_t dst_len = (uint16_t)msg_length;

            bool ok;

            uint8_t stack_memory[MSG_CNT * msg_length + POLY_CNT * ecc_length * 2];
            this->memory = stack_memory;

            Poly *msg_in  = &polynoms[ID_MSG_IN];
            Poly *msg_out = &polynoms[ID_MSG_OUT];
            Poly *epos    = &polynoms[ID_ERASURES];

            msg_in->Set(src_ptr, msg_length);
            msg_in->Set(ecc_ptr, ecc_length, msg_length);
            msg_out->Copy(msg_in);

            if(erase_pos == NULL)
            {
                epos->length = 0;
            }
            else
            {
                epos->Set(erase_pos, (uint16_t)erase_count);
                for(uint16_t i = 0; i < epos->length; i++)
                {
                    msg_in->at(epos->at(i)) = 0;
                }
            }

            if(epos->length > ecc_length) return 1;

            Poly *synd   = &polynoms[ID_SYNDROMES];
            Poly *eloc   = &polynoms[ID_ERRORS_LOC];
            Poly *reloc  = &polynoms[ID_TPOLY1];
            Poly *err    = &polynoms[ID_ERRORS];
            Poly *forney = &polynoms[ID_FORNEY];

            CalcSyndromes(msg_in);

            bool has_errors = false;
            for(uint16_t i = 0; i < synd->length; i++)
            {
                if(synd->at((uint16_t)i) != 0)
                {
                    has_errors = true;
                    break;
                }
            }

            if(!has_errors) goto return_corrected_msg;

            CalcForneySyndromes(synd, epos, src_len);
            FindErrorLocator(forney, NULL, epos->length);

            reloc->length = eloc->length;
            for(int16_t i = (int16_t)eloc->length - 1, j = 0; i >= 0; i--, j++)
            {
                reloc->at((uint16_t)j) = eloc->at((uint16_t)i);
            }

            ok = FindErrors(reloc, src_len);
            if(!ok) return 1;
            if(err->length == 0) return 1;

            for(uint16_t i = 0; i < err->length; i++)
            {
                epos->Append(err->at((uint16_t)i));
            }

            CorrectErrata(synd, epos, msg_in);

            return_corrected_msg:
            msg_out->length = dst_len;
            memcpy(dst_ptr, msg_out->ptr(), msg_out->length * sizeof(uint8_t));
            return 0;
        }

        /// Decode wrapper for full codeword
        int Decode(const void* src, void* dst, uint8_t* erase_pos = NULL, size_t erase_count = 0)
        {
            const uint8_t *src_ptr = (const uint8_t*) src;
            const uint8_t *ecc_ptr = src_ptr + msg_length;
            return DecodeBlock(src, ecc_ptr, dst, erase_pos, erase_count);
        }

        /**
         * Non-systematic decoding (floating window / cyclic shifts) with decimation attempts.
         * First, try to find a quotient q(x) such that q*g matches received with Hamming distance ≤ t.
         * If that fails, try decimations/strict algebraic reconstruction and, importantly,
         * attempt algebraic correction treating rotated_code as a full codeword polynomial.
         */
        int DecodeNonSystematic(const void* src, void* dst,
                                uint8_t* erase_pos = NULL, size_t erase_count = 0,
                                uint8_t shift_step = 1, uint16_t max_shifts = 0,
                                const uint8_t* decimations = NULL, size_t decimations_count = 0)
        {
            assert(msg_length + ecc_length < 256);
            const uint16_t n = (uint16_t)msg_length + (uint16_t)ecc_length;
            if(max_shifts == 0) max_shifts = n;
            if(shift_step == 0) shift_step = 1;

            const uint8_t *src_ptr = (const uint8_t*) src;
            uint8_t *dst_ptr = (uint8_t*) dst;

            /* allocate memory on stack for polynomials */
            uint8_t stack_memory[MSG_CNT * msg_length + POLY_CNT * ecc_length * 2];
            this->memory = stack_memory;

            /* Pointers to polynomials we will use */
            Poly *gen = &polynoms[ID_GENERATOR];
            Poly *pPoly = &polynoms[ID_MSG_IN];   // general purpose
            Poly *quot = &polynoms[ID_MSG_E];     // store quotient (candidate message)
            Poly *mulres = &polynoms[ID_TPOLY2];  // multiplication result buffer

            /* Build generator polynomial (placed into polynoms[ID_GENERATOR]) */
            gen->Reset();
            GeneratorPoly();

            /* local temporary arrays */
            uint8_t rotated_code[512];
            uint8_t reconstructed_code[512];

            uint8_t t_correct = ecc_length / 2; // number of symbol errors allowed

            // 1) For each cyclic shift compute quotient q and c' = q * g; if Hamming(c', r) <= t -> accept q
            for(uint16_t shift = 0; shift < n; ++shift)
            {
                // rotated_code[i] = src[(i + shift) % n]
                for(uint8_t i = 0; i < n; ++i)
                    rotated_code[i] = src_ptr[(i + shift) % n];

                // set pPoly to rotated_code (polynomial p)
                pPoly->Set(rotated_code, n, 0);

                // compute quotient (floor division) -> quot
                PolyDivGetQuotient(pPoly, gen, quot);

                // multiply quotient by generator to get candidate codeword
                mulres->Reset();
                gf::poly_mul(quot, gen, mulres); // mulres now has length = quot.len + gen.len - 1

                // pad/truncate mulres to length n
                memset(reconstructed_code, 0, n);
                for(uint8_t i = 0; i < (mulres->length && mulres->length <= n ? mulres->length : (mulres->length)); ++i)
                {
                    if(i < n) reconstructed_code[i] = mulres->at(i);
                }
                // If mulres->length > n, we still only compare first n coefficients (practically shouldn't happen)

                // compute Hamming distance between reconstructed_code and rotated_code
                uint8_t diff = 0;
                for(uint8_t i = 0; i < n; ++i)
                {
                    if(reconstructed_code[i] != rotated_code[i])
                    {
                        ++diff;
                        if(diff > t_correct) break;
                    }
                }

                if(diff <= t_correct)
                {
                    // accept quotient as recovered message (quot->ptr())
                    uint8_t copy_len = (quot->length < msg_length) ? quot->length : msg_length;
                    memcpy(dst_ptr, quot->ptr(), copy_len);
                    if(copy_len < msg_length) memset(dst_ptr + copy_len, 0, msg_length - copy_len);
                    return 0;
                }

                // --- NEW: If simple Hamming check fails, try to algebraically correct rotated_code as full codeword ---
                // This attempts full syndrome-based correction over the whole polynomial of length n.
                {
                    Poly *msg_in  = &polynoms[ID_MSG_IN];
                    Poly *msg_out = &polynoms[ID_MSG_OUT];
                    Poly *epos    = &polynoms[ID_ERASURES];
                    Poly *synd    = &polynoms[ID_SYNDROMES];
                    Poly *forney  = &polynoms[ID_FORNEY];
                    Poly *errloc  = &polynoms[ID_ERRORS_LOC];
                    Poly *reloc   = &polynoms[ID_TPOLY3];
                    Poly *err     = &polynoms[ID_ERRORS];

                    // place rotated_code as a full polynomial of length n (no split msg/ecc)
                    msg_in->Set(rotated_code, (uint16_t)n, 0);
                    msg_out->Copy(msg_in);

                    if(erase_pos == NULL) epos->length = 0;
                    else {
                        epos->Set(erase_pos, (uint16_t)erase_count);
                        for(uint16_t i = 0; i < epos->length; i++) msg_in->at(epos->at((uint16_t)i)) = 0;
                    }
                    if(epos->length > ecc_length) {
                        // too many erasures, skip algebraic attempt for this shift
                    } else {
                        // compute syndromes treating msg_in as full length-n codeword
                        CalcSyndromes(msg_in);

                        bool has_errors = false;
                        for(uint8_t i = 0; i < synd->length; i++) if(synd->at(i) != 0) { has_errors = true; break; }

                        if(!has_errors)
                        {
                            // already clean (no errors) -> get quotient directly
                            msg_out->length = (uint16_t)n;
                        }
                        else
                        {
                            // attempt to correct
                            CalcForneySyndromes(synd, epos, n);
                            FindErrorLocator(forney, NULL, epos->length);

                            reloc->length = errloc->length;
                            for(int8_t i = errloc->length-1, j = 0; i >= 0; i--, j++)
                            {
                                reloc->at((uint16_t)j) = errloc->at((uint16_t)i);
                            }

                            if(FindErrors(reloc, n))
                            {
                                if(err->length != 0)
                                {
                                    for(uint8_t i = 0; i < err->length; i++) epos->Append(err->at(i));
                                    CorrectErrata(synd, epos, msg_in);
                                }
                                else
                                {
                                    // no error positions found
                                }
                            }
                            else
                            {
                                // unable to find error positions for this rotated_code
                                // continue to next shift
                                goto skip_algebraic_division;
                            }
                        }

                        // now msg_in/msg_out contain corrected full codeword polynomial,
                        // compute quotient = corrected_code / g(x) to extract original message polynomial
                        msg_out->length = (uint16_t)n;
                        // set pPoly to corrected full codeword
                        pPoly->Copy(msg_out);
                        PolyDivGetQuotient(pPoly, gen, quot);
                        // quot now contains candidate message coefficients
                        uint16_t copy_len = (quot->length < msg_length) ? quot->length : msg_length;
                        memcpy(dst_ptr, quot->ptr(), copy_len);
                        if(copy_len < msg_length) memset(dst_ptr + copy_len, 0, msg_length - copy_len);
                        return 0;
                    }
                }
                skip_algebraic_division: ; // label target for continue
            }

            /* 2) If strict method did not give result, try decimations/heuristics/fallback (unchanged) */
            // Decimation attempts
            if(decimations && decimations_count)
            {
                for(size_t di = 0; di < decimations_count; ++di)
                {
                    uint16_t eta = decimations[di];
                    if(eta == 0) continue;
                    for(uint16_t sigma = 0; sigma < eta; ++sigma)
                    {
                        uint8_t out_c[msg_length];
                        if(ApplyPrecomputedDecimation(src_ptr, (uint8_t)eta, sigma, out_c) == 0)
                        {
                            memcpy(dst_ptr, out_c, msg_length);
                            return 0;
                        }
                        if(DecodeDualBasisDecimation(src_ptr, dst_ptr, (uint8_t)eta, (uint16_t)sigma, 255) == 0)
                        {
                            return 0;
                        }

                        // heuristic sampling
                        uint8_t samp[msg_length];
                        for(uint8_t j = 0; j < msg_length; ++j)
                            samp[j] = src_ptr[(sigma + (uint32_t)j * eta) % 255];

                        // encode candidate and compare to input
                        Poly *msg_out = &polynoms[ID_MSG_OUT];
                        msg_out->Reset();
                        msg_out->Set(samp, msg_length);
                        msg_out->length = msg_length + ecc_length;
                        for(uint8_t i = 0; i < msg_length; i++)
                        {
                            uint8_t coef = msg_out->at(i);
                            if(coef != 0)
                            {
                                for(uint16_t j = 1; j < gen->length; j++)
                                    msg_out->at(i+j) ^= gf::mul(gen->at(j), coef);
                            }
                        }
                        // compare
                        uint16_t differ = 0;
                        for(uint8_t k = 0; k < n; ++k)
                        {
                            if(msg_out->ptr()[k] != src_ptr[k]) { if(++differ > ecc_length) break; }
                        }
                        if(differ <= ecc_length/2)
                        {
                            memcpy(dst_ptr, samp, msg_length);
                            return 0;
                        }
                    }
                }
            }

            // Fallback: try floating-window decode (systematic decode on rotated view)
            {
                Poly *msg_in  = &polynoms[ID_MSG_IN];
                Poly *msg_out = &polynoms[ID_MSG_OUT];
                Poly *epos    = &polynoms[ID_ERASURES];
                Poly *synd    = &polynoms[ID_SYNDROMES];
                Poly *forney  = &polynoms[ID_FORNEY];
                Poly *errloc  = &polynoms[ID_ERRORS_LOC];
                Poly *reloc   = &polynoms[ID_TPOLY3];
                Poly *err     = &polynoms[ID_ERRORS];

                uint8_t rotated_erase_buf[512];
                uint8_t reconstructed_code2[512];
                uint8_t corrected_code_orig[512];

                for(uint16_t shift = 0, tried = 0; shift < n && tried < max_shifts; shift += shift_step, ++tried)
                {
                    for(uint8_t i = 0; i < n; ++i)
                        rotated_code[i] = src_ptr[(i + shift) % n];

                    msg_in->Set(rotated_code, (uint16_t)msg_length);
                    msg_in->Set(rotated_code + msg_length, ecc_length, msg_length);
                    msg_out->Copy(msg_in);

                    if(erase_pos == NULL) epos->length = 0;
                    else {
                        epos->Set(erase_pos, (uint16_t)erase_count);
                        for(uint16_t i = 0; i < epos->length; i++) msg_in->at(epos->at((uint16_t)i)) = 0;
                    }
                    if(epos->length > ecc_length) continue;

                    CalcSyndromes(msg_in);
                    bool has_errors = false;
                    for(uint16_t i = 0; i < synd->length; i++) if(synd->at((uint16_t)i) != 0) { has_errors = true; break; }

                    if(!has_errors)
                    {
                        msg_out->length = (uint16_t)msg_length;
                        memcpy(rotated_code, msg_out->ptr(), msg_out->length * sizeof(uint8_t));
                    }
                    else
                    {
                        CalcForneySyndromes(synd, epos, n);
                        FindErrorLocator(forney, NULL, epos->length);

                        reloc->length = errloc->length;
                        for(int16_t i = errloc->length-1, j = 0; i >= 0; i--, j++) reloc->at((uint16_t)j) = errloc->at((uint16_t)i);

                        if(!FindErrors(reloc, n)) continue;
                        if(err->length == 0) continue;

                        for(uint16_t i = 0; i < err->length; i++) epos->Append(err->at((uint16_t)i));
                        CorrectErrata(synd, epos, msg_in);

                        msg_out->length = (uint16_t)msg_length;
                        memcpy(rotated_code, msg_out->ptr(), msg_out->length * sizeof(uint8_t));
                    }

                    // reconstruct full code and rotate back
                    msg_out->Reset();
                    msg_out->Set(rotated_code, msg_length);
                    msg_out->length = (uint16_t)(msg_length + ecc_length);
                    for(uint16_t i = 0; i < msg_length; i++)
                    {
                        uint8_t coef = msg_out->at((uint16_t)i);
                        if(coef != 0)
                        {
                            for(uint16_t j = 1; j < gen->length; j++)
                                msg_out->at((uint16_t)(i + j)) ^= gf::mul(gen->at((uint16_t)j), coef);
                        }
                    }
                    for(uint16_t i = 0; i < n; i++) reconstructed_code2[i] = msg_out->ptr()[i];
                    for(uint16_t i = 0; i < n; i++) corrected_code_orig[(i + shift) % n] = reconstructed_code2[i];

                    pPoly->Set(corrected_code_orig, n, 0);
                    PolyDivGetQuotient(pPoly, gen, quot);
                    uint8_t copy_len = (quot->length < msg_length) ? quot->length : msg_length;
                    memcpy(dst_ptr, quot->ptr(), copy_len);
                    if(copy_len < msg_length) memset(dst_ptr + copy_len, 0, msg_length - copy_len);
                    return 0;
                }
            }

            // nothing found
            return 1;
        }

        /**
         * Algebraic reconstruction for decimated subsequence using dual basis (Gaussian elimination).
         * Solves A * C = gamma where A[r][j] = (epsilon^{2^j})^{pos_r}.
         */
        int DecodeDualBasisDecimation(const void* src, void* dst, uint8_t eta, uint16_t sigma, size_t period = 255)
        {
            assert(msg_length > 0 && msg_length < 256);
            const uint8_t m = msg_length;

            const uint8_t* src_ptr = (const uint8_t*) src;
            uint8_t* dst_ptr = (uint8_t*) dst;

            if(eta == 0) return 1;
            if(period == 0) return 1;

            /* Precompute e_j = epsilon^{2^j} */
            uint8_t e[256];
            for(uint8_t j = 0; j < m; ++j)
            {
                int exp = 1;
                for(uint8_t t = 0; t < j; ++t) exp = (exp * 2) % 255;
                e[j] = gf::pow(2, exp);
            }

            /* Build augmented matrix mat[m][m+1] */
            uint8_t mat[256][257];
            memset(mat, 0, sizeof(mat));

            for(uint8_t r = 0; r < m; ++r)
            {
                size_t pos = (sigma + (size_t)r * eta) % period;
                for(uint8_t j = 0; j < m; ++j)
                {
                    mat[r][j] = gf::pow(e[j], (intmax_t)pos);
                }
                mat[r][m] = src_ptr[pos % period];
            }

            /* Gaussian elimination over GF to RREF */
            for(uint8_t col = 0, row = 0; col < m && row < m; ++col)
            {
                int pivot = -1;
                for(uint8_t r = row; r < m; ++r) if(mat[r][col] != 0) { pivot = r; break; }
                if(pivot == -1) return 1;
                if((uint8_t)pivot != row)
                {
                    for(uint8_t c = col; c <= m; ++c) { uint8_t tmp = mat[row][c]; mat[row][c] = mat[pivot][c]; mat[pivot][c] = tmp; }
                }
                uint8_t inv_pivot = gf::inverse(mat[row][col]);
                for(uint8_t c = col; c <= m; ++c) mat[row][c] = gf::mul(mat[row][c], inv_pivot);
                for(uint8_t r = 0; r < m; ++r)
                {
                    if(r == row) continue;
                    uint8_t factor = mat[r][col];
                    if(factor != 0)
                    {
                        for(uint8_t c = col; c <= m; ++c)
                        {
                            uint8_t prod = gf::mul(factor, mat[row][c]);
                            mat[r][c] ^= prod;
                        }
                    }
                }
                ++row;
            }

            for(uint8_t j = 0; j < m; ++j) dst_ptr[j] = mat[j][m];
            return 0;
        }

        /**
         * Precompute inverse matrices for decimations (A^{-1}) to speed reconstruction to O(m^2).
         */
        void PrecomputeDecimationMatrices(const uint8_t* etas, size_t etas_count, bool all_sigmas = true, size_t period = 255)
        {
            assert(msg_length > 0 && msg_length < 256);
            const uint8_t m = msg_length;

            precomp_matrices.clear();

            std::vector<uint8_t> e(m);
            for(uint8_t j = 0; j < m; ++j)
            {
                int exp = 1;
                for(uint8_t t = 0; t < j; ++t) exp = (exp * 2) % 255;
                e[j] = gf::pow(2, exp);
            }

            for(size_t i = 0; i < etas_count; ++i)
            {
                uint16_t eta = etas[i];
                if(eta == 0) continue;
                uint16_t sigma_max = all_sigmas ? eta : 1;
                for(uint16_t sigma = 0; sigma < sigma_max; ++sigma)
                {
                    std::vector<uint8_t> A(m * m);
                    for(uint8_t r = 0; r < m; ++r)
                    {
                        size_t pos = (sigma + (size_t)r * eta) % period;
                        for(uint8_t j = 0; j < m; ++j)
                        {
                            A[r * m + j] = gf::pow(e[j], (intmax_t)pos);
                        }
                    }

                    std::vector<uint8_t> inv(m * m);
                    if(invertMatrixGF(A.data(), inv.data(), m))
                    {
                        PrecompEntry ent;
                        ent.eta = (uint8_t)eta;
                        ent.sigma = (uint16_t)sigma;
                        ent.inv = std::move(inv);
                        precomp_matrices.push_back(std::move(ent));
                    }
                }
            }
        }

        /**
         * Apply precomputed inverse matrix: out = inv * gamma
         */
        int ApplyPrecomputedDecimation(const uint8_t* src, uint8_t eta, uint16_t sigma, uint8_t* out, size_t period = 255) const
        {
            assert(msg_length > 0 && msg_length < 256);
            const uint8_t m = msg_length;

            for(size_t i = 0; i < precomp_matrices.size(); ++i)
            {
                if(precomp_matrices[i].eta == eta && precomp_matrices[i].sigma == sigma)
                {
                    const std::vector<uint8_t>& inv = precomp_matrices[i].inv;
                    uint8_t gamma[256];
                    for(uint8_t r = 0; r < m; ++r)
                    {
                        size_t pos = (sigma + (size_t)r * eta) % period;
                        gamma[r] = src[pos % period];
                    }
                    for(uint8_t row = 0; row < m; ++row)
                    {
                        uint8_t accum = 0;
                        for(uint8_t col = 0; col < m; ++col)
                        {
                            uint8_t a = inv[row * m + col];
                            if(a != 0 && gamma[col] != 0)
                                accum ^= gf::mul(a, gamma[col]);
                        }
                        out[row] = accum;
                    }
                    return 0;
                }
            }
            return 1;
        }

#ifndef RS_DEBUG
        private:
#endif

        enum POLY_ID
        {
            ID_MSG_IN = 0,
            ID_MSG_OUT,
            ID_GENERATOR,
            ID_TPOLY1,
            ID_TPOLY2,

            ID_MSG_E,

            ID_TPOLY3,
            ID_TPOLY4,

            ID_SYNDROMES,
            ID_FORNEY,

            ID_ERASURES_LOC,
            ID_ERRORS_LOC,

            ID_ERASURES,
            ID_ERRORS,

            ID_COEF_POS,
            ID_ERR_EVAL
        };

        uint8_t* memory; ///< pointer to polynomial memory (stack)
        Poly polynoms[MSG_CNT + POLY_CNT]; ///< polynomial workspace

        /* Construct generator polynomial g(x) */
        void GeneratorPoly()
        {
            Poly *gen = polynoms + ID_GENERATOR;
            gen->at(0) = 1;
            gen->length = 1;

            Poly *mulp = polynoms + ID_TPOLY1;
            Poly *temp = polynoms + ID_TPOLY2;
            mulp->length = 2;

            for(int16_t i = 0; i < ecc_length; i++)
            {
                mulp->at(0) = 1;
                mulp->at(1) = gf::pow(2, i);

                gf::poly_mul(gen, mulp, temp);
                gen->Copy(temp);
            }
        }

        void CalcSyndromes(const Poly *msg)
        {
            Poly *synd = &polynoms[ID_SYNDROMES];
            synd->length = (uint16_t)(ecc_length+1);
            synd->at(0) = 0;
            for(uint16_t i = 1; i < ecc_length+1; i++)
            {
                synd->at((uint16_t)i) = gf::poly_eval(msg, gf::pow(2, (uint8_t)(i-1)));
            }
        }

        void FindErrataLocator(const Poly *epos)
        {
            Poly *errata_loc = &polynoms[ID_ERASURES_LOC];
            Poly *mulp = &polynoms[ID_TPOLY1];
            Poly *addp = &polynoms[ID_TPOLY2];
            Poly *apol = &polynoms[ID_TPOLY3];
            Poly *temp = &polynoms[ID_TPOLY4];

            errata_loc->length = 1;
            errata_loc->at(0) = 1;

            mulp->length = 1;
            addp->length = 2;

            for(uint16_t i = 0; i < epos->length; i++)
            {
                mulp->at(0) = 1;
                addp->at(0) = gf::pow(2, epos->at((uint16_t)i));
                addp->at(1) = 0;

                gf::poly_add(mulp, addp, apol);
                gf::poly_mul(errata_loc, apol, temp);

                errata_loc->Copy(temp);
            }
        }

        void FindErrorEvaluator(const Poly *synd, const Poly *errata_loc, Poly *dst, uint8_t ecclen)
        {
            Poly *mulp = &polynoms[ID_TPOLY1];
            gf::poly_mul(synd, errata_loc, mulp);

            Poly *divisor = &polynoms[ID_TPOLY2];
            divisor->length = (uint16_t)(ecclen + 2);

            divisor->Reset();
            divisor->at(0) = 1;

            gf::poly_div(mulp, divisor, dst);
        }

        void CorrectErrata(const Poly *synd, const Poly *err_pos, const Poly *msg_in)
        {
            Poly *c_pos = &polynoms[ID_COEF_POS];
            Poly *corrected = &polynoms[ID_MSG_OUT];
            c_pos->length = err_pos->length;

            for(uint16_t i = 0; i < err_pos->length; i++)
                c_pos->at((uint16_t)i) = (uint8_t)(msg_in->length - 1 - err_pos->at((uint16_t)i));

            FindErrataLocator(c_pos);
            Poly *errata_loc = &polynoms[ID_ERASURES_LOC];

            Poly *rsynd = &polynoms[ID_TPOLY3];
            rsynd->length = synd->length;

            for(int16_t i = (int16_t)synd->length - 1, j = 0; i >= 0; i--, j++)
                rsynd->at((uint16_t)j) = synd->at((uint16_t)i);

            Poly *re_eval = &polynoms[ID_TPOLY4];
            FindErrorEvaluator(rsynd, errata_loc, re_eval, (uint8_t)(errata_loc->length - 1));

            Poly *e_eval = &polynoms[ID_ERR_EVAL];
            e_eval->length = re_eval->length;
            for(int16_t i = (int16_t)re_eval->length - 1, j = 0; i >= 0; i--, j++)
                e_eval->at((uint16_t)j) = re_eval->at((uint16_t)i);

            Poly *X = &polynoms[ID_TPOLY1];
            X->length = 0;

            for(uint16_t i = 0; i < c_pos->length; i++)
                X->Append(gf::pow(2, (uint8_t)(255 - c_pos->at((uint16_t)i))));

            Poly *E = &polynoms[ID_MSG_E];
            E->Reset();
            E->length = msg_in->length;

            Poly *err_loc_prime_temp = &polynoms[ID_TPOLY2];

            for(uint16_t i = 0; i < X->length; i++)
            {
                uint8_t Xi_inv = gf::inverse(X->at((uint16_t)i));

                err_loc_prime_temp->length = 0;
                for(uint16_t j = 0; j < X->length; j++)
                {
                    if(j != i)
                        err_loc_prime_temp->Append(gf::sub(1, gf::mul(Xi_inv, X->at((uint16_t)j))));
                }

                uint8_t err_loc_prime = 1;
                for(uint16_t j = 0; j < err_loc_prime_temp->length; j++)
                    err_loc_prime = gf::mul(err_loc_prime, err_loc_prime_temp->at((uint16_t)j));

                uint8_t y = gf::poly_eval(re_eval, Xi_inv);
                y = gf::mul(gf::pow(X->at((uint16_t)i), 1), y);

                E->at((uint16_t)err_pos->at((uint16_t)i)) = gf::div(y, err_loc_prime);
            }

            gf::poly_add(msg_in, E, corrected);
        }

        bool FindErrorLocator(const Poly *synd, Poly *erase_loc = NULL, size_t erase_count = 0)
        {
            Poly *error_loc = &polynoms[ID_ERRORS_LOC];
            Poly *err_loc   = &polynoms[ID_TPOLY1];
            Poly *old_loc   = &polynoms[ID_TPOLY2];
            Poly *temp      = &polynoms[ID_TPOLY3];
            Poly *temp2     = &polynoms[ID_TPOLY4];

            if(erase_loc != NULL)
            {
                err_loc->Copy(erase_loc);
                old_loc->Copy(erase_loc);
            }
            else
            {
                err_loc->length = 1;
                old_loc->length = 1;
                err_loc->at(0)  = 1;
                old_loc->at(0)  = 1;
            }

            uint16_t synd_shift = 0;
            if(synd->length > ecc_length)
                synd_shift = (uint16_t)(synd->length - ecc_length);

            for(uint16_t i = 0; i < (uint16_t)(ecc_length - erase_count); i++)
            {
                uint16_t K;
                if(erase_loc != NULL)
                    K = (uint16_t)(erase_count + i + synd_shift);
                else
                    K = (uint16_t)(i + synd_shift);

                uint8_t delta = synd->at(K);
                for(uint16_t j = 1; j < err_loc->length; j++)
                {
                    uint16_t index = (uint16_t)(err_loc->length - j - 1);
                    delta ^= gf::mul(err_loc->at(index), synd->at((uint16_t)(K-j)));
                }

                old_loc->Append(0);

                if(delta != 0)
                {
                    if(old_loc->length > err_loc->length)
                    {
                        gf::poly_scale(old_loc, temp, delta);
                        gf::poly_scale(err_loc, old_loc, gf::inverse(delta));
                        err_loc->Copy(temp);
                    }

                    gf::poly_scale(old_loc, temp, delta);
                    gf::poly_add(err_loc, temp, temp2);
                    err_loc->Copy(temp2);
                }
            }

            uint32_t shift = 0;
            while(err_loc->length && err_loc->at((uint16_t)shift) == 0) ++shift;

            uint32_t errs = err_loc->length - shift - 1;
            if(((errs - erase_count) * 2 + erase_count) > ecc_length)
                return false;

            memcpy(error_loc->ptr(), err_loc->ptr() + shift, (err_loc->length - shift) * sizeof(uint8_t));
            error_loc->length = (uint16_t)(err_loc->length - shift);
            return true;
        }

        bool FindErrors(const Poly *error_loc, size_t msg_in_size)
        {
            Poly *err = &polynoms[ID_ERRORS];

            uint16_t errs = (uint16_t)(error_loc->length - 1);
            err->length = 0;

            for(uint16_t i = 0; i < msg_in_size; i++)
            {
                if(gf::poly_eval(error_loc, gf::pow(2, (uint8_t)i)) == 0)
                {
                    err->Append((uint8_t)(msg_in_size - 1 - i));
                }
            }

            if(err->length != errs) return false;
            return true;
        }

        void CalcForneySyndromes(const Poly *synd, const Poly *erasures_pos, size_t msg_in_size)
        {
            Poly *erase_pos_reversed = &polynoms[ID_TPOLY1];
            Poly *forney_synd = &polynoms[ID_FORNEY];
            erase_pos_reversed->length = 0;

            for(uint16_t i = 0; i < erasures_pos->length; i++)
                erase_pos_reversed->Append((uint8_t)(msg_in_size - 1 - erasures_pos->at((uint16_t)i)));

            forney_synd->Reset();
            forney_synd->Set(synd->ptr()+1, (uint16_t)(synd->length - 1));

            for(uint16_t i = 0; i < erasures_pos->length; i++)
            {
                uint8_t x = gf::pow(2, erase_pos_reversed->at((uint16_t)i));
                for(int16_t j = 0; j < (int16_t)forney_synd->length - 1; j++)
                {
                    forney_synd->at((uint16_t)j) = gf::mul(forney_synd->at((uint16_t)j), x) ^ forney_synd->at((uint16_t)(j+1));
                }
            }
        }

        /**
         * Compute quotient of polynomial division p / q, store quotient in dst (dst->ptr()).
         */
        void PolyDivGetQuotient(const Poly *p, const Poly *q, Poly *dst)
        {
            Poly *work = &polynoms[ID_MSG_OUT];

            if(work->ptr() != p->ptr())
            {
                memcpy(work->ptr(), p->ptr(), p->length * sizeof(uint8_t));
            }
            work->length = p->length;

            uint8_t qlen = q->length;
            int out_len = (int)p->length - (int)qlen + 1;
            if(out_len < 0) out_len = 0;

            dst->length = (uint8_t)out_len;

            for(int i = 0; i < out_len; ++i)
            {
                uint8_t coef = work->at((uint16_t)i);
                if(q->at(0) != 1)
                    coef = gf::div(coef, q->at(0));
                dst->at((uint16_t)i) = coef;

                if(coef != 0)
                {
                    for(uint16_t j = 1; j < qlen; j++)
                    {
                        if(q->at((uint16_t)j) != 0 && (i + j) < work->length)
                            work->at((uint16_t)(i + j)) ^= gf::mul(q->at((uint16_t)j), coef);
                    }
                }
            }
        }

        /**
         * Invert m x m matrix A over GF(2^8); returns true if invertible and writes inverse to inv_out (row-major).
         */
        bool invertMatrixGF(const uint8_t* A_in, uint8_t* inv_out, uint8_t m)
        {
            uint8_t aug[256][512];
            memset(aug, 0, sizeof(aug));
            for(uint8_t r = 0; r < m; ++r)
            {
                for(uint8_t c = 0; c < m; ++c)
                    aug[r][c] = A_in[r * m + c];
                for(uint8_t c = 0; c < m; ++c)
                    aug[r][m + c] = (r == c) ? 1 : 0;
            }

            for(uint8_t col = 0, row = 0; col < m && row < m; ++col)
            {
                int pivot = -1;
                for(uint8_t r = row; r < m; ++r) if(aug[r][col] != 0) { pivot = r; break; }
                if(pivot == -1) return false;
                if((uint8_t)pivot != row)
                {
                    for(uint8_t c = 0; c < 2*m; ++c) { uint8_t tmp = aug[row][c]; aug[row][c] = aug[pivot][c]; aug[pivot][c] = tmp; }
                }
                uint8_t inv_pivot = gf::inverse(aug[row][col]);
                for(uint8_t c = 0; c < 2*m; ++c) aug[row][c] = gf::mul(aug[row][c], inv_pivot);

                for(uint8_t r = 0; r < m; ++r)
                {
                    if(r == row) continue;
                    uint8_t factor = aug[r][col];
                    if(factor != 0)
                    {
                        for(uint8_t c = 0; c < 2*m; ++c)
                        {
                            uint8_t prod = gf::mul(factor, aug[row][c]);
                            aug[r][c] ^= prod;
                        }
                    }
                }
                ++row;
            }

            for(uint8_t r = 0; r < m; ++r)
            {
                for(uint8_t c = 0; c < m; ++c)
                    inv_out[r * m + c] = aug[r][m + c];
            }
            return true;
        }

        struct PrecompEntry
        {
            uint8_t eta;
            uint16_t sigma;
            std::vector<uint8_t> inv; // inverse matrix row-major
        };

        std::vector<PrecompEntry> precomp_matrices;

    };

} // namespace RS

#endif // RS_HPP
