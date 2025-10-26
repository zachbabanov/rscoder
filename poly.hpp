#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include <string.h>

#if !defined RS_DEBUG && !defined __CC_ARM && !defined RS_NO_ASSERT
#include <assert.h>
#else
#define assert(dummy)
#endif

namespace RS
{

    struct Poly
    {
        Poly() : length(0), _memory(NULL) { }

        Poly(uint8_t id, uint16_t offset, uint16_t size) : length(0), _id(id), _size(size), _offset(offset), _memory(NULL) { }

        /// Append number at the end of polynomial
        /// \param num - number to append
        /// \return false if polynomial can't be stretched
        inline bool Append(uint8_t num)
        {
            assert((uint16_t)(length + 1) <= _size);
            ptr()[length++] = num;
            return true;
        }

        /// Polynomial initialization
        inline void Init(uint8_t id, uint16_t offset, uint16_t size, uint8_t** memory_ptr)
        {
            this->_id     = id;
            this->_offset = offset;
            this->_size   = size;
            this->length  = 0;
            this->_memory = memory_ptr;
        }

        /// Polynomial memory zeroing
        inline void Reset()
        {
            memset((void*)ptr(), 0, this->_size);
        }

        /// Copy polynomial to memory
        /// \param src - source byte-sequence
        /// \param len - size of polynomial
        /// \param offset - write offset
        inline void Set(const uint8_t* src, uint16_t len, uint16_t offset = 0)
        {
            assert(src && (uint16_t)len <= (uint16_t)(_size - offset));
            memcpy(ptr()+offset, src, len * sizeof(uint8_t));
            length = (uint16_t)(len + offset);
        }

#define poly_max(a, b) ((a > b) ? (a) : (b))

        inline void Copy(const Poly* src)
        {
            length = poly_max(length, src->length);
            Set(src->ptr(), length);
        }

        inline uint8_t& at(uint16_t i) const
        {
            assert(i < _size);
            return ptr()[i];
        }

        inline uint8_t id() const
        {
            return _id;
        }

        inline uint16_t size() const
        {
            return _size;
        }

        /// Returns pointer to memory of this polynomial
        inline uint8_t* ptr() const
        {
            assert(_memory && *_memory);
            return (*_memory) + _offset;
        }

        uint16_t length;

    protected:

        uint8_t   _id;
        uint16_t  _size;    ///< Size of reserved memory for this polynomial (was uint8_t, now uint16_t)
        uint16_t  _offset;  ///< Offset in memory
        uint8_t** _memory;  ///< Pointer to pointer to memory
    };
}

#endif // POLY_H
