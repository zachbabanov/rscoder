#pragma once

#include <Epsilon.h>

#include <algorithm>
#include <iterator>
#include <vector>

namespace rscoder {
    constexpr Byte Decode(const std::vector<Byte> &encodedPacket, int bytePosition) {
        auto summ = static_cast<Byte>(0);
#pragma unroll
        for (int currentByte = 0; currentByte < 8; ++currentByte) {
            const Byte byte = encodedPacket[currentByte];
            if (byte != 0) {
                const int lambdaIndex = Lambda[(currentByte + (bytePosition << 3)) & 0x3F];
                const int index = (IndexOf(byte) + lambdaIndex) % 255;
                summ ^= Epsilon[index];
            }
        }
        return summ;
    }

    constexpr static std::vector<Byte> Decode(std::vector<std::vector<Byte>> &&encodedSequence) {
        size_t numRows = encodedSequence.size();
        std::vector<Byte> decodedData(8 * encodedSequence.size());

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < 8; j++) {
                decodedData[(i << 3) + j] = Decode(encodedSequence[i], j);
            }
        }

        return decodedData;
    }

    constexpr Byte Decode(const std::vector<Byte> &encodedPacket, uint32_t bytePosition, int packetPosition) {
        auto summ = static_cast<Byte>(0);
        const int shift = ((packetPosition % 255 + 255) % 255) * (1 << bytePosition);
#pragma unroll
        for (int currentByte = 0; currentByte < 8; ++currentByte) {
            const Byte byte = encodedPacket[currentByte];
            if (byte != 0) {
                const int lambdaIndex = Lambda[(currentByte + (bytePosition << 3)) & 0x3F];
                int index = (IndexOf(byte) + lambdaIndex - shift) % 255;
                if (index < 0) {
                    index += 255;
                }
                summ ^= Epsilon[index];
            }
        }
        return summ;
    }

    constexpr static std::vector<Byte> Decode(std::vector<std::vector<Byte>> &&decimationPackets, int packetPosition, int decimationIndex) {
        std::vector<Byte> encodedSequence;
        encodedSequence.reserve(decimationPackets.size());
        std::vector<Byte> decodedData(decimationPackets.size() * decimationPackets[0].size());

        for (size_t decodedOctet = 0; decodedOctet < decimationPackets[0].size(); ++decodedOctet) {
            encodedSequence.clear();
            std::transform(decimationPackets.begin(), decimationPackets.end(),
                           std::back_inserter(encodedSequence),
                           [decodedOctet](const std::vector<Byte>& v) { return v[decodedOctet]; });

            for (int currentByte = 0; currentByte < 8; ++currentByte) {
                decodedData[(decodedOctet << 3) + ((currentByte - decimationIndex) & 7)] = Decode(encodedSequence, currentByte, packetPosition);
            }
        }

        return decodedData;
    }
} // Decoder