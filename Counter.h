#pragma once

#include <Decoder.h>

#include <optional>
#include <bitset>
#include <map>

namespace rscoder {
    struct CheckResult {
        uint8_t firstPacketIndex{};
        uint8_t decimationIndex{};
    };

    class Counter {
    public:
        std::optional<std::vector<uint8_t>> Check(uint32_t packetIndex, std::vector<uint8_t> &&packet) {
            _packets.emplace(packetIndex, std::move(packet));

            auto isValidDecimation = CheckDecimation();
            if (isValidDecimation) {
                std::vector<std::vector<uint8_t>> packets;
                const uint8_t difference = 1 << isValidDecimation->decimationIndex;

                for (int i = 0; i < 8; i++) {
                    packets.push_back(_packets.at(isValidDecimation->firstPacketIndex + i * difference));
                }

                return rscoder::Decode(std::move(packets), isValidDecimation->firstPacketIndex, isValidDecimation->decimationIndex);
            }

            return std::nullopt;
        }

    private:
        std::optional<CheckResult> CheckDecimation() {
            if (_packets.size() < 8) {
                return std::nullopt;
            }

            std::bitset<256> packetsPresent;
            std::vector<uint8_t> packetsIndexes;
            for (const auto &pair : _packets) {
                packetsPresent.set(pair.first);
                packetsIndexes.push_back(pair.first);
            }

            for (uint8_t startIndex : packetsIndexes) {
                for (uint8_t decimationCount = 0; decimationCount < 9; decimationCount++) {
                    const uint8_t difference = 1 << decimationCount;
                    bool isValidDecimation = true;

                    for (uint8_t i = 0; i < 8; ++i) {
                        uint8_t packetIndex = (startIndex + i * difference) % 255;
                        if (!packetsPresent.test(packetIndex)) {
                            isValidDecimation = false;
                            break;
                        }
                    }

                    if (isValidDecimation) {
                        return CheckResult{startIndex, decimationCount};
                    }
                }
            }

            return std::nullopt;
        }

        std::map<uint8_t, std::vector<uint8_t>> _packets;
    };
}