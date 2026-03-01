#pragma once

#include <Counter.h>

#include <set>

namespace rscoder::Decimation {
    class Decoder {
    public:
        std::optional<std::vector<uint8_t>> Decode(uint32_t blockIndex, uint32_t packetIndex, std::vector<uint8_t> &&packet) {
            if (_decodedBlocks.contains(blockIndex)) {
                return std::nullopt;
            }

            if (!_decimationsCounters.contains(blockIndex)) {
                _decimationsCounters.emplace(blockIndex, Counter());
            }

            auto isDecoded = _decimationsCounters.find(blockIndex)->second.Check(packetIndex, std::move(packet));

            if (isDecoded) {
                _decimationsCounters.erase(blockIndex);
                _decodedBlocks.emplace(blockIndex);
            }

            return isDecoded;
        }

    private:
        std::unordered_map<uint32_t, Counter> _decimationsCounters;
        std::set<uint32_t> _decodedBlocks;
    };
};