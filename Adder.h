#pragma once

#include <Cell.h>

namespace rscoder::Encoder {
    class Adder {
    public:
        constexpr static Byte Summarize(const std::vector<Cell>& cells) noexcept {
            Byte result = 0;
            for (const auto& cell : cells) {
                result ^= cell.GetOutput();
            }
            return result;
        }
    };
} // namespace rscoder::Encoder
