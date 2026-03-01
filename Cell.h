#pragma once

#include <Epsilon.h>

#include <iostream>
#include <vector>
#include <cmath>

namespace rscoder::Encoder {
    class Cell {
        Byte _output = 0;
        Byte _position = 0;
        Byte _shiftValue = 0;
        bool _initialized = false;

    public:
        explicit Cell() = default;
        explicit Cell(Byte pos) : _position(pos), _shiftValue((1 << pos) % 255) {}

        constexpr void StartHinge(Byte data) noexcept {
            _output = data;
            _initialized = true;
        }

        constexpr void Step() {
            const int index = (IndexOf(_output) + _shiftValue) % 255;
            _output = Epsilon[index];
        }

        [[nodiscard]] constexpr Byte GetOutput() const noexcept { return _output; }
        constexpr void Reset() noexcept { _initialized = false; _output = 0; }
        [[nodiscard]]  constexpr bool IsInitialized() const noexcept { return _initialized; }
    };
} // namespace rscoder::Encoder
