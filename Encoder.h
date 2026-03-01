#pragma once

#include <Adder.h>

#include <execution>
#include <numeric>

namespace rscoder::Encoder {
    class Encoder {
    private:
        size_t _numberOfCells;
        size_t _maxOffset;
        std::vector<std::vector<Byte>> _splitBuffer;
        std::vector<std::vector<Byte>> _encodedSequence;

        constexpr std::vector<std::vector<Byte>> SplitReadBuffer(const std::vector<Byte>& byteArray) const {
            std::vector<std::vector<Byte>> splitReadBuffer;
            splitReadBuffer.assign(_maxOffset, std::vector<Byte>(_numberOfCells, 0));

            for (size_t currentByte = 0; currentByte < byteArray.size(); ++currentByte) {
                size_t row = currentByte / _numberOfCells;
                size_t col = currentByte % _numberOfCells;
                splitReadBuffer[row][col] = byteArray[currentByte];
            }
            return splitReadBuffer;
        }

    public:
        explicit Encoder(size_t numberOfCells)
                : _numberOfCells(numberOfCells), _maxOffset(0) {
        }

        constexpr void Read(const std::vector<Byte>& byteArray) {
            _maxOffset = (byteArray.size() + _numberOfCells - 1) / _numberOfCells;
            _encodedSequence.assign(_maxOffset, std::vector<Byte>(255, 0));
            _splitBuffer = SplitReadBuffer(byteArray);
        }

        void Encode() {
            std::vector<size_t> offsets(_maxOffset);
            std::iota(offsets.begin(), offsets.end(), 0);
            
            std::for_each(std::execution::par_unseq, offsets.begin(), offsets.end(),
                [this](size_t offset) {
                    thread_local std::vector<Cell> localCells;
                    
                    if (localCells.size() != this->_numberOfCells) {
                        localCells.clear();
                        for (size_t i = 0; i < this->_numberOfCells; ++i) {
                            localCells.emplace_back(static_cast<Byte>(i));
                        }
                    }
                    
                    for (size_t currentByte = 0; currentByte < 255; ++currentByte) {
                        for (size_t currentCell = 0; currentCell < this->_numberOfCells; ++currentCell) {
                            if (!localCells[currentCell].IsInitialized()) {
                                localCells[currentCell].StartHinge(this->_splitBuffer[offset][currentCell]);
                            } else {
                                localCells[currentCell].Step();
                            }
                        }
                        
                        this->_encodedSequence[offset][currentByte] = Adder::Summarize(localCells);
                    }
                    
                    for (auto& cell : localCells) {
                        cell.Reset();
                    }
                });
        }

        [[nodiscard]] const std::vector<std::vector<Byte>>& GetEncodedSequence() const {
            return _encodedSequence;
        }

        [[nodiscard]] std::vector<std::vector<Byte>> Encode(const std::vector<Byte>& byteArray) {
            Read(byteArray);
            Encode();

            const size_t rows = _encodedSequence.size();
            const size_t cols = _encodedSequence[0].size();

            std::vector<std::vector<Byte>> realEncodedSequence(cols);

            for (auto& encodedPacket : realEncodedSequence) {
                encodedPacket.reserve(rows);
            }

            for (size_t row = 0; row < rows; ++row) {
                const auto& encodedSequenceRow = _encodedSequence[row];

                for (size_t col = 0; col < cols; ++col) {
                    realEncodedSequence[col].push_back(encodedSequenceRow[col]);
                }
            }

            return realEncodedSequence;
        }
    };
} // namespace rscoder::Encoder
