# RS Coder

Header-only C++17 library for Reed-Solomon erasure coding with **decimation** recovery.  
It splits input data into blocks, encodes each block into 255 symbols, and arranges them into 255 packets.  
If any 8 packets whose indices form an arithmetic progression with step 1, 2, 4, … are received,  
the original block can be decoded - ideal for lossy networks (e.g., UAV video streaming).

## Features

- **GF(256) operations** - precomputed `Epsilon` (powers of generator) and `Lambda` (log tables).
- **Parallel encoding** - uses `std::execution::par_unseq` (C++17).
- **Decimation decoding** - recover data from any 8 packets with indices `start + i * 2^d`.
- **Header‑only** - just copy the `.h` files and `#include` them.

## Requirements

- C++17 compiler (with `<execution>` support for parallel encoding)
- Standard C++ library

## Integration

Copy all `.h` files into your project and include the required ones:

```cpp
#include "Encoder.h"
#include "Decoder.h"
#include "Decimation.h"
```

## Usage Overview

- **Encoder** (`rscoder::Encoder::Encoder`) - encodes a byte vector into 255 packets.
- **Decoder free functions** (`rscoder::Decode`) - decode when all packets or a specific decimation are available.
- **Decimation::Decoder** - accumulates packets and automatically recovers blocks when a valid decimation appears.

## Example

The following example encodes a long string, takes the first 8 packets, and uses `Decimation::Decoder` to recover the original message.

```cpp
#include <Decimation.h>
#include <Encoder.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

int main() {
    std::string msg = "When managing first-person view unmanned aerial vehicles in scenarios involving hybrid "
                      "terrestrial-orbital communication networks for data transmission, a critical challenge is "
                      "the low intensity of video streams caused by the limited bandwidth of such communication "
                      "channels and long frames decoding time with neural codecs. In order to increase the video "
                      "stream intensity, frames interpolation can be used. In this article, five interpolation "
                      "techniques to enhance the intensity of FPV video streams are reviewed: linear "
                      "interpolation, morphing, optical flow, unidirectional optical flow, and spline "
                      "interpolation. A comparative analysis of their effectiveness is conducted using an author "
                      "dataset. The research highlights optical flow's superior quality and linear "
                      "interpolation's high processing speed. The work concludes with practical recommendations "
                      "for selecting optimal methods for UAVs operating under low network bandwidth and "
                      "constrained computational resources. It applies to the video stream intensity controls.";

    std::vector<uint8_t> byteArray(msg.size());
    std::transform(msg.begin(), msg.end(), byteArray.begin(),
                   [](unsigned char c) { return static_cast<uint8_t>(c); });

    rscoder::Encoder::Encoder encoder(8);
    rscoder::Decimation::Decoder decoder;

    auto packets = encoder.Encode(byteArray);
    std::vector<uint8_t> decodedData;

    for (int i = 0; i < 8; ++i) {
        auto result = decoder.Decode(0, i, std::move(packets[i]));

        if (result) {
            decodedData = std::move(*result);
        }
    }

    for (uint8_t b : decodedData) {
        std::cout << static_cast<char>(b);
    }
    std::cout << std::endl;

    return 0;
}
```
