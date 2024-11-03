# Reed-Solomon
Reed Solomon encoder and decoder library

## Overview

This RS implementation was designed for embedded purposes, so all the memory allocations performed on the stack.

## Build

There is no need in building RS library, cause all the implementation is in headers.

## Usage

You just need to include header <b>rs.hpp</b>

Template class ReedSolomon accepts two template arguments: message length and ecc length. <br>

Simple example: <br>
```

    char initialData[] = "There you store information which will be transefered and may need error correcting";
    const int dataLength = sizeof(message);
    const int parityLength = sizeof(message);
    
    char repairedData[dataLength];
    char encodedData[dataLength + parityLength];

    RS::ReedSolomon<dataLength, parityLength> coder;

    coder.Encode(initialData, encodedData);

    /// Corrupting bytes both of message and parity
    for(uint i = 0; i < ecclen / 4; i++) 
    {
        encodedData[rand() % dataLength] = 'E';
        encodedData[rand() % dataLength + rand() % parityLength] = 'E';
    }

    coder.Decode(encodedData, repairedData);

    std::cout << "Original:  " << initialData  << std::endl;
    std::cout << "Corrupted: " << encodedData  << std::endl;
    std::cout << "Repaired:  " << repairedData << std::endl;

```
