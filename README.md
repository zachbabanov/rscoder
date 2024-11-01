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

    char message[] = "Some very important message ought to be delivered";
    const int msglen = sizeof(message);
    const int ecclen = 8;
    
    char repaired[msglen];
    char encoded[msglen + ecclen];

    RS::ReedSolomon<msglen, ecclen> rs;

    rs.Encode(message, encoded);

    // Corrupting first 8 bytes of message (any 4 bytes can be repaired, no more)
    for(uint i = 0; i < ecclen / 2; i++) 
    {
        encoded[i] = 'E';
    }

    rs.Decode(encoded, repaired);

    std::cout << "Original:  " << message  << std::endl;
    std::cout << "Corrupted: " << encoded  << std::endl;
    std::cout << "Repaired:  " << repaired << std::endl;

```
