
# WebUSB and WebAssembly RTL-SDR Example
# Copyright 2019 Ahmet Inan <inan@aicodix.de>

CXXFLAGS = -ffreestanding -fvisibility=hidden --target=wasm32
LDFLAGS = -nostdlib -Wl,--allow-undefined,--export-dynamic,--export=__wasm_call_ctors

CXXFLAGS += -std=c++11 -W -Wall -Ofast -fno-exceptions -fno-rtti
CXXFLAGS += -I$(HOME)/dsp

example.wasm.gz: example.wasm
	gzip -f -n example.wasm

example.wasm: example.cc example.hh
	clang++ $(CXXFLAGS) $(LDFLAGS) -o example.wasm example.cc

.PHONY: clean

clean:
	rm -f example.wasm

