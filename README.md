---
Title: README.md
Author: Makeda Harris Morris (30113200, makeda.harris@ucalgary.ca)
Date: 05 December 2023
---

# Convolution Reverb Program

The goal of this assignment is to optimize, in stages, the performance of a convolution reverb
program. Convolution reverb is an audio digital signal processing technique where a “dry”
recording of an instrument (i.e., a recording without reverberation) is convolved with the impulse
response of an acoustical space such as a concert hall. The result of the convolution is a sound file
where the instrument sounds as if it were playing in the actual concert hall. This is a commonly
used, but computationally intensive, technique for adding natural-sounding reverberation to
recorded sounds.

## Usage

```bash
convolve inputfile IRfile outputfile
```

- **convolve:** A command-line program that takes in a dry recording and an impulse response,
and produces the convolved signal.
- **inputfile:** The dry recording in .wav format (monophonic, 16-bit, 44.1 kHz).
- **IRfile:** The impulse response in .wav format (monophonic, 16-bit, 44.1 kHz).
- **outputfile:** The convolved signal in .wav format.

### Example Usage

```bash
convolve dry_recording.wav impulse_response.wav convolved_output.wav
```

This example convolves the `dry_recording.wav` file with the `impulse_response.wav` file, producing the `convolved_output.wav`.

## Building the Program

To build the program, use the following command:

```bash
make
```

This command will also clean the program before building.

## Debugging the Program

To build the program with debugging symbols, use the following command:

```bash
make debug
```

## Cleaning the Program

To clean the program (remove obj and bin directories), use the following command:

```bash
make clean
```

## Requirements

- C++ compiler (g++)
- make utility

## Notes

- All .wav files should be monophonic, 16-bit, 44.1 kHz sound files.