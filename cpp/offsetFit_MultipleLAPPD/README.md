# offsetFit_MultipleLAPPD

This script performs LAPPD timing offset fitting and takes `LAPPDTree.root` (produced by the `LAPPD_EB` toolchain) as input.

## Requirements

* All LAPPD events must appear in sequence.
* All CTC events must appear in sequence.

## Procedure

1. Load `LAPPDTree.root`.
2. Load the `TimeStamp` tree and obtain the run number and part file number.
3. For each unique run number and part file number, perform the fit.
4. For each LAPPD ID:

   a. Loop over the `TimeStamp` tree and load:

      * Data timestamp
      * Data beamgate
      * PPS for Board 0 and Board 1
        
   b. Loop over the `Trig` (or `GTrig`) tree and load:

      * CTC PPS (Trigger Word = 32)
      * Target trigger word
        
5. Find the number of resets in the LAPPD PPS stream:

   * There must be at least one event in each reset.
   * Based on the event index order, only fit before the reset that does not contain a data event.
   * Save the data event and PPS index ordering.
     
6. Use the target trigger word to perform the offset fit.
7. Save the offset for:

   * LAPPD ID
   * Run Number
   * Part File Number
   * Event Index
   * Reset Number
     
8. Print information to output files.

## Input Arguments

The script takes five input arguments:

| Argument               | Description                                                  |
| ---------------------- | ------------------------------------------------------------ |
| `filename`             | Input file (`LAPPDTree.root`)                                |
| `fitTargetTriggerWord` | Undelayed Beam Trigger (14) for beam runs, 47 for laser runs |
| `triggerGrouped`       | 1 for beam runs, 0 for laser runs                            |
| `intervalInSeconds`    | See table below                                              |
| `processPfNumber`      | 0                                                            |

### intervalInSeconds

| Run Range                | Value |
| ------------------------ | ----- |
| Run < 5140               | 10    |
| Run = 6139 (Beam run)    | 10    |
| Run ≥ 5140 (except 6139) | 5     |

## Usage

### Beam Runs

#### Runs < 5140

```bash
root -l -q 'offsetFit_MultipleLAPPD.cpp("LAPPDTree.root", 14, 1, 10, 0)'
```

#### Runs ≥ 5140 (except Run 6139)

```bash
root -l -q 'offsetFit_MultipleLAPPD.cpp("LAPPDTree.root", 14, 1, 5, 0)'
```

#### Run 6139

```bash
root -l -q 'offsetFit_MultipleLAPPD.cpp("LAPPDTree.root", 14, 1, 10, 0)'
```

### Laser Runs

#### Runs < 5140

```bash
root -l -q 'offsetFit_MultipleLAPPD.cpp("LAPPDTree.root", 47, 0, 10, 0)'
```

#### Runs ≥ 5140 (except Run 6139)

```bash
root -l -q 'offsetFit_MultipleLAPPD.cpp("LAPPDTree.root", 47, 0, 5, 0)'
```
