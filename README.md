# PASTA
Pattern Analysis vs Sequence-based Tree Alignment

This Python 3 package was developed for aligning sequences derived from binary tree structures, specifically but not limited to neuronal trees. This package, and results from alignment using it were first published in:

Gillette TA, Hosseini P, Ascoli GA (2015) Topological characterization of neuronal arbor morphology via sequence representation. II. Global alignment. BMC Bioinformatics (submitted).

The main files are:
- spaghetti.py: Pairwise global alignment between sequences in a single fasta file, or between sequences in two separate fasta files. The alignment is based on the Needleman-Wunsch algorithm modified to respect three specific bifurcation types to which each character must correspond. Matching and gapping rules are modified based on these bifurcation types. A second mode enables conversion from raw scores to per-character scores based on the shorter sequence length.
- penne.py: Multiple alignment utilizing the same underlying methods as in spaghetti. The process by default iterates using a position-specific score matrix (PSSM) based on the initial alignment.
- ditalini.py: A small program which produces measures and statistics of a multiple alignment.
- orzo.py: Extracts potential domains from a multiple alignment. This module was not used in the aforementioned paper and has not been thoroughly tested.

Data-type files:
- sequence.py: Contains classes NeuriteSequence, ConsensusSequence, and MultipleSequenceAlignment, and various methods for dealing with score matrices and tree-type character mapping.
- matrix.py: Contains several classes for holding matrices and transposing indices. Also contains Iterative Proportional Fitting (IPF) logic designed to determine most significant domains using multiple separate metrics.

Heavy-lifting processing:
- pairwise.py: Performs pairwise global alignment (local alignment code in progress), with wrappers for running multiple pairwise alignments given a set of queries and targets.
- msa.py: Performs multiple sequence alignment, including iteration using PSSM.
- domain.py: Analyze a consensus sequence and extract out domains which satisfy user-provided arguments.

Other processing:
- buffer.py: A high-level module which performs functionality having to do with writing results produced from analysis.
- parameter.py: Handles and validates program arguments (argparse) and produces InputStateWrapper object.
- score_converter.py: Converts raw alignment scores to per-character scores.
- validate_tree.py: Tests input sequence for whether it is a valid tree sequence.
