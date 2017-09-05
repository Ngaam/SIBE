# Sibe

[Sibe](http://godzilla.uchicago.edu/pages/ngaam/sibe/) is a biological engine for both protein dynamics & analyses and statistical analysis in genomics. The code of Sibe is at a prototype. Hope that it can be useful in its current version, but it is likely that it contains bugs or other errors that impact its ability to give meaningful results. Sibe is released under the [BSD 2-Clause license](https://github.com/BVLC/caffe/blob/master/LICENSE).

## Functions

```
sibe: usage command line
 sibe <command> <args>

 Commands:
 basic information:
   about                -about sibe
   device_query         -show GPU diagnostic information 
  
 work on the metaomics (see also: sibe help metaomics)
   metaomics            -systematicallyanalyze sequencing data (reserved)
   omics_sff2fastq      -convert an SFF file from the 454 genome sequencer
                         to FASTQ including the sequences and quality scores
   omics_stats          -basic statistical analysis on omics data
   sequence_stats       -statistical analysis on a protein mutiple sequence alignment
   sequence_design      -design a protein sequence
   sequence_trim        -trim a multiple protein sequence alignment
   sequence_energy      -calculate energy for a protein sequence or multiple sequences
   sequence_potential   -estimate a multiple protein sequence alignment
   point_mutation       -point mutation for a given protein sequence
   residue_coupling     -coupling relationship betwen pairwise protein residues
   
 work on statistics tools (see also: sibe help statistics)
   stats_pca             -principle component analysis on a given matrix
   stats_ica             -independent component analysis on a given matrix
 
 work on optimization tools (see also: sibe help optimization)
   ai_chpso             -convergent heterogeneous particle swarm optimizer
 
 work on protein folding (see also: sibe help folding)
   fold_protein         -predict tertiary structure
   protein_folding      -predict folding pathways & tertiary structure
   residue_contact      -contacts between pairwise residues
   pdb_parser           -parser PDB file
 
 work on deep learning (see also: sibe help learning)
   learning             -learn a deep neural network from a given data-set 
```

## Installation on Unix-like system
### Dependencies 
CMake (2.8+) 
```
wget https://cmake.org/files/v3.8/cmake-3.5.1.tar.gz 
./configure 
gmake && make install 
cmake --version 
```
Boost (1.58.0+)
```
sudo apt-get install libboost-all-de
```
or 
```
https://dl.bintray.com/boostorg/release/1.58.0/source/boost_1_58_0.tar.gz
./configure 
./b2 --libdir=/usr/local/lib --includedir=/usr/local/include 
./b2 install 
```
Eigen3
```
sudo apt install libeigen3-dev
```
GFlags  
```
wget https://github.com/schuhschuh/gflags/archive/master.zip
unzip master.zip
cd gflags-master
mkdir build && cd build
export CXXFLAGS="-fPIC" && cmake .. && make VERBOSE=1
make && make install
```
Glog  
```
wget https://github.com/google/glog/archive/v0.3.4.tar.gz
tar zvxf glog-0.3.*.tar.gz
cd glog-0.3.*/
./configure
make && make install
```
### Install Sibe
```
git clone https://github.com/Ngaam/Sibe.git
mkdir build
cd build
cmake ..
make -j && make install
```


## Examples
### Sibe on protein sequence
**Trim a multiple sequence alignment**
```
sibe sequence_trim -msa=example/1hrd.msa
```
**Create a potential from a given multiple sequence alignment**
```
sibe sequence_potential -msa=example/1hrd_trimed.aln
```
**Statistical analysis on a given multiple sequence alignment**
```
sibe sequence_statistics -msa=example/1hrd.msa
```
**Design a protein sequence**
```
sibe sequence_design -seq=input_seq.fasta -mat=potential_matrix.mat -dseq=output_seq.fasta
```
**Calculate energy for a protein sequence or multiple ones**
```
sibe sequence_energy -mat=potential_matrix.mat -msa=input_seq.fasta
```
**Analyze coupling relationship betwen pairwise protein residues**
```
sibe residue_coupling -msa=example/1hrd.msa
```
### Sibe on protein folding simulation
**Start from protein amino acid sequence**
```
sibe protein_folding -fasta=example/test.fasta -param=cfg/psibe.cfg
```
**Or start from protein PDB**
```
sibe protein_folding -pdb=example/test.pdb -param=cfg/psibe.cfg
```
**Protein residue contacts prediction**
```
sibe residue_contact -msa=example/test.aln
```
### Sibe on omcis data
**Convert SFF to FASTQ**
```
sibe omics_sff2fastq -sff=test.sff
```
**Basic statistial analysis**
```
sibe omics_statistics -fastx=test_with_quality scores.fastq
```
### Statistical tools and optimization methods
**Principle component analysis on a given matrix**
```
stats_pca -mat=test_matrix.csv
```
**Independent component analysis on a given matrix**
```
stats_ica -mat=test_matrix.csv
```
**Convergent heterogeneous particle swarm optimizer**
```
sibe ai_chpso -ntrial=50 -max_iter=2000 -ngrp=8
```
### Sibe on fMRI data
**Learning**
```
sibe learning -train=*.mat -test=*.mat
```


## Citation
Sibe: a biological engine from protein sequence statistics to its folding and design.

