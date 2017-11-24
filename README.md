# Sibe

[Sibe](http://godzilla.uchicago.edu/pages/ngaam/sibe/) is a powerful digital engine, it can be used for biological science, particularly in protein dynamics & analyses and statistical analysis in genomics. The code of Sibe is at a prototype. Hope that it can be useful in its current version, but it is likely that it contains bugs or other errors that impact its ability to give meaningful results.

- [What is Sibe for](#what-is-sibe-for)
- [What are the commands](#what-are-the-commands)
- [Where to find Sibe web-server](#where-to-find-sibe-web-server)
- [How to install Sibe](#how-to-install)
    - [Dependencies](#dependency)
    - [Install Sibe](#install-sibe)
- [How to apply Sibe to your data](#how-to-apply-sibe)
    - [Sibe on protein sequence](#sibe-on-protein-sequence)
    - [Sibe on protein folding](#sibe-on-protein-folding)
    - [Sibe on omcis data](#sibe-on-omics)
    - [Tools in Sibe](#sibe-on-stats-opt)
    - [Sibe on learning](#sibe-on-learning)
- [Acknowledgements](#acknowlegements)
- [Reference](#how-to-cite-sibe)


<a name="what-is-sibe-for"></a>
## What is Sibe for?
Sibe is an analytical and computational framework, and it aims to provide a powerful tool for biological science, such as sequence data analysis, <i>in silico</i> protein folding and design. Though much of the software suite is oriented toward basic research on protein sequence analysis, folding and design, Sibe is also designed for extracting meaningful information hidden behind 'big data' based on machine learning. With the help of statistical analysis methods, Sibe can infer co-evolutionary information encoded in protein amino acids sequences for protein folding and design. Now, Sibe includes seven  easy-interfaced modules, several physical- & chemical-principles and statistical analysis methods, as well as different optimization solvers. The advantages of Sibe are:

- Particular development for bioinformatics
- Data-driven analyses and simulations
- Expressive architecture \& extensible code
- Scalable computation with high performance
- Good for research calculations \& industry deployment


<a name="what-are-the-commands"></a>
## What are the commands?
```
sibe: usage command line
 sibe <command> <args>

 Commands:
 basic information:
   about                -about sibe
   device_query         -show GPU diagnostic information 
  
 work on the metaomics (see also: sibe help metaomics)
  metaomics            -systematically analyze sequencing data (reserved)
  omics_sff2fastq      -convert an SFF file from the 454 genome sequencer
                        to FASTQ including the sequences and quality scores
  omics_stats          -statistical analysis on omics data (TODO)
  quality_control      -quality control for ChIP-seq/RNA-seq data (TODO)
  peak_calling         -peak calling for ChIP-seq/RNA-seq data (TODO)
 
 work on the proteomics (see also: sibe help proteomics)
  sequence_stats       -statistical analysis on a protein mutiple sequence alignment
  sequence_design      -design a protein sequence
  sequence_trim        -trim a multiple protein sequence alignment
  sequence_energy      -calculate energy for a protein sequence or multiple sequences
  sequence_potential   -estimate a multiple protein sequence alignment
  point_mutation       -point mutation for a given protein sequence
  residue_coupling     -coupling relationship betwen pairwise protein residues
  dna                  -calculations for DNA
 
 work on statistics tools (see also: sibe help statistics)
  stats_pca            -principle component analysis on a given matrix
  stats_ica            -independent component analysis on a given matrix
  stats_tsne           -t-distributed stochastic neighbor embedding (t-SNE) method
 
 work on optimization tools (see also: sibe help optimization)
  opt_chpso            -the convergent heterogeneous particle swarm optimizer
  opt_lbfgs            -the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm with limited-memory
  opt_gd               -the Gradient Descent method
 
 work on protein folding (see also: sibe help folding)
  fold_protein         -predict tertiary structure
  protein_folding      -predict folding pathways & tertiary structure
  residue_contact      -contacts between pairwise residues
  pdb_parser           -parser a PDB file and write it out in Sibe format
 
 work on deep learning (see also: sibe help learning)
  img_learning         -learn a deep neural network from a given image data-set
  fmri_learning        -learn a deep neural network from a given fMRI data-set
  phsior               -protein torsion angle predictor based on regressive CNN deep network

```

<a name="where-to-find-sibe-web-server"></a>
## Where to find Sibe web-server?
A web-sever of Sibe is in developing, please find more details at [Sibe web-server](http://godzilla.uchicago.edu/pages/ngaam/sibe2/index.html).

<a name="how-to-install"></a>
## How to install Sibe on Unix-like system?
<a name="dependency"></a>
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
<a name="install-sibe"></a>
### Install Sibe
```
git clone https://github.com/Ngaam/Sibe.git
mkdir build
cd build
cmake ..
make -j && make install
```

<a name="how-to-apply-sibe"></a>
## How to apply Sibe to your data?
<a name="sibe-on-protein-sequence"></a>
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
<a name="sibe-on-protein-folding"></a>
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
<a name="sibe-on-omics"></a>
### Sibe on omcis data
**Convert SFF to FASTQ**
```
sibe omics_sff2fastq -sff=test.sff
```
**Basic statistical analysis**
```
sibe omics_statistics -fastx=test_with_quality scores.fastq
```
<a name="sibe-on-stats-opt"></a>
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
<a name="sibe-on-learning"></a>
### Sibe on fMRI data
**Learning**
```
sibe fmri_learning -train=*.mat -test=*.mat
```

### Authors
[N. J. Cheung](#nyaam.ch@gmail.com)

### Contributors

<a name="acknowlegements"></a>
### Acknowledgements

<a name="how-to-cite-sibe"></a>
### How to cite Sibe?
Sibe: a biological engine from protein sequence statistics to its folding and design.


