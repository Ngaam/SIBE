# Sibe

[Sibe](http://godzilla.uchicago.edu/pages/ngaam/sibe/) is a biological engine for both protein dynamics & analyses and statistical analysis in genomics. Sibe is released under the [BSD 2-Clause license](https://github.com/BVLC/caffe/blob/master/LICENSE).

## Functions

```
sibe: usage command line
 sibe <command> <args>

 Commands:
  about                -about sibe
  device_query         -show GPU diagnostic information
  metaomics_statistics -statistical analysis on metaomisc (TODO)
  sequence_statistics  -statistical analysis on a protein multiple sequence alignment
  sequence_design      -design a protein sequence
  sequence_trim        -trim a multiple protein sequence alignment
  sequence_energy      -calculate energy for a protein sequence or multiple sequences
  sequence_potential   -estimate a multiple protein sequence alignment
  fold_protein         -predict tertiary structure
  protein_folding      -predict folding pathways & tertiary structure
  residue_contact      -contacts between pairwise residues
  residue_coupling     -coupling relationship between pairwise protein residues
  point_mutation       -point mutation for a given protein sequence
  pdb_parser           -read PDB file
  metaomics            -analyze sequencing data
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
./sibe sequence_trim -msa=example/1hrd.msa
```
**Create a potential from a given multiple sequence alignment**
```
./sibe sequence_potential -msa=example/1hrd_trimed.aln
```
**Statistical analysis on a given multiple sequence alignment**
```
./sibe sequence_statistics -msa=example/1hrd.msa
```
**Design a protein sequence**
```
./sibe sequence_design -seq=input_seq.fasta -mat=potential_matrix.mat -dseq=output_seq.fasta
```
**Calculate energy for a protein sequence or multiple ones**
```
./sibe sequence_energy -mat=potential_matrix.mat -msa=input_seq.fasta
```
**Analyze coupling relationship betwen pairwise protein residues**
```
./sibe residue_coupling -msa=example/1hrd.msa
```
### Sibe on protein folding simulation
**Start from protein amino acid sequence**
```
./sibe protein_folding -fasta=example/test.fasta -param=cfg/psibe.cfg
```
**Or start from protein PDB**
```
./sibe protein_folding -pdb=example/test.pdb -param=cfg/psibe.cfg
```
**Protein residue contacts prediction**
```
./sibe residue_contact -msa=example/test.aln
```
### Sibe on fMRI data
**Learning**
```
./sibe learning -train=*.mat -test=*.mat
```


## Citation
Sibe: a biological engine from protein sequence statistics to its folding and design.

