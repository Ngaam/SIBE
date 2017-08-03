/**
 * The C++ Package of SIBE (Statistical Informatics for Biological Engine).
 * Copyright (C) 2015- All Rights Reserved - Ngaam J. Cheung.
 * Contact: ngaam.ch 'AT' gmail.com
 *
 * This library was partially written at the James Franck Institute at 
 * the University of Chicago. Gordon Center for Integrated Science, 
 * E226. 929 East 57th Street, Chicago, IL 60637, USA, 
 * and at Department of Brain & Cognitive Sciences at Daegu Gyeongbuk 
 * Institute of Science & Technology. DGIST E4-505, Daegu 42988, Korea 
 * Homepage: http://godzilla.uchicago.edu/pages/ngaam/sibe/index.html
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 *
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

/**
 * This contains main().
 *
 */
#include <gflags/gflags.h>
#include <glog/logging.h>

#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "boost/algorithm/string.hpp"
//#include "boost/lexical_cast.hpp"

#include "sibe/sibe.hpp"

// #include "sibe/test.hpp"



//using sibe::Blob;
using sibe::Sibe;
using sibe::MT19937;
using sibe::Sequence;
// using sibe::Genome
// using sibe::Metaomics;
// using sibe::BAM;
using sibe::DNA;
using sibe::Atom;
using sibe::AminoAcid;
using sibe::Chain;
using sibe::Protein;
// using sibe::Solver;
// using sibe::Learner;
using sibe::Layer;
using sibe::Network;
using sibe::shared_ptr;

using sibe::Simulation;


using sibe::string;
//using sibe::Timer;
using sibe::vector;
//using std::ifstream;
using sibe::ostringstream;
//using std::ios;

// DEFINE_bool: boolean
// DEFINE_int32: 32-bit integer
// DEFINE_int64: 64-bit integer
// DEFINE_uint64: unsigned 64-bit integer
// DEFINE_double: double
// DEFINE_string: C++ string
DEFINE_string(gpu, "0",
    "Optional; run in GPU mode on given device IDs separated by ',', "
    "e.g. '-gpu=1,2'. '-gpu all' is to run it on all available GPUs.");
DEFINE_string(param, "",
    "'-param=*.par', the paramter settings text file.");
DEFINE_string(fasta, "",
    "'-fasta=*.fasta', the FASTA file.");
DEFINE_string(fastq, "",
    "'-fastq=*.fastq', the FASTQ file.");
DEFINE_string(msa, "",
    "'-msa=*.aln', the MSA text file.");
DEFINE_string(pdb, "",
    "'-pdb=*.pdb', the PDB file.");
DEFINE_string(mat, "",
    "'-mat=*.mat', the matrix text file.");
DEFINE_string(seq, "",
    "'-sep=*.fasta', the input FASTA file for the designed sequence.");
DEFINE_string(dseq, "",
    "'-dseq=*.fasta', the FASTA file of the designed sequence.");
DEFINE_string(model, "",
    "'-model=*.model', the model definition protocol buffer text file.");
DEFINE_string(output, "",
    "'-output=directory_path', the path for the output of Sibe (optional).");
DEFINE_int32(iterations, 200,
    "'-iterations=500', the number of iterations to run (optional).");
DEFINE_int32(iter_lsch, 50,
    "'-iter_lsch=50', the maximum number of iterations in linear search (optional).");
DEFINE_double(temperature, 1.0,
    "'-temperature=1.0', temperature for simulation (optional).");
DEFINE_double(lsbf, 0.5,
    "'-lsbf=0.5', line search backtracking factor in (0,1) (optional).");
DEFINE_double(tlssdc, 0.0001,
    "'-tlssdc=0.0001', tolerance for line search sufficient decrease criterion (optional).");
DEFINE_double(tlscc, 0.1,
    "'-tlscc=0.1', tolerance for line search curvature criterion (optional).");
DEFINE_double(threshold, 0.8,
    "'-threshold=0.8', threshold (optional).");

DEFINE_string(train, "",
    "'-train=*', training samples with labels starting from 0.");
DEFINE_string(test, "",
    "'-test=*', test samples with labels starting from 0.");
DEFINE_int32(echo, 20,
    "'-echo=20', number of echoes (optional), its default value is 20.");

// A simple registry for sibe commands.
typedef int (*BrewFunction)();
typedef std::map<sibe::string, BrewFunction> BrewMap;
BrewMap g_brew_map;

#define RegisterBrewFunction(func) \
namespace { \
class __Registerer_##func { \
 public: /* NOLINT */ \
  __Registerer_##func() { \
    g_brew_map[#func] = &func; \
  } \
}; \
__Registerer_##func g_registerer_##func; \
}

static BrewFunction GetBrewFunction(const sibe::string& name);
static void get_gpus(vector<int>* gpus); // Parse GPU ids or use all available devices
int device_query(); // Device Query: show diagnostic information for a GPU device.
RegisterBrewFunction(device_query);


int metaomics_statistics();
RegisterBrewFunction(metaomics_statistics);

int metaomics();
RegisterBrewFunction(metaomics);




int sequence_statistics(); // Handle Multiple Sequence Alignment
RegisterBrewFunction(sequence_statistics);

int sequence_design(); // Handle Multiple Sequence Alignment
RegisterBrewFunction(sequence_design);

int sequence_trim();
RegisterBrewFunction(sequence_trim);

int sequence_energy();
RegisterBrewFunction(sequence_energy);

int point_mutation();
RegisterBrewFunction(point_mutation);

// int dca();
// RegisterBrewFunction(dca);

int sequence_potential();
RegisterBrewFunction(sequence_potential);

int residue_contact();
RegisterBrewFunction(residue_contact)
// int fasta();
// RegisterBrewFunction(fasta);
int pdb_parser(); 
RegisterBrewFunction(pdb_parser);

int fold_protein(); 
RegisterBrewFunction(fold_protein);
int protein_folding(); 
RegisterBrewFunction(protein_folding);

int learning(); 
RegisterBrewFunction(learning);

int about();
RegisterBrewFunction(about);
int usage();
RegisterBrewFunction(usage);



int test(); // TODO
RegisterBrewFunction(test);







/*===== Main =====*/
int main(int argc, char *argv[]) {
//sibe::test();
  FLAGS_alsologtostderr = 1;                              // Print output to stderr (while still logging).
  gflags::SetVersionString(AS_STRING(SIBE_VERSION));      // Set version
  gflags::SetUsageMessage("usage command line\n"          // Usage message.
      " sibe <command> <args>\n"
      " \n"
      " These are common Sibe commands used in various situations:\n"
      " \n"
      " get sibe and system info. (see also: sibe help tutorial)\n"
      "  about                -about sibe\n"
      "  help                 -sibe help\n"
      "  device_query         -show GPU diagnostic information\n"
      "  time                 -execution time\n"
      " \n"
      " work on the metaomics (see also: sibe help metaomics)\n"
      "  metaomics            -analyze sequencing data (TODO: working on)\n"
      "  metaomics_statistics -statistical analysis on metaomics (TODO)\n"
      "  sequence_statistics  -statistical analysis on a protein mutiple sequence alignment\n"
      "  sequence_design      -design a protein sequence\n"
      "  sequence_trim        -trim a multiple protein sequence alignment\n"
      "  sequence_energy      -calculate energy for a protein sequence or multiple sequences\n"
      "  sequence_potential   -estimate a multiple protein sequence alignment\n"
      "  point_mutation       -point mutation for a given protein sequence\n"
      "  residue_coupling     -coupling relationship betwen pairwise protein residues\n"
      "  dna                  -calculations for DNA\n"
      " \n"
      " work on protein folding (see also: sibe help folding)\n"
      "  fold_protein         -predict tertiary structure (TODO)\n"
      "  protein_folding      -predict folding pathways & tertiary structure (TODO: working on)\n"
      "  residue_contact      -contacts between pairwise residues\n"
      "  pdb_parser           -read PDB file\n"
      // "  fasta                -amino acid sequnce (FASTA)\n"
      " \n"
      " work on deep learning (see also: sibe help learning)\n"
      "  learning             -learn a deep neural network from a given data-set (now only fMRI & MNIST)\n"
      " \n"
      " 'sibe help -a' and 'sibe help -g' list available subcommands and some concept guides.\n"
      "  See 'sibe help <command>' or 'sibe help <concept>' to read about a specific subcommand or concept."
      " \n"
      );
  sibe::GlobalInit(&argc, &argv);  // Run tool or show usage.

  if(argc < 2) {
    //GetBrewFunction(sibe::string(about))();
    if (argc == 1) {
      //GetBrewFunction(sibe::string(usage))();
    }
    gflags::ShowUsageWithFlagsRestrict(argv[0], "tool/sibe");
    return 0;
  }

  Timer tim;
  bool time_disp = true;
  //double t_load;
  double t_pass;
  if(string(argv[1]) != "about" && string(argv[1]) != "usage" && string(argv[1]) != "help") {

    //t_load = t.elapsed();
    tim.restart();

    printf("\n");
    printf("-- Sibe is running: %s\n", argv[1]);
    if(!FLAGS_output.size()) {
      FLAGS_output = "example/results/";
      printf("--  Sibe will save all output files in default directory: %s\n", FLAGS_output.c_str());
    }
    else {
      printf("--  Sibe will save all output files in the directory: %s\n", FLAGS_output.c_str());
    }
    string path = FLAGS_output;
    int found = path.find_last_of("/\\");
    int len = path.length();
    // if( path.compare(path.size()-1, 1, "/") ) { path += "/"; }
    if( found != len-1 ) { path += "/"; }
    FLAGS_output = path;
  }

  GetBrewFunction( sibe::string(argv[1]) )();

  if(string(argv[1]) != "about" && string(argv[1]) != "usage" && string(argv[1]) != "help") {
    t_pass = tim.elapsed();
    printf("--  Status Report (Step %i / Maximum iter. %i)\n", 1, 1);
    printf("\033[1;34m");
    //printf("--  Energy: %f\n", -1000.0);
    printf("\033[0m");



    if(time_disp) {
      printf("\033[1;34m");
      printf("--  Time Status Report (seconds)\n");
      printf("--    Loading Data:  %f\n", t_pass);
      printf("--  Total: %f\n", t_pass);
      printf("\033[0m");
    }
    tim.restart();
    printf("-- Sibe is completed: %s\n\n", argv[1]);
  }
  return EXIT_SUCCESS;
}



/*===== GetBrewFunction =====*/
static BrewFunction GetBrewFunction(const sibe::string& name) {
  if (g_brew_map.count(name)) {
    return g_brew_map[name];
  }
  else {
    LOG(ERROR) << "Available commands:";
    for (BrewMap::iterator it = g_brew_map.begin();
         it != g_brew_map.end(); ++it) {
      LOG(ERROR) << "\t" << it->first;
    }
    LOG(FATAL) << "Unknown command: " << name;
    return NULL;  // not reachable, just to suppress old compiler warnings.
  }
}
/*===== Get GPUs =====*/
// Parse GPU ids or use all available devices
static void get_gpus(vector<int>* gpus) {
  if (FLAGS_gpu == "all") {
    int count = 0;
#ifndef CPU_ONLY
    CUDA_CHECK(cudaGetDeviceCount(&count));
#else
    NO_GPU;
#endif
    for (int i = 0; i < count; ++i) {
      gpus->push_back(i);
    }
  }
  else if (FLAGS_gpu.size()) {
    // vector<string> strings;
    // boost::split(strings, FLAGS_gpu, boost::is_any_of(","));
    // for (int i = 0; i < strings.size(); ++i) {
    //   gpus->push_back(boost::lexical_cast<int>(strings[i]));
    // }
  }
  else {
    CHECK_EQ(gpus->size(), 0);
  }
}

/*===== Device Query =====*/
// sibe commands to call by
//     sibe <command> <args>
//
// To add a command, define a function "int command()" and register it with
// RegisterBrewFunction(action);
// Device Query: show diagnostic information for a GPU device.
int device_query() {
  CHECK_GT(FLAGS_gpu.size(), 0) << "Please input GPU id to query!";
  LOG(INFO) << "Querying GPUs " << FLAGS_gpu;
  vector<int> gpus;
  get_gpus(&gpus);
  for (int i = 0; i < gpus.size(); ++i) {
    sibe::Sibe::SetDevice(gpus[i]);
    sibe::Sibe::DeviceQuery();
  }
  return 0;
}


int metaomics_statistics() {

  return 0;
}

int metaomics() {

  // BAM<double> reader;
  // string fil = "/home/ngaam/learning_code/bamquality/example.bam";
  // reader.read_bam(fil);
  // reader.close_bam();
  // string bamHeader = reader.GetHeaderText();
  // LOG(INFO)<<bamHeader;
  

  return 0;
}






/*===== Protein Sequence Statistics =======*/
int sequence_statistics() {
  CHECK_GT(FLAGS_msa.size(), 0) << "Use '-msa=*' -output=* (optional) to add a Multiple Sequence Alignment (MSA)!";
  string postfix;

  Sequence<double> sequence;
  int nrow, ncol;
  vector<vector<int> > msa;
  //vector<double> _weight;

  // Read MSA
  sequence.read_msa(msa, &nrow, &ncol, FLAGS_msa);
  vector<vector<double> > _similarity_matrix;
  int nrow_scale = 5e3;
  if(nrow <= nrow_scale) {
    nrow_scale = nrow;
  } 
  else {
    // nrow_scale = 5e1;
    printf("--  The size of the MSA is too large, Sibe will compute sequence similarities among the first %i sequences!\n", nrow_scale);
  }
  printf("--  Number of sequences in the analyzed MSA: %i\n", nrow_scale);
  // double sim = 0.0;
  // _similarity_matrix.resize(nrow_scale, vector<double>(nrow_scale, 0.0));
  // for(int i = 0; i < nrow_scale; i++) {
  //   for(int j = i + 1; j < nrow_scale; j++) {
  //     sequence.sequence_similarity(&sim, msa[i], msa[j]);
  //     _similarity_matrix[i][j] = sim;
  //     _similarity_matrix[j][i] = sim;
  //   }
  // }
  // // Write similarities into a file.
  string matrix_name;
  // postfix = "similarity_matrix.txt";
  // sibe::define_output_file_name(matrix_name, FLAGS_msa, FLAGS_output, postfix);
  // sibe::write_matrix<double>(matrix_name, _similarity_matrix, nrow_scale, nrow_scale);

  //double _identity = 0.2; //0.8
  //sequence.sequence_weight(msa, _weight, _identity);
  // Positional conservation: Kullback-Leibler relative entropy
  int n_aa = 20;
  vector<vector<double> > relative_entropy;
  relative_entropy.clear();
  relative_entropy.resize(ncol, vector<double>(1, 0.0));

  vector<vector<double> > entropy_matrix;
  entropy_matrix.clear();
  entropy_matrix.resize(n_aa, vector<double>(ncol, 0.0));
  sequence.positional_conservation(relative_entropy, entropy_matrix, msa, nrow, ncol);

  // Write the matrix of relative entropies of each amino acid at each position
  // matrix_name = FLAGS_output + nam + "_entropy_matrix.txt";
  postfix = "entropy_matrix.txt";
  sibe::define_output_file_name(matrix_name, FLAGS_msa, FLAGS_output, postfix);
  sibe::write_matrix<double>(matrix_name, entropy_matrix, n_aa, ncol);
  postfix = "positional_conservation.txt";
  sibe::define_output_file_name(matrix_name, FLAGS_msa, FLAGS_output, postfix);
  sibe::write_matrix<double>(matrix_name, relative_entropy, ncol, 1);
  
  sibe::figure(FLAGS_msa, ncol, FLAGS_output);
  
  return 0;
}



// /*===== Direct Coupling Analysis ============*/
// int dca() {
//   CHECK_GT(FLAGS_msa.size(), 0) <<
//     "Please check the multiple sequence alignment (MSA)!";
//   if(!FLAGS_msa.size()) {
//     LOG(ERROR) << "Usage: ./sibe dca -msa=* ";
//   }
//   return 0;
// }
/*===== Sequence Potential ============*/
int sequence_potential() {
  CHECK_GT(FLAGS_msa.size(), 0)
    << "Please check the multiple sequence alignment (MSA)!";
  if(!FLAGS_msa.size()) {
    LOG(ERROR) << "Usage: ./sibe potentials -msa=* -output=* (optional) ";
  }

  double reweight_threshold = FLAGS_threshold;//F08; // 0.8
  Sequence<double> sequence;
  int nrow, ncol;
  vector<vector<int> > msa;
  vector<double> weight;
  sequence.read_msa(msa, &nrow, &ncol, FLAGS_msa);
  printf("--  Number of sequences in the analyzed MSA: %i\n", nrow);
  //LOG(INFO) << "Size of MSA is: " << nrow << "x" <<ncol;
  //LOG(INFO) << msa[1][2];
  //LOG(INFO)<<FLAGS_msa;
  FILE* msa_fid = fopen(FLAGS_msa.c_str(), "r");
  if(msa_fid == NULL) {
    LOG(ERROR) << "Cannot open " << FLAGS_msa.c_str();
  }

  int nsingle = ncol * (N_COMMON_AA - 1);
  //int nvar = nsingle + ncol * ncol * N_COMMON_AA * N_COMMON_AA;
  // int nsingle = nsingle;// + N_COMMON_AA - (nsingle % N_COMMON_AA);
  int nvar = nsingle + ncol*ncol*N_COMMON_AA*N_COMMON_AA;

  vector<double> var(nvar, 0.0);
  sequence.potential_bias_init(var, nrow, ncol, msa);
  vector<double> weights(nrow, 0.0);
  sequence.potential_initialization(msa, weights, reweight_threshold, nrow, ncol);

  double val;
  int para_int[] = {FLAGS_iterations, FLAGS_iter_lsch, nrow, ncol};//, nvar, nsingle};
  double param[] = {FLAGS_lsbf, FLAGS_tlssdc, FLAGS_tlscc};
  sequence.potential_conjugate_gradient(msa, para_int, var, &val, param, weights);
  string mtx_name;// = FLAGS_output + nam + "_potential.mat";

  string postfix = "potential.mat";
  sibe::define_output_file_name(mtx_name, FLAGS_msa, FLAGS_output, postfix);

  FILE* rawfile = fopen(mtx_name.c_str(), "w");
  sibe::matrix_write(rawfile, var, ncol);
  fclose(rawfile);
  return 0;
}
/*===== Protein Sequence Design ============*/
int sequence_design() {
  if(!FLAGS_fasta.size() || !FLAGS_mat.size() || !FLAGS_dseq.size()) {
    LOG(ERROR) << "Usage: sibe design -fasta=* -mat=* -dseq=*";
  }
  CHECK_GT(FLAGS_mat.size(), 0) << "Use '-mat=*.mat' to add a matrix file, then you can calculate sequence 'energy' and do protein sequence design!";
  Sequence<double> sequence;
//  vector<string> _msa_aln;
//  vector<string> _raw_mat;
  

  // Design a new sequence
  //LOG(INFO)<<FLAGS_msa.size();
  // LOG(INFO) << "Initial sequence: " << FLAGS_fasta;
  // LOG(INFO) << "Potential: " << FLAGS_mat;
  // LOG(INFO) << "Max_iteration: " << FLAGS_iterations;
  // LOG(INFO) << "Sequence trajectory: " << FLAGS_dseq;
  int pos = 0;
  pos = FLAGS_fasta.find_last_of("/");
  printf("--  WT sequence: %s, ", (FLAGS_fasta.substr(pos+1)).c_str());
  pos = FLAGS_dseq.find_last_of("/");
  printf("sequence trajactory: %s, ", (FLAGS_dseq.substr(pos+1)).c_str());
  printf("iterations: %i, ", FLAGS_iterations);
  printf("temperature: %f\n", FLAGS_temperature);
  sequence.sequence_design(FLAGS_fasta, FLAGS_dseq, FLAGS_mat, FLAGS_iterations, FLAGS_temperature);
  return 0;
}
/*===== Protein Sequence Trim ============*/
int sequence_trim() {
  if( !FLAGS_msa.size() ) {
    LOG(ERROR) << "Usage: sibe sequence_trim -msa=* -output=* (optional) ";
  }
  CHECK_GT(FLAGS_msa.size(), 0) << "Use '-msa=*.aln' to load a MSA file, then you can trim the sequence alignemt!";
  Sequence<double> sequence;
  if(FLAGS_msa.size()) {
    sequence.sequence_trim(FLAGS_msa, FLAGS_output);
  }
  return 0;
}
/*===== Protein Sequence Energy ============*/
int sequence_energy() {
  // TODO: path is big bug!!
  if( !FLAGS_msa.size() || !FLAGS_mat.size() ) {
    LOG(ERROR) << "Usage: sibe sequence_energy -msa=* -mat=* -output=* (optional) ";
  }
  CHECK_GT(FLAGS_msa.size(), 0) << "Use '-msa=*.aln' to load a MSA file, then you can calculate sequence 'energy'!";
  CHECK_GT(FLAGS_mat.size(), 0) << "Use '-mat=*.mat' to load a Mat file (potentials), then you can calculate sequence 'energy'!";
  Sequence<double> sequence;
//  vector<string> _msa_aln;
//  vector<string> _raw_mat;
  //LOG(INFO) << "Here we are at Protein Sequence Energy: " << FLAGS_msa.size();
  // Calculate energy for multiple sequences
  if(FLAGS_msa.size()) {
    //LOG(INFO) << "Here we are at Protein Sequence Energy";
    //sequence.sequence_trim(FLAGS_msa);
    sequence.calculate_sequence_energy(FLAGS_msa, FLAGS_mat, FLAGS_output);
  }
  return 0;
}
/*===== Protein Point Mutation ============*/
int point_mutation() {
  // ddG = dG_mut-dG_wt
  if( !FLAGS_fasta.size() || !FLAGS_mat.size() ) {
    LOG(ERROR) << "Usage: sibe point_mutation -fasta=* -mat=* -output=* (optional) ";
  }
  CHECK_GT(FLAGS_fasta.size(), 0) << "Use '-fasta=*.fasta' to load a FASTA file, then you can try point mutation!";
  CHECK_GT(FLAGS_mat.size(), 0) << "Use '-mat=*.mat' to load a Mat file (potentials), then you can try point mutation!";
  Sequence<double> sequence;
  if(FLAGS_fasta.size()) {
    // string fasta_wt;
    // sequence.read_fasta(fasta_wt, FLAGS_fasta);
    // int aa_num = fasta_wt.length();
    sequence.point_mutation(FLAGS_fasta, FLAGS_mat, FLAGS_output);
  }
  return 0;
}
/*===== Residue-Contact ============*/
int residue_contact() {
  CHECK_GT(FLAGS_msa.size(), 0)
    << "Please check the multiple sequence alignment (MSA)!";
  if(!FLAGS_msa.size()) {
    LOG(ERROR) << "Usage: ./sibe residue_contact -msa=* -output=* (optional) ";
  }

  double reweight_threshold = FLAGS_threshold; 
  Sequence<double> sequence;
  int nrow, ncol;
  vector<vector<int> > msa;
  vector<double> weight;
  sequence.read_msa(msa, &nrow, &ncol, FLAGS_msa);
  printf("--  Number of sequences in the given MSA: %i\n", nrow);
  FILE* msa_fid = fopen(FLAGS_msa.c_str(), "r");
  if(msa_fid == NULL) {
    LOG(ERROR) << "Cannot open " << FLAGS_msa.c_str();
  }

  int nsingle = ncol * (N_COMMON_AA - 1);
  int nvar = nsingle + ncol * ncol * N_COMMON_AA * N_COMMON_AA;

  vector<double> var(nvar, 0.0);
  sequence.potential_bias_init(var, nrow, ncol, msa);
  vector<double> weights(nrow, 0.0);
  sequence.potential_initialization(msa, weights, reweight_threshold, nrow, ncol);

  double val;
  int para_int[] = {FLAGS_iterations, FLAGS_iter_lsch, nrow, ncol};//, nvar, nsingle};
  double param[] = {FLAGS_lsbf, FLAGS_tlssdc, FLAGS_tlscc};
  sequence.potential_conjugate_gradient(msa, para_int, var, &val, param, weights);

  vector<vector<double> > mat;
  mat.resize(ncol, vector<double>(ncol, 0.0));
  sequence.residue_contact(var, mat, ncol);

  string mtx_name;
  string postfix = "contact.mat";
  sibe::define_output_file_name(mtx_name, FLAGS_msa, FLAGS_output, postfix);

  FILE* rawfile = fopen(mtx_name.c_str(), "w");
  sibe::write_matrix(rawfile, mat, ncol, ncol);
  fclose(rawfile);

  return 0;
}

/*===== PDB =======*/
int pdb_parser() {
  CHECK_GT(FLAGS_pdb.size(), 0)
    << "Please check the PDB file!";
  if(!FLAGS_pdb.size()) {
    LOG(ERROR) << "Usage: ./sibe pdb_parser -pdb=* -output=* (optional) ";
  }

  // sibe::ProteinParameter protein_param;
  // Protein<double> protein(protein_param);
  // Protein<double> protein();
  // Protein<double> protein();
  // protein.pdb_parser("/home/ngaam/sibe/example/test0.pdb");
  // protein.pdb_parser(FLAGS_pdb);
  // // protein.
  // // LOG(INFO)<<"at num: " << protein.protein_get_num_of_atoms() 
  // // << ", aa num: " << protein.protein_get_num_of_aminoacids() 
  // // << ", ch num: " << protein.protein_get_num_of_chains();
  // printf("--  PDB: %s\n", FLAGS_pdb.c_str());
  // for(int k = 0; k < protein.chains_.size(); k++) {
  //   printf("--  PDB: chain %i includes %i amino acids\n", (int)protein.chains_.size(), (int)protein.chains_[k].amino_acids_.size()); 
  // }
  // Atom <double> atom;
  // AminoAcid<double> amino_acid;
  // Chain<double> chain;

  // vector<vector<double> > avg_contact_matrix;
  // protein.residue_contact(avg_contact_matrix, "A", protein.chains_[0], 7.5, 1.0, 3, "CA", FLAGS_output);

  return 0;
}

/*===== Fold =====*/
int fold_protein() {
  // CHECK_GT(FLAGS_param.size(), 0) << "Need a -param=*.par of parameter settings to fold!";

  
  return 0;
}
/*===== Folding =====*/
int protein_folding() {
  CHECK_GT(FLAGS_fasta.size(), 0)
    << "Please check the sequence file!";
  CHECK_GT(FLAGS_param.size(), 0)
    << "Please check the paramter file!";
  CHECK_GT(FLAGS_param.size(), 0) 
    << "Use '-fasta=*.fasta' '-param=*.par' -output=* (optional) to run protein folding!";
  // Set device id and mode
//   vector<int> gpus;
//   get_gpus(&gpus);
//   if (gpus.size() != 0) {
//     LOG(INFO) << "Use GPU with device ID " << gpus[0];
// #ifndef CPU_ONLY
//     cudaDeviceProp device_prop;
//     cudaGetDeviceProperties(&device_prop, gpus[0]);
//     LOG(INFO) << "GPU device name: " << device_prop.name;
// #endif
//     Sibe::SetDevice(gpus[0]);
//     Sibe::set_mode(Sibe::GPU);
//   }
//   else {
//     LOG(INFO) << "Use CPU.";
//     Sibe::set_mode(Sibe::CPU);
//   }



  if(FLAGS_fasta.size() > 0) { // Start from sequence
    // Sequence<double> sequence;
    // string fasta_wt;
    // sequence.read_fasta(fasta_wt, FLAGS_fasta);
    // int num_aa = fasta_wt.length();
    // int num_atom = 332;
    // int num_chain = 1;
  }
  if(FLAGS_pdb.size() > 0) {

  }

  Simulation<double> simulation;
  // LOG(INFO)<<FLAGS_param.c_str()<<"\t"<<FLAGS_output;
  // simulation.module;

  simulation.select_driver_by_index(0); // There is 1
  simulation.run_simulation(FLAGS_fasta, FLAGS_param.c_str(), FLAGS_output); // in simulation.cpp

  // LOG(INFO)<<"DBL_EPSILON = "<<DBL_EPSILON;

// sibe::TestStatic();



  return 0;
}
/*===== Learning =====*/
int learning() {
// Labels MUST start from 0!!!!
  CHECK_GT(FLAGS_train.size(), 0)
    << "Please check the trainning samples file with full path!";
  CHECK_GT(FLAGS_test.size(), 0)
    << "Please check the test samples file with full path!";
  Network<double> network;
  int width, high;
  // Load fMRI data
  network.load_data_from_mat(network.train_label, network.train_data, FLAGS_train);
  network.load_data_from_mat(network.test_label, network.test_data, FLAGS_test);
  width = 486;
  high  = 522;
  // MNIST image: 6000(samples)x32x32x1
  // network.load_given_data();
  // width = 32;
  // high  = 32;

  int output_num = *max_element(network.train_label.begin(), network.train_label.end())+1;
  network.set_epoch_num(FLAGS_echo);
  //set callback funtion being called everytime when an epoch finished
  network.set_callback_func(sibe::learning_epoch<double>);
  // // load trained network
  // // network.load_network("weights.txt");

  sibe::Layer<double> *pl = NULL;
  //add input layer with 32 X 32 size and 1 channel
  // pl = network.insert_input_layer(32, 32, 1); // User's definition for width and high
  pl = network.insert_input_layer(width, high, 1); // User's definition for width and high
  //pl = network.insert_partial_weight_sharing_conv_layer(5, 5, 5, 5, 6);
  //add a full wights sharing convolutional layer with 5 X 5 filter size and 6 feature maps
  pl = network.insert_full_weight_sharing_conv_layer(pl, 5, 5, 6);
  // //add a maxpooling layer with 2 X 2 filter size
  pl = network.insert_max_pooling_layer(pl, 2, 2);
  // //add a full connected layer with 120 units output, using tanh activation function.
  pl = network.insert_full_connect_layer(pl, 120, TANH);
  // //add a softmax layer(output layer) with #odim units output
  network.insert_softmax_layer(pl, output_num);

  printf("--  Sibe's NiNet is being trained on the given data-set ...\n");
  // string weight_mat_name;
  if(network.train(network.train_data, network.train_label, FLAGS_output)) {
    printf("--  Sibe completed to train the NiNet network!\n");
  }
  // //save network
  // string postfix = ".wgt";
  // LOG(INFO)<<FLAGS_output;
  // sibe::define_output_file_name(weight_mat_name, FLAGS_msa, FLAGS_output, postfix);
  // network.save_network("model.txt");



  return 0;
}



/*===== About =====*/
int about() {
  //system("clear");
  printf( "===========================================================\n");
  printf( "|  SIBE: Statistical Information for Biological Engine    |\n");
  printf( "===========================================================\n");
  printf( "|    Author        : Ngaam J. Cheung                      |\n");
  printf( "|    Email         : ngaam.ch@gmail.com                   |\n");
  printf( "|    Release       : Version 1.2.1                        |\n");
  printf( "|    Release Date  : Jun. 29, 2016.                       |\n");
  printf( "|    Last Modified : Aug. 15, 2015.                       |\n");
  printf( "|    http://godzilla.uchicago.edu/pages/ngaam/index.html  |\n");
  // printf( " PSIBE: Protein Statistical Information for Folding        \n");
  // printf( " PSIBE Incentive Product - SIBE Executable Build           \n");
  printf( "|     Copyright (C) 2015-   by Ngaam J. Cheung.           |\n");
  printf( "-----------------------------------------------------------\n");
//  printf( " James Franck Institute, The University of Chicago         \n");
//  printf( " 929 East 57 Street, Chicago, IL 60637.                    \n");
  printf( "\n");

  printf( " Permission to use, copy, modify, and distribute this      \n");
  printf( " software and its documentation with or without            \n");
  printf( " smodifications and for any purpose and without fee is     \n");
  printf( " hereby granted, provided that any copyright notices       \n");
  printf( " appear in all copies and that both those copyright        \n");
  printf( " notices and this permission notice appear in supporting   \n");
  printf( " documentation, and that the names of the contributors     \n");
  printf( " or copyright holders not be used in advertising or        \n");
  printf( " publicity pertaining to distribution of the software      \n");
  printf( " without specific prior permission.                        \n");
  printf( "-----------------------------------------------------------\n");
  printf("\n");
  return 0;
}
int usage() {
  printf(" _________________________________________ \n");
  printf("/  ______||___   ___|| ______ \\| _______|\n");
  printf("| |           | |    | |    | || |        \n");
  printf("| \\______     | |    | |____| /| |_____ \n");
  printf("\\______  \\    | |    | ______ \\|  _____|\n");
  printf("       | |    | |    | |    | || |       \n");
  printf(" ______/ | ___| |___ | |____/ || |______ \n");
  printf("|_______/_|_________||________/|________|\n");
  return 0;
}



int test() {
  Sequence<double> sequence;
  int nrow, ncol;
  vector<vector<int> > msa;
  vector<double> weight;

  // double identity = 0.0;
  // Read MSA
  sequence.readMSA(msa, &nrow, &ncol, FLAGS_msa);
  // for(int i = 0; i < 10; i++)LOG(INFO)<<msa[0][0];
  sequence.sequence_reweight(msa, weight);
  sequence.count_marginals_in_msa(msa, weight);

  return 0;
}


