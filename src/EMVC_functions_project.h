#ifndef EMVC_FUNCTIONS_PROJECT
#define EMVC_FUNCTIONS_PROJECT

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define P 1000        //s x-axis length. Number of nucleotides that with be saved in one sweep
#define M 7         //number of learners
#define L 4         //number of different labels
#define K 10        //number of different classes
#define DISPLAY_PERIOD 100000	//how often the mpileup status message is updated
#define LOG(x) log2(x)
#define EXP(x) pow(2, x)

#define BED_REGIONS 10000000

unsigned int count_lines(FILE* file_ptr);
unsigned long * read_truth_bed_file (FILE * truth_bed_file, unsigned long bed_regions);

int decide_type(int pos, char ref, int *len, int *fwd_ref, int *fwd_alt, int *rev_ref, int *rev_alt, char *seq, char *qs, unsigned long * truth_bed);
void clean_seq(char *seq, char *qs, int *len, int *fwd_ref, int *fwd_alt, int *rev_ref, int *rev_alt);

void add_to_nobv(int pos, int len, int fwd_ref, int fwd_alt, int rev_ref, int rev_alt, char ref, char *seq, char *qs, FILE *nobv_file, int *s, int *nobv_indxs, int *nobv_cnt, char *nobv_file_name);
void add_to_obv_var(int pos, int len, int fwd_ref, int fwd_alt, int rev_ref, int rev_alt, char ref, char *seq, FILE *obv_var_indx_file, int *obv_var_indxs, int *obv_var_cnt, char *obv_var_file_name);
void add_to_obv_invar(int pos, char ref, FILE *obv_invar_indx_file, int *obv_invar_indxs, int *obv_invar_cnt, char *obv_invar_file_name);
void decide_m_l(int *m, int *l, int pos, char nuc, int qs);

int* update_matrix(int *matrix, int cnt, int size);
int* trim_matrix(int *matrix, int cnt, int size);

void initialize_pi(int T, double log_pi[K][T+1]);
void initialize_gamma(int T, double log_gamma[L][K][M][T+1]);

void run_EMVC(int N, int T, double *alpha, double log_pi[K][T+1], double log_gamma[L][K][M][T+1], int *s, int *nobv_indxs);

void compare_to_ref(int N, int *nobv_var_indxs, int *nobv_invar_indxs, int *nobv_var_cnt, int *nobv_invar_cnt, int *nobv_indxs);
void compare_to_ground_truth(int *tp, int *fp, int *fn, int obv_var_cnt, int nobv_var_cnt, int *obv_var_indxs, int *nobv_var_indxs, FILE *gt_file);
int read_ground_truth(int *position, FILE *gt_file);
int class_to_num(char a1, char a2);
int nuc_to_num(char nuc);
double compute_precision(int tp, int fp);
double compute_sensitivity(int tp, int fn);

void generate_vcf(char *chr, int obv_var_cnt, int nobv_var_cnt, int *obv_var_indxs, int *nobv_var_indxs, double *alpha, char *file_name);
void calculate_paramVCF(FILE *f, int *var_indxs, int i, int *pos, char *ref, int *len, int *fwd_ref, int *fwd_alt, int *rev_ref, int *rev_alt, int *qual, char *alt, char *gt, int ms);

void print_init();
void print_s(int *s, int cnt);
void print_s_file(int *s, int cnt, char *s_file_name);
void print_indxs(int *indxs, int cnt, int col);
void print_gamma(int T, double log_gamma[L][K][M][T+1], int t);
void print_pi(int T, double log_pi[K][T+1], int t);
void print_info(int tp, int fp, int fn, int obv_var_cnt, int obv_invar_cnt, int nobv_cnt, int nobv_var_cnt, int nobv_invar_cnt, int disc_cnt, double prec, double sens, char *gt_file_name, char *results_file_name, FILE *results_file);

int offset_3d(int x, int y, int z);
int offset_2d(int x, int y, int col);

#endif
