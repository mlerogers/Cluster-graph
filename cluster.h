#include <stdio.h>

//size of all hash tables made
#define HASHSIZE 31
//size of char array holding helix id
#define ARRAYSIZE 20
//default frequency threshold in percentage
#define DEF_THRESH_FREQ 10.0
//default frequency for a helix to be a given
#define DEF_THRESH_COMMON 99.0
//default total number of freq helices
#define DEF_NUMFREQ 7
//default profile frequency threshold in percentage
#define DEF_PROF_FREQ 2.0
//default percentage of structures represented by final profiles
#define DEF_THRESH_STRUCT 85.0
//initial size of memory allocated = ARRAYSIZE * SIZE
#define INIT_SIZE 2
//default output file name
#define DEF_OUTPUT "cluster.dot"

extern HASHTBL *bp;
extern HASHTBL *marginals;
extern HASHTBL *idhash;
extern HASHTBL *freq;
extern HASHTBL *cluster;
extern HASHTBL *binary;
extern HASHTBL *bracket;
extern char *OUTPUT;
extern char *INPUT;
//if filter is 1, keep number of freq helices under 10
extern int FILTER; 
//if 1, prints additional info
extern int VERBOSE;
//if freq > thresh_freq, is significant (in percentage)
extern double THRESH_FREQ;
//if freq > thresh_common, assume is in everything
extern double THRESH_COMMON;
//limit to number of frequent helices
extern int NUMFREQ;
//limit to number of profiles made
extern int NUMPROF;
//profiles must have freq at least PROF_FREQ
extern double PROF_FREQ;
//number of structures processed; sfold generates 1000
extern int NUMSTRUCTS;
//lower bound for percent of structures represented by final profiles
extern double THRESH_STRUCT;
//sets minimum length a freq helix must be
extern int LENGTH;
//makes triplet stats in process_structs
extern int STATS;

typedef struct linkedlist {
  int data;
  struct linkedlist *next;
} LIST;

//in cluster.c
char* input_seq(char *seqfile);
int process_structs(char *seqfile,char *name);
void longest_possible(int i,int j,int k,int id);
int match(int i,int j);
int print_all_helices(int total);

char** find_freq(int total);
double calc_entropy(int marg);
int charcompare(const void *v1, const void *v2);
void filter(char *key,char **mostfreq, int count);
int freqcompare(const void *v1, const void *v2);
int binsearch(char **mostfreq, char *key);
void freq_insert(char *key,int marg,int length);
//int make_graph(HASHTBL *marg, HASHTBL *max,HASHTBL *idhash,int total, char *name, FILE *fp);
int make_profiles(char *name);
int select_profiles(char **mostfreq,int notcommon);
int profcompare(const void *v1, const void *v2);
char* process_profile(HASHTBL *halfbrac,char *profile,int numhelix,int *size,int *most);
void make_brackets(HASHTBL *brac, int i, int j, int id);
void make_bracket_rep(HASHTBL *brac,char *profile);
char* resize(int *size,int total,char *s);
char* quicksort(char *profile,char *dup);
int compare(const void *v1, const void *v2);
char* strcat_front(char *s, char *ct);
void prune_profiles(char **mostfreq);
char *delete_helix(char *origprof, char *least,char *modprofile, int *m);
int print_cluster(FILE *fp);

/*in graph.c
int insert_graph(FILE *fp);
int insert_vertices(HASHTBL **graph,char *key, int count,FILE *fp);
int insert_vertex(HASHTBL *hash, char *key, int count);
int BUpowerset(HASHTBL **graph, char **array, int count, int length,FILE *fp) ;
int print_graph(HASHTBL **graph, int count,FILE *fp);
int print_edges(HASHTBL *hash,char **array,int length,FILE *fp);
*/

//in condensed.c
int insert_graph(FILE *fp);
int insert_and_binary(char *key,char *profile,int freq);
void make_key();
struct hashnode_s* insert_LCAs(FILE *fp,char **profileID,int k);
char* convert_binary(char *profile,int binary,int *count);
void process_input(FILE *fp);
HASHTBL* process_input_profile(FILE *fp,HASHTBL *brac,char *fullprofile, int fullnum,char *profile,int numhelix,char *diff,int prob);
KEY* find_parents(char *profile);
char* insert_diff(HASHTBL *temp,char *diff);
char* sort_input(char *profile,int length);
void find_centroid_edges(FILE *fp);
int print_vertices(FILE *fp);
int find_edges(FILE *fp,char *profile, int *found, int count);
void find_LCA_edges(FILE *fp,struct hashnode_s *node,int k);
struct hashnode_s* insert_LCA(char *profile,int num,int s,int frq);
void check_insert_edge(FILE *fp,char *profile,char *origprof);
char* find_diff(HASHTBL *hash,char *profile, char *origprof, int *k1, int *k2);
char* edge_label(HASHTBL *hash,char *profile, char *origprof,int k);
int make_binary(char *profile,int *k);
char* print_edge(KEY *node,char **table,int v,int* sum);
//char* print_edge(HASHTBL *binary,KEY *node,char **table,int v,int* sum);
//int print_condensed_graph(HASHTBL *binary,KEY *nodekey,LIST *nodelist);
