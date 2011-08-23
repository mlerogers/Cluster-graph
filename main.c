#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashtbl.h"
#include "cluster.h"

HASHTBL *bp;
HASHTBL *marginals;
HASHTBL *idhash;
HASHTBL *freq;
HASHTBL *cluster;
HASHTBL *binary;
HASHTBL *bracket;
char *OUTPUT;
char *INPUT;
int FILTER;
int VERBOSE;
int THRESH_FREQ;
int THRESH_COMMON;
int NUMFREQ;
int PRUNE;
int NUMPROF;
int PROF_FREQ;
int NUMSTRUCTS;
int LENGTH;
int STATS;

//input first the fasta file, then the sample_1000.out file run on the fasta
int main(int argc, char *argv[]) {
  int i,total;
  char **mostfreq;
  FILE *fp;

  OUTPUT = DEF_OUTPUT;
  INPUT = NULL;
  FILTER = 0;
  PRUNE = 0;
  THRESH_FREQ = DEF_THRESH_FREQ;
  THRESH_COMMON = DEF_THRESH_COMMON;
  VERBOSE = 0;
  NUMFREQ = 0;
  NUMPROF = DEF_NUMPROF;
  PROF_FREQ = 0;
  NUMSTRUCTS = 0;
  LENGTH = 0;
  STATS = 0;

  if (argc < 3) {
    fprintf(stderr,"Not enough arguments\n");
    exit(EXIT_FAILURE);
  }
  for (i = 3; i < argc; i++) {
    //printf("argv[%d] is %s\n",i,argv[i]);
    if (!strcmp(argv[i],"-f")) {
      FILTER = 1;
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMFREQ)) {
	//	sscanf(argv[i+1],"%s",val);
	//printf("val is %s and argv %s\n",val,argv[i+1]);
	NUMFREQ = atoi(argv[i+1]);
	i++;
      }
      else
	NUMFREQ = DEF_NUMFREQ;
    }
    else if (!strcmp(argv[i],"-t")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&THRESH_FREQ)) {
	THRESH_FREQ = atoi(argv[i+1]);
	if (THRESH_FREQ < 1 || THRESH_FREQ > 100) {
	  fprintf(stderr,"Error: invalid input for frequency threshold\n");
	  THRESH_FREQ = DEF_THRESH_FREQ;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-c")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&THRESH_COMMON)) {
	THRESH_COMMON = atoi(argv[i+1]);
	if (THRESH_COMMON < 1 || THRESH_COMMON > 100) {
	  fprintf(stderr,"Error: invalid input for common threshold\n");
	  THRESH_COMMON = DEF_THRESH_COMMON;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-p")) {
      //PRUNE = 1;
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMPROF)) {
	NUMPROF = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-q")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&PROF_FREQ)) {
	PROF_FREQ = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-l")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&LENGTH)) {
	LENGTH = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-s")) {
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMSTRUCTS)) {
	NUMSTRUCTS = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-o")) {
      if (i + 1 <= argc - 1) {
	OUTPUT = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-i")) {
      if (i + 1 <= argc - 1) {
	INPUT = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-v"))
      VERBOSE = 1;
    else if (!strcmp(argv[i],"-s"))
      STATS = 1;
  }
  
  if (!(bp = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for bp failed");
    exit(EXIT_FAILURE);
  }
  
  if (!(marginals = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for marginals failed");
    exit(EXIT_FAILURE);
  }

  if (!(idhash = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for idhash failed");
    exit(EXIT_FAILURE);
  }

  if (!(binary = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for binary failed");
    exit(EXIT_FAILURE);
  }

  total = process_structs(argv[1],argv[2]);
  if (VERBOSE) {
    printf("Threshold to find frequent helices: %d\%\n",THRESH_FREQ);
    printf("Maximum number of frequent helices: ");
    if (NUMFREQ == 0)
      puts("no limit");
    else
      printf("%d\n",NUMFREQ);
    printf("Maximum number of profiles: %d\n", NUMPROF);
    printf("Number of structures processed: %d\n",NUMSTRUCTS);
 }
  printf("Total number of equivalence helix classes: %d\n",total-1);
  if (VERBOSE)
    print_all_helices(total);

  //  make_graph(marginals,max,id,total,argv[1],fp);

  mostfreq = find_freq(total);
  printf("Total number of frequent helices: %d\n",hashtbl_numkeys(freq));
  cluster = make_cluster(argv[2],mostfreq);
  if (VERBOSE) {
    fp = fopen("cluster.out","w");
    fprintf(fp,"Processing %s\n",argv[2]);
    print_cluster(fp);
    fclose(fp);
  }
  printf("Total number of clusters: %d\n",hashtbl_numkeys(cluster)-1);
  
  fp = fopen(OUTPUT,"w");
  insert_graph(fp);  
  fputs("}",fp);
  fclose(fp);
  
  hashtbl_destroy(bp);
  hashtbl_destroy(marginals);
  hashtbl_destroy(idhash);
  hashtbl_destroy(binary);
  hashtbl_destroy(freq);
  hashtbl_destroy(cluster);
  hashtbl_destroy(bracket);

  return 0;
}
