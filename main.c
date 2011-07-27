#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashtbl.h"
#include "cluster.h"

HASHTBL *max;
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
int NUMPROF;
int NUMSTRUCTS;

int main(int argc, char *argv[]) {
  int i,total;
  char **mostfreq;
  FILE *fp;

  OUTPUT = DEF_OUTPUT;
  INPUT = NULL;
  FILTER = 1;
  NUMFREQ = DEF_NUMFREQ;
  NUMPROF = DEF_NUMPROF;
  THRESH_FREQ = DEF_THRESH_FREQ;
  THRESH_COMMON = DEF_THRESH_COMMON;
  VERBOSE = 0;
  NUMSTRUCTS = 0;

  if (argc < 2) {
    fprintf(stderr,"Not enough arguments\n");
    exit(EXIT_FAILURE);
  }
  for (i = 2; i < argc; i++) {
    //printf("argv[%d] is %s\n",i,argv[i]);
    if (!strcmp(argv[i],"-f")) {
      FILTER = 1;
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMFREQ)) {
	//	sscanf(argv[i+1],"%s",val);
	//printf("val is %s and argv %s\n",val,argv[i+1]);
	NUMFREQ = atoi(argv[i+1]);
	i++;
      }
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
      if ((i + 1 <= argc - 1) && sscanf(argv[i+1],"%d",&NUMPROF)) {
	NUMPROF = atoi(argv[i+1]);
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
  }
  if (!(max = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }
  if (!(marginals = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }

  if (!(idhash = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }

  if (!(binary = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() of binary failed");
    exit(EXIT_FAILURE);
  }

  
  total = process_structs(argv[1]);
  if (VERBOSE) {
    printf("Threshold to find frequent helices: %d\%\n",THRESH_FREQ);
    printf("Maximum number of frequent helices: %d\n",NUMFREQ);
    printf("Maximum number of profiles: %d\n", NUMPROF);
    printf("Number of structures processed: %d\n",NUMSTRUCTS);
 }
  printf("Total number of equivalence helix classes: %d\n",total-1);
  if (VERBOSE)
    print_all_helices(total);

  //  make_graph(marginals,max,id,total,argv[1],fp);

  mostfreq = find_freq(total);
  printf("Total number of frequent helices: %d\n",hashtbl_numkeys(freq));
  make_cluster(argv[1],mostfreq);
  
  if (VERBOSE) {
    fp = fopen("cluster.out","w");
    fprintf(fp,"Processing %s\n",argv[1]);
    print_cluster(fp);
    fclose(fp);
  }
  printf("Total number of clusters: %d\n",hashtbl_numkeys(cluster)-1);
  
  fp = fopen(OUTPUT,"w");
  insert_graph(fp);  
  fputs("}",fp);
  
  hashtbl_destroy(max);
  hashtbl_destroy(marginals);
  hashtbl_destroy(idhash);
  hashtbl_destroy(binary);
  hashtbl_destroy(freq);
  hashtbl_destroy(cluster);
  hashtbl_destroy(bracket);
  fclose(fp);
  return 0;
}
