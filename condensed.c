/*Used with main.c,cluster.c in cluster graph
Finds freq helices in 1000 structs; finds their LCA
Creates Hasse diagram with these vertices
Condensed version of graph.c
*/

#include "hashtbl.h"
#include "cluster.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static HASHTBL *edges;
static HASHTBL **graph;
static HASHTBL *input;
static char **table;
static char **modprofileID;
static int *sums;
static int most;
static int graphsize;

//runs through nodes in cluster, finding their LCA's
int insert_graph(FILE *fp) {
  int *count,i,*freq,k = 0,numkeys,total;
  char **profileID,*profile;
  KEY *node;
  struct hashnode_s *begin;

  if (!(edges = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }

  fputs("digraph G {\n",fp);
  fprintf(fp,"\tlabel = \"%s\";\n",OUTPUT);
  fprintf(fp,"\tpad = 0.5;\n");
  fprintf(fp,"\tnodesep = 0.5;\n");
  //fputs("{ node [shape = plaintext]; ",fp);
 //first entry is longest profile length
  node = hashtbl_getkeys(cluster);
  count = hashtbl_get(cluster,node->data);
  graphsize = *count;
  //printf("graphsize is %d for %s\n",*count,node->data);
  graph = malloc(sizeof(HASHTBL*)*(graphsize+1));
  for (i = 0; i <= graphsize; i++) graph[i] = NULL;
  /*print ranks
  for (i = 0; i < *count; i++) {
    fprintf(fp,"%d",i+1);
    if (i == *count-1)
      fprintf(fp,"; }\n");
    else
      fprintf(fp,"->");
  }
  */
  numkeys = hashtbl_numkeys(cluster)-1;
  sums = malloc(sizeof(int)*numkeys);
  k = numkeys-1;
  profileID = malloc(sizeof(char*)*numkeys);
  modprofileID = malloc(sizeof(char*)*numkeys);
  most = 0;
  //for each profile, insert with freq
  
  for (node = node->next; node; node = node->next) {
    freq = hashtbl_get(cluster,node->data);
    //need to insert into graph
    //printf("node data is %s\n",node->data);
    if (freq) {
      if (most < *freq)
	most = *freq;
      profile = malloc(strlen(node->data)+1);
      sums[k] = insert_and_binary(node->data,profile,*freq);      
      profileID[k] = node->data;
      modprofileID[k--] = profile;
      fprintf(fp,"\"%s\" [shape = box];\n",profile);
    }
    else
      fprintf(stderr,"No entry for %s\n",node->data);
    //printf("for profile %s binary rep is %d with ID %d\n",node->data,sums[k+1],k+1);
  }
  make_key();
  
  begin = insert_LCAs(fp,profileID,numkeys);
  //printf("inserted LCA's with numkeys %d\n",numkeys);
  find_LCA_edges(fp,begin,numkeys);
  if (INPUT)
    process_input(fp);
  total = print_vertices(fp);
  printf("Total number of vertices: %d\n",total);
  if (INPUT)
    hashtbl_destroy(input);
  for (i = 0; i < *count; i++) {
    if (graph[i])
      hashtbl_destroy(graph[i]);
  }
  free(sums);
  free(table);
  hashtbl_destroy(edges);
  return 0;
}

//inserts key into graph and converts from rep struct to binary notation
int insert_and_binary(char *key,char *profile,int freq) {
  char *blank = " ",*helix,*copy = strdup(key);
  int sum = 0,i,length = atoi(strtok(copy,blank)),*val;
  HASHTBL *hash = graph[length-1];

  profile[0] = '\0';
  for (i = 0; i < length; i++) {
    helix = strtok(NULL,blank);
    strcat(profile,helix);
    strcat(profile," ");
    sum += *((int*) hashtbl_get(binary,helix));
  }
  //printf("profile is %s\n",profile);
  if (!hash) {
    if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() failed");
      exit(EXIT_FAILURE);
    }
    graph[length-1] = hash;
  }
  val = malloc(sizeof(int));
  *val = freq;
  hashtbl_insert(hash,profile,val);
  //printf("profile %sinserted with length %d and freq %d\n",profile,length,*val);
  return sum;
}

//finds the LCA of powerset of set of profiles and insert if necessary
//graph holds the final graph, table[freq ID] = ID, profileID[profileID] = profile
//sums[profileID] = binary rep, k = number of profiles
//returns beginning of linked list containing all unique LCA's found
struct hashnode_s* insert_LCAs(FILE *fp,char **profileID,int k) {
  int num,count,frq,*val,bit,s=0,j,*found,size = INIT_SIZE;
  unsigned int i;
  char *profile;
  struct hashnode_s *node = NULL,*temp = NULL;

  found = malloc(sizeof(int)*k);

  //for every subset of minimal elements...
  //checking i = 1,2 is redundant (only 1 element in bin rep), so skip
  //printf("k is %d\n",k);
  for (i = 3; i < (1<<k); i++) {
    count = 0;
    frq = 0;
    j = hashtbl_numkeys(binary);
    //printf("number of freq helices %d\n",j);
    num = (1<<j)-1;
    //    printf("original num is %d with j %d\n",num,j);
    //...take their LCA...
    for (bit = 0, j = i; j > 0 && bit < k ; j >>= 1, bit++)
      if ((1 & j) == 1) {
	found[count++] = bit;
	//printf("found profile %s for bit %d\n",profileID[bit],bit);
	//sprintf(key,"%s",table[bit]);
	if (!(val = hashtbl_get(cluster,profileID[bit]))) 
	  fprintf(stderr, "No such %s in cluster\n",profileID[bit]);
	frq += *val;
	num &= sums[bit];
	//printf("num is %d after & with %d from bit %d: %s\n",num,sums[bit],bit,profileID[bit]);
      }
    //printf("binary profile is %d, k is %d,frq is %d, and binary helices is %d\n",i,k,frq,num);

    //...and insert
    if (num != 0) {
      profile = malloc(sizeof(char)*ARRAYSIZE*size);
      profile = convert_binary(profile,num,&s);
      //printf("profile is now %s\n",profile);
      find_edges(fp,profile,found,count);
      //printf("after edges with LCA profile %s and i %d\n",profile,i);
      temp = insert_LCA(profile,num,s,frq);
      if (temp) {
	temp->next = node;
	node = temp;
      } else
	free(profile);
    }
  }
  free(found);
  //  puts("after free found");
  return node;
}

//inserts profile into graph, updates freq if already exists
//saves LCA into linked list with num = binary rep
//s = number of helices in profile
struct hashnode_s* insert_LCA(char *profile,int num,int s,int frq) {
  int size,*val,*data;
  char *key;
  HASHTBL *hash = NULL;
  struct hashnode_s *temp = NULL;
  //KEY *node;

  hash = graph[s-1];
  if (!hash) {
    if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
      fprintf(stderr, "ERROR: hashtbl_create() failed");
      exit(EXIT_FAILURE);
    }
    graph[s-1] = hash;
  }

  val = hashtbl_get(hash,profile);
  if (!val) {
    val = malloc(sizeof(int));
    *val = frq;
    hashtbl_insert(hash,profile,val);
    //    printf("inserted %s into graph\n",profile);
    //if LCA isn't an original, save; will have to compare them later
    //below assumes s is double digits or less
    for (size = 1; size < strlen(profile)+3; size++) ;
    key = malloc(sizeof(char)*ARRAYSIZE*size);
    sprintf(key,"%d %s",s,profile);
    if (!hashtbl_get(cluster,key)) {
      //printf("%s as new vertex with bit %d\n",profile,num);
      temp = malloc(sizeof(struct hashnode_s));

      data = malloc(sizeof(int));
      *data = num;
      temp->key = profile;
      temp->data = data;

      //for (;temp;temp = temp->next)
      //printf("node is %s with bit %d\n",temp->key,*((int*)temp->data));
    }
    free(key);
  } else if (*val < frq) {
    //printf("profile %s inserted with freq %d, old freq %d\n\n",profile,frq,*val);
    *val = frq;
  }
  return temp;
}


//table[freq id] = helix id
//produces table,where ith element is val
//where (1<<i) = hash_get(binary,val)
void make_key() {
  int count;
  KEY *node =  hashtbl_getkeys(binary);

  count = hashtbl_numkeys(binary);
  //for (count = 0; (*num & (1<<count)) == 0; count++) ;
  table = malloc(sizeof(char*)*count);
  for (count--; count >= 0; count--) {
    //printf("table[%d] id is %s\n",count,node->data);
    table[count] = node->data;
    node = node->next;
  }
}

//converts binary rep to string of helices (profile)
//returns number of helices in profile
//change count to reflect number of helices in profile
char* convert_binary(char *profile,int binary, int *count) {
  int k,size = INIT_SIZE,num=0;

  profile[0] = '\0';
  for (k = 0; binary > 0; binary >>= 1, k++) 
    if ((binary & 1) == 1) {
      if (strlen(profile)+strlen(table[k]) > ARRAYSIZE*size-2) 
	profile = resize(&size,strlen(profile)+strlen(table[k])+2,profile);
      strcat(profile,table[k]);
      strcat(profile," ");
      num++;
    }
  *count = num;
  //printf("convert binary's profile is %s for binary %d\n",profile,binary);
  return profile;
}

//processes input file containing structures of interest, in triplet form
//similar code to make_cluster: converts triplet to profile,
void process_input(FILE *fp) {
  HASHTBL *halfbrac;
  FILE *file;
  char temp[100],tmp[ARRAYSIZE],*profile,*fullprofile,*diff;
  int i,j,k,id,*longest,last,lastprob;
  int numhelix = 0,fullnum = 0,size = INIT_SIZE,size2 = INIT_SIZE,size3 = INIT_SIZE;

  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  if (!(input = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for input failed");
    exit(EXIT_FAILURE);
  }
  longest = hashtbl_get(max,"longest");
  profile = malloc(sizeof(char)*ARRAYSIZE*size);
  fullprofile = malloc(sizeof(char)*ARRAYSIZE*size2);
  diff = malloc(sizeof(char)*ARRAYSIZE*size3);
  profile[0] = '\0';
  fullprofile[0] = '\0';
  diff[0] = '\0';
  if (!(file = fopen(INPUT,"r")))
    fprintf(stderr,"Cannot open %s\n",INPUT);
  while (fgets(temp,100,file)) {
    if (sscanf(temp,"%d %d %d",&i,&j,&k) == 3) {
      id = testwithin(i,j,k,*longest);
      //printf("id is %d for %d %d %d\n",id,i,j,k);
      sprintf(tmp,"%d",id);
      if (id != -1) {
	if (id != last) {
	  if (hashtbl_get(freq,tmp)) {   //is a freq helix, so save
	    numhelix++;
	    if (strlen(profile)+strlen(tmp) > (ARRAYSIZE*size-2)) 
	      profile = resize(&size,strlen(profile)+strlen(tmp)+1,profile);
	    //printf("adding %d to profile\n",id);
	    sprintf(profile,"%s%s ",profile,tmp);
	  }
	  else {
	    if (strlen(diff)+strlen(tmp)+2 > ARRAYSIZE*size3)
	      diff = resize(&size3,strlen(diff)+strlen(tmp)+2,diff);	  
	    sprintf(diff,"%s%s ",diff,tmp);
	    
	    //printf("printing diff %s\n",diff);
	  }
	  fullnum++;
	  if (strlen(fullprofile)+strlen(tmp) > (ARRAYSIZE*size-2)) 
	    fullprofile = resize(&size2,strlen(fullprofile)+strlen(tmp)+2,fullprofile);
	  sprintf(fullprofile,"%s%s ",fullprofile,tmp);
	  //printf("helix %d added is %s\n",fullnum,tmp);
	  last = id;
	  make_brackets(halfbrac,i,j,id);
	}
	else //if (id == -1)
	  printf("helix %d %d %d is duplicate with id %d\n",i,j,k,id);
      } else 
	printf("helix %d %d %d doesn't exist\n",i,j,k);
      
    }
    else if (sscanf(temp,"Structure %d (%d)",&i,&j) == 2) {
      if (strlen(fullprofile) == 0) {
	lastprob = j;
	continue;
      } 
      //printf("profile is %s, fullprofile %s with diff %s\n\n",profile,fullprofile,diff);
      halfbrac = process_input_profile(fp,halfbrac,fullprofile,fullnum,profile,numhelix,diff,lastprob);
      numhelix = 0;
      fullnum = 0;
      profile[0] = '\0';
      fullprofile[0] = '\0';
      diff[0] = '\0';
      lastprob = j;
    }
  }
  //printf("input profile is %s with fullprofile %s and diff %s\n",profile,fullprofile,diff);
  halfbrac = process_input_profile(fp,halfbrac,fullprofile,fullnum,profile,numhelix,diff,lastprob);
  free(profile);
  hashtbl_destroy(halfbrac);

  //finds edges between centroids
  find_centroid_edges(fp);

}

//takes input helix and checks if in graph
//if so, changes node shape to hexagon
//if not, insert new vertex
HASHTBL* process_input_profile(FILE *fp,HASHTBL *brac,char *fullprofile, int fullnum,char *profile,int numhelix, char *diff, int prob) {
  HASHTBL *hash, *temp=NULL;
  char *diff1,*bracket,*difftrip;
  int *val,k1=0,k2=0;
  KEY *parent,*next;

  if (!(temp = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }

  make_bracket_rep(brac,fullprofile);
  hashtbl_destroy(brac);
  
  profile = sort_input(profile,numhelix);
  //printf("sorted profile is %s with fullprofile %s and diff %s\n",profile,fullprofile,diff);

  if ((hash = graph[numhelix-1]) && (hashtbl_get(hash,profile))) {
    if (numhelix == fullnum) {
      //printf("case 1\n");
      fullprofile = sort_input(fullprofile,numhelix);
    } else {
    //cannot use find_diff because fullprofile has helices not in table[]
      //puts("case 2");
      diff = insert_diff(temp,diff);
      bracket = edge_label(temp,profile,fullprofile,fullnum);      
      fprintf(fp,"\"%s\"-> \"%s\" [label =\" %s\\n%s \",fontsize=8,style = dotted];\n",profile,fullprofile,diff,bracket); 
    }
  }
  else {
    /*          
    if (numhelix == fullnum)
      puts("case 3");
    else
      puts("case 4");
    */
    for (parent = find_parents(profile); parent; parent = next) {
      diff1 = find_diff(temp,parent->data,profile,&k1,&k2);      
      if (numhelix != fullnum) {
	difftrip = insert_diff(temp,diff);
	diff1 = realloc(diff1,strlen(diff1)+strlen(difftrip)+4);
	//printf("for parent %s, diff1 is now %s, diff is %s and difftrip is %s\n",parent->data,diff1,diff,difftrip);
	sprintf(diff1,"%s\\n%s",diff1,difftrip);
	//printf("Diff is %s for parent %s of profile %s; diff %s for %s\n",diff1,parent->data,profile,difftrip,fullprofile);	
      }
      bracket = edge_label(temp,parent->data,fullprofile,fullnum);
      fprintf(fp,"\"%s\"-> \"%s\" [label =\" %s\\n%s \",fontsize=8,style = dotted];\n",parent->data,fullprofile,bracket,diff1); 
      next = parent->next;
      free(parent);
    }
  }
  //printf("%s has size %d and prob %d\n",fullprofile,fullnum,prob);
  fprintf(fp,"\"%s\" [shape = hexagon];\n",fullprofile);
  val = malloc(sizeof(int)*2);
  val[0] = fullnum;
  val[1] = prob;
  hashtbl_insert(input,fullprofile,val);
  hashtbl_destroy(temp);  

  if (!(brac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }
  return brac;
}

//find all parents of profile in graph; profile has only freq helices
//transitive reduction to the rescue again
KEY* find_parents(char *profile) {
  int i,b1,b2,N,num = 0,k = 0;
  HASHTBL *hash;
  KEY *node,*parent,*begin = NULL;

  begin = malloc(sizeof(KEY));
  begin->data = "";
  begin->next = NULL;
  b2 = make_binary(profile,&num);
  //printf("finding parents for %s of length %d and binary %d\n",profile,num,b2);
  for (i = 0; i < num-1; i++) {
    if (!(hash = graph[i])) continue;
    for (node = hashtbl_getkeys(hash); node; node = node->next) {
      b1 = make_binary(node->data,&k);
      //printf("investigating %sof size %d with binary %d\n",node->data,i+1,b1);
      N = b1 & b2;
      if (N == b1) {
	//printf("found parent %sof size %d with binary %d\n",node->data,i+1,b1);
	parent = malloc(sizeof(KEY));
	parent->data = node->data;
	parent->next = begin;
	begin = parent;
      }
    }
  }
  return begin;
}

//diff contains helix nums that are inserted into a hash
//to be used before edge_label
//add triplet info to diff
char* insert_diff(HASHTBL *temp,char *diff) {
  int size = INIT_SIZE;
  char *blank = " ", *copy,*k,*val;

  if (strlen(diff) == 0) return "";
  copy = malloc(sizeof(char)*ARRAYSIZE*size);
  copy[0] = '\0';
  for (k = strtok(diff,blank); k; k = strtok(NULL,blank)) {
    hashtbl_insert(temp,k,"1");
    //printf("insert_diff: inserting %s into hash\n",k);
    val = hashtbl_get(idhash,k);
    if (strlen(copy)+strlen(k)+strlen(val)+5 > ARRAYSIZE*size)
      copy = resize(&size,strlen(copy)+strlen(k)+strlen(val)+5,copy);	  
    if (strlen(copy)> 0) 
      strcat(copy,"\\n");
    sprintf(copy,"%s%s: %s",copy,k,val);
  }
  //free(diff);
  return copy;
}

//sorts input profile; code similar to quicksort()
char* sort_input(char *profile,int length) {
  int *array,i=0;
  char *blank = " ", *k,*copy = mystrdup(profile);

  array = malloc(sizeof(int)*length);
  for (k = strtok(copy,blank); k ; k = strtok(NULL,blank)) 
    array[i++] = atoi(k); 
  qsort(array,length,sizeof(int),compare);
  profile[0] = '\0';
  for (i = 0; i < length; i++) {
    sprintf(profile,"%s%d ",profile,array[i]);
  }
  //  printf("input profile sorted now %s\n",profile);
  free(copy);
  free(array);
  return profile;
}

//inserts centroids into graph
//finds and prints edges between them
void find_centroid_edges(FILE *fp) {
  int *i,*zero;
  KEY *node;
  HASHTBL *hash = NULL;

  zero = malloc(sizeof(int));
  *zero = 0;
  for (node = hashtbl_getkeys(input); node; node = node->next) {
    i = hashtbl_get(input,node->data);
    if (*i-1 > graphsize)
      *i = graphsize+1;
    //printf("inserting %s into graph[%d]\n",node->data,*i-1);
    if (!(hash = graph[*i-1])) {
      if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
	fprintf(stderr, "ERROR: hashtbl_create() for input failed");
	exit(EXIT_FAILURE);
      }
      graph[*i-1] = hash;
    }
    if (!hashtbl_get(hash,node->data))
      hashtbl_insert(hash,node->data,zero);    
  }
}

int print_vertices(FILE *fp) {
  int i,*val,size = 5,total = 0,size2=INIT_SIZE,*frq = NULL,zero = 0,start,end;
  char *rank,*v;;
  HASHTBL *hash;
  KEY *node = NULL;

  rank = malloc(sizeof(char)*ARRAYSIZE*size);
  v = malloc(sizeof(char)*ARRAYSIZE*size2);
  for (i = 0; !graph[i]; i++) ;
  start = i;
  for (node = hashtbl_getkeys(graph[i]); node; node = node->next)
    check_insert_edge(fp,"",node->data);
  //print ranks
  fputs("{ node [shape = plaintext]; ",fp);
  if (graph[graphsize])
    end = graphsize;
  else
    end = graphsize-1;
  for ( ; i <= end; i++) {
    fprintf(fp,"%d",i+1);
    if (i == end)
      fprintf(fp,"; }\n");
    else
      fprintf(fp,"->");
  }

  for (i = start; i <= end; i++) {
    if (!(hash = graph[i])) continue;
    //printf("printing level %d\n",i+1);
    node = hashtbl_getkeys(hash);
    sprintf(rank,"{ rank = same; %d;",i+1);
    for (; node; node = node->next) {
      if (strlen(node->data)+ 4 > ARRAYSIZE*size2-1) {
	v = resize(&size2,strlen(node->data)+5,v);
	//printf("resizing v of size %d to %d\n",strlen(v),size2);
      }
      sprintf(v,"%d %s",i+1,node->data);
      frq = hashtbl_get(cluster,v);
      if (!frq) 
	frq = &zero;
      sprintf(v," \"%s\";",node->data);
      val = hashtbl_get(hash,node->data);
      //      printf("found %s for %s\n",hashtbl_get(bracket,node->data),node->data);
      if (VERBOSE)
	printf("Vertex %swith frequency %d, originally %d\n",node->data,*val,*frq);
      if (strlen(rank)+strlen(v) > ARRAYSIZE*size-1) {
	//	printf("resizing rank %s and v %s of size %d to %d\n",rank,v,strlen(rank)+strlen(v),size);
	rank = resize(&size,strlen(rank)+strlen(v)+1,rank);
      }
      strcat(rank,v);
      //fprintf(fp,"\"%s\" [label = \"%s ",node->data,hashtbl_get(bracket,node->data));
      fprintf(fp,"\"%s\" [label = \"",node->data);
      if (*frq == most)
	fprintf(fp,"**");
      //fprintf(fp,"%s",hashtbl_get(bracket,node->data));
      fprintf(fp,"(%d/%d)",*frq,*val);
      if (*frq == most)
	fprintf(fp,"**");
      if (INPUT && (val = hashtbl_get(input,node->data)))
	fprintf(fp,"\\n(%d)",val[1]);
      fprintf(fp,"\",style=filled,color=black,fillcolor=grey%d]\n",(1000-*frq)/20+49);
      //      fprintf(fp,"\"%s\" [shape = box, label = \"%s (%d)\"];\n",node->data,*val);
    }
    fprintf(fp,"%s }\n",rank);
    total += hashtbl_numkeys(hash);
    //v]0] = '\0';
  }
  return total;
}

//edges[LCA ID] = profile ID
//profile is LCA found, found[] is array of contributing originals, count their number
int find_edges(FILE *fp,char *profile, int *found, int count) {
  int i;
  char *origprof;

  //for each contributing original profile...
  for (i = 0; i < count; i++) {
    origprof = modprofileID[found[i]];
    //printf("original profile is %s, and subprofile is %swith total length %d\n",origprof,profile,strlen(profile)+strlen(origprof));
    if (!strcmp(profile,origprof)) continue;
    check_insert_edge(fp,profile,origprof);
  }

  return 0;
}

//find all edges originating from generated LCA vertices
//node is beginning of LCA list, k is number of original profiles
void find_LCA_edges(FILE *fp,struct hashnode_s *node,int k) {
  int LCA,size = INIT_SIZE,i,count=0;
  char *profile;
  struct hashnode_s *last = NULL,*temp = NULL; 

  if (!node) return;
  profile = malloc(sizeof(char)*ARRAYSIZE*size);
  for (; node; node = node->next) {
    if (last) {
        free(last->data);
        free(last->key);
        free(last);
    }
    last = node;
    //printf("considering LCA %s with bin %d\n",(char*)node->key,*((int*)node->data));
    //check LCA against original profiles for edges
    for (i = 0; i < k; i++) {
      LCA = *((int*)node->data) & sums[i];
      if (LCA == 0) continue;
      profile = convert_binary(profile,LCA,&count);
      if (LCA != *((int*)node->data))
	check_insert_edge(fp,profile,node->key);
    }
    //    count = 0;
    if (!node->next) continue;
    //check LCA against other LCA's
    for (temp = node->next; temp; temp = temp->next) {
      //printf("considering %s and %s\n",node->key,temp->key);
      LCA = *((int*)node->data) & *((int*)temp->data);
      if (LCA == 0) continue;
      profile = convert_binary(profile,LCA,&count);
      if (LCA != *((int*)node->data)) 
	check_insert_edge(fp,profile,node->key);
      else if (LCA != *((int*)temp->data)) 
	check_insert_edge(fp,profile,temp->key);
    }
  }
  free(last->data);
  free(last->key);
  free(last);
  free(profile);
}

//checks whether edge exists and inserts if not
//ie insert if profile -> origprof doesn't exist yet
void check_insert_edge(FILE *fp,char *profile,char *origprof) {
  int k1=0,k2=0,*f1,*f2;
  double ratio;
  char *diff,*copy,*brac;
  HASHTBL *hash = NULL;

  /*
  hash = graph[k1-1];
  f1 = hashtbl_get(hash,profile);
  hash = graph[k2-1];
  f2 = hashtbl_get(hash,origprof);
  ratio = ((double)*f2)/((double)*f1);
  //  printf("ratio for %s and %s is %d/%d = %.2f\n",origprof,profile,*f2,*f1,ratio);
  */

  copy = malloc(sizeof(char)*(strlen(profile)+strlen(origprof)+4));
  //  if (strlen(profile)+strlen(origprof) > (ARRAYSIZE*size-4)) 
  //copy = resize(&size,strlen(profile)+strlen(origprof)+4,copy);
  sprintf(copy,"%s-> %s",profile,origprof);

  if (hashtbl_get(edges,copy)) return;

  hashtbl_insert(edges,copy,"1");
  if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() failed");
    exit(EXIT_FAILURE);
  }
  diff = find_diff(hash,profile,origprof,&k1,&k2);
  brac = edge_label(hash,profile,origprof,k2);
  //color=gray%.0f,100-(ratio*100);
  if (strlen(brac) > 1)
    fprintf(fp,"\"%s\"-> \"%s\" [label=\"  %s\\n%s  \",fontsize=8];\n",profile,origprof,brac,diff);
  //if (VERBOSE)
  //printf("inserting edge %s\n",copy);

  free(copy);
  free(diff);
  free(brac);
  hashtbl_destroy(hash);
}  

//finds helical difference between profile -> origprof
//k1 and k2 are num of helices in profile and origprof
//return difference with triplets
char* find_diff(HASHTBL *hash,char *profile, char *origprof, int *k1, int *k2) {
  int size = INIT_SIZE,b1,b2,xor,i;
  char *diff,*id;

  diff = malloc(sizeof(char)*ARRAYSIZE*size);
  diff[0] = '\0';
  b1 = make_binary(profile,k1);
  b2 = make_binary(origprof,k2);
  xor = b1 ^ b2;
  //printf("b1 is %d, b2 is %d, and xor is %d\n",b1,b2,xor);
  for (i = 0; xor > 0; xor>>=1,i++) {
    if ((xor & 1)==1) {
      hashtbl_insert(hash,table[i],"1");
      id = hashtbl_get(idhash,table[i]);
      if (strlen(diff)+strlen(table[i])+strlen(id)+3 > ARRAYSIZE*size)
	diff = resize(&size,strlen(diff)+strlen(table[i])+strlen(id)+3,diff);
      if (strlen(diff)>1)
	strcat(diff,"\\n");
      sprintf(diff,"%s%s: %s",diff,table[i],id);
      
      //strcat(diff,table[i]);
      //strcat(diff," ");
    }
  }
  //printf("diff is %s between %s and %s\n",diff,origprof,profile);

  return diff;
}

//finds the bracket notation for profile->origprof based on originating profile of size k
//difference between profile and origprof is in hash
//make sure origprof has a bracket rep
char* edge_label(HASHTBL *hash,char *profile, char *origprof,int k) {
  int i=0,count = 0,j=0,num=0,*diff,*save,ind=0,m=0;
  char *origbrac,*brac,*blank = "[]",*copy,*val;
  char **array;

  if (!(origbrac = hashtbl_get(bracket,origprof))) return "";
  //  if (strlen(profile) == 0) 
  copy = mystrdup(origbrac);
  brac = malloc(strlen(origbrac));
  //printf("finding edge label between %s and %s\n",origbrac,profile);
  array = malloc(sizeof(char*)*k);
  diff = malloc(sizeof(int)*hashtbl_numkeys(hash));
  save = malloc(sizeof(int)*hashtbl_numkeys(hash));
  //put helices of origprof into array; array in chron order
  for (val = strtok(copy,blank); val; val = strtok(NULL,blank)) {
    //printf("val is %s\n",val);
    if (i >= k) fprintf(stderr,"mismatch between k=%d and number of helices %d for %s: edge_label\n",k,i,origprof);
    array[i++] = val;
    //printf("val %s is at %d\n",val,i-1);
  }
  //save index in origprof of all diff helices; in ascending order
  for (count = 0; count < i; count++) {
    if (hashtbl_get(hash,array[count])) {
      //printf("saving %s at %d to index %d\n",array[count],count,j);
      diff[j++] = count;
    }
  }
  copy[0] = '\0';
  brac[0] = '\0';
  count = 0;
  j = -1;
  //i is index for origbrac, j is index for what level we need to match ']'
  //ind is index for level that increases for '[' and decreases for ']'
  //num is number of '[' encountered
  //count is index of different helix being proecessed
  //m is index for array of origprof helices, to print out ones not in diff
  for (i = 0; origbrac[i] != '\0'; i++) {
    //keep track of how many '['s
    if (origbrac[i] == '[') {
      ind++;
      if (diff[count] == num++) {
	count++;
	save[++j] = ind;
	//printf("\nsaving %d to j=%d\n",ind,j);
	strcat(copy,"{");
      }
      else {
	strcat(brac,"[");
	strcat(brac,array[m]);
	strcat(copy,"[");
      }
      m++;
    }
    //keep track of which level brackets are at
    else if (origbrac[i] == ']')
      //printf("\nchecking ind %d against %d at %d\n",ind,save[j],j);
      if (j >= 0 && save[j] == ind--) {
	j--;
	strcat(copy,"}");
      }
      else {
	strcat(brac,"]");
	strcat(copy,"]");
      }
  }
  if (!(hashtbl_get(bracket,profile))) {
    hashtbl_insert(bracket,profile,brac);
    //printf("new brac is %s for %s and %s with copy %s\n",brac,origprof,profile,copy);
  }
  free(array);
  free(diff);
  free(save);
  return copy;
}

int make_binary(char *profile,int *k) {
  int sum=0,*bin = NULL;
  char *blank = " ";
  char *copy = strdup(profile);
  char *helix;
  for (helix = strtok(copy,blank); helix; helix = strtok(NULL,blank)) {
    bin = hashtbl_get(binary,helix);
    if (!bin) continue;
    sum += *bin;
    (*k)++;
  }
  free(copy);
  return sum;
}

char* print_edge(KEY *node,char **table,int v,int* sum) {
  int count = 0,val;
  char *v1 = malloc(sizeof(char)*30);
  
  v1[0] = '\0';
  if (v == 0) {
    printf("%s-- %s\n",node->data,"E ");
    return v1;
  }

  for (val = (v & 1); v > 0; v >>= 1, count++, val = (v & 1))
    if (val == 1) {
      *sum += *((int*) hashtbl_get(binary,table[count]));
      strcat(v1,table[count]);
      strcat(v1," ");
    }
  
  printf("%s-- %s\n",node->data,v1);
  return v1;
}

/*
int main() {
  int i,*bin;
  HASHTBL *binary = hashtbl_create(16,NULL);
  KEY *node,*temp,*begin;
  char **vals = NULL,*val;
  LIST *start=NULL,*curr = NULL,*tmp = NULL;

  vals = malloc(sizeof(char*)*4);
  for(i = 0; i<4; i++) {
    val = malloc(sizeof(char));
    sprintf(val,"%d",i+1);
    vals[i] = val;
    bin = malloc(sizeof(int));
    *bin = (1<<i);
    hashtbl_insert(binary,val,bin);
  }

  node = malloc(sizeof(KEY));
  node->data = "1 3 4 ";
  curr = malloc(sizeof(LIST));
  curr->data = 13;
  start = curr;
  begin = node;

  temp = malloc(sizeof(KEY));
  temp->data = "1 2 ";
  tmp = malloc(sizeof(LIST));
  tmp->data = 3;
  curr->next = tmp;
  curr = tmp;
  node->next = temp;
  node = temp;

temp = malloc(sizeof(KEY));
  temp->data = "2 3 4 ";
  tmp = malloc(sizeof(LIST));
  tmp->data = 14;
  curr->next = tmp;
  curr = tmp;
  node->next = temp;
  node = temp;

  temp = malloc(sizeof(KEY));
  temp->data = "2 4 ";
  tmp = malloc(sizeof(LIST));
  tmp->data = 10;
  curr->next = tmp;
  curr = tmp;
  curr->next = NULL;
  node->next = temp;
  temp->next = NULL;

  print_condensed_graph(binary,begin,start);

  return 0;
}

//nodekey is current in outer key chain, nodelist for outer list chain
int print_condensed_graph_list(HASHTBL *binary,KEY *nodekey,LIST *nodelist) {
  int subset,*sum;
  char **table = NULL,*sub;
  //subkey is new key for subset, temp is current in inner key chain, end is last in key chain, start is begin of new chain, newpos is current in new chain
  KEY *subkey,*tempkey,*lastkey=NULL,*startkey,*newkey=NULL;
  //nxt is second in list chain, last is lst in list chain, sublist is new list
  LIST *sublist,*templist = NULL,*lastlist=NULL,*startlist,*newlist;

  table = make_key(binary);
  for (; nodelist->next; nodelist = nodelist->next, nodekey = nodekey->next) {
    startkey = NULL;
    startlist = NULL;
    if (lastlist)
      free(lastlist);
    if (lastkey) 
      free(lastkey);
    for (templist = nodelist->next,tempkey = nodekey->next; templist; templist = templist->next,tempkey = tempkey->next) {
      subset = nodelist->data & templist->data;
      sub = NULL;
      printf("considering %sand %swith subset %d\n",nodekey->data,tempkey->data,subset);
      sum = malloc(sizeof(int));
      if (subset != nodelist->data)
	sub = print_edge(binary,nodekey,table,subset,sum);
      if (subset != templist->data)
	sub = print_edge(binary,tempkey,table,subset,sum);

      lastkey = tempkey;
      lastlist = templist;

      if (!sub || sub[0] == '\0') continue;

      //make a new key and attach
      subkey = malloc(sizeof(KEY));
      subkey->data = sub;
      if (!startkey) 
	startkey = subkey;
      else 
	newkey->next = subkey;
      newkey = subkey;

      //make a new list and attach
      sublist = malloc(sizeof(LIST));
      sublist->data = subset;
      if (!startlist)
	startlist = sublist;
      else
	newlist->next = sublist;
      newlist = sublist;

    }

    //hook up end of old chain to beginning of new
    lastkey->next = startkey;
    lastlist->next = startlist;
    //save to free
    lastkey = nodekey;
    lastlist = nodelist;
  }
  free(table);
  return 0;
}
*/
