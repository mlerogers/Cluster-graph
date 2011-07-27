/* Part of cluster graph
Program assigns all helices characterized by a maximal helix an ID and counts them
Finds all freq helices and clusters 1000 structs based on them
*/

#include "hashtbl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cluster.h"
#include <math.h>

static HASHTBL *triplet;
static HASHTBL *common;
//static char *comm;

/*
name is the file name containing structures to process
max is the hash of a hash: i->j->[k,id]
marginals is the hash: id->frequency of id
return idcount = num of helix equivalence classes + 1
*/
int process_structs(char *name) {
  FILE *fp;
  int i,j,k,helixid,idcount=1,m,n,o,longest=1,*lg,*id=NULL,last = 0;
  char tmp[100],key[ARRAYSIZE],*val,dbl[ARRAYSIZE];
  HASHTBL *temp;

  fp = fopen(name,"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",name);
    return 0;
  }
  triplet = hashtbl_create(HASHSIZE,NULL);
  while (fgets(tmp,100,fp) != NULL) {
    if (sscanf(tmp,"%d %d %d",&i,&j,&k) == 3) {
      //printf("%d %d %d\n",i,j,k);
      if (i == 0) continue;
      if (longest < k) 
	longest = k;
      helixid = testwithin(i,j,k,longest);
      //printf("id %d for %d %d %d\n",helixid,i,j,k);
      if (helixid == -1) { 
	insert(i,j,k,idcount);
	insertid(i,j,k,idcount);
	helixid = idcount++;
	//increment marginal
	sprintf(key,"%d",helixid);
	id = malloc(sizeof(int));
	*id = 1;
	hashtbl_insert(marginals,key,id);
	//triplet stats
	if (VERBOSE) {
	  temp = hashtbl_create(HASHSIZE,NULL);
	  sprintf(dbl,"%d %d",i,k);
	  lg = malloc(sizeof(int));
	  *lg = 1;
	  hashtbl_insert(temp,dbl,lg);
	  hashtbl_insert(triplet,key,temp);
	}
      }
      else {
	//printf("Found %d %d with id %d\n",i,j,helixid);
	sprintf(key,"%d",helixid);
	val = hashtbl_get(idhash,key);
	sscanf(val,"%d %d %d",&m,&n,&o);
	if (m > i || (m == i && o < k) ) {
	  //printf("replacing %d %d %d with %d %d %d\n",m,n,o,i,j,k);
	  insertid(i,j,k,helixid);
	}
	//increment marginal
	if (last != helixid) {
	  id = hashtbl_get(marginals,key);
	  ++*id;
	}
	//triplet stats
	if (VERBOSE) {
	  temp = hashtbl_get(triplet,key);
	  sprintf(dbl,"%d %d",i,k);
	  lg = hashtbl_get(temp,dbl);
	  if (lg)
	    (*lg)++;
	  else {
	    lg = malloc(sizeof(int));
	    *lg = 1;
	    hashtbl_insert(temp,dbl,lg);
	  }
	}
      }
      last = helixid;
    }
    else if (sscanf(tmp,"Structure %d",&i) == 1)
      NUMSTRUCTS++;
  }
  //hash doesn't alloc memory for vals, just for keys
  lg = malloc(sizeof(int));
  *lg = longest;
  hashtbl_insert(max,"longest",lg);
  
  //  printf("\n");
  if (fclose(fp))
    fprintf(stderr, "File %s not closed successfully\n",name);
  return idcount;   //idcount augmented, so can use < in loops
}

//checks if a longer or shorter helix already exists with ID
//if so, returns that ID
//inserts itself into hash, so future lookups wont have to check
int testwithin(int i, int j, int k,int longest) {
  int *id,m;
  if ((id = check(i,j)) != NULL) 
    return id[1];
  //puts("searching");
  for (m = 1; m <= longest; m++) {
    //if a shorter is found
    //printf("checking for m = %d\n",m);
    id = check(i+m,j-m);
    if (id != NULL && k >= m) {
      //      printf("Found shorter %d %d for %d %d\n",i+m,j-m,i,j);
      insert(i,j,k,id[1]);
      //insertid(i,j,k,id[1],idhash);
      return id[1];
    }
    id = check(i-m,j+m);
    //if a longer is found
    if (id != NULL && id[0] > m) {
      insert(i,j,k,id[1]);
      return id[1];
    }
  }
  return -1;
}

int* check (int i, int j) {
  HASHTBL *hash;
  char key[ARRAYSIZE],key2[ARRAYSIZE];
  int *val;

  sprintf(key,"%d",i);
  sprintf(key2,"%d",j);
  if ((hash = hashtbl_get(max,key)) != NULL) {
    if ((val = hashtbl_get(hash,key2)) != NULL)
      return val;
  }
  return NULL;
}

//This function populates hash with i->j->[k, idcount]
int insert (int i, int j, int k, int idcount) {
  char str[ARRAYSIZE], key[ARRAYSIZE];
  int *val;
  HASHTBL *tmp;
  //inserting into max
  val = malloc(sizeof(int)*2);
  val[0] = k;
  val[1] = idcount;
  sprintf(key,"%d",i);
  sprintf(str,"%d",j);
  if ((tmp = hashtbl_get(max,key)) != NULL) {
    hashtbl_insert(tmp,str,val);
    //    printf("found %d %d inserting %d\n",i,j,val[1]);
  }
  else {
    tmp = hashtbl_create(HASHSIZE,NULL);
    hashtbl_insert(tmp,str,val);
    hashtbl_insert(max,key,tmp);
    //    printf("inserting %d %d %d with id: %d\n",i,j,k,val[1]);
  }      
  return 0;
}

//inserts into hash with idcount->"i j k" if it is the maximal helix
//if found idcount in hash, assumes itself is larger and overwrites
int insertid (int i, int j, int k, int idcount) {
  char key[ARRAYSIZE], *get, *val;

  sprintf(key,"%d",idcount);
  if ((get = hashtbl_get(idhash,key)) != NULL) {
    sprintf(get,"%d %d %d",i,j,k);
    /*
    strcpy(get,newval);
    get = hashtbl_get(id,key);
    printf("inserting %s for id %d\n",get,idcount);
    */
  }
  else {
    val = malloc(sizeof(char)*ARRAYSIZE);
    sprintf(val,"%d %d %d",i,j,k);
    hashtbl_insert(idhash,key,val);
    //printf("inserting %s for id %d\n",val,idcount);
  }
  //  printf("inserting %d %d %d for id %d\n",i,j,k,idcount);

  return 0;
}

//looks up and prints actual helices for all id's
int print_all_helices(int total) {
  FILE *fp;
  HASHTBL *temp;
  KEY *node;
  char key[ARRAYSIZE],*val;
  int i,*m;

  fp = fopen("triplet.out","w");
  for (i = 1; i < total; i++) {
    sprintf(key,"%d",i);
    val = hashtbl_get(idhash,key);
    m = hashtbl_get(marginals,key);
    if (val != NULL)
      printf("Helix %d is %s with freq %d\n",i,val,*m);
    else
      printf("No entry for %d\n",i);

    //printing triplet info
    val = hashtbl_get(idhash,key);
    fprintf(fp,"For id %d with frequency %d, represented by %s:\n",i,*m,val);
    temp = hashtbl_get(triplet,key);
    if (!temp)
      fprintf(stderr,"error: no entry for %d: print_all_helices\n",i);
    else {
      for (node = hashtbl_getkeys(temp); node; node = node->next) {
	m = hashtbl_get(temp,node->data);
	fprintf(fp, "\t%s %d\n",node->data,*m);
      }
      
    }
      
  }
  fclose(fp);
  return 0;
}

//finds all frequent helices > threshold and prints them out
//inserts into hash with key = id, val = 1
char** find_freq(int total) {
  int *marg = NULL,length = 0,i,count = 0,percent,size = INIT_SIZE;
  double entropy;
  char key[ARRAYSIZE],**mostfreq;
  KEY *node = NULL;

  if (!(freq = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() of freq failed");
    exit(EXIT_FAILURE);
  }
  if (!(common = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() of common failed");
    exit(EXIT_FAILURE);
  }

  mostfreq = malloc(sizeof(char*)*NUMFREQ);
  //comm = malloc(sizeof(char)*ARRAYSIZE*size);
  //comm[0] = '\0';
  for (i = 1; i < total; i++) {
    sprintf(key,"%d",i);
    marg = hashtbl_get(marginals,key);
    //entropy = calc_entropy(*marg);
    //if (VERBOSE)
    //printf("entropy for %d is %.2f\n",i,entropy);
    if (!marg)
      fprintf(stderr,"error: marginal in find_freq is null");
    percent = *marg*100/NUMSTRUCTS;
    if (percent >= THRESH_COMMON) {
      if (VERBOSE)
	printf("Common helix %d with freq %d\n",i,*marg);
      hashtbl_insert(common,key,"1");
      //if (strlen(comm)+strlen(key) > (ARRAYSIZE*size-2)) 
      //comm = resize(&size,strlen(comm)+strlen(key)+1,comm);
      //sprintf(comm,"%s%s ",comm,key);
      //freq_insert(key,*marg,length++);
    }
    else if (percent >= THRESH_FREQ) {
      //printf("freq helix %s: with freq %d\n",key,*marg);
      if (FILTER) {
	filter(key,mostfreq,count++);
      } else {
	//if allow this option, need to sort
	freq_insert(key,*marg,length++);
      }
    }
  }
  //insert most freq helices (up to NUMFREQ)
  if (FILTER) {
    //add in common helices to array of most common; 
    if (count > NUMFREQ)
      count = NUMFREQ;

    //printf("count is now %d with numfreq %d\n",count,NUMFREQ);
    //for (percent = 0; percent < count; percent++)
    //printf("at %d is %s\n",percent,mostfreq[percent]);

    if ((i = hashtbl_numkeys(common)) > 0) {
      if (count+i > NUMFREQ)
	mostfreq = realloc(mostfreq,(count+i)*sizeof(char*));
      //printf("new size %d\n",count+i);
      
      i = count;
      count += hashtbl_numkeys(common);
      for (node = hashtbl_getkeys(common); node; node = node->next)
	mostfreq[i++] = node->data;
      //for (percent = 0; percent < count; percent++)
      //printf("at %d is %s\n",percent,mostfreq[percent]);
    }
    qsort(mostfreq,count,sizeof(char*),charcompare);
    for (i = 0; i < count; i++) {
      //printf("inserting %s\n",mostfreq[i]);
      marg = hashtbl_get(marginals,mostfreq[i]);
      freq_insert(mostfreq[i],*marg,length++);
    }
  }

    //sort it back for pruning purposes
  if (count == NUMFREQ)
    qsort(mostfreq,NUMFREQ,sizeof(char*),freqcompare);
  else
    for (i = count; i < NUMFREQ; i++)
      mostfreq[i] = "0";
  return mostfreq;
}

double calc_entropy(int marg) {
  double e = 0;
  double p = ((double)marg)/((double)NUMSTRUCTS);
  
  if (p == 1.0) return 0.0;
  //printf("marg is %d, total is %d and p is %.1f\n",marg,NUMSTRUCTS,p);
  e -=  p*log(p);
  p = 1-p;
  e -= p * log(p);
  return e;
}

int charcompare(const void *v1, const void *v2) {
  return (atoi(*((char**)v1))-atoi(*((char**)v2)));
}
//maintains most frequent helices, at most NUMFREQ
//opted for this rather than keep a long list of all freq helices, then sort and take top
//because long sequences may have very long list of frequent helices
void filter(char *key,char **mostfreq, int count) {
  int i,k,*least,*keyfreq;
  char *temp,*last;

  if (count < NUMFREQ) {
    mostfreq[count] = strdup(key);
    //printf("in mostfreq: %s as %d\n",mostfreq[count],count);
    if (count == NUMFREQ -1) {
      qsort(mostfreq,NUMFREQ,sizeof(char*),freqcompare);
    }
  } else {
    //if not among top 10 most freq, return
    least = hashtbl_get(marginals,mostfreq[NUMFREQ-1]);
    keyfreq = hashtbl_get(marginals,key);
    if (*keyfreq <= *least) return;
    //else, find where to insert...
    //    printf("testing key %s with freq %d\n",key,*((int*)hashtbl_get(marginals,key)));
    k = binsearch(mostfreq,key);
    if (*((int*)hashtbl_get(marginals,mostfreq[k])) > *((int*)hashtbl_get(marginals,key))) 
      k++;
    //printf("in mostfreq replacing %s with %s at %d\n",mostfreq[k],key,k);
    last = strdup(key);
    //..and bump everything down one
    for (i = k; i < NUMFREQ; i++ ) {
      temp = mostfreq[i];
      mostfreq[i] = last;
      last = temp;
    }
    free(last);
  }
}

//will sort to have descending freq
int freqcompare(const void *v1, const void *v2) {
  int *i1 =  hashtbl_get(marginals,*(char**)v1);
  if (!i1) fprintf(stderr,"%s not found in marginals\n",*(char**)v1);
  int *i2 =  hashtbl_get(marginals,*(char**)v2);
  if (!i2) fprintf(stderr,"%s not found in marginals\n",*(char**)v2);
  //  printf("i1 is %d and i2 is %d\n",*i1,*i2);
  return (*i2 - *i1);
}
//returns the index i next to where key should be (either side possible)
int binsearch(char **mostfreq, char *key) {
  int left = 0;
  int right = NUMFREQ - 1;
  int mid,val;
  char **l,**r;

  while (left < right) {
    mid = (left+right)/2;
    //printf("left is %d, right is %d, and mid %d\n",left,right,mid);
    l = &(mostfreq[mid]);
    r = &key;
    val = freqcompare(l,r);
    //    printf("val is %d\n",val);
    if ( val > 0)
      right = mid - 1;
    else if (val < 0)
      left = mid + 1;
    else
      return mid;
  }
  return left;
}

//inserts frequent helix and makes binary rep
void freq_insert(char *key,int marg,int length) {
  int *bin;
  char *val;

  //  printf("inserting into freq, key %s\n",key);
  hashtbl_insert(freq,key,"1");
  bin = malloc(sizeof(int));
  *bin = (1<<length);
  hashtbl_insert(binary,key,bin);
  //printf("inserting %s with value %d\n",key,(1<<(length-1)));
  val = hashtbl_get(idhash,key);
  if (val) {
    if (VERBOSE)
      printf("freq helix %s: %s with freq %d\n",key,val,marg);
  }
  else 
    printf("error: no triplet found for helix ID %s\n",key);
}

//like cluster.pl
//hash freq: key = helix ID. val = 1 (indicator)
//hash cluster: key = (num helices) list of helices. val = freq
HASHTBL* make_cluster(char *name,char **mostfreq) {
  FILE *fp;
  int num,id = 0,i,j,k,*longest=NULL,last = -1, iscommon = 0;
  int *count=NULL,most=0,numhelix = 0,size=INIT_SIZE;
  char temp[100],val[ARRAYSIZE],*l=NULL,*profile=NULL;
  HASHTBL *halfbrac;
  //KEY *node;

  longest = hashtbl_get(max,"longest");
  //  dup = malloc(sizeof(char)*ARRAYSIZE*size);
  profile = malloc(sizeof(char)*ARRAYSIZE*size);
//printf("length of profile is %d\n",(int)strlen(profile));
  profile[0] = '\0';

  if (!(cluster = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for cluster failed");
    exit(EXIT_FAILURE);
  }
  if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
    exit(EXIT_FAILURE);
  }
  if (!(bracket = hashtbl_create(HASHSIZE,NULL))) {
    fprintf(stderr, "ERROR: hashtbl_create() for bracket failed");
    exit(EXIT_FAILURE);
  }
  fp = fopen(name,"r");
  if (fp == NULL) {
    fprintf(stderr, "can't open %s\n",name);
    return 0;
  }
  while (fgets(temp,100,fp) != NULL) {
    if (sscanf(temp,"Structure %d",&num) == 1) {
      if (last == -1)
	continue;
      if (iscommon != hashtbl_numkeys(common)) {
	if (VERBOSE) 
	  printf("Found profile %snot having common helices\n",profile);
      }
      else
	profile = process_profile(halfbrac,profile,numhelix,&size,&most);
      if (VERBOSE && (count = hashtbl_get(cluster,profile)) && (*count == 1))
	printf("First struct %d with profile %s\n",num,profile);
      if (!(halfbrac = hashtbl_create(HASHSIZE,NULL))) {
	fprintf(stderr, "ERROR: hashtbl_create() for halfbrac failed");
	exit(EXIT_FAILURE);
      }
      last = 0;
      numhelix = 0;
      profile[0] = '\0';
      iscommon = 0;
    }
    else if (sscanf(temp,"%d%d%d",&i,&j,&k)!= EOF) {
      id = testwithin(i,j,k,*longest);
      if (id != -1 && id != last) {
	sprintf(val,"%d",id);
	if (hashtbl_get(common,val)) {
	  iscommon++;
	}
	l = hashtbl_get(freq,val);
	if (l != NULL) {   //is a freq helix, so save
	  numhelix++;
	  if (strlen(profile)+strlen(val) > (ARRAYSIZE*size-2)) 
	    profile = resize(&size,strlen(profile)+strlen(val)+1,profile);
	  //printf("adding %d to profile\n",id);
	  strcat(profile,val);
	  strcat(profile," ");
	  make_brackets(halfbrac,i,j,id);
	} 
	last = id;
      } 
    }
  }
  if (iscommon != hashtbl_numkeys(common)) {
    if (VERBOSE)
      printf("Found profile %snot having common helices\n",profile);
  }
  else
    profile = process_profile(halfbrac,profile,numhelix,&size,&most);
  for (; hashtbl_numkeys(cluster) > NUMPROF; prune_profiles(mostfreq))
    if (VERBOSE)
      printf("pruning profiles, now at %d profiles\n",hashtbl_numkeys(cluster));
    
  count = malloc(sizeof(int));
  *count = most;
  hashtbl_insert(cluster,"most",count);
  free(profile);
  //  hashtbl_destroy(halfbrac);
  return cluster;
}

//if profile is unique, insert into cluster, and make bracket representation
char* process_profile(HASHTBL *halfbrac,char *profile,int numhelix,int *size,int *most) {
  int *count,size2;
  char val[ARRAYSIZE],*dup;

  //printf("profile %s\n",profile);
  for (size2 = INIT_SIZE; strlen(profile) > (ARRAYSIZE *size2-1); size2++);  
  dup = malloc(sizeof(char)*ARRAYSIZE*size2);
  //numhelix += hashtbl_numkeys(common);
  sprintf(val,"%d ",numhelix);

  if (strlen(profile)+strlen(val)> ARRAYSIZE*(*size)-1)
    profile = resize(size,strlen(profile)+strlen(val)+1,profile);
  //  strcat(profile,comm);
  profile = strcat_front(profile,val);
  //printf("profile is %s, val is %s\n",profile,val);
  profile = quicksort(profile,dup);
  if ((count = hashtbl_get(cluster,profile)) == NULL) {
    count = malloc(sizeof(int));
    *count = 1;
    hashtbl_insert(cluster,profile,count);
    make_bracket_rep(halfbrac,dup);
    //printf("inserting %s for level %d\n",profile,numhelix-1);
  }
  else {
    ++*count;
    //printf("augmenting count of %s to %d\n",profile,*count-1);
  }
  if (*most < numhelix)
    *most = numhelix;
  free(dup);
  hashtbl_destroy(halfbrac);
  return profile;
}

//inserts bracket representation for i,j into a hash
void make_brackets(HASHTBL *brac, int i, int j, int id) {
  char key[ARRAYSIZE],*val;

  sprintf(key,"%d",i);
  val = malloc(sizeof(char)*ARRAYSIZE);
  sprintf(val,"[%d",id);
  //  printf("making bracket %s for %d\n",val,i);
  hashtbl_insert(brac,key,val);
  sprintf(key,"%d",j);
  val = malloc(sizeof(char)*2);
  val[0] = ']';
  val[1] = '\0';
  hashtbl_insert(brac,key,val);
}

//makes the bracket representation of dup, using values in hashtbl brac
//dup is a profile in graph
void make_bracket_rep(HASHTBL *brac,char *dup) {
  int num,*array,k=0,size = INIT_SIZE,total;
  char *profile,*val;
  KEY *node = NULL;

  num = hashtbl_numkeys(brac);
  array = malloc(sizeof(int)*num);
  for (node = hashtbl_getkeys(brac); node; node=node->next) 
    array[k++] = atoi(node->data);
  //sort by i,j position  
  qsort(array,num,sizeof(int),compare);
  profile = malloc(sizeof(char)*ARRAYSIZE*size);
  profile[0] = '\0';
  val = malloc(sizeof(char)*ARRAYSIZE);
  for (k = 0; k < num; k++) {
    sprintf(val,"%d",array[k]);
    val = hashtbl_get(brac,val);
    if ((total = strlen(profile)+strlen(val)) > ARRAYSIZE*size-1)
      profile = resize(&size,++total,profile);
    strcat(profile,val);
  }
  //if (VERBOSE)
  //printf("%s for %s\n",profile,dup);
  hashtbl_insert(bracket,dup,profile);
  free(val);
  free(array);
}

//resizes a dynamically allocated string s to appropriate size
char* resize(int *size,int total,char *s) {
  int old = *size;
  char *temp = NULL;
  for (; total > ARRAYSIZE * (*size);(*size)++) ;
  //if (old != *size)
  //printf("resizing to %d from %d for size %d\n",*size,old,total);
  temp = realloc(s,sizeof(char)*ARRAYSIZE*(*size));
  if (!temp)
    fprintf(stderr, "unable to resize %s\n",s);
  /*
  temp = malloc(sizeof(char)*ARRAYSIZE*(*size));
  strcpy(temp,s);
  free(s);
  s = temp;
  */
  return temp;
}

//implements quicksort to sort profile ID's into ascending order
char* quicksort(char *profile,char *dup) {
  int *array,i,length;
  char *blank = " ",*k=NULL,val[ARRAYSIZE],*copy = strdup(profile);

  k = strtok(copy,blank);
  length = atoi(k);
  array = malloc(sizeof(int)*length);
  for (i = 0; i < length; i++) {
    array[i] = atoi(strtok(NULL,blank));
    //printf("setting array[%d] to %d\n",i,array[i]);
  }
  qsort(array,length,sizeof(int),compare);
  dup[0] = '\0';
  for (i = 0; i < length; i++) {
    sprintf(val,"%d",array[i]);
    //printf("now adding %s\n",val);
    strcat(dup,val);
    strcat(dup," ");
  }
  //  strcat(k," ");
  //  profile = strcat_front(profile,k);
  sprintf(profile,"%s %s",k,dup);
  free(array);
  return profile;
}

int compare(const void *v1, const void *v2) {
  return (*(int*)v1 - *(int*)v2);
}

//concatenate ct to the front of s; return s
//assumes s has enough space to add ct
char* strcat_front(char *s, char *ct) {
  char *temp = mystrdup(s);
  strcpy(s,ct);
  strcat(s,temp);
  free(temp);
  return s;
  /*
  char *hold,*temp = malloc(sizeof(char)*(strlen(s)+strlen(ct)+1));
  //  printf("concatenating %s to %s\n",ct,s);
  strcpy(temp,ct);
  strcat(temp,s);
  hold = s;
  s = temp;
  printf("destroying %d and making new %d\n",hold,s);
  free(hold);
  return s;*/
}

//gets number of profiles to below NUMPROF
//eliminates profiles by removing least freq helix from consideration
void prune_profiles(char **mostfreq) {
  int k,*frq,*frq2,m;
  char *least,*profile;
  char *modprofile;
  KEY *node = NULL;
  HASHTBL *hash;

  for (k = NUMFREQ-1; !strcmp(mostfreq[k],"0"); k--) ;
  least = mostfreq[k];
  mostfreq[k] = "0";
  for (node = hashtbl_getkeys(cluster); node; node = node->next) {
    modprofile = malloc(strlen(node->data));
    profile = delete_helix(node->data,least,modprofile,&m);
    if (!strcmp(profile,node->data)) continue;
    frq = hashtbl_get(cluster,node->data);
    //modified profile exists or not
    if (!(frq2 = hashtbl_get(cluster,profile))) {
      //printf("adding %s t cluster\n",profile);
      hashtbl_insert(cluster,profile,frq);
      if (!(hash = hashtbl_create(HASHSIZE,NULL))) {
	fprintf(stderr, "ERROR: hashtbl_create() failed");
	exit(EXIT_FAILURE);
      }
      hashtbl_insert(hash,least,"1");
      edge_label(hash,modprofile,node->data,m);
    }
    else
      *frq2 += *frq;
    //printf("removing %s from cluster; num is %d\n",node->data,hashtbl_numkeys(cluster));
    hashtbl_remove(cluster,node->data);

    hashtbl_remove(bracket,modprofile);
  }
  
}

//takes cluster entry origprof and removes least if present
  char *delete_helix(char *origprof, char *least,char *modprofile,int *m) {
  int length,found = 0;
  char *k,*blank = " ";
  char *copy = mystrdup(origprof);
  char *profile = malloc(strlen(origprof));

  modprofile[0] = '\0';
  profile[0] = '\0';
  length = atoi(strtok(copy,blank));
  *m = length;
  for (k = strtok(NULL,blank); k; k = strtok(NULL,blank)) {
    sprintf(modprofile,"%s%s ",modprofile,k);
    if (strcmp(k,least))
      sprintf(profile,"%s%s ",profile,k);
    else
      found = 1;
  }
  if (found) {
    length--;
    //assuming num of helices in profile will never go past 4 digits
    k = malloc(sizeof(char)*5);
    sprintf(k,"%d ",length);
    profile = strcat_front(profile,k);
    if (VERBOSE)
      printf("profile %safter deleting %s is %s\n",origprof,least,profile);
    free(k);
    free(copy);
    return profile;
  }
  free(copy);
  return origprof;
}

int print_cluster(FILE *fp) {
  int *count, length,i;
  char *key, *blank = " ";
  KEY *node = hashtbl_getkeys(cluster);

  //first element is entry for most
  count = hashtbl_get(cluster,node->data);
  //printf("Longest chain of helices: %d\n",*count);
  fprintf(fp,"Using a threshold of %d:\n",THRESH_FREQ);
  for (node = node->next; node != NULL; node = node->next) {
    key = strdup(node->data);
    count = hashtbl_get(cluster,key);
    if (count)
      printf("Key %s has %d elements\n",key,*count);
    else
      printf("No entry for %s",key);

    length = atoi(strtok(key,blank));
    fprintf(fp,"Representative structure:\n");
    for (i = 0; i < length; i++) {
      fprintf(fp,"\t%s\n",(char*)hashtbl_get(idhash,strtok(NULL,blank)));
    }
    
  }
  return 0;
}
