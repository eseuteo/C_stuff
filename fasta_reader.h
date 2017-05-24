#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_SIZE 256
#define FAIL_FGETS 2

typedef int bool;

typedef struct fastafile_s {
	FILE *ffd;
	char buffer[MAX_SIZE];
} fastafile;

typedef struct freq_vector_s {
	char name[MAX_SIZE];
	int * vector;
} frequencies_vector;

int terror(int i);
bool is_valid(char c);
//void read_fasta(nuc_sequences * sequences, fastafile * fasta_f);
int get_index(char * k_mer, int k);
double euclidean_distance(int * vector1, int * vector2, int k);
