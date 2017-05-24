#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_SIZE 256

typedef int bool;
typedef struct fastafile_s {
	FILE *ffd;
	char buffer[MAX_SIZE];
} fastafile;

int terror(int i);
bool is_fasta_header(char *buffer);
bool is_valid(char c);

int main(int argc, char **argv) {
	fastafile *ff1;
	int alloc_multiplier = 1;
	char *aux_pchar;
	char *sequence;
	char *name;
	int used;

	ff1 = malloc(sizeof(fastafile));

	ff1->ffd = fopen(argv[1], "r");

	if (ff1->ffd == NULL) {
		terror(1);
	}
	int aux;
	aux = feof(ff1->ffd);
	printf("%d\n", aux);

	while (!aux){
		if ((aux_pchar = fgets(ff1->buffer, MAX_SIZE, ff1->ffd)) == NULL)
			terror(2);
		if (is_fasta_header(aux_pchar)){
			name = malloc(sizeof(char) * strlen(aux_pchar));
			strcpy(name, aux_pchar+1);
			sequence = malloc(sizeof(char) * MAX_SIZE * alloc_multiplier);
			used = 0;
		} else {
			for (aux_pchar = ff1->buffer; *aux_pchar != '\0'; aux_pchar++){
				if (is_valid(*aux_pchar)){
					sequence[used++] = *aux_pchar;
					if (used == alloc_multiplier * MAX_SIZE){
						if (realloc(sequence, MAX_SIZE * ++alloc_multiplier) == NULL)
							terror(3);
					}
				}
			}
		}
		sequence[used] = '\0';
		aux = feof(ff1->ffd);
	}
	// printf("Name of the sequence: %s\nSequence:\n", name);
	// printf("%s\n", sequence);

	free(sequence);
	free(name);
	fclose(ff1->ffd);
	free(ff1);
	return(0);
}

int terror(int i){
	fprintf(stderr, "Error: %d\n",i);
	return(i);
}

bool is_fasta_header(char *buffer){
	return buffer[0] == '>';
}

bool is_valid(char c){
	return c == 'A'|| c == 'C' || c == 'G' || c == 'T';
}
