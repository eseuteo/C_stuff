#include "fasta_reader.h"

sequences * fasta_reader(fastafile * ff){
  nuc_sequences[2] ans;
  char * name;
  char * sequence;
  char * aux_pchar;
  int i;
  int alloc_multiplier;
  int used;

  ans = malloc(sizeof(nuc_sequences));
  i = 0;

  while (!feof(ff->ffd)){
    if ((aux_pchar = fgets(ff->buffer, MAX_SIZE, ff->ffd)) == NULL)
			terror(FAIL_FGETS);
    if (*ff->buffer == '>'){
      name = malloc(sizeof(char) * strlen(aux_pchar));
			strcpy(name, aux_pchar+1);
      ans[i].name = name;

      alloc_multiplier = 1;
      sequence = malloc(sizeof(char) * MAX_SIZE * alloc_multiplier);
      ans[i].sequence = sequence;
			used = 0;
      i++;
    } else {
      for (aux_pchar = ff->buffer; *aux_pchar != '\0'; aux_pchar++){
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
  }
}

bool is_valid(char c){
	return c == 'A'|| c == 'C' || c == 'G' || c == 'T';
}

int terror(int i){
	fprintf(stderr, "Error: %d\n",i);
	return(i);
}
