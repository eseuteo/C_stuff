#include "fasta_reader.h"
#include <math.h>

#define K_MER_SIZE 15

int main(int argc, char ** argv){
  if (argc < 2)
    terror(10);

  fastafile * fasta_f;
  frequencies_vector frequencies_v[16];

  fasta_f = malloc(sizeof(fastafile));
  fasta_f->ffd = fopen(argv[1], "r");
  if (fasta_f->ffd == NULL){
    terror(15);
  }

  int num_of_sequence = 0;
  int index;
  char k_mer[K_MER_SIZE];
  char * aux_pchar;
  char * name;
  int aux = 0;

  while (aux_pchar = fgets(fasta_f->buffer, MAX_SIZE, fasta_f->ffd)){
  //  printf("%s\n", aux_pchar);
    if (feof(fasta_f->ffd)){
      printf("EOF reached\n");
      break;
    }
    if (*aux_pchar == '>'){
      frequencies_v[num_of_sequence].vector = calloc((int) pow(4, K_MER_SIZE), sizeof(int));
      // for (int i=0; i<(int) pow(4, K_MER_SIZE); i++){
      //   printf("%d\n",frequencies_v[0].vector[i]);
      // }
    //  strcpy(frequencies_v[num_of_sequence].name, aux_pchar+1);
    //  name = malloc(sizeof(char) * strlen(aux_pchar));
    //  strcpy(name, aux_pchar);
    //  frequencies_v[num_of_sequence].name = name;
      strcpy(frequencies_v[num_of_sequence].name, aux_pchar);
    //  printf("%s\n", frequencies_v[num_of_sequence].name);
      num_of_sequence++;
    } else {
      for (aux_pchar = fasta_f->buffer; *(aux_pchar+K_MER_SIZE) != '\0'; aux_pchar++){
        aux++;
      //  printf("%s\n", aux_pchar);
        strncpy(k_mer, aux_pchar, K_MER_SIZE);
        index = get_index(k_mer, K_MER_SIZE);
        frequencies_v[num_of_sequence-1].vector[index]++;
      }
    }
  }

  double euclidean_d = euclidean_distance(frequencies_v[1].vector, frequencies_v[2].vector, (int) pow(4, K_MER_SIZE));
  printf("%f\n", euclidean_d);
  // for (int i = 0; i<num_of_sequence; i++){
  //   for (int j=0; j<(int) pow(4, K_MER_SIZE); j++){
  //     printf("%d, ", frequencies_v[i].vector[j]);
  //   }
  //   printf("\n");
  // }
}

double euclidean_distance(int * vector1, int * vector2, int k){
  double res;
  for (int i=0; i<k; i++){
    res += pow((vector1[i]-vector2[i]), 2);
  }
  return sqrt(res);
}

int get_index(char * k_mer, int k){
  int res = 0;
  int i;
  int nucleotid;
  for (i = 0; i < k; i++){
    switch (k_mer[i]) {
      case 'A':
        nucleotid = 0;
        break;
      case 'C':
        nucleotid = 1;
        break;
      case 'G':
        nucleotid = 2;
        break;
      case 'T':
        nucleotid = 3;
        break;
      }
    res += pow(4,(k-1-i)) * nucleotid;
  }
    return res;
}

//   fastafile * fasta_f;
//   nuc_sequences sequences[MAX_SIZE];
//
//   fasta_f = malloc(sizeof(fastafile));
//   fasta_f->ffd = fopen(argv[1], "r");
//
//   if (fasta_f->ffd == NULL){
//     terror(15);
//   }
//
//   read_fasta(sequences, fasta_f);
// }
//
// void read_fasta(nuc_sequences * sequences, fastafile * fasta_f){
//   char * sequence;
//   char * name;
//   char * aux_pchar;
//   int different_sequences;
//   int used;
//   int alloc_multiplier = 1;
//
//   different_sequences = 0;
//
//   while (aux_pchar = fgets(fasta_f->buffer, MAX_SIZE, fasta_f->ffd)){
//     if (feof(fasta_f->ffd)){
//       printf("EOF reached\n");
//       break;
//     }
//     if (*aux_pchar == '>'){     // es el nombre
//     //  name = malloc(sizeof(char) * strlen(aux_pchar));
// 		//	strcpy(name, aux_pchar+1);
//       strcpy(sequences[different_sequences].name, aux_pchar+1);
//       printf("%s\n", sequences[different_sequences-1].name);
//       //printf("%s\n", sequences[different_sequences-1].name);
//       // alloc_multiplier = 1;
// 		  // used = 0;
//       // if (different_sequences){
//       //   sequences[different_sequences-1].sequence = sequence;
//       //   if (realloc(sequence, sizeof(char) * MAX_SIZE * alloc_multiplier) == NULL)
//       //      terror(3);
//       // }
//       // different_sequences++;
//     } else {
//       // for (aux_pchar = fasta_f->buffer; *aux_pchar != '\0'; aux_pchar++){
//       //   if (is_valid(*aux_pchar)){
//       //     sequence[used++] = *aux_pchar;
//       //     if (used == alloc_multiplier * MAX_SIZE){
//       //       if (realloc(sequence, MAX_SIZE * ++alloc_multiplier * sizeof(char)) == NULL)
//       // 			   terror(3);
//       //     }
//       //   }
//       // }
//       continue;
//     }
//     sequence[used] = '\0';
//     printf("%s\n", sequence);
//   }
// }
//
int terror(int i){
	fprintf(stderr, "Error: %d\n",i);
	return(i);
}
//
// bool is_valid(char c){
// 	return c == 'A'|| c == 'C' || c == 'G' || c == 'T';
// }
