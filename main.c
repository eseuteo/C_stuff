// g_prof
// medir tiempo
// linked list
#include "fasta_reader.h"
#include <math.h>

int four_powers[3][16] =
{
  {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824},
  {2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152, 8388608, 33554432, 134217728, 536870912, 2147483648},
  {3, 12, 48, 192, 768, 3072, 12288, 49152, 196608, 786432, 3145728, 12582912, 50331648, 201326592, 805306368, 3221225472}
};

int main(int argc, char ** argv){
  if (argc < 3)
    terror(10);

  int k_size = atoi(argv[1]);
  if (k_size < 2 || k_size > 15){
    terror(11);
  }

  fastafile * fasta_f;
  frequencies_vector frequencies_v[2];
  int num_of_sequence;
  int index;
  char * k_mer;
  int nucleotid;
  char * aux_pchar;
  char * name;
  int aux;
  int empezar;
  int i;

  k_mer = malloc(sizeof(char) * k_size);
  aux_pchar = malloc(sizeof(char) * k_size);
  fasta_f = malloc(sizeof(fastafile));

  for (i = 2; i <= 3 ; i++){
    fasta_f->ffd = fopen(argv[i], "r");
    if (fasta_f->ffd == NULL){
    terror(15);
    }

    num_of_sequence = i-2;
    empezar = 1;

    frequencies_v[num_of_sequence].vector = calloc(four_powers[1][k_size], sizeof(int));
    strcpy(frequencies_v[num_of_sequence].name, "initialization");

    do {
      nucleotid = fgetc(fasta_f->ffd);
      if (nucleotid == '\n')
        continue;
      if (nucleotid == '>'){
        int j = 0;
        while (nucleotid != '\n'){
          frequencies_v[num_of_sequence].name[j] = nucleotid;
          nucleotid = fgetc(fasta_f->ffd);
          j++;
        }
        frequencies_v[num_of_sequence].name[j] = '\0';
      } else {
        if (is_valid(nucleotid)){
          if (empezar){
            k_mer[0] = nucleotid;
            for (int j = 1; j <= k_size-1; j++){
              do {
                nucleotid = fgetc(fasta_f->ffd);
              } while (nucleotid == '\n');
              if (!is_valid(nucleotid)){
                printf("El k-mer es demasiado grande ?\n");
              }
              k_mer[j] = nucleotid;
            }
            empezar = 0;
          } else {
            strncpy(aux_pchar, k_mer+1, k_size-1);
            strncpy(k_mer, aux_pchar, k_size-1);
            k_mer[k_size-1] = nucleotid;
          }
          index = get_index(k_mer, k_size);
          frequencies_v[num_of_sequence].vector[index]++;
        } else {
          empezar = 1;
        }
      }
    } while (nucleotid != EOF);
    num_of_sequence++;
  }

  double euclidean_d = 0;
  euclidean_d = euclidean_distance(frequencies_v[0].vector, frequencies_v[1].vector, (int) pow(4, k_size));
  printf("%f\n\n", euclidean_d);

  for (int i = 0; i<num_of_sequence; i++){
    printf("%s\n", frequencies_v[i].name);
    for (int j=0; j<(int) pow(4, k_size); j++){
      printf("%d\n", frequencies_v[i].vector[j]);
    }
    printf("\n");
  }

  free(fasta_f);
  free(k_mer);
  free(aux_pchar);
  free(frequencies_v[0].vector);
  free(frequencies_v[1].vector);
}

double euclidean_distance(int * vector1, int * vector2, int k){
  double res = 0;
  for (int i=0; i<k; i++){
    res += pow((vector1[i]-vector2[i]), 2);
  }
  return sqrt(res);
}

bool is_valid(char c){
	return c == 'A'|| c == 'C' || c == 'G' || c == 'T';
}

int get_index(char * k_mer, int k){
  int res = 0;
  int i;
  int nucleotid = 0;
  for (i = 0; i < k; i++){
    switch (k_mer[i]) {
      case 'A':
        continue;
      case 'C':
        nucleotid = 0;
        break;
      case 'G':
        nucleotid = 1;
        break;
      case 'T':
        nucleotid = 2;
        break;
      }
      res += four_powers[nucleotid][k-1-i];
  }
    return res;
}

int terror(int i){
	fprintf(stderr, "Error: %d\n",i);
	return(i);
}
