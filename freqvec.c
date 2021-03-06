#include <string.h>
#include <stdio.h>
#include <math.h>

#define MAX_LINE 2048

int get_index(char * k_mer, int k);
void freq_vector(int * frequency_vector, char * sequence, int k);

/*
void main(){
  int k = 5;
  int frequency_vector[(int) pow(4, k)];
  char sequence[MAX_LINE] = "GGCAGATTCCCCCTAGACCCGCCCGCACCATGGTCAGGCATGCCCCTCCTCATCGCTGGGCACAGCCCAGAGGGTATAAACAGTGCTGGAGGCTGGCGGGGCAGGCCAGCTGAGTCCTGAGCAGCAGCCCAGCGCAGCCACCGAGACACCATGAGAGCCCTCACACTCCTCGCCCTATTGGCCCTGGCCGCACTTTGCATCGCTGGCCAGGCAGGTGAGTGCCCCCACCTCCCCTCAGGCCGCATTGCAGTGGGGGCTGAGAGGAGGAAGCACCATGGCCCACCTCTTCTCACCCCTTTGGCTGGCAGTCCCTTTGCAGTCTAACCACCTTGTTGCAGGCTCAATCCATTTGCCCCAGCTCTGCCCTTGCAGAGGGAGAGGAGGGAAGAGCAAGCTGCCCGAGACGCAGGGGAAGGAGGATGAGGGCCCTGGGGATGAGCTGGGGTGAACCAGGCTCCCTTTCCTTTGCAGGTGCGAAGCCCAGCGGTGCAGAGTCCAGCAAAGGTGCAGGTATGAGGATGGACCTGATGGGTTCCTGGACCCTCCCCTCTCACCCTGGTCCCTCAGTCTCATTCCCCCACTCCTGCCACCTCCTGTCTGGCCATCAGGAAGGCCAGCCTGCTCCCCACCTGATCCTCCCAAACCCAGAGCCACCTGATGCCTGCCCCTCTGCTCCACAGCCTTTGTGTCCAAGCAGGAGGGCAGCGAGGTAGTGAAGAGACCCAGGCGCTACCTGTATCAATGGCTGGGGTGAGAGAAAAGGCAGAGCTGGGCCAAGGCCCTGCCTCTCCGGGATGGTCTGTGGGGGAGCTGCAGCAGGGAGTGGCCTCTCTGGGTTGTGGTGGGGGTACAGGCAGCCTGCCCTGGTGGGCACCCTGGAGCCCCATGTGTAGGGAGAGGAGGGATGGGCATTTTGCACGGGGGCTGATGCCACCACGTCGGGTGTCTCAGAGCCCCAGTCCCCTACCCGGATCCCCTGGAGCCCAGGAGGGAGGTGTGTGAGCTCAATCCGGACTGTGACGAGTTGGCTGACCACATCGGCTTTCAGGAGGCCTATCGGCGCTTCTACGGCCCGGTCTAGGGTGTCGCTCTGCTGGCCTGGCCGGCAACCCCAGTTCTGCTCCTCTCCAGGCACCCTTCTTTCCTCTTCCCCTTGCCCTTGCCCTGACCTCCCAGCCCTATGGATGTGGGGTCCCCATCATCCCAGCTGCTCCCAAATAAACTCCAGAAG";
  freq_vector(frequency_vector, sequence, k);
  int i;
  for (i=0; i<(int) pow(4,k); i++){
    printf("%d\t", frequency_vector[i]);
  }
}
*/

void freq_vector(int * frequency_vector, char * sequence, int k) {
  int vector_size = (int) pow(4, k);
  int length_seq = strlen(sequence);
  int i = 0;
  int index;
  char k_mer[k];

  for (i=0; i<vector_size; i++)
    frequency_vector[i] = 0;

  for (i=0; i<length_seq-(k-1); i++){
    strncpy(k_mer, sequence+i, k);
    index = get_index(k_mer, k);
    frequency_vector[index]++;
  }
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
