// Takes alignments and writes the sequences in a FASTA file if Coverage and
// identity values are in the range set

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "commonFunctions.h"
#include "commonFunctions.c"
#define MAX_LENGTH 69
#define SEQ_BUFFER 100000
#define MIN_VALUE 50
#define MAX_VALUE 100
#define STATE_INIT 0
#define STATE_ID_0 1
#define STATE_ID_1 2
#define STATE_PERCENTAGE_COV 3
#define STATE_PERCENTAGE_ID 4
#define STATE_SEQS 5
#define STATE_SEQ0 6
#define STATE_SEQ1 7

FILE* fopen64(const char *filename, const char *type);
void init_args(int argc, char ** av, FILE ** input, FILE ** output, uint64_t * min_coverage, uint64_t * max_coverage, uint64_t * min_id, uint64_t * max_id);

int main(int argc, char ** av){
	FILE * input = NULL;
	FILE * output = NULL;
	uint64_t min_coverage = MIN_VALUE;
	uint64_t max_coverage = MAX_VALUE;
	uint64_t min_id = MIN_VALUE;
	uint64_t max_id = MAX_VALUE;

	init_args(argc, av, &input, &output, &min_coverage, &max_coverage, &min_id, &max_id);

	char * temp_seq_buffer = NULL;
	char * sequence0 = NULL;
	char * sequence1 = NULL;
	char * value = NULL;
	char c;
	uint64_t idx = 0, r = 0;
	uint64_t sequence0_length = 0, sequence1_length = 0;
	uint64_t coverage_value = 0, id_value = 0;
	uint64_t state = STATE_INIT;
	uint64_t value_length = 0;
	uint64_t line_count_seq_0 = 0, line_count_seq_1 = 0;

	if ((temp_seq_buffer = calloc(READBUF, sizeof(char))) == NULL) {
		terror("Could not allocate memory for read buffer");
	}

	if ((sequence0 = calloc(SEQ_BUFFER, sizeof(char))) == NULL) {
		terror("Could not allocate memory for read buffer");
	}

	if ((sequence1 = calloc(SEQ_BUFFER, sizeof(char))) == NULL) {
		terror("Could not allocate memory for read buffer");
  }

	if ((value = calloc(4, sizeof(char))) == NULL) {
		terror("Could not allocate memory for read buffer");
	}

	idx = READBUF + 1;
	sequence0[0] = '>';
	sequence1[0] = '>';

	c = buffered_fgetc(temp_seq_buffer, &idx, &r, input);

	while((!feof(input) || (feof(input) && idx < r))){
		switch (state) {
			case STATE_INIT :
				if (c == '(') {
					sequence0_length = 1;
					sequence1_length = 1;
					line_count_seq_0 = 0;
					line_count_seq_1 = 0;
					state = STATE_ID_0;
				}
				break;
			case STATE_ID_0 :
				if (isdigit(c)) {
					sequence0[sequence0_length++] = c;
				} else if (c == ',') {
					sequence0[sequence0_length++] = '\n';
					state = STATE_ID_1;
				}
				break;
			case STATE_ID_1 :
				if (isdigit(c)) {
					sequence1[sequence1_length++] = c;
				} else if (c == ')') {
					sequence1[sequence1_length++] = '\n';
					state = STATE_PERCENTAGE_COV;
				}
				break;
			case STATE_PERCENTAGE_COV :
				if (isdigit(c)) {
					value[value_length++] = c;
				} else if (c == '%') {
				//	printf("%" PRIu64 "\n", value_length);
					value[value_length++] = '\0';
					coverage_value = atoi(value);
					value_length = 0;
					// printf("percentage cov: %" PRIu64 "\n", coverage_value);
					// printf("decision (1: continue, 0: restart) --- %d\n", (coverage_value >= min_coverage && coverage_value <= max_coverage));
					state = (coverage_value >= min_coverage && coverage_value <= max_coverage) ? STATE_PERCENTAGE_ID : STATE_INIT;
				}
				break;
			case STATE_PERCENTAGE_ID :
				if (isdigit(c)) {
					value[value_length++] = c;
				} else if (c == '%') {
					value[value_length++] = '\0';
					id_value = atoi(value);
					value_length = 0;
					// printf("percentage id: %" PRIu64 "\n", id_value);
					// printf("decision (1: continue, 0: restart) --- %d\n", (id_value >= min_id && id_value <= max_id));
					// getchar();
					state = (id_value >= min_id && id_value <= max_id) ? STATE_SEQS : STATE_INIT;
				}
				break;
			case STATE_SEQS :
				if (c == 'D') {
					state = STATE_SEQ0;
				} else if (c == 'Q') {
					state = STATE_SEQ1;
				} else if (c == '(') {
					sequence0[sequence0_length] = '\0';
					sequence1[sequence1_length] = '\0';
					fprintf(output, "%s\n%s\n", sequence0, sequence1);
					sequence0_length = 1;
					sequence1_length = 1;
					line_count_seq_0 = 0;
					line_count_seq_1 = 0;
					state = STATE_ID_0;
				}
				break;
			case STATE_SEQ0 :
				if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') {
					if (line_count_seq_0 == MAX_LENGTH) {
						sequence0[sequence0_length++] = '\n';
						line_count_seq_0 = 0;
					}
					sequence0[sequence0_length++] = c;
					line_count_seq_0++;
				} else if (c == '\n') {
					state = STATE_SEQS;
				}
				break;
			case STATE_SEQ1 :
				if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') {
					if (line_count_seq_1 == MAX_LENGTH) {
						sequence1[sequence1_length++] = '\n';
						line_count_seq_1 = 0;
					}
					sequence1[sequence1_length++] = c;
					line_count_seq_1++;
				} else if (c == '\n') {
					state = STATE_SEQS;
				}
				break;
		}
		c = buffered_fgetc(temp_seq_buffer, &idx, &r, input);
  }

	if (coverage_value >= min_coverage && coverage_value <= max_coverage && id_value >= min_id && id_value <= max_id) {
		sequence0[sequence0_length] = '\0';
		sequence1[sequence1_length] = '\0';
		fprintf(output, "%s\n%s\n", sequence0, sequence1);
	}

	fclose(output);
	fclose(input);
	free(sequence0);
	free(sequence1);
	free(temp_seq_buffer);
	free(value);
}

void init_args(int argc, char ** av, FILE ** input, FILE ** output, uint64_t * min_coverage, uint64_t * max_coverage, uint64_t * min_id, uint64_t * max_id) {
    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "\tprograma0 -in input_file -out output_file [-mincov min_coverage] -maxcov max_coverage [-minid min_id] -maxid max_id\n");
            exit(1);
        }
        if(strcmp(av[pNum], "-in") == 0){
            *input = fopen64(av[pNum+1], "rt");
            if(input==NULL) terror("Could not open input file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *output = fopen64(av[pNum+1], "wt");
            if(output==NULL) terror("Could not open output file");
        }
        if(strcmp(av[pNum], "-mincov") == 0) {
            *min_coverage = (uint64_t) atoi(av[pNum+1]);
            if(*min_coverage < MIN_VALUE || *min_coverage > MAX_VALUE) terror("min_coverage must be between [50,100]");
        }
				if(strcmp(av[pNum], "-maxcov") == 0) {
            *max_coverage = (uint64_t) atoi(av[pNum+1]);
            if(*max_coverage < 52 || *max_coverage > MAX_VALUE) terror("max_coverage must be between [50,100]");
        }
				if(strcmp(av[pNum], "-minid") == 0) {
            *min_id = (uint64_t) atoi(av[pNum+1]);
            if(*min_id < MIN_VALUE || *min_id > MAX_VALUE) terror("min_id must be between [50,100]");
        }
				if(strcmp(av[pNum], "-maxid") == 0) {
            *max_id = (uint64_t) atoi(av[pNum+1]);
            if(*max_id < 52 || *max_coverage > MAX_VALUE) terror("max_id must be between [50,100]");
        }
        pNum++;
    }

    if(*input==NULL || *output==NULL ) terror("An input file, max coverage and max id values are required");
}
