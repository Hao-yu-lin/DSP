#include "../inc/hmm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double *Pi;
double (*A)[MAX_STATE];
double (*B)[MAX_STATE];

static void cal_delta(HMM *hmm, double delta[MAX_LINE][MAX_STATE], int observ[MAX_LINE], int num_time, int num_state);

void test(HMM *hmm, int num_model, FILE *test_fp, FILE *result_fp) {
    char test_seq[MAX_LINE];
    int observ[MAX_LINE];
    double delta[MAX_LINE][MAX_STATE];
    double now_prob;
    double max_prob;
    int num_state;
    int num_time;
    int len_seq;
    int model_idx;


    while (fscanf(test_fp, "%s", test_seq) > 0) {
        model_idx = 0;
        max_prob = 0;

        for (int i = 0; i < num_model; i++) {
            num_state = hmm[i].state_num;
            num_time = strlen(test_seq);
            len_seq = strlen(test_seq);

            Pi = hmm[i].initial;
            A = hmm[i].transition;
            B = hmm[i].observation;
            memset(delta, 0, sizeof(delta));
            

            for (int t = 0; t < len_seq; t++) {
                observ[t] = test_seq[t] - 'A';
            }
            
            // init delta
            for (int j = 0 ; j < num_state ; j++) {
                delta[0][j] = Pi[j] * B[observ[0]][j];
            }

            cal_delta(&hmm[i], delta, observ, num_time, num_state);
            // prob reset
            now_prob = 0;
            
            for (int j = 0 ; j < num_state ; j++) {
                if (delta[num_time - 1][j] > now_prob) {
                    now_prob = delta[num_time - 1][j];
                }
            }
            if (now_prob > max_prob) {
                max_prob = now_prob;
                model_idx = i;
            }
            
        }
        fprintf(result_fp, "%s %e\n", hmm[model_idx].model_name, max_prob);
    }
}

int main(int argc, char *argv[]){
    char *models_list_path = argv[1];
	char *seq_path = argv[2];
	char *output_result_path = argv[3];

    HMM hmm[5];
	load_models(models_list_path , hmm, 5); 

	FILE *fp = fopen(seq_path, "r");    
	FILE *ans = fopen(output_result_path, "w");
    
	test(hmm, 5, fp, ans);
    
    fclose(ans);
    fclose(fp);
    

    return 0;
}

static void cal_delta(HMM *hmm, double delta[MAX_LINE][MAX_STATE], int observ[MAX_LINE], int num_time, int num_state){
    int max_idx;
    double max;
    for (int t = 1 ; t < num_time ; t++) {
        for (int j = 0 ; j < num_state ; j++) {
            max_idx = 0;
            max = 0;
            for (int k = 0 ; k < num_state ; k++) {
                if (delta[t-1][k] * A[k][j] > max) {
                    max = delta[t-1][k] * A[k][j];
                    max_idx = k;
                }
            }
            delta[t][j] = delta[t-1][max_idx] * A[max_idx][j] * B[observ[t]][j];
        }
    }
}
