#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../inc/hmm.h"


void train(HMM *hmm, FILE *train_fp, int iter);
void cal_alpha(double alpha[MAX_LINE][MAX_STATE], int observ[MAX_LINE], int num_times, int num_state);
void cal_beta(double beta[MAX_LINE][MAX_STATE], int observ[MAX_LINE], int T, int num_state);
void cal_gamma(double sum_gamma[MAX_STATE], double sum_gamma_observ[MAX_OBSERV][MAX_STATE], double alpha[MAX_LINE][MAX_STATE], double beta[MAX_LINE][MAX_STATE], int observ[MAX_LINE], double pi_reest[MAX_STATE], int num_times, int num_state);
void cal_epsilon(double sum_epsilon[MAX_STATE][MAX_STATE],  double alpha[MAX_LINE][MAX_STATE], double beta[MAX_LINE][MAX_STATE], int observ[MAX_LINE], int num_times, int num_state);

double *pi;
double (*A)[MAX_STATE];
double (*B)[MAX_STATE];

int main(int argc, char *argv[]){
	int iter = atoi(argv[1]);
	char *model_init_path = argv[2];
	char *seq_path = argv[3];
	char *output_model_path = argv[4];

	HMM hmm;
	loadHMM(&hmm, model_init_path);
	FILE *filp = open_or_die(seq_path, "r");


	train(&hmm, filp, iter);

	FILE* output_file = fopen(output_model_path, "w");
	dumpHMM(output_file, &hmm);
	fclose(filp);
	fclose(output_file);
	return 0;
}

void train(HMM *hmm, FILE *train_fp, int iter){
	// HMM_PARAM parameter;
	const int num_state = hmm->state_num;
	const int num_observ = hmm->observ_num;

	// // build pi & a & b parameter
	pi = hmm->initial;
	A = hmm->transition;
	B = hmm->observation;

	// read train_sequence
    char train_seq[MAX_LINE];
	int observ[MAX_LINE];
	
	// parameter
	double alpha[MAX_LINE][MAX_STATE];	
	double beta[MAX_LINE][MAX_STATE];
	double sum_gamma[MAX_STATE];
	double sum_gamma_observ[MAX_OBSERV][MAX_STATE];
	double sum_epsilon[MAX_STATE][MAX_STATE];
	double pi_reest[MAX_STATE];	

	int count = 0;
	while(count < iter) {
		int num_samples = 0;
        int num_times = 0;
        int len_seq = 0;
		int T = 0;
		memset(pi_reest, 0, sizeof(double) * MAX_STATE);
		memset(sum_gamma, 0, sizeof(double) * MAX_STATE);
		memset(sum_gamma_observ, 0, sizeof(double) * MAX_OBSERV * MAX_STATE);
		memset(sum_epsilon, 0, sizeof(double) * MAX_STATE * MAX_STATE);
		fseek(train_fp,  0,  SEEK_SET);
		
		while (fscanf(train_fp, "%s", train_seq) != EOF) {
			len_seq = strlen(train_seq);
			if (num_samples == 0)
				num_times = len_seq;
			
			for (int t = 0; t < len_seq; t++){
                observ[t] = train_seq[t] - 'A';
            }
			T = num_times - 1;
			// init alpha beta
			for (int j = 0 ; j < num_state ; j++){
    			alpha[0][j] = pi[j] * B[observ[0]][j];
				beta[T][j] = 1;
			}
			
			cal_alpha(alpha, observ, num_times, num_state);
			cal_beta(beta, observ, T, num_state);
			cal_gamma(sum_gamma, sum_gamma_observ, alpha, beta, observ, pi_reest, num_times, num_state);
			cal_epsilon(sum_epsilon, alpha, beta, observ, num_times, num_state);
			
			num_samples++;
			memset(alpha, 0, sizeof(alpha));
			memset(beta, 0, sizeof(beta));
		}
		// update parameters
		for (int j = 0 ; j < num_state ; j++){
		    pi[j] = pi_reest[j] / num_samples;
		}
		for (int j = 0 ; j < num_state ; j++) {
		    for (int k = 0 ; k < num_state ; k++) {
				A[j][k] = sum_epsilon[j][k] / sum_gamma[j];
		    }
	    }
		for (int ob = 0 ; ob < num_observ ; ob++) {
		    for (int j = 0 ; j < num_state ; j++) {
			    B[ob][j] = sum_gamma_observ[ob][j] / sum_gamma[j];
		    }
	    }


		count ++;
	}
}

void cal_alpha(double alpha[MAX_LINE][MAX_STATE], 
				int observ[MAX_LINE], 
				int num_times, 
				int num_state){
	double sum_alpha;
	for (int t = 1 ; t <= num_times ; t++) {
		for (int j = 0 ; j < num_state ; j++) {
			sum_alpha = 0;
			for (int i = 0 ; i < num_state ; i++)
				sum_alpha += alpha[t-1][i] * A[i][j];
			alpha[t][j] = sum_alpha * B[observ[t]][j];
		}
	}

}

void cal_beta(double beta[MAX_LINE][MAX_STATE], 
				int observ[MAX_LINE], 
				int T, 
				int num_state){
	double beta_sum;
	for(int t = T-1 ; t >= 0 ; t--) {
		for (int j = 0 ; j < num_state ; j++) {
			beta_sum = 0;
			for (int k = 0 ; k < num_state ; k++){
				beta_sum = beta_sum + A[j][k] * B[observ[t+1]][k] * beta[t+1][k];
			}
			beta[t][j] = beta_sum;
		}
	}
}

void cal_gamma(double sum_gamma[MAX_STATE], 
				double sum_gamma_observ[MAX_OBSERV][MAX_STATE], 
				double alpha[MAX_LINE][MAX_STATE], 
				double beta[MAX_LINE][MAX_STATE],
				int observ[MAX_LINE],
				double pi_reest[MAX_STATE],
				int num_times, 
				int num_state){
	double gamma[num_state];
	for (int t = 0 ; t < num_times ; t++) {
		double gamma_sum = 0;
		for (int j = 0 ; j < num_state ; j++){
			gamma[j] = alpha[t][j] * beta[t][j];			
			gamma_sum = gamma_sum +  gamma[j];			
		}
		for (int j = 0 ; j < num_state ; j++){
			gamma[j] = gamma[j]/gamma_sum;					
			sum_gamma[j] = sum_gamma[j] + gamma[j];		
			sum_gamma_observ[observ[t]][j] = sum_gamma_observ[observ[t]][j]+ gamma[j];
		}
		
		if (t == 0) {
			for (int j = 0 ; j < num_state ; j++) {
				pi_reest[j] = pi_reest[j] + gamma[j];
			}
		}
	}
}

void cal_epsilon(double sum_epsilon[MAX_STATE][MAX_STATE],
				double alpha[MAX_LINE][MAX_STATE], 
				double beta[MAX_LINE][MAX_STATE], 
				int observ[MAX_LINE], 
				int num_times, 
				int num_state){
	double epsilon[num_state][num_state];
	for (int t = 0 ; t < num_times ; t++) {
		double epsil_sum = 0;
		for (int j = 0 ; j < num_state ; j++) {
			for (int k = 0 ; k < num_state ; k++) {
				epsilon[j][k] = alpha[t][j] * A[j][k] * B[observ[t+1]][k] * beta[t+1][k];
				epsil_sum = epsil_sum + epsilon[j][k];									
			}
		}
		for (int j = 0 ; j < num_state ; j++) {
			for (int k = 0 ; k < num_state ; k++) {
				if (epsil_sum != 0) {
					epsilon[j][k] = epsilon[j][k]/epsil_sum;
				}		
				sum_epsilon[j][k] = sum_epsilon[j][k] + epsilon[j][k];
			}
		}
	}
}


