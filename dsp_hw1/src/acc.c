#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
	FILE * acc_fp = fopen("acc.txt", "w");

	FILE * res_fp = fopen(argv[1], "r");
	FILE * ans_fp = fopen("./data/test_lbl.txt", "r");

	char res[20];
	char ans[20];

	int all = 0;
	int same = 0;

	while ( (fscanf(res_fp, "%s", res) != EOF) &&
			 fscanf(ans_fp, "%s", ans) != EOF) {
		if (strcmp(res, ans) == 0) {
			same++;
		}
		// printf("%s, %s\n", res, ans);
		all++;
		fscanf(res_fp, "%s", res);
	}
	double acc = (double)same/all;
	printf("%d, %d, %.6f\n", same, all, acc);
	fprintf(acc_fp, "%.6f", acc);
	fclose(acc_fp);
	fclose(res_fp);
	fclose(ans_fp);

	return 0;
}