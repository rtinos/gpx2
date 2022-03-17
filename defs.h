/******************************************************************************\
*								 Definitions							 *
\******************************************************************************/
#include <iostream>
using namespace std;

const int max_dm_size = 20000; // ~1.5G
enum { WORST, BEST };

/* Global variables */
extern double *coord_x, *coord_y;				// coordinates of the cities
extern int **W;									// weight matrix for the ATSP
extern int n_cities;							// number of cities
extern char *prob_name;							// name of the file for the weight matrix

/* Function declaration */
void gpx(int **solution_blue, int **solution_red, double *fitness_offspring_blue, double *fitness_offspring_red);
int *aloc_vectori(int lines);
double *aloc_vectord(int lines);
int **aloc_matrixi(int lines , int collums);
void desaloc_matrixi(int **Matrix , int lines);
void read_problem(char* filename);
void rand_perm(int *inp, int *out, int size);
int weight(int i, int j);
