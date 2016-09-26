/* Permutation Test Test
 *
 * A test to test the permutation test by testing many permutations of the permutation test.
 *
 * Compilation and execution:
 * gcc -DDEBUG=0 -O2 -Wall permutation_test.c -o permutation_test -lgsl -lgslcblas -lpthread -lm && time ./permutation_test < input > output
 *
 * Mac hints:
 * You need GSL. Installed easily with Brew (http://brew.sh/).
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <malloc/malloc.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

// How many threads should we spawn.
#define NUMTHREADS 5

// How many samples should we generate per test.
#define NUMSAMPLES 1000

// How many tests should we run.
#define NUMTESTS 1000

// Maximum number of samples to expect.
#define MAX_SAMPLES 10000

// Return the smallest of two integers.
int min(int a, int b) {
	return a < b ? a : b;
}

// Try to get random seed from /dev/random if possible. Uses time as a fallback.
unsigned long int random_seed()
{
    unsigned int seed;
    struct timeval tv;
    FILE *devrandom;

    if ((devrandom = fopen("/dev/random","r")) == NULL) {
        gettimeofday(&tv,0);
        seed = tv.tv_sec + tv.tv_usec;
    } else {
        fread(&seed,sizeof(seed),1,devrandom);
        fclose(devrandom);
    }

    return(seed);
}

// SimulationData struct will hold all data needed for a series of simulation runs.
typedef struct {
	float observed_delta;		// What difference in means was observed.
	int num_observations;		// How many values were observed.
	int num_observations_base;	// How many values were observed in base.
	float * observations;		// What values were observed.
	
	int start, end;				// Which indexes of tests should this series run.
	unsigned long int seed;		// Seed to seed the random number generator.
} SimulationData;

// SimulationThread will run a series of simulation runs on a separate thread.
void *SimulationThread(void *args) {
	SimulationData * sd = (SimulationData*) args;
	float * observations = sd->observations;
	float observed_delta = sd->observed_delta;
	int num_observations = sd->num_observations;
	int num_observations_base = sd->num_observations_base;
	
	// Default random number generator is good enough for simulations and FAST.
	gsl_rng * r = gsl_rng_alloc (gsl_rng_default);
    gsl_rng_set(r, sd->seed);

	// Allocate memory for repeated sampling without replacement.
	float * o = malloc( num_observations_base * sizeof( float ) );

	// Calculate grand total so we can efficiently calculate means later.
	float grand_total = 0, base_total = 0;
	for (int i = 0; i < num_observations; i++) {
	    grand_total += observations[i];
	    if (i < num_observations_base) { base_total += observations[i]; }
	}

    // For each test
    for (int i = sd->start; i < sd->end; i++)
    {
    	// We will count how often permutated mean difference was larger than observed. 
    	int permutated_was_bigger = 0;
    	
	    // For each sample
	    for (int j = 0; j < NUMSAMPLES; j++)
	    {
		    // Choose random observations without replacement as new base.
		    // Equivalent to shuffling the entire list, but at least twice as fast.
		    gsl_ran_choose(r, o, num_observations_base, observations, num_observations, sizeof (*observations));
	        
		    // Calculate total for new base.
		    float base_perm = 0;
		    for (int i = 0; i < num_observations_base; i++) {
		        base_perm += o[i];
		    }
		    
			// Mean of new base is total over observations in base.
		    float mean_a = base_perm / num_observations_base;
		    // Mean of new variant is grand total minus new base total over observations in variant.
		    float mean_b = (grand_total - base_perm) / (num_observations - num_observations_base);
		    // (See what I did there? Saved at least half of the summing needed per iteration!)
	    	
	    	float simulated_delta = fabs(mean_a - mean_b);

			if (observed_delta < simulated_delta) permutated_was_bigger++;
	    }
	    
	    if (permutated_was_bigger == 0) {
	    	printf("%f\n", 0.0); // TODO Think about this edge case.
	    } else {
	    	// p-value is the fraction of times permutated mean was more extreme than observed.
	    	printf("%f\n", (float) permutated_was_bigger / NUMSAMPLES);
	    }
    }

	// Destroy stuff
    gsl_rng_free(r); free(o);
    
    return NULL;
}

int main (void)
{
    // Create a generator using defaults.
    gsl_rng_env_setup();
    gsl_rng * r = gsl_rng_alloc (gsl_rng_default);

#ifndef DEBUG
	// This is super duper slow; enable only when needed.
    gsl_rng_set(r, random_seed());
#endif
	
	float * a_samples = malloc( MAX_SAMPLES * sizeof( float ) );
	float * b_samples = malloc( MAX_SAMPLES * sizeof( float ) );
	for (int i = 0; i < MAX_SAMPLES; i++) { a_samples[i] = 0; b_samples[i] = 0; }
	int a_size = 0, b_size = 0;
	float a_total = 0, b_total = 0;

	// Read exactly two lines of float data
	float float_in;
	while( scanf("%f",&float_in) ) {
		a_samples[a_size] = float_in;
		a_total += float_in;
		a_size++;

		if ( getchar() == '\n' ) { break; }
	}	
	while( scanf("%f",&float_in) ) {
		b_samples[b_size] = float_in;
		b_total += float_in;
		b_size++;		
		
		int c = getchar();
		if ( c == '\n' || c == EOF ) { break; }
	}

	// TODO Ensure base is always the variant with the LOWEST number of observations.
	// (Because optimisation when sampling without replacement.)
	int num_observations = a_size + b_size;
	int num_observations_base = a_size;
	float * observations = malloc( num_observations * sizeof( float ) );
	for (int i = 0; i < num_observations; i++) {
		if (i < num_observations_base) {
			observations[i] = a_samples[i];
		}
		else {
			observations[i] = b_samples[i - a_size];
		}
	}
	float observed_delta = fabs((a_total / a_size) - (b_total / b_size));
	
	// Multithreaded from here.
	pthread_t * thread = malloc( NUMTHREADS * sizeof( pthread_t ) );
    SimulationData * sd = malloc( NUMTHREADS * sizeof( SimulationData ) );
	int i, tasks_per_thread = (NUMTESTS + NUMTHREADS - 1) / NUMTHREADS;

	// Divide the work.
    for (i = 0; i < NUMTHREADS; i++) {
        sd[i].start = i * tasks_per_thread;
        sd[i].end = (i+1) * tasks_per_thread;
        sd[i].seed = gsl_rng_get(r);
        
        sd[i].observed_delta = observed_delta;
		sd[i].num_observations = num_observations;
		sd[i].num_observations_base = num_observations_base;
		sd[i].observations = observations;
    }
    sd[NUMTHREADS-1].end = NUMTESTS;
    
    // Kick off threads.
    for (i = 0; i < NUMTHREADS; i++) {
        pthread_create(&thread[i], NULL, SimulationThread, &sd[i]);
    }
    
    // Wait for all to complete.
    for (i = 0; i < NUMTHREADS; i++) {
        pthread_join(thread[i], NULL);
    }
        
    free(thread); free(sd);
    // End of multithreading
    
    gsl_rng_free(r);
    
    // REPORT HUGE SUCCES
    return 0;
}