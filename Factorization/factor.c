#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "factorutil.h"

/*
Factorization using trial division, Pollard's rho method (Brent's cycle detection) and Pollard's p-1 method
p_max: bound for trial division
B: bound for power of 2 to be tested
B1, B2: bounds for Pollard's p-1
run ./factor n p_max B B1 B2
*/

int main(int argc, char* argv[]){
    if (argc == 6){
        mpz_t n, p_max, B1, B2;
        mpz_inits(n, p_max, B1, B2, NULL);
        long long int B;

        mpz_set_str(n, argv[1], 10);
        mpz_set_str(p_max, argv[2], 10);
        sscanf(argv[3], "%lli", &B);
        mpz_set_str(B1, argv[4], 10);
        mpz_set_str(B2, argv[5], 10);

        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_printf("Computing the factorization of %Zd by\nTrial division up to %Zd\nPollard's rho method (with Brent's cycle detection), testing with power of 2 up to %lli iterations\nPollard's p-1 method with stage 1 smoothness bound %Zd and stage 2 smoothness bound %Zd\n\n", n, p_max, B, B1, B2);
        clock_t start = clock();

        int k = factorization(&f, n, p_max, B, B1, B2);

        clock_t end = clock();
        if (k == 2) gmp_printf("%Zd is a prime number\n\n", n);
        else if (k != -1){
            if (k == 1) printf("Fully factored !\n\n");
            else printf("Partially factored !\n\n");
            gmp_printf("%Zd = ", n);
            for (int i = 0; i < f.omega; i++){
                if (f.exp[i] != 1) gmp_printf("%Zd^%d", f.prime_factors[i], f.exp[i]);
                else gmp_printf("%Zd", f.prime_factors[i], f.exp[i]);
                if (i != f.omega - 1) printf(" * ");
            }
            printf("\n\n");
        } else gmp_printf("Failed\n\n");
        
        printf("Process finished in %.9f secs \n\n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        free(f.prime_factors);
        free(f.exp);

        mpz_clears(n, p_max, B1, B2, NULL);
    }
    return 0;
}