#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "factorutil.h"

/*
Pollard's rho method with Brent's cycle detection algorithm, test power of 2 up to bound B

run ./factor_pollard_rho n B
*/

int main(int argc, char* argv[]){
    if (argc == 3){
        mpz_t n;
        mpz_init(n);
        long long int B;

        mpz_set_str(n, argv[1], 10);
        sscanf(argv[2], "%lli", &B);


        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_printf("Computing the factorization of %Zd by Pollard's rho using Brent's cycle detection\nTesting power of 2 up to %lli\n\n", n, B);
        clock_t start = clock();

        int k = factorization_pollard_rho(&f, n, B);

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
        } else gmp_printf("Failed to factor %Zd by Pollard's rho method up to %lli iterations\n\n", n, B);
        
        printf("Process finished in %.9f secs \n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        free(f.prime_factors);
        free(f.exp);

        mpz_clear(n);
    }
    return 0;
}