#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "factorutil.h"

/*
Pollard's p-1 method in 2 stages with smoothness bounds B1 and B2

run ./factor_p-1 n B1 B2 
*/

int main(int argc, char* argv[]){
    if (argc == 4){
        mpz_t n, B1, B2;
        mpz_inits(n, B1, B2);

        mpz_set_str(n, argv[1], 10);
        mpz_set_str(B1, argv[2], 10);
        mpz_set_str(B2, argv[3], 10);


        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_printf("Computing the factorization of %Zd by Pollard's p - 1 method\nStage 1 smoothness bound: %Zd\nStage 2 smoothness bound: %Zd\n\n", n, B1, B2);
        clock_t start = clock();

        int k = factorization_p_minus_1(&f, n, B1, B2);

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
        } else gmp_printf("Failed to factor %Zd by Pollard's p-1 method\n\n", n);
        
        printf("Process finished in %.9f secs \n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        free(f.prime_factors);
        free(f.exp);

        mpz_clears(n, B1, B2, NULL);
    }
    return 0;
}