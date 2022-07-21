#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "factorutil.h"

/*
Trial division method up to bound p_max

run ./factor_division n p_max 
*/

int main(int argc, char* argv[]){
    if (argc == 3){
        mpz_t n, p_max;
        mpz_inits(n, p_max, NULL);

        mpz_set_str(n, argv[1], 10);
        mpz_set_str(p_max, argv[2], 10);


        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_printf("Computing the factorization of %Zd by trial division up to %Zd\n\n", n, p_max);
        clock_t start = clock();

        int k = factorization_division(&f, n, p_max);

        clock_t end = clock();
        if (k == 2) gmp_printf("%Zd is a prime number\n\n", n);
        else if (k != -1){
            if (k == 1) printf("Fully factored !\n\n");
            else printf("Partially factored !\n\n");
            gmp_printf("%Zd = ", n);
            for (int i = 0; i < f.omega; i++){
                if (f.exp[i] != 1) gmp_printf("%Zd^%lu", f.prime_factors[i], f.exp[i]);
                else gmp_printf("%Zd", f.prime_factors[i], f.exp[i]);
                if (i != f.omega - 1) printf(" * ");
            }
            printf("\n\n");
        } else gmp_printf("Failed to factor %Zd by trial division up to %Zd\n\n", n, p_max);
        
        printf("Process finished in %.9f secs \n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        printf("\n");
        free(f.prime_factors);
        free(f.exp);

        mpz_clears(n, p_max, NULL);
    }
    return 0;
}