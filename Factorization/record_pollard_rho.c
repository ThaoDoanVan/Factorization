#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gmp.h"
#include "factormethod.h"

/*
Finding factor of a number n with at least k digits, using Pollard's rho method

run ./record_p-1 n c B k

B: maximum power of 2 to be tested

c: specify the choice of function f(x) = x^2 + c mod n

recommend to choose c != 0, -2
*/


int main(int argc, char *argv[]){
    if ((argc == 5)){
        mpz_t n, d, c;
        long long B;
        int k;
        mpz_inits(n, d, c, NULL);

        mpz_set_str(n, argv[1], 10);
        mpz_set_str(c, argv[2], 10);
        sscanf(argv[3], "%lli", &B);
        k = atoi(argv[4]);

        gmp_printf("Finding factor of n = %Zd, using Pollard's rho method with Brent's cycle detection\nPolynomial f(x) = x^2 + %Zd\nTesting power of 2 up to %lli\n\n", n, c, B);

        mpz_t bound;
        mpz_init(bound);
        mpz_ui_pow_ui(bound, 10, k);
        clock_t start = clock();
        while(mpz_cmp(d, bound) < 0){
            
            int f = pollard_rho_brent_opt(d, n, c, B, 100);
            

            if (f == 0) gmp_printf("Found divisor: %Zd\n", d);
            else{
                printf("Failed !");
                break;
            }
            if (mpz_cmp(d, bound) > 0){
                int l = mpz_sizeinbase(d, 10);
                printf("Divisor has %d digits !\n", l);
                break;
            }
            mpz_divexact(n, n, d);
        }
        clock_t end  = clock();
        printf("Process finished in %f secs\n", (double)(end-start)/CLOCKS_PER_SEC);
        mpz_clears(n, d, c, bound, NULL);
    }

    return 0;
}