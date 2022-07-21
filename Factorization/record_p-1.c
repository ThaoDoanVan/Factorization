#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gmp.h"
#include "factormethod.h"

/*
Finding factor of a number n with at least k digits, using Pollard's p-1 method

run ./record_p-1 n B1 B2 k

B1, B2: bounds for each stage

using this script we found:

- 32 digit factor of 2^977 - 1
- 34 digit factor of 575th Fibonacci number 
- 66 digit factor of 960^119 - 1
*/

int main(int argc, char *argv[]){
    if ((argc == 4) || (argc == 5)){
        mpz_t n, d, B1, B2;
        mpz_inits(n, d, B1, B2, NULL);
        int k;
        if (argc == 4){

            mpz_set_str(n, argv[1], 10);
            mpz_set_str(B1, argv[2], 10);
            k = atoi(argv[3]);
            gmp_printf("Finding factor of n = %Zd with Pollard's p-1 method\nSmoothness bound B = %Zd\n\n", n, B1);
        }
        else {


            mpz_set_str(n, argv[1], 10);
            mpz_set_str(B1, argv[2], 10);
            mpz_set_str(B2, argv[3], 10);
            k = atoi(argv[4]);
            gmp_printf("Finding factor of n = %Zd with Pollard p-1 method in 2 phases\nPhase 1 smoothness bound B1 = %Zd\nPhase 2 smoothness bound B2 = %Zd\n\n", n, B1, B2);
        }
        mpz_t bound;
        mpz_init(bound);
        mpz_ui_pow_ui(bound, 10, k);
        clock_t start = clock();
        while(mpz_cmp(d, bound) < 0){
            
            int f = p_minus_1(d, n, B1, B2);
            

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
        mpz_clears(n, d, B1, B2, bound, NULL);
    }

    return 0;
}