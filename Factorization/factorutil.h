#include "factormethod.h"

/*
prime_factors is the list of prime factors and exp is the list of the corresponding exponents in the factorization
omega is the number of prime factors (the function which counts the number of prime factors of a number n is denoted by omega)
*/

typedef struct factor {
    mpz_t *prime_factors;
    unsigned long *exp;
    int omega;
} factor;

/*
set b and e to be the values satisfying b^e = n and such that e is maximum
only use when n is a perfect power, which is tested using function mpz_perfect_power_p
*/

void perfect_power(mpz_t b, unsigned long int *e, mpz_t n){
    mpz_t x;
    mpz_init(x);
    unsigned long int k = mpz_sizeinbase(n, 2);
    for(unsigned long int i = 2; i <= k; i++){
        int r = mpz_root(x, n, i);
        if (r != 0){
            *e = i;
            mpz_set(b, x);
        }
    }
    mpz_clear(x);
    return;
}


/*
compute the factorization of a number n by trial division method
return 2 if n is a prime
return 1 if n is fully factored,
return 0 if n is partially factored, in this case the final element of f->divisor is not prime
return -1 if n cannot be factored using bound p_max
*/

int factorization_division(factor *f, mpz_t n, mpz_t p_max){
    // primality testing
    if (mpz_probab_prime_p(n, 10) != 0) return 2;

    if (f == NULL) f = (factor*) malloc(sizeof(factor));
    f->prime_factors = (mpz_t*) malloc(sizeof(mpz_t));
    f->exp = (unsigned long*) malloc(sizeof(unsigned long));
    f->exp[0] = 1;
    f->omega = 1;

    int r = 0;
    mpz_t N, p, p_min;

    mpz_init(N);

    // perfect power testing
    if (mpz_perfect_power_p(n)){
        perfect_power(N, &(f->exp[0]), n);
        mpz_init_set(f->prime_factors[0], N);
        // primality testing for base N
        if (mpz_probab_prime_p(N, 10) != 0){
            mpz_clear(N);
            return 1;
        }
    } else {
        mpz_set(N, n);
        mpz_init_set(f->prime_factors[0], N);
    }

    mpz_inits(p, p_min, NULL);
    mpz_set_ui(p_min, 1);

    while (mpz_cmp_ui(N, 1) > 0){
        // trial division by primes from p_min to p_max
        r = successive_division(p, N, p_min, p_max);
        if (r != 0) break; // fail

        // store the found prime factor p
        unsigned long exp_N = f->exp[f->omega-1];
        mpz_set(f->prime_factors[f->omega-1], p);
        f->exp[f->omega-1] = 0;

        // find the corresponding exponent
        while(mpz_divisible_p(N, p)){
            mpz_divexact(N, N, p);
            f->exp[f->omega-1]++;
        }

        f->exp[f->omega-1] *= exp_N;

        // store the remaining factor N/p^e only when different than 1
        if (mpz_cmp_ui(N, 1) != 0){
            f->omega++;
            f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->omega)*(sizeof(mpz_t)));
            f->exp = (unsigned long*) realloc(f->exp, (f->omega)*(sizeof(unsigned long)));
            mpz_init_set(f->prime_factors[f->omega-1], N);
            f->exp[f->omega-1] = exp_N;
        }

        // primality testing for the remaining factor of n
        if (mpz_probab_prime_p(N, 10) != 0){
            r = 1;
            break;
        }

        // start the trial division from the primes > p
        mpz_set(p_min, p);
    }
    
    // remaining factor is 1, complete factorization
    if (mpz_cmp_ui(N, 1) == 0) r = 1;
    // remaining factor is nontrivial and not prime, incomplete factorization
    else if ((r != 1) && (mpz_cmp(N, n) < 0)) r = 0;

    mpz_clears(N, p, p_min, NULL);
    return r;
}


/*
compute the factorization of a number n by Pollard's rho method
return 2 if n is a prime
return 1 if n is fully factored,
return 0 if n is partially factored, in this case the final element of f->divisor is not prime
return -1 if n cannot be factored using Pollard's rho method up to B iterations
*/
int factorization_pollard_rho(factor *f, mpz_t n, long long int B){
    
    // primality testing
    if (mpz_probab_prime_p(n, 10) != 0) return 2;

    if (f == NULL) f = (factor*) malloc(sizeof(factor));
    f->prime_factors = (mpz_t*) malloc(sizeof(mpz_t));
    f->exp = (unsigned long*) malloc(sizeof(unsigned long));
    f->omega = 1;


    mpz_t N, d, p, c;
    mpz_init_set(N, n);

    // factor out power of 2
    if (mpz_divisible_2exp_p(N, 1)){
        mpz_init_set_ui(f->prime_factors[0], 2);
        f->exp[0] = 0;
        while(mpz_divisible_2exp_p(N, 1)){
            mpz_divexact_ui(N, N, 2);
            f->exp[0]++;
        }
        // store the remaining factor N/2^e only when different than 1
        if (mpz_cmp_ui(N, 1) != 0){
            f->omega++;
            f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->omega)*(sizeof(mpz_t)));
            mpz_init_set(f->prime_factors[f->omega-1], N);
            f->exp = (unsigned long*) realloc(f->exp, (f->omega)*(sizeof(unsigned long)));
            f->exp[f->omega-1] = 1;
            // remaining factor is prime, complete factorization
            if (mpz_probab_prime_p(N, 10) != 0){
                mpz_clear(N);
                return 1;
            }
        }
        // case N is a power of 2 
        else {
            mpz_clear(N);
            return 1;
        }
    }
    // perfect power test
    if (mpz_perfect_power_p(N)){
        perfect_power(N, &(f->exp[f->omega-1]), N);
        mpz_init_set(f->prime_factors[f->omega-1], N);
        // primality testing for base N
        if (mpz_probab_prime_p(N, 10) != 0){
            mpz_clear(N); 
            return 1;
        }
    } else mpz_init_set(f->prime_factors[0], N); f->exp[0] = 1;

    mpz_inits(d, p, c, NULL);
    mpz_set(d, N);
    mpz_set_ui(p, 1);
    int r = 0;

    while (mpz_cmp_ui(N, 1) > 0){
        // start with x^2 + 1
        mpz_set_ui(c, 1);
        r = pollard_rho_brent(p, d, c, B);

        // if fail, start with x^2 + 2
        if (r != 0){
            mpz_set_ui(c, 2);
            r = pollard_rho_brent(p, d, c, B);
        }
        
        // fail for both choices of polynomial
        if (r != 0) break;

        // a factor p is found, perform primality test
        if (mpz_probab_prime_p(p, 10) != 0) {
            // store the prime factor
            unsigned long exp_N = f->exp[f->omega-1];
            mpz_set(f->prime_factors[f->omega-1], p);
            f->exp[f->omega-1] = 0;
            
            // compute the exponent
            while(mpz_divisible_p(N, p)){
                mpz_divexact(N, N, p);
                f->exp[f->omega-1]++;
            }

            f->exp[f->omega-1] *= exp_N;

            // store the remaining factor N/p^e
            f->omega++;
            f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->omega)*(sizeof(mpz_t)));
            f->exp = (unsigned long*) realloc(f->exp, (f->omega)*(sizeof(unsigned long)));
            mpz_init_set(f->prime_factors[f->omega-1], N);
            f->exp[f->omega-1] = exp_N;

            // set d to the remaining factor N/p^e
            mpz_set(d, N);

            // primality test for d
            if (mpz_probab_prime_p(d, 10) != 0){
                // remaining factor is prime, hence complete factorization
                r = 1;
                break;  
            }          
        } else mpz_set(d, p); // in case factor p is not prime, apply Pollard's rho for p to find a prime factor of n
    }

     
    // remaining factor is 1, complete factorization
    if (mpz_cmp_ui(N, 1) == 0) r = 1;
    // remaining factor is nontrivial and not prime, incomplete factorization
    else if ((r != 1) && (mpz_cmp(N, n) < 0)) r = 0;

    mpz_clears(N, d, p, c, NULL);
    return r;
}

/*
compute the factorization of a number n by Pollard's p-1 method
return 2 if n is a prime
return 1 if n is fully factored,
return 0 if n is partially factored, in this case the final element of f->divisor is not prime
return -1 if n cannot be factored using Pollard's p-1 method with smoothness bounds B1 and B2
if B2 <= B1, only stage 1 is performed, else both stages are performed
*/
int factorization_p_minus_1(factor *f, mpz_t n, mpz_t B1, mpz_t B2){
    // primality testing
    if (mpz_probab_prime_p(n, 10) != 0) return 2;


    if (f == NULL) f = (factor*) malloc(sizeof(factor));
    f->prime_factors = (mpz_t*) malloc(sizeof(mpz_t));
    f->exp = (unsigned long*) malloc(sizeof(unsigned long));
    f->exp[0] = 1;
    f->omega = 1;


    mpz_t N, d, p, c;
    mpz_init_set(N, n);

    // factor out power of 2
    if (mpz_divisible_2exp_p(N, 1)){
        mpz_init_set_ui(f->prime_factors[0], 2);
        f->exp[0] = 0;
        while(mpz_divisible_2exp_p(N, 1)){
            mpz_divexact_ui(N, N, 2);
            f->exp[0]++;
        }
        // store the remaining factor N/2^e only when different than 1
        if (mpz_cmp_ui(N, 1) != 0){
            f->omega++;
            f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->omega)*(sizeof(mpz_t)));
            mpz_init_set(f->prime_factors[f->omega-1], N);
            f->exp = (unsigned long*) realloc(f->exp, (f->omega)*(sizeof(unsigned long)));
            f->exp[f->omega-1] = 1;
            // remaining factor is prime, complete factorization
            if (mpz_probab_prime_p(N, 10) != 0){
                mpz_clear(N);
                return 1;
            }
        }
        // case N is a power of 2 
        else {
            mpz_clear(N);
            return 1;
        }
    }

    // perfect power test
    if (mpz_perfect_power_p(N)){
        perfect_power(N, &(f->exp[f->omega-1]), N);
        mpz_init_set(f->prime_factors[f->omega-1], N);
        // primality testing for base N
        if (mpz_probab_prime_p(N, 10) != 0){
            mpz_clear(N);
            return 1;
        }
    } else mpz_init_set(f->prime_factors[0], N); f->exp[0] = 1;
    
    mpz_inits(d, p, c, NULL);
    mpz_set(d, N);
    mpz_set_ui(p, 1);

    int r = 0;

    while (mpz_cmp_ui(N, 1) > 0){
        r = p_minus_1(p, d, B1, B2);
        if (r != 0) break;
        
        // a factor p is found, perform primality test
        if (mpz_probab_prime_p(p, 10) != 0){
            // store the prime factor
            unsigned long exp_N = f->exp[f->omega-1];
            mpz_set(f->prime_factors[f->omega-1], p);
            f->exp[f->omega-1] = 0;
            
            // compute the exponent
            while(mpz_divisible_p(N, p)){
                mpz_divexact(N, N, p);
                f->exp[f->omega-1]++;
            }

            f->exp[f->omega-1] *= exp_N;

            // store the remaining factor N/p^e
            f->omega++;
            f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->omega)*(sizeof(mpz_t)));
            f->exp = (unsigned long*) realloc(f->exp, (f->omega)*(sizeof(unsigned long)));
            mpz_init_set(f->prime_factors[f->omega-1], N);
            f->exp[f->omega-1] = exp_N;

            // set d to the remaining factor N/p^e
            mpz_set(d, N);

            // primality test for d
            if (mpz_probab_prime_p(d, 10) != 0){
                // remaining factor is prime, hence complete factorization
                r = 1;
                break;  
            }          
        } else mpz_set(d, p);  // in case factor p is not prime, apply Pollard's p-1 for p to find a prime factor of n
    }

    // remaining factor is 1, complete factorization
    if (mpz_cmp_ui(N, 1) == 0) r = 1;
    // remaining factor is nontrivial and not prime, incomplete factorization
    else if ((r != 1) && (mpz_cmp(N, n) < 0)) r = 0;
    // else r = -1 and N = n indicate that algorithm fails to find a nontrivial factor of n -> failure 

    mpz_clears(N, d, p, c, NULL);
    return r;
}


/*
Factorization by trial division, Pollard's rho method and Pollard's p-1 method combined
*/
int factorization(factor *f, mpz_t n, mpz_t p_max, long long int B, mpz_t B1, mpz_t B2){
    mpz_t N, d, p, c;
    
    int r = factorization_division(f, n, p_max);
    
    if ((r == 1) || (r == 2)) return r;
    else{
        mpz_init_set(N, f->prime_factors[f->omega-1]);
        // printf("Trial division failed ! Apply Pollard's rho method\n\n");
    }


    mpz_inits(d, p, c, NULL);
    mpz_set(d, N);
    mpz_set_ui(p, 1);

    while (mpz_cmp_ui(N, 1) > 0){
        // Pollard's rho method with x^2 + 1
        mpz_set_ui(c, 1);
        r = pollard_rho_brent_opt(p, d, c, B, 100);
        
        // fail 
        if (r != 0) break;
        
        // a factor p is found, perform primality test
        if (mpz_probab_prime_p(p, 10) != 0) {
            unsigned long int exp_N = f->exp[f->omega-1];

            // store the prime factor
            mpz_set(f->prime_factors[f->omega-1], p);
            f->exp[f->omega-1] = 0;
            // compute the exponent
            while(mpz_divisible_p(N, p)){
                mpz_divexact(N, N, p);
                f->exp[f->omega-1]++;
            }

            f->exp[f->omega-1] *= exp_N;

            // store the remaining factor n1/p^e only when different than 1
            if (mpz_cmp_ui(N, 1) != 0){
                f->omega++;
                f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->omega)*(sizeof(mpz_t)));
                f->exp = (unsigned long*) realloc(f->exp, (f->omega)*(sizeof(unsigned long)));
                mpz_init_set(f->prime_factors[f->omega-1], N);
                f->exp[f->omega-1] = exp_N;
            }

            // set d to the remaining factor N/p^e
            mpz_set(d, N);

            // primality test for d
            if (mpz_probab_prime_p(d, 10) != 0){
                // remaining factor is prime, hence complete factorization
                r = 1;
                break;  
            }          
        } else mpz_set(d, p); // in case factor p is not prime, apply Pollard's rho for p to find a prime factor of n
    }

    if ((r == 1) || (mpz_cmp_ui(N, 1) == 0)){
        mpz_clears(N, p, d, c, NULL);
        return 1;
    } else {
        // printf("Pollard's rho failed ! Apply Pollard's p-1 method \n");
        mpz_set(d, N);
        while (mpz_cmp_ui(N, 1) > 0){
            r = p_minus_1(p, d, B1, B2);
            if (r != 0) break;
            
            // a factor p is found, perform primality test
            if (mpz_probab_prime_p(p, 10) != 0){
                unsigned long int exp_N = f->exp[f->omega-1];

                // store the prime factor
                mpz_set(f->prime_factors[f->omega-1], p);
                f->exp[f->omega-1] = 0;


                // compute the exponent
                while(mpz_divisible_p(N, p)){
                    mpz_divexact(N, N, p);
                    f->exp[f->omega-1]++;
                }

                f->exp[f->omega-1] *= exp_N;

                // store the remaining factor N/p^e only when different than 1
                if (mpz_cmp_ui(N, 1) != 0){
                    f->omega++;
                    f->prime_factors = (mpz_t*) realloc(f->prime_factors, (f->omega)*(sizeof(mpz_t)));
                    f->exp = (unsigned long*) realloc(f->exp, (f->omega)*(sizeof(unsigned long)));
                    mpz_init_set(f->prime_factors[f->omega-1], N);
                    f->exp[f->omega-1] = exp_N;
                }

                // set d to the remaining factor N/p^e
                mpz_set(d, N); 

                // primality test for d
                if (mpz_probab_prime_p(d, 10) != 0){
                    // remaining factor is prime, hence complete factorization
                    r = 1;
                    break;  
                }          
            } else mpz_set(d, p); // in case factor p is not prime, apply Pollard's rho for p to find a prime factor of n
        }
        
        // remaining factor is 1, complete factorization
        if (mpz_cmp_ui(N, 1) == 0) r = 1;
        // remaining factor is nontrivial and not prime, incomplete factorization
        else if ((r != 1) && (mpz_cmp(N, n) < 0)) r = 0;
        // else r = -1 and N = n indicate that algorithm fails to find a nontrivial factor of n -> failure 
        
        mpz_clears(N, p, d, c, NULL);
        return r;
    }     
}