/*
find a prime divisor p of a composite number n by trial division by prime numbers 
from p_min to p_max
*/
int successive_division(mpz_t p, mpz_t n, mpz_t p_min, mpz_t p_max){

    mpz_nextprime(p, p_min);
    int r = -1;

    while (mpz_cmp(p, p_max) <= 0){
        if (mpz_divisible_p(n, p)) {
            r = 0;
            break;
        }
        mpz_nextprime(p, p);
    }

    return r;
}

/*
finding a divisor d of a composite number n using Pollard's rho method
given n, B, c; a divisor is found by finding a collision of the function x^2 + c mod n
to find a collision, Floyd's cycle detection method is used
B is the maximum number of iterations
*/

int pollard_rho(mpz_t d, mpz_t n, mpz_t c, long long int B){
    mpz_t x, y, t;
    mpz_inits(x, y, t, NULL);

    /*
    // this is to test if randomize the initial x0 results in better running time on average
    gmp_randstate_t state;
    unsigned long int seed = time(NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
    mpz_urandomm(x, state, n);
    gmp_randclear(state);
    mpz_set(y, x);
    */

    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);

    long long int i = 0;
    while( (mpz_cmp_ui(d, 1) == 0) || (mpz_cmp(d, n) == 0)){
        // maximum number of iterations reached
        if (i > B) break;

        // compute f(x) = (x^2 + c) mod n
        mpz_mul(x, x, x);
        mpz_add(x, x, c);
        mpz_mod(x, x, n);
        i++;

        // compute f(f(y))
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);
        
        // compute gcd(|x-y|, n)
        mpz_sub(t, x, y);
        mpz_abs(t, t);
        mpz_gcd(d, t, n);
        
    }
    mpz_clears(x, y, t, NULL);
    if ((mpz_cmp_ui(d, 1) == 0) || (mpz_cmp(d, n) == 0)) return -1;
    else return 0;
}

/*
an optimized version of Pollard's rho method
*/
int pollard_rho_opt(mpz_t d, mpz_t n, mpz_t c, long long int B, int m){
    mpz_t x, y, s, t;
    mpz_inits(x, y, s, t, NULL);
    
    /*
    gmp_randstate_t state;
    unsigned long int seed = time(NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
    mpz_urandomm(x, state, n);
    gmp_randclear(state);
    mpz_set(y, x);
    */

    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);
    mpz_set_ui(s, 1);

    long long int i = 0;
    int g = 0;
    while( (mpz_cmp_ui(d, 1) == 0) || (mpz_cmp(d, n) == 0)){
        if (i > B) break;

        // compute f(x) = (x^2 + c) mod n
        mpz_mul(x, x, x);
        mpz_add(x, x, c);
        mpz_mod(x, x, n);
        i++;

        // compute f(f(y))
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mul(y, y, y);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);
        

        mpz_sub(t, x, y);
        mpz_abs(t, t);
        mpz_mul(s, s, t); 
        mpz_mod(s, s, n);
        g++;

        // compute gcd after iterated for m times
        // s = |x_1 - y_1| * |x_2 - y_2| * ... * |x_100 - y_100|
        if (g == m){
            mpz_gcd(d, s, n);
            g = 0;
        }
    }
    mpz_clears(x, y, s, t, NULL);
    if ((mpz_cmp_ui(d, 1) == 0) || (mpz_cmp(d, n) == 0)) return -1;
    else return 0;  
}

/*
Pollard's rho method using Brent's cycle detection
*/
int pollard_rho_brent(mpz_t d, mpz_t n, mpz_t c, long long int B){
    mpz_t x, y, t;
    mpz_inits(x, y, t, NULL);
        
    // initialize y and p
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);

    // r is a power of 2
    long long int r = 1;
    
    while ((mpz_cmp_ui(d, 1) == 0) || (mpz_cmp(d, n) == 0)){
        
        // x to the current position of y, namely x_i
        mpz_set(x, y);

        long long int j = 0;
        long long int k = 0;

        // y to x_(i+r)
        while (j < r){
            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);
            j++;
        }

        // test gcd(|x-y|, n); where y ranges from x_(i+r) to x_(i+2r)
        while ((k < r) && (mpz_cmp_ui(d, 1) == 0)){
            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);
            mpz_sub(t, x, y);
            mpz_abs(t, t);
            mpz_mod(t, t, n);                  
            mpz_gcd(d, t, n);
            k++;
        }

        // going to the next power of 2
        r = r*2;
        if (r >= B) break;
    }

    mpz_clears(x, y, t, NULL);
    if ((mpz_cmp(d, n) == 0) || (mpz_cmp_ui(d, 1) == 0)) return -1;
    else return 0;    
}

/*
optimized version of Pollard's rho method using Brent's cycle detection
*/
int pollard_rho_brent_opt(mpz_t d, mpz_t n, mpz_t c, long long int B, int m){
    mpz_t x, y, ys, s, t;
    mpz_inits(x, y, ys, s, t, NULL);
        
    // initialize y and p
    mpz_set_ui(y, 2);
    mpz_set_ui(d, 1);
    mpz_set_ui(s, 1);

    // r is a power of 2
    long long int r = 1;
    

    while ((mpz_cmp(d, n) == 0) || (mpz_cmp_ui(d, 1) == 0)) {
        // x to the current position of y, namely x_i
        mpz_set(x, y);
        long long int j = 0;
        long long int k = 0;

        // y to x_(i+r)
        while (j < r){
            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);
            j++;
        }

        // test gcd(|x-y|, n); where y ranges from x_(i+r) to x_(i+2r)
        while ((k < r) && (mpz_cmp_ui(d, 1) == 0)){
            mpz_set(ys, y);
            // compute the gcd after every m iterations
            for(int i = 0; ( (i < m) || (i < r-k) ); i++) {
                mpz_mul(y, y, y);
                mpz_add(y, y, c);
                mpz_mod(y, y, n);
                mpz_sub(t, x, y);
                mpz_abs(t, t);
                mpz_mul(s, s, t);
                mpz_mod(s, s, n);           
            }
            // s = x_k * x_(k+1) * ... * x_(k + m -1)
            mpz_gcd(d, s, n);
            k = k + m;
        }

        // go to the next power of 2
        r = r*2;
        if (r > B) break;   
    }

    /*
    // backtrack to the state before the last gcd computation
    // the original version of Brent backtracks to compute the factor in case d = n
    // for optimization purpose, we do not perform the backtracking step

    long long int l = 0;
    if (mpz_cmp(d, n) == 0) printf("Backtracking !\n");
    while (mpz_cmp(d, n) == 0){
        mpz_mul(ys, ys, ys);
        mpz_add(ys, ys, c);
        mpz_mod(ys, ys, n);
        mpz_sub(t, x, ys);
        mpz_abs(t, t);
        mpz_gcd(d, t, n);
        l++;
        if (l > B) break;
    }

    */

    mpz_clears(x, y, ys, s, t, NULL);
    if ((mpz_cmp(d, n) == 0) || (mpz_cmp_ui(d, 1) == 0)) return -1;
    else return 0; 
}


int p_minus_1(mpz_t d, mpz_t n, mpz_t B1, mpz_t B2){
    mpz_t a, p, q, t;
    mpz_init_set_ui(a, 1);

    gmp_randstate_t state;
    unsigned long int seed = time(NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
    mpz_urandomm(a, state, n);

    // randomize the base a
    while (mpz_cmp_ui(a, 1) <= 0) mpz_urandomm(a, state, n);
    
    // return gcd(a, n) if nontrivial
    mpz_init(t);
    mpz_gcd(t, a, n);
    if (mpz_cmp_ui(t, 1) > 0){
        mpz_set(d, t);
        mpz_clears(a, t, NULL);
        return 0;
    }

    mpz_inits(p, q, NULL);
    mpz_set_ui(p, 2);
    mpz_set(d, n);

    int r = -1;
    
    while(mpz_cmp(d, n) == 0){
        while (mpz_cmp(p, B1) <= 0){

            mpz_set_ui(q, 1);
            while(mpz_cmp(q, B1) <= 0){
                mpz_powm(a, a, p, n); // a <- a^p mod n
                mpz_mul(q, q, p); // q <- q*p
            }

            mpz_sub_ui(t, a, 1);
            mpz_gcd(d, t, n);

            // a divisor > 1 found
            if (mpz_cmp_ui(d, 1) > 0){
                r = 0;
                break;
            }

            // go to the next prime
            mpz_nextprime(p, p);
        }

        // all prime powers are tested, no divisor found
        if (mpz_cmp_ui(d, 1) == 0) r = -1;

        // found divisor d = n, rerandomize the base a and restart the calculation
        else if (mpz_cmp(d, n) == 0){
            while (mpz_cmp_ui(a, 1) <= 0) mpz_urandomm(a, state, n);
            mpz_gcd(t, a, n);

            if (mpz_cmp_ui(t, 1) != 0){
                mpz_set(d, t);
                r = 0;
                break;
            }
            mpz_set_ui(p, 2);
            mpz_set_ui(q, 1);
        }
    }

    gmp_randclear(state);

    // stage 2 of the algorithm, only started when B2 > B1 and no factor found after stage 1
    if ((mpz_cmp(B2, B1) > 0) && (r != 0)){
        // precompute the prime gap table
        // gap between primes less than 10^15 is less than 1000

        mpz_t p_next, g, gap[1000];
        mpz_inits(p_next, g, NULL);

        // at the end of stage 1, a = base^(product of prime powers <= B1) mod n
        // compute a^(2k) for 0 < 2k < 1000 and store in gap[2k-1] 
        for(int i = 1; i < 1000; i = i+2){
            mpz_init(gap[i]);
            mpz_powm_ui(gap[i], a, i+1, n);
        }

        // p is the first prime > B1, then first we must compute a^p mod n
        mpz_powm(a, a, p, n);
        

        // compute a = a^p * a^(next prime to p - p) for all primes p <= B2 
        while(mpz_cmp(p, B2) <= 0){
            mpz_nextprime(p_next, p);
            mpz_sub(g, p_next, p); // compute gap = next prime p_next - current prime p
            int i = mpz_get_ui(g);
            mpz_mul(a, a, gap[i-1]); // compute a^p * a^gap
            mpz_mod(a, a, n);
            
            // compute gcd(a - 1, n)
            mpz_sub_ui(t, a, 1);
            mpz_gcd(d, t, n);
            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0){
                r = 0;
                break;
            }

            // go to the next prime
            mpz_set(p, p_next);
        }

        // clear data;
        mpz_clears(p_next, g, NULL);
        for(int i = 1; i < 1000; i = i+2) mpz_clear(gap[i]);
    }
    
    mpz_clears(a, p, q, t, NULL);
    return r;
}