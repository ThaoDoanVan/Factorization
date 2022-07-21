/*
gcc -fPIC -shared shareLib.c -lgmp -o shareLib.so
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "factorutil.h"

void m_trival_division(char *fileName)
{
    FILE *file;
    file = fopen("result.tmp", "wt");

	mpz_t n, p_max;
	mpz_init(n);
	mpz_init(p_max);

	FILE *fp;
	fp = fopen(fileName, "r");
    int input_count = 0;
	input_count = input_count + gmp_fscanf(fp, "n = %Zu\n", n); 
	input_count = input_count + gmp_fscanf(fp, "p_max = %Zu\n", p_max); 
  	fclose(fp);

	if (input_count == 2)
    {
        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_fprintf(file,"Computing the factorization of %Zd by trial division up to %Zd\n\n", n, p_max);
        clock_t start = clock();

        int k = factorization_division(&f, n, p_max);
        clock_t end = clock();
        if (k == 2) gmp_fprintf(file,"%Zd is a prime number\n\n", n);
        else if (k != -1){
            if (k == 1) fprintf(file,"Fully factored !\n\n");
            else fprintf(file,"Partially factored !\n\n");
            gmp_fprintf(file,"%Zd = ", n);
            for (int i = 0; i < f.omega; i++){
                if (f.exp[i] != 1) gmp_fprintf(file,"%Zd^%lu", f.prime_factors[i], f.exp[i]);
                else gmp_fprintf(file,"%Zd", f.prime_factors[i], f.exp[i]);
                if (i != f.omega - 1) fprintf(file," * ");
            }
            fprintf(file,"\n\n");
        } else gmp_fprintf(file,"Failed to factor %Zd by trial division up to %Zd\n\n", n, p_max);
            
        fprintf(file,"Process finished in %.9f secs \n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        fprintf(file,"\n");
            
        free(f.prime_factors);
        free(f.exp);
    }
	else
    {
        fprintf(file,"Wrong input!");
    }

	mpz_clears(n,p_max, NULL);	
    fclose(file);
}

void m_pollard_rho(char *fileName)
{
    FILE *file;
    file = fopen("result.tmp", "wt");

    mpz_t n;
    mpz_init(n);
    long long int B;

    FILE *fp;
    fp = fopen(fileName, "r");
    int input_count = 0;
    input_count = input_count + gmp_fscanf(fp, "n = %Zu\n", n); 
    input_count = input_count + fscanf(fp, "B = %lli\n", &B); 
    fclose(fp);

    if (input_count == 2)
    {
        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_fprintf(file, "Computing the factorization of %Zd by Pollard's rho using Brent's cycle detection\nTesting power of 2 up to %lli\n\n", n, B);
        clock_t start = clock();

        int k = factorization_pollard_rho(&f, n, B);

        clock_t end = clock();
        if (k == 2) gmp_fprintf(file,"%Zd is a prime number\n\n", n);
        else if (k != -1){
            if (k == 1) fprintf(file,"Fully factored !\n\n");
            else fprintf(file,"Partially factored !\n\n");
            gmp_fprintf(file,"%Zd = ", n);
            for (int i = 0; i < f.omega; i++){
                if (f.exp[i] != 1) gmp_fprintf(file,"%Zd^%d", f.prime_factors[i], f.exp[i]);
                else gmp_fprintf(file,"%Zd", f.prime_factors[i], f.exp[i]);
                if (i != f.omega - 1) fprintf(file," * ");
            }
            fprintf(file,"\n\n");
        } else gmp_fprintf(file,"Failed to factor %Zd by Pollard's rho method up to %lli iterations\n\n", n, B);
        
        fprintf(file,"Process finished in %.9f secs \n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        free(f.prime_factors);
        free(f.exp);
    }
    else
    {
        fprintf(file,"Wrong input!");
    }

    mpz_clear(n);
    fclose(file);
}

void m_pollard_p1(char *fileName)
{
    FILE *file;
    file = fopen("result.tmp", "wt");

    mpz_t n, B1, B2;
    mpz_inits(n, B1, B2);

    FILE *fp;
    fp = fopen(fileName, "r");
    int input_count = 0;
    input_count = input_count + gmp_fscanf(fp, "n = %Zu\n", n); 
    input_count = input_count + gmp_fscanf(fp, "B1 = %Zu\n", B1); 
    input_count = input_count + gmp_fscanf(fp, "B2 = %Zu\n", B2);
    fclose(fp);

    if (input_count == 3)
    {
        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_fprintf(file,"Computing the factorization of %Zd by Pollard's p- 1 method\nStage 1 smoothness bound: %Zd\nStage 2 smoothness bound: %Zd\n\n", n, B1, B2);
        clock_t start = clock();

        int k = factorization_p_minus_1(&f, n, B1, B2);

        clock_t end = clock();
        if (k == 2) gmp_fprintf(file,"%Zd is a prime number\n\n", n);
        else if (k != -1){
            if (k == 1) fprintf(file,"Fully factored !\n\n");
            else fprintf(file,"Partially factored !\n\n");
            gmp_fprintf(file,"%Zd = ", n);
            for (int i = 0; i < f.omega; i++){
                if (f.exp[i] != 1) gmp_fprintf(file,"%Zd^%d", f.prime_factors[i], f.exp[i]);
                else gmp_fprintf(file,"%Zd", f.prime_factors[i], f.exp[i]);
                if (i != f.omega - 1) fprintf(file," * ");
            }
            fprintf(file,"\n\n");
        } else gmp_fprintf(file,"Failed to factor %Zd by Pollard's p-1 method\n\n", n);
        
        fprintf(file,"Process finished in %.9f secs \n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        free(f.prime_factors);
        free(f.exp);
    }
    else
    {
        fprintf(file,"Wrong input!");
    }

    mpz_clears(n, B1, B2, NULL);
    fclose(file);
}

void m_full_factorization(char *fileName)
{
    FILE *file;
    file = fopen("result.tmp", "wt");

    mpz_t n, p_max, B1, B2;
    mpz_inits(n, p_max, B1, B2, NULL);
    long long int B;

    FILE *fp;
    fp = fopen(fileName, "r");
    int input_count = 0;
    input_count = input_count + gmp_fscanf(fp, "n = %Zu\n", n); 
    input_count = input_count + gmp_fscanf(fp, "p_max = %Zu\n", p_max); 
    input_count = input_count + fscanf(fp, "B = %lli\n", &B); 
    input_count = input_count + gmp_fscanf(fp, "B1 = %Zu\n", B1); 
    input_count = input_count + gmp_fscanf(fp, "B2 = %Zu\n", B2);
    fclose(fp);

    if (input_count == 5)
    {
        factor f;
        f.prime_factors = NULL;
        f.exp = NULL;
        f.omega = 0;
        gmp_fprintf(file, "Computing the factorization of %Zd by\nTrial division up to %Zd\nPollard's rho method (with Brent's cycle detection), testing with power of 2 up to %lli iterations\nPollard's p-1 method with stage 1 smoothness bound %Zd and stage 2 smoothness bound %Zd\n\n", n, p_max, B, B1, B2);
        clock_t start = clock();

        int k = factorization(&f, n, p_max, B, B1, B2);

        clock_t end = clock();
        if (k == 2) gmp_fprintf(file,"%Zd is a prime number\n\n", n);
        else if (k != -1){
            if (k == 1) fprintf(file,"Fully factored !\n\n");
            else fprintf(file,"Partially factored !\n\n");
            gmp_fprintf(file,"%Zd = ", n);
            for (int i = 0; i < f.omega; i++){
                if (f.exp[i] != 1) gmp_fprintf(file,"%Zd^%d", f.prime_factors[i], f.exp[i]);
                else gmp_fprintf(file,"%Zd", f.prime_factors[i], f.exp[i]);
                if (i != f.omega - 1) fprintf(file," * ");
            }
            fprintf(file,"\n\n");
        } else gmp_fprintf(file,"Failed\n\n");
            
        fprintf(file,"Process finished in %.9f secs \n\n", (double)(end-start)/CLOCKS_PER_SEC);

        for (int i = 0; i < f.omega; i++) mpz_clear(f.prime_factors[i]);
        free(f.prime_factors);
        free(f.exp);
    }
    else
    {
        fprintf(file,"Wrong input!");
    }

    mpz_clears(n, p_max, B1, B2, NULL);
    fclose(file);
}