#include <stdio.h>
#include <time.h>
#include <gmp.h>

void random(mpz_t , mpz_t, mpz_t);
int isPrime(mpz_t, int);
void generatePrimeNum(mpz_t);
void fastExponent(mpz_t, mpz_t, mpz_t, mpz_t);
void gcd(mpz_t, mpz_t, mpz_t);
void invMulti(mpz_t, mpz_t, mpz_t);

int main(){
	FILE *f = fopen("out.txt", "w");
	srand(time(NULL));
	printf("Key generation:\n");
	// generate 2 random prime numbers p and q
	mpz_t p, q; 
	mpz_init(p); mpz_init(q);
	generatePrimeNum(p);
	generatePrimeNum(q);
	gmp_printf("p = %Zd\n\n", p);
	mpz_out_str(f, 10, p);
	fputc('\n', f);
	gmp_printf("q = %Zd\n\n", q);
	mpz_out_str(f, 10, q);
	fputc('\n', f);
	// calculate n = p * q
	mpz_t n; mpz_init(n);
	mpz_mul(n, p, q);
	gmp_printf("n = %Zd\n\n", n);
	// calculate phi(n)
	mpz_t phi; mpz_init(phi);
	mpz_sub_ui(p, p, 1);
	mpz_sub_ui(q, q, 1);
	mpz_mul(phi, p, q);
	mpz_out_str(f, 10, phi);
	fputc('\n', f);
	mpz_clear(p); mpz_clear(q);
	// generate publickey e
	mpz_t e, Gcd, t; 
	mpz_init(e); mpz_init(Gcd); mpz_init(t);
	mpz_set_ui(t, 2);
	do{
		random(t, phi, e);
		gcd(e, phi, Gcd);
	}
	while(mpz_cmp_ui(Gcd, 1) != 0);
	gmp_printf("e = %Zd\n\n", e);
	mpz_out_str(f, 10, e);
	fputc('\n', f);
	// generate privatekey d
	mpz_t d; mpz_init(d);
	invMulti(e,phi,d);
	mpz_clear(phi);
	gmp_printf("d = %Zd\n\n", d);
	// plaintext
	mpz_t m; mpz_init(m);
	mpz_set_str(m,"18111997", 10);
	gmp_printf("Plaintext: %Zd\n\n", m);
	// ciphertext
	mpz_t c; mpz_init(c);
	fastExponent(m, e, n, c);
	gmp_printf("Ciphertext: %Zd\n\n", c);
	mpz_out_str(f, 10, c);
	// recover plaintext
	mpz_t r; mpz_init(r);
	fastExponent(c, d, n, r);
	gmp_printf("Plaintext (recover): %Zd", r);
	// Clear temporary variables
	mpz_clear(n); mpz_clear(e); mpz_clear(d);
	mpz_clear(m); mpz_clear(c); mpz_clear(r);
	fclose(f);
	return 0;
}

void random(mpz_t start, mpz_t end, mpz_t res){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, rand() * rand());
	mpz_t t1, t2; mpz_init(t1); mpz_init(t2);
	mpz_set(t1, start); mpz_set(t2, end);
	mpz_sub(t2, t2, t1);
	mpz_urandomm(res, state, t2);
	mpz_add(res, res, t1);
	gmp_randclear(state);
}

int isPrime(mpz_t n, int test){
	if(mpz_cmp_ui(n, 2) < 0)
		return 0;
	int ok = 1;
	int i;
	mpz_t v, t, _n, p;
	mpz_init(v); mpz_init(t); mpz_init(_n); mpz_init(p);
	mpz_set_ui(t, 2);
	mpz_sub_ui(_n, n, 1);		// _n = n - 1
	for(i = 0; i < test; i++){
		random(t, _n, v);
		fastExponent(v, _n, n, p);
		if(mpz_cmp_ui(p, 1) != 0){
			ok = 0;
			break;
		}
	}
	// Clear temporary variables
	mpz_clear(v); mpz_clear(t); mpz_clear(_n); mpz_clear(p);
	// Get result
	return ok;
}

void generatePrimeNum(mpz_t res){
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, rand() * rand());
	while(!isPrime(res, 5000)){
		mpz_urandomb(res, state, 512);
	}
	gmp_randclear(state);
}

void fastExponent(mpz_t a, mpz_t m, mpz_t n, mpz_t res){
	mpz_t b; mpz_init(b);
	mpz_set_ui(b, 1);
	while(mpz_cmp(b, m) <= 0){
		mpz_mul_ui(b, b, 2);
	}
	mpz_div_ui(b, b, 2);
	mpz_t p; mpz_init(p);
	mpz_set_ui(p, 1);			// init p = 1
	while(mpz_cmp_ui(b, 0) > 0){
		// p = (p * p) mod n
		mpz_mul(p, p, p);		
		mpz_mod(p, p, n);
		mpz_t t; mpz_init(t);
		mpz_and(t, b, m);
		if(mpz_cmp_ui(t, 0) != 0){		// bit == 1 ?
			// p = (p * a) mod n
			mpz_mul(p, p, a);	
			mpz_mod(p, p, n);
		}
		mpz_div_ui(b, b, 2);
	}
	// Get result
	mpz_set(res, p);
	// Clear temporary variables
	mpz_clear(p); mpz_clear(b);
}

void gcd(mpz_t a, mpz_t b, mpz_t res){		// Euclidean
	mpz_t r1, r2, q, r;
	mpz_init(r1); mpz_init(r2);
	mpz_init(q); mpz_init(r);
	mpz_set(r1, a); mpz_set(r2, b);
	while(mpz_cmp_ui(r2, 0) > 0){
		mpz_div(q, r1, r2);
		mpz_mod(r, r1, r2);
		mpz_set(r1, r2);
		mpz_set(r2, r);
	}
	mpz_set(res, r1);
	// Clear temporary variables
	mpz_clear(r1); mpz_clear(r2); mpz_clear(q); mpz_clear(r);
}

void invMulti(mpz_t a, mpz_t n, mpz_t res){	// extendEuclidean
	mpz_t q, r1, r2, r, t1, t2, t;
	mpz_init(q); mpz_init(r); mpz_init(r1); mpz_init(r2);
	mpz_init(t1); mpz_init(t2); mpz_init(t);
	mpz_set(r1, n); mpz_set(r2, a);
	mpz_set_ui(t1, 0); mpz_set_ui(t2, 1);
	while(mpz_cmp_ui(r2, 0) > 0){
		mpz_div(q, r1, r2);     // q = r1 div r2
		mpz_mod(r, r1, r2);		// r = r1 mod r2
		mpz_set(r1, r2); mpz_set(r2, r);
		mpz_mul(q, q, t2);		
		mpz_sub(t, t1, q);
		mpz_set(t1, t2);
		mpz_set(t2, t);
	}
	if(mpz_cmp_ui(t1, 0) < 0)
		mpz_add(t1, t1, n);
	mpz_set(res, t1);
	// Clear temporary variables
	mpz_clear(q); mpz_clear(r); mpz_clear(r1); mpz_clear(r2);
	mpz_clear(t); mpz_clear(t1); mpz_clear(t2);
}
