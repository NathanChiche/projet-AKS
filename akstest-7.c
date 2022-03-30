#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h> 
//la librairie gmp est comprise dans mpfr





// étape 1 de AKS : renvoie 0 si n = a^b avec b>1, 1 sinon

int e1(mpz_t n)
{
	if (mpz_perfect_power_p(n)==0)
	{
		return 1;
	}
	else
	{
		return 0;
	} 
}

// étape 2&3 : recherche du r

void e23(mpfr_t *l, mpz_t *r, mpz_t n)
{	
	mpfr_t N;
	mpfr_init2(N, 200);
	/* On crée une variable N de type mpfr qui contient
	 * la valeur de n afin de calculer log(n).
	* Le 200 signifie que la précision est à 200bits près.*/
	
	mpfr_set_z(N, n, MPFR_RNDN);
	
	mpfr_log2(*l, N, MPFR_RNDN);// l =log(n) en base 2
	
	
	/* on affiche l sur le canal stdout
	 * ie le canal de sortie ie le terminal*/
	gmp_printf("log %Zd =\n", n);
	mpfr_out_str (stdout, 10, 0, *l, MPFR_RNDN);
	printf("\n\n");
	
	
	mpfr_mul(*l, *l, *l, MPFR_RNDN);
	
	// on affiche l avec sa nouvelle valeur
	gmp_printf("log²%Zd =\n", n);
	mpfr_out_str (stdout, 10, 0, *l, MPFR_RNDN);
	printf("\n\n");
	
	mpz_t d;// pgcd(r,n)
	mpz_init(d);
	
	mpz_t x;
	mpz_init(x);
	/* x prendra la valeur n afin de calculer 
	 * les puissances successives de n mod r.*/
	
	
	unsigned long int k; // compteur dans une boucle.
	int b = 0; 
	/* break le while lorsque b != 0
	 * Si b = 1, alors le r recherché est trouvé
	 * Si b = 2, alors n est composé */
	
	
	/* Cette boucle while se termine forcement car
	 * il existe r dans [2, log⁵n] qui satisfait 
	 * le if de la ligne 104*/
	while (b == 0)
	{
		mpz_add_ui(*r, *r, 1); 
		
		mpz_gcd(d, *r, n); // d prends la valeur pgcd(r, n)
		
		/*étape 3 d'AKS : si 1 < d < n, alors d | n
		 * donc n est composé*/
		if(mpz_cmp_ui(d, 1) != 0 && mpz_cmp(d, n) != 0) 
		{
			gmp_printf("COMPOSITE : %Zd\n", d);
			b = 2;
		}
		
		/* pour que n admette un ordre mod r,
		 * il faut s'assurer que d = 1*/
		if(mpz_cmp_ui(d, 1) == 0)
		{
			mpz_set(x, n); 
			// x est réinitialisé à n avant chaque boucle
			k = 1; 
			
			while (mpz_cmp_ui(x, 1) != 0 && mpfr_cmp_ui(*l, k) >= 0)
			{
				// calcul de n^k mod r
				mpz_mul(x, x, n);
				mpz_mod(x, x, *r);
				k++;
			}
			/* si l < k, alors n^a != 1 pour tout a<l
			 * donc ord(n) > l dans Z/rZ, donc r correspond*/
			if(mpfr_cmp_ui(*l, k) < 0) 
			{
				b = 1;
			}		
			
		}
	}

	if(b == 2)
	{
		mpz_set_ui(*r, 0); 	// r >= 2 lorsque b != 2
	}
	gmp_printf("r = %Zd\n\n", *r);
	
	mpfr_clear(N);

	mpz_clear(d);
	mpz_clear(x);
}


/* Calcul de phi(r) qui nous sera utile pour trouver la
 * borne l de l'étape 5.
 * il s'agit simplement d'identifier les élements inversibles
 * de Z/rZ en verifiant "pgcd(i, r)=1?" pour 0 < i < r */

void phi(mpz_t *k, mpz_t r)
{
	mpz_t d;	// pgcd(i,r)
	mpz_init(d);
	mpz_t i;
	mpz_init_set_ui(i, 2); 
	// on sait que 0 n'est pas inversible et 1 l'est.
	 
	
	while(mpz_cmp(r,i) > 0)
	{
		mpz_gcd(d, i, r);
		if(mpz_cmp_ui(d, 1) == 0)
		//on teste pgcd(i,r) pour tout 2<i<r
		{
			mpz_add_ui(*k, *k, 1);	
		}
		mpz_add_ui(i, i, 1);
	}
	
	mpz_clear(d);
	mpz_clear(i);
}


// étape 4 de l'algorithme
int e4(mpz_t n, mpz_t r)
{
	if(mpz_cmp(r, n) >= 0)
	{
		return 1;
	}
	return 0;
}




typedef struct 
{
    mpz_t* coef;
  	unsigned int deg;
}Polynome;



/* Nous allons à présent nous attarder sur l'étape 5
 * qui demande plus de travail que les précedentes.
 *  
 * 
 * On considère que le degré est entier car le degré d'un 
 * polynome est toujours inférieur ou égal à r-1 dans 
 * A = (Z/nZ)/((X^r)-1), avec r assez petit pour être 
 * transformé en unsigned int.
 * 
 * Dans l'étape 2, r est un mpz_t. R est une variable
 * de type unsigned int qui prends la valeur de r */



// Constructeur

void poly_initialise(Polynome** p, unsigned int R)
	{
  	*p = (Polynome*)malloc(sizeof(Polynome));
  	(*p)->coef = (mpz_t*)malloc((R) * sizeof(mpz_t));
  	(*p)->deg = 0;
  	unsigned int i;
  	for (int i = 0; i < R; i++)
  	{
  		mpz_init_set_ui((*p)->coef[i], 0);
  	}
}


// Destructeur

void detruire_polynome(Polynome **p, unsigned int R)
{
	unsigned int i;
  	for (i = 0; i < R; i++) {
   		mpz_clear((*p)->coef[i]);
  	}
  	free((*p)->coef); 
  	free(*p);			
}


/* On reinitialise le polynome p au polynome nul
 * Cette fonction correspond à la fonction précente,
 * sans réallouer la mémoire*/

void poly_zero(Polynome **p)
{
	unsigned int i;
	for(i = 0; i <= (*p)->deg; i++)
	{
		mpz_set_ui((*p)->coef[i], 0);
	}
	(*p)->deg = 0;
}



// teste si p_1 = p_2 coef par coef

int poly_egalite(Polynome *p_1, Polynome *p_2)
{
	/* On teste l'égalité des degrés afin de n'avoir qu'à
	 * verifier les coefficients inferieurs au degré*/	
	if(p_1->deg != p_2->deg) 				
	{
		return 0;
	}
	unsigned int i;
	for(i = 0; i <= p_1->deg; i++)
	{
		if(mpz_cmp(p_1->coef[i],p_2->coef[i]) != 0)
		{
			return 0;
		}
	}
	return 1;
}

// Affiche les coeffs non nuls du polynome p.

void poly_affiche(Polynome *p)
{
	int i = p->deg;
	gmp_printf("%Zdx^%d", p->coef[i], i);
	i--;
	while(i >= 0)
	{
		// N'affiche que les coefficients <= deg
		if(mpz_cmp_ui(p->coef[i], 0) > 0)
		{
			gmp_printf(" + %Zdx^%d", p->coef[i], i);
		}
		i--;
	}
	printf("\n");
}


/* On copie les coefficients de p sur res.
 * On considère que res est un polynôme déjà initialisé. */

void poly_copy(Polynome **res, Polynome *p)
{
	unsigned int i;
	if((*res)->deg > p->deg)
	{
		for(i = p->deg+1; i <= (*res)->deg; i++)
		{
			mpz_set_ui((*res)->coef[i], 0);
		}
	}
	(*res)->deg = p->deg;
	
	for(i = 0; i <= (*res)->deg; i++)
	{
		mpz_set((*res)->coef[i], p->coef[i]);
	}
}


/* Corrige le degré du polynôme si son degré est erroné.
 * Celà peut arriver lors de reductions modulo n*/

void poly_degreduit(Polynome **p)
{
	unsigned int i = (*p)->deg;
	// Tant que i>0 & coef[i] = 0
	while (i > 0 && mpz_cmp_ui((*p)->coef[i], 0) == 0)
	{
		i--;
	}
	(*p)->deg = i;
}

// Ajoute le terme X^k de coefficient a au polynome res.

void poly_addterme(Polynome **res, mpz_t a, unsigned int k)
{
	if ((*res)->deg < k)
	{
		(*res)->deg = k;
	}
	mpz_add((*res)->coef[k], (*res)->coef[k], a);
}


/* On additionne les polynomes p_1 et p_2 et on stocke 
 * le résultat dans res. 
 * Les coefficients sont également reduits modulo n */

void addition_poly(Polynome **res, Polynome *p_1, Polynome *p_2, mpz_t n)
{
	unsigned int degmax;
	unsigned int i;
	
	if (p_1->deg >= p_2->deg)
	{
		degmax = p_1->deg;
	}
	else
	{
		degmax = p_2->deg;
	}
	
	if((*res)->deg < degmax)
	{
		(*res)->deg = degmax;
	}
	//deg res = max(deg p_1,deg p_2)
	
	for(i = 0; i <= degmax; i++)
	{
		mpz_add((*res)->coef[i], p_1->coef[i], p_2->coef[i]);
		mpz_mod((*res)->coef[i], (*res)->coef[i], n);
	}
}

/* On multiplie res par un scalaire(= poly de deg0) a
 * les coefficients sont encore réduits modulo n.*/

void multiplication_scalaire(Polynome **res, mpz_t a, mpz_t n)
{
	unsigned int i;
	
	for(i = 0; i <= (*res)->deg; i++)
	{
		mpz_mul((*res)->coef[i], (*res)->coef[i], a);
		mpz_mod((*res)->coef[i], (*res)->coef[i], n);
	}
}


/* Ici, nous effectuons une multiplication de res par X^k.
 * Puisque X^r = 1 dans A, il suffit d'effectuer une rotation
 * des coefficients : le i-ème coefficient prend la valeur
 * du (i-k)-ème coefficient.*/

void multiplication_X_g(Polynome **res, unsigned int k, unsigned int R)
{
	/* Si k=0, il n'y a rien à faire.
	 * Si k >= r, alors le resultat est inchangé si l'on prend
	 * le reste de la division euclidienne de k par r.
	 * Cependant, cette précaution est inutile dans notre code
	 * puisque l'on effectue au plus une multiplication par X^{r-1}*/
	 
	if(k != 0)
	{
		mpz_t x;	// stocke la valeur du coef precedent
		mpz_init(x);
		mpz_t y;	// stocke la valeur du coef present
		mpz_init(y);
		mpz_t pgcd;
		mpz_init(pgcd);
		mpz_t r;
		mpz_init_set_ui(r, R); 
		// On a besoin d'un mpz pour calculer le pgcd
	
	
		unsigned int i;
		unsigned int j;
		unsigned int d = R-1;	// degré du polynôme
	
		mpz_gcd_ui(pgcd, r, k);
		
		// on effectue d rotations.
		for(i = 0; mpz_cmp_ui(pgcd, i) > 0; i++)
		{
			mpz_set(x, (*res)->coef[i]);
			j = i+k;
		
			while(j != i)	// fin de la (i+1)ème rotation
			{
				mpz_set(y, (*res)->coef[j]);
				// On a x = coef[j-k].
				mpz_set((*res)->coef[j], x);
				mpz_set(x, y);
				j = j+k;
				if(j >= R)
				{
					j = j-R;// X^j = X^(j-R)
				}
			}
			mpz_set((*res)->coef[i], x);
		}
		
		/* Si le degré >= R-k, le coefficient dominant
		 * le coefficient dominant change. 
		 * Il faut donc retrouver le bon degré. */
		if ((*res)->deg >= R-k)
		{
			while(mpz_cmp_ui((*res)->coef[d], 0) == 0 && d != 0)
			{
				d--;
			}
			(*res)->deg = d;
		}
		else
		{
			(*res)->deg = (*res)->deg+k;
		}
	
		mpz_clear(x);
		mpz_clear(y);
		mpz_clear(pgcd);
		mpz_clear(r);
	}
}


/* La fonction square du square and multiply
 * Si f = a_d*X^d + ..... + a_1*X + a_0, alors :
 * f² = a_d*(f*X^d) + ....+ a_1*(f*X) + a_0*f
 * f*X^i peut être calculé grâce à multiplication_X_g
 * ai*f peut être calculé grâce à multiplication_scalaire*/

void poly_square(Polynome **res, mpz_t n, unsigned int R)
{
	unsigned int i;
	Polynome *p_1, *p_2;
	
	poly_initialise(&p_1, R);
	poly_initialise(&p_2, R);
	poly_copy(&p_1, *res);	// p_1 <- res = f
	poly_zero(res);			// res <- 0
	
	for(i = 0; i <= p_1->deg; i++)
	{
		if(mpz_cmp_ui(p_1->coef[i], 0) > 0)
		{
			// p_2 <- p_1 = f
			poly_copy(&p_2, p_1);	
			multiplication_X_g(&p_2, i, R);
			multiplication_scalaire(&p_2, p_1->coef[i], n);
			// res <- res + a_i(f*X^i)
			addition_poly(res, *res, p_2, n);
		}
	}
	
	detruire_polynome(&p_1, R);
	detruire_polynome(&p_2, R);
}


/* Voici la fonction multiply du square and multiply.
 * Dans l'algorithme AKS, nous allons exclusivement 
 * multiplier par un polynome de la forme X+a dans 
 * le square and multiply.
 * Nous avons donc crée une fonction se limitant à cet usage
 * en additionnant une rotation d'ordre 1 et une 
 * multiplication_scalaire */

void poly_multiply(Polynome **res, mpz_t a, mpz_t n, unsigned int R)
{
	Polynome *p_1, *p_2;
	poly_initialise(&p_1, R);
	poly_initialise(&p_2, R);
	poly_copy(&p_1, *res);	// p_1 <- f
	poly_copy(&p_2, *res);	// p_2 <- f
	poly_zero(res);
	
	// La multiplication par 0 est inutile
	if(mpz_cmp_ui(a, 0) > 0)
	{
		multiplication_scalaire(&p_2, a, n);
		addition_poly(res, *res, p_2, n); // res <- a*f
	}
	multiplication_X_g(&p_1, 1, R);
	addition_poly(res, *res, p_1, n); 
	// res <- a*f + X*f = (X+a)*f
	
	detruire_polynome(&p_1, R);
	detruire_polynome(&p_2, R);
}

// La fonction square and multiply classique

void exponentiation_rapide(Polynome **res, mpz_t a, mpz_t n, unsigned int R)
{
	int i;
	mpz_t x;
	mpz_init_set_ui(x, 1);
	
	poly_zero(res);
	poly_addterme(res, x, 0);
	for(i = mpz_sizeinbase(n, 2) ; i>=0; i--)
	{
		poly_square(res, n, R); // res <- res²
		if(mpz_tstbit(n,i))
		{
			poly_multiply(res, a, n, R); // res <-
		}
	}
	
	/* De nombreux calculs modulaires ont été effectués
	 * après l'exponentiation. Nous réajustons donc le
	 * degré de res en appelant degreduit:
	 * Nous avons besoin d'avoir le bon degré pour 
	 * que poly_egalite fonctionne correctement.*/
	poly_degreduit(res);
	mpz_clear(x);
}


/* L'étape 5 de l'algorithme
 * Nous vérifions que (X+a)^n = (X+a) dans (Z/nZ)/((X^r)-1)
 * pour tout a < l
 * S'il existe a tel que (X+a)^n != X+a, l'algorithme s'arrete
 * et renvoie 0. Sinon, l'algorithme renvoie 1*/

int e5(Polynome **res, mpz_t n, unsigned int R, unsigned int l)
{
	Polynome *p_1;	// le polynôme X^n+a
	poly_initialise(&p_1, R);
	
	unsigned int m;
	
	mpz_t a;
	mpz_init_set_ui(a, 1);
	
	mpz_t x;
	mpz_init_set_ui(x, 1);
	
	mpz_t N;
	mpz_init(N);
	mpz_mod_ui(N, n, R);
	m = mpz_get_ui(N);	
	/* p_1 ne passe pas par l'exponentiation rapide,
	 * donc il n'y a pas de reduction du degré si N > r
	 * (ce qui est forcement le cas puisque nous sommes
	 * à l'étape 5. Il faut donc réduire le degré "à la main"
	 * avant afin que la fonction poly_affiche fonctionne
	 * correctement.*/
	
	while(mpz_cmp_ui(a, l) <= 0)
	{
		poly_zero(&p_1);
		exponentiation_rapide(res, a, n, R);// res = (X+a)^n
		poly_addterme(&p_1, x, m);
		poly_addterme(&p_1, a, 0);	// p_1 = X^m+a (= X^n+a)
		
		// teste si res = p_1 ie (X+a)^n = X^n+a
		if(poly_egalite(*res, p_1) == 0)
		{
			gmp_printf("a = %Zd\n", a);
			detruire_polynome(&p_1, R);
			mpz_clear(a);
			mpz_clear(x);
			mpz_clear(N);
			return 0;
		}
		mpz_add_ui(a, a, 1);
	}
	detruire_polynome(&p_1, R);
	mpz_clear(a);
	mpz_clear(x);
	mpz_clear(N);
	
	return 1;
}

/* l'algorithme AKS : il prend en entrée un mpz_t n et renvoie
 * 1 si n est premier, et 0 si n est composé.*/

int ISPRIME(mpz_t n)
{
	if(e1(n)==0)
	{
		gmp_printf("%Zd est composé : ETAPE 1\n\n", n);
		return 0;
	}
	
	mpz_t r;
	mpz_init_set_ui(r, 2);
	
	mpfr_t l_0;
	mpfr_init2(l_0, 200);	// 200bits de precision pour l_0
	
	e23(&l_0, &r, n); 
	
	/* A la sortie de e23, il y 2 possibilités:
	 * Si r = 0, alors la boucle s'est terminée car b=2,
	 * donc n est composé.
	 * Sinon, la boucle s'arrete forcement pour b=1, donc
	 * r prends la valeur correcte.
	 * De plus, l_0 resort avec la valeur log²n.*/
	
	
	if(mpz_cmp_ui(r, 0) == 0)
	{
		gmp_printf("%Zd est composé : ETAPE 3\n\n", n);
		return 0;
	}
	
	mpz_t k;
	mpz_init_set_ui(k, 1);
	/* k est initialisé à 1 car l'on ne teste pas pgcd(1,r)
	 * dans notre fonction phi*/
	
	phi(&k, r);	// k prends la valeur phi(r).
	
	if(e4(n, r))// si n <= r
	{
		gmp_printf("%Zd est premier : ETAPE 4\n\n", n);
		return 1;
	}
	
	unsigned int l;
	
	// l_0 <- l_0*phi(r) = log²(n)phi(r)
	mpfr_mul_z(l_0, l_0, k, MPFR_RNDN);	
	// l_0 <- sqrt(l_0) = sqrt(phi(r))*log(n)
	mpfr_sqrt(l_0, l_0, MPFR_RNDN);
	// l = int(l_0).
	l = mpfr_get_ui(l_0, MPFR_RNDD);
	
	printf("l = %u\n\n", l);
	
	
	unsigned int R;
	R = mpz_get_ui(r);
	
	Polynome *res;
	poly_initialise(&res, R);
	
	int p;
	
	p = e5(&res, n, R, l); //teste l'étape 5
	
	detruire_polynome(&res, R);
	mpfr_clear(l_0);
	mpz_clear(r);
	mpz_clear(k);
	
	/* On ne peut pas renvoyer e5 puisqu'il faut
	 * désallouer la mémoire avant.*/
	if(p == 0)
	{
		gmp_printf("%Zd est composé : ETAPE 5\n\n", n);
		return 0;
	}
	else
	{
		gmp_printf("%Zd est premier : ETAPE 6\n\n", n);
		return 1;
	}
}



int main(int argc, char **argv)
{
	mpz_t n;
	
	// prends en arg le terminal
	mpz_init_set_str(n, argv[1], 10);
	
	ISPRIME(n);
	
	mpz_clear(n);
	
	return 0;
}
