#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdbool.h>
#include <chrono>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <iterator>

#define MATRIX_INT std::vector<std::vector<int>>
#define MATRIX_FLOAT std::vector<std::vector<float>>

#define k_d 0.0
//#define k_on 0.56 //608.478 //5*1000000/(0.62*(4/3)*M_PI*(pow(7.5,3)))
//#define k_off 1.0
//#define k_on 0.56
#define k_off 10.0
#define k_on 0.04

#define NBrow 10


//VALEURS POUR RAN2
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)



# define L 100.0 //longueur de virus
# define diameter 60.0 //diamÃ¨tre du virus
# define D 10.0 //distance maximale pour les rÃ©actions
// # define MU 0.65 //viscositÃ© de l'eau 
// # define KbT 1.38*pow(10,-23)*308 //constante de Boltzmann et tempÃ©rature (K)
# define D_parallel 100000//3.38*pow(10,6) //coefficients de diffusion
# define D_orthogonal 100000//1.75*pow(10,6)
# define D_THETA 100
# define DT 0.00001 //time step
# define N_BRWN 100001 //500000 //nb de steps
# define Taille_Systeme 1000.0
# define N_R 100000//100000 //nombre de rÃ©cepteurs
# define N_L 200

#define R_Verlet 15.0 //rayon pour la liste de Verlet 




bool isodd(int num) //vÃ©rifie la paritÃ© d'un entier
{
	if(num % 2 == 0)
		return false;
	else
		return true;
}



float ran2(long *idum) //idum doit être un nombre négatif
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}



float RGV(long *g) //trouve le nombre alÃ©atoire pour le mvt brownien
{
	std::vector<float> random_numbers;
	float p;

    for(int j=0; j<2; j++)
    {
        p = ran2(g);
        random_numbers.push_back(p);
    }
 	float r1 = random_numbers[0];
	float r2 = random_numbers[1];

  	float r3 = sqrt(-2*(log(1-r1)));
  	float r = cos(r2*2*M_PI)*r3;

  	return r;
}



/*void initialize_IAV(std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<int> *T_l, std::vector<int> *L_to_R)
{//Initialise la position des ligands dans le repÃ¨re du CM du virus
	float b = D;
	float d = (sqrt(3)*b)/2;
	std::vector<float> y0;
	std::vector<float> y1;
	std::vector<float> y2;



	for (int i = 0; i <= int(L/(2*D)); i++)
	{
		y0.push_back(-L/2 + i*2*D);
	}

	for (int j = 0; j <= int(L/(2*D)+1); j++)
	{
		y1.push_back(-L/2 - D + j*2*D);
	}

	for (int k = 0; k <= int(L/(2*D)+2); k++)
	{
		y2.push_back(-L/2 - 2*D + k*2*D);
	}

   

	for (int a = 0; a < 9; a++)
	{
		if(a==0)
		{
			for (int l = 0; l < y0.size(); l++)
			{
				X_l_x->push_back(-2*b);
				X_l_y->push_back(y0.at(l));
			}
		}
		else if(a==8)
		{
			for (int m = 0; m < y0.size(); m++)
			{
				X_l_x->push_back(2*b);
				X_l_y->push_back(y0.at(m));
			}
		}
		else if(isodd(a)==true)
		{
			for (int n = 0; n < y1.size(); n++)
			{
				X_l_x->push_back((a-4)*b/2);
				X_l_y->push_back(y1.at(n));
			}
		}
		else if(isodd(a)==false)
		{
			for (int o = 0; o < y2.size(); o++)
			{
				X_l_x->push_back((a-4)*b/2);
				X_l_y->push_back(y2.at(o));
			}
		}	
	}

	for (int p = 0; p < X_l_y->size(); p++)
	{
		if(X_l_y->at(p) < (-2*L/5))	
			T_l->push_back(0);
		else
			T_l->push_back(1);
	
	    L_to_R->push_back(-1);
	}
}*/

void initialize_Virus(std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<int> *T_l, std::vector<int> *L_to_R, long *m)
{
	float x, y;

	while(X_l_x->size() < N_L)
	{
		x = ran2(m)*L - L/2;
		y = ran2(m)*L - L/2;

		if(x*x + y*y <= L/2 * L/2)
		{
			X_l_x->push_back(x);
			X_l_y->push_back(y);		
		}
	}

	for (int j = 0; j < X_l_x->size(); j++)
	{
		T_l->push_back(1);
		L_to_R->push_back(-1);
	}
}


void initialize_SA(std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_r, long *m)
{//initialise les rÃ©cpteurs de maniÃ¨re alÃ©atoire dans le plan
	float x;
	float y;
	int check;
	
 	x = ran2(m)*Taille_Systeme - Taille_Systeme/2;
	y = ran2(m)*Taille_Systeme - Taille_Systeme/2;
		
	X_r_x->push_back(x);
    X_r_y->push_back(y);	

	while(X_r_x->size() < N_R) 
	{
		x = ran2(m)*Taille_Systeme - Taille_Systeme/2;
		y = ran2(m)*Taille_Systeme - Taille_Systeme/2;
	    check = 0;

		X_r_x->push_back(x);
		X_r_y->push_back(y);  
	}

	for (int j = 0; j < X_r_x->size(); j++)
	{
		T_r->push_back(1);
	}
}





bool check_breakage(std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, float DX, float DY, float DTHETA)
{//vÃ©rifie si la diffusion ne brise pas un pont
	bool bk = false;
	int receptor = 0;
    float theta;
    float x;
    float y;

    theta = DTHETA + X_CM_theta->at((X_CM_theta->size()) - 1);

	for (int i = 0; i < X_l_x->size(); ++i)
	{
		if (T_l->at(i) == 2)
		{
			receptor = L_to_R->at(i);
		    
		    x = (X_l_x->at(i))*cos(theta) - (X_l_y->at(i))*sin(theta) + X_CM_x->at((X_CM_x->size()) -1) ;
		    y = (X_l_x->at(i))*sin(theta) + (X_l_y->at(i))*cos(theta) + X_CM_y->at((X_CM_y->size()) -1) ;
		    
		    
			if (pow(x + DX - X_r_x->at(receptor - 1),2) + pow(y + DY - X_r_y->at(receptor - 1),2) > D*D)
			{
				//std::cout << pow(x + DX - X_r_x->at(receptor - 1),2) + pow(y + DY - X_r_y->at(receptor - 1),2) << '\n';
				bk = true;
			}
		}
	}
	
	return bk;
}



void diffusion(std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_l, std::vector<int> *L_to_R, long *idum, float delta_t)
{//effectue la diffusion aprÃ¨s avoir vÃ©rifiÃ© si elle est autorisÃ©e
    std::vector<float> random_numbers;

    for(int k = 0; k<3; k++)
    {
        random_numbers.push_back(RGV(idum));
    }
    
    float old_theta = X_CM_theta->at((X_CM_theta->size()) - 1);
	
	float DX_Brwn = sqrt(2*D_parallel*delta_t)*random_numbers[0];    
	//float DX_Brwn = sqrt(2*D_parallel*delta_t)*cos(old_theta)*random_numbers[0] - sqrt(2*D_orthogonal*delta_t)*sin(old_theta)*random_numbers[1];
	// version de Fletcher: (sqrt(2*D_parallel*delta_t)*pow(cos(old_theta),2) + sqrt(2*D_orthogonal*delta_t)*pow(sin(old_theta),2))*random_numbers[0] + ((sqrt(2*D_parallel*delta_t) - sqrt(2*D_orthogonal*delta_t))*cos(old_theta)*sin(old_theta))*random_numbers[1];
	float DY_Brwn = sqrt(2*D_parallel*delta_t)*random_numbers[1];   
	//float DY_Brwn = sqrt(2*D_parallel*delta_t)*sin(old_theta)*random_numbers[0] + sqrt(2*D_orthogonal*delta_t)*cos(old_theta)*random_numbers[1]; 
	// version de Fletcher: ((sqrt(2*D_parallel*delta_t) - sqrt(2*D_orthogonal*delta_t))*cos(old_theta)*sin(old_theta))*random_numbers[0] + (sqrt(2*D_parallel*delta_t)*pow(sin(old_theta),2) + sqrt(2*D_orthogonal*delta_t)*pow(cos(old_theta),2))*random_numbers[1];
	float DTHETA_Brwn = sqrt(2*D_THETA*delta_t)*random_numbers[2];
	
	bool bk;
	

    
    bk = check_breakage(X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, T_l, L_to_R, DX_Brwn, DY_Brwn, DTHETA_Brwn);
        
        if(bk == false) //on effectue la diffusion
        {
        	X_CM_x->push_back(X_CM_x->at((X_CM_x->size()) -1) + DX_Brwn);
        	X_CM_y->push_back(X_CM_y->at((X_CM_y->size()) -1) + DY_Brwn);
        	X_CM_theta->push_back(X_CM_theta->at((X_CM_theta->size()) -1) + DTHETA_Brwn);
            //std::cout << "check breakage false " << '\n'; 
        }   
        
        if(bk == true) //le systÃ¨me ne bouge pas
        {
            X_CM_x->push_back(X_CM_x->at(X_CM_x->size() - 1));
        	X_CM_y->push_back(X_CM_y->at(X_CM_y->size() - 1));
        	X_CM_theta->push_back(X_CM_theta->at(X_CM_theta->size() - 1));
            //std::cout << "check breakage true" << '\n'; 
        }
}



int Cell_to_Index(int N_x, int N_y, int ix, int iy) //reÃ§oit les coordonnÃ©es de la cellule et renvoie l'indice de la cellule
{
	int ixp;
	int iyp;
	int nc;

	ixp= ix;
	iyp = iy;

	if (ix<1)
	{
		ixp = N_x;
	}
	if (ix>N_x)
	{
		ixp = 1;
	}
	if (iy<1)
	{
		iyp = N_y;
	}
	if (iy>N_y)
	{
		iyp = 1;
	}
	nc = (iyp - 1)*N_x + ixp;

	return nc;
}



void Initialize_L_L(MATRIX_INT *L_L, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> ll_R, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, float R_Cell, int step)
{//Initialise L_L en utilisant les listes de voisinage
    float theta;
    float x, y;
	for (int alpha = 0; alpha < N_L; ++alpha)
	{ 
		theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

		int ix = int((x + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((y + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        
		for (int j = 0; j < 9; ++j)
		{
			int i = Neigh[index-1][j];
			
			if(hoc_R[i] != 0)
			{
	            if(pow(x - X_r_x[hoc_R[i] - 1],2) + pow(y - X_r_y[hoc_R[i] - 1],2) < R_Verlet*R_Verlet)
    				(L_L->at(alpha)).push_back(hoc_R[i]); 
				int g = ll_R[hoc_R[i]];
				
				while(g != 0)
				{
	                if(pow(x - X_r_x[g - 1],2) + pow(y - X_r_y[g - 1],2) < R_Verlet*R_Verlet)
					    (L_L->at(alpha)).push_back(g);
					int k = g;	
					g = ll_R[k];
				}
			}
		}
		if(T_l[alpha] == 2 && pow(x - X_r_x[L_to_R[alpha] - 1],2) + pow(y - X_r_y[L_to_R[alpha] - 1],2) > R_Verlet*R_Verlet) //gère l'exception due à l'imprécision dans check_breakage
			(L_L->at(alpha)).push_back(L_to_R[alpha]);

		int N = (L_L->at(alpha)).size();
		(L_L->at(alpha)).insert((L_L->at(alpha)).begin(), N);
	}
}



void Initialize_L_R(MATRIX_INT *L_R, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_L, std::vector<int> ll_L, std::vector<int> T_r, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, float R_Cell, int step)
{//initialise L_R en utilisant les listes de voisinage
    float theta;
    float x, y;
    float count = 0;

    //std::cout << "un" << '\n'; 

	for (int alpha = 0; alpha < N_R; ++alpha)
	{ 
		int ix = int((X_r_x[alpha] + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((X_r_y[alpha] + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        
		for (int j = 0; j < 9; ++j)
		{
			int i = Neigh[index - 1][j];
			
			if(hoc_L[i] != 0)
			{
				theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    		x = (X_l_x->at(hoc_L[i] - 1))*cos(theta) - (X_l_y->at(hoc_L[i] - 1))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
				y = (X_l_x->at(hoc_L[i] - 1))*sin(theta) + (X_l_y->at(hoc_L[i] - 1))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

	            if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) < R_Verlet*R_Verlet)
    				(L_R->at(alpha)).push_back(hoc_L[i]); 
		    	
		    	//std::cout << alpha << '\n'; 			
				
				int g = ll_L[hoc_L[i]];
				
				count = 0;
				while(g != 0)
				{
					//std::cout << count << " " << ll_L.size() << '\n';
					theta = X_CM_theta.at((X_CM_theta.size()) - 1);
				    x = (X_l_x->at(g-1))*cos(theta) - (X_l_y->at(g-1))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
					y = (X_l_x->at(g-1))*sin(theta) + (X_l_y->at(g-1))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

	                if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) < R_Verlet*R_Verlet)
					    (L_R->at(alpha)).push_back(g);
					int k = g;	
					g = ll_L[k];
					count += 1;

					/*if(count >= ll_L.size())
					{
						std::cout << "ERROR: " << "receptor: " << alpha << " cellule: " << i << '\n' << "[ ";
						for (int w = 0; w<ll_L.size(); w++)
						{
							std::cout << ll_L[w] << " ";
						}
						std::cout << "]" << '\n' << "[";
						for (int w = 0; w < hoc_L.size(); w++)
						{
							std::cout << hoc_L[w] << " ";
						}
						std::cout << " ]" << '\n';
					}*/
				}

			}
		}

		
		if(T_r[alpha] == 2) //corrige l'exception due à l'imprécision numérique dans check_breakage
		{
			for(int z = 0; z < 199; z++)
			{
				if(T_l[z] == 2 && L_to_R[z] == alpha + 1)
				{
					theta = X_CM_theta.at((X_CM_theta.size()) - 1);
					x = (X_l_x->at(z))*cos(theta) - (X_l_y->at(z))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
					y = (X_l_x->at(z))*sin(theta) + (X_l_y->at(z))*cos(theta) + X_CM_y[X_CM_y.size() - 1];
					if(pow(x - X_r_x[alpha],2) + pow(y - X_r_y[alpha],2) > R_Verlet*R_Verlet)
						(L_R->at(alpha)).push_back(z+1);
				}
			}
		}

		int N = (L_R->at(alpha)).size();
		(L_R->at(alpha)).insert((L_R->at(alpha)).begin(), N);
	}
}


void Initialize_Affinity(MATRIX_FLOAT *a, MATRIX_INT *L_L, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *L_to_R, std::vector<int> T_l, std::vector<int> T_r, float k, int t1, int t2)
{//initialise les matrices a_d, a_on et a_off
	float a_sum = 0;
	float theta, x, y;

	for (int alpha = 0; alpha < N_L; ++alpha)
	{
		theta = X_CM_theta->at((X_CM_theta->size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
	
		a_sum = 0;
		for (int i = 1; i < (L_L->at(alpha)).size(); ++i)
		{
			if((pow(x - X_r_x->at(L_L->at(alpha)[i] - 1),2) + pow(y - X_r_y->at(L_L->at(alpha)[i] - 1),2) < D*D))
			{
		        if(T_l[alpha] == t1 && T_r[((L_L->at(alpha)).at(i)) - 1] == t2)
		        {
	    	       	if(t1 == 2 && t2 == 2) //gère le cas où on a un pont
	    	       	{
	    	       	    if(L_to_R->at(alpha) == L_L->at(alpha)[i])
	    	       	    {
	    	       	        a_sum += k;
			                (a->at(alpha)).push_back(k);
	    	       	    }
	    	       	    else
	    	       	    {
		                    (a->at(alpha)).push_back(0);
	    	       	    }
	    	       	}
	    	       	else //cas où on a pas de pont
	    	       	{
			            a_sum += k;
			            (a->at(alpha)).push_back(k);
	    	       	}
		        }
		        else
		        {
		            (a->at(alpha)).push_back(0);
			    }
			}
			else
			{
				(a->at(alpha)).push_back(0);
			}    
		}
		(a->at(alpha)).insert((a->at(alpha)).begin(), a_sum);
	}
}



std::vector<int> Index_to_Cell(int N_x, int N_y, int nc) //reÃ§oit l'indice de la cellule et renvoie ses coordonnÃ©es
{
	std::vector<int> ixy;
	int ix;
	int iy;

	iy = (nc-1)/N_x + 1;
	ix = nc - (iy - 1)*N_x;

	ixy.push_back(ix);
	ixy.push_back(iy);

	return ixy;
}



MATRIX_INT New_List(float R_Cell, int N_Cell, int N_part, int N_x, int N_y, std::vector<float> *X_part_x, std::vector<float> *X_part_y)
{//crÃ©e les listes hoc, ll et bl
	std::vector<int> hoc(N_Cell + 1);
	std::vector<int> ll(N_part + 1);
	std::vector<int> bl(N_part + 1);
	MATRIX_INT v;

	for (int icell = 1; icell < N_Cell + 1; ++icell)
	{
		hoc[icell] = 0;
	}
	for (int i = 1; i < N_part + 1; ++i)
	{
		int ix = int((X_part_x->at(i-1) + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((X_part_y->at(i-1) + Taille_Systeme/2)/R_Cell) + 1;
		int nc = Cell_to_Index(N_x,N_y,ix,iy);

		ll[i] = hoc[nc];
		if (hoc[nc] != 0)
		{
			bl[hoc[nc]] = i;
		}
		hoc[nc] = i;
	}
	
	v.push_back(ll);
	v.push_back(bl);
	v.push_back(hoc);
	
	return v;
}


MATRIX_INT Set_Neigh(int N_Cell, int N_x, int N_y) //crÃ©e les listes de voisinage
{
	int NNeigh = 9;
	std::vector<int> ixy(2);

	MATRIX_INT Neigh;
	Neigh.resize(N_Cell, std::vector<int>(NNeigh));

	for (int i = 0; i < N_Cell; ++i)
	{
		ixy = Index_to_Cell(N_x, N_y, i+1);
		int ix = ixy[0];
		int iy = ixy[1];
		int cnt = 0;
		for (int icx = ix-1; icx < ix + 2; ++icx)
		{
			for (int icy = iy-1 ; icy < iy + 2; ++icy)
			{
				int index = Cell_to_Index(N_x,N_y,icx,icy);
				Neigh[i][cnt] = index;
				cnt=cnt+1;
			}
		}
	}
	return Neigh;
}





void UpGrade(int Old_Cell, int New_Cell, int ipar, std::vector<int> *ll, std::vector<int> *bl, std::vector<int> *hoc)
{//met Ã  jour les cellules si des particules rentrent ou sortent
	//supprimer particule de l'ancienne cellule
	if (ll->at(ipar) != 0)
	{
		if (bl->at(ipar) != 0) //particule au milieu de la cellule
		{
			ll->at(bl->at(ipar)) = ll->at(ipar);
			bl->at(ll->at(ipar)) = bl->at(ipar);
		}
		else //particule est la premiÃ¨re de la cellule
		{
			bl->at(ll->at(ipar)) = 0;
			hoc->at(Old_Cell) = ll->at(ipar);
		}
	}
	else
	{
		if (bl->at(ipar) != 0) //particule est la derniÃ¨re de la celulle
		{
			ll->at(bl->at(ipar)) = 0;
		}
		else //particule est la seule dans la cellule
		{
			hoc->at(Old_Cell) = 0;
		}
	}


	//ajouter particule Ã  la nouvelle cellule
	bl->at(ipar) = 0;
	ll->at(ipar) = hoc->at(New_Cell);
	if (hoc->at(New_Cell) !=0)
	{
		bl->at(hoc->at(New_Cell)) = ipar;
	}
	hoc->at(New_Cell) = ipar;
}




int Choose(std::vector<float> V, long *g) //choisit une rÃ©action, un ligand ou un rÃ©cepteur
{
	float r = ran2(g); 
	float sumpar = 0.0;
	int i = 0;

        if(V[0] == 0)
        return -1;
        
	while(sumpar <= r && i<V.size()-1) 
        { 
		sumpar += (V[i+1]/V[0]);
                i += 1;
	}
	
	return i;
}


int Get_Id(std::vector<int> L_L_line, int i) //donne la position du rÃ©cepteur i dans le vecteur L_L[alpha] (alpha in [0,N_L-1], i in [1,N_R])
{
    int id = 1;
    for(int w = 0; w < (L_L_line.size()); w++)
    {
        if(L_L_line[w + 1] == i)
        {
            return id;
        }
        else
        {
            id += 1;   
        }
    }
    return -1; //si erreur
}


std::vector<float> Get_a_X_tot(MATRIX_FLOAT *a_X)  //crÃ©e le vecteur nÃ©cessaire pour que Choose puisse choisir un ligand une fois la rÃ©action choisie
{
	std::vector<float> a_X_tot;
	a_X_tot.clear();
	float a_X_sum = 0.0;

	for (int i = 0; i < N_L; ++i)
	{
		a_X_tot.push_back((a_X->at(i)).at(0));
		a_X_sum += (a_X->at(i)).at(0);
	}
	a_X_tot.insert(a_X_tot.begin(), a_X_sum);
	return a_X_tot;
}


std::vector<float> Get_a_tot(std::vector<float> *a_d_tot, std::vector<float> *a_on_tot, std::vector<float> *a_off_tot) 
{//crÃ©e le vecteur nÃ©cessaire pour choisir une rÃ©action
	float a_tot;
	float a_d = 0.0;
	float a_on = 0.0;
	float a_off = 0.0;
	float epsilon = 0.0000001;
	std::vector<float> a;
	a.clear();

	for (int i = 0; i < N_L; ++i)
	{
		a_d += a_d_tot->at(i+1);
		if(a_d < epsilon)
		{
		    a_d = 0.0;
		}
	}

	for (int j = 0; j < N_L; ++j)
	{
		a_on += a_on_tot->at(j+1);
		if(a_on < epsilon) //Ã  remplacer par k_on quand on remet k_on !=0
		{
		    a_on = 0.0;
		}
	}

	for (int k = 0; k < N_L; ++k)
	{
		a_off += a_off_tot->at(k+1);
		if(a_off < epsilon)
		{
		    a_off = 0.0;
		}
	}

	a_tot = a_d + a_on + a_off;
	a.push_back(a_tot);
	a.push_back(a_d);
	a.push_back(a_on);
	a.push_back(a_off);

	return a;
}


float Get_time_reaction(float a_tot, long *idum)
{
    float t;
    t =  (-1)*(log(1 - ran2(idum)))/a_tot;
    
    return t;
}



//ATTENTION ! MAINTENANT ON DOIT MODIFIER LES MÀJ CAR LES a_X ONT CHANGÉ DE FORMAT !!!!!!
void Make_reaction(MATRIX_FLOAT *a_d, MATRIX_FLOAT *a_on, MATRIX_FLOAT *a_off, std::vector<float> a_d_tot, std::vector<float> a_on_tot, std::vector<float> a_off_tot, MATRIX_INT L_L, MATRIX_INT L_R, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *T_l, std::vector<int> *T_r, std::vector<int> *L_to_R, int r, long *idum, int step)
{//choisit le ligand et le rÃ©cepteur et effectue la rÃ©action choisie auparavant 
    float a_sum;
    float theta, x, y;
	if (r == 1) //destruction NA
	{
		int alpha = Choose(a_d_tot, idum); //choix du ligand qui rÃ©agit

		int id = Choose((a_d->at(alpha-1)), idum); //choix du rÃ©cepteur qui rÃ©agit
		int i = L_L[alpha-1][id];
		

		//MISE A JOUR
		(a_d->at(alpha-1))[0] -= (a_d->at(alpha-1))[id];
		(a_d->at(alpha-1))[id] = 0;

		for (int j = 1; j < L_R[i-1].size(); j++) 
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(L_R[i-1][j] - 1))*cos(theta) - (X_l_y->at(L_R[i-1][j] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(L_R[i-1][j] - 1))*sin(theta) + (X_l_y->at(L_R[i-1][j] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);

			if((pow(x - X_r_x->at(i-1),2) + pow(y - X_r_y->at(i-1),2) < D*D))
			{
				//clock_t startTime = clock();
				int id2 = Get_Id(L_L[(L_R[i-1][j]) - 1], i);
				//std::cout << "test L_R id - r1" << '\n';
				//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
			    if(T_l->at(L_R[i-1][j] - 1) == 0)
				{
					(a_d->at(L_R[i-1][j] - 1))[0] -= (a_d->at(L_R[i-1][j] - 1))[id2];
				    (a_d->at(L_R[i-1][j] - 1))[id2] = 0.0;
				}
				if(T_l->at(L_R[i-1][j] - 1) == 1)
				{
					(a_on->at(L_R[i-1][j] - 1))[0] -= (a_on->at(L_R[i-1][j] - 1))[id2];
					(a_on->at(L_R[i-1][j] - 1))[id2] = 0.0;
				}
			}
		}
	    T_r->at(i-1) = 0;
	}


	else if(r == 2) //crÃ©ation pont
	{
		int alpha = Choose(a_on_tot, idum); //choix du ligand qui rÃ©agit
		int id = Choose((a_on->at(alpha - 1)), idum); //choix du rÃ©cepteur qui rÃ©agit
		int i = L_L[alpha - 1][id];

		//MISE A JOUR
		(a_off->at(alpha - 1))[0] += k_off;
		(a_off->at(alpha - 1))[id] = k_off; 
	
		for (int k = 1; k < L_R[i-1].size(); k++)
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(L_R[i-1][k] - 1))*cos(theta) - (X_l_y->at(L_R[i-1][k] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(L_R[i-1][k] - 1))*sin(theta) + (X_l_y->at(L_R[i-1][k] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);

			if((pow(x - X_r_x->at(i-1),2) + pow(y - X_r_y->at(i-1),2) < D*D))
			{
				//clock_t startTime = clock();
				int id2 = Get_Id(L_L[(L_R[i-1][k]) - 1], i);
				//std::cout << "test L_R id - r2" << '\n';
				//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';			    
				if (T_l->at(L_R[i-1][k] - 1) == 1)
				{
					(a_on->at(L_R[i-1][k] - 1))[0] -= (a_on->at(L_R[i-1][k] - 1))[id2];
					(a_on->at(L_R[i-1][k] - 1))[id2] = 0.0;				
				}
				else if (T_l->at(L_R[i-1][k] - 1) == 0)
				{
					(a_d->at(L_R[i-1][k] - 1))[0] -= (a_d->at(L_R[i-1][k] - 1))[id2];
					(a_d->at(L_R[i-1][k] - 1))[id2] = 0.0;						
				}
			}
		}
		
		for (int k = 0; k < (a_on->at(alpha-1)).size(); k++)
		{
		    a_on->at(alpha-1).at(k) = 0.0;
		}
		
		T_l->at(alpha - 1) = 2;
		T_r->at(i - 1) = 2;
		L_to_R->at(alpha - 1) = i;	
	}


	else if(r == 3) //destruction pont
	{
		int alpha = Choose(a_off_tot, idum); //choix du ligand qui rÃ©agit

		int id = Choose((a_off->at(alpha-1)), idum); //choix de la position du rÃ©cepteur qui rÃ©agit
		int i = L_L[alpha-1][id]; //indice du rÃ©cepteur qui rÃ©agit


		//MISE A JOUR
		(a_off->at(alpha-1))[0] = 0.0;
		(a_off->at(alpha-1))[id] = 0.0; 

		T_l->at(alpha-1) = 1;
		T_r->at(i - 1) = 1;

		for (int m = 1; m < L_R[i-1].size(); m++)
		{

			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(L_R[i-1][m] - 1))*cos(theta) - (X_l_y->at(L_R[i-1][m] - 1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(L_R[i-1][m] - 1))*sin(theta) + (X_l_y->at(L_R[i-1][m] - 1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);

			if((pow(x - X_r_x->at(i-1),2) + pow(y - X_r_y->at(i-1),2) < D*D))
			{
				//clock_t startTime = clock();
				int id2 = Get_Id(L_L[(L_R[i-1][m]) - 1], i);
				//std::cout << "test L_R id - r3" << '\n';
				//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';		

				if (T_l->at(L_R[i-1][m] - 1) == 1 && L_R[i-1][m] != alpha)	
				{
					(a_on->at(L_R[i-1][m] - 1))[0] += k_on;
					(a_on->at(L_R[i-1][m] - 1))[id2] = k_on;				
				}
				else if (T_l->at(L_R[i-1][m] - 1) == 0)
				{
					(a_d->at(L_R[i-1][m] - 1))[0] += k_d;
					(a_d->at(L_R[i-1][m] - 1))[id2] = k_d;						
				}
			}
		}
		
		a_sum = 0;
		for (int k = 1; k < L_L[alpha-1].size(); k++)
		{
			theta = X_CM_theta->at((X_CM_theta->size()) - 1);
		    x = (X_l_x->at(alpha-1))*cos(theta) - (X_l_y->at(alpha-1))*sin(theta) + X_CM_x->at(X_CM_x->size() - 1);
			y = (X_l_x->at(alpha-1))*sin(theta) + (X_l_y->at(alpha-1))*cos(theta) + X_CM_y->at(X_CM_y->size() - 1);
		
		    if ((T_r->at(L_L[alpha-1][k] - 1) == 1)  && (pow(x - X_r_x->at(L_L[alpha-1][k] - 1),2) + pow(y - X_r_y->at(L_L[alpha-1][k] - 1),2) < D*D))
		    {
		        a_on->at(alpha-1)[k] += k_on;
		        a_sum += k_on;
		    }
		    else
		    {
		        a_on->at(alpha-1)[k] = 0.0;
		    }
		}
		a_on->at(alpha-1)[0] = a_sum; 
		L_to_R->at(alpha - 1) = -1;
		
	}    
	
	else
	{
	    std::cout << "error in choice of reaction: no reaction has been chosen" << '\n';
	}
}



void Update_Cell_L(std::vector<int> *ll, std::vector<int> *bl, std::vector<int> *hoc, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<int> *check, int N_x, int N_y, float R_Cell, float X_CM_x_init, float X_CM_y_init, float X_CM_theta_init)
{//met à jour toutes les cellules après déplacement du virus
	float theta = X_CM_theta->at(X_CM_theta->size() - 1);
	float otheta = X_CM_theta_init;


	for (int i = 1; i < N_L + 1; i++)
	{

				float x = (X_l_x->at(i-1))*cos(theta) - (X_l_y->at(i-1))*sin(theta);
				float y = (X_l_x->at(i-1))*sin(theta) + (X_l_y->at(i-1))*cos(theta);
				    
				float x_p = x + X_CM_x->at(X_CM_x->size() - 1);
				float y_p = y + X_CM_y->at(X_CM_y->size() - 1);

				int ix = int((x_p + Taille_Systeme/2)/R_Cell) + 1;
				int iy = int((y_p + Taille_Systeme/2)/R_Cell) + 1;
				int nc = Cell_to_Index(N_x,N_y,ix,iy);



				float ox = (X_l_x->at(i-1))*cos(otheta) - (X_l_y->at(i-1))*sin(otheta);
				float oy = (X_l_x->at(i-1))*sin(otheta) + (X_l_y->at(i-1))*cos(otheta);
				    
				float ox_p = ox + X_CM_x_init;
				float oy_p = oy + X_CM_y_init;

				int oix = int((ox_p + Taille_Systeme/2)/R_Cell) + 1;
				int oiy = int((oy_p + Taille_Systeme/2)/R_Cell) + 1;
				int oc = Cell_to_Index(N_x,N_y,oix,oiy);

				if(oc != nc)
				{
					//std::cout << "test a" << "\n";
					UpGrade(oc, nc, i, ll, bl, hoc);
					//std::cout << "test b" << '\n';
				}
			
		
	}

}



void Update_Lists(MATRIX_INT *L_L, MATRIX_INT *L_R, MATRIX_FLOAT *a_d, MATRIX_FLOAT *a_on, MATRIX_FLOAT *a_off, std::vector<float> *a_d_tot, std::vector<float> *a_on_tot, std::vector<float> *a_off_tot, std::vector<float> *affinity, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> hoc_L, std::vector<int> ll_L, std::vector<int> ll_R, std::vector<float> X_l_x, std::vector<float> X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, int N_x, int N_y, float R_Cell, std::vector<int> T_l, std::vector<int> T_r, std::vector<int> *L_to_R, int step)
{//update L_L, L_R, a_d, a_on, a_off en les vidant, puis en les remplissant Ã  nouveau avec les fonctions Initialize
    for(int i=0; i<N_L; i++)
    {
        (L_L->at(i)).clear();
        (a_d->at(i)).clear();
        (a_on->at(i)).clear();
        (a_off->at(i)).clear();
    }

    for(int i=0; i<N_R; i++)
    {
        (L_R->at(i)).clear();
    }

    a_d_tot->clear();
    a_on_tot->clear();
    a_off_tot->clear();
    affinity->clear();
    
	Initialize_L_L(L_L, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, *L_to_R, N_x, N_y, R_Cell, step);
	Initialize_L_R(L_R, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_L, ll_L, T_r, T_l, *L_to_R, N_x, N_y, R_Cell, step);


	Initialize_Affinity(a_d, L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_d, 0, 1);
	Initialize_Affinity(a_on, L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_on, 1, 1);
    Initialize_Affinity(a_off, L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, L_to_R, T_l, T_r, k_off, 2, 2);

    *a_d_tot = Get_a_X_tot(a_d);
	*a_on_tot = Get_a_X_tot(a_on);
	*a_off_tot = Get_a_X_tot(a_off);

	*affinity = Get_a_tot(a_d_tot, a_on_tot, a_off_tot);
}


void Update_L_L(MATRIX_INT *L_L, MATRIX_INT *L_R, MATRIX_INT *L_L_old, std::vector<int> *check, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> X_r_x, std::vector<float> X_r_y, std::vector<float> X_CM_x, std::vector<float> X_CM_y, std::vector<float> X_CM_theta, MATRIX_INT Neigh, std::vector<int> hoc_R, std::vector<int> ll_R, std::vector<int> T_l, std::vector<int> L_to_R, int N_x, int N_y, float R_Cell)
{
	for(int k = 0; k < N_L; k++)
	{	
		(L_L_old->at(k)).clear();
		L_L_old->at(k) = L_L->at(k);
	}

	float theta;
    float x, y;
    int alpha;

	//std::cout << "test 1bis" << '\n';

	for (int beta = 0; beta < check->size(); ++beta)
	{ 
		//std::cout << check->size() << '\n';
		//std::cout << beta << '\n';
		alpha = check->at(beta) - 1;
		(L_L->at(alpha)).clear();
		theta = X_CM_theta.at((X_CM_theta.size()) - 1);
	    x = (X_l_x->at(alpha))*cos(theta) - (X_l_y->at(alpha))*sin(theta) + X_CM_x[X_CM_x.size() - 1];
		y = (X_l_x->at(alpha))*sin(theta) + (X_l_y->at(alpha))*cos(theta) + X_CM_y[X_CM_y.size() - 1];

		int ix = int((x + Taille_Systeme/2)/R_Cell) + 1;
		int iy = int((y + Taille_Systeme/2)/R_Cell) + 1;
		int index = Cell_to_Index(N_x, N_y, ix, iy);
        
		for (int j = 0; j < 9; ++j)
		{
			int i = Neigh[index-1][j];
			
			if(hoc_R[i] != 0)
			{
	            if(pow(x - X_r_x[hoc_R[i] - 1],2) + pow(y - X_r_y[hoc_R[i] - 1],2) < R_Verlet*R_Verlet)
    				(L_L->at(alpha)).push_back(hoc_R[i]); 
				int g = ll_R[hoc_R[i]];
				
				while(g != 0)
				{
	                if(pow(x - X_r_x[g - 1],2) + pow(y - X_r_y[g - 1],2) < R_Verlet*R_Verlet)
					    (L_L->at(alpha)).push_back(g);
					int k = g;	
					g = ll_R[k];
				}
			}
		}
		if(T_l[alpha] == 2 && pow(x - X_r_x[L_to_R[alpha] - 1],2) + pow(y - X_r_y[L_to_R[alpha] - 1],2) > R_Verlet*R_Verlet) //gère l'exception due à l'imprécision dans check_breakage
			(L_L->at(alpha)).push_back(L_to_R[alpha]);

		int N = (L_L->at(alpha)).size();
		(L_L->at(alpha)).insert((L_L->at(alpha)).begin(), N);
	}
}

bool check_intersection(int element, std::vector<int> v)
{
	for(int j = 0; j < v.size(); j++)
	{
		if(element == v[j])
		{
			return true;
		}		
	}
	return false;
}


void Update_L_R(MATRIX_INT *L_L_old, MATRIX_INT *L_L, MATRIX_INT *L_R, std::vector<int> *check)
{
	clock_t startTime = clock();
	int alpha;
	bool check_inter;
	std::vector<int> v_old;
	std::vector<int> v_new;

	for (int beta = 0; beta < check->size(); ++beta)
	{ 
		alpha = check->at(beta);

		v_old = L_L_old->at(alpha - 1);
		v_new = L_L->at(alpha - 1);

		v_old.erase(v_old.begin());
		v_new.erase(v_new.begin());

	    std::sort(v_old.begin(), v_old.end());
	    std::sort(v_new.begin(), v_new.end());
	 
	    std::vector<int> v_intersection;
	 
	    std::set_intersection(v_old.begin(), v_old.end(),
	                          v_new.begin(), v_new.end(),
	                          std::back_inserter(v_intersection));


	    /*std::cout << '\n' << "v_new: ";
	    for(int w = 0; w < v_new.size(); w++)
	    {
	    	std::cout << v_new[w] << " ";
	    }
	    std::cout << '\n';

	    std::cout << '\n' << "v_old: ";
	    for(int w = 0; w < v_old.size(); w++)
	    {
	    	std::cout << v_old[w] << " ";
	    }
	    std::cout << '\n';

	    std::cout << '\n' << "v_intersection: ";
	    for(int w = 0; w < v_intersection.size(); w++)
	    {
	    	std::cout << v_intersection[w] << " ";
	    }
	    std::cout << '\n';*/

	    //ENLEVER LIGAND POUR LES RECEPTEURS SORTANTS
		//clock_t startTime2 = clock();
	    for(int j = 0; j < v_old.size(); j++)
	    {
	    	if(v_intersection.size() != 0)
	    	{
	    		//clock_t startTime6 = clock();
	    		check_inter = check_intersection(v_old[j],v_intersection); //REGARDER LA DIFFERENCE ENTRE CHECk_INTER ET GET_ID, PQ 20x LE TEMPS ??? 
				//std::cout << "test L_R inter" << '\n';
				//std::cout << double( clock() - startTime6 ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
	    		if(check_inter == false)
	    		{
	    				//clock_t startTime3 = clock();
					    int i = v_old[j];
	    				int id = Get_Id(L_R->at(i-1), alpha);

						//std::cout << "test L_R id1" << '\n';
						//std::cout << double( clock() - startTime3 ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';

	    				(L_R->at(i-1)).erase(L_R->at(i-1).begin() + id);
	    				(L_R->at(i-1).at(0)) -= 1;
	    		}
	    	}
	    	else
	    	{
	    		clock_t startTime4 = clock();

		    	int i = v_old[j];
				int id = Get_Id(L_R->at(i-1), alpha);

				//std::cout << "test L_R id2" << '\n';
				//std::cout << double( clock() - startTime4 ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
				
				(L_R->at(i-1)).erase(L_R->at(i-1).begin() + id);	
				(L_R->at(i-1).at(0)) -= 1;
	    	}
	    }
		//std::cout << "test L_R OUT" << '\n';
		//std::cout << double( clock() - startTime2 ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';


	    //AJOUTER LIGAND POUR LES RECEPTEURS ENTRANTS
		//clock_t startTime5 = clock();
	    for(int j = 0; j < v_new.size(); j++)
	    {
	    	if(v_intersection.size() != 0)
	    	{
	    		check_inter = check_intersection(v_new[j],v_intersection);
	    		if(check_inter == false)
	    		{
	    				int i = v_new[j];
	    				(L_R->at(i-1)).push_back(alpha);
	    				L_R->at(i-1).at(0) += 1;
	    		}
	    	}
	    	else
	    	{
		    	int i = v_new[j];
				(L_R->at(i-1)).push_back(alpha);	
				L_R->at(i-1).at(0) += 1;
	    	}
	    }
		//std::cout << "test L_R IN" << '\n';
		//std::cout << double( clock() - startTime5 ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
	}

	
	v_old.clear();
	v_new.clear();
	//std::cout << "test L_R" << '\n';
	//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
}


void Update_Affinity(MATRIX_FLOAT *a_on, MATRIX_FLOAT *a_off, MATRIX_FLOAT *a_d, std::vector<float> *a_d_tot, std::vector<float> *a_on_tot, std::vector<float> *a_off_tot, std::vector<float> *affinity, MATRIX_INT *L_L, std::vector<float> *X_CM_x, std::vector<float> *X_CM_y, std::vector<float> *X_CM_theta, std::vector<float> *X_l_x, std::vector<float> *X_l_y, std::vector<float> *X_r_x, std::vector<float> *X_r_y, std::vector<int> *L_to_R, std::vector<int> T_l, std::vector<int> T_r)
{
	for(int i=0; i<N_L; i++)
    {
        (a_d->at(i)).clear();
        (a_on->at(i)).clear();
        (a_off->at(i)).clear();
    }

    a_d_tot->clear();
    a_on_tot->clear();
    a_off_tot->clear();
    affinity->clear();

	Initialize_Affinity(a_d, L_L, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_d, 0, 1);
	Initialize_Affinity(a_on, L_L, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_on, 1, 1);
    Initialize_Affinity(a_off, L_L, X_CM_x, X_CM_y, X_CM_theta, X_l_x, X_l_y, X_r_x, X_r_y, L_to_R, T_l, T_r, k_off, 2, 2);

    *a_d_tot = Get_a_X_tot(a_d);
	*a_on_tot = Get_a_X_tot(a_on);
	*a_off_tot = Get_a_X_tot(a_off);

	*affinity = Get_a_tot(a_d_tot, a_on_tot, a_off_tot);


}




int main()
{
	clock_t startTime = clock();
    std::ofstream trajectory;
    trajectory.open ("trajectory.txt");

    trajectory << "X_CM,Y_CM,THETA_CM,N_P" << '\n';

    std::ofstream trajectory_2;
    trajectory_2.open ("trajectory_2.txt");

    std::ofstream receptor;
    receptor.open ("recetpor.txt");

    std::ofstream ligand;
    ligand.open ("ligand.txt");

    std::ofstream run_info;
    run_info.open ("run_info.txt");

    std::ofstream Liste_Ligand;
    Liste_Ligand.open ("L_L.txt");

    /*std::ofstream Liste_Receptor;
    Liste_Receptor.open ("L_R.txt");*/

    std::ofstream Affinity;
    Affinity.open ("Affinity.txt");
  
	
    srand (time(NULL));
    long idum = -1*(rand() % 888888888 + 111111111); //-18225; //-111113882;

    run_info << "idum: " << idum << '\n';
    run_info << "k_d: " << k_d << "  k_off: " << k_off << "  k_on: " << k_on << "  DT: " << DT << "  N_BRWN: " << N_BRWN << " NBrow: " << NBrow << "  Taille: " << Taille_Systeme << "  N_R: " << N_R << "  R_Verlet: " << R_Verlet << '\n';
    
    std::cout << "seed du GNPA: " << idum << '\n';
    
  	int pct = 0;

	std::vector<float> a_d_tot;
	std::vector<float> a_on_tot;
	std::vector<float> a_off_tot;

	std::vector<float> affinity;

	MATRIX_INT v_L;
	std::vector<int> ll_L;
	std::vector<int> bl_L;

	MATRIX_INT v_R;
	std::vector<int> ll_R;
	std::vector<int> bl_R;

	std::vector<int> hoc_L;
	std::vector<int> hoc_R;

	float R_Cell;
	int N_x;
	int N_y;
	int N_Cell;
	int nc;
	int ix;
	int iy;

	MATRIX_INT Neigh;

	std::vector<float> X_l_x;  // dans le repÃ¨re du CM
	std::vector<float> X_l_y;
	
	std::vector<float> X_CM_x; // dans le repÃ¨re du laboratoire
	std::vector<float> X_CM_y;
	std::vector<float> X_CM_theta; 

	std::vector<float> X_r_x; // dans le repÃ¨re du laboratoire
	std::vector<float> X_r_y; 

	std::vector<int> T_l;
	std::vector<int> T_r;
	std::vector<int> L_to_R;
	
	float time_react;
	float time_tot;
	std::vector<float> time_vect;
	time_vect.push_back(0);
	
	int N_ponts;
    std::vector<int> ponts_vect;
    ponts_vect.push_back(0.0);
    int count;

	initialize_Virus(&X_l_x, &X_l_y, &T_l, &L_to_R, &idum);


	ligand << "X_l_x = [";
	for(int i = 0; i<X_l_x.size();i++)
	{
		ligand << X_l_x[i] << ',';
	}
	ligand << "]" << '\n' << "X_l_y =[";
	for(int i = 0; i<X_l_y.size();i++)
	{
		ligand << X_l_y[i] << ',';
	}	
	ligand << "]";

	
	MATRIX_INT L_L; 
	L_L.resize(N_L, std::vector<int>(0));
	MATRIX_INT L_R;
	L_R.resize(N_R, std::vector<int>(0));
	MATRIX_INT L_L_old; 
	L_L_old.resize(N_L, std::vector<int>(0));

	MATRIX_FLOAT a_d;
	a_d.resize(N_L, std::vector<float>(0));
	MATRIX_FLOAT a_on;
	a_on.resize(N_L, std::vector<float>(0));
	MATRIX_FLOAT a_off;
	a_off.resize(N_L, std::vector<float>(0));

    X_CM_x.push_back(0.0);
    X_CM_y.push_back(0.0);
    X_CM_theta.push_back(0.0);

    initialize_SA(&X_r_x, &X_r_y, &T_r, &idum);

	receptor << "X_r_x = [";
	for(int i = 0; i<X_r_x.size();i++)
	{
		receptor << X_r_x[i] << ',';
	}
	receptor << "]" << '\n' << "X_r_y =[";
	for(int i = 0; i<X_r_y.size();i++)
	{
		receptor << X_r_y[i] << ',';
	}	
	receptor << "]";


	N_x = int(Taille_Systeme/R_Verlet); //D->R_Verlet
	N_y = N_x;
	N_Cell = N_x*N_y;

	R_Cell = Taille_Systeme/N_x;
	
	//std::cout << R_Cell << '\n';
	//std::cout << "test 1" << '\n';
	//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

	Neigh = Set_Neigh(N_Cell, N_x, N_y);
    
	v_L = New_List(R_Cell, N_Cell, N_L, N_x, N_y, &X_l_x, &X_l_y);
	ll_L = v_L[0];
	bl_L = v_L[1];
	hoc_L = v_L[2];

	v_R = New_List(R_Cell, N_Cell, N_R, N_x, N_y, &X_r_x, &X_r_y);
	ll_R = v_R[0];
	bl_R = v_R[1];
	hoc_R = v_R[2];

	//std::cout << "size ll_L: " << ll_L.size() << " size bl_L: " << bl_L.size() << " size hoc_L: " << hoc_L.size();


	//std::cout << "test 2" << '\n';
	//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

	Initialize_L_L(&L_L, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, L_to_R, N_x, N_y, R_Cell, 0);
	Initialize_L_R(&L_R, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_L, ll_L, T_r, T_l, L_to_R, N_x, N_y, R_Cell,0);

	//std::cout << "test 3" << '\n';
	//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

	Initialize_Affinity(&a_d, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_d, 0, 1);
	Initialize_Affinity(&a_on, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_on, 1, 1);
    Initialize_Affinity(&a_off, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r, k_off, 2, 2);

	//std::cout << "test 4" << '\n';
	//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

    a_d_tot = Get_a_X_tot(&a_d);
	a_on_tot = Get_a_X_tot(&a_on);
	a_off_tot = Get_a_X_tot(&a_off);
	
	std::cout << '\n';

	affinity = Get_a_tot(&a_d_tot, &a_on_tot, &a_off_tot);

  	//std::vector<float> X_init;
  	std::vector<float> displacement;
  	std::vector<float> x;
  	std::vector<float> y;
  	std::vector<float> x_init;
  	std::vector<float> y_init;
  	for(int w=0; w<N_L;w++)
  	{
  		x_init.push_back(0.0);
  		y_init.push_back(0.0);
  		x.push_back(0.0);
  		y.push_back(0.0);
  		displacement.push_back(0.0);
  	}
  	float X_CM_x_init = 0.0;
  	float X_CM_y_init = 0.0;
  	float X_CM_theta_init = 0.0;
  	std::vector<int> check;
  	float theta;

  	for(int w=0; w < N_L; w++)
	{
		theta = X_CM_theta_init;
		x_init[w] = X_l_x[w]*cos(theta) - X_l_y[w]*sin(theta) + X_CM_x_init;
	    y_init[w] = X_l_x[w]*sin(theta) + X_l_y[w]*cos(theta) + X_CM_y_init;

		//X_init[w] = sqrt(pow(x,2) + pow(y,2));
	}




/*
	Liste_Ligand << '\n' << "L_L: " << '\n';
    for(int i = 0; i < N_L; i++)
    {
        Liste_Ligand << "Ligne de Ligand " << i+1 << ": ";
        for(int k = 0; k < L_L[i].size(); k++)
        {
            Liste_Ligand << L_L[i][k] << ' ';
        }
        Liste_Ligand << '\n';
    }
    
    Liste_Ligand << '\n';*/
	

	/*Liste_Receptor << '\n' << "L_R: " << '\n';
    for(int i = 0; i < N_R; i++)
    {
    	Liste_Receptor << "Ligne de Recepteur " << i+1 << ": ";
        for(int k = 0; k < L_R[i].size(); k++)
        {
        Liste_Receptor << L_R[i][k] << ' ';
        }
    	Liste_Receptor << '\n';
    }
    
	Liste_Receptor << '\n';*/


	/*Affinity << '\n' << "a_on: " << '\n';
    for(int i = 0; i < N_L; i++)
    {
        Affinity << "Ligne de Ligand " << i+1 << ": ";
        for(int k = 0; k < a_on[i].size(); k++)
        {
            Affinity << a_on[i][k] << ' ';
        }
        Affinity << '\n';
    }
    
    Affinity << '\n';


    Affinity << '\n' << "a_off: " << '\n';
    for(int i = 0; i < N_L; i++)
    {
        Affinity << "Ligne de Ligand " << i+1 << ": ";
        for(int k = 0; k < a_off[i].size(); k++)
        {
            Affinity << a_off[i][k] << ' ';
        }
        Affinity << '\n';
    }
    
    Affinity << '\n';


    Affinity << '\n' << "a_d: " << '\n';
    for(int i = 0; i < N_L; i++)
    {
        Affinity << "Ligne de Ligand " << i+1 << ": ";
        for(int k = 0; k < a_d[i].size(); k++)
        {
            Affinity << a_d[i][k] << ' ';
        }
        Affinity << '\n';
    }
    
    Affinity << '\n';*/


    for (int i = 0; i < N_BRWN; ++i)
	{	
		//Liste_Ligand << '\n' << '\n' << "STEP " << i <<'\n';
		/*Liste_Receptor << '\n' << '\n' << "STEP " << i <<'\n';*/	
		//trajectory << '\n' << '\n' << "STEP " << i << '\n';
		//Affinity << '\n' << '\n' << "STEP " << i << '\n';

		time_tot = 0.0;
		count = 0;
		N_ponts = 0;		

		float DT1;

		DT1=DT;

		 if(i==0)
		{
		    DT1=1.0;
		}

		//std::cout << "test 5" << '\n';
		//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';
		/*if(hoc_L.at(3452)==157)
			std::cout << '\n' << "ici c'est déjà faux 1 (étape brwn " << i << ')' << '\n';
		*/
		while(time_tot < DT1)
    	{
    		/*std::cout << "DT1: " << DT1 << '\n';
    		std::cout << time_tot << '\n';*/
    		int r = Choose(affinity, &idum); //choix du type de rÃ©action
    	    time_react = Get_time_reaction(affinity[0], &idum);
    	    time_tot += time_react;
    	    if(time_tot < DT1)
    	    {
           		//std::cout << "test " << time_tot << '\n';
        		//std::cout << r << ' ';
        	    Make_reaction(&a_d, &a_on, &a_off, a_d_tot, a_on_tot, a_off_tot, L_L, L_R, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &T_r, &L_to_R, r, &idum, i);
        	    
                a_d_tot.clear();
                a_on_tot.clear();
                a_off_tot.clear();
                affinity.clear();
        	    
        	    a_d_tot = Get_a_X_tot(&a_d);
                a_on_tot = Get_a_X_tot(&a_on);
            	a_off_tot = Get_a_X_tot(&a_off);
            	affinity = Get_a_tot(&a_d_tot, &a_on_tot, &a_off_tot);

        	    count += 1;
    	    }
    	}

    	if(i==0)
    	{
   			run_info << '\n' << "thermalisation terminée en ";
			run_info << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
    	}

    	if(i==0)
    	{
   			std::cout << '\n' << "thermalisation terminée en ";
			std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " secondes." << '\n';
    	}


    	//std::cout << '\n' << '\n' << "Nombre de réactions: " << count << '\n' << '\n';

		/*if(hoc_L.at(3452)==157)
			std::cout << '\n' << "ici c'est déjà faux 2 (étape brwn " << i << ')' << '\n';
     	
     	if(i==0){

		}*/
		//std::cout << "test 6" << '\n';		
		//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

		for(int iBrow=0; iBrow < NBrow; iBrow++)
		{
		   diffusion(&X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &T_l, &L_to_R, &idum, DT/NBrow);
		}

		//std::cout << "test 7" << '\n';		
		//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

		//X_CM_current = sqrt(pow(X_CM_x[X_CM_x.size() - 1] + X_CM_y[X_CM_x.size() - 1] ,2));
		//std::cout << pow((X_init[w] - X_current[w]),2) << " " << pow(((D - R_Cell)*0.99),2) << '\n';

		/*if(hoc_L.at(3452)==157)
			std::cout << '\n' << "ici c'est déjà faux 3 (étape brwn " << i << ')' << '\n';
		*/
		check.clear();
		//std::cout << check.size() << '\n';


		for(int w=0; w < N_L; w++)
		{
			//Liste_Ligand << sqrt(pow(x_init[w],2)+pow(y_init[w],2)) << " " <<'\n';

			theta = X_CM_theta[X_CM_theta.size() - 1];
			x[w] = X_l_x[w]*cos(theta) - X_l_y[w]*sin(theta) + X_CM_x[X_CM_x.size() -1];
		    y[w] = X_l_x[w]*sin(theta) + X_l_y[w]*cos(theta) + X_CM_y[X_CM_y.size() -1];
			displacement[w] = pow(x[w]-x_init[w],2) + pow(y[w]-y_init[w],2);

			//Liste_Ligand << displacement[w] << " ";
			//Liste_Ligand << '\n' << "Déplacement % X_init: " << sqrt(displacement[w]) << " --- " << sqrt(pow(((D - R_Cell)*0.99),2)) << '\n';

			if( displacement[w] > pow(((D - R_Cell)*0.99),2) )
				check.push_back(w + 1);
		}

		//std::cout << "test 8" << '\n';		
		//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

		//std::cout << "check size after calculation: " << check.size() << '\n' << '\n';
		/*if(hoc_L.at(3452)==157)
			std::cout << '\n' << "ici ça devient faux (étape brwn " << i << ')' << '\n';
		*/
		
		if(check.size() != 0) //pow((X_init - X_current),2) >= pow(((D - R_Cell)*0.99),2) )
	    {
    		/*clock_t startTime3 = clock();
	    	std::cout << "test 9a" << '\n';		
			std::cout << double( clock() - startTime3 ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';*/

	    	/*std::cout << "check = [";
	    	for (int p = 0; p < check.size(); p++)
	    	{
	    		std::cout << check[p] << " ";
	    	}
	    	std::cout << "]" << '\n';*/


			/*std::cout << "size ll_L: " << ll_L.size() << " size bl_L: " << bl_L.size() << " size hoc_L: " << hoc_L.size();

			if(hoc_L.at(3452)==157)
				std::cout << '\n' << "ici ça devient faux (upgrade ligand " << ')' << '\n';
    		std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';

			for (int w = 0; w < ll_L.size() + 1; w++)
			{
				std::cout << ll_L[w] << " ";
			}
			std::cout << "]" << '\n' << "[";
			for (int w = 0; w < hoc_L.size() + 1; w++)
			{
				std::cout << hoc_L[w] << " ";
			}
			std::cout << " ]" << '\n';*/

	    	Update_Cell_L(&ll_L, &bl_L, &hoc_L, &X_l_x, &X_l_y, &X_CM_x, &X_CM_y, &X_CM_theta, &check, N_x, N_y, R_Cell, X_CM_x_init, X_CM_y_init, X_CM_theta_init);

	    	/*std::cout << "test 10a" << '\n';
	    	std::cout << check.size() << '\n';		
			std::cout << double( clock() - startTime3 ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';*/

    		Update_L_L(&L_L, &L_R, &L_L_old, &check, &X_l_x, &X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, Neigh, hoc_R, ll_R, T_l, L_to_R, N_x, N_y, R_Cell);

	    	/*std::cout << "test 10b" << '\n';
	    	std::cout << check.size() << '\n';		
			std::cout << double( clock() - startTime3 ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';*/

    		Update_L_R(&L_L_old, &L_L, &L_R, &check);

 	    	/*std::cout << "test 10c" << '\n';
	    	std::cout << check.size() << '\n';		
			std::cout << double( clock() - startTime3 ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';*/

			for(int w=0; w < N_L + 1; w++)
			{	
				for(int v=0; v < check.size(); v++)
				{
					if(w + 1 == check[v])
					{	
						x_init[w] = x[w];
						y_init[w] = y[w];	
					}
				}	
	    	}

	    	X_CM_x_init = X_CM_x[X_CM_x.size() - 1];
	    	X_CM_y_init = X_CM_y[X_CM_y.size() - 1];
	    	X_CM_theta_init = X_CM_theta[X_CM_theta.size() - 1];

	    	//std::cout << "test 11a" << '\n';		
			//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';
	    	//Liste_Ligand << '\n' << "on update la liste de Verlet ici ! " << '\n';
	    	/*Liste_Receptor << '\n' << "on update la liste de Verlet ici ! " << '\n';*/
	    	/*Liste_Ligand << "check = [";
	    	for (int p = 0; p < check.size(); p++)
	    	{
	    		Liste_Ligand << check[p] << " ";
	    	}
	    	Liste_Ligand << "]" << '\n';*/
	    	//trajectory << '\n' << "on update la liste de Verlet ici ! " << '\n';
	    }


    	if( (X_CM_x.at(X_CM_x.size() - 1) != X_CM_x.at(X_CM_x.size() - 1-NBrow)) && (X_CM_y.at(X_CM_y.size() - 1) != X_CM_y.at(X_CM_y.size() - 1-NBrow)) && (X_CM_theta.at(X_CM_theta.size() - 1) != X_CM_theta.at(X_CM_theta.size() - 1-NBrow)) )
    	{	
	    	//std::cout << "test 9b" << '\n';		
			//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';
			//clock_t startTime2 = clock();

			Update_Affinity(&a_on, &a_off, &a_d, &a_d_tot, &a_on_tot, &a_off_tot, &affinity, &L_L, &X_CM_x, &X_CM_y, &X_CM_theta, &X_l_x, &X_l_y, &X_r_x, &X_r_y, &L_to_R, T_l, T_r);
			//Update_Lists(&L_L, &L_R, &a_d, &a_on, &a_off, &a_d_tot, &a_on_tot, &a_off_tot, &affinity, Neigh, hoc_R, hoc_L, ll_L, ll_R, X_l_x, X_l_y, X_r_x, X_r_y, X_CM_x, X_CM_y, X_CM_theta, N_x, N_y, R_Cell, N_L, T_l, T_r, &L_to_R, i);

	    	//std::cout << "test 10b" << '\n';		
			//std::cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';
 	    	//std::cout << "test 11" << '\n';
	    	//std::cout << check.size() << '\n';		
			//std::cout << double( clock() - startTime2 ) / (double)CLOCKS_PER_SEC<< " seconds." << '\n';
    	}


		for(int k = 0; k<T_l.size(); k++) //calcul du nombre de ponts
		{
		    if(T_l[k] == 2)
		    {
		        N_ponts += 1;
		    }
		}
		
		ponts_vect.push_back(N_ponts);		
    

    //PRINTS
	    if(i % 100 == 0)
	    {

			/*Liste_Ligand << '\n' << "L_L: " << '\n';
		    for(int i = 0; i < N_L; i++)
		    {
		        Liste_Ligand << "Ligne de Ligand " << i+1 << ": ";
		        for(int k = 0; k < L_L[i].size(); k++)
		        {
		            Liste_Ligand << L_L[i][k] << ' ';
		        }
		        Liste_Ligand << '\n';
		    }
		    
		    Liste_Ligand << '\n';*/

			/*Liste_Receptor << '\n' << "L_R: " << '\n';
		    for(int i = 0; i < N_R; i++)
		    {
		    	if(L_R[i][0] != 0)
		    	{
			    	Liste_Receptor << "Ligne de Recepteur " << i+1 << ": ";
			        for(int k = 0; k < L_R[i].size(); k++)
			        {
			        	Liste_Receptor << L_R[i][k] << ' ';
			        }
			    	Liste_Receptor << '\n';
			    }
		    }
		    
			Liste_Receptor << '\n';*/

			/*Affinity << '\n' << "a_on: " << '\n';
		    for(int i = 0; i < N_L; i++)
		    {
		        Affinity << "Ligne de Ligand " << i+1 << ": ";
		        for(int k = 0; k < a_on[i].size(); k++)
		        {
		            Affinity << a_on[i][k] << ' ';
		        }
		        Affinity << '\n';
		    }
		    
		    Affinity << '\n';


		    Affinity << '\n' << "a_off: " << '\n';
		    for(int i = 0; i < N_L; i++)
		    {
		        Affinity << "Ligne de Ligand " << i+1 << ": ";
		        for(int k = 0; k < a_off[i].size(); k++)
		        {
		            Affinity << a_off[i][k] << ' ';
		        }
		        Affinity << '\n';
		    }
		    
		    Affinity << '\n';


		    Affinity << '\n' << "a_d: " << '\n';
		    for(int i = 0; i < N_L; i++)
		    {
		        Affinity << "Ligne de Ligand " << i+1 << ": ";
		        for(int k = 0; k < a_d[i].size(); k++)
		        {
		            Affinity << a_d[i][k] << ' ';
		        }
		        Affinity << '\n';
		    }
		    
		    Affinity << '\n';*/

	    	trajectory << X_CM_x.at(X_CM_x.size() - 1) << "," << X_CM_y.at(X_CM_y.size() - 1) << "," << X_CM_theta.at(X_CM_theta.size() - 1) << "," << ponts_vect.at(ponts_vect.size() - 1) << '\n';
	      	trajectory.flush();

	      	/*trajectory_2 << '[';
	      	for(int k=0; k<T_r.size(); k++)
	      	{
	      		if(T_r[k] == 0)
	      		{
	      			trajectory_2 << k << ',' ;
	      		}
	      	}
	      	trajectory_2 << ']' << ',';
	      	trajectory_2.flush();*/	     	 	
	    }		

	    if(i % ((N_BRWN-1)/100) == 0)
	    {
	    	std::cout << pct << " %" << '\n';
	    	pct += 1;
	    }
		
	}	

   run_info << "run terminé en " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds.";
   run_info.flush();


   trajectory.close();
   trajectory_2.close();
   receptor.close();
   ligand.close();
   run_info.close();
   Liste_Ligand.close();
   Affinity.close();
     
    return 0;
}
