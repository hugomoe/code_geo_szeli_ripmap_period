//SECTION 5.5 Ripmap
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//construit le ripmap (qui est une image 4 fois plus grande, toujours de profondeur 3)
//On a un ripmap de type float*
//Les deux premiers entiers correspondent a la taille du rectangle
//on a une fonction qui donne les coordonnées dans le ripmap en fonction de la taille des rectangles


//a inclure dans main :
//	int w = e->iw;
//on va construire le ripmap de img.
//	float *ripmap=malloc(12*w*w*sizeof(float));
//	int logw = (int) log2(e->iw);
//	e->img = build_ripmap(e->img,ripmap,w,e->ih,e->pd,logw);

//int coord(int i,int j,int u,int v,int w,int l){
//	return (2*(1-1/pow(2,i))*w + 4*(1-1/pow(2,j))*w*w + u + v*2*w)*3+l;
///}


//des problèmes avec les statiques.

/*float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}*/

/*int good_modulus(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	return r;
}*/


// périodisation :


int calc_rot(double pt[2],int w2,int h2,double A[6]){
	double ptbis[2] ;
        double det_A=A[0]*A[4]-A[1]*A[3] ;
        //double AA[6] = {A[4]/det_A,-A[1]/det_A,(A[1]*A[5]-A[4]*A[2])/det_A,-A[3]/det_A,A[0]/det_A,(-A[0]*A[5]+A[3]*A[2])/det_A};
        double AA[6] = {A[0],A[1],A[2],A[3],A[4],A[5]};
        ptbis[0]=AA[0]*pt[0]+AA[1]*pt[1]+AA[2]; // terme de translation pour faire une rotation par rapport au centre de l'image
        ptbis[1]=AA[3]*pt[0]+AA[4]*pt[1]+AA[5];
        pt[0]=ptbis[0];
        pt[1]=ptbis[1];

 return 0 ;
}



int calc_vect(double vect1[2],double vect2[2],double centre[2],int w2,int h2 ,int dim1,int dim2,double A[6]){
	double pt1[2]={0,0};
        double pt2[2]={0,dim2};
        double pt3[2]={dim1,0};
        //double pt4[2]={dim1,dim2}; useless...
        calc_rot(pt1,dim1,dim2,A);
        calc_rot(pt2,dim1,dim2,A);
	calc_rot(pt3,dim1,dim2,A);
	//calc_rot(pt4,dim1,dim2,A);
        // attention : prendre dim1=w1, dim2=h1 ...
        
        
        vect1[0]=-(pt1[0]-pt2[0]) ;
        vect1[1]=-(pt1[1]-pt2[1]) ;
        vect2[0]=-(pt1[0]-pt3[0]) ;
	vect2[1]=-(pt1[1]-pt3[1]) ;

        double vect1_test[2] ={vect1[0],vect1[1]};
        double vect2_test[2] ={vect2[0],vect2[1]};

        /*prendre un base directe (à priori elle l'était déjà avant) :
        
        if ((vect1[0]*vect2[1]-vect1[1]*vect2[0])<0){
        	vect1[0] =vect2_test[0];
		vect1[1] =vect2_test[1];
		vect2[0] =vect1_test[0];         
		vect2[1] =vect1_test[1]; }       
        fin de prendre une base directe*/ 
        
        //calcul du centre :
        centre[0]=w2/2-(1.0/2.0)*vect1[0]-(1.0/2.0)*vect2[0] ;
        centre[1]=h2/2-(1.0/2.0)*vect1[1]-(1.0/2.0)*vect2[1] ;





 return 0 ;
}

int periodisation(double uv[2],int w2,int h2 ,int w1,int h1,double A[6]){

double vect1[2];
double vect2[2];
double centre[2];
double AAA[6]={A[0],A[1],-w1/2.0,A[3],A[4],-h1/2.0};
calc_vect(vect2,vect1,centre,h2,w2,h1,w1,AAA) ;

double invcentre[2]={centre[1],centre[0]};
centre[1]=invcentre[1] ;
centre[0]=invcentre[0] ;






double wx[2]={uv[1],uv[0]};
uv[1]=wx[0] ;
uv[0]=wx[1] ;



//Changement d'origine :
double wi=uv[0]-centre[0] ;
double wj=uv[1]-centre[1] ;

//Changement de base :
double det=vect1[0]*vect2[1]-vect1[1]*vect2[0] ;
double wa=(vect2[1]*wj-vect2[0]*wi)/det;
double wb=(-vect1[1]*wj+vect1[0]*wi)/det;

//reduction modulo ce qui va bien :
double ra=wa-floor(wa);
double rb=wb-floor(wb);

//changement de base inverse :
 wi=ra*vect1[1]+rb*vect2[1] ;
 wj=ra*vect1[0]+rb*vect2[0] ;
//retour à l’origine :
  wx[0]=wi+centre[0] ;
  wx[1]=wj+centre[1] ;



uv[1]=wx[1] ;
uv[0]=wx[0] ;
	return 0 ;
}



int coord_period(int i,int j,int u,int v,int w,int l,double A[6],int w1,int h1){
	/*double uv[2]  ;
        uv[0]=(double)u ;
        uv[1]=(double)v ;
        periodisation(uv,w/pow(2,i),w/pow(2,j),w1/pow(2,i),w1/pow(2,j),A);
        
        int x = good_modulus(uv[0],w/pow(2,i));
	int y = good_modulus(uv[1],w/pow(2,j));*/

	int x = good_modulus(u,w/pow(2,i));
	int y = good_modulus(v,w/pow(2,j));

	return (2*(1-1/pow(2,i))*w + 4*(1-1/pow(2,j))*w*w + x + y*2*w)*3+l;
}




//fin périodisation


int coord(int i,int j,int u,int v,int w,int l){
	int x = good_modulus(u,w/pow(2,i));
	int y = good_modulus(v,w/pow(2,j));
	return (2*(1-1/pow(2,i))*w + 4*(1-1/pow(2,j))*w*w + x + y*2*w)*3+l;
}


//On construit le ripmap. On a ici fait un flitre gaussien dans la direction de la compression

#define TAPSR 5     //nombre de coefficients non nuls du filtre gaussien (doit être impaire)
#define SIG 0.6    //paramètre du filtre gaussien


//deux fonctions de flitrage verticale et horizontale

float gaussian_filter_h(float *r,int i,int j,int u,int v,int w,int l,float *g){
	float a,b;
	a = b = 0;
	for(int f=0;f<TAPSR;f++){a += r[coord(i-1,j,2*u+f-(TAPSR-1)/2,v,w,l)]*g[f];}
	for(int f=0;f<TAPSR;f++){b += r[coord(i-1,j,2*u+1+f-(TAPSR-1)/2,v,w,l)]*g[f];}
	return (a+b)/2;
}

float gaussian_filter_v(float *r,int i,int j,int u,int v,int w,int l,float *g){
	float a,b;
	a = b = 0;
	for(int f=0;f<TAPSR;f++){a += r[coord(i,j-1,u,2*v+f-(TAPSR-1)/2,w,l)]*g[f];}
	for(int f=0;f<TAPSR;f++){b += r[coord(i,j-1,u,2*v+1+f-(TAPSR-1)/2,w,l)]*g[f];}
	return (a+b)/2;
}


void build_ripmap(float *img,float *r,int w,int h,int pd,int logw){
	int l,ll,i,j,u,v,w1,w2;
	
	//on déclare un filtre gaussien 1D;
	float gauss1D[TAPSR];
        float M_PI=3.141592653 ;
	for(int f=0;f<TAPSR;f++){
		gauss1D[f]=exp(-pow((TAPSR-1)/2-f,2)/(2*pow(SIG,2)))/(sqrt(2*M_PI)*SIG);
	}
	float total = 0;
	for(int f=0;f<TAPSR;f++){total = total + gauss1D[f];}
	for(int f=0;f<TAPSR;f++){gauss1D[f] = gauss1D[f]/total;}
	
	for(u=0;u<w;u++){
		for(v=0;v<h;v++){
			for(l=0;l<3;l++){
				if(l>pd-1){ll=pd-1;}else{ll=l;}
				r[coord(0,0,u,v,w,l)]=img[(u+v*w)*pd+ll];
			}
		}
	}
//on construit l'image d'origine en (0,0)	

	for(i=1;pow(2,i)<=w;i++){
		w1=w/pow(2,i);
		for(u=0;u<w1;u++){
			for(v=0;v<w;v++){
				for(l=0;l<3;l++){    
					r[coord(i,0,u,v,w,l)] = gaussian_filter_h(r,i,0,u,v,w,l,gauss1D);
				}
			}
		}
	}
//ici on moyenne simplement pour avoir tous les rectangles aplatis en largeur	

	for(i=0;pow(2,i)<=w;i++){
		w1=w/pow(2,i);
		for(j=1;pow(2,j)<=w;j++){
			w2=w/pow(2,j);
			for(u=0;u<w1;u++){
				for(v=0;v<w2;v++){
					for(l=0;l<3;l++){
						r[coord(i,j,u,v,w,l)] = gaussian_filter_v(r,i,j,u,v,w,l,gauss1D);
					}
				}
			}
		}
	}
//ici on obtient les rectangles aplatis en hauteur a partir de ceux calculés avant 
	//return r;
}


//Pour la fonction de distance, on prend dh = |du/dx|+|du/dy| et dl = |dv/dx| + |dv/dy|
//On a ainsi le plus petit rectangle qui contient l'affinité localement


// on a utilisé cela comme ref pour écrire les fonctions
// c = D[9]x+D[10]y+D[11]
// du/dx = D[1]y+D[2]/c
// du/dy = D[5]x+D[6]/c
// dv/dx = D[3]y+D[4]/c
// dv/dy = D[7]x+D[8]/c


void precal_D(double H[3][3],double *D){
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0];
	D[1]=a[0]*a[7]-a[6]*a[1];
	D[2]=a[0]*a[8]-a[6]*a[2];
	D[3]=a[3]*a[7]-a[6]*a[4];
	D[4]=a[3]*a[8]-a[6]*a[5];
	D[5]=-D[1];
	D[6]=a[1]*a[8]-a[7]*a[2];
	D[7]=-D[3];
	D[8]=a[4]*a[8]-a[7]*a[5];
	D[9]=a[6];
	D[10]=a[7];
	D[11]=a[8];
}


/*double absd(double p){if(p>0){return p;}{return -p;}} *///une valeur absolu pour les double


//un variable globale que l'on additionne a la distance
#define D_BIAS 0.
#define D_COEFF 1.

//la fonction de calcul de D. On a pris le plus petit rectangle qui encadre le parallélogramme
void cal_D(int x,int y,double *D,double *d,double *coo){
	double p[2]={x,y};
	double a;
	a = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	
//On declare les dérivé partielles	
	double dudx,dudy,dvdx,dvdy;
	dudx = (D[1]*p[1]+D[2])/a;
	dudy = (D[5]*p[0]+D[6])/a;
	dvdx = (D[3]*p[1]+D[4])/a;
	dvdy = (D[7]*p[0]+D[8])/a;
	
//On decale l'origine du rectangle pour que le parallélogramme soit à l'intérieur	
	if(dudx<0){coo[0] += dudx*D_COEFF;}
	if(dudy<0){coo[0] += dudy*D_COEFF;}
	if(dvdx<0){coo[1] += dvdx*D_COEFF;}
	if(dvdy<0){coo[1] += dvdy*D_COEFF;}
	
	d[0] = ( absd(dudx) + absd(dudy) + D_BIAS )*D_COEFF;
	d[1] = ( absd(dvdx) + absd(dvdy) + D_BIAS )*D_COEFF;
}


//l'interpolation bilinéaire
static float bilinear_ripmap(float *x, int w,
		float p, float q, int l, int d1, int d2,double A[6],int w1,int h1){	
	float pp = p - 1/2*(pow(2,d1)-1)/pow(2,d1);
	float qq = q - 1/2*(pow(2,d2)-1)/pow(2,d2);   //on decale car les rectangles représente plusieurs pixels
	int ip = floor(p);
	int iq = floor(q);
	float a = x[coord_period(d1,d2,ip,iq,w,l,A,w1,h1)];
	float b = x[coord_period(d1,d2,ip+1,iq,w,l,A,w1,h1)];
	float c = x[coord_period(d1,d2,ip,iq+1,w,l,A,w1,h1)];
	float dd = x[coord_period(d1,d2,ip+1,iq+1,w,l,A,w1,h1)];
	return evaluate_bilinear_cell(a, b, c, dd, pp-ip, qq-iq);
}


//la fonction principale. On a distingué beaucoup de cas suivant dh<=1, dh>=logw, etc...

static float ripmap_interpolation_at(float *r, int w, int h,
		float x, float y, int l,double d[2],double A[6],int w1,int h1){
		int logw = (int) log2(w);
		int dh = floor(log2(d[0])) + 1;
		int dl = floor(log2(d[1])) + 1;
		float a,b,c,dd,u,v;
		if((pow(2,dh))>w && (pow(2,dl)>w)){return r[coord(logw,logw,0,0,w,l)];}
		if(dh<=1 && dl<=1){return bilinear_ripmap(r,w,x,y,l,0,0,A,w1,h1);}
		if(pow(2,dh)>w){
			if(dl<=1){return bilinear_ripmap(r,w,0,y,l,logw,0,A,w1,h1);}
			{
				a=bilinear_ripmap(r,w,0,y/pow(2,dl-1),l,logw,dl-1,A,w1,h1);
				b=bilinear_ripmap(r,w,0,y/pow(2,dl-2),l,logw,dl-2,A,w1,h1);
				return (dl-log2(d[1]))*b + (log2(d[1])-dl+1)*a;
			}
		}
		if(pow(2,dl)>w){
			if(dh<=1){return bilinear_ripmap(r,w,x,0,l,0,logw,A,w1,h1);}
			{
				a=bilinear_ripmap(r,w,x/pow(2,dh-1),0,l,dh-1,logw,A,w1,h1);
				b=bilinear_ripmap(r,w,x/pow(2,dh-2),0,l,dh-2,logw,A,w1,h1);
				return (dh-log2(d[0]))*b + (log2(d[0])-dh+1)*a;
			}
		}
		if(dh<=1){
			a=bilinear_ripmap(r,w,x,y/pow(2,dl-1),l,0,dl-1,A,w1,h1);
			b=bilinear_ripmap(r,w,x,y/pow(2,dl-2),l,0,dl-2,A,w1,h1);
			return (dl-log2(d[1]))*b + (log2(d[1])-dl+1)*a;}
		if(dl<=1){
			a=bilinear_ripmap(r,w,x/pow(2,dh-1),y,l,dh-1,0,A,w1,h1);
			b=bilinear_ripmap(r,w,x/pow(2,dh-2),y,l,dh-2,0,A,w1,h1);
			return (dh-log2(d[0]))*b + (log2(d[0])-dh+1)*a;
		
		}
		a=bilinear_ripmap(r,w,x/pow(2,dh-1),y/pow(2,dl-1),l,dh-1,dl-1,A,w1,h1);
		b=bilinear_ripmap(r,w,x/pow(2,dh-2),y/pow(2,dl-1),l,dh-2,dl-1,A,w1,h1);
		c=bilinear_ripmap(r,w,x/pow(2,dh-1),y/pow(2,dl-2),l,dh-1,dl-2,A,w1,h1);
		dd=bilinear_ripmap(r,w,x/pow(2,dh-2),y/pow(2,dl-2),l,dh-2,dl-2,A,w1,h1);
		u=(pow(2,dh)-d[0])/pow(2,dh);
		v=(pow(2,dl)-d[1])/pow(2,dl);
		return u*v*dd + (1-u)*v*b + (1-v)*u*c + (1-u)*(1-v)*a;
	}
