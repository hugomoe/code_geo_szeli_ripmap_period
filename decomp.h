#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "affine.h"
#include "ripmap.h"

//Pour les homographie
// 0 1 2
// 3 4 5
// 6 7 8

//Pour les affinité
// 0 1 2
// 3 4 5

//Et pour dH
// lambda  N


void decomp(double H[9],double A[6],double H0[9],double B[6]){
/**
  * @param
  *     H,A,H0,B : H = A H0 B
  */
    double a=H[0], b=H[1], p=H[2], c=H[3], d=H[4], q=H[5], r=H[6], s=H[7], t=H[8];
    double t0, t1;
    //on suppose r,s != 0,0



    /*
    techniquement, il existe une infinité de couples (t0,t1) tels que
    TH soit géométriquement décomposable (toute une droite affine en vérité)
    mais comme on ne fait que des rotations, le résultat et les étapes intermédiaires
    sont toutes de même difficulté, quelles que soient les valeurs choisies de t0,t1.
    On prend donc l'un des deux nul, arbitrairement
    */
    if(fabs(a*s-b*r)<fabs(c*s-d*r)){
        t0 = 0.;
        t1 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(c*s-d*r);
    }else{
        t0 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(a*s-b*r);;
        t1 = 0.;
    }



    //on change les notations en plein milieu : a,b,p,... désigneront les coef de TH
        //T=1,0,t0, 0,1,t1, 0,0,1
        //en gros, on translate les coefficients
    a += t0*r;
    b += t0*s;
    p += t0*t;
    c += t1*r;
    d += t1*s;
    q += t1*t;







    double Nphi = sqrt(pow(r,2.)+pow(s,2.));
    B[0] = r/Nphi;
	B[1] = s/Nphi;
	B[2] = 0.;
	B[3] = -s/Nphi;
	B[4] = r/Nphi;
	B[5] = 0.;



    double Npsi = sqrt(pow(a*s-b*r,2.)+pow(c*s-d*r,2.));
	A[0] = (c*s-d*r)/Npsi;
	A[1] = (a*s-b*r)/Npsi;
	A[2] = -t0;
	A[3] = -(a*s-b*r)/Npsi;
	A[4] = (c*s-d*r)/Npsi;
    A[5] = -t1;

    H0[0] = -(a*d-b*c)*Nphi/Npsi;
    H0[1] = 0.;
    H0[2] = (p*(c*s-d*r)-q*(a*s-b*r))/Npsi;
    H0[3] = 0.;
    H0[4] = -Npsi/Nphi;
    H0[5] = (p*(a*s-b*r)+q*(c*s-d*r))/Npsi;
    H0[6] = Nphi;
    H0[7] = 0.;
    H0[8] = t;

//printf("\nA :\n%f %f %f\n%f %f %f\n",A[0],A[1],A[2],A[3],A[4],A[5]);
//printf("H0 :\n%f %f %f\n%f %f %f\n%f %f %f\n",H0[0],H0[1],H0[2],H0[3],H0[4],H0[5],H0[6],H0[7],H0[8]);
//printf("B :\n%f %f %f\n%f %f %f\n",B[0],B[1],B[2],B[3],B[4],B[5]);
}



int apply_homography_ripmap(float *img,float *img_f,int w,int h,int w_f,int h_f,int mu2,int nu2,int mu3,int nu3,double H0[9],double A[6],int w1,int h1){

	double H[3][3];
	H[0][0]=H0[0];
	H[0][1]=H0[1];
	H[0][2]=H0[2];
	H[1][0]=H0[3];
	H[1][1]=H0[4];
	H[1][2]=H0[5];
	H[2][0]=H0[6];
	H[2][1]=H0[7];
	H[2][2]=H0[8];

	double D[12];
	precal_D(H,D);
//on va construire le ripmap de img.
	float *ripmap=malloc(12*w*h*sizeof(float));
	int logw = (int) log2(w);
	build_ripmap(img,ripmap,w,h,3,logw);

	for (int j = 0; j < h_f; j++)
	for (int i = 0; i < w_f; i++)
	{
		double p[2] = {i+mu3, j+nu3};
		apply_homography(p, H, p);
		p[0]=p[0]-mu2;
		p[1]=p[1]-nu2;
		p[0] = (p[0] - 0.5) * w / (w - 1.0);
		p[1] = (p[1] - 0.5) * h / (h - 1.0);

		double z[2]={i+mu3,j+nu3};
		double a;
		a = pow(D[9]*z[0]+D[10]*z[1]+D[11],2);
		double dudx,dudy,dvdx,dvdy;
		dudx = (D[1]*z[1]+D[2])/a;
		dudy = (D[5]*z[0]+D[6])/a;
		dvdx = (D[3]*z[1]+D[4])/a;
		dvdy = (D[7]*z[0]+D[8])/a;

		double det = absd( dudy*dvdx - dvdy*dudx );
		double d[2];
		d[0] = absd(dudx) + absd(dudy);
		d[1] = absd(dvdx) + absd(dvdy);

		//on doit faire une copie de p, sinon il ne veut pas dŽcaler...
		double coo[2];
		coo[0]=p[0]; coo[1]=p[1];
		/*
                //essai de périodisation : en fait c'est normal que ça n'a pas fonctionné.
        double uv[2]  ;
        uv[0]=(double)coo[0] ;
        uv[1]=(double)coo[1] ;
        periodisation(uv,w,h,w1,h1,A);

        coo[0] = uv[0];
        coo[1] = uv[1];
                //fin essai de périodisation
        */

		//On decale l'origine du rectangle pour que le parallŽlogramme soit ˆ l'intŽrieur
		if(dudx<0){coo[0] += dudx;}
		if(dudy<0){coo[0] += dudy;}
		if(dvdx<0){coo[1] += dvdx;}
		if(dvdy<0){coo[1] += dvdy;}

		for (int l = 0; l < 3; l++)
		{

                //essai de périodisation : en fait c'est normal que ça n'a pas fonctionné.
        double uv[2]  ;
        uv[0]=(double)coo[0] ;
        uv[1]=(double)coo[1] ;
        periodisation(uv,w,h,w1,h1,A);

        coo[0] = uv[0];
        coo[1] = uv[1];
                //fin essai de périodisation
        int idx = l + 3 * (w_f * j + i);
        img_f[idx]=ripmap_interpolation_at(ripmap, w, h, coo[0], coo[1], l, d,A,w1,h1);
		}
	}
	return 0;
}



void apply_homo_final(float *img,float *img_f,int w,int h,int w_f,int h_f,double H[3][3]){
/**
  * @param
  *     img : image d'entrée
  *     img_f : image de sortie
  *     w,h : dimension de l'image d'entrée
  *     H : matrice de l'homographie img_f(x,y)=img(H(x,y))
  */
//LES BLOCS DOCUMENTATION PLUS BAS NE SERVIRONT QUE QUAND TOUT FONCTIONNERA ET SERA MODIFIÉ
//printf("w x h : %d x %d\n",w,h);printf("w_f x h_f : %d x %d\n",w_f,h_f);
	double *a = *H;
	//faire un cas où a[8]=0 (avec une translation d'un pixel
        //i.e. agrandir l'image d'un pixel puis la rogner)
    //dans le cas de la décomposimon, on ne devrait pas avoir à diviser par a[8]
        //a condition que l'homo particuliere accepte le cas a[8]=0
	for(int i=0;i<8;i++){a[i]=a[i]/a[8];}
	a[8]=1.;
printf("------------------- debut de apply_homo_final -------------------");
//code impression
printf("\nH :\n");
printf("%f %f %f\n",H[0][0],H[0][1],H[0][2]);
printf("%f %f %f\n",H[1][0],H[1][1],H[1][2]);
printf("%f %f %f\n",H[2][0],H[2][1],H[2][2]);
//fin code impression
	if(a[6]==0 && a[7]==0){
    //cas d'une affinité
		double A[6];
		for(int i=0;i<6;i++){A[i]=a[i];}
		apply_affinite(img,img_f,w,h,w_f,h_f,A);
		//il faudrait rogner, si jamais A translate
	}
	else{
		double A[6];
		double H0[9];
		double B[6];

		decomp(a,A,H0,B);









        //calcul de w1,h1,mu1,nu1
		double x4[4] = {0,w_f,w_f,0};
        double y4[4] = {0,0,h_f,h_f};
        //x1,y1 = H(x4,y4)
        double x1[4];
        double y1[4];
        for(int k=0;k<4;k++){
            //si (x4[k]*a[6]+y4[k]*a[7]+a[8]), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!prévoir un cas pour simuler l'infini
            x1[k]=(x4[k]*a[0]+y4[k]*a[1]+a[2])/(x4[k]*a[6]+y4[k]*a[7]+a[8]);
            y1[k]=(x4[k]*a[3]+y4[k]*a[4]+a[5])/(x4[k]*a[6]+y4[k]*a[7]+a[8]);
        }
        double temp_min_x1 = x1[0];
        double temp_min_y1 = y1[0];
        double temp_max_x1 = x1[0];
        double temp_max_y1 = y1[0];
        for(int k=1;k<4;k++){
            if(x1[k]<temp_min_x1){temp_min_x1 = x1[k];}
            if(y1[k]<temp_min_y1){temp_min_y1 = y1[k];}
            if(x1[k]>temp_max_x1){temp_max_x1 = x1[k];}
            if(y1[k]>temp_max_y1){temp_max_y1 = y1[k];}
        }
        int min_x1 = floor(temp_min_x1);
        int min_y1 = floor(temp_min_y1);
        int max_x1 = ceil(temp_max_x1);
        int max_y1 = ceil(temp_max_y1);
        int w1 = max_x1-min_x1, h1 = max_y1-min_y1, mu1 = min_x1, nu1 = min_y1;
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!à faire
        //tronquer l'image si on ne périodise pas
        //périodiser puis tronquer la périodisée en ne gardant au plus qu'une période si on périodise
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pour l'instant, on oublie ça
        w1=w;
        h1=h;
        mu1=0;
        nu1=0;
        //x1,y1 sont tjrs les sommets du quadrilatère utile
        //si ils dépassaient de l'image, peut-être vaut-il mieux leur donner maintenant la valeur "coins de img1"
		x1[0]=0; x1[1]=w; x1[2]=w; x1[3]=0;
		y1[0]=0; y1[1]=0; y1[2]=h; y1[3]=h;



		//x2,y2 = AA(x1,y1) où AA = A^(-1)
		double det_A = A[0]*A[4]-A[3]*A[1];
        double AA[6] = {
            A[4]/det_A,
            -A[1]/det_A,
            (A[1]*A[5]-A[4]*A[2])/det_A,
            -A[3]/det_A,
            A[0]/det_A,
            (-A[0]*A[5]+A[3]*A[2])/det_A
        };

        double x2[4];
        double y2[4];
        for(int k=0;k<4;k++){
            x2[k]=x1[k]*AA[0]+y1[k]*AA[1]+AA[2];
            y2[k]=x1[k]*AA[3]+y1[k]*AA[4]+AA[5];
        }
        double temp_min_x2 = x2[0];
        double temp_min_y2 = y2[0];
        double temp_max_x2 = x2[0];
        double temp_max_y2 = y2[0];
        for(int k=1;k<4;k++){
            if(x2[k]<temp_min_x2){temp_min_x2 = x2[k];}
            if(y2[k]<temp_min_y2){temp_min_y2 = y2[k];}
            if(x2[k]>temp_max_x2){temp_max_x2 = x2[k];}
            if(y2[k]>temp_max_y2){temp_max_y2 = y2[k];}
        }

        int min_x2 = floor(temp_min_x2);
        int min_y2 = floor(temp_min_y2);
        int max_x2 = ceil(temp_max_x2);
        int max_y2 = ceil(temp_max_y2);
        int w2 = max_x2-min_x2, h2 = max_y2-min_y2, mu2 = min_x2, nu2 = min_y2;
        //augmenter w2, h2 pour avoir même parité que w1, h1
        if((w2-w1)%2!=0){w2++;}
        if((h2-h1)%2!=0){h2++;}
/*int ns2 = 512; //ns la new size
mu2 = mu2+(w2-ns2)/2; nu2 = nu2+(h2-ns2)/2;
w2 = ns2; h2 = ns2;*/
printf("\nw2 x h2 : %d x %d\n",w2,h2);
printf("mu2 x nu2 : %d x %d\n",mu2,nu2);



        //x3,y3 = B(x4,y4)
        double x3[4];
        double y3[4];
        for(int k=0;k<4;k++){
            x3[k]=x4[k]*B[0]+y4[k]*B[1]+B[2];
            y3[k]=x4[k]*B[3]+y4[k]*B[4]+B[5];
        }
        double temp_min_x3 = x3[0];
        double temp_min_y3 = y3[0];
        double temp_max_x3 = x3[0];
        double temp_max_y3 = y3[0];
        for(int k=1;k<4;k++){
            if(x3[k]<temp_min_x3){temp_min_x3 = x3[k];}
            if(y3[k]<temp_min_y3){temp_min_y3 = y3[k];}
            if(x3[k]>temp_max_x3){temp_max_x3 = x3[k];}
            if(y3[k]>temp_max_y3){temp_max_y3 = y3[k];}
        }
        int min_x3 = floor(temp_min_x3);
        int min_y3 = floor(temp_min_y3);
        int max_x3 = ceil(temp_max_x3);
        int max_y3 = ceil(temp_max_y3);
        int w3 = max_x3-min_x3, h3 = max_y3-min_y3, mu3 = min_x3, nu3 = min_y3;
        //augmenter w3, h3 pour avoir même parité que w_f, h_f
        if((w3-w_f)%2!=0){w3++;}
        if((h3-h_f)%2!=0){h3++;}
/*int ns3 = 512; //ns la new size
mu3 = mu3+(w3-ns3)/2; nu3 = nu3+(h3-ns3)/2;
w3 = ns3; h3 = ns3;*/
printf("w3 x h3 : %d x %d\n",w3,h3);
printf("mu3 x nu3 : %d x %d\n",mu3,nu3);



        int w4 = w_f, h4 = h_f;
        int mu4 = 0, nu4 = 0;
printf("w4 x h4 : %d x %d\n",w4,h4);
printf("mu4 x nu4 : %d x %d\n",mu4,nu4);



        //à ce stade, on a pris en compte une grande partie des translations de A et , via les mu,nu
        //on modifie donc A et B en conséquence

        A[2] = mu2*A[0] + nu2*A[1] + A[2] - mu1;
        A[5] = mu2*A[3] + nu2*A[4] + A[5] - nu1;

        B[2] = mu4*B[0] + nu4*B[1] + B[2] - mu3;
        B[5] = mu4*B[3] + nu4*B[4] + B[5] - nu3;









printf("\ntransformation de l'image en cours\n");
//        /**
        float *img1 = malloc(3*w1*h1*sizeof(float));
		for(int i=0;i<3*w1*h1;i++){img1[i]=img[i];}

        float *img2 = malloc(3*w2*h2*sizeof(float));
		apply_affinite(img1,img2,w1,h1,w2,h2,A);
		//free(img1);
printf("application de A  : done\n");

		float *img3 = malloc(3*w3*h3*sizeof(float));
		apply_homography_ripmap(img2,img3,w2,h2,w3,h3,mu2,nu2,mu3,nu3,H0,A,w1,h1);
		//free(img2);
printf("application de H0 : done\n");

		float *img4 = malloc(3*w_f*h_f*sizeof(float));
		apply_affinite(img3,img4,w3,h3,w4,h4,B);
		//free(img3);
printf("application de B  : done\n");

		for(int l=0;l<3;l++)
            for(int i=0;i<w_f;i++)
                for(int j=0;j<h_f;j++)
                    {img_f[(i+w_f*j)*3+l]=img4[(i+w4*j)*3+l];}
        //free(img4);
	}
printf("-------------------- fin de apply_homo_final --------------------\n");
}
