#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1
#define N 2
#define epssm 0.0001

float *xp,**yp,dxsav,a,eps;
int kmax,kount,*iup,*jup,*idown,*jdown;

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_vector(v,nl,nh)
float *v;
long nh,nl;
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void derivs(x,y,dydx)
float dydx[],x,y[];
{
        dydx[1]=y[2];
        dydx[2]=(y[2]*(1-y[2])*(y[2]-a)-y[1])/eps;
}

void rkck(y,dydx,n,x,h,yout,yerr,derivs)
float dydx[],h,x,y[],yerr[],yout[];
int n;
void (*derivs)();
{
        int i;
        static float a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
                b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
                b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
                b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
                b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
                c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
                dc5 = -277.0/14336.0;
        float dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
                dc4=c4-13525.0/55296.0,dc6=c6-0.25;
        float *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

        ak2=vector(1,n);
        ak3=vector(1,n);
        ak4=vector(1,n);
        ak5=vector(1,n);
        ak6=vector(1,n);
        ytemp=vector(1,n);
        for (i=1;i<=n;i++)
                ytemp[i]=y[i]+b21*h*dydx[i];
        (*derivs)(x+a2*h,ytemp,ak2);
        for (i=1;i<=n;i++)
                ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
        (*derivs)(x+a3*h,ytemp,ak3);
        for (i=1;i<=n;i++)
                ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
        (*derivs)(x+a4*h,ytemp,ak4);
        for (i=1;i<=n;i++)
                ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
        (*derivs)(x+a5*h,ytemp,ak5);
        for (i=1;i<=n;i++)
                ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
        (*derivs)(x+a6*h,ytemp,ak6);
        for (i=1;i<=n;i++)
                yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
        for (i=1;i<=n;i++)
                yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
        free_vector(ytemp,1,n);
        free_vector(ak6,1,n);
        free_vector(ak5,1,n);
        free_vector(ak4,1,n);
        free_vector(ak3,1,n);
        free_vector(ak2,1,n);
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
float *hdid,*hnext,*x,dydx[],eps,htry,y[],yscal[];
int n;
void (*derivs)();
{
        void rkck();
        int i;
        float errmax,h,xnew,*yerr,*ytemp;

        yerr=vector(1,n);
        ytemp=vector(1,n);
        h=htry;
        for (;;) {
                rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
                errmax=0.0;
                for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
                errmax /= eps;
                if (errmax > 1.0) {
                        h=SAFETY*h*pow(errmax,PSHRNK);
                        if (h < 0.1*h) h *= 0.1;
                        xnew=(*x)+h;
                        continue;
                } else {
                        if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
                        else *hnext=5.0*h;
                        *x += (*hdid=h);
                        for (i=1;i<=n;i++) y[i]=ytemp[i];
                        break;
                }
        }
        free_vector(ytemp,1,n);
        free_vector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

#define MAXSTP 10000
#define TINY 1.0e-30

void odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
float eps,h1,hmin,x1,x2,ystart[];
int *nbad,*nok,nvar;
void (*derivs)(),(*rkqs)();
{
        int nstp,i;
        float xsav,x,hnext,hdid,h;
        float *yscal,*y,*dydx;

        yscal=vector(1,nvar);
        y=vector(1,nvar);
        dydx=vector(1,nvar);
        x=x1;
        h=SIGN(h1,x2-x1);
        *nok = (*nbad) = kount = 0;
        for (i=1;i<=nvar;i++) y[i]=ystart[i];
        if (kmax > 0) xsav=x-dxsav*2.0;
        for (nstp=1;nstp<=MAXSTP;nstp++) {
                (*derivs)(x,y,dydx);
                for (i=1;i<=nvar;i++)
                        yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
                if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
                        xp[++kount]=x;
                        for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
                        xsav=x;
                }
                if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
                (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
                if (hdid == h) ++(*nok); else ++(*nbad);
                if ((x-x2)*(x2-x1) >= 0.0) {
                        for (i=1;i<=nvar;i++) ystart[i]=y[i];
                        if (kmax) {
                                xp[++kount]=x;
                                for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
                        }
                        free_vector(dydx,1,nvar);
                        free_vector(y,1,nvar);
                        free_vector(yscal,1,nvar);
                        return;
                }
                h=hnext;
        }
}
#undef MAXSTP
#undef TINY

float gasdev(idum)
long *idum;
{
        float ran3();
        static int iset=0;
        static float gset;
        float fac,rsq,v1,v2;

        if  (iset == 0) {
                do {
                        v1=2.0*ran3(idum)-1.0;
                        v2=2.0*ran3(idum)-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                    if (mk < MZ) mk += MBIG;
                              mj=ma[ii];
                      }
                      for (k=1;k<=4;k++)
                              for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                              }
                      inext=0;
                      inextp=31;
                      *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

main()
{    FILE *fp1;
     int nbad,nok,t,lattice_size_x,i;
     float *sol,h1=0.1,hmin=0.0,x1,x2,*lap,deltax,fac;
     float *u,*v,timestep,kappa;

     fp1=fopen("spiralmodel1d","w"); 
     a=0.2;
     eps=0.05;
     deltax=0.2;
     timestep=deltax*deltax/5;
     kappa=timestep/(deltax*deltax);
     lattice_size_x=512;
     u=vector(1,lattice_size_x);
     v=vector(1,lattice_size_x);
     lap=vector(1,lattice_size_x);
     sol=vector(1,2);
     iup=ivector(1,lattice_size_x);
     idown=ivector(1,lattice_size_x);
     for (i=1;i<=lattice_size_x;i++)
      {iup[i]=i+1;
       idown[i]=i-1;}
     iup[lattice_size_x]=1;
     idown[1]=lattice_size_x;
     for (i=1;i<=lattice_size_x;i++)
      if ((i>lattice_size_x/2-5)&&(i<lattice_size_x/2-2)) 
       {u[i]=0.25;v[i]=0.0;} 
      else 
       {u[i]=0;v[i]=0;}
     /* plot initial condition first */
     for (i=1;i<=lattice_size_x;i++)
      fprintf(fp1,"%f %f %f\n",(i-lattice_size_x/2)/10.0,u[i],-v[i]);
     fac=timestep*kappa/deltax/deltax;
     x1=0;x2=timestep;
	 for (t=1;t<=1000;t++)
      {for (i=1;i<=lattice_size_x;i++)
        lap[i]=0;
       for (i=1;i<=lattice_size_x;i++)
        {sol[1]=v[i];sol[2]=u[i];
		 odeint(sol,N,x1,x2,epssm,h1,hmin,&nok,&nbad,derivs,rkqs);
		 v[i]=sol[1];u[i]=sol[2];
         if (fabs(u[i])>0.00001) {lap[i]=-2*u[i];lap[iup[i]]=u[i];lap[idown[i]]=u[i];}}
       for (i=1;i<=lattice_size_x;i++)
        if (fabs(lap[i])>0) u[i]+=fac*lap[i];}
     for (i=1;i<=lattice_size_x;i++)
      fprintf(fp1,"%f %f %f\n",(i-lattice_size_x/2)/10.0,u[i],-v[i]);
     fclose(fp1);
}  
