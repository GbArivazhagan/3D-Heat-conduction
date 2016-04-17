//Code for 3D heat conduction using flux based approach
#include<stdio.h>
#include<iostream>
#include<time.h>
#include<string.h>
using namespace std;
double ***x,***y,***z,***xx,***xy,***xz,***Qxold,***Qyold,***Qzold;
double ***qx,***qy,***qz,***T,***q,***qold,***Told,***unstead;
void printit(int);
int Imax=66; //Max grid points
int Jmax=66;
int Kmax=258;
int imax=Imax-1; //To denote the flux indices
int jmax=Jmax-1;
int kmax=Kmax-1;
int i,j,k;
int main()
{
	int ii,ij,ik,iter;//have used two indices; imax,jmax for flux and Imax and Jmax for Temperature for convenience
	double dx,dy,dz,dt,Lx,Ly,Lz; //grid properties
	
	Lx=1.0;
	Ly=1.0;
	Lz=1.0;
	iter=0; //To keep track of iteration
	imax=Imax-1; //To denote the flux indices
	jmax=Jmax-1;
	kmax=Kmax-1;
	double rho=7750; //Physical parameters
	double cp=500;
	double kk=16.2;
	double qvolx=0; //Volumetric heat generation in x and y direction
	double qvoly=0;
	double qvolz=0;
	clock_t t1,t2;
	t1=clock();
//	double unstead[kmax-1][imax-1][jmax-1]; //Unsteadiness used for convergence
	double alpha; //Courant number
	double unsteadi=1;
	x=new double**[Kmax];
	y=new double**[Kmax];
	z=new double**[Kmax];
	xx=new double**[kmax];
	xy=new double**[kmax];
	xz=new double**[kmax];
	Qxold=new double**[kmax];
	Qyold=new double**[kmax];
	Qzold=new double**[kmax];
	qold=new double**[kmax];
	q=new double**[kmax];
	unstead=new double**[kmax-1];
	for(k=0;k<Kmax;k++)
	{
		x[k]=new double*[Jmax];
		y[k]=new double*[Jmax];
		z[k]=new double*[Jmax];
		for(j=0;j<Jmax;j++)
		{
			x[k][j]=new double[Imax];
			y[k][j]=new double[Imax];
			z[k][j]=new double[Imax];
			
		}
	}
	for(k=0;k<kmax;k++)
	{
		xx[k]=new double*[jmax];
		xy[k]=new double*[jmax];
		xz[k]=new double*[jmax];
		Qxold[k]=new double*[jmax];
		Qyold[k]=new double*[jmax];
		Qzold[k]=new double*[jmax];
		qold[k]=new double*[jmax];
		q[k]=new double*[jmax];
		for(j=0;j<jmax;j++)
		{
			xx[k][j]=new double[imax];
			xy[k][j]=new double[imax];
			xz[k][j]=new double[imax];
			Qxold[k][j]=new double[imax];
			Qyold[k][j]=new double[imax];
			Qzold[k][j]=new double[imax];
			qold[k][j]=new double[imax];
			q[k][j]=new double[imax];
		}
	}
	for(k=0;k<kmax-1;k++)
	{
		unstead[k]=new double*[jmax-1];
		for(j=0;j<jmax-1;j++)
		{
			unstead[k][j]=new double[imax-1];
		}
	}
	T = new double **[Kmax];
	double **Ti = new double *[Jmax*Kmax];
	double *Tii = new double[Imax*Jmax*Kmax];
	for(k=0; k<Kmax; k++, Ti += Jmax)
	{
		T[k] = Ti;
		for(j=0; j<Jmax; j++, Tii += Imax)
			T[k][j] = Tii;
	}
	Told = new double **[Kmax];
	double **Toldi = new double *[Jmax*Kmax];
	double *Toldii = new double[Imax*Jmax*Kmax];
	for(k=0; k<Kmax; k++, Toldi += Jmax)
	{
		Told[k] = Toldi;
		for(j=0; j<Jmax; j++, Toldii += Imax)
			Told[k][j] = Toldii;
	}
	qx = new double **[kmax];
	double **qxi = new double *[kmax*jmax];
	double *qxii = new double[imax*jmax*kmax];
	for(k=0; k<kmax; k++, qxi += jmax)
	{
		qx[k] = qxi;
		for(j=0; j<jmax; j++, qxii += imax)
			qx[k][j] = qxii;
	}
	qy = new double **[kmax];
	double **qyi = new double *[jmax*kmax];
	double *qyii = new double[imax*jmax*kmax];
	for(k=0; k<kmax; k++, qyi += jmax)
	{
		qy[k] = qyi;
		for(j=0; j<jmax; j++, qyii += imax)
			qy[k][j] = qyii;
	}
	qz = new double **[kmax];
	double **qzi = new double *[jmax*kmax];
	double *qzii = new double[imax*jmax*kmax];
	for(k=0; k<kmax; k++, qzi += jmax)
	{
		qz[k] = qzi;
		for(j=0; j<jmax; j++, qzii += imax)
			qz[k][j] = qzii;
	}
/*	double x[Kmax][Imax][Jmax]; //x coordinate for Temperature
	double y[Kmax][Imax][Jmax]; //y coordinate for Temperature
	double z[Kmax][Imax][Jmax];
	double T[Kmax][Imax][Jmax]; //Temperature
	double Tol[Kmax][Imax][Jmax]; //To store temperature values
	double xx[Kmax-1][Imax-1][Jmax-1]; //x coordinate for flux
	double xy[kmax][imax][jmax]; //y coordinate for flux
	double xz[kmax][imax][jmax];
	double Qxold[kmax][imax][jmax]; //To store conduction heat transfer in x direction and y direction
	double Qyold[kmax][imax][jmax];
	double Qzold[kmax][imax][jmax];
	double qold[kmax][imax][jmax];
	double qx[kmax][imax][jmax]; //Conduction flux
	double qy[kmax][imax][jmax];
	double qz[kmax][imax][jmax];*/
	double Qgen;
	
	FILE *fp,*fp1;
	char *fname1="final.dat";
	char *fname2="initial.dat";
	//Calculating required parameters
	dx=Lx/(imax-1);
	dy=Ly/(jmax-1);
	dz=Lz/(kmax-1);
	dt=0.15;
	alpha=k/(rho*cp);
	fp=fopen(fname1,"w");
	fp1=fopen(fname2,"w");
	//Domain
	//Defining coordinates for flux
	for(k=0;k<kmax;k++)
	{
		for(i=0;i<imax;i++)
	{
		for(j=0;j<jmax;j++)
		{
				xx[k][i][j]=j*dx;
				xy[k][i][j]=i*dy;
				xz[k][i][j]=k*dz;
//				printf("\n xx[%d][%d][%d]=%f",k,i,j,xy[k][i][j]);
		}
	}	
	}

	//Defining coordinates for temperature
	for(k=1;k<kmax;k++)
	{
		for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
				x[k][i][j]=(xx[k][i][j-1]+xx[k][i][j])/2;
				y[k][i][j]=(xy[k][i-1][j]+xy[k][i][j])/2;
				z[k][i][j]=(xz[k-1][i][j]+xz[k][i][j])/2;
//			printf("\n x[%d][%d][%d]=%f",k,i,j,z[k][i][j]);
		}
	}	
	}
	
	for(j=1;j<jmax;j++)
	{
		x[0][0][j]=x[1][1][j];
		x[0][Imax-1][j]=x[1][1][j];
		y[0][0][j]=0;
		y[0][Imax-1][j]=Ly;
		z[0][0][j]=0;
		z[0][Imax-1][j]=0;
		x[Kmax-1][0][j]=x[1][1][j];
		y[Kmax-1][0][j]=0;
		z[Kmax-1][0][j]=Lz;
		x[Kmax-1][Imax-1][j]=x[1][1][j];
		y[Kmax-1][Imax-1][j]=Ly;
		z[Kmax-1][Imax-1][j]=Lz;
//		printf("x[0][0][%d] = %f\n",j,x[0][Imax-1][j]);
	}
	for(i=1;i<imax;i++)
	{
		x[0][i][0]=0;
		y[0][i][0]=y[1][i][1];
		z[0][i][0]=0;
		x[0][i][Jmax-1]=Lx;
		y[0][i][Jmax-1]=y[1][i][1];
		z[0][i][Jmax-1]=0;
		x[Kmax-1][i][0]=0;
		y[Kmax-1][i][0]=y[1][i][1];
		z[Kmax-1][i][0]=Lz;
		x[Kmax-1][i][Jmax-1]=Lx;
		y[Kmax-1][i][Jmax-1]=y[1][i][1];
		z[Kmax-1][i][Jmax-1]=Lz;
	}
	for(k=1;k<kmax;k++)
	{
		x[k][0][0]=0;
		x[k][Imax-1][0]=0;
		y[k][0][0]=0;
		y[k][Imax-1][0]=Ly;
		z[k][0][0]=z[k][1][1];
		z[k][Imax-1][0]=z[k][1][1];
		x[k][0][Jmax-1]=Lx;
		y[k][0][Jmax-1]=0;
		z[k][0][Jmax-1]=z[k][1][1];
		x[k][Imax-1][Jmax-1]=Lx;
		y[k][Imax-1][Jmax-1]=Ly;
		z[k][Imax-1][Jmax-1]=z[k][1][1];
//		printf("x[0][0][%d] = %f\n",j,x[0][Imax-1][j]);
	}
	for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
			x[0][i][j]=x[1][1][j];
			y[0][i][j]=y[1][i][1];
			z[0][i][j]=0;
			x[Kmax-1][i][j]=x[1][1][j];
			y[Kmax-1][i][j]=y[1][i][1];
			z[Kmax-1][i][j]=Lz;
		}
	}
		for(k=1;k<kmax;k++)
	{
		for(j=1;j<jmax;j++)
		{
			x[k][Imax-1][j]=x[1][1][j];
			y[k][Imax-1][j]=Ly;
			z[k][Imax-1][j]=z[k][1][1];
			x[k][0][j]=x[1][1][j];
			y[k][0][j]=0;
			z[k][0][j]=z[k][1][1];
		}
	}
		for(i=1;i<imax;i++)
	{
		for(k=1;k<kmax;k++)
		{
			x[k][i][Jmax-1]=Lx;
			y[k][i][Jmax-1]=y[1][i][1];
			z[k][i][Jmax-1]=z[k][1][1];
			x[k][i][0]=0;
			y[k][i][0]=y[1][i][1];
			z[k][i][0]=z[k][1][1];
		}
	}
	x[0][0][0]=0;
	x[Kmax-1][Imax-1][Jmax-1]=Lx;
	x[Kmax-1][0][0]=0;
	x[Kmax-1][Imax-1][0]=0;
	x[Kmax-1][0][Jmax-1]=Lx;
	x[0][0][Jmax-1]=Lx;
	x[0][Imax-1][0]=0;
	x[0][Imax-1][Jmax-1]=Lx;
	
	y[0][0][0]=0;
	y[Kmax-1][Imax-1][Jmax-1]=Ly;
	y[Kmax-1][0][0]=0;
	y[Kmax-1][Imax-1][0]=Ly;
	y[Kmax-1][0][Jmax-1]=0;
	y[0][0][Jmax-1]=0;
	y[0][Imax-1][0]=Ly;
	y[0][Imax-1][Jmax-1]=Ly;
	
	z[0][0][0]=0;
	z[Kmax-1][Imax-1][Imax-1]=Lz;
	z[Kmax-1][0][0]=Lz;
	z[Kmax-1][Imax-1][0]=Lz;
	z[Kmax-1][0][Imax-1]=Lz;
	z[0][0][Jmax-1]=0;
	z[0][Imax-1][0]=0;
	z[0][Imax-1][Jmax-1]=0;

/*/Printing the domain	
	for(k=0;k<Kmax;k++)
	{
		for(i=0;i<Imax;i++)
		{
			for(j=0;j<Jmax;j++)
			{
				fprintf(fp,"\n z[%d][%d][%d]=%f",k,i,j,z[k][i][j]);
			}
		}
	}
*/	//Initial condition
	for(k=1;k<kmax;k++)
	{
		for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
		T[k][i][j]=30;	
		}
	}
		
	}
//Boundary condition
for(j=0;j<Jmax;j++)
{
	for(i=0;i<Imax;i++)
	{	
		T[0][i][j]=100;
		T[Kmax-1][i][j]=300;
	}
}
for(k=0;k<Kmax;k++)
{
	for(j=0;j<Jmax;j++)
	{
		T[k][0][j]=400;
		T[k][Imax-1][j]=200;
	}
}
for(k=0;k<Kmax;k++)
{
	for(i=0;i<Imax;i++)
	{
		T[k][i][0]=500;
		T[k][i][Jmax-1]=600;
	}
}
	//Print Initial condition
	printf("\n");
	for(k=0;k<Kmax;k++)
	{
	for(i=0;i<Imax;i++)
	{
		for(j=0;j<Jmax;j++)
		{
			if(i==0&&j==0&&k==0)
			{
			fprintf(fp1,"VARIABLES = \"Z\", \"X\", \"Y\", \"T\"\n");
			fprintf(fp1,"ZONE I=%d, J=%d, K=%d, F=POINT",Imax,Jmax,Kmax);	
			}
			fprintf(fp1,"\n%0.6f %0.6f %0.6f %0.6f ",z[k][i][j],x[k][i][j],y[k][i][j],T[k][i][j]);
		}
	}
	}
	fclose(fp1);
	//Calculating volumetric heat generation
	Qgen=(qvolx*dx)+(qvoly*dy)+(qvolz*dz);
	//Iteration begins here
	while(iter<=383045)
	{
	iter+=1;
	unsteadi=0;
	for(k=0;k<Kmax;k++)
	{
	for(i=0;i<Imax;i++)
	{
		for(j=0;j<Jmax;j++)
		{
			Told[k][i][j]=T[k][i][j];
		}
	}
	}
	//Calculating flux in x direction
	for(k=0;k<kmax;k++)
	{
	for(i=0;i<imax;i++)
	{
	for(j=0;j<jmax;j++)
	{
	if((j==0)||(j==jmax-1))
		qx[k][i][j]=-(2*kk)*((Told[k][i][j+1]-Told[k][i][j])/dx);
	else
		qx[k][i][j]=-kk*((Told[k][i][j+1]-Told[k][i][j])/dx);
		
	//	printf("\n qx[%d][%d]=%f",i,j,qx[i][j]);
	}
	}
	}
	//Caculating flux in y direction
	for(k=0;k<kmax;k++)
	{
	for(i=0;i<imax;i++)
	{
	for(j=0;j<jmax;j++)
	{
		if((i==0)||(i==imax-1))
		qy[k][i][j]=-(2*kk)*((Told[k][i+1][j]-Told[k][i][j])/dy);
		else
		qy[k][i][j]=-kk*((Told[k][i+1][j]-Told[k][i][j])/dy);
		
	//	printf("\n qy[%d][%d]=%f",i,j,qy[i][j]);
	}	
	}
	}
	//Caculating flux in y direction
	for(k=0;k<kmax;k++)
	{
	for(i=0;i<imax;i++)
	{
	for(j=0;j<jmax;j++)
	{
		if((k==0)||(k==kmax-1))
		qz[k][i][j]=-(2*kk)*((Told[k+1][i][j]-Told[k][i][j])/dz);
		else
		qz[k][i][j]=-kk*((Told[k+1][i][j]-Told[k][i][j])/dz);
		
	//	printf("\n qy[%d][%d]=%f",i,j,qy[i][j]);
	}
	}
	}
	//Calculating temperature from flux
	for(k=1;k<kmax;k++)
	{
	for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
			Qxold[k][i][j]=qx[k][i][j-1]-qx[k][i][j];
		}
	}
	}
	for(k=1;k<kmax;k++)
	{
		for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
			Qyold[k][i][j]=qy[k][i-1][j]-qy[k][i][j];
		}
	}
	}
		for(k=1;k<kmax;k++)
	{
		for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
			Qzold[k][i][j]=qz[k-1][i][j]-qz[k][i][j];
		}
	}
	}
	for(k=1;k<kmax;k++)
	{
	for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
			T[k][i][j]=Told[k][i][j]+(dt/(rho*cp*dx))*(Qxold[k][i][j])+(dt/(rho*cp*dy))*(Qyold[k][i][j])+(Qzold[k][i][j])*(dt/(rho*cp*dz))+Qgen;
		//	printf("T[%d][%d]=%f \t",i,j,T[i][j]);
			
		}
		//printf("\n");
	}
	}
	if((iter%2000)==0)
	{
	printit(iter);
	}
	//Calculating unsteadiness - assuming Lc=Lx
	
	/*for(i=1;i<imax;i++)
	{
		for(j=1;j<jmax;j++)
		{
			if((T[i+1][j]-T[i][j-1])!=0)
			unstead[i-1][j-1]=((Lx*Lx)*(T[i][j]-Told[i][j]))/(alpha*(T[i+1][j]-T[i][j-1])*dt);
			else
			unstead[i-1][j-1]=0;
		//	printf("\n %d %d %f",i,j,unstead[i-1][j-1]);
		}
	}
	for(i=0;i<imax-1;i++)
	{
		for(j=0;j<jmax-1;j++)
		{
			if(unstead[i][j]>unsteadi)
			{
				unsteadi=unstead[i][j];
			//	printf("\n %f",unsteadi);
			}
		}
	}
	//printf("\n Iteration : %d and unsteadiness = %f",iter,unsteadi);
	*/
		//Calculating unsteadiness - with only temperature difference
/*	for(k=1;k<kmax-1;k++)
	{
	for(i=1;i<imax-1;i++)
	{
		for(j=1;j<jmax-1;j++)
		{
			unstead[k][i][j]=T[k][i][j]-Told[k][i][j];
			if(unstead[k][i][j]<0)
			{
				unstead[k][i][j]=-unstead[k][i][j];
			}
		//	printf("\n %d %d %f",i,j,unstead[i-1][j-1]);
		}
	}
	}*/
/*	for(k=1;k<kmax-1;k++)
	{
	for(i=1;i<imax-1;i++)
	{
		for(j=1;j<jmax-1;j++)
		{
			if(unstead[k][i][j]>unsteadi)
			{
				unsteadi=unstead[k][i][j];
			//	printf("\n %f",unsteadi);
			}
		}
	}
	}*/
//	printf("\n Iteration : %d unsteadiness: %0.6f",iter,unsteadi);
	//printf("\n %d",iter);
	}
	for(k=0;k<Kmax;k++)
	{
	for(i=0;i<Imax;i++)
	{
		for(j=0;j<Jmax;j++)
		{
			if(i==0&&j==0&&k==0)
			{
			fprintf(fp,"VARIABLES = \"Z\", \"X\", \"Y\", \"T\"\n");
			fprintf(fp,"ZONE I=%d, J=%d, K=%d, F=POINT",Imax,Jmax,Kmax);	
			}
			fprintf(fp,"\n%0.6f %0.6f %0.6f %0.6f ",z[k][i][j],x[k][i][j],y[k][i][j],T[k][i][j]);
		}
	}
	}
	fclose(fp);
	t2=clock();
	printf("\n Time taken is : %0.6f",(double)((t2-t1)*1.0/CLOCKS_PER_SEC));
	return 0;
}
void printit(int iter)
{
	char tem[80]={0};
	int aa;
	char ch[10];
	aa = sprintf(ch, "%d",iter);
	strcat(tem,ch);
	strcat(tem,".dat");
	FILE *fs;
	fs=fopen(tem,"w");
	for(k=0;k<Kmax;k++)
	{
	for(i=0;i<Imax;i++)
	{
		for(j=0;j<Jmax;j++)
		{
			if(i==0&&j==0&&k==0)
			{
			fprintf(fs,"VARIABLES = \"Z\", \"X\", \"Y\", \"T\"\n");
			fprintf(fs,"ZONE I=%d, J=%d, K=%d, F=POINT",Imax,Jmax,Kmax);	
			}
			fprintf(fs,"\n%0.6f %0.6f %0.6f %0.6f ",z[k][i][j],x[k][i][j],y[k][i][j],T[k][i][j]);
		}
	}
	}
	fclose(fs);
}
