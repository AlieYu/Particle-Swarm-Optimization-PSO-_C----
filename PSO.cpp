/************************************************************************/
/* Author:      yuzelin       -  https://github.com/AlieYu               */
/* Email:       zelinyu@qq.com                                            */
/* FileName:    PSO.cpp                                                  */         
/* LastChange:  2015-4-28                                                 */
/************************************************************************/ 

  

#include<stdlib.h>
#include<time.h>
#include<math.h>
#include"PSO.H"

//微粒的无参构造函数
PARTICLE::PARTICLE()
{
	X=0;
	V=0;
	XBest=0;
	Dim=0;
}
//微粒的有参构造函数
PARTICLE::PARTICLE(int n)
{
	Dim=n;
	X=new double[Dim];
	V=new double[Dim];
	XBest=new double[Dim];
}
//微粒的析构函数
PARTICLE::~PARTICLE()
{
	if(Dim)
	{
		delete[] X;
		delete[] V;
		delete[] XBest;
	}
}
//设置微粒的维数
void PARTICLE::SetDim(int d)
{
	if(X)
	{
		delete[] X;
	}
	if(V)
	{
		delete[] V;
	}
	if(XBest)
	{
		delete[] XBest;
	}
	Dim=d;
	X=new double[Dim];
	V=new double[Dim];
	XBest=new double[Dim];
}

//PSO的无参构造函数
PSO::PSO()
{
	Particle=0;
	PNum=0;
	GBestIndex=0;
	Xup=0;
	Xdown=0;
//	W=1;
	W_max=1;
	W_min=0.6;

	C1=2;
	C2=2;
	Com=0;
}
//PSO的有参构造函数
PSO::PSO(int dim,int n)
{
	Particle=new PARTICLE[n];
	for(int i=0;i<n;i++)
	{
		Particle[i].SetDim(dim);
	}
	PNum=n;
	GBestIndex=0;
	Xup=new double[dim];
	Xdown=new double[dim];
	Vmax=new double[dim];
//	W=1;
	W_max=1;
	W_min=0.6;

	C1=2;
	C2=2;
	Com=0;
}
//PSO的析构函数
PSO::~PSO()
{
	if(Particle)
	{
		delete[] Particle;
	}
	if(Xup)
	{
		delete[] Xup;
	}
	if(Xdown)
	{
		delete[] Xdown;
	}
	if(Vmax)
	{
		delete[] Vmax;
	}
}
//设置坐标上界
void PSO::SetXup(double *up)
{
	if(!Particle)
	{
		return;
	}
	for(int i=0;i<Particle[0].Dim;i++)
	{
		Xup[i]=up[i];
	}
}
//设置坐标下界
void PSO::SetXdown(double *d)
{
	if(!Particle)
	{
		return;
	}
	for(int i=0;i<Particle[0].Dim;i++)
	{
		Xdown[i]=d[i];
	}
}
//设置最大速度
void PSO::SetVmax(double *max)
{
	if(!Particle)
	{
		return;
	}
	for(int i=0;i<Particle[0].Dim;i++)
	{
		Vmax[i]=max[i];
	}
}
void PSO::SetVmax(double p)
{
	if(!Particle)
	{
		return;
	}
	for(int i=0;i<Particle[0].Dim;i++)
	{
		Vmax[i]=(Xup[i]-Xdown[i])*p;
	}
}
//初始化群体
void PSO::Initialize()
{
	if(!Particle)
	{
		return;
	}
	static int kk=(unsigned)time(NULL);
	srand((unsigned)time(NULL)+kk++);

	GBestIndex=0;
	
	//初始化所有粒子的个体
	for(int i=0;i<PNum;i++)
	{
		for(int j=0;j<Particle[i].Dim;j++)
		{
			Particle[i].X[j]=rand()/(double)RAND_MAX*(Xup[j]-Xdown[j])+Xdown[j];//随机初始化坐标
			Particle[i].XBest[j]=Particle[i].X[j];
			Particle[i].V[j]=(rand()/(double)RAND_MAX-0.5)* Vmax[j];//随机初始化速度
		}
		Particle[i].Fit=GetFit(Particle[i]);//计算每个微粒适合度
		Particle[i].FitBest=Particle[i].Fit;//计算最优适合度值
		if(Particle[i].Fit>Particle[GBestIndex].Fit)
		{
			//如果这个鸟的适合度大于群体的最大适合度的话，记录下查找群体的最优微粒
			GBestIndex=i;
		}
	}
}
//计算群体各个微粒的适合度
void PSO::CalFit()
{
	if(!Particle)
	{
		return;
	}
	for(int i=0;i<PNum;i++)
	{
		Particle[i].Fit=GetFit(Particle[i]);
	}
}
//微粒飞翔，产生新一代微粒
void PSO::ParticleFly()
{
	static double FitBak[100];//用来存放备份的适合度值
	if(!Particle)
	{
		return;
	}
	static int tt=(unsigned)time(NULL);
	srand((unsigned)time(NULL)*tt++);

	static int kk=2;//迭代次数
	double W;
	W=W_max-kk*(W_max-W_min)/IteorMax;//由于W由IteorMax计算，Run前需要先设置IteorMax。
	kk++;

	//整个群体飞向新的位置
	for(int i=0;i<PNum;i++)
	{
		for(int j=0;j<Particle[i].Dim;j++)
		{
			Particle[i].V[j]=W*Particle[i].V[j]+
							 rand()/(double)RAND_MAX*C1*(Particle[i].XBest[j]-Particle[i].X[j])+
							 rand()/(double)RAND_MAX*C2*(Particle[GBestIndex].XBest[j]-Particle[i].X[j]);

		}
		for(int j=0;j<Particle[i].Dim;j++)
		{
			if(Particle[i].V[j]>Vmax[j])
			{
				Particle[i].V[j]=Vmax[j];
			}
			if(Particle[i].V[j]<-Vmax[j])
			{
				Particle[i].V[j]=-Vmax[j];
			}
		}
		for(int j=0;j<Particle[i].Dim;j++)
		{
			Particle[i].X[j]+=Particle[i].V[j];//修改坐标
			if(Particle[i].X[j]>Xup[j])
			{
				Particle[i].X[j]=Xup[j];
			}
			if(Particle[i].X[j]<Xdown[j])
			{
				Particle[i].X[j]=Xdown[j];
			}
		}
	}
	
	//计算各微粒的适合度
	CalFit();
	for(i=0;i<PNum;i++)
	{
		FitBak[i]=Particle[i].Fit;
	}
	//设置个体的最好位置
	for(int i=0;i<PNum;i++)
	{
		if(Particle[i].Fit>=Particle[i].FitBest)
		{
			Particle[i].FitBest=Particle[i].Fit;
			for(int j=0;j<Particle[i].Dim;j++)
			{
				Particle[i].XBest[j]=Particle[i].X[j];
			}
		}
	}

	//设置群体中新的最优个体
	GBestIndex=0;
	for(i=0;i<PNum;i++)
	{
		if( (Particle[i].FitBest>=Particle[GBestIndex].FitBest) && i!=GBestIndex)
		{
			GBestIndex=i;
		}
	}
}

//按最多运行次数运行群粒算法，返回最优粒子
PARTICLE& PSO::Run(int n)
{
	Initialize();
	double *opt_p=new double[Particle[0].Dim];//通讯用数组，最优点坐标
	double **opt_a=new double*[PNum];		  //通讯用数组，所有点坐标

	for(int i=0;i<n;i++)
	{
		ParticleFly();
		if(Com)			//通讯函数存在，完成通讯
		{
			for(int k=0;k<Particle[0].Dim;k++)
			{
				opt_p[k]=Particle[GBestIndex].XBest[k];//拷贝当前最优点坐标
			}
			for(int k=0;k<PNum;k++)
			{
				opt_a[k]=Particle[k].X;//指向所有点坐标
			}
			if(!Com(Particle[GBestIndex].FitBest,opt_p,opt_a,GBestIndex))
			{
				break;
			}
		}
	}
	delete[] opt_p;
	delete[] opt_a;
	return Particle[GBestIndex];
}
//按最佳适合度运行群粒算法
PARTICLE& PSO::Run(double fit)
{
	Initialize();
	double *opt_p=new double[Particle[0].Dim];//通讯用数组，最优点坐标
	double **opt_a=new double*[PNum];		  //通讯用数组，所有点坐标

	do
	{
		ParticleFly();
		if(Com)			//通讯函数存在，完成通讯
		{
			for(int k=0;k<Particle[0].Dim;k++)
			{
				opt_p[k]=Particle[GBestIndex].XBest[k];//拷贝最优点坐标
			}
			for(int k=0;k<PNum;k++)
			{
				opt_a[k]=Particle[k].X;//指向所有点坐标
			}
			if(!Com(Particle[GBestIndex].FitBest,opt_p,opt_a,GBestIndex))
			{
				break;
			}
		}
	}while(Particle[GBestIndex].FitBest<fit);
	delete[] opt_p;
	delete[] opt_a;
	return Particle[GBestIndex];
}
//返回最佳个体
double PSO::GetBest(double *r)
{
	for(int i=0;i<Particle[GBestIndex].Dim;i++)
	{
		r[i]=Particle[GBestIndex].XBest[i];
	}
	return Particle[GBestIndex].FitBest;
}