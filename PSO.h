/************************************************************************/
/* Author:      yuzelin       -  https://github.com/AlieYu               */
/* Email:       zelinyu@qq.com                                            */
/* FileName:    PSO.h                                                       */           
/* LastChange:  2015-4-28                                                 */
/************************************************************************/   
#ifndef _PSO_H
#define _PSO_H

//微粒类
class PARTICLE
{
public:
	double *X;			//微粒的坐标数组
	double *V;			//微粒的速度数组
	double *XBest;		//微粒的最好位置数组
	int    Dim;			//微粒的维数
	double Fit;			//微粒的适合度
	double FitBest;		//微粒最好位置的适合度
	
	PARTICLE();		    //空构造函数
	PARTICLE(int n);    //维数为参数的构造函数

	~PARTICLE();	    //析构函数

	void SetDim(int d); //设置微粒的维数

};

//群粒子类
class PSO
{
protected:
	PARTICLE *Particle; //微粒群数组
	int    PNum;	    //微粒个数
	int    GBestIndex;	//最好微粒索引
	double W_max;		//惯性权重的最大值
	double W_min;		//惯性权重的最小值
	int    IteorMax;	//最大迭代次数

	double C1;			//加速度系数1
	double C2;			//加速度系数2
	double *Xup;		//微粒坐标上界数组
	double *Xdown;		//微粒坐标下界数组
	double *Vmax;		//微粒最大速度数组

	void Initialize();			//初始化群体
	void CalFit();				//计算全体适合度
	virtual void ParticleFly();	//微粒飞翔，产生新一代微粒

	
	//通讯函数，返回值为false时，系统停止优化

	bool (*Com)(double, /*最优微粒适合度*/ double*, /*最优微粒坐标数组*/
				double**, /*所有微粒坐标指针数组*/ int /*当前最优微粒索引*/ );

public:
	PSO();						//空构造函数
	PSO(int dim,int n);			//dim为微粒维数，n为微粒个数

	~PSO();						//析构函数

	void SetXup(double*);		//设置微粒坐标上界
	void SetXdown(double*);		//设置微粒坐标下界
	void SetVmax(double*);		//设置微粒最大速度,以数组为参数
	void SetVmax(double);		//设置微粒最大速度，以坐标的上下界的百分比为参数
	void SetC1(double c){C1=c;}	//设置C1
	void SetC2(double c){C2=c;}	//设置C2
	void SetCom(void *p)		//设置通讯函数
	{
		Com=(bool(*)(double,double*,double**,int))p;
	}

	void SetIteorMax(int iteor)
	{
		IteorMax=iteor;			//设置最大迭代次数
	}
	//计算特定微粒坐标所对应的适合度，必须由派生类的实际PSO类定义,以便计算适合度
	virtual float GetFit(PARTICLE&)=0;

	PARTICLE& Run(int max);		//按最多次数限制运行PSO
	PARTICLE& Run(double fit);	//按最佳适合度目标运行PSO
	double GetBest(double*);	//获得最佳微粒适合度和坐标

};

#endif
