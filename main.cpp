#include <glut\glut.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <math.h>
#include <time.h>

using namespace std;

float crossX;
float crossY;
int winW = 1280;
int winH = 720;
int mousePosX;
int mousePosY;

float SLOMO_K = 1;

class PointF {
public:
	float x;
	float y;

	void setPos(float x, float y) {
		this->x = x;
		this->y = y;
	}
};

class Sector
{
public:
	int num1;
	int num2;
	Sector()
	{
		num1=0; num2=0;
	}
};

Sector sector[100][100];

int curDescI, curDescJ;

class Point3d
{
public:
	float x,y,z;	
	Point3d(float x,float y,float z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
};

//bool DEBUG_MODE = false;
float dT=1;
float shadow_volume = 0.9;
bool SH1 = false;
float ANGLE = 0.0;
float UDANG = 0.0;
int WIDTH = 100;
int ZFAR = 10000;
int viewAngle = 40;
int DISTANCETOSCENE = 1200;
bool pauseMode = false;


#define		_torad	0.017453293 //3.1415926/180;
#define		_toa	57.29577952 //180/3.1415926;
#define MAX(a,b) ((a < b) ? (b) : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SIGN(x) ((x) < (0) ? (-1) : (1))
#define SQR(x)  ((x) * (x))

int TOTALCOUNT =100;


#define MAXPRIMCOUNT 20000
#define MAXENEMYCOUNT 1000
#define MAXDOORCOUNT 100
#define CHECKPOINT_COUNT 10
#define MAXITEMS 30

#define DEFAULT_GRAVITY 0.16


//Line *line[3000]; 


#define GOOD 0
#define EVIL 1

#define MAX_INTERVALS 8


GLboolean lightingSaveState;

void pushLighting()
{
	lightingSaveState = glIsEnabled(GL_LIGHTING);	
}

void popLighting()
{	
	if (lightingSaveState == GL_TRUE)
		glEnable(GL_LIGHTING);
	else 
		glDisable(GL_LIGHTING);
}

class Timers
{
	int intervals[MAX_INTERVALS];
public:
	int starttime[MAX_INTERVALS];
	int time[MAX_INTERVALS];
	int typeOfImpact;
	//	vector<int> intervals;
	//	vector<int>	::iterator pInter;

	//	int intervalMilliSec;
			 
	
	Timers()
	{
		typeOfImpact = GOOD;
		restart();

		for(int i=0;i<MAX_INTERVALS;i++)
			intervals[i] = 0;
	}

	virtual void 	actionOnTime(int numAction) { }

	 
	void setInterval(int index, int intervalsMilliSec)
	{
		if(index >= 0 && index < MAX_INTERVALS)
			intervals[index] = intervalsMilliSec;

		starttime[index] = glutGet(GLUT_ELAPSED_TIME) - rand()% (1 + glutGet(GLUT_ELAPSED_TIME));
		restart();
	}


	bool timer()
	{	
		for(int i=0;i<MAX_INTERVALS;i++)
			time[i] = glutGet(GLUT_ELAPSED_TIME)-starttime[i]; 

		bool flag = false;
		for(int i=0;i<MAX_INTERVALS;i++)
		{
			if ( intervals[i]!=0 && (time[i] >= intervals[i] || starttime[i] == glutGet(GLUT_ELAPSED_TIME)) )
			{	starttime[i] += time[i];
				actionOnTime(i);
				flag = true;
			}
		}
		return flag;
	}

	void restart()
	{	
		for(int i=0;i<MAX_INTERVALS;i++)
			starttime[i] = glutGet(GLUT_ELAPSED_TIME);
	}
};

class FVector
{
public: float x0,y0,z0, // откуда смотрит вектор
			x1,y1,z1; // куда смотрит вектор

	float r;
	float R;
	static float gravity;
	
	FVector()
	{
		//
		// TODO: Add constructor logic here
		//
	}
	void calc_kas(float oldx,  
		float oldy,  
		// flolat oldz,
		float x,  
		float y//,  
		// flolat z
		) 
	{ 
		x0=x;
		y0=y;
		//	z0=z;

		x1=x+(x-oldx);
		y1=y+(y-oldy);
		//z1=z+(z-oldz);
	}

	void calc_gravity()
	{
		y1-=gravity/SLOMO_K;
	}

	void calc_antigravity()
	{
		y1+=gravity/SLOMO_K;
	}

	void calc_in_water( )
	{	
		y1+= (gravity/SLOMO_K + gravity/SLOMO_K*0.4);
	}

	void calc_veter(float massa, float xa,float ya )
	{
		x1+=gravity/SLOMO_K;
	}

	void add(FVector *vector)
	{

	}

	void addImpulse(float impx, float impy)
	{
		this->x1 += impx; 
		this->y1 += impy;
	}
};

void ResultofAddVect(FVector *FRES,FVector A,FVector B)
{

	FRES->x0=A.x0;
	FRES->y0=A.y0;
	FRES->x1=A.x1+B.x1-FRES->x0;
	FRES->y1=A.y1+B.y1-FRES->y0;
}

class Unit : public Timers, public PointF
{
public:

	float dx,dy;
	float oldx,oldy;
	bool inWater;
	bool onHard;
	bool onAlienUnit;
	float radius;
	short isMassive;
	FVector FMove;
	bool friction;
	bool spirit;
	float inertion;
	
	Unit()
	{
		x=0;
		y=0;
		oldx=x;
		oldy=y;
		dx=0;
		dy=0;
		radius = 0;
		isMassive = 1;
		friction = true;
		spirit = false;
		inWater = false;
		inertion = 1;

	}
	
	void setPos(float x, float y)
	{
		this->x = x; 
		this->y = y;
		oldx=x;
		oldy=y;
	}

	void setRadius(float r)
	{
		this->radius = r; 
	}
	void setDxDy(float dx, float dy)
	{
		this->dx = dx/SLOMO_K; 
		this->dy = dy/SLOMO_K;
	}
	void setMassive(short value)
	{
		isMassive = value;
	}
	void setFriction(bool fri)
	{
		friction = fri;
	}
	void setSpirit(bool spi)
	{
		spirit = spi;
	}
	void setInertion(float ine)
	{
		inertion= ine;
	}
	void CalcVector()
	{
		if(inertion)
			FMove.calc_kas(oldx,oldy,x,y);

		if(isMassive==1)
		{
			FMove.calc_gravity();

			if(inWater) 
				FMove.calc_in_water(); 
		}

		if(isMassive == -1)
		{
			FMove.calc_antigravity();

			if(inWater) 
				FMove.calc_in_water(); 
		}
	}

	void CalcDxDy()
	{	
		dx = inertion*(FMove.x1 - FMove.x0)/SLOMO_K;
		dy = inertion*(FMove.y1 - FMove.y0)/SLOMO_K ;

		if(isMassive)
			if(inWater) 
				Trenie(0.985);

	}

	void Trenie(float trenie)
	{	
		dx*=trenie;
		dy*=trenie;
	}
	
	void Move()
	{
		x+=dx;
		y+=dy ;

		oldx=x-dx*SLOMO_K;
		oldy=y-dy*SLOMO_K;				
	}
};

class RopeElement: public Unit
{
public: //FVector FSh;
	//FVector FMove;

	bool fix;

	//	 float old_dx,old_dy;		

	bool blue ;
	bool empty ;
	bool visible;
	bool enable ;
	bool moveble;

	int massa;


	float oldrad ;

	bool inprocess ;

	bool Inert;

	int iter;
	float cred,
		cgreen,
		cblue;

	bool cc;
	int coltype ;

	int extype;


	int NUM;
	
	RopeElement(float x, float y, float rad,int NUM)
	{	cc =false;
		iter=0;
		massa = 1;		
		blue = true;
		empty = true;
		visible = true;
		enable = true;
		moveble = false;

		Inert = true;
		fix = false;
		this->NUM = NUM;		
		inprocess = true;

		cred = 0.7;
		cgreen =0.7;
		cblue = 0.7;

		oldx=this->x = x;
		oldy=this->y = y;

		oldrad = this->radius = rad;

		dx = dy = 0;
	}
	void actionOnTime(int intervalMilliSec)
	{

	}

	int soot(int a)
	{
		if(a<18) return a+18;
		else return a-18;
	}

	void CalcVectorInWater()
	{

		FMove.calc_kas(oldx,oldy,x,y);
		FMove.calc_in_water();

	}

	void CalcDxDyInWater()
	{
		///	old_dx = dx;
		//	old_dy = dy;		
		dx = (FMove.x1 - FMove.x0);
		dy = (FMove.y1 - FMove.y0);		
		dx*=0.985;
		dy*=0.985;
	}

	void Show(float addx,float addy)
	{
		if(visible)
		{	
			glPushMatrix();

			glTranslatef(x,y,0);
			glColor3f(cred,cgreen,cblue); 
			glutSolidSphere(radius,8,8);

			glPopMatrix();
		}
	}
};

inline void tr(float ax,float ay,float az,
			   float bx,float by,float bz,
			   float cx,float cy,float cz
			   )
{

	static float v1[3],v2[3],v3[3];
	v1[0] = ax;
	v1[1] = ay;
	v1[2] = az;

	v2[0] = bx;
	v2[1] = by;
	v2[2] = bz;

	v3[0] = cx;
	v3[1] = cy;
	v3[2] = cz;

	float vector1X = ax-bx;
	float vector1Y = ay-by;
	float vector1Z = az-bz;

	float vector2X = bx-cx;
	float vector2Y = by-cy;
	float vector2Z = bz-cz;

	float normalX = vector1Y * vector2Z - vector1Z * vector2Y;     
	float normalY = vector1Z * vector2X - vector1X * vector2Z;
	float normalZ = vector1X * vector2Y - vector1Y * vector2X;

	glNormal3f(normalX,normalY,normalZ);
	glVertex3fv( v1);
	glVertex3fv( v2);
	glVertex3fv( v3);

	glNormal3f(-normalX,-normalY,-normalZ);
	glVertex3fv( v3);
	glVertex3fv( v2);
	glVertex3fv( v1);
}

class Line
{
public:
	PointF a;
	PointF b;

	PointF med;

	PointF KK ;
	float Length ;

	int rotMode ;	
	bool Rotable ;
	bool Bad ;		
	bool Lift ;		

	PointF pnt[4] ;
	float ang;


	float Angle;

	float OldAngle;

	float dA;
	int NUM;
	float Width;

	bool invertNormal;

	Line () {}

	void set(float x1, float y1, float x2, float y2 ,float Width, int num)
	{
		this->Width=Width;
		invertNormal = false;
		Angle = ang=0;
		Rotable = false;
		Bad = false;		
		Lift = false;	
		NUM = num;
		a.x = x1;
		a.y = y1;
		b.x = x2;
		b.y = y2;
		med.x = (a.x+b.x)*0.5f;
		med.y = (a.y+b.y)*0.5f;		
		KK = normal();
		Length = (float)sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y)*(a.y - b.y));
	}

	void SetRotable(bool f1,int f2)
	{
		Rotable = f1;
		rotMode = f2;
	}

	void SetBad(bool f)
	{
		Bad = f;
	}		

	void SetLift(bool f)
	{
		Lift = f;
	}			


	void Rotate (float angle)
	{
		return;
		// PointF rotPnt;

		//	if (rotMode == 0) 
		{	a.x = med.x - (float)cos (angle*_torad)*Length*0.5f;
			a.y = med.y - (float)sin (angle*_torad)*Length*0.5f;

			b.x = med.x + med.x - a.x;
			b.y = med.y + med.y - a.y;

		//		if(NUM<0)
		//	{
		//		b.x = med.x ;
		//		b.y = med.y ;
		//	}
		}

		/*if (rotMode == 1) 
		{
			b.x = a.x + (float)cos (angle*_torad)*Length;
			b.y = a.y + (float)sin (angle*_torad)*Length;
		}
		*/
		if(b.x<a.x)
		{	
			PointF tmpp = a;	a = b; b = tmpp;
		}
	}

	PointF normal()
	{
		PointF rezz;

		float nx = a.x-b.x;
		float ny = a.y-b.y;
		float L = 1/sqrt (nx*nx + ny*ny);
		nx*=L;
		ny*=L; 

		if(!invertNormal)
		{
			rezz.x = ny ;
			rezz.y = -nx;
		}
		else
		{
			rezz.x = -ny;
			rezz.y = nx;
		}

		return rezz; 

	}

	void ShowLite ()
	{	
		tr(a.x,a.y,Width,
			b.x,b.y,Width,
			b.x,b.y,-Width);

		tr(b.x,b.y,-Width,
			a.x,a.y,-Width,
			a.x,a.y,Width);
	}

	void ShowF ()
	{	
		glBegin(GL_TRIANGLES);
		tr(a.x,a.y,Width,
			b.x,b.y,Width,
			b.x,b.y,-30);

		tr(b.x,b.y,-30,
			a.x,a.y,-30,
			a.x,a.y,Width);
		glEnd();
	}

	void ShowShadow (bool flag, PointF *light, float *onlgt)		
	{
		//	 glDisable(GL_DEPTH_TEST);


		static int numLights=(((int)sizeof(light))>>1);


		int numONLights = 0;
		for (int i=0; i<numLights;i++ )
			if(onlgt[i]) numONLights++;


		for(int i=0; i<numLights; i++)
		{
			if(!onlgt[i]) 
				continue;


			//float mtx = (med.x-light[i].x);
			//float mty = (med.y-light[i].y);	

			//float kk = (float)Math.Sqrt(mtx*mtx +mty*mty );	

			float tmp2x = (b.x-light[i].x);
			float tmp2y = (b.y-light[i].y);
			float tmp3x = (a.x-light[i].x);
			float tmp3y = (a.y-light[i].y);





			/*		
			tmp2x /= kk;
			tmp2y /= kk;
			tmp3x /= kk;
			tmp3y /= kk;
			*/
			//	int ck=220-flag*50;
			//	float dk = 0.2f;



			pnt[0].x = a.x;
			pnt[0].y = a.y;			

			pnt[1].x = b.x;
			pnt[1].y = 	b.y;


			pnt[2].x =b.x+tmp2x*1000;
			pnt[2].y =b.y+tmp2y*1000 ;
			pnt[3].x =a.x+tmp3x*1000 ;
			pnt[3].y =a.y+tmp3y*1000 ;


			//glDisable(GL_BLEND);

			glDisable (GL_CULL_FACE);
			glBegin(GL_QUADS);

			//	glColor4f(0,0,0,shadow_volume/numONLights);

			//	glColor4f(0.1,0.1,0.1,0.2);	
			glColor3f(0,0,0);

			if(SH1)
			{
				glVertex3f(pnt[0].x,pnt[0].y,Width+3);
				glVertex3f(pnt[1].x,pnt[1].y,Width+3);
				//	glColor4f(0,0,0,1);	

				glVertex3f(pnt[2].x,pnt[2].y,Width+3);
				glVertex3f(pnt[3].x,pnt[3].y,Width+3);
			}

			glVertex3f(pnt[0].x,pnt[0].y,10);
			glVertex3f(pnt[1].x,pnt[1].y,10);

			glVertex3f(pnt[2].x,pnt[2].y,10);
			glVertex3f(pnt[3].x,pnt[3].y,10);

			/*
			glVertex3f(pnt[2].x,pnt[2].y,-30);
			glVertex3f(pnt[3].x,pnt[3].y,-30);
			glVertex3f(pnt[3].x+tmp3x*500,pnt[3].y+tmp3y*500,-30);
			glVertex3f(pnt[2].x+tmp2x*500,pnt[2].y+tmp2y*500,-30);
			*/

			glEnd();


			//glEnable(GL_BLEND);
			glEnable(GL_CULL_FACE);

			//	if(UDANG>-50) continue;

			/*
			glDisable (GL_CULL_FACE);
			glBegin(GL_QUADS);
			glColor4f(0,0,0,shadow_volume/numONLights);	
			glVertex3f(pnt[0].x,pnt[0].y,150);
			glVertex3f(pnt[1].x,pnt[1].y,150);
			glColor4f(0,0,0,shadow_volume);
			glVertex3f(pnt[2].x,pnt[2].y,150);
			glVertex3f(pnt[3].x,pnt[3].y,150);
			glEnd();
			glEnable(GL_CULL_FACE);
			*/

		}


		//linegr.DrawString(KK.y.ToString(),new Font("arial",10),Brushes.Black,0.5f*(a.x+b.x)-5f,0.5f*(a.y+b.y)+5f);
		//linegr.DrawString(KK.x.ToString(),new Font("arial",10),Brushes.Black,0.5f*(a.x+b.x)-5f,0.5f*(a.y+b.y)-15f);			
		//	glEnable(GL_DEPTH_TEST);
	}
};






float minirand()
{ return (float)(rand()%10000)/10000;
}


int interrRand(int a, int b)
{ return a + rand()%(b-a);
}

double FSIN(float x)
{float xx=x*x;//2
float xxx=xx*x;//3
float xxxxxx=xxx*xxx; //6

return (double)(x-0.166666667*xxx
				+0.008333333*xxx*xx
				-0.000198413*xxxxxx*x
				+0.000002756*xxxxxx*xxx
				-0.000000025*xxxxxx*xxx*xx);

}



double FCOS(float x)
{ 
	float xx=x*x;
	float xxxx=xx*xx;
	float xxxxxxxx=xxxx*xxxx;

	return (double)
		(1.00000000
		-0.500000000*xx
		+0.041666667*xxxx
		-0.001388889*xxxx*xx
		+0.000024802*xxxxxxxx
		-0.000000276*xxxxxxxx*xx
		+0.000000002*xxxxxxxx*xxxx);
}

/*
bool isGraphActual(PointF a,Unit b)
{

}
*/

PointF ToPointF(float x,float y)
{
	PointF rezp;
	rezp.x = x;
	rezp.y = y;
	return rezp;
}


float getAngle(float X1,float Y1,float X2,float Y2)
{
	float result_ugol;
	float b;  b=fabs(Y1-Y2);
	float c;  c=fabs(X1-X2);


	if(Y2<=Y1 && X2>=X1) result_ugol=360-atan(b/c)*_toa;
	if(Y2<=Y1 && X2<X1)  result_ugol=180+atan(b/c)*_toa;
	if(Y2>Y1 && X2<=X1)  result_ugol=180-atan(b/c)*_toa;
	if(Y2>Y1 && X2>X1)   result_ugol=atan(b/c)*_toa;  

	return result_ugol;
}

int numMax3(float a,float b,float c , int *num2 , int *num3)
{
	if(a>b && a>c) {*num2 = 1; *num3 = 2; return 0;}
	if(b>a && b>c) {*num2 = 0; *num3 = 2; return 1;}
	if(c>a && c>b) {*num2 = 0; *num3 = 1; return 2;}
	return 0;
}

//	PointF medium();


#define ITNUM 4

float Sqroot(float x) {
	int expo,i;
	float a,b;
	if(x<=0.F) return(0.F);
	/* decompose x into mantisse ranged [0.5,1) and exponent. Machine-
	dependent operation is presented here as a function call. */
	x=frexp(x,&expo);
	/* odd exponent: multiply mantisse by 2 and decrease exponent 
	making it even. Now the mantisse is ranged [0.5,2.) */
	if(expo&1) {x*=2.F;expo--;}
	/* initial approximation */
	a=1.F;
	/* process ITNUM Newtonian iterations with the mantisse as 
	an argument */
	for(i=ITNUM;i>0;i--) {
		b=x/a; a+=b; a*=0.5F;
	}
	/* divide the exponent by 2 and combine the result.
	Function ldexp() is opposite to frexp. */
	a=ldexp(a,expo/2);
	return(a);
}

void ch(int *a,int *b)
{ 
	int t = *a;
	*a= *b;
	*b = t;
}
void ch(float *a,float *b)
{ 
	float t = *a;
	*a= *b;
	*b = t;
}


bool CircleIntersects(float *x, float *y, float R, float L, 
					  PointF A,PointF  B)
{
	// единичный вектор отрезка AB
	float overL = 1/ L;
	float Xv = (B.x - A.x) * overL ;
	float Yv = (B.y - A.y) * overL;
	float Xd = (A.x - *x);
	float Yd = (A.y - *y);
	float b = 2 * (Xd * Xv + Yd * Yv);
	float c = Xd * Xd + Yd * Yd - (R) * (R);
	float c4 = c + c; c4 += c4;
	float D = b * b - c4;
	if (D < 0) return false; // нет корней, нет пересечений

	D = Sqroot(D);
	float l1 = ( - b + D) * 0.5f;
	float l2 = ( - b - D) * 0.5f;
	bool intersects1 = ((l1 >= 0.0) && (l1 <= L));
	bool intersects2 = ((l2 >= 0.0) && (l2 <= L));
	bool intersects = intersects1 || intersects2;

	return intersects;
}



bool  simple_rasst (Unit *u1,Unit *u2, float R)
{

	return (u1->x< u2->x+R && u1->x > u2->x - R &&
		u1->y< u2->y+R && u1->y > u2->y - R );

	//	return (fabs(X1-X2)<R && fabs(Y1-Y2)<R && fabs(Z1-Z2)<R);

} ;


bool  simple_rasst (PointF p1,PointF p2, float R)
{

	return (p1.x< p2.x+R && p1.x > p2.x - R &&
		p1.y< p2.y+R && p1.y > p2.y - R );

	//	return (fabs(X1-X2)<R && fabs(Y1-Y2)<R && fabs(Z1-Z2)<R);

} ;


bool  simple_rasst (float x1,float y1, float x2,float y2, float R)
{

	return (x1< x2+R && x1 > x2 - R &&
		y1< y2+R && y1 > y2 - R );

	//	return (fabs(X1-X2)<R && fabs(Y1-Y2)<R && fabs(Z1-Z2)<R);

} ;



float getC2 (float A,float B)
{
	static float tmax,b;		
	tmax = MAX(A,B);
	b =  MIN(A,B)/tmax; 
	return  tmax+ 0.425*b*b*tmax;
} 

float Dist(float x1, float y1, float x2, float y2)
{
	//   return  getC2(fabs(x1 - x2),fabs(y1 - y2)); 
	return sqrtf((x1 - x2) * (x1 - x2)+ (y1 - y2) * (y1 - y2)); 
}

float Dist(float *x1, float *y1, float *x2, float *y2)
{
	//   return  getC2(fabs(x1 - x2),fabs(y1 - y2)); 
	return sqrtf((*x1 - *x2) * (*x1 - *x2)+ (*y1 - *y2) * (*y1 - *y2)); 
}

float Dist(PointF *p1,PointF  *p2)
{
	return Dist(&p1->x, &p1->y, &p2->x ,&p2->y);
}

float Dist(PointF p1,PointF  p2)
{
	return Dist(p1.x, p1.y, p2.x ,p2.y);
}

float Dist(Unit *elem1, PointF *p)
{

	return Dist(&elem1->x,&elem1->y, &p->x ,&p->y );

}

float Dist(Unit *u1,Unit  *u2)
{
	return Dist(&u1->x,&u1->y,  &u2->x,&u2->y);		
}

float Dist(Unit *u, float x, float y)
{
	return Dist(&u->x,&u->y, &x, &y);	    
}

float distance(PointF a,PointF  b,PointF c)
{																																																																												
	float dx = a.x - b.x;
	float dy = a.y - b.y;
	float D = dx * (c.y - a.y) - dy * (c.x - a.x);
	return (float)(D / sqrt(dx * dx + dy * dy));
}

float distance(PointF a,PointF  b,Unit *c)
{																																																																												
	float dx = a.x - b.x;
	float dy = a.y - b.y;
	float D = dx * (c->y - a.y) - dy * (c->x - a.x);
	return (float)(D / sqrt(dx * dx + dy * dy));
}


class TrianglePlatform;
bool pointInPlatform(PointF point, TrianglePlatform  *platform);


int SIZE1 = 120 ,SIZE2 = 115; 

float limitX=SIZE1*100;
float limitY=SIZE2*100; 


#define UN   -1
#define GEO   0
#define HYDRO 1
#define ICE	  2
#define SAND  3
#define LAVA  4


struct PlatformView
{
	float color[4];
	float colorDark[4];
};

class TrianglePlatform 
{
public:

	int type; 

	PointF A;
	PointF B;
	PointF C;
	Line line[3];

	float Width;
	int Count;
	PointF boxMIN,boxMAX;
	int NUM;
	bool isHalfBox;
	float dx,dy;
	bool isLegal;
	float medx;
	float medy;	
	float medz;	
	bool show;
	bool isDepth;
	bool uppest;
	int fgh;
	float bmx,bmy;
	PlatformView view;

	PointF st_medium;

	TrianglePlatform()
	{
		this->type = UN;
		//	set(GEO, ToPointF(0,0)	,ToPointF(0,0)	,ToPointF(0,0), 0,true,num);
	}

	void set(int type, PointF a,PointF b,PointF c , float Width ,bool ishb, int NUM  )
	{

		this->type = type;
		this->NUM = NUM;
		isHalfBox = ishb;

		isLegal = false;

		static int cc=0;
		Count = cc;
		cc++;

		fgh = 1100+minirand()*50;

		dx = 0;
		dy=0;

		if(type!=HYDRO  && type!=LAVA)
			Width = 100 ;

		this->Width=Width;

		switch (type)
		{
		case SAND:
			this->setColor(1,1,.5,.5); break;
		case ICE:
			this->setColor(0.75,1,1,0.9); break;
		case GEO:
			this->setColor(.5,.5,.5, 1); break;
			break;
		case HYDRO:
			this->setColor(0,0.7,0.6,0.4); break;
		case LAVA:
			this->setColor(0.5,0,0,0.6); break;		
		default:
			break;
		}

		if(type == ICE) 
		{	
			a.x+=1.0;
			b.x+=1.0;
			c.x+=1.0;
			a.y+=1.0;
			b.y+=1.0;
			c.y+=1.0;
		}

		setPos (a,b,c);
	}

	void setPos(PointF a,PointF b,PointF c)
	{
		this->A = a;
		this->B = b;
		this->C = c;

		bmx = 	(A.x + C.x)*0.5;
		bmy = 	(A.y + C.y)*0.5;

		isDepth=false;

		uppest= false;

		boxMIN.x = MIN(a.x, MIN(b.x,c.x));
		boxMIN.y = MIN(a.y, MIN(b.y,c.y));
		boxMAX.x = MAX(a.x, MAX(b.x,c.x)); 
		boxMAX.y = MAX(a.y, MAX(b.y,c.y)); 
		medx= (A.x+B.x+C.x)/3;
		medy= (A.y+B.y+C.y)/3;
		medy= Width/2;		

		line[0].set(a.x,a.y,b.x,b.y,Width,0);
		line[1].set(b.x,b.y,c.x,c.y,Width,1);
		line[2].set(c.x,c.y,a.x,a.y,Width,2);

		if(pointInPlatform(ToPointF((A.x+B.x)/2+line[0].normal().x,  (A.y+B.y)/2 + line[0].normal().y), this))
			line[0].invertNormal = true;

		if(pointInPlatform(ToPointF((C.x+B.x)/2+line[1].normal().x,  (C.y+B.y)/2 + line[1].normal().y), this))
			line[1].invertNormal = true;

		if(pointInPlatform(ToPointF((A.x+C.x)/2+line[2].normal().x,  (A.y+C.y)/2 + line[2].normal().y), this))
			line[2].invertNormal = true;
	}

	~TrianglePlatform()
	{
		//	delete 	&line[0] ;

		//	line[0] = NULL;
		//	line[1]	= NULL;
		//	line[2] = NULL;
	}

	void actionOnTime(int intervalMilliSec)
	{
	}
	
	void setColor(float r,float g,float b,float a)
	{
		view.color[0]=r;
		view.color[1]=g;
		view.color[2]=b;
		view.color[3]=a;

		view.colorDark[0]=r*0.5;
		view.colorDark[1]=g*0.5;
		view.colorDark[2]=b*0.5;
		view.colorDark[3]=a*0.5;
	}

	int numMin3(float a,float b,float c )
	{
		if(a<=b && a<=c) return 0;
		if(b<=a && b<=c) return 1;
		if(c<=a && c<=b) return 2;

		return 0;
	}

	PointF GetLineNormal(int i)
	{	return line[i].normal();
	}

	Line GetLine(int i)
	{	return line[i];
	}

	PointF GetLineA(int i) 
	{	return line[i].a;
	}

	PointF GetLineB(int i) 
	{	return line[i].b;
	}

	float GetLineLength(int i)
	{	return line[i].Length;
	}

	int numNearLine (PointF P)
	{
		return 	numMin3(fabs(distance(line[0].a, line[0].b, P)),
			fabs(distance(line[1].a, line[1].b, P)),
			fabs(distance(line[2].a, line[2].b, P))
			);
	}

	PointF medium()
	{

		return st_medium;
	}

	bool HalfCheck(PointF P,bool pr)
	{
		PointF a,b,c ;
		if(pr)
		{
			a = this->C;
			b = this->B;
			c = this->A;
		}
		else
		{
			a = this->A;
			b = this->B;
			c = this->C;
		}

		if ((P.x-a.x)*(a.y-b.y) - (P.y-a.y)*(a.x-b.x) >= 0) 
		{
			if ((P.x-b.x)*(b.y-c.y) - (P.y-b.y)*(b.x-c.x) >= 0) 
			{
				if ((P.x-c.x)*(c.y-a.y) - (P.y-c.y)*(c.x-a.x) >= 0) 
					return true;

			}
		}	  

		return false;
	}
	
	void Trenie(RopeElement *b)
	{

		b->dx = 0;
		b->dy = 0;
		b->Move();	
	}

	bool inBox(Unit *u)
	{
		if(    u->x >= boxMIN.x 
			&& u->x <= boxMAX.x
			&& u->y >= boxMIN.y
			&& u->y <= boxMAX.y)
			return true;
		return false;
	}


	void Move()
	{
		A.x+=dx;
		B.x+=dx;
		C.x+=dx;

		A.y+=dy;
		B.y+=dy;
		C.y+=dy;

		boxMIN.x = MIN(A.x, MIN(B.x,C.x));
		boxMIN.y = MIN(A.y, MIN(B.y,C.y));
		boxMAX.x = MAX(A.x, MAX(B.x,C.x)); 
		boxMAX.y = MAX(A.y, MAX(B.y,C.y)); 
		medx= (A.x+B.x+C.x)/3;
		medy= (A.y+B.y+C.y)/3;

		st_medium.x = (boxMIN.x + boxMAX.x) /2;
		st_medium.y = (boxMIN.y + boxMAX.y) /2;	
	}

	void Show()
	{}
};

bool pointInPlatform(PointF point, TrianglePlatform  *platform)
{
	PointF a,b,c ;


	a = (*platform).C;
	b = (*platform).B;
	c = (*platform).A;

	if ((point.x-a.x)*(a.y-b.y) - (point.y-a.y)*(a.x-b.x) >= 0) 
	{
		if ((point.x-b.x)*(b.y-c.y) - (point.y-b.y)*(b.x-c.x) >= 0) 
		{
			if ((point.x-c.x)*(c.y-a.y) - (point.y-c.y)*(c.x-a.x) >= 0) 
				return true;

		}
	}	  


	a = (*platform).A;
	b = (*platform).B;
	c = (*platform).C;


	if ((point.x-a.x)*(a.y-b.y) - (point.y-a.y)*(a.x-b.x) >= 0) 
	{
		if ((point.x-b.x)*(b.y-c.y) - (point.y-b.y)*(b.x-c.x) >= 0) 
		{
			if ((point.x-c.x)*(c.y-a.y) - (point.y-c.y)*(c.x-a.x) >= 0) 
				return true;

		}
	}	  

	return false;
}


void Away(Unit* b1, Unit* b2,float dist)
{  
	if (fabs(b1->x - b2->x) > dist)
		return;

	if (fabs(b1->y - b2->y) > dist)
		return;
 
	float d = Dist(b1, b2);
///	float k = (b1->radius / b2->radius) ;
  
	float antipath_dist = fabs(d - dist)  * 0.5f ; 
 
  
	float vec_x = (b2->x - b1->x) / d;
	float vec_y = (b2->y - b1->y) / d;
  
	float dx =  vec_x * antipath_dist;
	float dy =  vec_y * antipath_dist; 
    
	if (d < dist) {
		b1->x -= dx;
		b1->y -= dy;
		b2->x += dx;
		b2->y += dy;
	} 
}

bool  AwayUnitOne(Unit *u1, Unit *u2,float dist)
{	
	bool res;
	if(Dist(u1,u2)-(dist)<0) res=true;
	else  res=false;

	for(int c=0;Dist(u1,u2)-(dist)<0 && c<1000;c++)
	{ 	
			u1->x+= 0.02f*(u1->x - u2->x);
			u1->y+= 0.02f*(u1->y - u2->y);
	}
	return res;
}

void  AwayM(Unit *b1, Unit *b2,float dist,int k)
{
	for(int c=0;Dist(b1,b2)-(dist)<0 && c<k;c++)
	{ 
		b1->x+= 0.02f*(b1->x - b2->x);
		b1->y+= 0.02f*(b1->y - b2->y);
		b2->x+= 0.02f*(b2->x - b1->x);	
		b2->y+= 0.02f*(b2->y - b1->y);
	}
}

bool  Away(Unit *b1, PointF *p,float dist)
{
	bool res;
	if(Dist(b1,p)-(dist)<=0) res=true;
	else  res=false;
		
	for(int c=0;Dist(b1,p)-(dist)<=0 && c<1000;c++)
	{ 
		b1->x+= 0.01f*(b1->x - p->x);
		b1->y+= 0.01f*(b1->y - p->y);
				
	}
	return res;
	
}



void SVIAZKA(Unit *u1,Unit *u2, float dist,int K)
{

	if(u1->spirit || u2->spirit)
		return ;

	for(int c=0;Dist(u1,u2)-(dist)>0 && c<K;c++)
	{

		u1->y-= 0.001*(u1->y - u2->y);
		u1->x-= 0.001*(u1->x - u2->x);
		u2->y-= 0.001*(u2->y - u1->y);
		u2->x-= 0.001*(u2->x - u1->x);	
	}		

	for(int c=0;Dist(u1,u2)-(dist)<0 && c<1000;c++)
	{

		u1->y+= 0.001*(u1->y - u2->y);
		u1->x+= 0.001*(u1->x - u2->x);
		u2->y+= 0.001*(u2->y - u1->y);
		u2->x+= 0.001*(u2->x - u1->x);							
	}				

}


void S_VIAZKA(Unit *u1,Unit *u2, float dist,int K)
{	

	if(u1->spirit || u2->spirit)
		return ;


	int c;
	for( c=0;Dist(u1,u2)-(dist)>0 && c<K;c++)
	{

		u1->y-= 0.01*(u1->y - u2->y);
		u1->x-= 0.01*(u1->x - u2->x);
		u2->y-= 0.01*(u2->y - u1->y);
		u2->x-= 0.01*(u2->x - u1->x);	
	}	


	for( c=0;Dist(u1,u2)-(dist)<0 && c<K;c++)
	{

		u1->y+= 0.01*(u1->y - u2->y);
		u1->x+= 0.01*(u1->x - u2->x);
		u2->y+= 0.01*(u2->y - u1->y);
		u2->x+= 0.01*(u2->x - u1->x);							
	}




}

void SVIAZKA_UnitOne(Unit *u1,Unit *u2, float dist,int K)
{

	if(u1->spirit || u2->spirit)
		return ;

	for(int c=0;Dist(u1,u2)-(dist)>0 && c<K;c++)
	{

		u1->y-= 0.01*(u1->y - u2->y);
		u1->x-= 0.01*(u1->x - u2->x);

	}		

	for(int c=0;Dist(u1,u2)-(dist)<0 && c<1000;c++)
	{

		u1->y+= 0.01*(u1->y - u2->y);
		u1->x+= 0.01*(u1->x - u2->x);						
	}				

}

void TROS(Unit **u,int n1,int n2,float dist,int k)
{

	for(int i=n1;i<n2;i++)			
		SVIAZKA(u[i],u[i+1], dist ,k);	

	for(int i=n1+1;i<n2+1;i++)			
		SVIAZKA(u[i],u[i-1], dist ,k);	

}



void TROS(vector<RopeElement> *ball,float dist,int k)
{
	vector<RopeElement>::iterator pBall;
	for(pBall=ball->begin(); pBall!=ball->end()-1; pBall++)
	{

		S_VIAZKA(&(*pBall),&(*(pBall+1)), dist ,k);	
	}

}



void PROCH(Unit *u1,Unit *u2, float dist,int K)
{	

	for(int c=0; Dist(u1,u2)-(dist)<0 && c<K;c++)
	{


		u1->y+=  0.02*(u1->y - u2->y)/Dist(u1,u2);
		u1->x+=  0.02*(u1->x - u2->x)/Dist(u1,u2);	
		u2->y+=  0.02*(u2->y - u1->y)/Dist(u1,u2);
		u2->x+=  0.02*(u2->x - u1->x)/Dist(u1,u2);			


	}

}


void  Away(Unit *b, Line *line)
{
	if(Dist(b,&line->a)<b->radius)  Away(b,&line->a,b->radius);
	if(Dist(b,&line->b)<b->radius)  Away(b,&line->b,b->radius);

	/*	float dd1;
	float dd2 ;		

	dd1 = distance (line->a,
	line->b,
	ToPointF(b->x,b->y));

	dd2 = distance (line->a,
	line->b,
	ToPointF(b->oldx,b->oldy));*/	


	/*	if(SIGN(dd1)!= SIGN(dd2) 
	&& CircleIntersects(&b->x,&b->y,b->radius, line->Length, line->a,line->b))

	{
	b->x-=b->dx;
	b->y-=b->dy;

	}*/


	//	for(int o=0;o<10;o++)
	while(CircleIntersects(&b->x,&b->y,b->radius, line->Length, line->a,line->b))
	{
		b->x+=0.1f * line->normal().x ;
		b->y+=0.1f * line->normal().y ;	
	}
}	

int  Collision(Unit *u, TrianglePlatform *platform)
{

	int result = -1;

	if(u->spirit)
		return -1;

	int c=0;
	int numNL = 0;
	if(platform == NULL)
	{	//MessageBox(NULL,"platform == NULL","platform == NULL",MB_OK);
		return -1;
	}
	//	cout<<platform->NUM<<' ';
	if(platform->type < 0 || platform->type >4 )
	{	//MessageBox(NULL,"platform->type < 0 || platform->type >3","platform->type < 0 || platform->type >3",MB_OK);
		return -1;
	}

	//	int flag=0;
	switch (platform->type)
	{

	case GEO:
	case ICE:			

		if(platform->inBox(u)) 
		{
			numNL = platform->numNearLine(ToPointF(u->x,u->y));
			//	u->onHard = false;
			if(platform->isHalfBox)
			{

				PointF normal = platform->line[numNL].normal();

				u->x += normal.x * 0.5 * SIGN(distance(platform->line[numNL].a,platform->line[numNL].b , ToPointF(u->x,u->y)));
				u->y += normal.y * 0.5 * SIGN(distance(platform->line[numNL].a,platform->line[numNL].b , ToPointF(u->x,u->y)));

				if(platform->type!=ICE)
				{

					if(u->friction)
					{	u->dx = 0;
					u->dy = 0;
					u->Move();
					}
				}
				u->onHard = true;

			}

			//	u->onHard = false;
			for(c=0; platform->inBox(u) && pointInPlatform(ToPointF(u->x,u->y ),platform) ;)	
			{ 

				u->x += platform->line[numNL].normal().x * 0.5 * SIGN(distance(platform->line[numNL].a,platform->line[numNL].b , ToPointF(u->x,u->y)));
				u->y += platform->line[numNL].normal().y * 0.5 * SIGN(distance(platform->line[numNL].a,platform->line[numNL].b , ToPointF(u->x,u->y)));

				if(platform->type!=ICE)
				{
					if(u->friction)
					{
						u->dx = 0;
						u->dy = 0;
						u->Move();
					}
				}

				u->onHard = true;				
			}

			/*	if(platform->inBox(u) && platform->isHalfBox  && platform->type==ICE)
			{
			if((u->y - u->dy) < platform->boxMAX.y 
			&&  fabs(u->y - platform->boxMAX.y ) < fabs(u->y - platform->boxMIN.y ))
			u->y =platform->boxMAX.y;


			}*/

		}

		if(u->radius > 0)
		{

			for( c=0;c<3;c++)
			{


				Away(u,&platform->line[c].a,u->radius);
				Away(u,&platform->line[c].b,u->radius);

				while(CircleIntersects(&u->x,&u->y,u->radius, 
					platform->line[c].Length, 
					platform->line[c].a,
					platform->line[c].b))
				{
					u->x+=0.1f * platform->line[c].normal().x ;
					u->y+=0.1f * platform->line[c].normal().y ;
					u->onHard = true;	
					result = platform->NUM;
					u->onHard = true;

					if(platform->type==ICE || !u->friction) {}
					else
					{

						u->dx = 0;
						u->dy = 0;
						u->Move();


					}

				}



			}

		}


		//	Away(u, &platform->line[numNL]);			

		break;
	case SAND:
		if( platform->inBox(u))
		{	if(pointInPlatform(ToPointF(u->x,u->y ),platform)	)
		{ 
			u->dx = 0;
			u->dy = 0;
			u->Move();
			result = platform->NUM;
		}
		}
		break;
	case HYDRO:
	case LAVA:
		if(platform->inBox(u)) 
		{
			if(platform->isHalfBox  )
			{	 
				u->inWater = true;
				result = platform->NUM;
				break;
			}

			if( pointInPlatform(ToPointF(u->x,u->y),platform))
			{	 
				u->inWater = true;
				result = platform->NUM;
			}

		}

		break;
	default: 

		result = -1;
		break;

	}


	return result;		
}

enum ImpactType 
{BALL,
LINE,
_MAX_IMPACT_TYPE_
};

class Impact: public Unit
{
public:
	float H;
	int type;
	float ammoDist; 
	Impact()
	{	
	}

	void set(int type)
	{
		H=100;
		this->type = type;
		radius = 0;
		ammoDist = 2000;
		typeOfImpact = EVIL;
	}


};


enum WeaponTypes{
	LAMP			,
	PISTOL			,
	GUN				,
	AUTOMAT			,
	GRENADE			,
	MACHINE_GUN		,
	BLASTER			,
	ROCKET_LAUNCHER ,
	SUPER_WEAPON	,
	MAXWEAPONTYPES	
};

class Weapon: public Unit
{
public:
	Impact impact[10];
	int impactCounter;
	int type;
	float angle;
	bool fire,re;
	float impactSpeed;
	float impactVectorX,impactVectorY;
	float health;
	bool active;
	void set(int type)
	{
		this->type = type;
		//	setPos(x, y);
		impactCounter = 0;
		angle = 0;
		fire = re = false;
		impactVectorX = 0;
		impactVectorY = 1;
		impactSpeed = 0;

		radius=20;
		health = 100;

		setInertion(0.99);


		active=true;
		switch (type)
		{

		case LAMP:	
			radius=15;
			setFriction(false);
			setMassive(true);
			setInterval(0, 50);
			break;

		case PISTOL:
			setInterval(0, 500);
			impactSpeed   = 2000;
			break;

		case GUN:
			setInterval(0, 900);
			setInterval(2, 900);
			impactSpeed   = 2000;
			break;

		case AUTOMAT:
			setInterval(0, 100);
			setInterval(1, 50);
			impactSpeed   = 2000;
			break;

		case MACHINE_GUN:
			setInterval(0, 50);
			setInterval(1, 25);
			impactSpeed   = 2000;
			break;


		case GRENADE:
			setInterval(0, 500);
			impactSpeed   = 12.5;
			break;

		case BLASTER:
			setInterval(0, 500);
			impactSpeed   = 10;
			break;
		}

		for(int i=0;i<10; i++)
			impact[0].H = 0;
		setFriction(false);

		restart();
	}



	void startFire()
	{
		//	restart();
		fire = true;
		actionOnTime(0);

		for(int i=0;i<10; i++)
			impact[0].setPos(x,y);
	}
	void endFire()
	{
		fire = false;
	}

	void setImpactVector(float vx,float vy)
	{
		impactVectorX = vx;
		impactVectorY = vy;
	}


	/*	float duloLength = 80;
	float ammorad = 2000;//mygish->impact[0].ammoDist;

	float dist = Dist(mygish->medium(),ToPointF(crossX,crossY));
	float k =	dist / ammorad ;
	float k2 =  dist / duloLength ;


	float ammox = mygish->medium().x + ( crossX - mygish->medium().x ) / k;
	float ammoy = mygish->medium().y + ( crossY - mygish->medium().y ) / k;

	float dulox = mygish->medium().x + ( crossX - mygish->medium().x ) / k2;
	float duloy = mygish->medium().y + ( crossY - mygish->medium().y ) / k2;
	*/


	void actionOnTime(int numInterval)
	{

		if(health < 0)
			health = 0;

		if(type == LAMP  && numInterval == 0)
		{
			//if(!inWater) 
			y-=0.05;

		}

		if(re)
		{
			switch (type)
			{

			case GUN:
				switch(numInterval)
				{
				case 2:
					re = false; 
					break;

				}
				break;
			}

		}


		if(fire)
		{
			switch (type)
			{

			case PISTOL:
			case AUTOMAT:
			case MACHINE_GUN:				
				switch(numInterval)
				{
				case 0:
					impact[0].H = 100;
					impact[0].setPos(x,y);
					impact[1].setPos(x+impactVectorX*impactSpeed,y+impactVectorY*impactSpeed);
					break;

				}
				break;

			case GUN:

				switch(numInterval)
				{
				case 0:
					if(!re)
					{
						impact[0].H = 100;
						impact[0].setPos(x,y);
						impact[1].setPos(x+impactVectorX*impactSpeed,y+impactVectorY*impactSpeed);
						impact[2].setPos(x+impactVectorX*impactSpeed+(minirand()*500-250),y+impactVectorY*impactSpeed+(minirand()*500-250));					
						impact[3].setPos(x+impactVectorX*impactSpeed+(minirand()*500-250),y+impactVectorY*impactSpeed+(minirand()*500-250));
						impact[4].setPos(x+impactVectorX*impactSpeed+(minirand()*500-250),y+impactVectorY*impactSpeed+(minirand()*500-250));

						impact[5].setPos(x+impactVectorX*impactSpeed,y+impactVectorY*impactSpeed);
						impact[6].setPos(x+impactVectorX*impactSpeed+(minirand()*500-250),y+impactVectorY*impactSpeed+(minirand()*500-250));					
						impact[7].setPos(x+impactVectorX*impactSpeed+(minirand()*500-250),y+impactVectorY*impactSpeed+(minirand()*500-250));
						impact[8].setPos(x+impactVectorX*impactSpeed+(minirand()*500-250),y+impactVectorY*impactSpeed+(minirand()*500-250));

						for(int i=1; i<9; i++)impact[i].H = 100;


						re = true;
					}
					break;

				}
				break;

			case GRENADE:

				switch(numInterval)
				{
				case 0:
					impact[impactCounter].setPos(x,y);
					//	impact[impactCounter].setDxDy(cos(angle*_torad)*impactSpeed,sin(angle*_torad)*impactSpeed);
					impact[impactCounter].setDxDy(impactVectorX*impactSpeed,impactVectorY*impactSpeed);						
					impact[impactCounter].Move();
					impact[impactCounter].H=100;
					impactCounter++; if(impactCounter>=10) impactCounter=0;	
					break;

				}
				break;


			case BLASTER:

				switch(numInterval)
				{
				case 0:
					impact[impactCounter].setPos(x,y);
					//	impact[impactCounter].setDxDy(cos(angle*_torad)*impactSpeed,sin(angle*_torad)*impactSpeed);
					impact[impactCounter].setDxDy(impactVectorX*impactSpeed,impactVectorY*impactSpeed);						
					impact[impactCounter].Move();
					impact[impactCounter].H=100;
					impact[impactCounter].setMassive(false);
					impactCounter++; if(impactCounter>=10) impactCounter=0;	
					break;

				}
				break;


			}
		}
	}
};


#define YELLOW  0
#define BLUE	1
#define RED		2
#define GREEN	3

class Door: public TrianglePlatform, public Timers
{
public :
	bool isOpen;
	Door()
	{
		isOpen=false;
		setInterval(0,10);
	}

	void close()
	{
		isOpen=false;
	}

	void open()
	{
		isOpen=true;
	}

	void actionOnTime(int numInterval)
	{
		switch(numInterval)
		{
		case 0:
			if(isOpen)
			{
				if(A.y < C.y - 10) 
					A.y +=10;
				setPos(A,B,C);
				//	MessageBox(NULL,"rty","ty3",MB_OK);
			}

			//	

			break;
		}
	}
};

#define YELLOW  0
#define BLUE	1
#define RED		2
#define GREEN	3

class Key: public Unit
{
public:
	int NUM;
	Key () {}

	void actionOnTime(int numInterval) 
	{
	}

	void set(int num, float x, float y )
	{
		NUM = num;
		radius = 20;
		setPos(x,y);
	}

};

enum ItemTypes 
{
	HEALTH,
	AMMO,
	MAXITEMTYPES
};

class Item:public Unit
{
public:
	Item()
	{
		radius = 5+rand()%15;
	}
	int type;

	void actionOnTime(int intervalMilliSec)
	{

	}

	void setType(int type)
	{
		this->type = type;

	}
};

class Hero: public Timers 
{
public:
	bool LIP   ;
	bool SKOL  ;
	bool TVER ;	
	bool FIRE ;
	RopeElement *elem[100];
	int N ;
	float H ,oldH;	
	float Ox ;
	int Score ; 
	PointF oldMed;
	PointF med;
	PointF boxMIN, boxMAX;
	float color[4];
	float tmpcolor[4];
	bool Red;
	bool DelOx;
	float diameter;
	vector<TrianglePlatform>::iterator p;
	///	Impact impact[1];
	bool SLOMO;

	int currentWeapon;

	Weapon weapon;

	Hero()
	{	DelOx = false;
		Red = false;
		Score= 0;
		N = 24;
		LIP  = false;
		SKOL = false;
		TVER = false;
		FIRE = false;
		SLOMO=false;
		H = 100;
		Ox = 10;
		tmpcolor[0] = color[0] = 1;
		tmpcolor[1] = color[1] = 1;
		tmpcolor[2] = color[2] = 1;
		tmpcolor[3] = color[3] = 1;
		SetPos(limitX/2,limitY/2);
		currentWeapon = 1;		
		//	impact
	}

	void SetPos(float x ,float y)
	{
		for(int i=0;i<N;i++)
			elem[i] = new RopeElement(x+70*sin((float)i*_torad*15),y+70*cos((float)i*_torad*15),6.2,i);

		//		starttime = glutGet(GLUT_ELAPSED_TIME);
	}


	void actionOnTime(int numInterval){};
	
	void calcMedium()
	{
		float XMed=0;
		float YMed=0;	

		int kkk=0;

		for(int i=0;i<N;i++)
		{
			if(elem[i]->visible)
			{
				XMed+=elem[i]->x;
				YMed+=elem[i]->y;
				kkk++;
			}
		}

		XMed/=kkk;
		YMed/=kkk;

		med.x = XMed;
		med.y = YMed;

		boxMIN.x = XMed - 200;
		boxMIN.y = YMed - 200;
		boxMAX.x = XMed + 200;
		boxMAX.y = YMed + 200;
	}

	PointF medium()
	{
		calcMedium();
		return med;
	}

	/*	PointF boxMedium ()
	{

	//	st_medium.x = (boxMIN.x + boxMAX.x) /2;
	//	st_medium.y = (boxMIN.y + boxMAX.y) /2;	
	return medium();
	}*/

	int soot(int a)
	{
		if(a<N/2) 
			return a+N/2;
		else 
			return a-N/2;
	}


	int numMedX()
	{
		int res=0;
		for(int i=1;i<N;i++)
		{
			if(fabs(elem[i]->x -medium().x ) > fabs(elem[i-1]->x -medium().x )) res = i-1;
			else res = i;
		}
		return res;
	}


	void Kaplya()
	{	
		int K=25;
		int L=90;

		if(H>0) 
		{	
			TROS ((Unit**)elem,0,N-1,12.4,1000);	
			SVIAZKA (elem[0],elem[N-1],12.4,1000);
		}

		for(int i=0;i<N;i++)
		{
			if(TVER)
			{K=2000;
			L=100;

			elem[i]->y-=0.1;
			}

			K=180000/SQR(diameter=Dist(elem[i],elem[soot(i)])) ;

			PROCH(elem[i],elem[soot(i)],L+12.4,K);

		}	


	}



	void Slipanie(TrianglePlatform *prim, 
		int NUM_TR , 
		float limitDist,
		bool moveUP,
		bool moveDOWN,		 
		bool moveLEFT,
		bool moveRIGHT)
	{


		if(prim->type ==  HYDRO || prim->type ==  LAVA  || prim->type ==  ICE || prim->type ==  SAND )
			return;

		//	for(int k=0;k<3;k++)
		for(int i=0;i<N;i++)
		{



			int numNL = prim->numNearLine(ToPointF(elem[i]->x,elem[i]->y));

			if(CircleIntersects(
				&elem[i]->x,
				&elem[i]->y,
				elem[i]->radius+limitDist,
				prim->GetLineLength(numNL), 	
				prim->GetLineA(numNL),	
				prim->GetLineB(numNL)))
			{

				//	gravity = 0;	
				float dd = distance (prim->GetLineA(numNL),
					prim->GetLineB(numNL),
					ToPointF(elem[i]->x,elem[i]->y));

				elem[i]->x +=prim->GetLineNormal(numNL).x * SIGN(dd);
				elem[i]->y +=prim->GetLineNormal(numNL).y * SIGN(dd);	

				int N = soot(i);


				if(H>0)
				{
					if(moveUP)		elem[N]->y+=0.30;	
					if(moveDOWN)	elem[N]->y-=0.30;	
				}


			}
			//	else 	gravity = 0.1;	
		}
		//		}

	}

	void AdditionColliz(TrianglePlatform *platform,
		int NUM_TR, 
		bool moveDOWN,
		bool moveUP,
		bool moveLEFT,
		bool moveRIGHT,
		bool jump, 
		int i)
	{	 
		static float f = 1; 
		
		if(platform->type ==  HYDRO )
		{	
			return;
		}

		if(platform->type ==  LAVA )
		{	
			for(int i=0;i<N;i++)
				if(pointInPlatform(ToPointF(elem[i]->x,elem[i]->y) , platform))
				{ 
					this->H -=0.01;
				}
				return;
		}

		if(platform->type ==  GEO 
			|| platform->type ==  ICE 
			||  platform->type ==  SAND)
		{
			//				int M=i;
			int N =soot(i);
			int numNL = platform->numNearLine(
				ToPointF(elem[i]->x,elem[i]->y) );

			//		if(pointInPlatform(ToPointF(elem[i]->x,elem[i]->y) , platform))
			//		if(elem[i]->onHard)
			//		if(distance(platform->line[numNL].a,platform->line[numNL].b , elem[i])<=elem[i]->radius)
			//{ 

			if(H>0 && platform->type !=  ICE)
			{					
				if(moveLEFT)	elem[N]->x-=0.5;	
				if(moveRIGHT)	elem[N]->x+=0.5;		
			}

			//	if(moveUP)		elem[N]->y+=0.15;	
			//	if(moveDOWN)	elem[N]->y-=0.15;
			//	if(platform->type == ICE)
			//	{
			//				if(moveLEFT)	elem[M]->x+=0.7;	
			//				if(moveRIGHT)	elem[M]->x-=0.7;	
			//	}

			if(H > 0 && jump)				
			{	
				f+=6;
				elem[N]->x += f * platform->GetLineNormal(numNL).x ;
				elem[N]->y += f * platform->GetLineNormal(numNL).y ;
			}
					
			f=1;

			//	if(!SKOL)
			if(platform->type != ICE)
			{
				/*	elem[M]->dx = platform->dx;
				elem[M]->dy = platform->dy;
				elem[M]->Move();*/
			}
		//}
		} 
	}


	bool inBox(Unit *u)
	{
		if(    u->x >= boxMIN.x 
			&& u->x <= boxMAX.x
			&& u->y >= boxMIN.y
			&& u->y <= boxMAX.y)
			return true;
		return false;
	}

	/*	void damaging(float value)
	{
	this->H -= value;
	color[0] = 1;
	color[1] = 0;
	color[2] = 0;
	}
	*/
};

void HeroCollision(Hero *mygish,Unit *unit)
{

	for(int i=0;i<mygish->N;i++)	
	{

		if(mygish->LIP && Dist(mygish->elem[i],unit)<=unit->radius+mygish->elem[i]->radius)
		{	if (unit->typeOfImpact == EVIL ) 
			SVIAZKA_UnitOne(mygish->elem[i],unit,unit->radius,5000);
			else
			{	
				if((i%2)==0)
					SVIAZKA(mygish->elem[i],unit,unit->radius,5000);
			}

		}
		else 
		{	
			if (unit->typeOfImpact == EVIL ) 
				AwayUnitOne(	mygish->elem[i], unit, unit->radius);
			else 
				Away(	mygish->elem[i], unit, unit->radius);
		}
	}



	/*	RopeElement *unit2;
	unit2 = new RopeElement(mygish->medium().x, mygish->medium().y,1,0);

	if(Dist(unit2,unit)<unit2->radius)	
	{//	Away(	unit2, unit, 60 );
	unit->y+=45;
	unit->setDxDy(0,0);
	unit->Move();
	}

	delete unit2;*/

}

struct EnemyView
{
	float angleHorizontal;
	float angleVertical;
	float x;
	float y;
};

enum EnemyType 
{JUMPER,
ANTI_JUMPER,
TURRET,
FLY,
_MAX_ENEMY_TYPE_
};

class Enemy: public Unit
{


public:

	int type;
	float dspeed;

	bool moveUp;
	bool moveDown;
	bool moveLeft;
	bool moveRight;
	bool move;
	short rotate;
	float angle;
	float dangle;
	bool ACTUAL;
	Impact impact[10];
	int impactCounter;
	EnemyView enemyView;
	float goalx,goaly;
	float health,maxHealth;
	float impactSpeed;
	bool impactMassive;
	int tmpInterval;
	int respawnInterval;

	bool tmpFireType;
	float oldDistToGoal;

	PointF startPosition;

	bool sleeping;
	Enemy()
	{
	}

	void set(int type,float x, float y)
	{
		this->type = type;
		typeOfImpact = EVIL;
		setPos(x, y);

		startPosition = ToPointF(x,y);
		moveUp = false;
		moveDown = false;
		moveLeft = false;
		moveRight = false;
		move=false;
		rotate = 0;
		angle = rand()%360;
		dangle = 0.1;
		dspeed = 2;
		ACTUAL = false;
		restart();
		impactCounter=0;
		goalx = x+1;
		goaly = y;
		oldDistToGoal=0;

		impactMassive = 0;
		tmpFireType = 0;

		tmpInterval = 100+rand()%700;


		respawnInterval = 60000;

		setType(type);

		sleeping = true;


		maxHealth = health = radius;
		impactSpeed = 5+minirand()*10;
		dx=dy=0;

		//setInterval(3, respawnInterval);


		// minirand()>0.5;

		//	if(!impactMassive)
		//		impactSpeed = 3+minirand()*3;

	}

	void setType(int type)
	{
		switch(type)
		{	case ANTI_JUMPER:
		setMassive(-1);
		//	setFriction(false);
		impactMassive = 1;  
		tmpFireType = 1;
		setInterval(1, 2000+rand()%1000);
		case JUMPER: 

			dspeed = 4;	
			//	dangle = 0.8;
			radius = 20+rand()%20;
			tmpInterval = 200+rand()%1500;
			setInterval(0, tmpInterval);
			//	setFriction(false);	

			break;
		case TURRET: 

			dspeed = 0;	
			dangle = minirand()*5;
			radius = 20+rand()%20;	
			setInterval(0, 500 + rand()%1500);
			setInterval(2, 100);
			setMassive(false);
			setFriction(false);
			//	tmpFireType=true;
			impactMassive = 0;
			break;

		case FLY: 
			dspeed = 4;	
			dangle = minirand()*5;
			rotate = 1;
			radius = 15+rand()%35;
			setInterval(0, 1000+rand()%500);
			setInterval(1, 2000+rand()%1000);
			setMassive(false);	
			setFriction(false);

			impactMassive = 0;  
			tmpFireType = 0;

			break;

		}
	}
	void setLive()
	{
		set(type, x,  y);

	}

	void setDeath()
	{
		health = 0;
		//setInertion(0.97);
		if(!spirit)setMassive(1);
		typeOfImpact = GOOD;
	}

	void actionOnTime(int numInterval)
	{




		//	if(numInterval == 3)
		//	{



		//	}


		if(health <=0)
			return;	

		if(sleeping)
			return;

		switch (type)
		{
		case ANTI_JUMPER:
		case JUMPER:
			switch(numInterval)
			{
			case 0:
				//	setInterval(0, rand()%700);
				//		move = !move;
				//		dspeed = 2+minirand() * 3;


				if(onHard )
				{

					FMove.addImpulse((2*(goalx > x)-1)*2, 0);

					if(fabs(x - goalx) > radius*3)
						FMove.addImpulse(0,  (2*(goaly > y)-1) * (tmpInterval/100));

					CalcDxDy();	
					Move();
				}

				if(inWater ) 
				{

					if(goaly > y)
						FMove.addImpulse((2*(goalx > x)-1)*2, 2);
					else 
						FMove.addImpulse((2*(goalx > x)-1)*2, -5);
					CalcDxDy();	
					Move();
				}



				/*	float xx=1,gxx=1;
				float currentDistToGoal = Dist(&x,&xx,&goalx,&gxx);
				if(oldDistToGoal < currentDistToGoal)
				setInertion(0.97);	
				else setInertion(1);

				oldDistToGoal = currentDistToGoal;*/
				break;

			}
			break;

		case TURRET:

			switch(numInterval)
			{
				//	case 2:
				////		angle = getAngle(x,y,goalx,goaly);
				//	break;

			case 0:
				rotate = 1 - rotate;
				if(rotate==0)
				{
					impact[impactCounter].setPos(x,y);
					/*

					impact[impactCounter].setDxDy(cos(angle*_torad)*impactSpeed,sin(angle*_torad)*impactSpeed);

					impact[impactCounter].setFriction(!impactMassive);
					impact[impactCounter].setMassive(impactMassive);
					impact[impactCounter].Move();
					impact[impactCounter].H=100;
					impactCounter++; if(impactCounter>=10) impactCounter=0;

					*/


					if(impactMassive)
						angle = getAngle(x,y,goalx,goaly + fabs(goaly - y));
					else
						angle = getAngle(x,y,goalx,goaly);

					impact[impactCounter].setPos(x,y);
					impact[impactCounter].setDxDy(cos(angle*_torad)*impactSpeed,sin(angle*_torad)*impactSpeed);

					impact[impactCounter].setFriction(!impactMassive);
					impact[impactCounter].setMassive(impactMassive);
					impact[impactCounter].Move();
					impact[impactCounter].H=100;
					impactCounter++; if(impactCounter>=10) impactCounter=0;
				}
				break;

			}

			break;

		case FLY:
			switch(numInterval)
			{
			case 0:

				float currentDistToGoal = Dist(&x,&y,&goalx,&goaly);

				if(oldDistToGoal < currentDistToGoal)
					setInertion(0.97);	
				else setInertion(1);


				oldDistToGoal = currentDistToGoal;


				angle = getAngle(x,y,goalx,goaly+radius*0.5);

				//	CalcVector();

				float vimpx = cos(angle*_torad);
				float vimpy = sin(angle*_torad);
				if( (vimpx  + vimpy) <= 2)
				{//	FMove.addImpulse(vimpx*0.05,vimpy*0.05);
					//	CalcDxDy();	
					//	Move();	
				}


				moveRight = moveUp = moveDown = moveLeft = false;

				if(goalx > x) moveRight = true;
				else		moveLeft = true;

				if(goaly > y) moveUp = true;
				else		moveDown = true;


				break;

			}
			break;
		}


		if(numInterval == 1)
		{
			if(tmpFireType)
			{
				if(impactMassive)
				{	angle = getAngle(x,y,goalx,goaly + fabs(goaly - y));

				if(goaly < y && fabs(goalx - x) < radius*3)
					for(impactCounter=0; impactCounter <5; impactCounter++)
					{
						impactSpeed=minirand()*4-2;
						impact[impactCounter].setPos(x,y);
						impact[impactCounter].setDxDy(cos(angle*_torad)*impactSpeed,sin(angle*_torad)*impactSpeed);

						impact[impactCounter].setFriction(!impactMassive);
						impact[impactCounter].setMassive(impactMassive);
						impact[impactCounter].Move();
						impact[impactCounter].H=100;
					}
				}
				else
				{
					angle = getAngle(x,y,goalx,goaly);

					impact[impactCounter].setPos(x,y);
					impact[impactCounter].setDxDy(cos(angle*_torad)*impactSpeed,sin(angle*_torad)*impactSpeed);

					impact[impactCounter].setFriction(!impactMassive);
					impact[impactCounter].setMassive(impactMassive);
					impact[impactCounter].Move();
					impact[impactCounter].H=100;
					impactCounter++; if(impactCounter>=10) impactCounter=0;
				}
			}

		}






	}

	void setGoal(float gx,float gy)
	{

		goalx = gx;
		goaly = gy;

	}

	void Controller()
	{


		/*	if(move)
		{
		if(x < goalx) 
		moveRight = true;
		else 		
		moveLeft = true;

		if(y < goaly) 
		moveUp = true;
		else 		
		moveDown = true;

		}
		else moveRight = moveLeft =moveDown = moveUp = false;

		if(moveUp) y+=dspeed;
		if(moveDown) y-=dspeed;
		if(moveRight) x+=dspeed;
		if(moveLeft) x-=dspeed;
		*/
		//	if(rotate == 1) angle+=dangle;


		if(health >0)
			angle = getAngle(x,y,goalx,goaly);

		if(health <= 0)
		{
			setDeath();	
		}	


	}

	void InvalidateView()
	{
		enemyView.x = x;
		enemyView.y = y;
	}


};

/*
void SetDif(float r,float g ,float b)
{
	GLfloat	diffuse[4]={r,g,b,1};
	glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);

}
*/



void ShowHero(Hero *mygish)
{
	glColor4f(0.5,0.5,0.5,1); 

	for(int i=0;i<mygish->N;i++)
	{
		float xii = (mygish->elem[i]->x + mygish->medium().x)/2; // seredina
		float yii = (mygish->elem[i]->y + mygish->medium().y)/2;

		float xiii = (xii + mygish->elem[i]->x )/2; // krai i seredina
		float yiii = (yii + mygish->elem[i]->y )/2;

		//			float xiiii = (xiii + mygish->elem[i]->x)/2;
		//			float yiiii = (yiii + mygish->elem[i]->y)/2;

		glPushMatrix();	
			//glColor4f(0.5,0.5,0.5,1);
			glColor3f(0.5,0.5*mygish->H/100,0.5*mygish->H/100 );

			if(mygish->H<=0)
			glColor4f(0.5,0,0,1);

			glTranslatef(xiii,yiii,0);

			if(mygish->LIP )
			{
				if(i%2)	{ 
					glColor3f(0,1,0); 	
					glPushMatrix(); 
					glTranslatef(0,0,10); 
					glRotatef(90+getAngle(
					mygish->elem[i]->x,
					mygish->elem[i]->y,
					mygish->medium().x,
					mygish->medium().y) 
					,0,0,1); 
					glRotatef(90,1,0,0);
					glTranslatef(0,0,-18); 
					glutSolidCone(11,10,8,8); 
					glPopMatrix(); 
				}
			}  
			else {
				glutSolidSphere(18,7,7); 
			}

			if(mygish->TVER)
			{
				glColor3f(0.5,mygish->H/200,mygish->H/200 );
				glPushMatrix(); 
					glRotatef(-90+getAngle(
						mygish->elem[i]->x,
						mygish->elem[i]->y,
						mygish->medium().x,
						mygish->medium().y) 
						,0,0,1); 
					glRotatef(90,1,0,0);
					glutSolidCone(5,35,8,8);
				glPopMatrix(); 
			} 

			glColor3f(1,mygish->H/100,mygish->H/100 );
			glTranslatef(0,0,10);
			glutSolidSphere(11,8,8);
		glPopMatrix();
	}
}


void ShowEnemy(Enemy *enemy)
{
	//	glDisable(GL_DEPTH_TEST);

	if(enemy->moveLeft)
	{
		if(enemy->enemyView.angleHorizontal>-70)  enemy->enemyView.angleHorizontal -= 3;
	}
	if (enemy->moveRight)
	{
		if(enemy->enemyView.angleHorizontal<70)  enemy->enemyView.angleHorizontal += 3;		
	}

	if(enemy->moveUp)
	{
		if(enemy->enemyView.angleVertical>-20)  enemy->enemyView.angleVertical -= 2;
	}
	if (enemy->moveDown)
	{
		if(enemy->enemyView.angleVertical<20)  enemy->enemyView.angleVertical += 2;		
	}

	if(!enemy->moveLeft 
		&& !enemy->moveRight) 
	{
		if(enemy->enemyView.angleHorizontal<0)  enemy->enemyView.angleHorizontal += 1;	
		if(enemy->enemyView.angleHorizontal>0)  enemy->enemyView.angleHorizontal -= 1;	
	}

	if( !enemy->moveDown
		&& !enemy->moveUp) 
	{
		if(enemy->enemyView.angleVertical<0)  enemy->enemyView.angleVertical += 1;	
		if(enemy->enemyView.angleVertical>0)  enemy->enemyView.angleVertical -= 1;	
	}


	glPushMatrix();


	//	glDisable(GL_DEPTH_TEST);
	glTranslatef(enemy->x,enemy->y, 0);
	if(enemy->health <= 0 )
		glTranslatef(0,-enemy->radius*0.1, 0);	

	glRotatef(enemy->enemyView.angleHorizontal,0,1,0);
	glRotatef(enemy->enemyView.angleVertical,1,0,0);

	if(enemy->health <= 0 )		
		glScalef(1,0.95,1);

	float colorLive[3]  = {1,enemy->health/enemy->maxHealth,enemy->health/enemy->maxHealth};

	if(enemy->health<=0) colorLive[0] = 0.33;

	if(enemy->spirit) 
	{
		glDisable(GL_DEPTH_TEST);
		float colorLive1[4]  = {1,1,1,0.7};
		glColor4fv(colorLive1);
		glutSolidSphere(enemy->radius,16,16);
		glEnable(GL_DEPTH_TEST);
	}
	else
	{
		switch(enemy->type)
		{

		case ANTI_JUMPER:
			glRotatef(180,1,0,0);
		case JUMPER:
			glColor3fv(colorLive);
			glRotatef(-90,1,0,0);
			glTranslatef(0,1, -enemy->radius);
			glutSolidCone(enemy->radius, enemy->radius*2, 8, 8);
			glRotatef(-180,1,0,0);
			glutSolidCone(enemy->radius, 0, 8, 8);
			break;
		case FLY:
			glColor3fv(colorLive);
			glutSolidSphere(enemy->radius,16,16);

			//	glTranslatef(0,0, 8);
			//	glColor3f(0,0.5,1);
			//	glutSolidSphere(39,16,16);

		 
			pushLighting();

				glTranslatef(0,0, enemy->radius-enemy->radius*0.5);
				glColor3f(cos( enemy->angle * _torad),0,0);
				glutSolidSphere(enemy->radius*0.6,16,16);


				glTranslatef(0,0, enemy->radius*0.35);
				glColor3f(sin( enemy->angle * _torad),0,0);
				glutSolidSphere(enemy->radius*0.3,16,16);			
				 
			popLighting();
 
			 
			break;

		case TURRET:

			glColor3fv(colorLive);

			glRotatef(enemy->angle,0,0,1);
			glRotatef(90,0,1,0);
			//	glCullFace(GL_FRONT_AND_BACK);

			//	glScalef(45,45,45);
			//	glutSolidCone(enemy->radius,enemy->radius*2,10,10);
			glutSolidSphere(enemy->radius,16,16);
			glRotatef(180,0,1,0);
			glTranslatef(0,0,-enemy->radius*1.5);	
			glutSolidCone(enemy->radius*0.4,enemy->radius*2,10,10);
			glutSolidTorus(enemy->radius*0.05,enemy->radius*0.4,10,10);
			break;
		default:
			glColor3fv(colorLive);
			glutSolidSphere(enemy->radius,16,16);
			break;

		}

	}
	//	glEnable(GL_DEPTH_TEST);
	glPopMatrix();



	//	glEnable(GL_DEPTH_TEST);	
}


void ShowEnemyImpacts(Enemy *enemy)
{
	static float scale = 10;
	static bool scaleUp = true;
	if(scaleUp) scale+=0.5;
	else scale-=0.5;
	if(scale > 10) scaleUp=false;
	if(scale < 4) scaleUp=true;
	glColor3f(1,1,0);
	for(int i=0;i<10;i++)
	{	
		if(enemy->impact[i].H > 0)
		{	
			pushLighting();
				if(!enemy->impact[i].isMassive) {	glColor3f(1,0,0);}
				glPushMatrix();	
				glTranslatef(enemy->impact[i].x,enemy->impact[i].y, 0);
				glRotatef(45,1,1,1);
				glScalef(scale,scale,scale);
				if(!enemy->impact[i].isMassive)
					glutSolidSphere(1,8,8);
				else glutSolidCube(1);
				glPopMatrix();
				//	if(!enemy->impact[i].isMassive) {	 }
			 
			popLighting();
		}
	}
}
 


void ShowLimit()
{
	/*	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glPushMatrix();
	glColor3f(.5,1,.5);
	glBegin(GL_QUADS);
	glNormal3f(1,0,0);
	for(int j=-30;j<150;j+=30)		
	for(int i=0;i<limitY;i+=30)
	{
	glVertex3f(0,i,  j);
	glVertex3f(0,i+30,j);
	glVertex3f(0,i+30,j+30);
	glVertex3f(0,i,j+30);
	}
	glEnd();
	glPopMatrix();

	glCullFace(GL_FRONT);
	glPushMatrix();
	glBegin(GL_QUADS);
	glNormal3f(-1,0,0);
	for( j=-30;j<150;j+=30)		
	for(int i=0;i<limitY;i+=30)
	{
	glVertex3f(limitX,i,  j);
	glVertex3f(limitX,i+30,j);
	glVertex3f(limitX,i+30,j+30);
	glVertex3f(limitX,i,j+30);
	}
	glEnd();
	glPopMatrix();*/

}

/*
void ShowShadowsStena()
{

int	NUM_PRIM_GR_ACTUAL = n_prim_graph_actual.size();

glPushMatrix();

for(int n = 0; n < NUM_PRIM_GR_ACTUAL; n++)
//	if(	prim[n_prim_graph_actual[n]].ACTUAL)
for(int j=0;j<3;j++)
{	

prim[n_prim_graph_actual[n]].line[j]->ShowShadow(0,light);
}
glPopMatrix();
}
*/
/*
void ShowShadowsPlatform(TrianglePlatform* platform)
{
	glPushMatrix();

	if(	platform->ACTUAL)
		for(int j=0;j<3;j++)
		{	
			platform->line[j]->ShowShadow(0,light);
		}
	glPopMatrix();
}
*/

float dddx[10000];			
float dddy[10000];
bool lx[10000];			
bool ly[10000];
TrianglePlatform		prim[MAXPRIMCOUNT];



void ShowDoor(Door *door)
{
	switch (door->NUM)
	{
	case YELLOW:  
		glColor3f(1,1,0);
		break;
	case BLUE:
		glColor3f(0,0,1);
		break;
	case RED:
		glColor3f(1,0,0);
		break;
	case GREEN:
		glColor3f(0,1,0);
		break;
	default:
		glColor3f(1,1,1);
		break;
	}


	door->line[0].ShowLite();
	door->line[1].ShowLite();
	door->line[2].ShowLite();

	tr(	door->A.x,door->A.y,door->Width,   
		door->C.x,door->C.y,door->Width,		   
		door->B.x,door->B.y,door->Width);
}

void ShowKey(Key *key)
{
	switch (key->NUM)
	{
	case YELLOW:  
		glColor3f(1,1,0);
		break;
	case BLUE:
		glColor3f(0,0,1);
		break;
	case RED:
		glColor3f(1,0,0);
		break;
	case GREEN:
		glColor3f(0,1,0);
		break;
	default:
		glColor3f(1,1,1);
		break;
	}

	glPushMatrix();
	glTranslatef(key->x, key->y, 0);
	//	glScalef(key->radius*1.2,key->radius*1.2,key->radius*1.2);
	glutSolidTorus(key->radius*0.3,key->radius-key->radius*0.3,16,16);

	glPopMatrix();

}


void ShowPlatform(TrianglePlatform *platform)
{/*
	if(sector[curDescI][curDescJ].num1 == platform->NUM
 	 || sector[curDescI][curDescJ].num2 == platform->NUM)
 	 	 return;*/

	float mx, my;

	mx= platform->medx ;
	my= platform->medy ;

//	platform->view.color[3]=1;

	switch (platform->type)
	{
	case SAND:
	 
		glColor4fv(platform->view.color);

		platform->line[0].ShowLite();
		platform->line[1].ShowLite();

		tr(platform->A.x,platform->A.y,platform->Width,
			platform->C.x,platform->C.y,platform->Width,
		platform->B.x,platform->B.y,platform->Width);
		break;

	case ICE:
 
		glColor4fv(platform->view.color);
		platform->line[0].ShowLite();
		platform->line[1].ShowLite();
		platform->line[2].ShowLite();
		tr(platform->A.x,platform->A.y,platform->Width,
		platform->C.x,platform->C.y,platform->Width,
		platform->B.x,platform->B.y,platform->Width);
		break;

	case GEO:	
 
		glColor4fv(platform->view.color); 

		if(  !platform->isDepth )
		{
			platform->line[0].ShowLite();
			platform->line[1].ShowLite();

			if(!platform->isHalfBox)
				platform->line[2].ShowLite();
		}
 
		float a;
		if((platform->NUM%3) == 0)
			a=0.95;
		else
			a=1.07;
			
		tr( mx,my,platform->Width*a,
		platform->B.x,platform->B.y,platform->Width,
		platform->A.x,platform->A.y,platform->Width);
 
		tr(platform->C.x,platform->C.y,platform->Width,
		platform->B.x,platform->B.y,platform->Width,
		mx,my,platform->Width*a);
 

		tr(platform->A.x,platform->A.y,platform->Width,
		platform->C.x,platform->C.y,platform->Width,
		mx,my,platform->Width*a);
		break;

	case HYDRO: 
	case LAVA:
		platform->view.color[3] = 0.6;
		glColor4fv(platform->view.color);

		if(dddx[platform->NUM] > 20 ) lx[platform->NUM] = false;
		if(dddx[platform->NUM] < -20 ) lx[platform->NUM] = true;

		if(dddy[platform->NUM] > 20 ) ly[platform->NUM] = false;
		if(dddy[platform->NUM] < -20 ) ly[platform->NUM] = true;

		if(lx[platform->NUM]) dddx[platform->NUM]+=0.5;
		else	dddx[platform->NUM]-=0.5;

		if(ly[platform->NUM]) dddy[platform->NUM]+=0.5;
		else    dddy[platform->NUM]-=0.5;	


		float mx= platform->medx + dddx[platform->NUM];
		//		float my= platform->medy + dddy[platform->NUM];
		float mz= platform->medz + dddy[platform->NUM];		
 
		tr(platform->A.x,platform->A.y,platform->Width-5,
		platform->C.x,platform->C.y,platform->Width-5,		 
		platform->B.x, platform->B.y,platform->Width-5);

		if(platform->uppest)
		{
			mz=0;
		 

			tr(platform->A.x,platform->C.y,platform->Width-5,
			platform->C.x,platform->C.y,platform->Width-5,		 
			mx,platform->C.y+dddy[platform->NUM],mz);
			 
			tr(platform->A.x,platform->C.y,platform->Width-5,
			mx,platform->C.y+dddy[platform->NUM],mz,
			platform->A.x,platform->C.y,-platform->Width	 
			);
		 
			tr(platform->C.x,platform->C.y,platform->Width-5,
			platform->A.x,platform->C.y,-platform->Width, 
			mx,platform->C.y+dddy[platform->NUM],mz	 	 
			);
 
			tr(platform->C.x,platform->C.y,-platform->Width,
			platform->C.x,platform->C.y,platform->Width-5,		 
			mx+50,platform->C.y+dddx[platform->NUM],mz);

			 
			tr(platform->C.x,platform->C.y,-platform->Width,
			mx+50,platform->C.y+dddx[platform->NUM],mz,
			platform->A.x,platform->C.y,-platform->Width	 
			);
	 
			tr(platform->C.x,platform->C.y,platform->Width-5,
			platform->A.x,platform->C.y,-platform->Width, 
			mx+50,platform->C.y+dddx[platform->NUM],mz	 	 
			);	
		}
//		tr(platform->A.x,platform->A.y,platform->Width-5,   platform->C.x,platform->C.y,platform->Width-5,		  platform->B.x, platform->B.y,platform->Width-5);

		break;
	}
}


void ShowMapBackground(int iii, int jjj, int area)
{
	glBegin(GL_TRIANGLES);
	glColor4f(.2,.2,.2,0.7);
	//	glColor3f(.5, 1,.5);
	int z =  -700; //- prim[0].Width;
	for (int n = iii-area; n<=iii+area; n++) 
	{
		for(int m = jjj-area; m<=jjj+area; m++)
		{  
			//	  if(n >=0 && m>=0 && n<100 && m<100)
			//	  if(prim[sector[n][m].num1].type != HYDRO)
 
			PointF a = ToPointF(n*SIZE1,(100-m)*SIZE2);
			PointF b = ToPointF(n*SIZE1+SIZE1,(100-m)*SIZE2);
			PointF c = ToPointF(n*SIZE1,(100-m)*SIZE2+SIZE2); 
			PointF d = ToPointF(n*SIZE1+SIZE1,(100-m)*SIZE2+SIZE2);  

			float bmx =	(a.x + d.x)*0.5;
			float bmy = (a.y + d.y)*0.5;

			tr(b.x,b.y,z,  
				bmx,bmy,z,  
				a.x,a.y,z);

			tr(c.x,c.y,z,
				bmx,bmy,z,  
				d.x,d.y,z);
		}
	}
	glEnd(); 
}

void ShowCheckPoint(PointF *checkPoint)
{
	glPushMatrix();
	glColor3f(1,1,0);
	glTranslatef(checkPoint->x, checkPoint->y, 0);
	glScalef(30,45,30);
	glutSolidOctahedron();
	glPopMatrix();

}


void ShowExit(PointF *exitPoint)
{
	glPushMatrix();
	pushLighting();
		glColor3f(0.2,0,0);
		glTranslatef(exitPoint->x, exitPoint->y, 0);
		glutSolidCube(120);
	popLighting();
	glPopMatrix();

}

void ShowImpactLine(Hero *mygish)
{
	//	float maxx =  MAX(fabs(mygish->medium().x - mygish->impactx), fabs(mygish->medium().x - mygish->impactx) );
	//	float minn =  MIN(fabs(mygish->medium().x - mygish->impactx), fabs(mygish->medium().x - mygish->impactx) );


	float duloLength = 80;
	float ammorad = 2000;//mygish->impact[0].ammoDist;

	float dist = Dist(mygish->medium(),ToPointF(crossX,crossY));
	float k =	dist / ammorad ;
	float k2 =  dist / duloLength ;


	float ammox = mygish->medium().x + ( crossX - mygish->medium().x ) / k;
	float ammoy = mygish->medium().y + ( crossY - mygish->medium().y ) / k;

	float dulox = mygish->medium().x + ( crossX - mygish->medium().x ) / k2;
	float duloy = mygish->medium().y + ( crossY - mygish->medium().y ) / k2;

	//	mygish->impact[0].x = ammox;
	//	mygish->impact[0].y = ammoy;

	static int c = 0;
	static int c2 = 0;
	//	static int c3 = 0;
	//	static int c4 = 0;

	glColor3f(1,0,0);
	glBegin(GL_LINES);
		glColor3f(0,0,0);	
		glVertex3f(mygish->medium().x, mygish->medium().y 	,50);
		glVertex3f(dulox, 	duloy,50);
	glEnd();

	if(mygish->FIRE)
	{

		if(c++%2 == 0) return;
		if(c2++%2 == 0) return;
		//	if(c3++%2 == 0) return;

		glLineWidth(3);

		pushLighting();
			glColor3f(1,0,0);
			//glTranslatef(exitPoint->x, exitPoint->y, 50);
			//	glutSolidCube(120);
			glBegin(GL_LINES);
			//	glColor3f(0.8,0.8,0.8);
			glColor3f(minirand(),minirand(),minirand());
			glVertex3f(mygish->medium().x, mygish->medium().y 	,40);
			glVertex3f(ammox, 	ammoy,40);
			glEnd();
		popLighting();
	}
}

void ShowWeaponWithImpacts(Weapon *weapon)
{		
	if(weapon->type == ROCKET_LAUNCHER)
		return;

	switch(weapon->type)
	{
	case 0: 
		glColor3f(1,1,1); 
		glPushMatrix();
		glTranslatef(weapon->x, weapon->y, 0);				
		glutSolidSphere(weapon->radius,16,16);
		glPopMatrix();
		return;
		break;
	case 1: 	glColor3f(.3,.3,0);
		break;
	case 2: 	glColor3f(0,.3,0);
		break;
	case 3: 	glColor3f(0,.3,.3);
		break;
	case 4: 	glColor3f(0,0,1);
		break;
	case 5: 	glColor3f(.3,0,.3);
		break;
	case 6: 	glColor3f(.3,.5,.3);
		break;			
	case 7: 	glColor3f(.5	,.3,.5);
		break;
	case 8: 	glColor3f(.5,.5,.5);
		break;	
	}

	glPushMatrix();
	glTranslatef(weapon->x, weapon->y, 0);
	glScalef(weapon->radius*1.2,weapon->radius*1.2,weapon->radius*1.2);
	glutSolidIcosahedron();
	glPopMatrix();

	for(int i=0;i<10;i++)
	{	
		if(weapon->impact[i].H > 0)
		{
			glPushMatrix();	
			glTranslatef(weapon->impact[i].x,weapon->impact[i].y, 0);
			glRotatef(45,1,1,1);
			glScalef(10,10,10);
			glutSolidCube(1);
			glPopMatrix();
		}
	}
 

	if(weapon->fire && weapon->impactSpeed == 2000)
	{

		//static int c = 0;
		//	static int c2 = 0;
		//	if(c++%2 == 0) return;
		//	if(c2++%2 == 0) return;
		if(weapon->impact[0].H > 0)	
		{
			for(int i=1;i<9;i++)
			{	if(weapon->impact[i].x != 0)	
			{
				glColor3f(0.5,0.5,0.5);
				glBegin(GL_LINES);
					glVertex3f(weapon->impact[0].x, weapon->impact[0].y ,0);
					glVertex3f(weapon->impact[i].x, weapon->impact[i].y	,0);
				glEnd();
			}
			weapon->impact[i].H = 0;

			}
			weapon->impact[0].H = 0;
		}
	}
}

void ShowItem(Item *item)
{

	static float angle = 0 ;
	static float angle2 = 0 ;
	if((angle+=0.1 )>= 360) angle = 0;
	if((angle2+=0.03 )>= 360) angle2 = 0;
	glPushMatrix();
	glTranslatef(item->x, item->y, 0);	
	glRotatef(angle2,0,1,0);
	glRotatef(angle,0,0,1);


	switch(item->type)
	{

	case HEALTH:
		glColor3f(1,0,0);
		glPushMatrix();
		glScalef(3,1,1);
		glutSolidCube(item->radius);
		glPopMatrix();

		glPushMatrix();
		glScalef(1,3,1);
		glutSolidCube(item->radius);
		glPopMatrix();
		break;			
	case AMMO: 	
		break;
	default:
		glutSolidCube(item->radius);
		break;

	}

	glPopMatrix();
}

class LevelHelper
{
	int linesCount;
public:
	int numCurrentLevel;
//	std::string currentLevel;

	LevelHelper()
	{
		numCurrentLevel = 1;
		linesCount = 0;
	}

	std::string nextLevel(std::string scriptFile)
	{
		std::ifstream in(scriptFile.c_str());
		std::string level;
		std::string curLine;		

		if (!in.is_open()) {
			printf("file '%s' not found\n", scriptFile.c_str());
			in.close();
			return "";
		}

	  
		if(numCurrentLevel  > linesCount)
			numCurrentLevel = 1;

		linesCount = 0;
		while (in >> curLine){ 
			linesCount++;

			if(linesCount == numCurrentLevel) { 
				level = curLine;
			}
		}
		

		numCurrentLevel ++; 
		  
		return level;
	}
/*
	std::string getCurrentLevel()
	{

		return currentLevel;

	}
	*/
};


bool onlgt[2];


bool moveUP=false;
bool moveDOWN=false;
bool moveLEFT=false;
bool moveRIGHT=false;


bool moveUP_2=false;
bool moveDOWN_2=false;
bool moveLEFT_2=false;
bool moveRIGHT_2=false;


bool jump=false;
float YROT= 0;
float XROT= 0;
float ZROT= 0;
bool cameraHeroOrientation = true;
char command[100];
bool  playing=true;
//PointF light[2];
PointF MEDGISH;
GLenum format;
GLenum type;




PointF  camera;

int NUMLINES;
int NUM_TR;
int RMX = limitX;
int LMX = 0;
int BMY = limitY;
int TMY = 0;

Hero gish;

Enemy			enemy[MAXENEMYCOUNT];
Weapon			weapon[11];
Door			door[MAXDOORCOUNT];
Key				key[MAXDOORCOUNT];
Item			item[MAXITEMS];

int enemyCount = 0;
int doorCount = 0;
int keyCount = 0;
int itemCount = 0;

vector<RopeElement>		ball;


vector<int> n_prim_graph_actual;
vector<int> n_prim_actual;
//vector<TrianglePlatform>		::iterator iter_prim;
vector<RopeElement>		::iterator pBall;
vector<Enemy>			::iterator iter_enemy;	
vector<int> n_enemy_actual;

PointF Hero;
PointF exitLevel;
PointF startLevel;
PointF respawnPoint;
PointF checkPoint[CHECKPOINT_COUNT];


LevelHelper levelHelper;

float	FVector::gravity = 0.16;

bool isMap = false;

bool noclipMode = false;


float colorGeo[3];
float colorWalk[3];

float lightColor[100][100][3];

bool researched[100][100];
char dat[101][101];



GLfloat	real_lightmodel[4]={0.1,0.1,0.1,1};

PointF forbombPoint;


void GetQuadCoord(int *I, int *J, PointF heroPos) {
	*I = heroPos.x / SIZE1;
	*J = 101 - heroPos.y / SIZE2;

	if (*I < 0 || *I > 99)
		*I = 0;
	if (*J < 0 || *J > 99)
		*J = 0;

}


bool underCameraUnit (PointF unit, PointF camera, int area)
{
	int iii,jjj;
	GetQuadCoord(&iii,&jjj, unit);

	int iii2,jjj2;
	GetQuadCoord(&iii2,&jjj2, camera);

	if(abs(iii-iii2) < area && abs(jjj-jjj2) < area)
		return true;

	return false;
}


bool researchedUnit (PointF unit)
{
	int iii,jjj;
	GetQuadCoord(&iii,&jjj, unit);

	if(iii >=0 && iii<100 && jjj>=0 && jjj<100)
		if(researched[iii][jjj])
			return true;

	return false;
}


void DescCollision(Unit *unit, int *primIndex) {

	int iii, jjj;
	GetQuadCoord(&iii, &jjj, *unit);
	int primi = -1;

	for (int n = iii - 2; n <= iii + 1; n++)
		for (int m = jjj - 1; m <= jjj + 1; m++) {
			if (n >= 0 && m >= 0 && n < 100 && m < 100) {

				if ((primi = Collision(unit, &prim[sector[n][m].num1])) != -1) {
					*primIndex = primi;
				}
				if ((primi = Collision(unit, &prim[sector[n][m].num2])) != -1) {
					*primIndex = primi;
				}
			}
		}

		for (int d = 0; d < doorCount; d++) {
			if (door[d].inBox(unit))
				Collision(unit, &door[d]);
		}

}

void SimpleCollision(Unit *unit) {
	int iii, jjj;
	GetQuadCoord(&iii, &jjj, *unit);
	Collision(unit, &prim[sector[iii][jjj].num1]);
	Collision(unit, &prim[sector[iii][jjj].num2]);
}

void Collisions()
{
	GetQuadCoord(&curDescI,&curDescJ,gish.medium());

	//curDescI = gish.medium().x / SIZE1;
	//curDescJ = 101-gish.medium().y / SIZE2;
	int i;

	if(!noclipMode)
	{

		if(gish.LIP && !gish.SKOL)	
		{


			for (int n = curDescI-1; n<=curDescI+1; n++)
				for(int m = curDescJ-1; m<=curDescJ+1; m++)
				{	
					if(n >=0 && m>=0 && n<100 && m<100)
					{	
						for(int k =0; k<2 ; k++)
						{	gish.Slipanie(&prim[sector[n][m].num1],NUM_TR,6,
						moveUP || moveUP_2 ,
						moveDOWN || moveDOWN_2,
						moveLEFT || moveLEFT_2,
						moveRIGHT || moveRIGHT_2
						);

						gish.Slipanie(&prim[sector[n][m].num2],NUM_TR,6,
							moveUP || moveUP_2 ,
							moveDOWN || moveDOWN_2,
							moveLEFT || moveLEFT_2,
							moveRIGHT || moveRIGHT_2
							);
						}
					}
				}

				for(int k=0;k < MAXDOORCOUNT;k++)
					gish.Slipanie(&door[k],NUM_TR,6.51,
					moveUP || moveUP_2 ,
					moveDOWN || moveDOWN_2,
					moveLEFT || moveLEFT_2,
					moveRIGHT || moveRIGHT_2
					);

		}




		for(  i=0;i<gish.N;i++)	
		{
			gish.elem[i]->inWater = false;
			gish.elem[i]->onHard = false;
			int primi=-1;
			DescCollision(gish.elem[i],&primi);

			if(primi >=0 && primi<MAXPRIMCOUNT)
			{


				gish.AdditionColliz(
					&prim[primi],
					NUM_TR,
					moveUP || moveUP_2 ,
					moveDOWN || moveDOWN_2,
					moveLEFT || moveLEFT_2,
					moveRIGHT || moveRIGHT_2,
					jump,
					i);

			}

		}
		FVector::gravity =DEFAULT_GRAVITY;


	}
	else FVector::gravity = 0;


	//float dx  = -(float)sin((-90+ZROT)*_torad) * 0.05f;
	//	float dy  = (float)cos((-90+ZROT)*_torad) * 0.05f;
	gish.oldMed = gish.medium();

	for( i=0;i<gish.N;i++)
	{

		if(!gish.elem[i]->moveble) 
			gish.elem[i]->CalcVector();




		if(moveDOWN || moveDOWN_2 ) gish.elem[i]->FMove.addImpulse(0, -0.04);
		if(moveUP|| moveUP_2 ) 	 gish.elem[i]->FMove.addImpulse(0, 0.04);


		if(YROT >= 90 &&  YROT < 270)
		{ if(moveRIGHT || moveRIGHT_2 ) gish.elem[i]->FMove.addImpulse(-0.04, 0);	
		if(moveLEFT || moveLEFT_2) gish.elem[i]->FMove.addImpulse(0.04, 0);
		}
		else
		{ if(moveLEFT || moveLEFT_2) gish.elem[i]->FMove.addImpulse(-0.04, 0);	
		if(moveRIGHT || moveRIGHT_2) gish.elem[i]->FMove.addImpulse(0.04, 0);
		}


		//	if(moveDOWN) gish.elem[i]->FMove.addImpulse(-dx, -dy);
		//	if(moveUP) 	 gish.elem[i]->FMove.addImpulse(dx, dy);



		if(!gish.elem[i]->moveble) 
		{	gish.elem[i]->CalcDxDy();	
		gish.elem[i]->Move();
		}

	}	

	gish.Kaplya();
	////////////
	/////////////
	////////////



	size_t ENEMY_ACTUAL_COUNT = n_enemy_actual.size();

	for(int d=0;d < MAXDOORCOUNT;d++) 
	{
		int primi;
		DescCollision(&key[d],&primi);
		HeroCollision(&gish,&key[d]);
		key[d].CalcVector();
		key[d].CalcDxDy();
		key[d].Move();

	}


	for(int k=0;k<11;k++)
	{	weapon[k].inWater = false;
	weapon[k].onHard = false;


	if(gish.currentWeapon != k) 
	{	HeroCollision(&gish,&weapon[k]);
	int primi;
	DescCollision(&weapon[k],&primi);	

	for ( int d  = 0; d < 10; d++)
		if(d!=k) 
			Away(&weapon[k],&weapon[d], weapon[d].radius + weapon[k].radius);



	for(int d=0; d < MAXDOORCOUNT;d++)	
	{
		Away(&weapon[k],&key[d], key[d].radius + weapon[k].radius);
	}


	//	int ENEMY_ACTUAL_COUNT = n_enemy_actual.size();
	for ( size_t j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
	{
		Away(&weapon[k],&enemy[n_enemy_actual[j]], enemy[n_enemy_actual[j]].radius + weapon[k].radius);
	}
	}

	weapon[k].CalcVector();
	weapon[k].CalcDxDy();
	weapon[k].Move();



	for(int i=0; i<10;i++) {
		weapon[k].impact[i].inWater = false;
	}

	for(int i=0; i<10; i++) {
		int iii,jjj;
		GetQuadCoord(&iii,&jjj,	weapon[k].impact[i]);	
		Collision( &weapon[k].impact[i],	&prim[sector[iii][jjj].num1]);
		Collision( &weapon[k].impact[i],	&prim[sector[iii][jjj].num2]);

		for(int d=0; d < MAXDOORCOUNT;d++)	
		{
			Collision( &weapon[k].impact[i],	&door[d]);
		}


		if(weapon[k].impact[i].H)
		{	weapon[k].impact[i].CalcVector();
		weapon[k].impact[i].CalcDxDy();
		weapon[k].impact[i].Move();
		}
	}
	}



	for (size_t  j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
	{	
		Enemy *tmpEnemy;  
		tmpEnemy = &enemy[n_enemy_actual[j]];
		tmpEnemy->inWater = false;
		tmpEnemy->onHard  = false;

		for ( size_t k  = 0; k < ENEMY_ACTUAL_COUNT; k++)
		{	if(n_enemy_actual[j]==n_enemy_actual[k]) continue;
		Enemy *tmpEnemy2; 
		tmpEnemy2 = &enemy[n_enemy_actual[k]];
		Away(	tmpEnemy2, tmpEnemy, tmpEnemy->radius + tmpEnemy2->radius);
		}

		HeroCollision(&gish,tmpEnemy); /////!
		int primi;
		DescCollision(&enemy[n_enemy_actual[j]],&primi);
		enemy[n_enemy_actual[j]].CalcVector();
		enemy[n_enemy_actual[j]].CalcDxDy();
		enemy[n_enemy_actual[j]].Move();
	}



	for (size_t  j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
	{	

		Enemy *tmpEnemy;  
		tmpEnemy = &enemy[n_enemy_actual[j]];


		for(int k=0; k<10;k++)
		{
			tmpEnemy->impact[k].inWater = false;
			tmpEnemy->impact[k].onHard = false;

			SimpleCollision(&tmpEnemy->impact[k]); 

			HeroCollision(&gish,&tmpEnemy->impact[k]); 


			if(tmpEnemy->impact[k].onHard)tmpEnemy->impact[k].H = 0;

			for (size_t  d  = 0; d < ENEMY_ACTUAL_COUNT; d++)
			{
				if(d!=j) 
				{
					if(enemy[n_enemy_actual[d]].health > 0)
						PROCH(&tmpEnemy->impact[k],&enemy[n_enemy_actual[d]], enemy[n_enemy_actual[d]].radius,5000);
				}
			}


			if(tmpEnemy->impact[k].H)
			{	tmpEnemy->impact[k].CalcVector();
			tmpEnemy->impact[k].CalcDxDy();
			tmpEnemy->impact[k].Move();

			if(tmpEnemy->impact[k].onHard )
			{	tmpEnemy->impact[k].H=0;
			tmpEnemy->impact[k].setPos(10000,1000000);
			}

			}

		}		

	}




	//		TROS ((Unit**)ball,0,ball.size()-1,12.4,2000);

	TROS(&ball,20,5000);
	//	SVIAZKA (&ball[0],&ball[ball.size()-1],12.4,2000);





	//		if(moveUP) weapon[GUN].y+=5;
	//		if(moveDOWN) weapon[GUN].y-=5;
	//		if(moveLEFT) weapon[GUN].x-=5;
	//		if(moveRIGHT) weapon[GUN].x+=5;

	//		if(moveUP_2) weapon[SUPER_WEAPON].y+=5;
	//		if(moveDOWN_2) weapon[SUPER_WEAPON].y-=5;
	//		if(moveLEFT_2) weapon[SUPER_WEAPON].x-=5;
	//		if(moveRIGHT_2) weapon[SUPER_WEAPON].x+=5;


	ball[0].x = weapon[SUPER_WEAPON].x;
	ball[0].y =	weapon[SUPER_WEAPON].y;

	//		ball[ball.size()-1].x = weapon[GUN].x;
	//		ball[ball.size()-1].y =	weapon[GUN].y;





	for(pBall=ball.begin(); pBall!=ball.end(); pBall++)
	{


		for (size_t  j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
		{	


			Enemy *tmpEnemy;  
			tmpEnemy = &enemy[n_enemy_actual[j]];
			if(tmpEnemy->health > 0)
				Away(tmpEnemy,&(*pBall), pBall->radius + tmpEnemy->radius);

		}


		int primi;
		DescCollision(&(*pBall),&primi);


		pBall->CalcVector();
		pBall->CalcDxDy();
		pBall->Move();
		/*	for( int i=0;i<gish.N;i++)	
		{

		if(gish.LIP && Dist(gish.elem[i],pBall)<gish.elem[i]->radius+pBall->radius+5)
		{
		SVIAZKA(gish.elem[i],pBall,14,5000);
		}
		else 	Away(pBall,	gish.elem[i],	6+pBall->radius);

		}*/
	}

}



void setDeath()
{

	/*
	glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.9);
	glEnable (GL_LIGHT1);
	glDisable (GL_LIGHT0);
	real_lightmodel[0]=0;
	real_lightmodel[1]=0;
	real_lightmodel[2]=0;
	real_lightmodel[3]=1;
	//	glClearColor(0.1,0,0,1);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,real_lightmodel);
	GLfloat	diffuse[4]={1,0.5,0.5,1};
	glLightfv(GL_LIGHT1,GL_DIFFUSE,diffuse);
	*/
}


void setDay()
{
	glEnable (GL_LIGHTING); 
	glEnable (GL_LIGHT0); 

	GLfloat	diffuse[4]={1,1,1,1};

	GLfloat commonColor[] = {.5, .55, .55, 1};
	
	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);

	glClearColor(commonColor[0],commonColor[1],commonColor[2],commonColor[3]);

	real_lightmodel[0]=0.0;
	real_lightmodel[1]=0.0;
	real_lightmodel[2]=0.0;
	real_lightmodel[3]=1;
//	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,real_lightmodel);

	glEnable(GL_FOG);
	glFogf(GL_FOG_MODE, GL_LINEAR);
	glFogf(GL_FOG_START, 1300.0f);
	glFogf(GL_FOG_END, 2000.0f);
	//	glFogf(GL_FOG_DENSITY, 0.0006f);
	glFogfv(GL_FOG_COLOR, commonColor);
}

void setLocalLight(int number,
				   int x,int y,float z,
				   float difR, float difG, float difB,
				   GLfloat k)
{	
	return;

	glDisable (GL_LIGHT0);
	glLightf(GL_LIGHT0+number, GL_QUADRATIC_ATTENUATION, k);
	glEnable (GL_LIGHT0+number);

	real_lightmodel[0]=0.02;
	real_lightmodel[1]=0.02;
	real_lightmodel[2]=0.02;
	real_lightmodel[3]=1;
	glClearColor(.5,.5,.5,1);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,real_lightmodel);

	//	int n,m;
	//	GetQuadCoord(&n,&m, ToPointF(x,y));

	//	int n = x/SIZE1;
	//	int m = 101 - y/SIZE2;

	//	if(number > 1)
	/*	{
	glLightfv(GL_LIGHT0+number,GL_DIFFUSE,diffuse);
	}
	else{*/
	GLfloat	diffuse[4]={difR,difG,difB,1};
	glLightfv(GL_LIGHT0+number,GL_DIFFUSE,diffuse);
	//	}

	GLfloat myLightPosition[] = {
		x,
		y,	
		z, 
		1
	};

	glLightfv(GL_LIGHT0+number, GL_POSITION, myLightPosition);	
}

int starttime1=0,time_1;
int starttime2=0,time_2;
int starttime3=0,time_3;
int tGT=0;

void LoadMap(char *path);

int deathPause = 0;
void idle()
{ 
	int n;

	tGT = glutGet(GLUT_ELAPSED_TIME);         
	time_1 = tGT-starttime1; 
	time_2 = tGT-starttime2; 
	time_3 = tGT-starttime3;

	if(time_3>1000/100)
	{
		glutPostRedisplay();			
		starttime3+=time_3;
	}

	if(pauseMode)
		return;

	gish.timer();

	for(int i=0; i<10; i++)
		weapon[i].timer();

	for(int d=0;d<MAXDOORCOUNT;d++)
		door[d].timer();

	//	if ( gish.time[0] >= 45000 )  
	//	{
	//		respawnPoint = gish.medium();
	//  		gish.starttime[0] += gish.time[0];	
	//	}

	// int numEnemyActual = n_enemy_actual.size();
	for( n=0; n < enemyCount; n++)
	{ //  if(	enemy[n_enemy_actual[n]].health>0)
		enemy[n].timer();
	}				

	if ( time_2 > 1000/5)  
	{ starttime2 += time_2; 

		/*   int numEnemyActual = n_enemy_actual.size();
		for( n=0; n < numEnemyActual; n++)
		{
		enemy[n_enemy_actual[n]].setGoal(gish.medium().x, gish.medium().y);

		}*/

		if( gish.H<= 0 )
		{
				gish.H = 0 ;
				deathPause++;
				setDeath();
				if(deathPause > 25)
				{	 //glClearColor(0,0,0,1);
					//		     setNight();
					gish.H = 100;
					PointF tmpP = respawnPoint;
				//	LoadMap((char*)levelHelper.getCurrentLevel().c_str());
					gish.SetPos(tmpP.x, tmpP.y);

					deathPause=0;
				}

		}



		if(Dist(gish.medium(), exitLevel) < 50)
			LoadMap((char*)levelHelper.nextLevel("mapscript.dat").c_str());

	for(int d=0; d < MAXDOORCOUNT; d++)
		for(int k=0; k < MAXDOORCOUNT; k++)
			if(key[k].NUM == door[d].NUM)
				if(key[k].x + key[k].radius > door[d].boxMIN.x 					  
					&& key[k].x - key[k].radius < door[d].boxMAX.x 
					&&key[k].y + key[k].radius > door[d].boxMIN.y 					  
					&& key[k].y - key[k].radius < door[d].boxMAX.y  )
				{
					door[d].open();
					key[k].radius = 0;
				}



	for(n = 0;  n<MAXPRIMCOUNT; 	n++)
		prim[n].Move();


	n_enemy_actual.clear();

	for (n = 0; n < enemyCount; n++)
	{	
		enemy[n].ACTUAL = 
			simple_rasst((PointF)enemy[n],gish.medium(),1700)  ;


		if( enemy[n].ACTUAL)
		{	n_enemy_actual.push_back(n); 
			enemy[n].setGoal(gish.medium().x, gish.medium().y);	 
		}

	}
	/////////////////////////////////////////////////////////////////
	gish.oldH = gish.H;


	for(n=0; n< CHECKPOINT_COUNT; n++) 
		if(Dist(gish.medium(), checkPoint[n]) < 70)
			respawnPoint =  checkPoint[n] ; 

 
 
	for(n=0; n< itemCount; n++)
	{
		if(Dist(gish.medium(), item[n]) < 70 )
		{
			if(gish.H < 100) 
			{	  gish.H += item[n].radius*2 ;
			item[n].setPos(100000,100000);

			}
			if(gish.H > 100) 
				gish.H = 100;
		}
	}

	/*	if(gish.FIRE)
	{
	PointF A = gish.medium();
	PointF B = ToPointF(gish.impact[0].x,gish.impact[0].y);


	int ENEMY_ACTUAL_COUNT = n_enemy_actual.size();
	float dist = 10000000000;
	int numnearEnemy = -1;
	for ( int j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
	{
	if(enemy[n_enemy_actual[j]].health>0)
	{

	float curdist = Dist( (PointF)enemy[n_enemy_actual[j]], gish.medium()) ;

	if(	 curdist < dist 
	&& CircleIntersects(&enemy[n_enemy_actual[j]].x, &enemy[n_enemy_actual[j]].y, enemy[n_enemy_actual[j]].radius, Dist(A, B), A, B))					
	{	dist = curdist;
	numnearEnemy = n_enemy_actual[j];
	}
	}
	}

	if(numnearEnemy > -1)
	{
	enemy[numnearEnemy].health -= 3;
	}
	}
	*/		

	size_t ENEMY_ACTUAL_COUNT = n_enemy_actual.size();
	for (size_t  j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
	{
		if(researchedUnit(	enemy[n_enemy_actual[j]]) )
			enemy[n_enemy_actual[j]].sleeping = false;


		if(	enemy[n_enemy_actual[j]].health>0)	
		{
			Enemy *tmpEnemy;  
			tmpEnemy = &enemy[n_enemy_actual[j]];
 

			for( int i=0;i<gish.N;i++)	
			{
				if(Dist(gish.elem[i],tmpEnemy)<tmpEnemy->radius)
				{
					gish.H -= minirand()*3;
				}

				if(gish.TVER && Dist(gish.elem[i],tmpEnemy)<tmpEnemy->radius+10)
				{
					tmpEnemy->health -=2;
				}

				/*	for(int i=0;i<10;i++)
				if(Dist(&gish.medium(), &tmpEnemy->impact[i])<tmpEnemy->radius + 20)
				{
				gish.H -= minirand()*10;
				tmpEnemy->impact[i].H=0;
				tmpEnemy->impact[i].setPos(tmpEnemy->x,tmpEnemy->y );
				//		tmpEnemy->impact[i].setDxDy(0, 0);
				}*/
			}


			for(int k=0;k<10;k++)
			{
				for( int i=0;i<gish.N;i++)
					if(Dist(gish.medium(), tmpEnemy->impact[k])<(tmpEnemy->radius + 20))
					{
						gish.H -= minirand()*10;
						tmpEnemy->impact[k].H=0;
						tmpEnemy->impact[k].setPos(tmpEnemy->x,tmpEnemy->y );
						//		tmpEnemy->impact[i].setDxDy(0, 0);
					}
 
			}
		}
	}
	
	gish.Red = false;

	if(gish.H < gish.oldH)
	{
		gish.color[0] = 1;		  
		gish.color[1] = 0;		  
		gish.color[2] = 0;		  
	}
	else
	{
		gish.color[0] = 1;		  
		gish.color[1] = 1;		  
		gish.color[2] = 1;
	}

	static float night_color = 0;
	static bool rassvet = true;
	if(night_color <= 0) rassvet = true;
	if(night_color >= 1) rassvet = false;

	if(rassvet) night_color+=0.0005;
	else night_color-=0.0005;
	// glClearColor(night_color*0.95,night_color,night_color,1);
	}//if
	/*
	if(glutGetModifiers() && GLUT_ACTIVE_SHIFT )
	{
		SLOMO_K=0.3;
	}
	else
		SLOMO_K=1.0;
	 */	
	// int dtimer = 1000/100;

	//if(gish.SLOMO)
		//dtimer = 1000/10;
	
	if (time_1 > 1000/100)  
	{
		starttime1 += time_1; 
		//		static float deltaOpen[4];
		static int napr=1;
		static int  TT = 0;

		if(TT++>300) {TT=0; napr=-napr;}

		MEDGISH = gish.medium();
		int k;

		//		for(int k=0; k<10; k++)
		//			weapon[k].setPos(0, 0);


		if(curDescI >=0 && curDescI <100 && curDescJ >=0 && curDescJ <100) 
		{
			int k=0;

			for(int i=curDescI; k<4; i++, k++ )
			{
				if(i>=0 && i<100)
					researched[i][curDescJ] = true;

				if( dat[i][curDescJ] == '0' || dat[i][curDescJ] == 'W') {}
				else 	  break;
			}
			k=0;

			for(int i=curDescI, k=0; k<4; i--, k++ )
			{
				if(i>=0 && i<100)
					researched[i][curDescJ] = true;

				if( dat[i][curDescJ] == '0' || dat[i][curDescJ] == 'W') {}
				else 	 	  break;
			}


			k=0;

			for(int j=curDescJ; k<4; j++, k++ )
			{
				if(j>=0 && j<100)
					researched[curDescI][j] = true;

				if( dat[curDescI][j] == '0' || dat[curDescI][j] == 'W') {}
				else 	 	  break;
			} 
			k=0;


			for(int j=curDescJ; k<4; j--, k++ )
			{
				if(j>=0 && j<100)
					researched[curDescI][j] = true;

				if( dat[curDescI][j] == '0' ||  dat[curDescI][j] == 'W') {}
				else 	 	  break;
			}

			for(int i=curDescI-1; i<curDescI+2; i++)
				for(int j=curDescJ-1; j<curDescJ+2; j++)
				{	 

					if( dat[i][j] != '0' )
						break;
					if(i>=0 && i<100 && j>=0 && j<100)
						researched[i][j] = true;
				}
		}

		//	if(gish.LIP)
		for(int i=0;i<gish.N;i++)
			for( k=0;k<MAXWEAPONTYPES;k++)
				if(k!=gish.currentWeapon 
					&& weapon[gish.currentWeapon].y > weapon[k].y 
					&& Dist(&weapon[gish.currentWeapon],&weapon[k]) < weapon[k].radius*2.2 )
				{	
					weapon[gish.currentWeapon].setPos(gish.medium().x, gish.medium().y+70);
					gish.currentWeapon = k;
					i=gish.N;
					k=MAXWEAPONTYPES;
				}

				
				weapon[gish.currentWeapon].setPos(gish.medium().x, gish.medium().y);

				for(int  i=0;i<gish.N;i++)	
				{
					if(Dist(gish.elem[i], gish.elem[gish.soot(i)]) < 20  )
						gish.H -= 1;

					for(int k=0; k < MAXDOORCOUNT; k++)
						if(Dist(gish.elem[i], &key[k]) < key[k].radius*1.1 )
						{ key[k].setSpirit(true);
					//	key[k].setMassive(false);

					}
				}

				for( k=0; k < MAXDOORCOUNT; k++)
				{
					if(key[k].spirit)
						key[k].setPos(gish.elem[key[k].NUM]->x , gish.elem[key[k].NUM]->y);

				}


				Collisions();

				size_t ENEMY_ACTUAL_COUNT = n_enemy_actual.size();
				for (size_t n = 0; n < ENEMY_ACTUAL_COUNT; n++)
				{

					enemy[n_enemy_actual[n]].Controller();
					enemy[n_enemy_actual[n]].InvalidateView();

					//	iter_enemy->Move();
				}

				float dist = Dist(gish.medium(),ToPointF(crossX,crossY));
				float impactVectorX = ( crossX - gish.medium().x ) / dist;
				float impactVectorY = ( crossY - gish.medium().y ) / dist;

				weapon[gish.currentWeapon].setImpactVector(impactVectorX,impactVectorY);

				if(weapon[gish.currentWeapon].fire)
				{	
				/*	for(int i=0;i<gish.N ;i++)
					{//	gish.elem[i]->CalcVector();
					gish.elem[i]->setDxDy(0,0);
					//	gish.elem[i]->CalcDxDy();
					gish.elem[i]->Move();
					}
					*/
				/*	for(int i=0;i<gish.N ;i++)
					{
						gish.elem[i]->oldx += weapon[gish.currentWeapon].impactVectorX*0.3;
						gish.elem[i]->oldy += weapon[gish.currentWeapon].impactVectorY*0.3;

					//	gish.elem[i]->dx *= 0.1;
					//	gish.elem[i]->dy *= 0.1;

						gish.elem[i]->CalcVector();
						gish.elem[i]->CalcDxDy();
						gish.elem[i]->Move();
					}
					*/
				 
				}


				for(pBall=ball.begin(); pBall!=ball.end(); pBall++)
				{
					for (  size_t j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
					{	
						Enemy *tmpEnemy;  
						tmpEnemy = &enemy[n_enemy_actual[j]];
						if(Dist(&(*pBall), tmpEnemy) < tmpEnemy->radius + pBall->radius)
						{
							tmpEnemy->health = 0;
						}
					}
				}

				for( k=0; k<10; k++)
				{
					if(k==SUPER_WEAPON) 
					{
						if(weapon[k].fire)	
							for(int i=0;i<20;i++)
								ball[i].setMassive(-1);
						else 	 
							for(int i=0;i<20;i++)
								ball[i].setMassive(1);
					}

					if(weapon[k].fire)
					{


						if(k == GRENADE || k==BLASTER || k==ROCKET_LAUNCHER)
						{			
							size_t ENEMY_ACTUAL_COUNT = n_enemy_actual.size();
							for ( size_t j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
							{
								if(	enemy[n_enemy_actual[j]].health>0)	
								{
									Enemy *tmpEnemy;  
									tmpEnemy = &enemy[n_enemy_actual[j]];

									for( int i=0;i<10;i++)	
									{					

										if(weapon[k].impact[i].H > 0 )
										{
											if(Dist(&weapon[k].impact[i],tmpEnemy)<tmpEnemy->radius)
											{
												tmpEnemy->health -= 25;
												weapon[k].impact[i].H=0;
											}
										}

										//	for(int i=0;i<10;i++)
										//	if(Dist(ToPointF(gish.elem[i]->x,gish.elem[i]->y) ,
										//			ToPointF(tmpEnemy->impact[i].x,tmpEnemy->impact[i].y)
										//			)<gish.elem[i]->radius+35)
										//	{
										//		gish.H -= 0.1;
										//		tmpEnemy->impact[i].H=0;
										//		tmpEnemy->impact[i].setPos(100000000, 100000000);
										//		tmpEnemy->impact[i].setDxDy(0, 0);
										//	}

									}


							}
						}
					}
					else if(weapon[k].impact[0].H > 0)
					{


						for(int i=0; i<10; i++)
						{
							if(weapon[k].impact[i].x == 0)
								continue;

							PointF A = ToPointF(
								weapon[k].impact[0].x, 
								weapon[k].impact[0].y);

							PointF B = ToPointF(
								weapon[k].impact[i].x, 
								weapon[k].impact[i].y);


							size_t ENEMY_ACTUAL_COUNT = n_enemy_actual.size();
							float dist = 10000000000;
							int numnearEnemy = -1;
							for ( size_t j  = 0; j < ENEMY_ACTUAL_COUNT; j++)
							{
								//	if(enemy[n_enemy_actual[j]].health>0)
								{

									float curdist = Dist( (PointF)enemy[n_enemy_actual[j]], A) ;

									if(	 curdist < dist 
										&& CircleIntersects(
										&enemy[n_enemy_actual[j]].x, 
										&enemy[n_enemy_actual[j]].y, 
										enemy[n_enemy_actual[j]].radius, 
										Dist(A, B), A, B))					
									{	dist = curdist;
									numnearEnemy = n_enemy_actual[j];
									}
								}
							}

							if(numnearEnemy > -1)
							{
								enemy[numnearEnemy].health -= 3;

								//	enemy[numnearEnemy].CalcVector();
								//	enemy[numnearEnemy].setDxDy(0,0);
								//	enemy[numnearEnemy].CalcDxDy();
								//	enemy[numnearEnemy].Move();


								enemy[numnearEnemy].CalcVector();
								//	enemy[numnearEnemy].setDxDy(weapon[k].impactVectorX*0.3 ,weapon[k].impactVectorY*0.3);
								enemy[numnearEnemy].FMove.addImpulse(weapon[k].impactVectorX*0.7 ,weapon[k].impactVectorY*0.7  );
								enemy[numnearEnemy].CalcDxDy();
								enemy[numnearEnemy].Move();

								//	weapon[k].impact[0].H = 0;
							}
						}
					}
				}

			}
	}
}

bool levelFreeSpace(char c)
{
	return (c=='0' || c=='W');
}

void LoadMap(char* path)
{
	ball.clear();

	NUM_TR = 0;
	FILE *fp = NULL;

	colorGeo[0] = minirand();
	colorGeo[1] = minirand();
	colorGeo[2] = minirand();

	colorWalk[0] = 1-minirand()*0.5;
	colorWalk[1] = 1-minirand()*0.5;
	colorWalk[2] = 1-minirand()*0.5;	

	for(int n=0;n<100;n++)
		for(int m=0;m<100; m++)
		{
			lightColor[n][m][0] = 1-minirand()*0.5;
			lightColor[n][m][1] = 1-minirand()*0.5;
			lightColor[n][m][2] = 1-minirand()*0.5;

			researched[n][m] = false;
		}

		for(int i=0;i<10;i++)
			weapon[i].set(i);

		for(int i=0; i< enemyCount; i++)
		{
			enemy[i].setPos(4000000,4000000);
		}
	 
		
		int minLevelI = 100;
		int minLevelJ = 100;
		int maxLevelI = 0;
		int maxLevelJ = 0;

		if(path[0]!=0)
		{
			fp=fopen(path,"r");
			
			int j;	
			
			for(int i = 0; i< MAXPRIMCOUNT; i++ ) {
				prim[i].type = UN;
			}



			for( j=0; j<100; j++)
				for(int i=0; i<100; i++)
				{ 
						dat[i][j] = '0';
						sector[i][j].num1 = 0;
						sector[i][j].num2 = 0;
				}
			 
			for( j=0; j<100; j++)
			{	  
				char str[101];
				fscanf(fp, "%s", str);

				for(int i=0; i<100; i++)
				{

					dat[i][j] =  str[i];//getc(fp);
				}
			}

			int counter = 0;
			int num = 0;
			float redcol=0.75,greencol=0.75,bluecol=0.75;

			int ienem=0;

			enemyCount = 0;
			doorCount = 0;
			keyCount = 0;
			itemCount = 0;

			for( j=0; j<100; j++)
			{
				for(int i=0; i<100; i++)
				{
					if((counter++)>1000)
					{
						counter = 0;
						redcol=0.4+minirand()*0.1;
						greencol=0.4+minirand()*0.1;
						bluecol=0.4+minirand()*0.1;
					}

					int I = i;
					int J = j;

					char c = '0';
					c= dat[i][j];

					PointF A = ToPointF(I*SIZE1,		 limitY - J*SIZE2);
					PointF B = ToPointF(I*SIZE1,		  limitY - J*SIZE2+SIZE2); 	
					PointF C = ToPointF(I*SIZE1 + SIZE1, limitY - J*SIZE2+SIZE2);
					PointF D = ToPointF(I*SIZE1 + SIZE1, limitY - J*SIZE2);

					PointF E;
					E.x = (A.x + D.x) / 2; 
					E.y = A.y-SIZE2/2;

					PointF A0 = ToPointF(I*SIZE1,		 limitY - J*SIZE2);
					PointF B0 = ToPointF(I*SIZE1,		  limitY - J*SIZE2+SIZE2); 	
					PointF C0 = ToPointF(I*SIZE1 + SIZE1, limitY - J*SIZE2+SIZE2);
					PointF D0 = ToPointF(I*SIZE1 + SIZE1, limitY - J*SIZE2); 


					if(c!='0')
					{
						if(I > maxLevelI) maxLevelI = I;
						if(J > maxLevelJ) maxLevelJ = J;

						if(I < minLevelI) minLevelI = I;
						if(J < minLevelJ) minLevelJ = J;

					}

					if(c == '1' || c == 'I'  || c == 's' || c== 'F')
					{

						sector[i][j].num1 = num;
						sector[i][j].num2 = num+1;

						if(c=='1')
						{
							prim[num].set(GEO, A	,B	,C,WIDTH,true,num);
							num++;
							prim[num].set(GEO,C	,D	,A,WIDTH,true,num);	
							num++;
						}
						else  	if(c=='I')
						{
							prim[num].set(ICE, A	,B	,C,WIDTH,true,num); 
							num++;
							prim[num].set(ICE,C	,D	,A,WIDTH,true,num);
							num++;
						}
						else  	if(c=='s')
						{
							prim[num].set(SAND, A	,B	,C,WIDTH,true,num);
							num++;
							prim[num].set(SAND, C	,D	,A,WIDTH,true,num);
							num++;
						}

						if(dat[i-1][j]==c && dat[i][j-1]==c &&
							dat[i+1][j]==c && dat[i][j+1]==c &&
							i>=0 && i<100 && j>=0 && j<100	)
						{
							prim[num-1].isDepth = true;
							prim[num-2].isDepth = true;
						}

					}

					if(c == 'W')
					{
						sector[i][j].num1 = num;
						sector[i][j].num2 = num+1;

						prim[num].set(HYDRO, A0	,B0	,C0,WIDTH,true, num);
						num++;
						prim[num].set(HYDRO, C0	,D0	,A0,WIDTH,true, num);
						num++;
						 
						prim[num-2].uppest = dat[i][j-1]!='W';
					}

					if(c == 'L')
					{	sector[i][j].num1 = num;

						sector[i][j].num2 = num+1;
						prim[num].set(LAVA,A0	,B0	,C0,WIDTH,true,num); 
						num++;
						prim[num].set(LAVA,C0	,D0	,A0,WIDTH,true,num); 
						num++;
						 
						prim[num-2].uppest = (dat[i][j-1]!='L');
					 
					}

					/*	if(c == 'X')
					{
					thpick.push_back(ThornPrickle(A	,B	,C,WIDTH,true,num++);		
					thpick.push_back(ThornPrickle(C	,D	,A,WIDTH,true,num++);		
					(thpick.end()-2)->setColor(0,0,0,.5,1); 
					(thpick.end()-1)->setColor(0,0,0,.5,1);
					}*/


					if(c == '2')
					{	sector[i][j].num1 = num;

					prim[num].set(GEO,A	,B	,C,WIDTH,false,num);
					num++;

						if(i<100 && dat[i+1][j]=='W') 
						{		
							prim[num].set(HYDRO, C0	,D0	,A0,WIDTH,false, num);
							sector[i][j].num2 = num;	
							num++;
						}

					}


					if(c == '3')
					{	sector[i][j].num1 = num;
						prim[num].set(GEO, A	,B	,D,WIDTH,false,num);
						num++;

						if(i<100 && dat[i+1][j]=='W') 
						{		
							prim[num].set(HYDRO,D0	,B0	,C0,WIDTH,false,num);		
							sector[i][j].num2 = num;
							num++;
						}
					}

					if(c == '4')
					{		
						prim[num].set(GEO, D,  B , C,WIDTH,false,num);
						sector[i][j].num1 = num;
						num++;

						if(i>0 && dat[i-1][j]=='W') 
						{
							sector[i][j].num1 = num;
							prim[num].set(HYDRO,A	,B	,D,WIDTH,false,num);		 
							num++;
						} 
					}

					if(c == '5')
					{
						prim[num].set(GEO,C	,D	,A,WIDTH,false,num);
						sector[i][j].num1 = num;
						num++;

						if(i>0 && dat[i-1][j]=='W') 
						{
							prim[num].set(HYDRO,A0	,B0	,C0,WIDTH,false,num);		 
							sector[i][j].num2 = num;
							num++;
						}
					}

					if(c == '6')
					{	sector[i][j].num1 = num;

						prim[num].set(ICE,A	,B	,C,WIDTH,false,num);
						num++;

						if(i<100 && dat[i+1][j]=='W') 
						{		
							prim[num].set(HYDRO, C0	,D0	,A0,WIDTH,false, num);
							sector[i][j].num2 = num;	
							num++;
						}
					}

					if(c == '7')
					{	sector[i][j].num1 = num;
						prim[num].set(ICE, A	,B	,D,WIDTH,false,num);
						num++;
						if(i<100 && dat[i+1][j]=='W') 
						{		
							prim[num].set(HYDRO,D0	,B0	,C0,WIDTH,false,num);		
							sector[i][j].num2 = num;
							num++;
						}
					}

					if(c == '8')
					{		
						if(i>0 && dat[i-1][j]=='W') 
						{
							sector[i][j].num1 = num;
							prim[num].set(HYDRO,A	,B	,D,WIDTH,false,num);		 
							num++;
						} 

						prim[num].set(ICE, D,  B , C,WIDTH,false,num++);
						sector[i][j].num2 = num;
						num++;
					}

					if(c == '9')
					{
						prim[num].set(ICE,C	,D	,A,WIDTH,false,num);
						sector[i][j].num1 = num;
						num++;

						if(i>0 && dat[i-1][j]=='W') 
						{
							prim[num].set(HYDRO,A0	,B0	,C0,WIDTH,false,num);		 
							sector[i][j].num2 = num;
							num++;
						}

					}

					if(c == 'R')	{	door[doorCount].set(GEO,E, B, C,WIDTH,false,RED); doorCount++;}
					if(c == 'G')	{	door[doorCount] .set(GEO,E, B, C,WIDTH,false,GREEN);doorCount++;}
					if(c == 'B')	{	door[doorCount].set(GEO,E, B, C,WIDTH,false,BLUE); doorCount++;}
					if(c == 'Y')	{	door[doorCount].set(GEO,E, B, C,WIDTH,false,YELLOW);doorCount++;}

					if(c == 'r')	{	key[keyCount].set(RED,I*SIZE1+75,limitY - J*SIZE2+75);keyCount++;}
					if(c == 'g')	{	key[keyCount]. set(GREEN,I*SIZE1+75,limitY - J*SIZE2+75);keyCount++;}
					if(c == 'b')	{	key[keyCount].set(BLUE,I*SIZE1+75,limitY - J*SIZE2+75);keyCount++;}
					if(c == 'y')	{	key[keyCount].set(YELLOW,I*SIZE1+75,limitY - J*SIZE2+75);keyCount++;}

					if(c == 'S')
					{
						startLevel.x = I*SIZE1+60;
						startLevel.y = limitY - J*SIZE2+60;
						respawnPoint = startLevel;
						gish.SetPos(I*SIZE1+60,limitY - J*SIZE2+60);

						item[itemCount].setPos(I*SIZE1+60, limitY - J*SIZE2+60);
						item[itemCount].setType(HEALTH);
						itemCount++;
					}

					if(c == 't' )
					{
						int enT = TURRET;
						enemy[ienem++].set(enT,I*SIZE1+SIZE1/2, limitY - J*SIZE2+SIZE2/2);enemyCount++;

						//	int count = rand()%5;
						//	if(enT != TURRET)
						//		for(int cn=0; cn<count; cn++)	
						//			enemy[ienem++].set(enT,I*SIZE1+SIZE1/2, limitY - J*SIZE2+SIZE2/2);

					} 

					if(c == 'e' )
					{
						int enT = rand()%_MAX_ENEMY_TYPE_;
						enemy[ienem++].set(enT,I*SIZE1+SIZE1/2, limitY - J*SIZE2+SIZE2/2);enemyCount++;
					}

					if(c == 'E')
					{	exitLevel.x = I*SIZE1+SIZE1/2;
						exitLevel.y = limitY - J*SIZE2+SIZE2/2;
					}
				}
			}
			fclose(fp);
			fp = NULL;
		}
		else
		{
			int	num=0;

			for(int i=0;i<50;i++)
			{		int x1= rand () %(int)limitX;
			int x2= x1+100+rand()%1000;
			int y1 = rand () %(int)limitY;
			int y2 = y1+100+rand ()%1000;

			prim[num].set(GEO, ToPointF(MIN(x1,x2),limitY-y1),
				ToPointF(MAX(x1,x2),limitY-y1+rand()%150-75),
				ToPointF((x1+x2)/2+rand()%700-350,limitY-y2),
				120,false,num++);
			prim[num].setColor(minirand(),
				minirand(),
				minirand(),1); 
			}
		}
	  
		for(int i=0;i<MAXWEAPONTYPES-1;i++)
		{
		 
				weapon[i].setPos(
					interrRand(minLevelI,maxLevelI)*SIZE1+SIZE1/2, 
				(100-interrRand(minLevelJ,maxLevelJ))*SIZE2+SIZE2/2);
		}

		for(int  icheck=0;icheck<CHECKPOINT_COUNT;icheck++)
		{
			checkPoint[icheck].setPos(
				interrRand(minLevelI,maxLevelI+10)*SIZE1+60, 
				(100-interrRand(minLevelJ,maxLevelJ+10))*SIZE2+60);
		}


		for( itemCount=0;itemCount<MAXITEMS;itemCount++)
		{

			item[itemCount].setPos(
				interrRand(minLevelI,maxLevelI+10)*SIZE1+60, 
				(100-interrRand(minLevelJ,maxLevelJ+10))*SIZE2+60);

			item[itemCount].setType(HEALTH);

		}
		
		for(int g=0;g<20;g++)
		{
			ball.push_back(RopeElement(weapon[SUPER_WEAPON].x+g,weapon[SUPER_WEAPON].y+10,10,g));

			ball[g].setInertion(0.99);

			//  if(g<12)ball[g].setMassive(-1);
			//  else 
			ball[g].setMassive(0);

			ball[g].setFriction(false);

			ball[g].cred = 0.1;
			ball[g].cgreen = 0;
			ball[g].cblue = 0;
		}
		
		ball[19].setMassive(1);
		
		playing = true;
}

bool console_mode=false;

float t=0;


//static float a = 0;


void ShowUnits()
{

	ShowHero(&gish);

	for(int i=0;i<MAXWEAPONTYPES;i++)
	{
		if(underCameraUnit (weapon[i],  camera,12) /*&& researchedUnit(weapon[i])*/)
			ShowWeaponWithImpacts(&weapon[i]);
	}

	//		if(researchedUnit(exitLevel))
	ShowExit(&exitLevel);

	for(int n=0;n<CHECKPOINT_COUNT; n++)
	{
		if(underCameraUnit (checkPoint[n],  camera,12) /*&& researchedUnit(checkPoint[n])*/)
			ShowCheckPoint(&checkPoint[n]);
	}

	for(int n=0;n<itemCount; n++)
	{
		if(underCameraUnit (item[n],  camera,12) /*&& researchedUnit(item[n])*/)
			ShowItem(&item[n]);
	}



	size_t ENEMY_ACTUAL_COUNT = n_enemy_actual.size();
	for (size_t n = 0; n < ENEMY_ACTUAL_COUNT; n++)
	{ //  if(	enemy[n_enemy_actual[n]].health>0)

		glPushMatrix();
		if(underCameraUnit (enemy[n_enemy_actual[n]],  camera,12) /*&& researchedUnit(enemy[n_enemy_actual[n]])*/)
			ShowEnemy(&enemy[n_enemy_actual[n]]);
		glPopMatrix();

		ShowEnemyImpacts(&enemy[n_enemy_actual[n]]);
	}
	
	glPushMatrix();
	for (int i = 0; i<keyCount; i++)
	{  
		ShowKey(&key[i]);
	}
	glPopMatrix();	
}




void ShowMapAllPlatforms(int iii, int jjj, int area)
{
	glBegin(GL_TRIANGLES);
	for (int n = iii-area; n<=iii+area; n++)
		for(int m = jjj-area; m<=jjj+area; m++)
		{  
			if(n >=0 && m>=0 && n<100 && m<100)
			{
				if(prim[sector[n][m].num1].type == GEO )
					ShowPlatform(&prim[sector[n][m].num1]);
				if(prim[sector[n][m].num2].type == GEO)
					ShowPlatform(&prim[sector[n][m].num2]);	
			}
		}
 	glEnd(); 

	glEnable(GL_BLEND);
	glBegin(GL_TRIANGLES);
	for (int n = iii-area; n<=iii+area; n++)
		for(int m = jjj-area; m<=jjj+area; m++)
		{  
			if(n >=0 && m>=0 && n<100 && m<100)
			{
				if(prim[sector[n][m].num1].type == LAVA )
					ShowPlatform(&prim[sector[n][m].num1]);
				if(prim[sector[n][m].num2].type == LAVA )
					ShowPlatform(&prim[sector[n][m].num2]);	
			}
		}

	for (int n = iii-area; n<=iii+area; n++)
		for(int m = jjj-area; m<=jjj+area; m++)
		{  
			if(n >=0 && m>=0 && n<100 && m<100)
			{
				if(prim[sector[n][m].num1].type == HYDRO)
					ShowPlatform(&prim[sector[n][m].num1]);
				if(prim[sector[n][m].num2].type == HYDRO)
					ShowPlatform(&prim[sector[n][m].num2]);	
			}
		}
 	glEnd(); 

	 
	glBegin(GL_TRIANGLES);
	for (int n = iii-area; n<=iii+area; n++)
		for(int m = jjj-area; m<=jjj+area; m++)
		{  
			if(n >=0 && m>=0 && n<100 && m<100)
			{
				if(prim[sector[n][m].num1].type == ICE || prim[sector[n][m].num1].type == SAND)
					ShowPlatform(&prim[sector[n][m].num1]);
				if(prim[sector[n][m].num2].type == ICE || prim[sector[n][m].num2].type == SAND)
					ShowPlatform(&prim[sector[n][m].num2]);	
			}

		}
 	glEnd(); 
 
	
	glDisable(GL_BLEND);
}


void displayByCamera(PointF cameraPos, int interval, int area)
{
	glLoadIdentity (); 

	glPushMatrix();

	YROT =  -(mousePosX - (winW>>1)) * 0.08;
	XROT =  -(mousePosY - (winH>>1)) * 0.1;

	glTranslatef(0,0 ,-(DISTANCETOSCENE - (DISTANCETOSCENE/3)));
	if(interval==0)
	{
		glRotatef(-XROT+5,1,0,0);	
		glRotatef(-YROT,0,1,0);	
	}
	glTranslatef(0,0 ,-(DISTANCETOSCENE/3));
	
	GLfloat myLightPosition[] = {-0.66, 1, 0.5, 0};
	glLightfv(GL_LIGHT0, GL_POSITION, myLightPosition);	

	if(isMap)
		glScalef(0.3,0.3,0.3);
	
	camera = cameraPos;

	glTranslatef(-cameraPos.x ,-cameraPos.y ,0);

	/*
	for(int nl=0; nl<8;nl ++)
		glDisable(GL_LIGHT0+nl);
		*/
	int lightNumber = 1;

	/*	setLocalLight(lightNumber, 	gish.medium().x, gish.medium().y ,0*150,
	1,weapon[BOMB].health/100,weapon[BOMB].health/100,
	0.000001 );

	lightNumber++;*/

/*
	if(prim[sector[curDescI][curDescJ].num2].type == HYDRO  || prim[sector[curDescI][curDescJ].num2].type == LAVA  )
	{//	ShowWeaponWithImpacts(&weapon[LAMP]);
		setLocalLight(lightNumber, gish.medium().x,gish.medium().y,180,		1,1,1, 0.0001);
		lightNumber++;  
	}
*/

	 
	int inerv = 25;

	//	if(prim[sector[curDescI][curDescJ].num2].type == HYDRO  )inerv=10;

 /*
	 	for(int n = curDescI-inerv; n<=curDescI+inerv; n++)
			for(int m = curDescJ-inerv; m<=curDescJ+inerv; m++)
			{
				//	if(n > 0 && m > 0 && n < 100 && m < 100 ) {}
				//	else continue;

				if(lightNumber<7)
				{
					if(((n%inerv) == 0 && (m%inerv) == 0 ) )
					{
						setLocalLight(lightNumber, n*SIZE1,(100-m)*SIZE2,200,
							lightColor[n][m][0],lightColor[n][m][1],lightColor[n][m][2], 0.000001);

						glPushMatrix();
						glTranslatef(n*SIZE1,(100-m)*SIZE2,10);
						glColor3f(1,1,1
							// (lightColor[n][m][0]*0.5 + 0.5),
							// (lightColor[n][m][1]*0.5 + 0.5),
							// (lightColor[n][m][2]*0.5 + 0.5)
							);
						glutSolidSphere(10,8,8);
						glPopMatrix();

						lightNumber++;
					}
				}
				else break;
			}
*/

			if(playing)
			{
				crossX = gish.medium().x + (mousePosX - (winW >> 1));
				crossY = gish.medium().y - (mousePosY - (winH >> 1));

				glColor3f(1,1,.9);
				//	glLineWidth(5);
				glLineWidth(1);
				glDisable(GL_DEPTH_TEST);
				glPushMatrix();
					glTranslatef(crossX,crossY,10);
					glutWireCube(20);
				glPopMatrix();
				glEnable(GL_DEPTH_TEST);


				::ShowImpactLine(&gish);

				int iii,jjj;
				GetQuadCoord(&iii,&jjj,cameraPos);

				glBegin(GL_TRIANGLES);
				for (int i = 0; i<doorCount; i++)
				{  
					//	if(researchedUnit(door[i].A))
					ShowDoor(&door[i]);	
				}
				glEnd();

				ShowUnits();
				/*

//		if(interval == 0)
//		{
//			glPushMatrix();
//			glTranslatef(0,0,-prim[0].Width*2);
//		//	glScalef(1,1,-1);
//			ShowUnits();
//			glPopMatrix();
//		}

				*/

				 

				ShowMapBackground(iii, jjj, area);
				ShowMapAllPlatforms(iii, jjj, area);


				for(pBall=ball.begin(); pBall!=ball.end(); pBall++)
				{
					pBall->Show(0,0);
				}
					 
		}
	glPopMatrix();
}

 
void display()
{ 
	PointF camHero  =  gish.medium();
	camHero.y += 50; 
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	displayByCamera(camHero,0,12); 

	if(gish.H>0)
	{
		glPushMatrix();
			glColor3f(0,0.75,0);
			glTranslatef(0,-3.4,-10);
			glPushMatrix();
				glScalef(gish.H*0.01,0.05,0.1);
				glutSolidCube(1);
			glPopMatrix();

			glPushMatrix();
				glScalef(1,0.01,0.01);
				glutSolidCube(1);
			glPopMatrix();
		glPopMatrix();
	} 
	//	if(pauseMode) { scene.fontBig.glDrawText(winH,winW/2,winH/2,"PAUSE");	}
	
	glutSwapBuffers();
	 
}

void mouse(int key,int state,int x,int y)
{ 
	if(key==GLUT_LEFT_BUTTON)
	{
		if(state==GLUT_DOWN)	 
		{
			weapon[gish.currentWeapon].startFire();
		}
		if(state==GLUT_UP)		
		{	
			weapon[gish.currentWeapon].endFire();
		}
	}	

	if(key==GLUT_RIGHT_BUTTON)
	{
		if(state==GLUT_UP)		
		{
			gish.LIP = false;
		}
		if(state==GLUT_DOWN)	 
		{
			gish.LIP = true;
		}
	}

	if(key==GLUT_MIDDLE_BUTTON)
	{
		if(state==GLUT_UP)		
		{gish.TVER = false;
		}
		if(state==GLUT_DOWN)
		{gish.TVER = true;
		}
	}

}

void keyboard(unsigned char key, int x, int y)
{	
	static int n=0;
	if(key == 27)	
	{
		exit(0);
	}

	if(!console_mode)
	{
		if(key == 'a' || key == 'ф')	moveLEFT_2 = true;
		if(key == 'd' || key == 'в')	moveRIGHT_2 = true;
		if(key == 'w' || key == 'ц')	moveUP_2 = true;
		if(key == 's' || key == 'ы')	moveDOWN_2= true;

		if(key ==  'r' || key ==  'к' ) 
		{
			cameraHeroOrientation = false;
			/*	for(int i=0;i <enemyCount;i++)
			{	enemy[i].setSpirit(true);
			enemy[i].setType(FLY);
			}*/

		}

		if(key ==  'p' || key ==  'з' ) 
		{
			pauseMode = !pauseMode;
		}

		if(key >= 48 && key <=57) 
			gish.currentWeapon = key - 48;
		
		if (key ==  32)	jump = true;

		/*	if(key ==  'l' || key ==  'z')	gish.LIP = true;
		if(key ==  'k')	gish.SKOL = true;
		if(key ==  'x')	gish.TVER = true;
		if(key ==  'h')	gish.FIRE = true;
		*/

		if( key == 'm' )
		{
			isMap = !isMap;
		}
		if(key == 'y') SH1 = !SH1;
	}
	else
	{
		if(key==13)
		{command[n] = '\0';n=0; 

		if(strstr(command,"maps"))
		{
			LoadMap(command);
			levelHelper.numCurrentLevel--;
		}

		if(strstr(command,"noclip"))
		{
			noclipMode = !noclipMode;
		}

		console_mode = false;
		}
		else if (key==8)
		{
			if(n>0) n--;
			for(int i=n;i<100;i++ )
				command[i] = '\0';
		}
		else	
		{
			command[n] = key;
			n++;
		}
	}

	if(key == '`')
	{
		console_mode = !console_mode;
		n=0;
		for(size_t i=0;i<strlen(command);i++)
			command[i]=0;
	}
}

void keyboardUp(unsigned char key, int x, int y)
{
	if(key == 27) exit(0);

	if(!console_mode)
	{
		if(key == 'a' || key == 'ф')	moveLEFT_2 = false;
		if(key == 'd' || key == 'в')	moveRIGHT_2 = false;
		if(key == 'w' || key == 'ц')	moveUP_2 = false;
		if(key == 's' || key == 'ы')	moveDOWN_2= false;

		if(key ==  'l' || key ==  'z')	gish.LIP = false;
		if(key ==  'k')	gish.SKOL = false;
		if(key ==  'x')	gish.TVER = false;
		if(key ==  'h')	gish.FIRE = false;

		//	if(key ==  'e' || key ==  'у' )	SLOMO_K = 1;

		if(key ==  'r' || key ==  'к' ) cameraHeroOrientation = true; 

		if(key ==  32)	jump = false;
	}
}

void specialUp(int key,int x, int y)
{
	if(key == GLUT_KEY_LEFT)	moveLEFT = false;
	if(key == GLUT_KEY_RIGHT)	moveRIGHT = false;
	if(key == GLUT_KEY_UP)		moveUP = false;
	if(key == GLUT_KEY_DOWN)	moveDOWN= false;
}

void special(int key,  int x, int y)
{
	if(key == GLUT_KEY_LEFT)	moveLEFT = true;
	if(key == GLUT_KEY_RIGHT)	moveRIGHT = true;
	if(key == GLUT_KEY_UP)		moveUP = true;
	if(key == GLUT_KEY_DOWN)	moveDOWN= true;
}

void init()
{
	srand ((unsigned)time(NULL));

#ifdef ANDROID
	LoadMap((char*)levelHelper.nextLevel("/sdcard/mapscript.dat").c_str());
#else
	LoadMap((char*)levelHelper.nextLevel("mapscript.dat").c_str());
#endif

	glBlendFunc	(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);	
	glEnable(GL_DEPTH_TEST);   
	glEnable(GL_CULL_FACE);

	starttime1=glutGet(GLUT_ELAPSED_TIME);
	starttime2=glutGet(GLUT_ELAPSED_TIME);
	starttime3=glutGet(GLUT_ELAPSED_TIME);
	tGT=glutGet(GLUT_ELAPSED_TIME);

	setDay();
}

void reshape(int width, int height)
{
	glViewport(0, 0, width, height);

	winW = width;
	winH = height;

	glMatrixMode(GL_PROJECTION);	
	glLoadIdentity();
	gluPerspective(viewAngle,    (float)width / (float)height, 5,  ZFAR);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void passiveMotion (int x, int y) {
	::mousePosX = x;
	::mousePosY = y;
}

void motion (int x, int y) {
	::mousePosX = x;
	::mousePosY = y;
}

int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

	glutInitWindowSize(winW, winH);	
	//glutInitWindowPosition(screenW/2-winW/2, screenH/2-winH/2);
	glutCreateWindow("");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutKeyboardUpFunc(keyboardUp);
	glutSpecialFunc(special);
	glutSpecialUpFunc(specialUp);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);	
	glutPassiveMotionFunc(passiveMotion);	
	glutIdleFunc(idle);

	glutSetCursor(GLUT_CURSOR_NONE);
	//	glutFullScreen();

	init();

	glutMainLoop();	

	return 0;
}

