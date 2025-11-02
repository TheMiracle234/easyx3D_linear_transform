#pragma once
#include <iostream>
#include <cstdio>
#include <conio.h>
#include <vector>
#include <string>
#include <map>
#include <Windows.h>
#include <functional>
#include <algorithm>
#include <memory>

#define ROW 670 //横向像素数 230 700
#define COL 500 //纵向像素数 120  500
#define PI 3.1415926
//#define THETA_A 3.1415926/3 //横向分量
#define THETA_B (3.1415926/5) //纵向分量 3.1415926/10
#define DST 0.7    //距离大圆的距离  0.48
//#define OFS 4//可接受射线与点的误差 2
#define delta_theta  (3.1415926*0.006) //控制转动步长 3.1415926*0.01
#define delta_dst  0.2 //控制移动步长 5
//线性变换参数
#define d_theta (3.1415926*0.006)
#define d_mui 0.01

//测试
static int a1 = 0;
static int a2 = 0;

void ADD(double obj1[3], double obj2[3]);
//检查按键连续按下
char CheckPress();
void onePress(char& ch,char key,bool& flag);
void preciseSleep(double seconds);

struct save_pix {
    std::vector<std::vector<int>> loc;
    double distance;
    int color;
};

class Surface;
class Pixel;
class Camera
{
private:
    std::vector<save_pix> pixels;
    double loc[3];//位置
    double vec[5][3];//方向
    double theta[5][2];//角度   0:主，1234：分
    double src_RC[2][3];//屏幕像素向量
    double alpha[3];
    double mui[3][3] = { {1,0,0},{0,1,0},{0,0,1}, };
    char positive = 1;
    bool showmuis = false;
    bool showlines = false;
    void showMuis();
public:
    void RenewTheta();
    void RenewVec();
    void RenewSrc();
    Camera();
    ~Camera();
    void spanChange(double p[3][3]);
    void Change(double p[3][3]);
    void PrintInfo();
    bool shoot(Pixel& pix, int* the_screen = nullptr);
    bool shoot(Surface& surf, int);
    void Show(int);
    void Command(char ch);
    void pushFunc(std::vector<Surface>& surf,std::function<double(double,double)> func, double the_max[3], double the_min[3], double step = 0.25);

    void Aplication(std::vector<Surface> &surf);
};

