#include <iostream>
#include <vector>
#include <thread>
#include <conio.h>
#include <windows.h>
#include <ctime>
#include <easyx.h>

#include "Camera.h"
#include "Surface.h"

clock_t the_time1 = 0;
clock_t the_time2 = 0;

int main()
{
    initgraph(ROW, COL);
    setbkmode(TRANSPARENT);

    std::vector<Surface> surf;
    Camera obj1;
    obj1.PrintInfo();
    std::cout << std::endl;
    int range = 10;
    int high = 1;
    //添加图像方法1（推荐）
    //y=sin(x)*sin(y)
    struct { double operator()(double x, double y) { return  sin(x) * sin(y); } } myFunc;
    double max[3] = { 3.1415*2,3.1415 * 3,1 };double min[3] = { 0,0,-1 };
    obj1.pushFunc(surf,myFunc,max,min,0.25);
    //抛物面
    struct { double operator()(double x, double y) { return  ((x+5)*(x+5) + y*y)/10; } } myFunc1;
    double max1[3] = { 0,5,2.5 };double min1[3] = { -10,-5,0 };
    obj1.pushFunc(surf, myFunc1, max1, min1,0.25);

    //马鞍面
    struct { double operator()(double x, double y) { return  (x*x-y*y) / 20+5; } } myFunc2;
    double max2[3] = { 10,10,100 };double min2[3] = { -10,-10,-100 };
    obj1.pushFunc(surf, myFunc2, max2, min2, 1);

    //球面
    //struct { double operator()(double x, double y) { return  sqrt(100-x*x-y*y); } } myFunc3;
    //double max3[3] = { 15,15,10 };double min3[3] = { -15,-15,-5 };
    //obj1.pushFunc(surf, myFunc3, max3, min3, 0.01);
    //struct { double operator()(double x, double y) { return  -sqrt(100 - x * x - y * y); } } myFunc4;
    //double max4[3] = { 15,15,5 };double min4[3] = { -15,-15,-10 };
    //obj1.pushFunc(surf, myFunc4, max4, min4, 0.01);

    //添加图像方法2（不便用函数表示）
    const double little = 0.01;
    const int range2 = 10;
    double point3[6][3][3] = {
    {{0,0,0},{0,0,0},{range2,0,0},},
    {{0,0,0},{0,0,0},{range2,0,0},},
    {{0,0,0},{0,0,0},{0,range2,0},},
    {{0,0,0},{0,0,0},{0,range2,0},},
    {{0,0,0},{0,0,0},{little,little,range2},},
    {{0,0,0},{0,0,0},{little,little,range2},},
    };

    int color[6][3] = {
        {0x00,0xAA,0x00,},
        {0x00,0xAA,0x00,},
        {0xBB,0x00,0x00,},
        {0xBB,0x00,0x00,},
        {0x00,0x00,0xAA,},
        {0x00,0x00,0xAA,},
    };
    int color2[3] = { 0xFF,0xFF,0xFF, };
    //固定
    for (int i = 0;i < 6;++i) {
        surf.emplace_back(Surface(point3[i], color2, "fixed_RGB", true, true));
    }
    for (int i = 0;i < 6;++i) {
        surf.emplace_back(Surface(point3[i], color[i], "fixed_RGB", false, true));
    }
    
    obj1.Aplication(surf);
    
    system("cls");
    std::cout << "结束，按任意键退出" << std::endl;
    return 0;
}
