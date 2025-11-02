#pragma once

#include <iostream>

class Pixel
{
private:
    double loc[3];
    double the_distance = 0;
    int color;
    bool has = false;
    friend class Camera;

public:
    Pixel(double a, double b, double c, int _color = 0xFFFFFF);
    Pixel(double _loc[3], int _color = 0xFFFFFF);
    void PrintInfo();
};


