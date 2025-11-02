#pragma once

#include <string>

#include <easyx.h>

class Camera;
class Surface
{
    friend class Camera;
private:
    std::string color_type;
    double the_max[3];
    double the_min[3];
    int RGB[3];
public:
    double p[3][3];//点索引，维度索引
    bool fix;
    bool line;
    Surface(double _p[3][3], int _RGB[3], std::string color_ = "fixed_RGB", bool _fix = false, bool _line = false);
    int color(double x, double y, double z);
};

