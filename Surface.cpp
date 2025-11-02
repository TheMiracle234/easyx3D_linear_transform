#include "Surface.h"

Surface::Surface(double _p[3][3], int _RGB[3], std::string color_, bool _fix, bool _line)
{
    fix = _fix;
    line = _line;
    color_type = color_;
    for (int i = 0; i < 3; ++i)
    {
        RGB[i] = _RGB[i];
        for (int j = 0; j < 3; ++j)
        {
            p[i][j] = _p[i][j];
        }
    }
}
int Surface::color(double x, double y, double z)
{
    if (color_type == "fixed_RGB")
    {
        int color_code = RGB(RGB[0], RGB[1], RGB[2]);
        return color_code;
    }
    return 0xFFFFFF;
}