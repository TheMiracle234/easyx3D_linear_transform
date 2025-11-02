#include "Pixel.h"
Pixel::Pixel(double a, double b, double c, int _color)
{
    loc[0] = a;
    loc[1] = b;
    loc[2] = c;
    color = _color;
}
Pixel::Pixel(double _loc[3], int _color)
{
    loc[0] = _loc[0];
    loc[1] = _loc[1];
    loc[2] = _loc[2];
    color = _color;
}
void Pixel::PrintInfo()
{
    std::cout << "loc: ";
    for (int i = 0; i < 3; ++i)
    {
        std::cout << loc[i] << " ";
    }
    std::cout << "\ncolor: " << color << std::endl;
    std::cout << "has: " << has << std::endl;
    std::cout << "the distance: " << the_distance << std::endl;
}
