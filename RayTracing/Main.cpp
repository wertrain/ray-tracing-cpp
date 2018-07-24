#include <cinttypes>
#include "PPM.h"

struct V
{
    double x, y, z;
    V(const double v) : x(v), y(v), z(v) {}
    V(const double x, const double y, const double z) : x(x), y(y), z(z) {}
    double &operator[](int i) const { (&x)[i]; }
};

int main()
{
    PPM ppm(400, 320);
    ppm.Fill(PPM::RGB(255,0,255));
    ppm.SaveAndPreview("result.ppm");
    return 0;
}