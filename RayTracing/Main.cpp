#include <cinttypes>
#include "PPM.h"
#include "Math.h"

struct Ray
{
    V o; ///< ���_
    V d; ///< ����
};

int main()
{
    PPM ppm(400, 320);
    ppm.Fill(PPM::RGB(255,0,255));
    ppm.SaveAndPreview("result.ppm");
    return 0;
}