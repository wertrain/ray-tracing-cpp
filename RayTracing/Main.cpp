#include <cinttypes>
#include "PPM.h"

int main()
{
    PPM ppm(400, 320);
    ppm.Fill(0xFF00FF);
    ppm.SaveAndPreview("result.ppm");
    return 0;
}