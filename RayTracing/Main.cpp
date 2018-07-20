#include <fstream>

int main()
{
    const float imageWidth = 1200.0f;
    const float imageHeight = 800.0f;
    
    std::ofstream ofs("result.ppm");
    ofs << "P3" << std::endl;
    ofs << static_cast<int>(imageWidth) << " " << static_cast<int>(imageHeight) << std::endl;
    ofs << "255" << std::endl;
    for (int y = 0; y < static_cast<int>(imageHeight); ++y)
    {
        for (int x = 0; x < static_cast<int>(imageWidth); ++x)
        {
            ofs << "255 0 255" << std::endl;
        }
    }
    ofs.close();

    system("PNMViewer\\PNMViewer.exe result.ppm");
    return 0;
}