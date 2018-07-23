/**
 * @file PPM.cpp
 * @author wertrain
 */
#include <cinttypes>
#include <fstream>
#include <string>

#include "PPM.h"

uint32_t PPM::RGB(const uint8_t r, const uint8_t g, const uint8_t b)
{
    return static_cast<uint32_t>((r << 16) | (g << 8) | (b << 0));
}

PPM::PPM(const uint32_t width, const uint32_t height)
    : width_(width)
    , height_(height)
{
    pixels_ = new uint32_t[width_ * height_];
}

PPM::~PPM()
{
    if (pixels_)
    {
        delete[] pixels_;
        pixels_ = nullptr;
    }
}

void PPM::Fill(const uint32_t color)
{
    for (int y = 0; y < static_cast<int>(height_); ++y)
    {
        for (int x = 0; x < static_cast<int>(width_); ++x)
        {
            pixels_[x + (y * width_)] = color;
        }
    }
}

void PPM::SetPixel(const uint32_t x, const uint32_t y, const uint32_t color)
{
    pixels_[x + (y * width_)] = color;
}

bool PPM::Save(const char* filename)
{
    std::ofstream ofs(filename);
    ofs << "P3" << std::endl;
    ofs << static_cast<int>(width_) << " " << static_cast<int>(height_) << std::endl;
    ofs << "255" << std::endl;
    for (int y = 0; y < static_cast<int>(height_); ++y)
    {
        for (int x = 0; x < static_cast<int>(width_); ++x)
        {
            const uint32_t color = pixels_[x + (y * width_)];
            uint8_t red = (color & 0xFF0000)   >> 16;
            uint8_t green = (color & 0x00FF00) >> 8;
            uint8_t blue = (color & 0x0000FF)  >> 0;
            ofs << std::to_string(red) << " " << std::to_string(green) << " " << std::to_string(blue) << std::endl;
        }
    }
    ofs.close();
    return true;
}

bool PPM::SaveAndPreview(const char* filename)
{
    if (Save(filename))
    {
        char buffer[256] = { 0 };
        strcat_s(buffer, "PNMViewer\\PNMViewer.exe ");
        strcat_s(buffer, filename);
        return system(buffer) == 0;
    }
    return false;
}

const uint32_t PPM::GetWidth() const
{
    return width_;
}

const uint32_t PPM::GetHeight() const
{
    return height_;
}

uint32_t &PPM::operator[](int i)
{
    return pixels_[i];
}
