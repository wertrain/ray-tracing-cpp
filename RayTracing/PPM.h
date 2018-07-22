/**
 * @file PPM.h
 * @author wertrain
 */
#pragma once

class PPM
{
public:
    PPM(const uint32_t width, const uint32_t height);
    ~PPM();

    void Fill(const uint32_t color);
    void SetPixel(const uint32_t x, const uint32_t y, const uint32_t color);
    bool Save(const char* filename);
    bool SaveAndPreview(const char* filename);

private:
    uint32_t* pixels_;
    uint32_t width_;
    uint32_t height_;
};