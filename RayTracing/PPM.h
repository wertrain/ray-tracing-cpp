/**
 * @file PPM.h
 * @author wertrain
 */
#pragma once

/**
 * PPM ‰æ‘œ‚Ìo—Í
 */
class PPM
{
public:
    static uint32_t RGB(const uint8_t r, const uint8_t g, const uint8_t b);

public:
    PPM(const uint32_t width, const uint32_t height);
    ~PPM();

    void Fill(const uint32_t color);
    void SetPixel(const uint32_t x, const uint32_t y, const uint32_t color);
    bool Save(const char* filename);
    bool SaveAndPreview(const char* filename);
    const uint32_t GetWidth() const;
    const uint32_t GetHeight() const;
    void Dump() const;
    uint32_t &operator[](int i);

private:
    uint32_t* pixels_;
    uint32_t width_;
    uint32_t height_;
};
