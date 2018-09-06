/**
 * @file Main.cpp
 * @author wertrain
 */
#define _USE_MATH_DEFINES
#include <cinttypes>
#include <vector>
#include <optional>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <random>
#include <tuple>
#include <omp.h>
#include "PPM.h"
#include "Math.h"

namespace
{
    struct Sphere;

    /// ���C
    struct Ray
    {
        V o; ///< ���_
        V d; ///< ����
    };

    /// �q�b�g
    struct Hit
    {
        double t;             ///< ���C�̌��_������������_�܂ł̋���
        const Sphere* sphere; ///< �q�b�g�������ւ̃|�C���^
        V p;                  ///< ���������_
        V n;                  
        ///< p �ł̖@��
    };

    /// �ގ��̃^�C�v
    enum class SurfaceType
    {
        Diffuse, ///< �g�U��
        Mirror,  ///< ����
        Fresnel, ///< ����
    };

    /// ��
    struct Sphere
    {
        V p;              ///< ���S�ʒu
        double r;         ///< ���a
        SurfaceType type; ///< �\�ʍގ��̎��
        V R;              ///< ���˗� [Reflectance]
        V Le;             ///< �Ɠx [Illuminance]
        double ior;       ///< ���ܗ�

        /** 
         * ��������
         */
        std::optional<Hit> Intersect(const Ray& ray, const double tmin, const double tmax) const
        {
            // ����̓_ x: |x - ���̒��S�ʒu| = ���̔��a
            // ���C��̓_x: x = ���C�̌��_ + ���������_�܂ł̋��� * ���C�̕���;
            const V op = p - ray.o;
            const double b = dot(op, ray.d);
            const double det = b * b - dot(op, op) + r * r;
            if (det < 0) { return {}; }
            const double t1 = b - sqrt(det);
            if (tmin < t1 && t1 < tmax) { return Hit{ t1, this }; }
            const double t2 = b + sqrt(det);
            if (tmin < t2 && t2 < tmax) { return Hit{ t2, this }; }
            return {};
        }
    };

    /// �V�[��
    struct Scene
    {
        std::vector<Sphere> spheres;
        std::optional<Hit> Intersect(const Ray& ray, double tmin, double tmax) const
        {
            std::optional<Hit> minh;
            for (const auto& sphere : spheres)
            {
                const auto h = sphere.Intersect(ray, tmin, tmax);
                
                if (!h)
                {
                    continue;
                }

                minh = h;
                tmax = minh->t;
            }
            if (minh)
            {
                const auto* s = minh->sphere;
                minh->p = ray.o + ray.d * minh->t;
                minh->n = (minh->p - s->p) / s->r;
            }
            return minh;
        }
    };

    /// �J����
    struct Camera
    {
        V eye;         ///< �J�����ʒu
        V center;      ///< �����_
        V up;          ///< �����
        double fov;    ///< ��p
        double aspect; ///< �A�X�y�N�g��

        /** 
         * �J�����̊��x�N�g�����v�Z
         */
        void Calculate(V& wE, V& uE, V& vE)
        {
            wE = normalize(eye - center);
            uE = normalize(cross(up, wE));
            vE = cross(wE, uE);
        }
    };
    
    /// ��������
    struct Random
    {
        std::mt19937 engine;
        std::uniform_real_distribution<double> dist;
        Random() {};
        Random(int seed)
        {
            engine.seed(seed);
            dist.reset();
        }
        double Next() { return dist(engine); }
    };

    /// �g�[���}�b�v
    uint8_t tonemap(const double v)
    {
        return static_cast<uint8_t>(std::min(
            std::max(int(std::pow(v, 1 / 2.2) * 255), 0), 255
        ));
    }

    V fromColor(const uint32_t color)
    {
        double red = static_cast<double>((color & 0xFF0000) >> 16);
        double green = static_cast<double>((color & 0x00FF00) >> 8);
        double blue = static_cast<double>((color & 0x0000FF) >> 0);
        return V(red, green, blue);
    }
}

// �`��p�����[�^
struct RenderParam
{
    int width;
    int height;
    int samplesPerPixel;
    int depth;
};

// ���킹���V�[�������
void createTestScene(const RenderParam& param, Scene& scene, Camera& camera)
{
    scene.spheres.push_back({ V(27,16.5,47), 16.5, SurfaceType::Diffuse, (.999) });
    scene.spheres.push_back({ V(73,16.5,78), 16.5, SurfaceType::Diffuse, V(.999) });
    scene.spheres.push_back({ V(1e5 + 1, 40.8, 81.6),   1e5, SurfaceType::Mirror, V(.75, .75, .75) });  // Left 
    scene.spheres.push_back({ V(-1e5 + 99, 40.8, 81.6), 1e5, SurfaceType::Mirror, V(.75, .75, .75) });  // Right 
    scene.spheres.push_back({ V(50, 40.8, 1e5),         1e5, SurfaceType::Diffuse, V(.25, .25, .75) });  // Back
                                                                                                         //scene.spheres.push_back({ V(50, 40.8, -1e5 + 170), 1e5, SurfaceType::Diffuse, V(.75, .75, .75) }); // Front
    scene.spheres.push_back({ V(50, 1e5, 81.6),         1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Bottom
    scene.spheres.push_back({ V(50, -1e5 + 81.6, 81.6), 1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Top
    scene.spheres.push_back({ V(50,681.6 - .27,81.6), 600, SurfaceType::Diffuse, V(), V(12) }); // Light

                                                                                                // �J�����p�����[�^�ݒ�   
    camera.eye = V(5.0, 52.0, 280.0);
    camera.center = camera.eye + V(0.6, -0.052612, -0.3);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = DegreeToRadian(30.0);
    camera.aspect = static_cast<double>(param.width) / static_cast<double>(param.height);
}

void chapter1(std::vector<V>& image, const RenderParam& param)
{
    Scene scene;
    // Spheres
    scene.spheres.push_back({ V(27,16.5,47), 16.5, SurfaceType::Diffuse, (.999) });
    scene.spheres.push_back({ V(73,16.5,78), 16.5, SurfaceType::Diffuse, V(.999) });
    // Cornell Box
    scene.spheres.push_back({ V(1e5 + 1, 40.8, 81.6),   1e5, SurfaceType::Diffuse, V(.75, .25, .25) });  // Left 
    scene.spheres.push_back({ V(-1e5 + 99, 40.8, 81.6), 1e5, SurfaceType::Diffuse, V(.25, .25, .75) });  // Right 
    scene.spheres.push_back({ V(50, 40.8, 1e5),         1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Back
    //scene.spheres.push_back({ V(50, 40.8, -1e5 + 170), 1e5, SurfaceType::Diffuse, V(.75, .75, .75) }); // Front
    scene.spheres.push_back({ V(50, 1e5, 81.6),         1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Bottom
    scene.spheres.push_back({ V(50, -1e5 + 81.6, 81.6), 1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Top
    // Light
    scene.spheres.push_back({ V(50,681.6 - .27,81.6), 600, SurfaceType::Diffuse, V(), V(12) });

    const int width = param.width;
    const int height = param.height;
    const int SamplesPerPixel = param.samplesPerPixel;

    // �J�����p�����[�^�ݒ�
    Camera camera;
    camera.eye = V(50.0, 52.0, 295.6);
    camera.center = camera.eye + V(0.0, -0.042612, -1.0);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = DegreeToRadian(30.0);
    camera.aspect = static_cast<double>(width) / static_cast<double>(height);

    V wE, uE, vE;
    camera.Calculate(wE, uE, vE);

    // ���Ԍv���̊J�n
    const std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); // �v���J�n����

    // �}���`�X���b�h�ݒ�
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < width * height; ++i)
    {
        thread_local Random rng(42 + omp_get_thread_num()); // �X���b�h���Ƃ̗���������

        for (int j = 0; j < SamplesPerPixel; ++j)
        {
            const int x = i % width;
            const int y = height - i / width; // �㉺���]
            Ray ray;
            ray.o = camera.eye;
            ray.d = [&]()
            {
                // �v���C�}�����C�̐���
                // 1. �J�������W�n�ł̕��������߂�
                // 2. ��������[���h���W�n�ɕϊ�
                const double tf = std::tan(camera.fov * 0.5);
                // �s�N�Z�����Ń����_���Ȉʒu�����߂�
                const double rpx = 2.0 * (x + rng.Next()) / width - 1;
                const double rpy = 2.0 * (y + rng.Next()) / height - 1;
                const V w = normalize(V(camera.aspect * tf * rpx, tf * rpy, -1.0));
                return uE * w.x + vE * w.y + wE * w.z;
            }();

            // ���C�̕��������܂����̂ŁA���̔��˂��v�Z���Ă���
            V L(0);  // �Ɠx�̈ꎞ�ϐ�
            V th(1); // ���̌o�H��ʂ����P�x�̕ω����L�^���邽�߂̕ϐ��ithroughput�j
            for (int depth = 0; depth < param.depth; ++depth)
            {
                // �^����ꂽ���C�ɑ΂����������
                if (const auto h = scene.Intersect(ray, 1e-4, 1e+10)) // �������g�ɂԂ���̂�����邽�߁A1e-4 �̂悤�ȏ������l������
                {
                    // ���ˉ񐔂��Ƃ̋P�x���v�Z���đ������킹��
                    L = L + th * h->sphere->Le;
                    // ���������ʒu�����Ƃɂ��Ď��̃��C�𐶐�����
                    ray.o = h->p;
                    // ���C�̕����̓����_���Ɍv�Z����
                    ray.d = [&]() {
                        // ���[�J�����W�̒��ł̕���

                        // 1. �ڋ�Ԃ̊��x�N�g�������߂�
                        // 2. 1. �̋�ԂŎ��̕����������_���ɐ���
                        // 3. 2. �̕��������[���h���W�n�ɕϊ�
                        const auto n = dot(h->n, -ray.d) > 0 ? h->n : -h->n;
                        const auto& tuple = tangentSpace(n);
                        const auto d = [&]() {
                            const double r = sqrt(rng.Next());
                            const double t = 2.0 * M_PI * rng.Next();
                            const double x = r * cos(t);
                            const double y = r * sin(t);
                            return V(x, y,
                                std::sqrt(
                                    std::max(.0, 1 - x*x - y*y)));
                        }();
                        // �����̃��[���h���W�ϊ�
                        return std::get<0>(tuple) * d.x + std::get<1>(tuple) * d.y + n * d.z;
                    }();
                    // throughput ���X�V
                    th = th * h->sphere->R;
                    if (std::max({ th.x, th.y, th.z }) == 0)
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            image[i] = image[i] + L / SamplesPerPixel;
        }
    }
    const int64_t elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count();
    std::cout << elapsed << " millisec." << std::endl;
}

void chapter2(std::vector<V>& image, const RenderParam& param)
{
    Scene scene;
    scene.spheres.push_back({ V(27,16.5,47), 16.5, SurfaceType::Mirror, (.999) });
    scene.spheres.push_back({ V(73,16.5,78), 16.5, SurfaceType::Fresnel, V(.999) });
    // Cornell Box
    scene.spheres.push_back({ V(1e5 + 1, 40.8, 81.6),   1e5, SurfaceType::Diffuse, V(.75, .25, .25) });  // Left 
    scene.spheres.push_back({ V(-1e5 + 99, 40.8, 81.6), 1e5, SurfaceType::Diffuse, V(.25, .25, .75) });  // Right 
    scene.spheres.push_back({ V(50, 40.8, 1e5),         1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Back
    //scene.spheres.push_back({ V(50, 40.8, -1e5 + 170), 1e5, SurfaceType::Diffuse, V(.75, .75, .75) }); // Front
    scene.spheres.push_back({ V(50, 1e5, 81.6),         1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Bottom
    scene.spheres.push_back({ V(50, -1e5 + 81.6, 81.6), 1e5, SurfaceType::Diffuse, V(.75, .75, .75) });  // Top
    // Light
    scene.spheres.push_back({ V(50,681.6 - .27,81.6), 600, SurfaceType::Diffuse, V(), V(12) });

    const int width = param.width;
    const int height = param.height;
    const int SamplesPerPixel = param.samplesPerPixel;

    // �J�����p�����[�^�ݒ�
    Camera camera;
    camera.eye = V(50.0, 52.0, 295.6);
    camera.center = camera.eye + V(0.0, -0.042612, -1.0);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = DegreeToRadian(30.0);
    camera.aspect = static_cast<double>(width) / static_cast<double>(height);

    // ���킹���V�[�����쐬
    //scene.spheres.clear();
    //createTestScene(param, scene, camera);

    V wE, uE, vE;
    camera.Calculate(wE, uE, vE);

    // ���Ԍv���̊J�n
    const std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); // �v���J�n����
    
    // �}���`�X���b�h�ݒ�
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < width * height; ++i)
    {
        thread_local Random rng(42 + omp_get_thread_num()); // �X���b�h���Ƃ̗���������

        for (int j = 0; j < SamplesPerPixel; ++j)
        {
            const int x = i % width;
            const int y = height - i / width; // �㉺���]
            Ray ray;
            ray.o = camera.eye;
            ray.d = [&]()
            {
                // �v���C�}�����C�̐���
                // 1. �J�������W�n�ł̕��������߂�
                // 2. ��������[���h���W�n�ɕϊ�
                const double tf = std::tan(camera.fov * 0.5);
                // �s�N�Z�����Ń����_���Ȉʒu�����߂�
                const double rpx = 2.0 * (x + rng.Next()) / width - 1;
                const double rpy = 2.0 * (y + rng.Next()) / height - 1;
                const V w = normalize(V(camera.aspect * tf * rpx, tf * rpy, -1.0));
                return uE * w.x + vE * w.y + wE * w.z;
            }();

            // ���C�̕��������܂����̂ŁA���̔��˂��v�Z���Ă���
            V L(0);  // �Ɠx�̈ꎞ�ϐ�
            V th(1); // ���̌o�H��ʂ����P�x�̕ω����L�^���邽�߂̕ϐ��ithroughput�j
            for (int depth = 0; depth < param.depth; ++depth)
            {
                // �^����ꂽ���C�ɑ΂����������
                if (const auto h = scene.Intersect(ray, 1e-4, 1e+10)) // �������g�ɂԂ���̂�����邽�߁A1e-4 �̂悤�ȏ������l������
                {
                    // ���ˉ񐔂��Ƃ̋P�x���v�Z���đ������킹��
                    L = L + th * h->sphere->Le;
                    // ���������ʒu�����Ƃɂ��Ď��̃��C�𐶐�����
                    ray.o = h->p;
                    // ���C�̕����̓����_���Ɍv�Z����
                    ray.d = [&]() {
                        // ���[�J�����W�̒��ł̕���

                        // �ގ����Ƃ̔��˂��v�Z
                        switch (h->sphere->type)
                        {
                            // �g�U����
                            case SurfaceType::Diffuse:
                            {
                                // 1. �ڋ�Ԃ̊��x�N�g�������߂�
                                // 2. 1. �̋�ԂŎ��̕����������_���ɐ���
                                // 3. 2. �̕��������[���h���W�n�ɕϊ�
                                const auto n = dot(h->n, -ray.d) > 0 ? h->n : -h->n;
                                const auto& tuple = tangentSpace(n);
                                const auto d = [&]() {
                                    const double r = sqrt(rng.Next());
                                    const double t = 2.0 * M_PI * rng.Next();
                                    const double x = r * cos(t);
                                    const double y = r * sin(t);
                                    return V(x, y,
                                        std::sqrt(
                                            std::max(.0, 1 - x*x - y*y)));
                                }();
                                // �����̃��[���h���W�ϊ�
                                return std::get<0>(tuple) * d.x + std::get<1>(tuple) * d.y + n * d.z;
                            }
                            // ���ʔ���
                            case SurfaceType::Mirror:
                            {
                                const auto wi = -ray.d;
                                return 2 * dot(wi, h->n) * h->n - wi;
                            }
                            // ����
                            case SurfaceType::Fresnel:
                            {
                                // ���˂��郌�C�̕����ɉ����āA�@������ܗ��̔���v�Z����
                                const auto wi = -ray.d;
                                const auto into = dot(wi, h->n) > 0;
                                const auto n = into ? h->n : -h->n;
                                const auto ior = h->sphere->ior;
                                const auto eta = into ? 1 / ior : ior;
                                // �X�l���̖@�����g���ċ��܂����ۂ̃x�N�g�����v�Z
                                const auto wt = [&]()->std::optional<V>
                                {
                                    const auto t = dot(wi, n);
                                    const auto t2 = 1 - eta * eta * (1 - t * t);
                                    if (t2 < 0) { return {}; } // �S���˂��N����P�[�X
                                    return eta * (n * t - wi) - n * sqrt(t2);
                                }();
                                if (!wt) return 2 * dot(wi, h->n) * h->n - wi;

                                // Schlick �̋ߎ����g���ăt���l���̎��̋ߎ����v�Z
                                const auto Fr = [&]()
                                {
                                    const auto cos = into ? dot(wi, h->n) : dot(*wt, h->n);
                                    const auto r = (1 - ior) / (1 + ior);
                                    return r * r + (1 - r * r * pow(1 - cos, 5));
                                }();
                                return rng.Next() < Fr ? 2 * dot(wi, h->n) * h->n - wi : *wt;
                            }
                        }
                        // �����ɂ͓��B���Ȃ��͂�
                        return V();
                    }();
                    // throughput ���X�V
                    th = th * h->sphere->R;
                    if (std::max({ th.x, th.y, th.z }) == 0)
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            image[i] = image[i] + L / SamplesPerPixel;
        }
    }
    const int64_t elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count();
    std::cout << elapsed << " millisec." << std::endl;
}

int main()
{
    const int width = 800;
    const int height = 600;
    const int pixels = width * height;
    const int samplesPerPixel = 1000;
    const int depth = 10;

    RenderParam param;
    param.width = width;
    param.height = height;
    param.samplesPerPixel = samplesPerPixel;
    param.depth = depth;

    std::vector<V> I(pixels);
    chapter2(I, param);

    PPM ppm(width, height);
    for (int i = 0; i < pixels; ++i)
    {
        ppm[i] = PPM::RGB(tonemap(I[i].x), tonemap(I[i].y), tonemap(I[i].z));
    }
    ppm.SaveAndPreview("result.ppm");

    return 0;
}