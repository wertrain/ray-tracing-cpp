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
        V n;                  ///< p �ł̖@��
    };

    /// ��
    struct Sphere
    {
        V p;       ///< ���S�ʒu
        double r;  ///< ���a
        V R;       ///< ���˗� [Reflectance]
        V Le;      ///< �Ɠx [Illuminance]

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

int main()
{
    Scene scene;
    //scene.spheres.push_back({ V(0), 1 });

    scene.spheres.push_back({ V(27,16.5,47), 16.5, V(.999) });
    scene.spheres.push_back({ V(73,16.5,78), 16.5, V(.999) });

    // Cornell Box
    scene.spheres.push_back({ V(1e5 + 1, 40.8, 81.6),   1e5, V(.75, .25, .25) });  // Left 
    scene.spheres.push_back({ V(-1e5 + 99, 40.8, 81.6), 1e5, V(.25, .25, .75) });  // Right 
    scene.spheres.push_back({ V(50, 40.8, 1e5),         1e5, V(.75, .75, .75) });  // Back
    //scene.spheres.push_back({ V(50, 40.8, -1e5 + 170), 1e5, V(.75, .75, .75) }); // Front
    scene.spheres.push_back({ V(50, 1e5, 81.6),         1e5, V(.75, .75, .75) });  // Bottom
    scene.spheres.push_back({ V(50, -1e5 + 81.6, 81.6), 1e5, V(.75, .75, .75) });  // Top

    // Light
    scene.spheres.push_back({ V(50,681.6 - .27,81.6), 600, V(), V(12) });

    const int width = 300;
    const int height = 300;
    const int SamplesPerPixel = 100;

    Camera camera;
    camera.eye = V(50.0, 52.0, 295.6);
    camera.center = camera.eye + V(0.0, -0.042612, -1.0);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = DegreeToRadian(30.0);
    camera.aspect = static_cast<double>(width) / static_cast<double>(height);
    
    V wE, uE, vE;
    camera.Calculate(wE, uE, vE);

    std::vector<V> I(width * height);
    const std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); // �v���J�n����

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
            for (int depth = 0; depth < 10; ++depth)
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
                            const double t = 2 * M_PI * rng.Next();
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
            I[i] = I[i] + L / SamplesPerPixel;
        }
    }
    const int64_t elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count();
    std::cout << elapsed << " millisec." << std::endl;

    PPM ppm(width, height);
    for (int i = 0; i < width * height; ++i)
    {
        ppm[i] = PPM::RGB(tonemap(I[i].x), tonemap(I[i].y), tonemap(I[i].z));
    }
    ppm.SaveAndPreview("result.ppm");

    return 0;
}