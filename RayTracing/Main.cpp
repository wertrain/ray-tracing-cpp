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

    uint8_t tonemap(const double v) 
    {
        return std::min<uint8_t>(
            std::max<uint8_t>(static_cast<uint8_t>(std::pow(v, 1 / 2.2) * 255), 0), 255
        );
    }
}

int main()
{
    Scene scene;
    //scene.spheres.push_back({ V(0), 1 });

    scene.spheres.push_back({ V(27,16.5,47), 16.5, V(.99) });
    scene.spheres.push_back({ V(73,16.5,78), 16.5, V(.99) });

    // Cornell Box
    scene.spheres.push_back({ V(1e5 + 1, 40.8, 81.6),   1e5, V(.75, .25, .25) });  // Left 
    scene.spheres.push_back({ V(-1e5 + 99, 40.8, 81.6), 1e5, V(.25, .25, .75) });  // Right 
    scene.spheres.push_back({ V(50, 40.8, 1e5),         1e5, V(.75, .75, .75) });  // Back
    //scene.spheres.push_back({ V(50, 40.8, -1e5 + 170), 1e5, V(.75, .75, .75) }); // Front
    scene.spheres.push_back({ V(50, 1e5, 81.6),         1e5, V(.75, .75, .75) });  // Bottom
    scene.spheres.push_back({ V(50, -1e5 + 81.6, 81.6), 1e5, V(.75, .75, .75) });  // Top

    // Light
    scene.spheres.push_back({ V(50,681.6 - .27,81.6), 600 });

    const int width = 400;
    const int height = 400;

    Camera camera;
    camera.eye = V(50.0, 52.0, 295.6);
    camera.center = camera.eye + V(0.0, -0.042612, -1.0);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = DegreeToRadian(30.0);
    camera.aspect = static_cast<double>(width) / static_cast<double>(height);
    
    V wE, uE, vE;
    camera.Calculate(wE, uE, vE);

    PPM ppm(width, height);

    //#pragma omp parallel for schedule(dynamic, 1)
    const std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); // �v���J�n����
    for (int i = 0; i < width * height; ++i)
    {
        const int x = i % width;
        //const int y = i / width;
        const int y = height - i / width; // �㉺���]
        Ray ray;
        //ray.o = V(2.0 * static_cast<double>(x) / width - 1, 2.0 * static_cast<double>(y) / height - 1, 5.0);
        //ray.d = V(0, 0, -1);
        ray.o = camera.eye;
        ray.d = [&]()
        {
            // �v���C�}�����C�̐���
            // 1. �J�������W�n�ł̕��������߂�
            // 2. ��������[���h���W�n�ɕϊ�
            const double tf = std::tan(camera.fov * 0.5);
            const double rpx = 2.0 * x / width - 1;
            const double rpy = 2.0 * y / height - 1;
            const V w = normalize(V(camera.aspect * tf * rpx, tf * rpy, -1.0));
            return uE * w.x + vE * w.y + wE * w.z;
        }();

        if (const auto h = scene.Intersect(ray, 0, 1e+10))
        {
            // �����o�[�g����
            // ���˗������ˌ��Ɩ@���̊Ԃ� cos �ɔ�Ⴗ��
            const auto n = h->sphere->R * dot(h->n, -ray.d); // �@���Ɠ��ˌ��i�J�������C�̕���������ˁj
            ppm[i] = PPM::RGB(tonemap(std::abs(n.x)), tonemap(std::abs(n.y)), tonemap(std::abs(n.z)));
        }
        else
        {
            ppm[i] = PPM::RGB(0, 0, 0);
        }
    }
    const int64_t elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count();
    std::cout << elapsed << " millisec." << std::endl;

    ppm.SaveAndPreview("result.ppm");
    return 0;
}