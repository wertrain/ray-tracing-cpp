#define _USE_MATH_DEFINES
#include <cinttypes>
#include <vector>
#include <optional>
#include <algorithm>
#include <cmath>
#include "PPM.h"
#include "Math.h"

namespace
{
    struct Sphere;

    struct Ray
    {
        V o; ///< ���_
        V d; ///< ����
    };

    struct Hit
    {
        double t;             ///< ���C�̌��_������������_�܂ł̋���
        const Sphere* sphere; ///< �q�b�g�������ւ̃|�C���^
        V p;                  ///< ���������_
        V n;                  ///< p �ł̖@��
    };

    struct Sphere
    {
        V p;      ///< ���S�ʒu
        double r; ///< ���a

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

    struct Camera
    {
        V eye;         ///< �J�����ʒu
        V center;      ///< �����_
        V up;          ///< �����
        double fov;    ///< ��p
        double aspect; ///< �A�X�y�N�g��
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
    scene.spheres.push_back({ V(0), 1 });

    const int width = 400;
    const int height = 320;

    Camera camera;
    camera.eye = V(5.0, 5.0, 5.0);
    camera.center = V(0.0, 0.0, 0.0);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = 30.0 * M_PI / 180.0;
    camera.aspect = static_cast<double>(width) / static_cast<double>(height);

    PPM ppm(width, height);

    for (int i = 0; i < width * height; ++i)
    {
        const int x = i % width;
        const int y = i / width;
        Ray ray;
        ray.o = V(2. * static_cast<double>(x) / width - 1, 2. * static_cast<double>(y) / height - 1, 5.);
        ray.d = V(0, 0, -1);

        if (const auto h = scene.Intersect(ray, 0, 1e+10))
        {
            const auto n = h->n;
            ppm[i] = PPM::RGB(tonemap(n.x), tonemap(n.y), tonemap(n.z));
        }
        else
        {
            ppm[i] = PPM::RGB(0, 0, 0);
        }
    }

    ppm.SaveAndPreview("result.ppm");
    return 0;
}