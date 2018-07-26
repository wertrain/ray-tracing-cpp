#include <cinttypes>
#include <vector>
#include <optional>
#include "PPM.h"
#include "Math.h"

namespace
{
    struct Sphere;

    struct Ray
    {
        V o; ///< 原点
        V d; ///< 方向
    };

    struct Hit
    {
        double t;             ///< レイの原点から交差した点までの距離
        const Sphere* sphere; ///< ヒットした球へのポインタ
    };

    struct Sphere
    {
        V p;      ///< 中心位置
        double r; ///< 半径

        std::optional<Hit> Intersect(const Ray& ray, const double tmin, const double tmax) const
        {
            // 球上の点 x: |x - 球の中心位置| = 球の半径
            // レイ上の点x: x = レイの原点 + 交差した点までの距離 * レイの方向;
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
            return minh;
        }
    };
}

int main()
{
    Scene scene;
    scene.spheres.push_back({ V(0), 1 });

    const int width = 400;
    const int height = 320;

    PPM ppm(width, height);

    for (int i = 0; i < width * height; ++i)
    {
        const int x = i % width;
        const int y = i / width;
        Ray ray;
        ray.o = V(2.0 * (static_cast<double>(x) / width) - 1, 2.0 * (static_cast<double>(y) / height) - 1, 5);
        ray.d = (0, 0, -1);

        if (const auto h = scene.Intersect(ray, 0, 1e+10))
        {
            ppm[i] = PPM::RGB(255, 0, 255);
        }
        else
        {
            ppm[i] = PPM::RGB(0, 0, 0);
        }
    }

    ppm.SaveAndPreview("result.ppm");
    return 0;
}