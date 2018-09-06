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

    /// レイ
    struct Ray
    {
        V o; ///< 原点
        V d; ///< 方向
    };

    /// ヒット
    struct Hit
    {
        double t;             ///< レイの原点から交差した点までの距離
        const Sphere* sphere; ///< ヒットした球へのポインタ
        V p;                  ///< 交差した点
        V n;                  
        ///< p での法線
    };

    /// 材質のタイプ
    enum class SurfaceType
    {
        Diffuse, ///< 拡散面
        Mirror,  ///< 鏡面
        Fresnel, ///< 屈折
    };

    /// 球
    struct Sphere
    {
        V p;              ///< 中心位置
        double r;         ///< 半径
        SurfaceType type; ///< 表面材質の種類
        V R;              ///< 反射率 [Reflectance]
        V Le;             ///< 照度 [Illuminance]
        double ior;       ///< 屈折率

        /** 
         * 交差判定
         */
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

    /// シーン
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

    /// カメラ
    struct Camera
    {
        V eye;         ///< カメラ位置
        V center;      ///< 注視点
        V up;          ///< 上方向
        double fov;    ///< 画角
        double aspect; ///< アスペクト比

        /** 
         * カメラの基底ベクトルを計算
         */
        void Calculate(V& wE, V& uE, V& vE)
        {
            wE = normalize(eye - center);
            uE = normalize(cross(up, wE));
            vE = cross(wE, uE);
        }
    };
    
    /// 乱数生成
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

    /// トーンマップ
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

// 描画パラメータ
struct RenderParam
{
    int width;
    int height;
    int samplesPerPixel;
    int depth;
};

// 合わせ鏡シーンを作る
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

                                                                                                // カメラパラメータ設定   
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

    // カメラパラメータ設定
    Camera camera;
    camera.eye = V(50.0, 52.0, 295.6);
    camera.center = camera.eye + V(0.0, -0.042612, -1.0);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = DegreeToRadian(30.0);
    camera.aspect = static_cast<double>(width) / static_cast<double>(height);

    V wE, uE, vE;
    camera.Calculate(wE, uE, vE);

    // 時間計測の開始
    const std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); // 計測開始時間

    // マルチスレッド設定
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < width * height; ++i)
    {
        thread_local Random rng(42 + omp_get_thread_num()); // スレッドごとの乱数生成器

        for (int j = 0; j < SamplesPerPixel; ++j)
        {
            const int x = i % width;
            const int y = height - i / width; // 上下反転
            Ray ray;
            ray.o = camera.eye;
            ray.d = [&]()
            {
                // プライマリレイの生成
                // 1. カメラ座標系での方向を求める
                // 2. それをワールド座標系に変換
                const double tf = std::tan(camera.fov * 0.5);
                // ピクセル内でランダムな位置を決める
                const double rpx = 2.0 * (x + rng.Next()) / width - 1;
                const double rpy = 2.0 * (y + rng.Next()) / height - 1;
                const V w = normalize(V(camera.aspect * tf * rpx, tf * rpy, -1.0));
                return uE * w.x + vE * w.y + wE * w.z;
            }();

            // レイの方向が決まったので、光の反射を計算していく
            V L(0);  // 照度の一時変数
            V th(1); // 光の経路を通した輝度の変化を記録するための変数（throughput）
            for (int depth = 0; depth < param.depth; ++depth)
            {
                // 与えられたレイに対する交差判定
                if (const auto h = scene.Intersect(ray, 1e-4, 1e+10)) // 自分自身にぶつかるのを避けるため、1e-4 のような小さい値を入れる
                {
                    // 反射回数ごとの輝度を計算して足し合わせる
                    L = L + th * h->sphere->Le;
                    // 交差した位置をもとにして次のレイを生成する
                    ray.o = h->p;
                    // レイの方向はランダムに計算する
                    ray.d = [&]() {
                        // ローカル座標の中での方向

                        // 1. 接空間の基底ベクトルを求める
                        // 2. 1. の空間で次の方向をランダムに生成
                        // 3. 2. の方向をワールド座標系に変換
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
                        // 方向のワールド座標変換
                        return std::get<0>(tuple) * d.x + std::get<1>(tuple) * d.y + n * d.z;
                    }();
                    // throughput を更新
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

    // カメラパラメータ設定
    Camera camera;
    camera.eye = V(50.0, 52.0, 295.6);
    camera.center = camera.eye + V(0.0, -0.042612, -1.0);
    camera.up = V(0.0, 1.0, 0.0);
    camera.fov = DegreeToRadian(30.0);
    camera.aspect = static_cast<double>(width) / static_cast<double>(height);

    // 合わせ鏡シーンを作成
    //scene.spheres.clear();
    //createTestScene(param, scene, camera);

    V wE, uE, vE;
    camera.Calculate(wE, uE, vE);

    // 時間計測の開始
    const std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); // 計測開始時間
    
    // マルチスレッド設定
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < width * height; ++i)
    {
        thread_local Random rng(42 + omp_get_thread_num()); // スレッドごとの乱数生成器

        for (int j = 0; j < SamplesPerPixel; ++j)
        {
            const int x = i % width;
            const int y = height - i / width; // 上下反転
            Ray ray;
            ray.o = camera.eye;
            ray.d = [&]()
            {
                // プライマリレイの生成
                // 1. カメラ座標系での方向を求める
                // 2. それをワールド座標系に変換
                const double tf = std::tan(camera.fov * 0.5);
                // ピクセル内でランダムな位置を決める
                const double rpx = 2.0 * (x + rng.Next()) / width - 1;
                const double rpy = 2.0 * (y + rng.Next()) / height - 1;
                const V w = normalize(V(camera.aspect * tf * rpx, tf * rpy, -1.0));
                return uE * w.x + vE * w.y + wE * w.z;
            }();

            // レイの方向が決まったので、光の反射を計算していく
            V L(0);  // 照度の一時変数
            V th(1); // 光の経路を通した輝度の変化を記録するための変数（throughput）
            for (int depth = 0; depth < param.depth; ++depth)
            {
                // 与えられたレイに対する交差判定
                if (const auto h = scene.Intersect(ray, 1e-4, 1e+10)) // 自分自身にぶつかるのを避けるため、1e-4 のような小さい値を入れる
                {
                    // 反射回数ごとの輝度を計算して足し合わせる
                    L = L + th * h->sphere->Le;
                    // 交差した位置をもとにして次のレイを生成する
                    ray.o = h->p;
                    // レイの方向はランダムに計算する
                    ray.d = [&]() {
                        // ローカル座標の中での方向

                        // 材質ごとの反射を計算
                        switch (h->sphere->type)
                        {
                            // 拡散反射
                            case SurfaceType::Diffuse:
                            {
                                // 1. 接空間の基底ベクトルを求める
                                // 2. 1. の空間で次の方向をランダムに生成
                                // 3. 2. の方向をワールド座標系に変換
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
                                // 方向のワールド座標変換
                                return std::get<0>(tuple) * d.x + std::get<1>(tuple) * d.y + n * d.z;
                            }
                            // 鏡面反射
                            case SurfaceType::Mirror:
                            {
                                const auto wi = -ray.d;
                                return 2 * dot(wi, h->n) * h->n - wi;
                            }
                            // 屈折
                            case SurfaceType::Fresnel:
                            {
                                // 入射するレイの方向に応じて、法線や屈折率の比を計算する
                                const auto wi = -ray.d;
                                const auto into = dot(wi, h->n) > 0;
                                const auto n = into ? h->n : -h->n;
                                const auto ior = h->sphere->ior;
                                const auto eta = into ? 1 / ior : ior;
                                // スネルの法則を使って屈折した際のベクトルを計算
                                const auto wt = [&]()->std::optional<V>
                                {
                                    const auto t = dot(wi, n);
                                    const auto t2 = 1 - eta * eta * (1 - t * t);
                                    if (t2 < 0) { return {}; } // 全反射が起きるケース
                                    return eta * (n * t - wi) - n * sqrt(t2);
                                }();
                                if (!wt) return 2 * dot(wi, h->n) * h->n - wi;

                                // Schlick の近似を使ってフレネルの式の近似を計算
                                const auto Fr = [&]()
                                {
                                    const auto cos = into ? dot(wi, h->n) : dot(*wt, h->n);
                                    const auto r = (1 - ior) / (1 + ior);
                                    return r * r + (1 - r * r * pow(1 - cos, 5));
                                }();
                                return rng.Next() < Fr ? 2 * dot(wi, h->n) * h->n - wi : *wt;
                            }
                        }
                        // ここには到達しないはず
                        return V();
                    }();
                    // throughput を更新
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