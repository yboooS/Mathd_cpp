//
//  Vector4d.hpp
//  XYZ
//  Vector4d常用类模板C++实现。封装了常用成员与方法。
//
//  Created by yangbo on 24-08-08.
#ifndef VECTOR4D_H
#define VECTOR4D_H

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip> 
#include "math_utils.hpp"
namespace XYZ {
    template <typename T>
    class Vector4d {
    public:
        T x, y, z, w;
        // 静态常量
        static const Vector4d<T> zero;
        static const Vector4d<T> one;
		// 构造函数
		Vector4d(T x = T(0), T y = T(0), T z = T(0), T w = T(0)) : x(x), y(y), z(z), w(w) {}
        // 析构函数
        ~Vector4d() {};
		// 索引器
        T& operator[](int index) {
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: throw std::out_of_range("Invalid Vector4d index!");
            }
        }
        const T& operator[](int index) const {
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: throw std::out_of_range("Invalid Vector4d index!");
            }
        }
    // 操作符重载
    public:
        // 赋值操作符
        Vector4d<T>& operator=(const Vector4d<T>& v) {
            if (this != &v) {
				x = v.x; y = v.y; z = v.z; w = v.w;
            }
            return *this;
        }
        Vector4d<T>& operator+=(const Vector4d<T>& v) {
			x += v.x; y += v.y; z += v.z; w += v.w;
            return *this;
        }
        Vector4d<T>& operator-=(const Vector4d<T>& v) {
			x -= v.x; y -= v.y; z -= v.z; w -= v.w;
            return *this;
        }
        Vector4d<T>& operator*=(T s) {
			x *= s; y *= s; z *= s; w *= s;
            return *this;
        }
        Vector4d<T>& operator/=(T s) {
            s = T(1) / s;
            x *= s;
            y *= s;
            z *= s;
			w *= s;
            return *this;
        }
        // 普通
        friend Vector4d<T> operator+(const Vector4d<T>& a, const Vector4d<T>& b) {
            return Vector4d<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
        }

        friend Vector4d<T> operator-(const Vector4d<T>& a) {
            return Vector4d<T>(-a.x, -a.y, -a.z, -a.w);
        }

        friend Vector4d<T> operator-(const Vector4d<T>& a, const Vector4d<T>& b) {
            return Vector4d<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
        }

        friend Vector4d<T> operator*(T d, const Vector4d<T>& a) {
            return Vector4d<T>(a.x * d, a.y * d, a.z * d, a.w * d);
        }

        friend Vector4d<T> operator*(const Vector4d<T>& a, T d) {
            return Vector4d<T>(a.x * d, a.y * d, a.z * d, a.w * d);
        }

        friend Vector4d<T> operator/(const Vector4d<T>& a, T d) {
            return Vector4d<T>(a.x / d, a.y / d, a.z / d, a.w / d);
        }

        friend bool operator==(const Vector4d<T>& lhs, const Vector4d<T>& rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
        }

        friend bool operator!=(const Vector4d<T>& lhs, const Vector4d<T>& rhs) {
            return !(lhs == rhs);
        }
		// 类型转换
        operator Vector2d<T>() const { return Vector2d<T>{x, y};}
        operator Vector3d<T>() const { return Vector3d<T>{x, y, z}; }
	public:
        // 赋值 （修改对象内部状态）
        void Set(T new_x = 0, T new_y = 0, T new_z = 0) {
            x = new_x;
            y = new_y;
            z = new_z;
        }
        // 归一化（修改对象内部状态）
        void Normalize() {
            T len = sqrt(SqrMagnitude());
            if (len > sTolerance)
            {
                x /= len;
                y /= len;
                z /= len;
                w /= len;
            }
        }
        // 模长
        T Magnitude() const {
            return sqrt(x * x + y * y + z * z + w * w);
        }
        // 模长平方
        T SqrMagnitude() const {
            return x * x + y * y + z * z + w * w;
        }
        // 缩放
        void Scale(const Vector4d<T>& scale) {
            x *= scale.x;
            y *= scale.y;
            z *= scale.z;
            w *= scale.w;
        }
        // 点乘
        T Dot(const Vector4d<T>& vec) const {
            return vec.x * x + vec.y * y + vec.z * z + vec.w * w;
        }
        // 距离
        T Distance(const Vector4d<T>& vec) const {
            return sqrt((x - vec.x) * (x - vec.x) + (y - vec.y) * (y - vec.y) + (z - vec.z) * (z - vec.z) + (w - vec.w) * (w - vec.w));
        }
        // 格式化字符串
        std::string ToString() {
            std::stringstream ss;
            ss << "(" << x << ", " << y << ", " << z << ", " << w << ")";
            return ss.str();
        }
        // 哈希值
        size_t GetHashCode() {
            size_t hash = 17;
            hash = hash * 31 + std::hash<double>{}(x);
            hash = hash * 31 + std::hash<double>{}(y);
            hash = hash * 31 + std::hash<double>{}(z);
            hash = hash * 31 + std::hash<double>{}(w);
            return hash;
        }
        // 判断是否相等
        bool Equals(const Vector4d<T>& other) {
            return this->x == other.x && this->y == other.y && this->z == other.z && this->w == other.w;
        }
        // 格式化字符串
        std::string ToString(const std::string& format) {
            std::stringstream ss;
            ss << "(";
            if (format == "F1") {
                ss << std::fixed << std::setprecision(1) << x << ", " << y << ", " << z << ", " << w;
            }
            else if (format == "F2") {
                ss << std::fixed << std::setprecision(2) << x << ", " << y << ", " << z << ", " << w;
            }
            else if (format == "F3") {
                ss << std::fixed << std::setprecision(3) << x << ", " << y << ", " << z << ", " << w;
            }
            ss << ")";
            return ss.str();
        }
    // 静态公共方法
    public:
        // 归一化
        static inline Vector4d<T> Normalize(const Vector4d<T>& vec) {
            T len = vec.SqrMagnitude();
            if (len > sTolerance) {
                len = sqrt(len);
                return Vector4d<T>(vec.x / len, vec.y / len, vec.z / len, vec.w / len);
            }
            return vec;
        }
        // 模长
        static inline T Magnitude(const Vector4d<T>& vec) {
            return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + vec.w * vec.w);
        }
        // 模长平方
        static inline T SqrMagnitude(const Vector4d<T>& v) {
            return v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
        }
        // 点乘
        static inline T Dot(const Vector4d<T>& vec1, const Vector4d<T>& vec2) {
            return T(vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z + vec1.w * vec2.w);
        }
        // 距离
        static inline T Distance(const Vector4d<T>& vec1, const Vector4d<T>& vec2) {
            return sqrt((vec1.x - vec2.x) * (vec1.x - vec2.x) + (vec1.y - vec2.y) * (vec1.y - vec2.y) + (vec1.z - vec2.z) * (vec1.z - vec2.z) + (vec1.w - vec2.w) * (vec1.w - vec2.w));
        }
        // 线性插值
        static inline Vector4d<T> Lerp(const Vector4d<T>& a, const Vector4d<T>& b, const T& t) {
            if (t <= sTolerance) {
                return a;
            }
            else if (t >= 1 - sTolerance) {
                return b;
            }
            return a + (b - a) * t;
        }
        // 线性插值（无限制）
        static inline Vector4d<T> LerpUnclamped(const Vector4d<T>& a, const Vector4d<T>& b, const T& t) {
            return a + (b - a) * t;
        }
        // 最大值（X,Y,Z均取最大）
        static inline Vector4d<T> Max(const Vector4d<T>& lhs, const Vector4d<T>& rhs) {
            Vector4d<T> temp;
            temp.x = std::max(lhs.x, rhs.x);
            temp.y = std::max(lhs.y, rhs.y);
            temp.z = std::max(lhs.z, rhs.z);
            temp.w = std::max(lhs.w, rhs.w);
            return temp;
        }
        // 最小值（X,Y,Z均取最小）
        static inline Vector4d<T> Min(const Vector4d<T>& lhs, const Vector4d<T>& rhs) {
            Vector4d<T> temp;
            temp.x = std::min(lhs.x, rhs.x);
            temp.y = std::min(lhs.y, rhs.y);
            temp.z = std::min(lhs.z, rhs.z);
            temp.w = std::min(lhs.w, rhs.w);
            return temp;
        }
        // 向目标点移动
        static inline Vector4d<T> MoveTowards(const Vector4d<T>& current, const Vector4d<T>& target, const T& maxDistanceDelta) {
            Vector4d<T> a = target - current;
            T magnitude = a.Magnitude();
            if (magnitude <= maxDistanceDelta || magnitude < sTolerance)
            {
                return target;
            }
            return current + a / magnitude * maxDistanceDelta;
        }
        // 向量投影
        static inline Vector4d<T> Project(const Vector4d<T>& vector, const Vector4d<T>& onNormal) {
            if (vector == zero || onNormal == zero) {
                return zero;
            }
            return Dot(vector, onNormal) / (onNormal.Magnitude() * onNormal.Magnitude()) * onNormal;
           
        }
        // 向量放缩
        static inline Vector4d<T> Scale(const Vector4d<T>& a, const Vector4d<T>& b) {
            return Vector4d<T>(a.x * b.x, a.y * b.y, a.z * b.z);
        }
    };
    template <typename T>
    const Vector4d<T> Vector4d<T>::zero{ 0, 0, 0 ,0 };

    template <typename T>
    const Vector4d<T> Vector4d<T>::one{ 1, 1, 1 ,1 };
    
} // namespace XYZW

#endif // VECTOR4D_H




