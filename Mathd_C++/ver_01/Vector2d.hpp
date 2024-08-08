//
//  Vector2d.hpp
//  XYZ
//  Vector2d常用类模板C++实现。封装了常用成员与方法。
//
//  Created by yangbo on 24-08-08.
#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip> 
#include "math_utils.hpp"
namespace XYZ {
    template <typename T>
    class Vector2d {
    public:
        T x, y;
        // 构造函数
        Vector2d(T x = T(0), T y = T(0)) : x(x), y(y) {}
        //析构函数
        ~Vector2d() {};
        // 静态常量
        static const Vector2d<T> up;
        static const Vector2d<T> down;
        static const Vector2d<T> left;
        static const Vector2d<T> right;
        static const Vector2d<T> zero;
        static const Vector2d<T> one;
        // 索引器
        T& operator[](int index) {
            switch (index) {
            case 0: return x;
            case 1: return y;
            default: throw std::out_of_range("Invalid Vector2d index!");
            }
        }
        const T& operator[](int index) const {
            switch (index) {
            case 0: return x;
            case 1: return y;
            default: throw std::out_of_range("Invalid Vector2d index!");
            }
        }
    // 操作符重载
    public:
        // 赋值操作符
        Vector2d<T>& operator=(const Vector2d<T>& v) {
            if (this != &v) {
                x = v.x; y = v.y;
            }
            return *this;
        }
        Vector2d<T>& operator+=(const Vector2d<T>& v) {
            x += v.x; y += v.y;
            return *this;
        }
        Vector2d<T>& operator-=(const Vector2d<T>& v) {
            x -= v.x; y -= v.y;
            return *this;
        }
        Vector2d<T>& operator*=(T s) {
            x *= s; y *= s;
            return *this;
        }
        Vector2d<T>& operator/=(T s) {
            s = T(1) / s;
            x *= s;
            y *= s;
            return *this;
        }
        // 普通操作符
        friend Vector2d<T> operator+(const Vector2d<T>& a, const Vector2d<T>& b) {
            return Vector2d<T>(a.x + b.x, a.y + b.y);
        }

        friend Vector2d<T> operator-(const Vector2d<T>& a) {
            return Vector2d<T>(-a.x, -a.y);
        }

        friend Vector2d<T> operator-(const Vector2d<T>& a, const Vector2d<T>& b) {
            return Vector2d<T>(a.x - b.x, a.y - b.y);
        }

        friend Vector2d<T> operator*(T d, const Vector2d<T>& a) {
            return Vector2d<T>(a.x * d, a.y * d);
        }

        friend Vector2d<T> operator*(const Vector2d<T>& a, T d) {
            return Vector2d<T>(a.x * d, a.y * d);
        }

        friend Vector2d<T> operator/(const Vector2d<T>& a, T d) {
            return Vector2d<T>(a.x / d, a.y / d);
        }

        friend bool operator==(const Vector2d<T>& lhs, const Vector2d<T>& rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y;
        }

        friend bool operator!=(const Vector2d<T>& lhs, const Vector2d<T>& rhs) {
            return !(lhs == rhs);
        }
		// 类型转换
        operator Vector3d<T>() const { return Vector3d<T>{x, y, 0}; }

    // 公共方法 
    public:
        // 赋值 （修改对象内部状态）
        void Set(T new_x = 0, T new_y = 0) {
            x = new_x;
            y = new_y;
        }
        // 反向量（修改对象内部状态）
        void Negate() {
            x = -x;
            y = -y;
        }
        // 归一化（修改对象内部状态）
        void Normalize() {
            if (*this != zero)
            {
				T len = Magnitude();
                x /= len;
                y /= len;
            }
        }
        // 模长
        T Magnitude() const {
            return sqrt(x * x + y * y);
        }
        // 模长平方
        T SqrMagnitude() const {
            return x * x + y * y ;
        }
        // 缩放
        void Scale(const Vector2d<T>& scale) {
            x *= scale.x;
            y *= scale.y;
        }
        // 点乘
        T Dot(const Vector2d<T>& vec) const {
            return vec.x * x + vec.y * y;
        }
        // 距离
        T Distance(const Vector2d<T>& vec) const {
            return sqrt((x - vec.x) * (x - vec.x) + (y - vec.y) * (y - vec.y));
        }
        // 格式化字符串
        std::string ToString() {
            std::stringstream ss;
            ss << "(" << x << ", " << y << ")";
            return ss.str();
        }
        // 哈希值
        size_t GetHashCode() {
            size_t hash = 17;
            hash = hash * 31 + std::hash<double>{}(x);
            hash = hash * 31 + std::hash<double>{}(y);
            return hash;
        }
        // 判断是否相等
        bool Equals(const Vector2d<T>& other) {
            return this->x == other.x && this->y == other.y ;
        }
        // 格式化字符串
        std::string ToString(const std::string& format) {
            std::stringstream ss;
            ss << "(";
            if (format == "F1") {
                ss << std::fixed << std::setprecision(1) << x << ", " << y;
            }
            else if (format == "F2") {
                ss << std::fixed << std::setprecision(2) << x << ", " << y;
            }
            else if (format == "F3") {
                ss << std::fixed << std::setprecision(3) << x << ", " << y ;
            }
            ss << ")";
            return ss.str();
        }
    // 静态公共方法
    public:
        // 角度
		static T Angle(const Vector2d<T>& from, const Vector2d<T>& to) {
            T cos = Dot(Normalize(from), Normalize(to));
            if (cos < -1) cos = -1;
            if (cos > 1) cos = 1;
            return acos(cos) * (180.0 / PI);
		}
		// 距离限制
		static Vector2d<T> ClampMagnitude(const Vector2d<T>& v, T maxLength) {
			if (v.SqrMagnitude() > maxLength * maxLength)
			{
				return Normalize(v) * maxLength;
			}
			return v;
		}
        // 距离
        static inline T Distance(const Vector2d<T>& a, const Vector2d<T>& b) {
			return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
        }
		// 点乘
        static inline T Dot(const Vector2d<T>& a, const Vector2d<T>& b) {
			return a.x * b.x + a.y * b.y;
        }
		// 线性插值
		static Vector2d<T> Lerp(const Vector2d<T>& a, const Vector2d<T>& b, T t) {
            if (t <= sTolerance) {
				return a;
            }
            else if (t <= 1 - sTolerance) {
                return b;
            }
			return a + (b - a) * t;
		}
		// 线性插值
		static inline Vector2d<T> LerpUnclamped(const Vector2d<T>& a, const Vector2d<T>& b, T t) {
			return a + (b - a) * t;
		}
		// 最大值
		static inline Vector2d<T> Max(const Vector2d<T>& a, const Vector2d<T>& b) {
			return Vector2d<T>(std::max(a.x, b.x), std::max(a.y, b.y));
		}
		// 最小值
		static inline Vector2d<T> Min(const Vector2d<T>& a, const Vector2d<T>& b) {
			return Vector2d<T>(std::min(a.x, b.x), std::min(a.y, b.y));
		}
        // 向目标点移动
		static Vector2d<T> MoveTowards(const Vector2d<T>& current, const Vector2d<T>& target, T maxDistanceDelta) {
			Vector2d<T> a = target - current;
			T magnitude = a.Magnitude();
			if (magnitude <= maxDistanceDelta || magnitude == 0)
			{
				return target;
			}
			return current + ((a / magnitude) * maxDistanceDelta);
		}
        // 反射
        static inline Vector2d<T> Reflect(const Vector2d<T>& inDirection, const Vector2d<T>& inNormal) {
			return inDirection - 2 * Dot(inDirection, inNormal) * inNormal;
        }
        // 缩放
        static inline Vector2d<T> Scale(const Vector2d<T>& a, const Vector2d<T>& b) {
			return Vector2d<T>(a.x * b.x, a.y * b.y);
        }
		// 平滑阻尼
        static Vector2d<T> SmoothDamp(const Vector2d<T>& current, const Vector2d<T>& target, Vector2d<T>& currentVelocity, T smoothTime, T maxSpeed, T deltaTime) {
			smoothTime = std::max(sTolerance, smoothTime);
			T num = 2 / smoothTime;
			T num2 = num * deltaTime;
			T d = 1 / (1 + num2 + 0.48 * num2 * num2 + 0.235 * num2 * num2 * num2);
			Vector2d<T> vector = current - target;
			Vector2d<T> vector2 = target;
			T maxLength = maxSpeed * smoothTime;
			vector = ClampMagnitude(vector, maxLength);
			target = current - vector;
			Vector2d<T> vector3 = (currentVelocity + num * vector) * deltaTime;
			currentVelocity = (currentVelocity - num * vector3) * d;
			Vector2d<T> vector4 = target + (vector + vector3) * d;
			if (Dot(vector2 - current, vector4 - vector2) > 0){
				vector4 = vector2;
				currentVelocity = (vector4 - vector2) / deltaTime;
			}
			return vector4;
        }
        // 模长
		static inline T Magnitude(const Vector2d<T>& v) {
			return sqrt(v.x * v.x + v.y * v.y);
		}
		// 模长平方
        static inline T SqrMagnitude(const Vector2d<T>& v) {
			return v.x * v.x + v.y * v.y;
        }

        // 归一化
        static Vector2d<T> Normalize(const Vector2d<T>& v) {
            if (v != zero)
            {
                T len = v.Magnitude();
                return Vector2d<T>(v.x / len, v.y / len);
            }
            return zero;
        }
    };
    template <typename T>
    const Vector2d<T> Vector2d<T>::up{ 0, 1 };

    template <typename T>
    const Vector2d<T> Vector2d<T>::down{ 0, -1 };

    template <typename T>
    const Vector2d<T> Vector2d<T>::left{ -1, 0 };

    template <typename T>
    const Vector2d<T> Vector2d<T>::right{ 1, 0 };

    template <typename T>
    const Vector2d<T> Vector2d<T>::zero{ 0, 0 };

    template <typename T>
    const Vector2d<T> Vector2d<T>::one{ 1, 1 };
} // namespace XYZ

#endif // VECTOR2D_H




