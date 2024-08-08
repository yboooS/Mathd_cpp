//
//  Vector3d.hpp
//  XYZ
//  Vector3d常用类模板C++实现。封装了常用成员与方法。
//
//  Created by yangbo on 24-08-08.
#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip> 
#include "Matrix4x4d.hpp"
#include "math_utils.hpp"


namespace XYZ {
    template <typename T>
    class Vector3d {
    public:
        T x, y, z;
		// 静态常量
        static const Vector3d<T> up;
        static const Vector3d<T> down;
        static const Vector3d<T> left;
        static const Vector3d<T> right;
        static const Vector3d<T> forward;
        static const Vector3d<T> back;
        static const Vector3d<T> zero;
        static const Vector3d<T> one;
        // 构造函数
        Vector3d(T x = T(0), T y = T(0), T z = T(0), bool isUniform = false) : x(x), y(isUniform ? x : y), z(isUniform ? x : z) { }
        Vector3d(const Vector3d& v) : x(v.x), y(v.y), z(v.z) { }
       
        // 析构函数
        ~Vector3d() {};
		// 索引器
        T& operator[](int index) {
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Invalid Vector3d index!");
            }
        }
        const T& operator[](int index) const {
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Invalid Vector3d index!");
            }
        }
	// 操作符重载
    public:
        // 赋值操作符
        Vector3d<T>& operator=(const Vector3d<T>& v) {
            if (this != &v) {
                x = v.x; y = v.y; z = v.z;
            }
            return *this;
        }
        Vector3d<T>& operator+=(const Vector3d<T>& v) {
            x += v.x; y += v.y; z += v.z;
            return *this;
        }
        Vector3d<T>& operator-=(const Vector3d<T>& v) {
            x -= v.x; y -= v.y; z -= v.z;
            return *this;
        }
        Vector3d<T>& operator*=(T s) {
            x *= s; y *= s; z *= s;
            return *this;
        }
        Vector3d<T>& operator/=(T s) {
            s = T(1) / s;
            x *= s;
            y *= s;
            z *= s;
            return *this;
        }
        // 普通操作符
        friend Vector3d<T> operator+(const Vector3d<T>& a, const Vector3d<T>& b) {
            return Vector3d<T>(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        friend Vector3d<T> operator-(const Vector3d<T>& a) {
            return Vector3d<T>(-a.x, -a.y, -a.z);
        }

        friend Vector3d<T> operator-(const Vector3d<T>& a, const Vector3d<T>& b) {
            return Vector3d<T>(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        friend Vector3d<T> operator*(T d, const Vector3d<T>& a) {
            return Vector3d<T>(a.x * d, a.y * d, a.z * d);
        }

        friend Vector3d<T> operator*(const Vector3d<T>& a, T d) {
            return Vector3d<T>(a.x * d, a.y * d, a.z * d);
        }

        friend Vector3d<T> operator/(const Vector3d<T>& a, T d) {
            return Vector3d<T>(a.x / d, a.y / d, a.z / d);
        }

        friend bool operator==(const Vector3d<T>& lhs, const Vector3d<T>& rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
        }

        friend bool operator!=(const Vector3d<T>& lhs, const Vector3d<T>& rhs) {
            return !(lhs == rhs);
        }
		// 类型转换
        operator Vector2d<T>() const { return Vector2d<T>{x, y}; }
        operator Vector4d<T>() const { return Vector4d<T>{x, y, z, 0}; }

    // 公共方法 
	public:
		// 赋值 （修改对象内部状态）
        void Set(T new_x = 0, T new_y = 0, T new_z = 0) {
            x = new_x;
            y = new_y;
            z = new_z;
        }
        // 反向量（修改对象内部状态）
        void Negate() {
            x = -x;
            y = -y;
            z = -z;
        }
		// 归一化（修改对象内部状态）
        void Normalize() {
            T len = sqrt(SqrMagnitude());
            if (len > sTolerance)
            {
                x /= len;
                y /= len;
                z /= len;
            }
        }
		// 模长
        T Magnitude() const {
            return sqrt(x * x + y * y + z * z);
        }
		// 模长平方
        T SqrMagnitude() const {
            return x * x + y * y + z * z;
        }
        // 缩放
		void Scale(const Vector3d<T>& scale) {
			x *= scale.x;
			y *= scale.y;
			z *= scale.z;
        }
        // 点乘
        T Dot(const Vector3d<T>& vec) const {
            return vec.x * x + vec.y * y + vec.z * z;
        }
        // 叉乘
        Vector3d<T> Cross(const Vector3d<T>& vec) const {
            return Vector3d<T>(y * vec.z - z * vec.y,
                z * vec.x - x * vec.z,
                x * vec.y - y * vec.x);
        }
        // 距离
        T Distance(const Vector3d<T>& vec) const {
            return sqrt((x - vec.x) * (x - vec.x) + (y - vec.y) * (y - vec.y) + (z - vec.z) * (z - vec.z));
        }
        // 格式化字符串
        std::string ToString() {
            std::stringstream ss;
            ss << "(" << x << ", " << y << ", " << z << ")";
            return ss.str();
        }
        // 哈希值
        size_t GetHashCode() {
            size_t hash = 17;
            hash = hash * 31 + std::hash<double>{}(x);
            hash = hash * 31 + std::hash<double>{}(y);
            hash = hash * 31 + std::hash<double>{}(z);
            return hash;
        }
        // 判断是否相等
        bool Equals(const Vector3d<T>& other) {
            return *this == other;
        }
        // 格式化字符串
        std::string ToString(const std::string& format) {
            std::stringstream ss;
            ss << "(";
			if (format == "F1") {
				ss << std::fixed << std::setprecision(1) << x << ", " << y << ", " << z ;
			}
			else if (format == "F2") {
				ss << std::fixed << std::setprecision(2) << x << ", " << y << ", " << z ;
			}
			else if (format == "F3") {
				ss << std::fixed << std::setprecision(3) << x << ", " << y << ", " << z ;
			}
            ss << ")";
            return ss.str();
        }
    // 静态公共方法
    public:
		// 归一化
        static inline Vector3d<T> Normalize(const Vector3d<T>& vec) {
            T len = vec.SqrMagnitude();
            if (len > sTolerance) { 
                len = sqrt(len); 
                return Vector3d<T>(vec.x / len, vec.y / len, vec.z / len);
            }
            return vec; 
        }
		// 模长
        static inline T Magnitude(const Vector3d<T>& vec) {
            return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
        }
		// 模长平方
        static inline T SqrMagnitude(const Vector3d<T>& v) {
            return v.x * v.x + v.y * v.y + v.z * v.z;
        }
		// 点乘
        static inline T Dot(const Vector3d<T>& vec1, const Vector3d<T>& vec2) {
            return T(vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
        }
		// 叉乘
        static inline Vector3d<T> Cross(const Vector3d<T>& vec1, const Vector3d<T>& vec2) {
            return Vector3d<T>(T(vec1.y * vec2.z - vec1.z * vec2.y),
                T(vec1.z * vec2.x - vec1.x * vec2.z),
                T(vec1.x * vec2.y - vec1.y * vec2.x));
        }
		// 距离
        static inline T Distance(const Vector3d<T>& vec1, const Vector3d<T>& vec2) {
            return sqrt((vec1.x - vec2.x) * (vec1.x - vec2.x) + (vec1.y - vec2.y) * (vec1.y - vec2.y) + (vec1.z - vec2.z) * (vec1.z - vec2.z));
        }
        // 夹角大小
        static inline T Angle(const Vector3d<T>& from, const Vector3d<T>& to) {
            T cos = Dot(Normalize(from), Normalize(to));
            if (cos < -1) cos = -1;
            if (cos > 1) cos = 1;
            return acos(cos) * (180.0 / PI);
        }
        // 夹角大小（弧度）
        static inline T AngleBetween(const Vector3d<T>& from, const Vector3d<T>& to) {
            T cos = Dot(Normalize(from), Normalize(to));
            if (cos < -1) cos = -1;
            if (cos > 1) cos = 1;
            return acos(cos);
        }
        // 距离限制
        static inline Vector3d<T> ClampMagnitude(const Vector3d<T>& vector, const T& maxLength){
            if (vector.SqrMagnitude() > maxLength * maxLength)
            {
                return Normalize(vector) * maxLength;
            }
            else
            {
                return vector;
            }
        }
        // 向量投影
        static inline Vector3d<T> Project(const Vector3d<T>& vector, const Vector3d<T>& onNormal) {
            if (vector == zero || onNormal == zero) {
                return zero;
            }
            return Dot(vector, onNormal) / (onNormal.Magnitude() * onNormal.Magnitude()) * onNormal;
        }
        // 去除
        static inline Vector3d<T> Exclude(const Vector3d<T>& excludeThis, const Vector3d<T>& fromThat) {
            return fromThat - Project(fromThat, excludeThis);
        }
        // 线性插值
        static inline Vector3d<T> Lerp(const Vector3d<T>& a, const Vector3d<T>& b, const T& t) {
            if (t <= sTolerance) {
				return a;
            }
            else if (t >= 1 - sTolerance) {
				return b;
            }
			return a + (b - a) * t;
        }
        // 线性插值（无限制）
        static inline Vector3d<T> LerpUnclamped(const Vector3d<T>& a, const Vector3d<T>& b, const T& t) {
            return a + (b - a) * t;
        }
        // 最大值（X,Y,Z均取最大）
		static inline Vector3d<T> Max(const Vector3d<T>& lhs, const Vector3d<T>& rhs) {
            Vector3d<T> temp;
            temp.x = std::max(lhs.x, rhs.x);
            temp.y = std::max(lhs.y, rhs.y);
            temp.z = std::max(lhs.z, rhs.z);
            return temp;
		}
		// 最小值（X,Y,Z均取最小）
		static inline Vector3d<T> Min(const Vector3d<T>& lhs, const Vector3d<T>& rhs) {
			Vector3d<T> temp;
			temp.x = std::min(lhs.x, rhs.x);
			temp.y = std::min(lhs.y, rhs.y);
			temp.z = std::min(lhs.z, rhs.z);
			return temp;
		}
        // 向目标点移动
		static inline Vector3d<T> MoveTowards(const Vector3d<T>& current, const Vector3d<T>& target, const T& maxDistanceDelta) {
			Vector3d<T> a = target - current;
			T magnitude = a.Magnitude();
			if (magnitude <= maxDistanceDelta || magnitude < sTolerance)
			{
				return target;
			}
			return current + a / magnitude * maxDistanceDelta;
		}
        // 正交法线
        void static OrthoNormalize( Vector3d<T>& normal,  Vector3d<T>& tangent) {
            T mag = normal.Magnitude();
            if (mag > sTolerance)
                normal /= mag;
            else
                normal = Vector3d<T>(1, 0, 0);

            T dot0 = Dot(normal, tangent);
            tangent -= dot0 * normal;
            mag = tangent.Magnitude();
            if (mag < sTolerance)
                tangent = OrthoNormalVectorFast(normal);
            else
                tangent /= mag;
        }
        // 正交法线
        void static OrthoNormalize(Vector3d<T>& normal, Vector3d<T>& tangent, Vector3d<T>& binormal) {
            T mag = normal.Magnitude();
            if (mag > 0)
                normal /= mag;
            else
                normal = Vector3d<T>(1, 0, 0);

            T dot0 = Dot(normal, tangent);
            tangent -= dot0 * normal;
            mag = tangent.Magnitude();
            if (mag > 0)
                tangent /= mag;
            else
                tangent = OrthoNormalVectorFast(normal);

            T dot1 = Dot(tangent, binormal);
            dot0 = Dot(normal, binormal);
            binormal -= dot0 * normal + dot1 * tangent;
            mag = binormal.Magnitude();
            if (mag > 0)
                binormal /= mag;
            else
                binormal = Cross(normal, tangent);
        }
        // 向量在平面上的投影
        static inline Vector3d<T> ProjectOnPlane(const Vector3d<T>& vector, const Vector3d<T>& planeNormal) {
            return vector - Project(vector, planeNormal);
        }
		// 反射
		static inline Vector3d<T> Reflect(const Vector3d<T>& inDirection, const Vector3d<T>& inNormal) {
			return -2 * Dot(inNormal, inDirection) * inNormal + inDirection;
		}
		// 向目标向量旋转
		static Vector3d<T> RotateTowards(const Vector3d<T>& current, const Vector3d<T>& target, const T& maxRadiansDelta, const T& maxMagnitudeDelta) {
            T currentMag = Magnitude(current);
            T targetMag = Magnitude(target);

            if (currentMag > sTolerance && targetMag > sTolerance)
            {
                Vector3d<T> currentNorm = current / currentMag;
                Vector3d<T> targetNorm = target / targetMag;

                T dot = Dot(currentNorm, targetNorm);

                if (dot > 1)
                {
                    return MoveTowards(current, target, maxMagnitudeDelta);
                }
                else if (dot < -1)
                {
                    Vector3d<T> axis = OrthoNormalVectorFast(currentNorm);
                    Matrix4x4d<T> m = SetAxisAngle(axis, maxRadiansDelta);
                    Vector3d<T> rotated = m * currentNorm;
                    rotated *= ClampedMove(currentMag, targetMag, maxMagnitudeDelta);
                    return rotated;
                }
                else
                {
                    T angle = acos(dot);
                    Vector3d<T> axis = Normalize(Cross(currentNorm, targetNorm));
                    Matrix4x4d<T> m = SetAxisAngle(axis, std::min(maxRadiansDelta, angle));
                    Vector3d<T> rotated = m * currentNorm;
                    rotated *= ClampedMove(currentMag, targetMag, maxMagnitudeDelta);
                    return rotated;
                }
            }
            else
            {
                return MoveTowards(current, target, maxMagnitudeDelta);
            }
		}
		// 向量放缩
		static inline Vector3d<T> Scale(const Vector3d<T>& a, const Vector3d<T>& b) {
			return Vector3d<T>(a.x * b.x, a.y * b.y, a.z * b.z);
		}
        // 球形插值
        static Vector3d<T> Slerp(const Vector3d<T>& lhs, const Vector3d<T>& rhs, const T& t) {
            if (t < 0) return SlerpUnclamped(lhs, rhs, 0);
            if (t > 1) return SlerpUnclamped(lhs, rhs, 1);
            return SlerpUnclamped(lhs, rhs, t);
        }
		// 球形插值（无限制）
		static Vector3d<T> SlerpUnclamped(const Vector3d<T>& lhs, const Vector3d<T>& rhs, const T& t) {
            T lhsMag = Magnitude(lhs);
            T rhsMag = Magnitude(rhs);

            if (lhsMag < sTolerance || rhsMag < sTolerance)
                return Lerp(lhs, rhs, t);

            T lerpedMagnitude = rhsMag * t + lhsMag * (1 - t);

            T dot = Dot(lhs, rhs) / (lhsMag * rhsMag);

            if (dot > 1)
            {
                return Lerp(lhs, rhs, t);
            }
            else if (dot < -1)
            {
                Vector3d<T> lhsNorm = lhs / lhsMag;
                Vector3d<T> axis = OrthoNormalVectorFast(lhsNorm);
                Matrix4x4d<T> m = SetAxisAngle(axis, PI * t);
                Vector3d<T> slerped = m * lhsNorm;
                slerped *= lerpedMagnitude;
                return slerped;
            }
            else
            {
                Vector3d<T> axis = Cross(lhs, rhs);
                Vector3d<T> lhsNorm = lhs / lhsMag;
                axis = Normalize(axis);
                double angle = acos(dot) * t;

                Matrix4x4d<T> m = SetAxisAngle(axis, angle);
                Vector3d<T> slerped = m * lhsNorm;
                slerped *= lerpedMagnitude;
                return slerped;
            }
		}
		
	// 私有方法
	private:
		static Vector3d<T> OrthoNormalVectorFast(const Vector3d<T>& n) {
			Vector3d<T> res;
			if (std::abs(n.z) > sqrt(0.5))
			{
				T a = n.y * n.y + n.z * n.z;
				T k = 1.0 / sqrt(a);
				res.x = 0;
				res.y = -n.z * k;
				res.z = n.y * k;
			}
			else
			{
				T a = n.x * n.x + n.y * n.y;
				T k = 1.0 / sqrt(a);
				res.x = -n.y * k;
				res.y = n.x * k;
				res.z = 0;
			}
			return res;
		}
        static inline T ClampedMove(T lhs, T rhs, T clampedDelta) {
            double delta = rhs - lhs;
            if (delta > 0.0F)
                return lhs + std::min(delta, clampedDelta);
            else
                return lhs - std::min(-delta, clampedDelta);
        }
        static Matrix4x4d<T> SetAxisAngle(const Vector3d<T>& rotationAxis, const T& radians) {

            Matrix4x4d<T> m;

            T s, c;
            T xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

            s = sin(radians);
            c = cos(radians);

            xx = rotationAxis.x * rotationAxis.x;
            yy = rotationAxis.y * rotationAxis.y;
            zz = rotationAxis.z * rotationAxis.z;
            xy = rotationAxis.x * rotationAxis.y;
            yz = rotationAxis.y * rotationAxis.z;
            zx = rotationAxis.z * rotationAxis.x;
            xs = rotationAxis.x * s;
            ys = rotationAxis.y * s;
            zs = rotationAxis.z * s;
            one_c = 1 - c;

            m(0, 0) = (one_c * xx) + c;
            m(0, 1) = (one_c * xy) - zs;
            m(0, 2) = (one_c * zx) + ys;

            m(1, 0) = (one_c * xy) + zs;
            m(1, 1) = (one_c * yy) + c;
            m(1, 2) = (one_c * yz) - xs;

            m(2, 0) = (one_c * zx) - ys;
            m(2, 1) = (one_c * yz) + xs;
            m(2, 2) = (one_c * zz) + c;

            return m;
        }
    };
    template <typename T>
    const Vector3d<T> Vector3d<T>::up{ 0, 1, 0 };

    template <typename T>
    const Vector3d<T> Vector3d<T>::down{ 0, -1, 0 };

    template <typename T>
    const Vector3d<T> Vector3d<T>::left{ -1, 0, 0 };

    template <typename T>
    const Vector3d<T> Vector3d<T>::right{ 1, 0, 0 };

    template <typename T>
    const Vector3d<T> Vector3d<T>::forward{ 0, 0, 1 };

    template <typename T>
    const Vector3d<T> Vector3d<T>::back{ 0, 0, -1 };

    template <typename T>
    const Vector3d<T> Vector3d<T>::zero{ 0, 0, 0 };

    template <typename T>
    const Vector3d<T> Vector3d<T>::one{ 1, 1, 1 };



} // namespace XYZ

#endif // VECTOR3_H


