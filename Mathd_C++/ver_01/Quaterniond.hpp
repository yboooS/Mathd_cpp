//
//  Quaterniond.hpp
//  XYZ
//  Quaterniond常用类模板C++实现。封装了常用成员与方法。
//
//  Created by yangbo on 24-08-08.
#ifndef QUATERNIOND_H
#define QUATERNIOND_H

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip> 
#include "math_utils.hpp"
#include "Vector3d.hpp"
#include "Matrix4x4d.hpp"

namespace XYZ {
    template <typename T>
    class Quaterniond{
    public:
        T x, y, z, w;
        // 静态常量
        static const Quaterniond<T> identity;
        // 构造函数
        Quaterniond(T x = T(0), T y = T(0), T z = T(0), T w = T(0)) : x(x), y(y), z(z), w(w) { }
        // 析构函数
        ~Quaterniond() {};
        // 索引器
        T& operator[](int index) {
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: throw std::out_of_range("Invalid Quaterniond index!");
            }
        }
        const T& operator[](int index) const {
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: throw std::out_of_range("Invalid Quaterniond index!");
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
        Vector3d<T>& operator*=(T s) {
            x *= s; y *= s; z *= s;
            return *this;
        }
        // 普通操作符
        friend Quaterniond<T> operator*(const Quaterniond<T>& lhs, const Quaterniond<T>& rhs) {
            return Vector3d<T>(lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y,
                lhs.w * rhs.y + lhs.y * rhs.w + lhs.z * rhs.x - lhs.x * rhs.z,
                lhs.w * rhs.z + lhs.z * rhs.w + lhs.x * rhs.y - lhs.y * rhs.x,
                lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z);
        }
        friend Vector3d<T> operator*(const Quaterniond<T>& rotation, const Vector3d<T>& point) {
            T num = rotation.x * 2;
            T num2 = rotation.y * 2;
            T num3 = rotation.z * 2;
            T num4 = rotation.x * num;
            T num5 = rotation.y * num2;
            T num6 = rotation.z * num3;
            T num7 = rotation.x * num2;
            T num8 = rotation.x * num3;
            T num9 = rotation.y * num3;
            T num10 = rotation.w * num;
            T num11 = rotation.w * num2;
            T num12 = rotation.w * num3;
            Vector3d<T> result;
            result.x = (1 - (num5 + num6)) * point.x + (num7 - num12) * point.y + (num8 + num11) * point.z;
            result.y = (num7 + num12) * point.x + (1 - (num4 + num6)) * point.y + (num9 - num10) * point.z;
            result.z = (num8 - num11) * point.x + (num9 + num10) * point.y + (1 - (num4 + num5)) * point.z;
            return result;
        }
        friend bool operator==(const Quaterniond<T>& lhs, const Quaterniond<T>& rhs) {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w;
        }
        friend bool operator!=(const Quaterniond<T>& lhs, const Quaterniond<T>& rhs) {
            return !(lhs == rhs);
        }


        // 公共方法 
    public:
        // 赋值 欧拉角（修改对象内部状态）
        void SetEulerAngles(Vector3d<T>& value) {
            *this = Euler(value);
        }
        // 设置四元数（修改对象内部状态）
        void Set(T new_x = 0, T new_y = 0, T new_z = 0, T new_w = 0) {
            x = new_x;
            y = new_y;
            z = new_z;
            w = new_w;
        }
        // 旋转的表示
        void SetFromToRotation(const Vector3d<T>& fromDirection, const Vector3d<T>& toDirection) {
            // 计算旋转的逻辑
            *this = FromToRotation(fromDirection, toDirection);
        }
        // 设置注视旋转
        void SetLookRotation(const Vector3d<T>& view) {
            // 计算旋转的逻辑
            *this = LookRotation(view);
        }
        // 设置注视旋转
        void SetLookRotation(const Vector3d<T>& view, const Vector3d<T>& up = Vector3d<T>::Up()) {
            // 计算旋转的逻辑
            *this = LookRotation(view, up);
        }
        // 转换为角轴
        void ToAngleAxis(T& angle, Vector3d<T>& axis) const {
            angle = 2.0 * std::acos(w);
            if (angle == 0) {
                axis = Vector3d<T>::Right();
                return;
            }

            T div = 1.0 / std::sqrt(1 - w * w);
            axis.Set(x * div, y * div, z * div);
            angle = angle * 180.0 / PI;  // 
        }
        // 点乘
        T Dot(const Quaterniond<T>& vec) const {
            return vec.x * x + vec.y * y + vec.z * z + vec.w * w;
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
        bool Equals(const Quaterniond<T>& other) {
            return *this == other;
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
        // 获取欧拉角
        static inline Vector3d<T> GetEulerAngles() {
            Matrix4x4d<T> m = QuaternionToMatrix(*this);
            return (MatrixToEuler(m) * 180 / PI);
        }
        static Matrix4x4d<T> QuaternionToMatrix(const Quaterniond<T>& quat) {
            Matrix4x4d<T> m;

            T x = quat.x * 2;
            T y = quat.y * 2;
            T z = quat.z * 2;
            T xx = quat.x * x;
            T yy = quat.y * y;
            T zz = quat.z * z;
            T xy = quat.x * y;
            T xz = quat.x * z;
            T yz = quat.y * z;
            T wx = quat.w * x;
            T wy = quat.w * y;
            T wz = quat.w * z;

            m[0] = 1.0 - (yy + zz);
            m[1] = xy + wz;
            m[2] = xz - wy;
            m[3] = 0.0;

            m[4] = xy - wz;
            m[5] = 1.0 - (xx + zz);
            m[6] = yz + wx;
            m[7] = 0.0;

            m[8] = xz + wy;
            m[9] = yz - wx;
            m[10] = 1.0 - (xx + yy);
            m[11] = 0.0;

            m[12] = 0.0;
            m[13] = 0.0;
            m[14] = 0.0;
            m[15] = 1.0;

            return m;
        }
        // 夹角大小
		static inline T Angle(const Quaterniond<T>& a, const Quaterniond<T>& b) {
			T f = Quaterniond<T>::Dot(a, b);
			return acos(std::min(std::abs(f), 1.0)) * 2 * (180 / PI);
		}
		// 轴向旋转
		static Quaterniond<T> AxisAngle(const T& angle, const Vector3d<T>& axis) {
			axis.Normalize();
			angle = angle / 180 * PI;

			Quaterniond<T> q;

			T halfAngle = angle * 0.5;
			T s = sin(halfAngle);
			q.x = axis.x * s;
			q.y = axis.y * s;
			q.z = axis.z * s;
			q.w = cos(halfAngle);

			return q;
		}
        // 点乘
		static inline T Dot(const Quaterniond<T>& a, const Quaterniond<T>& b)
		{
			return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
		}

        // 欧拉角转四元数
        static inline Quaterniond<T> Euler(const Vector3d<T>& euler)
        {
            return Euler(euler.x, euler.y, euler.z);
        }
		// 欧拉角转四元数
        static Quaterniond<T> Euler(const T& x, const T& y, const T& z)
        {
            T cX = cos(x * PI / 360);
            T sX = sin(x * PI / 360);

            T cY = cos(y * PI / 360);
            T sY = sin(y * PI / 360);

            T cZ = cos(z * PI / 360);
            T sZ = sin(z * PI / 360);

            Quaterniond<T> qX(sX, 0, 0, cX);
            Quaterniond<T> qY(0, sY, 0, cY);
            Quaterniond<T> qZ(0, 0, sZ, cZ);

            Quaterniond<T> q = (qY * qX) * qZ;

            return q;
            
        }
        // 向量间的角度
        static Quaterniond<T> FromToRotation(const Vector3d<T>& fromDirection, const Vector3d<T>& toDirection){
			fromDirection.Normalize();
			toDirection.Normalize();
			T dot = Vector3d<T>::Dot(fromDirection, toDirection);
			if (dot > 1 - sTolerance) {
				return Quaterniond<T>::identity;
			}
			if (dot < sTolerance - 1) {
				Vector3d<T> axis = Vector3d<T>::Cross(Vector3d<T>::right, fromDirection);
				if (axis.Magnitude() < sTolerance) {
					axis = Vector3d<T>::Cross(Vector3d<T>::up, fromDirection);
				}
				axis.Normalize();
				return Quaterniond<T>(axis, 180);
			}
			T angle = acos(dot) * 180 / PI;
			Vector3d<T> axis = Vector3d<T>::Cross(fromDirection, toDirection);
			axis.Normalize();
			return Quaterniond<T>(axis, angle);
        }
        // 四元数的逆
		static inline Quaterniond<T> Inverse(const Quaterniond<T>& rotation) {
			return Quaterniond<T>(-rotation.x, -rotation.y, -rotation.z, rotation.w);
		}
		// 线性插值
        static inline Quaterniond<T> Lerp(const Quaterniond<T>& a, const Quaterniond<T>& b, T t) {
            if (t > 1) {
				t = 1;
            }
            if (t < 0) {
				t = 0;
            }
			return LerpUnclamped(a, b, t);
        }
		// 线性插值（无限制）
		static Quaterniond<T> LerpUnclamped(const Quaterniond<T>& a, const Quaterniond<T>& b, T t) {
            Quaterniond<T> temQuat;
            if (Dot(a, b)<0) {
                temQuat.Set(
                    a.x + t * (-b.x - a.x),
					a.y + t * (-b.y - a.y), 
                    a.z + t * (-b.z - a.z),
					a.w + t * (-b.w - a.w));
            }
            else
			{
				temQuat.Set(
					a.x + t * (b.x - a.x),
					a.y + t * (b.y - a.y),
					a.z + t * (b.z - a.z),
                    a.w + t * (b.w - a.w));
			}
            T nor = sqrt(Dot(temQuat, temQuat));
			return Quaterniond<T>(temQuat.x / nor, temQuat.y / nor, temQuat.z / nor, temQuat.w / nor);
		}
        // 注视旋转
        static inline  Quaterniond<T> LookRotation(const Vector3d<T>& forward) {
			Vector3d<T>& up = Vector3d<T>::up;
			return LookRotation(forward, up);
        }
		// 注视旋转
        static inline  Quaterniond<T> LookRotation(const Vector3d<T>& forward, const Vector3d<T>& upwards = Vector3d<T>::Up()) {
			Matrix4x4d<T> m = LookRotationToMatrix(forward, upwards);
            return MatrixToQuaternion(m);
        }
        // 向目标角度旋转
		static inline Quaterniond<T> RotateTowards(const Quaterniond<T>& from, const Quaterniond<T>& to, T maxDegreesDelta) {
			T num = Quaterniond<T>::Angle(from, to);
			if (num == 0) {
				return to;
			}
		
			else {
				T t = std::min(1, maxDegreesDelta / num);
				return SlerpUnclamped(from, to, t);
				
			}
		}
        // 球形插值
        static inline Quaterniond<T> Slerp(const Quaterniond<T>& a, const Quaterniond<T>& b, T t) {
            if (t > 1) {
                t = 1;
			}
            if (t < 0) {
				t = 0;
            }
			return SlerpUnclamped(a, b, t);
        }
		// 球形插值（无限制）
        static Quaterniond<T> SlerpUnclamped(const Quaterniond<T>& q1, const Quaterniond<T>& q2, T t) {
			T dot = Quaterniond<T>::Dot(q1, q2);
			Quaterniond<T> tempQuat;
            if (dot < 0) {
				dot = -dot;
				tempQuat.Set(-q1.x, -q1.y, -q1.z, -q1.w);
			}
			else {
                tempQuat = q2;
            }

            if (dot < 1) {
				T angle = acos(dot);
				T sinadiv, sinat, sinaomt;
				sinadiv = 1 / sin(angle);
                sinat = sin(angle * t);
				sinaomt = sin(angle * (1 - t));
                tempQuat.Set(
					(q1.x * sinaomt + tempQuat.x * sinat) * sinadiv, 
					(q1.y * sinaomt + tempQuat.y * sinat) * sinadiv, 
					(q1.z * sinaomt + tempQuat.z * sinat) * sinadiv, 
					(q1.w * sinaomt + tempQuat.w * sinat) * sinadiv);
				return tempQuat;
            }
            else {
				return Lerp(q1, tempQuat, t);
            }
        }



        // 私有方法
    private:
        Vector3d<T> MatrixToEuler(const Matrix4x4d<T>& m)
        {
            Vector3d<T>  v;
            if (m(1, 2) < 1)
            {
                if (m(1, 2) > -1)
                {
                    v.x = asin(-m(1, 2));
                    v.y = atan2(m(0, 2), m(2, 2));
                    v.z = atan2(m(1, 0), m(1, 1));
                }
                else
                {
                    v.x = PI * 0.5;
                    v.y = atan2(m(0, 1), m(0, 0));
                    v.z = 0;
                }
            }
            else
            {
                v.x = -PI * 0.5;
                v.y = atan2(-m(0, 1), m(0, 0));
                v.z = 0;
            }

            for (int i = 0; i < 3; i++)
            {
                if (v[i] < 0)
                {
                    v[i] += 2 * PI;
                }
                else if (v[i] > 2 * PI)
                {
                    v[i] -= 2 * PI;
                }
            }

            return v;
        }
        static Quaterniond<T> MatrixToQuaternion(const Matrix4x4d<T>& m) {
			Quaterniond<T> quat;
            T fTrace = m(0, 0) + m(1, 1) + m(2, 2);
            T root;

            if (fTrace > 0)
            {
                root = sqrt(fTrace + 1);
                quat.w = 0.5 * root;
                root = 0.5 / root;
                quat.x = (m(2, 1) - m(1, 2)) * root;
                quat.y = (m(0, 2) - m(2, 0)) * root;
                quat.z = (m(1, 0) - m(0, 1)) * root;
            }
            else
            {
                std::array<int, 3> s_iNext = { 1, 2, 0 };
                int i = 0;
                if (m(1, 1) > m(0, 0))
                {
                    i = 1;
                }
                if (m(2, 2) > m(i, i))
                {
                    i = 2;
                }
                int j = s_iNext[i];
                int k = s_iNext[j];

                root = sqrt(m(i, i) - m(j, j) - m(k, k) + 1);
                if (root < 0)
                {
                    // 提示越界异常
					throw std::out_of_range("Quaterniond<T> MatrixToQuaternion(const Matrix4x4d<T>& m) root < 0");
                }
                quat[i] = 0.5 * root;
                root = 0.5 / root;
                quat.w = (m(k, j) - m(j, k)) * root;
                quat[j] = (m(j, i) + m(i, j)) * root;
                quat[k] = (m(k, i) + m(i, k)) * root;
            }
            T nor = sqrt(Dot(quat, quat));
            quat = Quaterniond<T>(quat.x / nor, quat.y / nor, quat.z / nor, quat.w / nor);

            return quat;
        }
        static Matrix4x4d<T> LookRotationToMatrix(const Vector3d<T>& viewVec, const Vector3d<T>& upVec) {
			Vector3d<T> z = viewVec;
			Matrix4x4d<T> m;

			T mag = Vector3d<T>::Magnitude(z);
            if (mag < sTolerance) {
				m = Matrix4x4d<T>::identity;
            }
			z /= mag;
			Vector3d<T> x = Vector3d<T>::Cross(upVec, z);
			mag = x.Magnitude();
            if (mag < sTolerance) {
				m = Matrix4x4d<T>::identity;
            }
			x /= mag;
			Vector3d<T> y = Vector3d<T>::Cross(z, x);
			m(0, 0) = x.x;
			m(1, 0) = x.y;
			m(2, 0) = x.z;
			m(0, 1) = y.x;
			m(1, 1) = y.y;
			m(2, 1) = y.z;
			m(0, 2) = z.x;
			m(1, 2) = z.y;
			m(2, 2) = z.z;

			return m;
        }
    };
    template <typename T>
    const Quaterniond<T> Quaterniond<T>::identity{ 0, 0, 0, 1 };




} // namespace XYZ

#endif // QUATERNIOND_H


