//
//  Matrix4x4d.hpp
//  XYZ
//  Unity Matrix4x4常用类模板C++实现。封装了常用成员与方法。
//
//  Created by yangbo on 24-08-08.
#ifndef MATRIX4X4D_H
#define MATRIX4X4D_H

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iomanip> 
#include "math_utils.hpp"
#include "Vector3d.hpp"
#include "Vector4d.hpp"
#include "Quaterniond.hpp"

namespace XYZ {

    template <typename T>
    class Matrix4x4d {
    public:
        T m00, m01, m02, m03;
        T m10, m11, m12, m13;
        T m20, m21, m22, m23;
        T m30, m31, m32, m33;
        // 构造函数
        Matrix4x4d(T m00 = T(0), T m01 = T(0), T m02 = T(0), T m03 = T(0),
            T m10 = T(0), T m11 = T(0), T m12 = T(0), T m13 = T(0),
            T m20 = T(0), T m21 = T(0), T m22 = T(0), T m23 = T(0),
            T m30 = T(0), T m31 = T(0), T m32 = T(0), T m33 = T(0))
            : m00(m00), m01(m01), m02(m02), m03(m03),
            m10(m10), m11(m11), m12(m12), m13(m13),
            m20(m20), m21(m21), m22(m22), m23(m23),
            m30(m30), m31(m31), m32(m32), m33(m33) {}
		// 析构函数
        ~Matrix4x4d() {}
        // 静态常量
        static const Matrix4x4d<T> identity;
        static const Matrix4x4d<T> zero;
        // 索引器
        T& operator[](int index) {
            switch (index) {
            case 0: return m00;
            case 1: return m10;
            case 2: return m20;
            case 3: return m30;
            case 4: return m01;
            case 5: return m11;
            case 6: return m21;
            case 7: return m31;
            case 8: return m02;
            case 9: return m12;
            case 10: return m22;
            case 11: return m32;
            case 12: return m03;
            case 13: return m13;
            case 14: return m23;
            case 15: return m33;
            default: throw std::out_of_range("Invalid matrixd index!");
            }
        }

        T& operator()(int row, int column) {
            return (*this)[row + column * 4];
        }
        // 设置一维索引值
        void set(int index, T value) {
            (*this)[index] = value;
        }

        // 设置二维索引值
        void set(int row, int column, T value) {
            (*this)(row, column) = value;
        }

        // 左乘 友元函数
        friend Vector4d<T> operator*(const Matrix4x4d<T>& lhs, const Vector4d<T>& v) {
            Vector4d<T> result;
            result.x = lhs.m00 * v.x + lhs.m01 * v.y + lhs.m02 * v.z + lhs.m03 * v.w;
            result.y = lhs.m10 * v.x + lhs.m11 * v.y + lhs.m12 * v.z + lhs.m13 * v.w;
            result.z = lhs.m20 * v.x + lhs.m21 * v.y + lhs.m22 * v.z + lhs.m23 * v.w;
            result.w = lhs.m30 * v.x + lhs.m31 * v.y + lhs.m32 * v.z + lhs.m33 * v.w;
            return result;
        }
        // 相乘 友元函数
        friend Matrix4x4d<T> operator*(const Matrix4x4d<T>& lhs, const Matrix4x4d<T>& rhs) {
            Matrix4x4d<T> result;
            result.m00 = lhs.m00 * rhs.m00 + lhs.m01 * rhs.m10 + lhs.m02 * rhs.m20 + lhs.m03 * rhs.m30;
            result.m01 = lhs.m00 * rhs.m01 + lhs.m01 * rhs.m11 + lhs.m02 * rhs.m21 + lhs.m03 * rhs.m31;
            result.m02 = lhs.m00 * rhs.m02 + lhs.m01 * rhs.m12 + lhs.m02 * rhs.m22 + lhs.m03 * rhs.m32;
            result.m03 = lhs.m00 * rhs.m03 + lhs.m01 * rhs.m13 + lhs.m02 * rhs.m23 + lhs.m03 * rhs.m33;
            result.m10 = lhs.m10 * rhs.m00 + lhs.m11 * rhs.m10 + lhs.m12 * rhs.m20 + lhs.m13 * rhs.m30;
            result.m11 = lhs.m10 * rhs.m01 + lhs.m11 * rhs.m11 + lhs.m12 * rhs.m21 + lhs.m13 * rhs.m31;
            result.m12 = lhs.m10 * rhs.m02 + lhs.m11 * rhs.m12 + lhs.m12 * rhs.m22 + lhs.m13 * rhs.m32;
            result.m13 = lhs.m10 * rhs.m03 + lhs.m11 * rhs.m13 + lhs.m12 * rhs.m23 + lhs.m13 * rhs.m33;
            result.m20 = lhs.m20 * rhs.m00 + lhs.m21 * rhs.m10 + lhs.m22 * rhs.m20 + lhs.m23 * rhs.m30;
            result.m21 = lhs.m20 * rhs.m01 + lhs.m21 * rhs.m11 + lhs.m22 * rhs.m21 + lhs.m23 * rhs.m31;
            result.m22 = lhs.m20 * rhs.m02 + lhs.m21 * rhs.m12 + lhs.m22 * rhs.m22 + lhs.m23 * rhs.m32;
            result.m23 = lhs.m20 * rhs.m03 + lhs.m21 * rhs.m13 + lhs.m22 * rhs.m23 + lhs.m23 * rhs.m33;
            result.m30 = lhs.m30 * rhs.m00 + lhs.m31 * rhs.m10 + lhs.m32 * rhs.m20 + lhs.m33 * rhs.m30;
            result.m31 = lhs.m30 * rhs.m01 + lhs.m31 * rhs.m11 + lhs.m32 * rhs.m21 + lhs.m33 * rhs.m31;
            result.m32 = lhs.m30 * rhs.m02 + lhs.m31 * rhs.m12 + lhs.m32 * rhs.m22 + lhs.m33 * rhs.m32;
            result.m33 = lhs.m30 * rhs.m03 + lhs.m31 * rhs.m13 + lhs.m32 * rhs.m23 + lhs.m33 * rhs.m33;
            return result;
        }
        // 相等
        friend bool operator==(const Matrix4x4d<T>& lhs, const Matrix4x4d<T>& rhs) {
            return lhs.m00 == rhs.m00 && lhs.m01 == rhs.m01 && lhs.m02 == rhs.m02 && lhs.m03 == rhs.m03 &&
				lhs.m10 == rhs.m10 && lhs.m11 == rhs.m11 && lhs.m12 == rhs.m12 && lhs.m13 == rhs.m13 &&
				lhs.m20 == rhs.m20 && lhs.m21 == rhs.m21 && lhs.m22 == rhs.m22 && lhs.m23 == rhs.m23 &&
				lhs.m30 == rhs.m30 && lhs.m31 == rhs.m31 && lhs.m32 == rhs.m32 && lhs.m33 == rhs.m33;
        }
        // 不等
        friend bool operator!=(const Matrix4x4d<T>& lhs, const Matrix4x4d<T>& rhs) {
			return !(lhs == rhs);
		}
    public:
        // 行列式
        T determinant() const {
            T result =
                m03 * m12 * m21 * m30 - m02 * m13 * m21 * m30 - m03 * m11 * m22 * m30 + m01 * m13 * m22 * m30 +
                m02 * m11 * m23 * m30 - m01 * m12 * m23 * m30 - m03 * m12 * m20 * m31 + m02 * m13 * m20 * m31 +
                m03 * m10 * m22 * m31 - m00 * m13 * m22 * m31 - m02 * m10 * m23 * m31 + m00 * m12 * m23 * m31 +
                m03 * m11 * m20 * m32 - m01 * m13 * m20 * m32 - m03 * m10 * m21 * m32 + m00 * m13 * m21 * m32 +
                m01 * m10 * m23 * m32 - m00 * m11 * m23 * m32 - m02 * m11 * m20 * m33 + m01 * m12 * m20 * m33 +
                m02 * m10 * m21 * m33 - m00 * m12 * m21 * m33 - m01 * m10 * m22 * m33 + m00 * m11 * m22 * m33;
            return result;
        }

        // 逆矩阵
        Matrix4x4d<T> inverse() const {
            int m = 4;
            int n = 4;
            std::array<std::array<T, 2 * 4 + 1>, 2 * 4 + 1> array{};

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    array[i][j] = (*this)(i, j);
                }
            }

            for (int k = 0; k < m; k++) {
                for (int t = n; t <= 2 * n; t++) {
                    array[k][t] = (t - k == m) ? T(1) : T(0);
                }
            }

            for (int k = 0; k < m; k++) {
                if (array[k][k] != T(1)) {
                    T bs = array[k][k];
                    array[k][k] = T(1);
                    for (int p = k + 1; p < 2 * n; p++) {
                        array[k][p] /= bs;
                    }
                }
                for (int q = 0; q < m; q++) {
                    if (q != k) {
                        T bs = array[q][k];
                        for (int p = 0; p < 2 * n; p++) {
                            array[q][p] -= bs * array[k][p];
                        }
                    }
                }
            }

            Matrix4x4d<T> result;
            for (int x = 0; x < m; x++) {
                for (int y = n; y < 2 * n; y++) {
                    result(x, y - n) = array[x][y];
                }
            }
            return result;
        }

        // 是否为单位矩阵
        bool isIdentity() const {
            return *this == identity;
        }

        // 转置矩阵
        Matrix4x4d<T> transpose() const {
            Matrix4x4d<T> temp;
            temp.m00 = m00; temp.m10 = m01; temp.m20 = m02; temp.m30 = m03;
            temp.m01 = m10; temp.m11 = m11; temp.m21 = m12; temp.m31 = m13;
            temp.m02 = m20; temp.m12 = m21; temp.m22 = m22; temp.m32 = m23;
            temp.m03 = m30; temp.m13 = m31; temp.m23 = m32; temp.m33 = m33;
            return temp;
        }
    public:
        // 行列式
        static T Determinant(const Matrix4x4d<T>& m) {
            return m.determinant();
        }
        // 逆矩阵
        static Matrix4x4d<T> Inverse(const Matrix4x4d<T>& m) {
            return m.inverse();
        }
        // 转置矩阵
        static Matrix4x4d<T> Transpose(const Matrix4x4d<T>& m) {
			return m.transpose();
		}
        // (未实现)
        static Matrix4x4d<T> LookAt(const Vector3d<T>& eye, const Vector3d<T>& target, const Vector3d<T>& up) {
			Vector3d<T> zaxis = (eye - target).normalize();
			Vector3d<T> xaxis = up.cross(zaxis).normalize();
			Vector3d<T> yaxis = zaxis.cross(xaxis);

			Matrix4x4d<T> result;
			result.m00 = xaxis.x; result.m01 = yaxis.x; result.m02 = zaxis.x; result.m03 = T(0);
			result.m10 = xaxis.y; result.m11 = yaxis.y; result.m12 = zaxis.y; result.m13 = T(0);
			result.m20 = xaxis.z; result.m21 = yaxis.z; result.m22 = zaxis.z; result.m23 = T(0);
			result.m30 = -xaxis.dot(eye); result.m31 = -yaxis.dot(eye); result.m32 = -zaxis.dot(eye); result.m33 = T(1);
			return result;
		}
        // 正交矩阵
        static Matrix4x4d<T> Ortho(const T left, const T right, const T bottom, const T top, const T zNear, const T zFar) {
            Matrix4x4d<T> result = identity;

            T deltax = right - left;
            T deltay = top - bottom;
            T deltaz = zFar - zNear;

            result(0, 0) = 2.0 / deltax;
            result(0, 3) = -(right + left) / deltax;
            result(1, 1) = 2.0 / deltay;
            result(1, 3) = -(top + bottom) / deltay;
            result(2, 2) = -2.0 / deltaz;
            result(2, 3) = -(zFar + zNear) / deltaz;
            return result;
        }
        // 透视矩阵
        static Matrix4x4d<T> Perspective(const T fov, const T aspect, const T zNear, const T zFar) {
			Matrix4x4d<T> result = zero;

            T cotangent, deltaZ;
            T radians = (fov / 2.0) * (PI / 180);
            cotangent = std::cos(radians) / std::sin(radians);
            deltaZ = zNear - zFar;

            result(0, 0) = cotangent / aspect; result(0, 1) = 0; result(0, 2) = 0; result(0, 3) = 0;
            result(1, 0) = 0; result(1, 1) = cotangent; result(1, 2) = 0; result(1, 3) = 0;
            result(2, 0) = 0; result(2, 1) = 0; result(2, 2) = (zFar + zNear) / deltaZ; result(2, 3) = 2 * zNear * zFar / deltaZ;
            result(3, 0) = 0; result(3, 1) = 0; result(3, 2) = -1; result(3, 3) = 0;

			return result;
		}
        // 缩放
        static Matrix4x4d<T> Scale(const Vector3d<T>& v) {
            Matrix4x4d<T> result;
            result.m00 = v.x; result.m01 = T(0); result.m02 = T(0); result.m03 = T(0);
            result.m10 = T(0); result.m11 = v.y; result.m12 = T(0); result.m13 = T(0);
            result.m20 = T(0); result.m21 = T(0); result.m22 = v.z; result.m23 = T(0);
            result.m30 = T(0); result.m31 = T(0); result.m32 = T(0); result.m33 = T(1);
			return result;
		}
        // 
        static Matrix4x4d<T> Translate(const Vector3d<T>& v) {
            Matrix4x4d<T> result;
            result.m00 = T(1); result.m01 = T(0); result.m02 = T(0); result.m03 = v.x;
            result.m10 = T(0); result.m11 = T(1); result.m12 = T(0); result.m13 = v.y;
            result.m20 = T(0); result.m21 = T(0); result.m22 = T(1); result.m23 = v.z;
            result.m30 = T(0); result.m31 = T(0); result.m32 = T(0); result.m33 = T(1);
            return result;
        }
        // 平移旋转缩放矩阵
        static Matrix4x4d<T> TRS(const Vector3d<T>& pos, const Quaterniond<T>& q, const Vector3d<T>& s) {
			Matrix4x4d<T> m = Quaterniond<T>::QuaternionToMatrix(q);
            m[0] *= s[0]; m[1] *= s[0]; m[2] *= s[0];
            m[4] *= s[1]; m[5] *= s[1]; m[6] *= s[1];
            m[8] *= s[2]; m[9] *= s[2]; m[10] *= s[2];
            m[12] = pos[0]; m[13] = pos[1]; m[14] = pos[2];
			return m;
		}
        // 获取列
        Vector4d<T> GetColumn(const T index) const {
			switch (index) {
			case 0: return Vector4d<T>(m00, m10, m20, m30);
			case 1: return Vector4d<T>(m01, m11, m21, m31);
			case 2: return Vector4d<T>(m02, m12, m22, m32);
			case 3: return Vector4d<T>(m03, m13, m23, m33);
			default: throw std::out_of_range("Invalid matrixd column index!");
			}
		}
        // 获取行
        Vector4d<T> GetRow(const T index) const {
            switch (index) {
            case 0: return Vector4d<T>(m00, m01, m02, m03);
            case 1: return Vector4d<T>(m10, m11, m12, m13);
            case 2: return Vector4d<T>(m20, m21, m22, m23);
            case 3: return Vector4d<T>(m30, m31, m32, m33);
            default: throw std::out_of_range("Invalid matrixd row index!");
            }
        }
        // 乘
        Vector3d<T> MultiplyPoint(const Vector3d<T>& v) {
            Vector3d<T> result;
            result.x = m00 * v.x + m01 * v.y + m02 * v.z + m03;
            result.y = m10 * v.x + m11 * v.y + m12 * v.z + m13;
            result.z = m20 * v.x + m21 * v.y + m22 * v.z + m23;
            T num = m30 * v.x + m31 * v.y + m32 * v.z + m33;
            num = 1.0 / num;
            result.x *= num;
            result.y *= num;
            result.z *= num;
            return result;
        }
        // 乘3x4矩阵
        Vector3d<T> MultiplyPoint3x4(const Vector3d<T>& v) {
            Vector3d<T> result;
            result.x = m00 * v.x + m01 * v.y + m02 * v.z + m03;
            result.y = m10 * v.x + m11 * v.y + m12 * v.z + m13;
            result.z = m20 * v.x + m21 * v.y + m22 * v.z + m23;
            return result;
        }
        // 乘向量
        Vector3d<T> MultiplyVector(const Vector3d<T>& v) {
            Vector3d<T> result;
            result.x = m00 * v.x + m01 * v.y + m02 * v.z;
            result.y = m10 * v.x + m11 * v.y + m12 * v.z;
            result.z = m20 * v.x + m21 * v.y + m22 * v.z;
            return result;
        }
        // 设置列
        void SetColumn(const T i, const Vector4d<T>& v) {
            *this(0, i) = v.x;
            *this(1, i) = v.x;
            *this(2, i) = v.x;
            *this(3, i) = v.x;
		}
        // 设置行
        void SetRow(const T i, const Vector4d<T>& v) {
            *this(i, 0) = v.x;
            *this(i, 1) = v.x;
            *this(i, 2) = v.x;
            *this(i, 3) = v.x;
        }
        // 设置平移旋转缩放矩阵
        void SetTRS(const Vector3d<T>& pos, const Quaterniond<T>& q, const Vector3d<T>& s) {
			*this = TRS(pos, q, s);
		}
        // 格式化字符串
        std::string ToString() {
            std::stringstream ss;
            ss << "(" << m00 << ", " << m01 << ", " << m02 << ", " << m03 << "," << m10 << ", " << m11 << ", " << m12 << ", " << m13 << "," << m20 << ", " << m21 << ", " << m22 << ", " << m23 << "," << m30 << ", " << m31 << ", " << m32 << ", " << m33 << ")";
            return ss.str();
        }
        // 哈希值
        size_t GetHashCode() {
            size_t hash = 17;
            hash = hash * 31 + std::hash<double>{}(m00);
            hash = hash * 31 + std::hash<double>{}(m01);
            hash = hash * 31 + std::hash<double>{}(m02);
            hash = hash * 31 + std::hash<double>{}(m03);
            hash = hash * 31 + std::hash<double>{}(m10);
            hash = hash * 31 + std::hash<double>{}(m11);
            hash = hash * 31 + std::hash<double>{}(m12);
            hash = hash * 31 + std::hash<double>{}(m13);
            hash = hash * 31 + std::hash<double>{}(m20);
            hash = hash * 31 + std::hash<double>{}(m21);
            hash = hash * 31 + std::hash<double>{}(m22);
            hash = hash * 31 + std::hash<double>{}(m23);
            hash = hash * 31 + std::hash<double>{}(m30);
            hash = hash * 31 + std::hash<double>{}(m31);
            hash = hash * 31 + std::hash<double>{}(m32);
            hash = hash * 31 + std::hash<double>{}(m33);

            return hash;
        }
        // 判断是否相等
        bool Equals(const Vector4d<T>& other) {
            return *this == other;
        }
        // 格式化字符串
        std::string ToString(const std::string& format) {
            std::stringstream ss;
            ss << "(";
            if (format == "F1") {
                ss << std::fixed << std::setprecision(1) << m00 << ", " << m01 << ", " << m02 << ", " << m03 << "," << m10 << ", " << m11 << ", " << m12 << ", " << m13 << "," << m20 << ", " << m21 << ", " << m22 << ", " << m23 << "," << m30 << ", " << m31 << ", " << m32 << ", " << m33;
            }
            else if (format == "F2") {
                ss << std::fixed << std::setprecision(2) << m00 << ", " << m01 << ", " << m02 << ", " << m03 << "," << m10 << ", " << m11 << ", " << m12 << ", " << m13 << "," << m20 << ", " << m21 << ", " << m22 << ", " << m23 << "," << m30 << ", " << m31 << ", " << m32 << ", " << m33;
            }
            else if (format == "F3") {
                ss << std::fixed << std::setprecision(3) << m00 << ", " << m01 << ", " << m02 << ", " << m03 << "," << m10 << ", " << m11 << ", " << m12 << ", " << m13 << "," << m20 << ", " << m21 << ", " << m22 << ", " << m23 << "," << m30 << ", " << m31 << ", " << m32 << ", " << m33;
            }
            ss << ")";
            return ss.str();
        }
    };
    template <typename T>
    const Matrix4x4d<T> Matrix4x4d<T>::identity = Matrix4x4d<T>(1, 0, 0, 0,
                                                                0, 1, 0, 0,
                                                                0, 0, 1, 0,
                                                                0, 0, 0, 1);

    template <typename T>
    const Matrix4x4d<T> Matrix4x4d<T>::zero = Matrix4x4d<T>(0, 0, 0, 0,
                                                            0, 0, 0, 0,
                                                            0, 0, 0, 0,
                                                            0, 0, 0, 0);
} // namespace MX4

#endif // MATRIX4X4D_H




