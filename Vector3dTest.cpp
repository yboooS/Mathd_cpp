// ConsoleApplication2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "Vector3d.h"
using namespace XYZ;

using namespace std;

int main()
{

    // 1. 测试默认构造函数
    Vector3d<double> v0;
    std::cout << "Test 1 - Default Constructor: v0 = (" << v0.x << ", " << v0.y << ", " << v0.z << ")\n";

    // 2. 测试参数化构造函数
    Vector3d<double> v1(3.0, 4.0, 5.0);
    std::cout << "Test 2 - Parameterized Constructor: v1 = (" << v1.x << ", " << v1.y << ", " << v1.z << ")\n";

    // 3. 测试拷贝构造函数
    Vector3d<double> v2 = v1;
    std::cout << "Test 3 - Copy Constructor: v2 = (" << v2.x << ", " << v2.y << ", " << v2.z << ")\n";

    // 4. 测试operator[]
    std::cout << "Test 4 - operator[]: v1[0] = " << v1[0] << ", v1[1] = " << v1[1] << ", v1[2] = " << v1[2] << "\n";

    // 5. 测试operator=
    Vector3d<double> v3;
    v3 = v1;
    std::cout << "Test 5 - operator=: v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 6. 测试operator+=
    v3 += v1;
    std::cout << "Test 6 - operator+=: v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 7. 测试operator-=
    v3 -= v1;
    std::cout << "Test 7 - operator-=: v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 8. 测试operator*=
    v3 *= 2;
    std::cout << "Test 8 - operator*=: v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 9. 测试operator/=
    v3 /= 2;
    std::cout << "Test 9 - operator/=: v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 10. 测试operator+
    Vector3d<double> v4 = v1 + v2;
    std::cout << "Test 10 - operator+: v4 = (" << v4.x << ", " << v4.y << ", " << v4.z << ")\n";

    // 11. 测试operator-
    Vector3d<double> v5 = v1 - v2;
    std::cout << "Test 11 - operator-: v5 = (" << v5.x << ", " << v5.y << ", " << v5.z << ")\n";

    // 12. 测试operator*
    Vector3d<double> v6 = v1 * 2;
    std::cout << "Test 12 - operator*: v6 = (" << v6.x << ", " << v6.y << ", " << v6.z << ")\n";

    // 13. 测试operator/
    Vector3d<double> v7 = v1 / 2;
    std::cout << "Test 13 - operator/: v7 = (" << v7.x << ", " << v7.y << ", " << v7.z << ")\n";

    // 14. 测试operator==
    bool isEqual = v1 == v2;
    std::cout << "Test 14 - operator==: isEqual = " << isEqual << "\n";

    // 15. 测试operator!=
    bool isNotEqual = v1 != v3;
    std::cout << "Test 15 - operator!=: isNotEqual = " << isNotEqual << "\n";



    // 18. 测试Set()
    v3.Set(7.0, 8.0, 9.0);
    std::cout << "Test 18 - Set(): v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 19. 测试Negate()
    v3.Negate();
    std::cout << "Test 19 - Negate(): v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 20. 测试Normalize()
    v3.Normalize();
    std::cout << "Test 20 - Normalize(): v3 = (" << v3.x << ", " << v3.y << ", " << v3.z << ")\n";

    // 21. 测试Magnitude()
    double mag = v1.Magnitude();
    std::cout << "Test 21 - Magnitude(): mag = " << mag << "\n";

    // 22. 测试SqrMagnitude()
    double sqrMag = v1.SqrMagnitude();
    std::cout << "Test 22 - SqrMagnitude(): sqrMag = " << sqrMag << "\n";

    // 23. 测试Scale()
    v1.Scale(v2);
    std::cout << "Test 23 - Scale(): v1 = (" << v1.x << ", " << v1.y << ", " << v1.z << ")\n";

    // 24. 测试Dot()
    double dot = Vector3d<double>::Dot(v1, v2);
    std::cout << "Test 24 - Dot(): dot = " << dot << "\n";

    // 25. 测试Cross()
    Vector3d<double> cross = Vector3d<double>::Cross(v1, v2);
    std::cout << "Test 25 - Cross(): cross = (" << cross.x << ", " << cross.y << ", " << cross.z << ")\n";

    // 26. 测试Distance()
    double dist = Vector3d<double>::Distance(v1, v2);
    std::cout << "Test 26 - Distance(): distance = " << dist << "\n";

    // 27. 测试Normalize() (static)
    Vector3d<double> norm = Vector3d<double>::Normalize(v1);
    std::cout << "Test 27 - Normalize() (static): norm = (" << norm.x << ", " << norm.y << ", " << norm.z << ")\n";

    // 28. 测试Magnitude() (static)
    double magStatic = Vector3d<double>::Magnitude(v1);
    std::cout << "Test 28 - Magnitude() (static): magStatic = " << magStatic << "\n";

    // 29. 测试SqrMagnitude() (static)
    double sqrMagStatic = Vector3d<double>::SqrMagnitude(v1);
    std::cout << "Test 29 - SqrMagnitude() (static): sqrMagStatic = " << sqrMagStatic << "\n";

    // 30. 测试Dot() (static)
    double dotStatic = Vector3d<double>::Dot(v1, v2);
    std::cout << "Test 30 - Dot() (static): dotStatic = " << dotStatic << "\n";

    // 31. 测试Cross() (static)
    Vector3d<double> crossStatic = Vector3d<double>::Cross(v1, v2);
    std::cout << "Test 31 - Cross() (static): crossStatic = (" << crossStatic.x << ", " << crossStatic.y << ", " << crossStatic.z << ")\n";

    // 32. 测试Distance() (static)
    double distStatic = Vector3d<double>::Distance(v1, v2);
    std::cout << "Test 32 - Distance() (static): distStatic = " << distStatic << "\n";

    // 33. 测试Angle()
    double angle = Vector3d<double>::Angle(v1, v2);
    std::cout << "Test 33 - Angle(): angle = " << angle << "\n";

    // 34. 测试AngleBetween()
    double angleBetween = Vector3d<double>::AngleBetween(v1, v2);
    std::cout << "Test 34 - AngleBetween(): angleBetween = " << angleBetween << "\n";

    // 35. 测试ClampMagnitude()
    Vector3d<double> clamped = Vector3d<double>::ClampMagnitude(v1, 1.0);
    std::cout << "Test 35 - ClampMagnitude(): clamped = (" << clamped.x << ", " << clamped.y << ", " << clamped.z << ")\n";

    // 36. 测试Project()
    Vector3d<double> projected = Vector3d<double>::Project(v1, v2);
    std::cout << "Test 36 - Project(): projected = (" << projected.x << ", " << projected.y << ", " << projected.z << ")\n";

    // 37. 测试Exclude()
    Vector3d<double> excluded = Vector3d<double>::Exclude(v1, v2);
    std::cout << "Test 37 - Exclude(): excluded = (" << excluded.x << ", " << excluded.y << ", " << excluded.z << ")\n";

    // 38. 测试Lerp()
    Vector3d<double> lerped = Vector3d<double>::Lerp(v1, v2, 0.5);
    std::cout << "Test 38 - Lerp(): lerped = (" << lerped.x << ", " << lerped.y << ", " << lerped.z << ")\n";

    // 39. 测试LerpUnclamped()
    Vector3d<double> lerpedUnclamped = Vector3d<double>::LerpUnclamped(v1, v2, 0.5);
    std::cout << "Test 39 - LerpUnclamped(): lerpedUnclamped = (" << lerpedUnclamped.x << ", " << lerpedUnclamped.y << ", " << lerpedUnclamped.z << ")\n";

    // 40. 测试Max()
    Vector3d<double> max = Vector3d<double>::Max(v1, v2);
    std::cout << "Test 40 - Max(): max = (" << max.x << ", " << max.y << ", " << max.z << ")\n";

    // 41. 测试Min()
    Vector3d<double> min = Vector3d<double>::Min(v1, v2);
    std::cout << "Test 41 - Min(): min = (" << min.x << ", " << min.y << ", " << min.z << ")\n";

    // 42. 测试MoveTowards()
    Vector3d<double> moved = Vector3d<double>::MoveTowards(v1, v2, 0.1);
    std::cout << "Test 42 - MoveTowards(): moved = (" << moved.x << ", " << moved.y << ", " << moved.z << ")\n";

    // 43. 测试OrthoNormalize()
    Vector3d<double> tangent(1, 1, 0);
    Vector3d<double> normalizedNormal = v1;
    Vector3d<double>::OrthoNormalize(normalizedNormal, tangent);
    std::cout << "Test 43 - OrthoNormalize(): normalizedNormal = (" << normalizedNormal.x << ", " << normalizedNormal.y << ", " << normalizedNormal.z << ")\n";
    std::cout << "Tangent after OrthoNormalize = (" << tangent.x << ", " << tangent.y << ", " << tangent.z << ")\n";

    // 44. 测试OrthoNormalize()
    Vector3d<double> tangent2(1, 1, 0);
    Vector3d<double> normalizedNormal2 = v1;
    Vector3d<double> binormal(2, 1, 0);;
    Vector3d<double>::OrthoNormalize(normalizedNormal2, tangent2, binormal);
    std::cout << "Test 44 - OrthoNormalize(): normalizedNormal = (" << normalizedNormal2.x << ", " << normalizedNormal2.y << ", " << normalizedNormal2.z << ")\n";
    std::cout << "Tangent after OrthoNormalize = (" << tangent2.x << ", " << tangent2.y << ", " << tangent2.z << ")\n";
    std::cout << "Tangent after OrthoNormalize = (" << binormal.x << ", " << binormal.y << ", " << binormal.z << ")\n";

    // 45. 测试ProjectOnPlane()
    Vector3d<double> projectedOnPlane = Vector3d<double>::ProjectOnPlane(v1, v2);
    std::cout << "Test 45 - ProjectOnPlane(): projectedOnPlane = (" << projectedOnPlane.x << ", " << projectedOnPlane.y << ", " << projectedOnPlane.z << ")\n";

    // 46. 测试Reflect()
    Vector3d<double> reflected = Vector3d<double>::Reflect(v1, v2);
    std::cout << "Test 46 - Reflect(): reflected = (" << reflected.x << ", " << reflected.y << ", " << reflected.z << ")\n";

    // 47. 测试RotateTowards()
    Vector3d<double> rotated = Vector3d<double>::RotateTowards(v1, v2, PI / 4, 10.0);
    std::cout << "Test 47 - RotateTowards(): rotated = (" << rotated.x << ", " << rotated.y << ", " << rotated.z << ")\n";

    // 48. 测试Scale() (static)
    Vector3d<double> scaled = Vector3d<double>::Scale(v1, v2);
    std::cout << "Test 48 - Scale() (static): scaled = (" << scaled.x << ", " << scaled.y << ", " << scaled.z << ")\n";

    // 49. 测试Slerp()
    Vector3d<double> slerped = Vector3d<double>::Slerp(v1, v2, 0.5);
    std::cout << "Test 49 - Slerp(): slerped = (" << slerped.x << ", " << slerped.y << ", " << slerped.z << ")\n";

    // 50. 测试SlerpUnclamped()
    Vector3d<double> slerpedUnclamped = Vector3d<double>::SlerpUnclamped(v1, v2, 0.5);
    std::cout << "Test 50 - SlerpUnclamped(): slerpedUnclamped = (" << slerpedUnclamped.x << ", " << slerpedUnclamped.y << ", " << slerpedUnclamped.z << ")\n";

    // 51. 测试ToString()
    std::string v1Str = v1.ToString();
    std::cout << "Test 51 - ToString(): " << v1Str << "\n";

    // 52. 测试GetHashCode()
    size_t hash = v1.GetHashCode();
    std::cout << "Test 52 - GetHashCode(): " << hash << "\n";

    // 53. 测试Equals()
    bool equals = v1.Equals(v2);
    std::cout << "Test 53 - Equals(): " << equals << "\n";

    // 54. 测试ToString() with format
    std::string formattedString = v1.ToString("F2");
    std::cout << "Test 54 - ToString() with format: " << formattedString << "\n";

    //// 55. 测试OrthoNormalVectorFast()
    //Vector3d<double> orthoNormal = Vector3d<double>::OrthoNormalVectorFast(v1);
    //std::cout << "Test 55 - OrthoNormalVectorFast(): orthoNormal = (" << orthoNormal.x << ", " << orthoNormal.y << ", " << orthoNormal.z << ")\n";

    //// 56. 测试SetAxisAngle()
    //Vector3d<double> axis(0, 1, 0);
    //Matrix4x4d<double> rotationMatrix = Matrix4x4d<double>::SetAxisAngle(axis, PI / 2);
    //std::cout << "Test 56 - SetAxisAngle(): rotationMatrix elements = ... (simplified output)\n";

    // 57. 未来可能添加的其他测试项
    // 确保所有函数都进行了测试
    std::cout << "Test 57 - Placeholder for future tests\n";
    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
