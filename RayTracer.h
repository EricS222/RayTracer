#pragma once
#include <vector>
#include <string>

using namespace std;

class RayTracer
{
	public:

		struct PIXEL
		{
			PIXEL(int r, int g, int b)
			{
				this->Red = r;
				this->Green = g;
				this->Blue = b;
			}
			float Red;
			float Green;
			float Blue;
		};

		struct Vector
		{
			Vector() {}
			Vector(float xVal, float yVal, float zVal)
			{
				this->x = xVal;
				this->y = yVal;
				this->z = zVal;
			}
			float x;
			float y;
			float z;
		};

		struct Ray
		{
			Ray() {}
			Ray(Vector& o, Vector& d, float size) : origin(o), direction(d), t(size) {}
			Vector origin;
			Vector direction;
			float t;
		};

		struct Sphere
		{
			Sphere(Vector o, float r)
			{
				this->origin = o;
				this->radius = r;
			}
			Vector origin;
			float radius;

		};

		struct Triangle
		{
			Triangle(Vector A, Vector B, Vector C)
			{
				this->a = A;
				this->b = B;
				this->c = C;
			}
			Vector a;
			Vector b;
			Vector c;
		};

		RayTracer();
		float sphereIntersection(Sphere sphere, Ray ray);
		float sphereShading(float diffuseCoefficient, float intensity, float t, float specular_coefficient, float phong_exponent, float ambient1, float ambient2, Vector source, Ray ray, Sphere sphere);
		float triangleShading(float diffuseCoefficient, float intensity, float t, float specular_coefficient, float phong_exponent, float ambient1, float ambient2, Vector source, Ray ray, Triangle tetrahedron);
		float triangleIntersection(Triangle tetrahedron, Ray ray);
		Vector crossProduct(Vector a, Vector b)
		{
			Vector cross_product(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
			return cross_product;
		}
		Vector normalVector(Vector a)
		{
			float one = 1;
			float two = -1;
			float x = (a.z * two - a.y * one) / a.x;
			Vector output(x, one, two);
			float output_length = sqrt(dotProduct(output, output));
			Vector output_normalized(output.x / output_length, output.y / output_length, output.z / output_length);
			return output_normalized;
		}
		float dotProduct(Vector a, Vector b);
		vector<float> getPixels();

	private:
		vector<vector<PIXEL>> screen;
		vector<float> pixels;
		unsigned char * RGB;




};