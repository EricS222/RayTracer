#include "RayTracer.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

RayTracer::RayTracer()
{
	//i,j viewing rectangle
	float width = 400;
	float height = 400;

	//set up objects
	Sphere sphere1(Vector(0.65, 0.25, 0), 0.14);
	Sphere sphere2(Vector(0.78, 0.52, 0), 0.12);
	Triangle tetrahedron1(Vector(0.30, 0.6, -1), Vector(0.30, 0.40, 0), Vector(0.10, 0.35, -1));
	Triangle tetrahedron2(Vector(0.50, 0.35, -1), Vector(0.30, 0.40, 0), Vector(0.30, 0.6, -1));
	Triangle tetrahedron3(Vector(0.30, 0.40, 0), Vector(0.50, 0.35, -1), Vector(0.10, 0.35, -1));
	Triangle tetrahedron4(Vector(0.5, 0.35, -1), Vector(0.3, 0.6, -1), Vector(0.10, 0.35, -1));

	for (int i = 0; i < height; i++)
	{
		screen.push_back(vector<PIXEL>());
		for (int j = 0; j < width; j++)
		{
			/*if (i > 200 && (j == 320 || j == 319 || j == 321))
			{
				PIXEL pixel(0, 0, 0);
				screen[i].push_back(pixel);
			}
			else if (i > 100 && (j == 260 || j == 261 || j == 259))
			{
				PIXEL pixel(0, 0, 0);
				screen[i].push_back(pixel);
			}
			else if (i > 200 && (j == 120 || j == 121 || j == 119))
			{
				PIXEL pixel(0, 0, 0);
				screen[i].push_back(pixel);
			}
			else
			{*/
				PIXEL pixel(150, 150, 150);
				screen[i].push_back(pixel);
			//}
		}
	}

	Vector LookAt(0, 0, -1);
	Vector Viewpoint;
	Vector light_source(1, 0, 1);
	float diffuse_coefficient = 0.8;
	float lightSource_intensity = 175;
	float specular_coefficient = 100;
	float phong_exponent = 50;
	float ambient_intensity = 1;
	float ambient_coefficient = 20;
	float mirror_coefficient = 0.9;

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{

			//compute viewing ray
			Ray viewingRay;
			/*float u = ((j + 1 / 2) / (width));
			float v = ((i + 1 / 2) / (height));
			viewingRay.origin.x = 0;
			viewingRay.origin.y = 0;
			viewingRay.origin.z = 1;
			viewingRay.direction.x = u;
			viewingRay.direction.y = v;
			viewingRay.direction.z = -1; */

			viewingRay.origin.x = ((j + 1 / 2) / (width));
			viewingRay.origin.y = ((i + 1/2) / (height));
			viewingRay.origin.z = 0;
			viewingRay.direction = LookAt;

			//intersect ray with scene

			float intersection1 = sphereIntersection(sphere1, viewingRay);
			float intersection2 = sphereIntersection(sphere2, viewingRay);
			float intersection3 = triangleIntersection(tetrahedron1, viewingRay);
			float intersection4 = triangleIntersection(tetrahedron2, viewingRay);
			float intersection5 = triangleIntersection(tetrahedron3, viewingRay);
			float intersection6 = triangleIntersection(tetrahedron4, viewingRay);

			if (intersection1 != -1)
			{
				float RGB = sphereShading(diffuse_coefficient, lightSource_intensity, intersection1, specular_coefficient, phong_exponent, ambient_intensity, ambient_coefficient, light_source, viewingRay, sphere1);
				Vector view_vector(viewingRay.origin.x + intersection1 * viewingRay.direction.x, viewingRay.origin.y + intersection1 * viewingRay.direction.y, viewingRay.origin.z + intersection1 * viewingRay.direction.z);
				Vector normal_vector(view_vector.x - sphere1.origin.x, view_vector.y - sphere1.origin.y, view_vector.z - sphere1.origin.z);
				float normal_length = sqrt(dotProduct(normal_vector, normal_vector));
				Vector normal_normalized(normal_vector.x / normal_length, normal_vector.y / normal_length, normal_vector.z / normal_length);
				Ray temp;
				float constant = (float)-2.0f * (dotProduct(viewingRay.direction, normal_normalized));
				Vector R(constant * normal_normalized.x + viewingRay.direction.x, constant * normal_normalized.y + viewingRay.direction.y, constant * normal_normalized.z + viewingRay.direction.z);
				temp.origin.x = view_vector.x;
				temp.origin.y = view_vector.y;
				temp.origin.z = view_vector.z;
				temp.direction.x = R.x;
				temp.direction.y = R.y;
				temp.direction.z = R.z;
				float intersect = sphereIntersection(sphere2, temp);
				if (intersect != -1)
				{
					if (intersect < 0)
					{
						float RGB2 = sphereShading(diffuse_coefficient, lightSource_intensity, intersect, specular_coefficient, phong_exponent, ambient_intensity, ambient_coefficient, light_source, viewingRay, sphere2);
						screen[i][j].Red = RGB;
						screen[i][j].Blue = 0;
						screen[i][j].Green = RGB2 * mirror_coefficient;
					}
					else
					{
						screen[i][j].Red = RGB;
						screen[i][j].Blue = 0;
						screen[i][j].Green = 0;
					}
				}
				else
				{
					screen[i][j].Red = RGB;
					screen[i][j].Blue = 0;
					screen[i][j].Green = 0;
				}
			}
			else if (intersection2 != -1)
			{
				float RGB = sphereShading(diffuse_coefficient, lightSource_intensity, intersection2, specular_coefficient, phong_exponent, ambient_intensity, ambient_coefficient, light_source, viewingRay, sphere2);
				screen[i][j].Red = 0;
				screen[i][j].Blue = 0;
				screen[i][j].Green = RGB;
			}
			else
			{
				Ray shadowRay;
				shadowRay.origin.x = viewingRay.origin.x;
				shadowRay.origin.y = viewingRay.origin.y;
				shadowRay.origin.z = viewingRay.origin.z - 0.2;
				shadowRay.direction.x = light_source.x;
				shadowRay.direction.y = light_source.y;
				shadowRay.direction.z = light_source.z;
				float shadowIntersection1 = sphereIntersection(sphere1, shadowRay);
				float shadowIntersection2 = sphereIntersection(sphere2, shadowRay);
				if (shadowIntersection1 < 0 && shadowIntersection1 != -1 && screen[i][j].Red != 0)
				{
					screen[i][j].Red = ambient_intensity * ambient_coefficient;
					screen[i][j].Blue = ambient_intensity * ambient_coefficient;
					screen[i][j].Green = ambient_intensity * ambient_coefficient;
				}
				if (shadowIntersection2 < 0 && shadowIntersection2 != -1 && screen[i][j].Red != 0)
				{
					screen[i][j].Red = ambient_intensity*ambient_coefficient;
					screen[i][j].Blue = ambient_intensity * ambient_coefficient;
					screen[i][j].Green = ambient_intensity * ambient_coefficient;
				}
			}

			if (intersection3 != -1)
			{
				float RGB = triangleShading(diffuse_coefficient, lightSource_intensity, intersection3, specular_coefficient, phong_exponent, ambient_intensity, ambient_coefficient, light_source, viewingRay, tetrahedron1);
				screen[i][j].Red = 0;
				screen[i][j].Blue = RGB + 20;
				screen[i][j].Green = 0;
			}

			else if (intersection4 != -1)
			{
				float RGB = triangleShading(diffuse_coefficient, lightSource_intensity, intersection4, specular_coefficient, phong_exponent, ambient_intensity, ambient_coefficient, light_source, viewingRay, tetrahedron2);
				screen[i][j].Red = 0;
				screen[i][j].Blue = RGB + 20;
				screen[i][j].Green = 0;
			}

			else if (intersection5 != -1)
			{
				float RGB = triangleShading(diffuse_coefficient, lightSource_intensity, intersection5, specular_coefficient, phong_exponent, ambient_intensity, ambient_coefficient, light_source, viewingRay, tetrahedron3);
				screen[i][j].Red = 0;
				screen[i][j].Blue = RGB + 20;
				screen[i][j].Green = 0;
			}
			else if (intersection6 != -1)
			{
				float RGB = triangleShading(diffuse_coefficient, lightSource_intensity, intersection6, specular_coefficient, phong_exponent, ambient_intensity, ambient_coefficient, light_source, viewingRay, tetrahedron4);
				screen[i][j].Red = 0;
				screen[i][j].Blue = RGB + 20;
				screen[i][j].Green = 0;
			}
		}
	}
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			pixels.push_back(screen[i][j].Red);
			pixels.push_back(screen[i][j].Green);
			pixels.push_back(screen[i][j].Blue);
		}
	}
}

float RayTracer::dotProduct(Vector a, Vector b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

float RayTracer::sphereIntersection(Sphere sphere, Ray ray)
{
	Vector center = Vector(ray.origin.x - sphere.origin.x, ray.origin.y - sphere.origin.y, ray.origin.z - sphere.origin.z);
	Vector two_d = Vector(2 * ray.direction.x, 2 * ray.direction.y, 2 * ray.direction.z);
	float discriminant = pow(dotProduct(two_d, center), 2) - 4 * 
		dotProduct(ray.direction, ray.direction) * (dotProduct(center, center) - pow(sphere.radius,2));

	float d = sqrt(discriminant);

	if (discriminant < 0)
	{
		return -1;
	}
	else if (discriminant == 0)
	{
		float t1 = (dotProduct(two_d, center) + d) / (2 * (dotProduct(ray.direction, ray.direction)));
		return t1;

	}
	else if (discriminant > 0)
	{
		float t1 = (dotProduct(two_d, center) + d) / (2 * (dotProduct(ray.direction, ray.direction)));
		float t2 = (dotProduct(two_d, center) - d) / (2 * (dotProduct(ray.direction, ray.direction)));
		if (t1 < t2)
		{
			return t1;
		}
		else
		{
			return t2;
		}
	}


}

float RayTracer::sphereShading(float diffuseCoefficient, float intensity, float t, float specular_coefficient, float phong_exponent,float ambient1, float ambient2, Vector source, Ray ray, Sphere sphere)
{
	Vector view_vector(ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, ray.origin.z + t * ray.direction.z);
	Vector normal_vector(view_vector.x - sphere.origin.x, view_vector.y - sphere.origin.y, view_vector.z - sphere.origin.z);
	float normal_length = sqrt(dotProduct(normal_vector, normal_vector));
	float light_length = sqrt(dotProduct(source, source));
	Vector normal_normalized(normal_vector.x / normal_length, normal_vector.y / normal_length, normal_vector.z / normal_length);
	Vector light_normalized(source.x / light_length, source.y / light_length, source.z / light_length);
	float lambert = max(0.f, dotProduct(normal_normalized, light_normalized));
	Vector h_vector(view_vector.x + source.x, view_vector.y + source.y, view_vector.z + source.z);
	float h_vector_length = sqrt(dotProduct(h_vector, h_vector));
	Vector h_normalized(h_vector.x / h_vector_length, h_vector.y / h_vector_length, h_vector.z / h_vector_length);
	float phong = pow(max(0.f, dotProduct(normal_normalized, h_normalized)),phong_exponent);
	return ((ambient1 * ambient2) + (diffuseCoefficient * intensity * lambert) + (specular_coefficient * phong));;
}

float RayTracer::triangleShading(float diffuseCoefficient, float intensity, float t, float specular_coefficient, float phong_exponent, float ambient1, float ambient2, Vector source, Ray ray, Triangle tetrahedron)
{
	Vector a_c(tetrahedron.a.x - tetrahedron.c.x, tetrahedron.a.y - tetrahedron.c.y, tetrahedron.a.z - tetrahedron.c.z);
	Vector b_c(tetrahedron.b.x - tetrahedron.c.x, tetrahedron.b.y - tetrahedron.c.y, tetrahedron.b.z - tetrahedron.c.z);
	Vector normal_vector(crossProduct(b_c, a_c).x / sqrt(dotProduct(b_c, a_c)), crossProduct(b_c, a_c).y / sqrt(dotProduct(b_c, a_c)), crossProduct(b_c, a_c).z / sqrt(dotProduct(b_c, a_c)));
	Vector view_vector(ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, ray.origin.z + t * ray.direction.z);
	float light_length = sqrt(dotProduct(source, source));
	Vector light_normalized(source.x / light_length, source.y / light_length, source.z / light_length);
	float lambert = max(0.f, dotProduct(normal_vector, light_normalized));
	Vector h_vector(view_vector.x + source.x, view_vector.y + source.y, view_vector.z + source.z);
	float h_vector_length = sqrt(dotProduct(h_vector, h_vector));
	Vector h_normalized(h_vector.x / h_vector_length, h_vector.y / h_vector_length, h_vector.z / h_vector_length);
	float phong = pow(max(0.f, dotProduct(normal_vector, h_normalized)), phong_exponent);
	return ((ambient1 * ambient2) + (diffuseCoefficient * intensity * lambert));
}


float RayTracer::triangleIntersection(Triangle tetrahedron, Ray ray)
{
	Vector a_b(tetrahedron.a.x - tetrahedron.b.x, tetrahedron.a.y - tetrahedron.b.y, tetrahedron.a.z - tetrahedron.b.z);
	Vector a_c(tetrahedron.a.x - tetrahedron.c.x, tetrahedron.a.y - tetrahedron.c.y, tetrahedron.a.z - tetrahedron.c.z);
	Vector b_c(tetrahedron.b.x - tetrahedron.c.x, tetrahedron.b.y - tetrahedron.c.y, tetrahedron.b.z - tetrahedron.c.z);
	Vector c_a(tetrahedron.c.x - tetrahedron.a.x, tetrahedron.c.y - tetrahedron.a.y, tetrahedron.c.z - tetrahedron.a.z);

	Vector normal_vector(crossProduct(b_c, a_c).x / sqrt(dotProduct(b_c, a_c)), crossProduct(b_c, a_c).y / sqrt(dotProduct(b_c, a_c)), crossProduct(b_c, a_c).z / sqrt(dotProduct(b_c, a_c)));
	float d = dotProduct(normal_vector, tetrahedron.b);
	float t = (d - (dotProduct(ray.origin, normal_vector))) / dotProduct(ray.direction, normal_vector);
	if (dotProduct(normal_vector, ray.direction) == 0)
	{
		return -1;
	}
	Vector x(ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, ray.origin.z + t * ray.direction.z);
	Vector x_a(x.x - tetrahedron.a.x, x.y - tetrahedron.a.y, x.z - tetrahedron.a.z);
	Vector x_b(x.x - tetrahedron.b.x, x.y - tetrahedron.b.y, x.z - tetrahedron.b.z);
	Vector x_c(x.x - tetrahedron.c.x, x.y - tetrahedron.c.y, x.z - tetrahedron.c.z);
	if (dotProduct(crossProduct(a_b, x_b), normal_vector) >= 0 && dotProduct(crossProduct(b_c, x_c), normal_vector) >= 0 && dotProduct(crossProduct(c_a, x_a), normal_vector) >= 0)
	{
		return t;
	}
	else
	{
		return -1;
	}


}






vector<float> RayTracer::getPixels()
{
	return pixels;
}
