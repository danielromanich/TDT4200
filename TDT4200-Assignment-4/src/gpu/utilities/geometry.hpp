#pragma once

#include <string>
#include <vector>
#include "mat4.hpp"
#include "cuda_runtime.h"

class Face;
class Material;
class Mesh;

class globalLight {
public:
	float3 direction;
	float3 colour;
	globalLight(float3 const vdirection, float3 const vcolour) : direction(vdirection), colour(vcolour) {}
};

class Material {
	public:
		static Material None;
		std::string name;
		float Ns;
		float3 Ka;
		float3 Kd;
		float3 Ks;
		float3 Ke;
		float Ni;
		float d;
		unsigned int illum;

		Material() : name("None"), Ns(100.0f), Ka(make_float3(1.0f, 1.0f, 1.0f)),
			Kd(make_float3(0.5f, 0.5f, 0.5f)), Ks(make_float3(0.5f, 0.5f, 0.5f)), Ke(make_float3(0.0f, 0.0f, 0.0f)),
			Ni(1.0), d(1.0), illum(2) {}

		Material(std::string vname) : name(vname), Ns(100.0f), Ka(make_float3(1.0f, 1.0f, 1.0f)),
			Kd(make_float3(0.5f, 0.5f, 0.5f)), Ks(make_float3(0.5f, 0.5f, 0.5f)), Ke(make_float3(0.0f, 0.0f, 0.0f)),
			Ni(1.0), d(1.0), illum(2) {}
};

float3 computeInterpolatedNormal(
        float3 const &normal0, float3 const &normal1, float3 const &normal2,
        float3 const &weights);

bool inPointInTriangle(float4 const &v0, float4 const &v1, float4 const &v2,
                       unsigned int const x, unsigned int const y,
                       float &u, float &v, float &w);

float computeDepth(
        float4 const &v0, float4 const &v1, float4 const &v2,
        float3 const &weights);


class GPUMesh {
public:
	float4* vertices;
	float3* normals;

	unsigned long vertexCount = 0;

	float3 objectDiffuseColour;

	bool hasNormals;

	GPUMesh clone();
};
