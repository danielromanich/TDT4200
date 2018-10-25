#include <iostream>
#include "geometry.hpp"

Material Material::None = Material();

inline float dot(float3 a, float3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float3 normalize(float3 v)
{
    float invLen = 1.0f / sqrtf(dot(v, v));
    v.x *= invLen;
    v.y *= invLen;
    v.z *= invLen;
    return v;
}

GPUMesh GPUMesh::clone() {
	GPUMesh clonedMesh;

	clonedMesh.vertices = new float4[this->vertexCount];
	std::copy(this->vertices, this->vertices + this->vertexCount, clonedMesh.vertices);

	clonedMesh.normals = new float3[this->vertexCount];
	std::copy(this->normals, this->normals + this->vertexCount, clonedMesh.normals);
	
	clonedMesh.objectDiffuseColour = this->objectDiffuseColour;

	clonedMesh.hasNormals = this->hasNormals;

	clonedMesh.vertexCount = this->vertexCount;

	return clonedMesh;
}


bool inPointInTriangle(
		float4 const &v0, float4 const &v1, float4 const &v2,
		unsigned int const x, unsigned int const y,
		float &u, float &v, float &w) {
		u = (((v1.y - v2.y) * (x    - v2.x)) + ((v2.x - v1.x) * (y    - v2.y))) /
				 	 (((v1.y - v2.y) * (v0.x - v2.x)) + ((v2.x - v1.x) * (v0.y - v2.y)));
		if (u < 0) {
			return false;
		}
		v = (((v2.y - v0.y) * (x    - v2.x)) + ((v0.x - v2.x) * (y    - v2.y))) /
					(((v1.y - v2.y) * (v0.x - v2.x)) + ((v2.x - v1.x) * (v0.y - v2.y)));
		if (v < 0) {
			return false;
		}
		w = 1 - u - v;
		if (w < 0) {
			return false;
		}
		return true;
}

float3 computeInterpolatedNormal(
		float3 const &normal0,
		float3 const &normal1,
		float3 const &normal2,
		float3 const &weights
	) {
	float3 weightedN0, weightedN1, weightedN2;

	weightedN0.x = (normal0.x * weights.x);
	weightedN0.y = (normal0.y * weights.x);
	weightedN0.z = (normal0.z * weights.x);

	weightedN1.x = (normal1.x * weights.y);
	weightedN1.y = (normal1.y * weights.y);
	weightedN1.z = (normal1.z * weights.y);

	weightedN2.x = (normal2.x * weights.z);
	weightedN2.y = (normal2.y * weights.z);
	weightedN2.z = (normal2.z * weights.z);

	float3 weightedNormal;

	weightedNormal.x = weightedN0.x + weightedN1.x + weightedN2.x;
	weightedNormal.y = weightedN0.y + weightedN1.y + weightedN2.y;
	weightedNormal.z = weightedN0.z + weightedN1.z + weightedN2.z;

	return normalize(weightedNormal);
}

float computeDepth(
		float4 const &v0, float4 const &v1, float4 const &v2,
		float3 const &weights) {
	return weights.x * v0.z + weights.y * v1.z + weights.z * v2.z;
}