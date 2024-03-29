// Spring 2019

#include "PhongMaterial.hpp"

PhongMaterial::PhongMaterial(
	const glm::vec3& kd, const glm::vec3& ks, double shininess, double reflection )
	: m_kd(kd)
	, m_ks(ks)
	, m_shininess(shininess)
	, reflection( reflection )
{}

PhongMaterial::~PhongMaterial()
{}
