// Spring 2019

#pragma once

#include <glm/glm.hpp>

#include "Material.hpp"

class PhongMaterial : public Material {
public:
  PhongMaterial(const glm::vec3& kd, const glm::vec3& ks, double shininess, double reflection = 0.0f );
  virtual ~PhongMaterial();

  glm::vec3 m_kd; // change this back
  glm::vec3 m_ks;
  double m_shininess;
  double reflection;

private:
//  glm::vec3 m_kd;
//  glm::vec3 m_ks;
//  double m_shininess;
};
