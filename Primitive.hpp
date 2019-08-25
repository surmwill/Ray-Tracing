// Spring 2019

#pragma once

#include <glm/glm.hpp>

class Primitive {

public:
  enum Type{ None, Sphere, Cube, NonhierSphere, NonhierBox };
  virtual ~Primitive();
  virtual Type getType() const;
  virtual bool getIntersection( const glm::vec3& eye, const glm::vec3& ray, double* t, glm::vec3& norm, const glm::mat4& trans = glm::mat4{ 1.0f } ) const;

protected:
	const double epsilon = 0.001f;
};

class Sphere : public Primitive {
public:
  virtual ~Sphere();
  Type getType() const;
  bool getIntersection( const glm::vec3& eye, const glm::vec3& ray, double* t, glm::vec3& norm, const glm::mat4& trans = glm::mat4{ 1.0f } ) const override;
};

class Cube : public Primitive {
public:
  virtual ~Cube();
  Type getType() const;
  bool getIntersection( const glm::vec3& eye, const glm::vec3& ray, double* t, glm::vec3& norm, const glm::mat4& trans = glm::mat4{ 1.0f } ) const override;
};

class NonhierSphere : public Primitive {
public:
  NonhierSphere(const glm::vec3& pos, double radius)
    : m_pos(pos), m_radius(radius)
  {
  }
  virtual ~NonhierSphere();
  Type getType() const;
  bool getIntersection( const glm::vec3& eye, const glm::vec3& ray, double* t, glm::vec3& norm, const glm::mat4& trans = glm::mat4{ 1.0f } ) const override;
  const glm::vec3& getPos() const;
  double getRadius() const;

private:
  glm::vec3 m_pos;
  double m_radius;
};

class NonhierBox : public Primitive {
public:
  NonhierBox(const glm::vec3& pos, double size)
    : m_pos(pos), m_size(size)
  {
  } 
  virtual ~NonhierBox();
  Type getType() const; 
  bool getIntersection( const glm::vec3& eye, const glm::vec3& ray, double* t, glm::vec3& norm, const glm::mat4& trans = glm::mat4{ 1.0f } ) const override;
  const glm::vec3& getPos() const;
  double getSize() const;

private:
  glm::vec3 m_pos;
  double m_size;
};;
