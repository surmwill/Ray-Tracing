// Spring 2019

#include "Primitive.hpp"

#include <glm/ext.hpp>
#include <glm/gtx/io.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "polyroots.hpp"
#include <iostream>
#include <limits>
#include <cstdlib>

using namespace glm;
using namespace std;

Primitive::~Primitive()
{
}

Primitive::Type Primitive::getType() const { return Type::None; }

bool Primitive::getIntersection( const vec3& eye, const vec3& ray, double* t, vec3& norm, const mat4& trans ) const {
	// cout << "hello" << endl;
	return false;
}

// Sphere
Sphere::~Sphere()
{
}
  
Primitive::Type Sphere::getType() const { return Type::Sphere; }

bool Sphere::getIntersection( const vec3& eye, const vec3& ray, double* t, vec3& norm, const mat4& trans ) const {
	vec3 m_pos = vec3( 0.0f, 0.0f, 0.0f );
	double m_radius = 1.0f;

	const vec3& v = eye - m_pos;
	double roots[2];

	const int num_roots = quadraticRoots( dot( ray, ray ), 2 * dot( ray, v ), dot( v, v ) - pow( m_radius, 2 ), roots ) ;
	if( num_roots == 0 ) return false;
	else if( num_roots == 1 ) *t = roots[0];
	else if( num_roots == 2 ) {
		if( roots[0] <= 0 ) *t = roots[1];
		else if( roots[1] <= 0 ) *t = roots[0];
		else *t = std::min( roots[0], roots[1] );
	}

	if( *t <= 0 || isnan( *t ) ) return false;
	else if( abs( *t ) < epsilon ) return getIntersection( eye + ( ray * ( *t + epsilon ) ), ray, t, norm );  

	norm = normalize( ( eye + ( ray * *t ) ) - m_pos ); 
	return true;
}

// Cube
Cube::~Cube()
{
}

Primitive::Type Cube::getType() const { return Type::Cube; }

bool Cube::getIntersection( const vec3& eye, const vec3& ray, double* t, vec3& norm, const mat4& trans ) const {
	vec3 m_pos = vec3( 0.0f, 0.0f, 0.0f );
	double m_size = 1.0f;
	
	const vec3 min_box = m_pos;  // diagonally opposite corners; (x, y, z)
    const vec3 max_box = min_box + vec3( m_size ); // (x+r, y+r, z+r)

    // dirty, but time -- epsilon check 
    if( ray.x == 0 && ( eye.x >= max_box.x || eye.x <= min_box.x ) ) return false;
    if( ray.y == 0 && ( eye.y >= max_box.y || eye.y <= min_box.y ) ) return false;
    if( ray.z == 0 && ( eye.z >= max_box.z || eye.z <= min_box.z ) ) return false;

    float tx_min = ( min_box.x - eye.x ) / ray.x;
    float tx_max = ( max_box.x - eye.x ) / ray.x;

    float t_min = tx_min;
    float t_max = tx_max;
    if( t_min > t_max ) swap( t_min, t_max );

    float ty_min = ( min_box.y - eye.y ) / ray.y;
    float ty_max = ( max_box.y - eye.y ) / ray.y;
    if( ty_min > ty_max ) swap( ty_min, ty_max );

    if( ( t_min > ty_max ) || ( ty_min > t_max ) ) return false;

    if( ty_min > t_min ) t_min = ty_min;
    if( ty_max < t_max ) t_max = ty_max;

    float tz_min = ( min_box.z - eye.z ) / ray.z;
    float tz_max = ( max_box.z - eye.z ) / ray.z;
    if( tz_min > tz_max ) swap( tz_min, tz_max );

    if( ( t_min > tz_max ) || ( tz_min > t_max ) ) return false;

    if( tz_min > t_min ) t_min = tz_min;
    if( tz_max < t_max ) t_max = tz_max;

    if( t_min <= 0 ) *t = t_max;
    else if( t_max <= 0 ) *t = t_min;
    else *t = std::min( t_min, t_max );

    if( *t <= 0 || isnan( *t ) ) return false;
    else if( abs( *t ) < epsilon ) return getIntersection( eye + ( ray * ( *t + epsilon ) ), ray, t, norm );

    if( *t == tx_max || *t == tx_min ) norm = ray.x > 0 ? vec3( -1.0f, 0.0f, 0.0f ) : vec3( 1.0f, 0.0f, 0.0f );
    else if( *t == ty_max || *t == ty_min ) norm = ray.y > 0 ? vec3( 0.0f, -1.0f, 0.0f ) : vec3( 0.0f, 1.0f, 0.0f );
    else if( *t == tz_max || *t == tz_min ) norm = ray.z > 0 ? vec3( 0.0f, 0.0f, -1.0f ) : vec3( 0.0f, 0.0f, 1.0f );

	return true;
}

// NonhierSphere
NonhierSphere::~NonhierSphere()
{
}

Primitive::Type NonhierSphere::getType() const { return Type::NonhierSphere; }

bool NonhierSphere::getIntersection( const vec3& eye, const vec3& ray, double* t, vec3& norm, const mat4& trans ) const {
	const vec3& v = eye - m_pos;
	double roots[2];

	const int num_roots = quadraticRoots( dot( ray, ray ), 2 * dot( ray, v ), dot( v, v ) - pow( m_radius, 2 ), roots ) ;
	if( num_roots == 0 ) return false;
	else if( num_roots == 1 ) *t = roots[0];
	else if( num_roots == 2 ) {
		if( roots[0] <= 0 ) *t = roots[1];
		else if( roots[1] <= 0 ) *t = roots[0];
		else *t = std::min( roots[0], roots[1] );
	}

	if( *t <= 0 || isnan( *t ) ) return false;
	else if( abs( *t ) < epsilon ) return getIntersection( eye + ( ray * ( *t + epsilon ) ), ray, t, norm );  

	norm = normalize( ( eye + ( ray * *t ) ) - m_pos ); 
	return true;
}

const vec3& NonhierSphere::getPos() const { return m_pos; }

double NonhierSphere::getRadius() const { return m_radius; }

// NonhierBox
NonhierBox::~NonhierBox()
{
}

Primitive::Type NonhierBox::getType() const { return Type::NonhierBox; }

bool NonhierBox::getIntersection( const vec3& eye, const vec3& ray, double* t, vec3& norm, const mat4& trans ) const {
	const vec3 min_box = m_pos;  // diagonally opposite corners; (x, y, z)
	const vec3 max_box = min_box + vec3( m_size ); // (x+r, y+r, z+r)

	// dirty, but time -- epsilon check 
	if( ray.x == 0 && ( eye.x >= max_box.x || eye.x <= min_box.x ) ) return false;
	if( ray.y == 0 && ( eye.y >= max_box.y || eye.y <= min_box.y ) ) return false;
	if( ray.z == 0 && ( eye.z >= max_box.z || eye.z <= min_box.z ) ) return false;

	float tx_min = ( min_box.x - eye.x ) / ray.x;
	float tx_max = ( max_box.x - eye.x ) / ray.x;

	float t_min = tx_min;
	float t_max = tx_max;
	if( t_min > t_max ) swap( t_min, t_max );

	float ty_min = ( min_box.y - eye.y ) / ray.y;
	float ty_max = ( max_box.y - eye.y ) / ray.y;
	if( ty_min > ty_max ) swap( ty_min, ty_max );

	if( ( t_min > ty_max ) || ( ty_min > t_max ) ) return false;

	if( ty_min > t_min ) t_min = ty_min;
	if( ty_max < t_max ) t_max = ty_max;

	float tz_min = ( min_box.z - eye.z ) / ray.z;
	float tz_max = ( max_box.z - eye.z ) / ray.z;
	if( tz_min > tz_max ) swap( tz_min, tz_max );

	if( ( t_min > tz_max ) || ( tz_min > t_max ) ) return false;

	if( tz_min > t_min ) t_min = tz_min;
	if( tz_max < t_max ) t_max = tz_max;

	if( t_min <= 0 ) *t = t_max;
	else if( t_max <= 0 ) *t = t_min;
	else *t = std::min( t_min, t_max );

	if( *t <= 0 || isnan( *t ) ) return false; 
	else if( abs( *t ) < epsilon ) return getIntersection( eye + ( ray * ( *t + epsilon ) ), ray, t, norm );  

	if( *t == tx_max || *t == tx_min ) norm = ray.x > 0 ? vec3( -1.0f, 0.0f, 0.0f ) : vec3( 1.0f, 0.0f, 0.0f );
	else if( *t == ty_max || *t == ty_min ) norm = ray.y > 0 ? vec3( 0.0f, -1.0f, 0.0f ) : vec3( 0.0f, 1.0f, 0.0f );
	else if( *t == tz_max || *t == tz_min ) norm = ray.z > 0 ? vec3( 0.0f, 0.0f, -1.0f ) : vec3( 0.0f, 0.0f, 1.0f );
	/*
	if( *t == tx_max ) norm = vec3( 1.0f, 0.0f, 0.0f );
	else if( *t == tx_min ) norm = vec3( 1.0f, 0.0f, 0.0f );
    else if( *t == ty_max ) norm = vec3( 0.0f, 1.0f, 0.0f );
	else if( *t == ty_min ) norm = vec3( 0.0f, 1.0f, 0.0f );
	else if( *t == tz_max ) norm = vec3( 0.0f, 0.0f, 1.0f );
	else if( *t == tz_min ) norm = vec3( 0.0f, 0.0f, 1.0f );
	*/

	// cout << "point: " << eye + ( ray * *t ) << endl; // " norm: " << norm <<  endl;/
	return true;
}

const vec3& NonhierBox::getPos() const { return m_pos; }

double NonhierBox::getSize() const { return m_size; }
