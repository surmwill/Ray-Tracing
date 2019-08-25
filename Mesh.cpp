// Spring 2019

#include <iostream>
#include <fstream>

#include <glm/ext.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/io.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "polyroots.hpp"
#include <limits>
#include <cstdlib>
#include "Primitive.hpp"

// #include "cs488-framework/ObjFileDecoder.hpp"
#include "Mesh.hpp"

using namespace std;
using namespace glm;

Mesh::Mesh( const std::string& fname )
	: m_vertices()
	, m_faces()
{
	std::string code;
	double vx, vy, vz;
	size_t s1, s2, s3;

	std::ifstream ifs( fname.c_str() );
	while( ifs >> code ) {
		if( code == "v" ) {
			ifs >> vx >> vy >> vz;
			vec3 vert( vx, vy, vz );
			if( !init ) {
				min = vert;
				max = vert;
				init = true;
			}
			else {
				if( vert.x < min.x ) min.x = vert.x;
				if( vert.y < min.y ) min.y = vert.y;
				if( vert.z < min.z ) min.z = vert.z;
				if( vert.x > max.x ) max.x = vert.x;
				if( vert.y > max.y ) max.y = vert.y;
				if( vert.z > max.z ) max.z = vert.z;
			}
			m_vertices.push_back( vert );
		}
		else if( code == "f" ) {
			ifs >> s1 >> s2 >> s3;
			m_faces.push_back( Triangle( s1 - 1, s2 - 1, s3 - 1 ) );
		}
	}
	m_pos = min;
	m_size = max.x - min.x;
	if( max.y - min.y > m_size ) m_size = max.y - min.y;
	if( max.z - min.z > m_size ) m_size = max.z - min.z;
}

bool Mesh::getIntersection( const vec3& eye, const vec3& ray, double* t, vec3& norm, const mat4& trans ) const {
	// cout << m_vertices.size() << endl;
	#ifdef RENDER_BOUNDING_VOLUMES 
		return cubeGetIntersection( eye, ray, t, norm, trans );
	#endif

	if( !cubeGetIntersection( eye, ray, t, norm, trans ) ) return false;

	*t = -1;
	for( const Triangle& face : m_faces ) {
		const vec3 p0 = m_vertices[ face.v1 ];
		const vec3 p1 = m_vertices[ face.v2 ];
		const vec3 p2 = m_vertices[ face.v3 ];
		
		vec3 edge1, edge2, h, s, q;
		float a, f, u, v;
		edge1 = p1 - p0;
		edge2 = p2 - p0;
		h = cross( ray, edge2 );
		a = dot( edge1, h );
		if( abs( a ) < epsilon ) continue;
		f = 1.0f / a;
		s = eye - p0;
		u = f * dot( s, h );
		if( u < 0.0f - epsilon || u > 1.0f + epsilon ) continue;
		q = cross( s, edge1 );
		v = f * dot( ray, q );
		if( v < 0.0f - epsilon || u + v > 1.0f + epsilon ) continue;
		float t1 = f * dot( edge2, q );
		if( t1 > epsilon  && ( *t < 0 || t1 < ( *t - epsilon ) ) ) {
			norm = normalize( cross( ( p1 - p0 ), ( p2 - p0 ) ) );
			*t = t1;
		}
	}

	if( *t < 0 || isnan( *t ) ) return false;
	return true;
}

bool Mesh::cubeGetIntersection( const vec3& eye, const vec3& ray, double* t, vec3& norm, const mat4& trans ) const {
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

std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  out << "mesh {";
  /*
  
  for( size_t idx = 0; idx < mesh.m_verts.size(); ++idx ) {
  	const MeshVertex& v = mesh.m_verts[idx];
  	out << glm::to_string( v.m_position );
	if( mesh.m_have_norm ) {
  	  out << " / " << glm::to_string( v.m_normal );
	}
	if( mesh.m_have_uv ) {
  	  out << " / " << glm::to_string( v.m_uv );
	}
  }

*/
  out << "}";
  return out;
}
