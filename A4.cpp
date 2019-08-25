// Spring 2019

#include <glm/ext.hpp>
#include <glm/gtx/io.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "A4.hpp"
#include "GeometryNode.hpp"
#include "PhongMaterial.hpp"
#include "polyroots.hpp"
#include "Mesh.hpp"
#include "scene_lua.hpp"
#include "Primitive.hpp"

#include <cfloat>
#include <algorithm>

using namespace std;
using namespace glm;

void A4_Render(
		// What to render  
		SceneNode * root,

		// Image to write to, set to a given width and height  
		Image & image,

		// Viewing parameters  
		const glm::vec3 & eye,
		const glm::vec3 & view,
		const glm::vec3 & up,
		double fovy,

		// Lighting parameters  
		const glm::vec3 & ambient,
		const std::list<Light *> & lights
) {
	std::cout << "Calling A4_Render(\n" <<
		  "\t" << *root <<
          "\t" << "Image(width:" << image.width() << ", height:" << image.height() << ")\n"
          "\t" << "eye:  " << glm::to_string(eye) << std::endl <<
		  "\t" << "view: " << glm::to_string(view) << std::endl <<
		  "\t" << "up:   " << glm::to_string(up) << std::endl <<
		  "\t" << "fovy: " << fovy << std::endl <<
          "\t" << "ambient: " << glm::to_string(ambient) << std::endl <<
		  "\t" << "lights{" << std::endl;

	for(const Light * light : lights) {
		std::cout << "\t\t" <<  *light << std::endl;
	}
	std::cout << "\t}" << std::endl;
	std:: cout <<")" << std::endl;

	//for( auto& pair : mesh_map ) {
	//	cout << pair.first << " " << pair.second << endl;
	//}

	size_t h = image.height();
	size_t w = image.width();

	// plane parameters
	const float d = 1.0f;     // assume d is 1
    const float plane_h = 2.0f * d * glm::tan( glm::radians( fovy / 2 ) );
    const float plane_w = (float)w / (float)h * plane_h;

    mat4 trans1{ 1.0f };
    trans1[3][0] = -(float)w / 2.0f;
    trans1[3][1] = -(float)h / 2.0f;
    trans1[3][2] = d;
    //std::cout << "T1" << trans1 << std::endl;

    mat4 scale2{ 1.0f };
    scale2[0][0] = -plane_h / (float)h;
    scale2[1][1] = plane_h / (float)h;
    scale2[2][2] = 1;
    //std::cout << "S2" << scale2 << std::endl;

    mat4 rot3{ 1.0f };
    // vec3 r = -( view - eye );   // corresponds to w in the notes ( should be view - eye but... )
    vec3 r = view;   // lookAt - lookFrom (eye towards to look at )
    r = glm::normalize( r );
	// cout << "r: " << r << endl;
    vec3 u = glm::cross( up, r );
    u = glm::normalize( u );
    vec3 v = -glm::cross( r, u );
    rot3[0][0] = u.x;
    rot3[0][1] = u.y;
    rot3[0][2] = u.z;
    rot3[1][0] = v.x;
    rot3[1][1] = v.y;
    rot3[1][2] = v.z;
    rot3[2][0] = r.x;
    rot3[2][1] = r.y;
    rot3[2][2] = r.z;
    //cout << "R3" << rot3 << endl;

    mat4 trans4{ 1.0f };
    trans4[3][0] = eye.x;
    trans4[3][1] = eye.y;
    trans4[3][2] = eye.z;
    //cout << "T4" << trans4 << endl;

    const mat4 pixel_to_world = trans4 * rot3 * scale2 * trans1;
	const int center_x = w / 2;
	const int center_y = h / 2;
	const double radius = ( (double)center_y / 2.5f );

	const vec3 dark_blue( 0.0f, 0.0f, 0.5f );
	const vec3 light_blue( 0.8f, 0.9f, 1.0f );
	const vec3 white( 1.0f, 1.0f, 1.0f );
	const vec3 light_brown( 0.7f, 0.6f, 0.4f );
	const vec3 dark_brown( 0.3f, 0.2f, 0.0f );
	const vec3 orange( 0.9f, 0.5f, 0.0f );
	const vec3 yellow( 1.0f, 0.9f, 0.2f );

	for (uint y = 0; y < h; ++y) {
		for (uint x = 0; x < w; ++x) {
			const vec4 p_world = pixel_to_world * vec4{ x, y, 0.0f, 1.0f };
			vec3 ray = vec3( p_world ) - eye;
			vec3 eye2 = eye;

			vec3 col = getCol( eye2, ray, root, lights, ambient );
			if( col.x == -1.0f && col.y == -1.0f && col.z == -1.0f ) {
				col = vec3( 1.0f );
				double t;
				if( y <= center_y ) {
					double dist_x = abs( (double)x - (double)center_x );
					double dist_y = abs( (double)y - (double)center_y );
					double dist_center = sqrt( pow( dist_x , 2 ) + pow( dist_y, 2 ) );
					
					if( dist_center <= radius ) {
						t = dist_center / radius;
						col = t * orange + ( 1 - t ) * yellow;
					}
					else {
						t = (double)y /  center_y;
						col = t * light_blue + ( 1 - t ) * dark_blue;
					}
				}
				else {
					t = (double)( y - center_y ) / (double)( h - center_y );
					col = t * dark_brown + ( 1 - t ) * light_brown;
				}
			}

			// Red: 
			image(x, y, 0) = col.x;
			// Green: 
			image(x, y, 1) = col.y;
			// Blue: 
			image(x, y, 2) = col.z;

			/*
			double t = DBL_MAX;
			Material* mat = nullptr;
			Material** p_mat = &mat;
			vec3 norm{ 1.0f };
		
			// collision, have point and material
			if( checkCollisions( eye, ray, root, &t, norm, p_mat ) ) {
				vec3 col = getCol( eye, ray, norm, *( static_cast<PhongMaterial*>( *p_mat ) ) );
				// cout << norm << endl;
				
				// Red: 
				image(x, y, 0) = col.x;
				// Green: 
				image(x, y, 1) = col.y;
				// Blue: 
				image(x, y, 2) = col.z;
			}
			else {
				// Red: 
				image(x, y, 0) = (double)1.0;
				// Green: 
				image(x, y, 1) = (double)1.0;
				// Blue: 
				image(x, y, 2) = (double)1.0;
			}
			*/
		}
	}
}

// if we return true t will contain the distance to the nearest object hit and p_mat will contain the material of what was hit
bool checkCollisions( vec3 point, vec3 ray, SceneNode* node, double* t, vec3& norm, vec3& hit_point, Material** p_mat, mat4 trans ) {
	bool collision = false;

	//vec3 p_world = -(ray - point);
	
	//cout << "ray before: " << vec4( ray, 1.0f ) << endl;
	ray = vec3( node->get_inverse() * vec4( ray, 0.0f ) );
	point = vec3( node->get_inverse() * vec4( point, 1.0f ) );
	trans *= node->get_transform();
	//trans = node->get_inverse() * trans;
	//ray = point - p_world;
	// cout << "ray after: " << ray << endl;
	
	// cout << node->get_inverse() << endl;

	// first test collision against us, the current node
	if( node->m_nodeType == NodeType::GeometryNode ) {
		GeometryNode* geo_node = static_cast<GeometryNode*>( node );
		double test_t;
		vec3 test_norm;

		if( geo_node->m_primitive->getIntersection( point, ray, &test_t, test_norm ) && test_t < *t ) {
			// cout << "intersection: " << ray << endl;
			*t = test_t;
			norm = test_norm;
			// norm = normalize( transpose( mat3( node->get_inverse() ) )  * norm );
			norm = normalize( transpose( inverse( mat3( trans ) ) ) * norm ); 
			*p_mat = geo_node->m_material;
			hit_point = point + ( ray * *t );
			collision = true;
			hit_point = vec3( node->get_transform() * vec4( hit_point, 1.0f ) );
		} 
	}

	// then test collision against all children
	for( SceneNode* child : node->children ) { 
		if( checkCollisions( point, ray, child, t, norm, hit_point, p_mat, trans ) ) {
			collision = true;
			hit_point = vec3( node->get_transform() * vec4( hit_point, 1.0f ) );
		}
	 }

	// hit_point = vec3( node->get_transform() * vec4( hit_point, 1.0f ) );

	// _point 
	//point = vec3( node->get_transform() * vec4( point, 1.0f ) );
	//ray = vec3( node->get_transform() * vec4( ray, 0.0f ) );

	return collision;
}

// returns true if the shadow ray hit an object
bool shadowRayHit( vec3 point, vec3 shadow_ray, SceneNode* node ) {
	double t_placeholder;
	vec3 n_placeholder;

	shadow_ray = vec3( node->get_inverse() * vec4( shadow_ray, 0.0f ) );
	point = vec3( node->get_inverse() * vec4( point, 1.0f ) );

	if( node->m_nodeType == NodeType::GeometryNode ) {
		GeometryNode* geo_node = static_cast<GeometryNode*>( node );
		if( geo_node->m_primitive->getIntersection( point, shadow_ray, &t_placeholder, n_placeholder ) ) {
			// cout << "hit: " << geo_node->m_name << " " << t_placeholder << " " << endl;
			return true;
		}
	}

	for( SceneNode* child : node->children ) {
		if( shadowRayHit( point, shadow_ray, child ) ) return true;	
	}

	shadow_ray = vec3( node->get_transform() * vec4( shadow_ray, 0.0f ) );

	return false;
}

vec3 getCol( vec3& point, vec3& ray, SceneNode* root, const list<Light*>& lights, const vec3& ambient, int hits_left )
{
	double t = DBL_MAX;
	Material* mat = nullptr;
	Material** p_mat = &mat;
	vec3 norm{ 1.0f };
	vec3 hit_point{ 1.0f };
	vec3 col;

	if( checkCollisions( point, ray, root, &t, norm, hit_point, p_mat, mat4{ 1.0f } ) ) {
		// col = vec3( 0.0f, 1.0f, 0.0f );
		PhongMaterial* phong_mat = static_cast<PhongMaterial*>( mat );
		col = ambient * phong_mat->m_kd;
	
		for( const Light* light : lights ) {
			vec3 shadow_ray = light->position - hit_point;
			vec3 refl_ray = reflect( ray, norm );
			if( !shadowRayHit( hit_point, shadow_ray, root ) ) {
				const double r = length( shadow_ray );
				const double* falloff = light->falloff;
				const double atten = 1.0f / ( falloff[0] + r * falloff[1] + pow( r, 2 ) * falloff[2] );

				const vec3 diffuse = dot( normalize( shadow_ray ), norm ) * phong_mat->m_kd;
				const vec3 spec = pow( dot( normalize( reflect( -shadow_ray, norm ) ), normalize( -ray ) ), phong_mat->m_shininess ) * phong_mat->m_ks;  
				col += atten * ( ( diffuse + spec ) * light->colour );	
			} 
			if( phong_mat->reflection > 0 && hits_left > 0) {
				col += phong_mat->reflection * getCol( hit_point, refl_ray, root, lights, ambient, hits_left - 1 );
			}
		}
	}
	// the first ray we sent didn't hit anything, set the pixel to the background
	else if( hits_left == MAX_HITS ) {
		// return vec3( 2.0f );
		col = vec3{ -1.0f, -1.0f, -1.0f };
		//col = vec3( 1.0f );
	}
	// we prev sent at least one ray which hit an object, and this ray can get to the light source; calc colour
	else col = vec3( 0.0f );
	return col;
}
