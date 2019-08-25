// Spring 2019

#pragma once

#include <glm/glm.hpp>

#include "SceneNode.hpp"
#include "Light.hpp"
#include "Image.hpp"

class GeometryNode;
class Material;
class PhongMaterial;

const int MAX_HITS = 3;

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
);

bool checkCollisions( glm::vec3 point, glm::vec3 ray, SceneNode* node, double* t, glm::vec3& norm, glm::vec3& hit_point, Material** p_mat, glm::mat4 trans );

bool shadowRayHit( glm::vec3 point, glm::vec3 ray, SceneNode* node );

glm::vec3 getCol( glm::vec3& point, glm::vec3& ray, SceneNode* root, const std::list<Light*>& lights, const glm::vec3& ambient, int hits_left = MAX_HITS );
