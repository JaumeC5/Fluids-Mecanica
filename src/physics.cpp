//////////////////////
//// LIFE IS GOOD ////
//Toni Ferrari Juan///
//////////////////////

#include <imgui\imgui.h>
#include <imgui\imgui_impl_glfw_gl3.h>
#include <iostream>
#include <glm\gtc\matrix_transform.hpp>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <../include/GLFW/glfw3.h>
using namespace std;
bool show_test_window = false;

// Helps to reset the system
bool reset = false;
float A; // gui
// Collision detectors
bool collisionLeft;
bool collisionRight;
bool collisionUp;
bool collisionDown;
bool collisionFront;
bool collisionBack;
int Rtime;
float distanceSprings;
float ellogation;

// Positions
glm::vec3 *initialMesh;
glm::vec3* actualMesh;
glm::vec3 * tempMesh;
glm::vec3 * finalMesh;

// External Force
glm::vec3 f;

// Sphere things
glm::vec3 centreSphere;
float RandomRadiusSphere;
float sphereX;
float sphereY;
float sphereZ;

time_t theTime = time(0);

// Constants
float Ke = 900;
float Kd = 50;

// Collision detector

bool hasCollision(glm::vec3 Pt, glm::vec3 n, float d, glm::vec3 PtPost, int plane) {

	float getPos;
	getPos = ((glm::dot(n, Pt) + d) * (glm::dot(n, PtPost) + d));
	if (getPos <= 0) { return true; }
	else { return false; }
}

// Bounce

void collidePlane(glm::vec3 n, float d, glm::vec3 &pos, glm::vec3 &prevPos) {

	pos = pos - (1 + 0.8f) * (glm::dot(pos, n) + d) * n;
	prevPos = prevPos - (1 + 0.8f) * (glm::dot(prevPos, n) + d) * n;
}

// Getting alpha - obsolet ( 1st method used to calculate contact point with the sphere )

float getTheAlpha(float a, float b, float c) {
	float res1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
	float res2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);

	if (res1 <= 1 && res1 >= 0) return res1;
	if (res2 <= 1 && res2 >= 0) return res2;

}

// Now that's good. Getting the point of contact between particle and the Sphere

glm::vec3 getQ(glm::vec3 position, glm::vec3 lastP, glm::vec3 cS, float r) {

	float a, b, c, alpha = 0;
	glm::vec3 dir = glm::normalize(position - lastP);
	a = glm::dot(dir, dir);
	b = 2 * glm::dot(dir, lastP - cS);
	c = glm::dot(lastP - cS, lastP - cS) - r*r;

	float res1 = (-b + sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);
	float res2 = (-b - sqrt(pow(b, 2) - 4 * a*c)) / (2 * a);

	if (res1 <= 1 && res1 >= 0) alpha = res1;
	if (res2 <= 1 && res2 >= 0) alpha = res2;

	glm::vec3 Q = lastP + (position - lastP)*alpha;
	return Q;
}

namespace Sphere {
	extern void setupSphere(glm::vec3 pos = glm::vec3(5.f, 1.f, 0.f), float radius = 1.f);
	extern void cleanupSphere();
	extern void updateSphere(glm::vec3 pos, float radius = 1.f);
	extern void drawSphere();
}
namespace ClothMesh {
	extern const int numCols;
	extern const int numRows;
	extern const int numVerts;
	extern void updateClothMesh(float* array_data);
}

void GUI() {
	{	//FrameRate
		ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
		if (ImGui::Button("Reset")) {
			reset = true;
		}
		
	}
	// ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
	if (show_test_window) {//
		ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
		ImGui::ShowTestWindow(&show_test_window);
	}
	// ImGui test window. Most of the sample code is in ImGui::ShowTestWindow()
	if (show_test_window) {//
		ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
		ImGui::ShowTestWindow(&show_test_window);
	}
}

class Particle {
public:

	glm::vec3 initialPos;
	glm::vec3 currPos;
	
	glm::vec3 Q;

	bool collided = true;
	bool calculated = false;
	float alpha;
	int mass = 1;
};

Particle* parVerts;

// Nice formula bro -- I know man, I know :v

void sphereCollision(glm::vec3 Q, glm::vec3 &actPos, glm::vec3 &lastPos, glm::vec3 cS) {

	glm::vec3 n = glm::normalize(Q - cS);
	float d = -glm::dot(n, Q);
	if (glm::length(actPos - centreSphere) <= RandomRadiusSphere) {
		collidePlane(n, d, actPos, lastPos);	// Execute bounce
	}
}

void PhysicsCleanup() {

	delete[] initialMesh;
	delete[] actualMesh;
	delete[] tempMesh;
	delete[] finalMesh;
	delete parVerts;
}

// Let's initiate all that stuff
float lambda; // --> modificar amb sa gui
float K;
float omega; // frequency GUI
glm::vec3 dir = glm::vec3(1.0f, 0.f, 0.0f);
void PhysicsInit() {

	reset = false;
	Rtime = 0;
	lambda = 2;
	K = (2 * 3.14159) / lambda;
	A = 2;
	omega = 6;

	initialMesh = new glm::vec3[14 * 18];
	actualMesh = new glm::vec3[14 * 18];
	parVerts = new Particle[14 * 18];

	for (int j = 0; j < ClothMesh::numRows; j++) {
		for (int i = 0; i < ClothMesh::numCols; i++) {

			initialMesh[(j * ClothMesh::numCols + i)] = glm::vec3(-4.5 + j*0.5f, 0.f, -4.5 + i*0.5f);
			actualMesh[(j * ClothMesh::numCols + i)] = initialMesh[(j * ClothMesh::numCols + i)];
		
			//The omnipresent gravity
			f = glm::vec3(0, -9.81, 0);

		
		}
	}
	sphereX = -4 + rand() % 8;
	sphereY = 0 + rand() % 9;
	sphereZ = -4 + rand() % 8;
	centreSphere = glm::vec3(sphereX, sphereY, sphereZ);
	RandomRadiusSphere = 0.5f + rand() % 3;
}


void resetAll()
{
	theTime++;
	if (theTime > 800) // 20 secs
	{
		PhysicsCleanup();
		PhysicsInit();
		Sphere::updateSphere(glm::vec3(sphereX, sphereY, sphereZ), RandomRadiusSphere);
		theTime = 0;
	}
	if (reset == true)
	{
		PhysicsCleanup();
		PhysicsInit();
		Sphere::updateSphere(glm::vec3(sphereX, sphereY, sphereZ), RandomRadiusSphere);
		theTime = 0;
	}
	centreSphere = glm::vec3(sphereX, sphereY, sphereZ);
}

void verletSprings(int i, int j, float dt)
{

		actualMesh[(j * ClothMesh::numCols + i)].x = initialMesh[(j * ClothMesh::numCols + i)].x - (1.f / 1.f) * A * sin(1.f* initialMesh[(j * ClothMesh::numCols + i)].x - 1.f);
		actualMesh[(j * ClothMesh::numCols + i)].y = A * cos(1 * initialMesh[(j * ClothMesh::numCols + i)].y - 1.f);



		collisionLeft = hasCollision(actualMesh[(j * ClothMesh::numCols + i)], glm::vec3(1, 0, 0), 5, initialMesh[(j * ClothMesh::numCols + i)], 1);
		//right
		collisionRight = hasCollision(actualMesh[(j * ClothMesh::numCols + i)], glm::vec3(-1, 0, 0), 5, initialMesh[(j * ClothMesh::numCols + i)], 2);
		//up
		collisionUp = hasCollision(actualMesh[(j * ClothMesh::numCols + i)], glm::vec3(0, -1, 0), 10, initialMesh[(j * ClothMesh::numCols + i)], 3);
		//down
		collisionDown = hasCollision(actualMesh[(j * ClothMesh::numCols + i)], glm::vec3(0, 1, 0), 0, initialMesh[(j * ClothMesh::numCols + i)], 4);
		//front
		collisionFront = hasCollision(actualMesh[(j * ClothMesh::numCols + i)], glm::vec3(0, 0, -1), 5, initialMesh[(j * ClothMesh::numCols + i)], 5);
		//back
		collisionBack = hasCollision(actualMesh[(j * ClothMesh::numCols + i)], glm::vec3(0, 0, 1), 5, initialMesh[(j * ClothMesh::numCols + i)], 6);

		if (collisionLeft) {
			collidePlane(glm::vec3(1, 0, 0), 5, initialMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)]);
		}
		if (collisionRight) {
			collidePlane(glm::vec3(-1, 1, 0), 5, initialMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)]);
		}
		if (collisionUp) {
			collidePlane(glm::vec3(0, -1, 0), 10, initialMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)]);
		}
		if (collisionDown) {
			collidePlane(glm::vec3(0, 1, 0), 0, initialMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)]);
		}
		if (collisionFront) {
			collidePlane(glm::vec3(0, 0, -1), 5, initialMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)]);
		}
		if (collisionBack) {
			collidePlane(glm::vec3(0, 0, 1), 5, initialMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)]);
		}

		parVerts[(j * ClothMesh::numCols + i)].Q = getQ(actualMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)], centreSphere, RandomRadiusSphere);
		sphereCollision(parVerts[(j * ClothMesh::numCols + i)].Q, initialMesh[(j * ClothMesh::numCols + i)], actualMesh[(j * ClothMesh::numCols + i)], centreSphere);
		parVerts[(j * ClothMesh::numCols + i)].currPos = actualMesh[(j * ClothMesh::numCols + i)];
		parVerts[(j * ClothMesh::numCols + i)].initialPos = initialMesh[(j * ClothMesh::numCols + i)];
}

void PhysicsUpdate(float dt) {
	
		for (int j = 0; j < ClothMesh::numRows; j++) {
			for (int i = 0; i < ClothMesh::numCols; i++) {
				actualMesh[j * ClothMesh::numCols + i] = initialMesh[j * ClothMesh::numCols + i] - (dir / K)*A*sin(glm::dot(dir, (initialMesh[j * ClothMesh::numCols + i]) - omega*dt));
				actualMesh[j * ClothMesh::numCols + i].y = A*cos(glm::dot(dir, actualMesh[j * ClothMesh::numCols + i]) - omega*dt);

				Rtime = glfwGetTime();
			}
		}
	
	resetAll();
	ClothMesh::updateClothMesh(&actualMesh[0].x);
}