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
float A;
// Collision detectors
int Rtime;

// Positions
glm::vec3* initialMesh;
glm::vec3* actualMesh;
glm::vec3 *temp;
glm::vec3 *lastMesh;

//Distance between nodes and sphere

glm::vec3 *q;

// External Force
glm::vec3 gravity;
glm::vec3 buoyancy;
glm::vec3 drag;
glm::vec3 lift;

//Total force

glm::vec3 total;

// Water density
float density = 10.0;

// Coeficients
float dragC = 0.47;
float massa = 1.f;


// Sphere things
glm::vec3 centreSphere;
float RandomRadiusSphere;
float sphereX;
float sphereY;
float sphereZ;
float subArea;

time_t theTime = time(0);
float totalTime = 0;

float lambda;
float K;
float omega;
glm::vec3 dir = glm::vec3(1.0f, 0.f, 0.0f);
float height;
float vsub = 0;
float diff;
glm::vec3 sphereSpeed;

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
	if (ImGui::CollapsingHeader("Spring parameters")) {
		ImGui::SliderFloat("Amplitude", &A, 0.1, 1.5);
		ImGui::SliderFloat("Frequency", &omega, 1, 10);
		ImGui::SliderFloat("Direction x", &dir.x, 0, 1);
		ImGui::SliderFloat("Direction y", &dir.y, 0, 1);
		ImGui::SliderFloat("Direction z", &dir.z, 0, 1);
		ImGui::SliderFloat("Height", &height, 0.5, 10);
		ImGui::SliderFloat("Mass", &massa, 0.1, 10);
		ImGui::SliderFloat("Liquid's density", &density, 1, 30);

	}

}

void calculateDrag() {
	subArea = abs(2 * 3.14150*RandomRadiusSphere*diff);
	drag.y = -0.5*density*subArea * glm::length(sphereSpeed)*sphereSpeed.y * 0.2;

}

bool hasCollision(glm::vec3 Pt, glm::vec3 n, float d, glm::vec3 PtPost, int plane) {

	float getPos;
	getPos = ((glm::dot(n, Pt) + d) * (glm::dot(n, PtPost) + d));
	if (getPos <= 0) { return true; }
	else { return false; }
}

void collidePlane(glm::vec3 n, float d, glm::vec3 &pos, glm::vec3 &prevPos) {

	pos = pos - (1 + 0.8f) * (glm::dot(pos, n) + d) * n;
	prevPos = prevPos - (1 + 0.8f) * (glm::dot(prevPos, n) + d) * n;
}

float getTheAlpha(float a, float b, float c) {
	float res1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
	float res2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);

	if (res1 <= 1 && res1 >= 0) return res1;
	if (res2 <= 1 && res2 >= 0) return res2;

}

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

class Particle {
public:

	glm::vec3 initialPos;
	glm::vec3 currPos;

	bool collided = true;
	bool calculated = false;
	float alpha;
	int mass = 1;
};

Particle* parVerts;


void sphereCollision(glm::vec3 Q, glm::vec3 &actPos, glm::vec3 &lastPos, glm::vec3 cS) {

	glm::vec3 n = glm::normalize(Q - cS);
	float d = -glm::dot(n, Q);
	if (glm::length(actPos - centreSphere) <= RandomRadiusSphere) {
		//collidePlane(n, d, actPos, lastPos);
		total = buoyancy;
	}



}

void PhysicsCleanup() {

	delete[] initialMesh;
	delete[] actualMesh;
	delete parVerts;
}


void PhysicsInit() {

	sphereSpeed = glm::vec3(0.f, 0.f, 0.f);
	reset = false;
	Rtime = 0;
	lambda = 2;
	K = (2 * 3.14159) / lambda;
	A = 0.4;
	omega = 6;
	height = 0.5;
	density = 1.7;
	vsub = 0;
	subArea = 0;
	gravity = glm::vec3(0, -9.81, 0);
	drag = glm::vec3(0, 0, 0);
	lift = glm::vec3(0, 0, 0);
	total = gravity;
	diff = 0;
	initialMesh = new glm::vec3[14 * 18];
	temp = new glm::vec3[14 * 18];
	lastMesh = new glm::vec3[14 * 18];

	actualMesh = new glm::vec3[14 * 18];
	parVerts = new Particle[14 * 18];
	q = new glm::vec3[14 * 18];

	for (int j = 0; j < ClothMesh::numRows; j++) {
		for (int i = 0; i < ClothMesh::numCols; i++) {
			temp[(j * ClothMesh::numCols + i)] = glm::vec3(0.f, 0.f, 0.f);
			lastMesh[(j * ClothMesh::numCols + i)] = glm::vec3(0.f, 0.f, 0.f);

			initialMesh[(j * ClothMesh::numCols + i)] = glm::vec3(-4.5 + j*0.5f, 0.f, -4.5 + i*0.5f);
			actualMesh[(j * ClothMesh::numCols + i)] = initialMesh[(j * ClothMesh::numCols + i)];

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

void PhysicsUpdate(float dt) {
	totalTime += dt;

	buoyancy = density*vsub*9.81f*glm::vec3(0.f, 1.f, 0.f);
	total = (gravity + buoyancy + drag)*massa;
	sphereSpeed.y = sphereSpeed.y + dt * (total.y / massa);

	cout << diff << endl;

	sphereY = sphereY + dt*sphereSpeed.y;
	for (int j = 0; j < ClothMesh::numRows; j++) {
		for (int i = 0; i < ClothMesh::numCols; i++) {
			temp[j * ClothMesh::numCols + i] = actualMesh[j * ClothMesh::numCols + i];
			actualMesh[j * ClothMesh::numCols + i] = initialMesh[j * ClothMesh::numCols + i] - (dir / K)*A*sin(glm::dot(dir, (initialMesh[j * ClothMesh::numCols + i]) - omega*totalTime));
			actualMesh[j * ClothMesh::numCols + i].y = height + A*cos(glm::dot(dir, actualMesh[j * ClothMesh::numCols + i]) - omega*totalTime);
			lastMesh[j * ClothMesh::numCols + i] = temp[j * ClothMesh::numCols + i];
			q[j * ClothMesh::numCols + i] = getQ(actualMesh[(j * ClothMesh::numCols + i)], lastMesh[(j * ClothMesh::numCols + i)], centreSphere, RandomRadiusSphere);
			diff = abs(sphereY - RandomRadiusSphere - actualMesh[j * ClothMesh::numCols + i].y);


			if (sphereY - RandomRadiusSphere <= actualMesh[j * ClothMesh::numCols + i].y) {
				vsub = (RandomRadiusSphere * 2) *diff;
				calculateDrag();
			}
			else {
				vsub = 0;
			}

			Rtime = glfwGetTime();
		}
	}

	Sphere::updateSphere(glm::vec3(sphereX, sphereY, sphereZ), RandomRadiusSphere);
	resetAll();
	ClothMesh::updateClothMesh(&actualMesh[0].x);
}

