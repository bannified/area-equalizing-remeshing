#include "mesh.h"

#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include "math.h"
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "mesh.h"
#include <map>
#include <queue>
#include <string>
#include <iomanip>
using namespace std;

#include "vec3.h"

void myObjType::draw() {

	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glPushMatrix();
	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	glScalef(4.0 / longestSide, 4.0 / longestSide, 4.0 / longestSide);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);
	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);
// uncomment the following after you computed the normals
		glNormal3dv(nlist[i]);    
		for (int j = 0; j < 3; j++)
			glVertex3dv(vlist[tlist[i][j]]);
		glEnd();
	
	}
	glDisable(GL_LIGHTING);

	glPopMatrix();
}

void myObjType::writeFile(char* filename)
{
	ofstream outFile;
	outFile.open(filename, ofstream::out);
	if (!outFile.is_open()) {
		cout << "File cannot be written to: " << filename << endl;
		exit(2);
	}

	// write vertices
	const string vertexPrefix = "v";
	for (int i = 1; i < vcount; i++) {
		double* vertex = vlist[i];
		outFile << vertexPrefix << " " << vertex[0] << " " << vertex[1] << " " << vertex[2] << endl;
	}

	const string facePrefix = "f";
	for (int i = 1; i < tcount; i++) {
		int* tri = tlist[i];
		outFile << facePrefix << " " << tri[0] << " " << tri[1] << " " << tri[2] << endl;
	}

	cout << "Writing to " << filename << " done.\n";
}

void myObjType::readFile(char* filename)
{
	cout << "Opening " << filename << endl;
	ifstream inFile;
	inFile.open(filename);
	if (!inFile.is_open()) {
		cout << "We cannot find your file " << filename << endl;
		exit(1);
	}

	string line;
	int i, j;
	bool firstVertex = 1;
	double currCood;

	while (getline(inFile, line))
	{
		if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
		{
			if (line[0] == 'v')
			{
				vcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ') j++;
					currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
					if (firstVertex) 
						lmin[k] = lmax[k] = currCood;
					else {
						if (lmin[k] > currCood)
							lmin[k] = currCood;
						if (lmax[k] < currCood)
							lmax[k] = currCood;
					}
					i = j;
				}

				firstVertex = 0;
			}
			if (line[0] == 'f')
			{
				tcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) {
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ' && linec[j] != '\\') j++;
					tlist[tcount][k] = atof(line.substr(i, j - i).c_str());
					i = j;
					fnlist[tcount][k] = 0;
					while (linec[j] != ' ') j++;

				}

			}


		}
	}

	// We suggest you to compute the normals here
    for (int i = 1; i <= tcount; i++) {
        int* vertices = tlist[i];
        vec3 v0 = vec3(vlist[vertices[0]][0], vlist[vertices[0]][1], vlist[vertices[0]][2]);
        vec3 v1 = vec3(vlist[vertices[1]][0], vlist[vertices[1]][1], vlist[vertices[1]][2]);
        vec3 v2 = vec3(vlist[vertices[2]][0], vlist[vertices[2]][1], vlist[vertices[2]][2]);

        vec3 v01 = v1 - v0;
        vec3 v02 = v2 - v0;
        vec3 crossProd = cross(v01, v02);
        crossProd.normalize();

        crossProd.copyToArray(nlist[i]);
    }

    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();
}



void myObjType::computeStat()
{
	int i;
    double minAngle = 180;
    double maxAngle = 0;

    for (int i = 1; i <= tcount; i++) {
        int* vertices = tlist[i];
        vec3 v0 = vec3(vlist[vertices[0]][0], vlist[vertices[0]][1], vlist[vertices[0]][2]);
        vec3 v1 = vec3(vlist[vertices[1]][0], vlist[vertices[1]][1], vlist[vertices[1]][2]);
        vec3 v2 = vec3(vlist[vertices[2]][0], vlist[vertices[2]][1], vlist[vertices[2]][2]);

        vec3 v01 = v1 - v0;
        vec3 v02 = v2 - v0;
        
        vec3 v12 = v2 - v1;
        vec3 v10 = -v01;

        double dot102 = dot(v01, v02);
        double angle102 = rad2Deg(acos(dot102 / (magnitude(v01) * magnitude(v02))));

        double dot210 = dot(v12, v10);
        double angle210 = rad2Deg(acos(dot210 / (magnitude(v12) * magnitude(v10))));

        double angle120 = 180.0f - angle102 - angle210;

        double tMinAngle = min(angle102, min(angle210, angle120));
        double tMaxAngle = max(angle102, max(angle210, angle120));

        if (tMinAngle < minAngle) {
            minAngle = tMinAngle;
        }
        int minBucketIndex = tMinAngle / 10.0f;
        ++(statMinAngle[minBucketIndex]);

        if (tMaxAngle > maxAngle) {
            maxAngle = tMaxAngle;
        }
        int maxBucketIndex = tMaxAngle / 10;
        ++(statMaxAngle[maxBucketIndex]);
    }

    cout << "Min. angle = " << minAngle << endl;
    cout << "Max. angle = " << maxAngle << endl;

	cout << "Statistics for Maximum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMaxAngle[i] << " ";
	cout << endl;
	cout << "Statistics for Minimum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMinAngle[i] << " ";
	cout << endl;
}

int myObjType::org(OrTri ot)
{
    int version = ver(ot);
    int* tri = tlist[idx(ot)];

    switch (version) {
    case 0:
    case 5:
        return tri[0];
    case 1:
    case 3:
        return tri[1];
    case 2:
    case 4:
        return tri[2];
    }
}

int myObjType::dest(OrTri ot)
{
    int version = ver(ot);
    int* tri = tlist[idx(ot)];

    switch (version) {
    case 2:
    case 3:
        return tri[0];
    case 0:
    case 4:
        return tri[1];
    case 1:
    case 5:
        return tri[2];
    }
}
