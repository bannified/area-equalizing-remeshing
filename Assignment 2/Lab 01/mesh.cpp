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
#include <algorithm>
#include <iterator>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "mesh.h"
#include "ScopedTimer.h"
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
    {
        ScopedTimer timer("File parsing");

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
                        while (linec[j] != ' ') j++;

                    }
                }
            }
        }
    }

    {
        ScopedTimer timer("fnext Population");

        std::fill(fnlist[0], fnlist[0] + MAXT * NUM_FNEXT, 0);

        // Lab 2 (main): Populating fnext list
        for (int i = 1; i <= tcount; i++) {
            //std::cout << "Populating fnext for triangle: " << i << "/" << tcount << std::endl;
            int* vertices = tlist[i];
            // Fill fnext list for face i.
            for (int k = 0; k < NUM_FNEXT; k++) {
                if (fnlist[i][k] != 0) {
                    // already filled.
                    continue;
                }

                OrTri ot = makeOrTri(i, k);
                int origin = org(ot);
                int destination = dest(ot);

                int fnextTriIdx = 0;
                // finding index of triangle that shares the same edge.
                for (int curr = 1; curr <= tcount; curr++) {
                    if (curr == i) continue;
                    int* currTriVerts = tlist[curr];

                    bool foundOrg = false, foundDest = false;

                    for (int vert = 0; vert < 3; vert++) {
                        int vertIdx = currTriVerts[vert];
                        if (!foundOrg && vertIdx == origin) {
                            foundOrg = true;
                            continue;
                        }

                        if (!foundDest && vertIdx == destination) {
                            foundDest = true;
                            continue;
                        }
                    }

                    if (foundDest && foundOrg) {
                        fnextTriIdx = curr;
                        break;
                    }
                }

                if (fnextTriIdx == 0) {
                    // no triangle found for this fnext
                    fnlist[i][k] = 0;
                    continue;
                }

                // Finding triangle version
                for (int version = 0; version < 6; version++) {
                    OrTri fnextOt = makeOrTri(fnextTriIdx, version);
                    if (dest(fnextOt) == destination && org(fnextOt) == origin) {
                        fnlist[i][k] = fnextOt;
                        fnlist[fnextTriIdx][version] = ot;
                        break;
                    }
                }
            }
        }
    }
    
    //printfnList();

    // Lab 2 Optional: Computing number of components
    {
        ScopedTimer timer("Computing number of components");

        numComponents = 0;
        int unvisitedTCount = tcount;
        bool visitedTris[MAXT] = { 0 };

        while (unvisitedTCount > 0) {
            int startingTri = 0;
            for (int i = 1; i <= tcount; i++) {
                if (!visitedTris[i]) {
                    startingTri = i;
                    break;
                }
            }

            ++numComponents;

            // start BFS, ignoring nodes already at frontier
            std::queue<int> frontier;
            frontier.push(startingTri);

            while (!frontier.empty()) {
                int curr = frontier.front();
                frontier.pop();

                if (visitedTris[curr]) {
                    continue;
                }

                int* fnexts = fnlist[curr];
                for (int f = 0; f < NUM_FNEXT; f++) {
                    int index = idx(fnexts[f]);
                    if (index > 0 && !visitedTris[index]) {
                        frontier.push(index);
                    }
                }
                triComponentNumber[curr] = numComponents;
                visitedTris[curr] = true;
                --unvisitedTCount;
                //std::cout << "triangles left: " << unvisitedTCount << "/" << tcount << std::endl;
            }
        }
    }

    std::cout << "Number of components: " << numComponents << std::endl;

    {
        ScopedTimer timer("Normal computation");

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
    }

    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    {
        ScopedTimer timer("Stats computation");
        computeStat();
    }
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

int myObjType::last(OrTri ot)
{
    int version = ver(ot);
    int* tri = tlist[idx(ot)];

    switch (version) {
    case 1:
    case 4:
        return tri[0];
    case 2:
    case 5:
        return tri[1];
    case 0:
    case 3:
        return tri[2];
    }
}

void myObjType::printfnList()
{
    for (int i = 1; i <= tcount; i++) {
        int* trifn = fnlist[i];

        std::cout << "Triangle " << i << "(";
        printOrTri(makeOrTri(i, 0));
        std::cout << ") : ";
        for (int k = 0; k < 6; k++) {
            printOrTri(trifn[k]);
            std::cout << "    ";
        }

        std::cout << endl;
    }
}

void myObjType::printOrTri(OrTri ot)
{
    std::cout << org(ot) << "|" << dest(ot) << "|" << last(ot);
}
