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
#include "limits"
#include <algorithm>
#include <iterator>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "ScopedTimer.h"
#include <map>
#include <queue>
#include <string>
#include <iomanip>
#include <utility>
#include <unordered_map>
using namespace std;

const string gPlyElementString = "element ";
const string gPlyVertexString = "vertex ";
const string gPlyPropertyString = "property ";
const string gPlyFaceString = "face ";
const string gPlyFloatString = "float ";
const string gPlayFloat32String = "float32 ";
const string gPlayEndHeaderString = "end_header";

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
        if (!smooth) glNormal3dv(nlist[i]);
        for (int j = 0; j < 3; j++) {
            if (smooth) glNormal3dv(vnlist[tlist[i][j]]);
            glMaterialfv(GL_FRONT, GL_AMBIENT, colorlist[tlist[i][j]]);
            glVertex3dv(vlist[tlist[i][j]]);
        }
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

    std::string fn(filename);
    std::string ext = fn.substr(fn.find_last_of('.'));

    string line;
    int i, j;
    bool firstVertex = 1;
    double currCood;
    {
        ScopedTimer timer("File parsing");

        if (ext == ".obj") {
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
        } else if (ext == ".ply") {
            // parse .ply

            int numVertexElements = 0;
            int numVertices = 0;
            int numFaces = 0;

            // header
            while (getline(inFile, line)) {
                size_t pos = line.find(gPlayEndHeaderString);
                if (pos != string::npos) {
                    break; // header ended
                }

                // handling elements
                pos = line.find(gPlyElementString);
                if (pos != string::npos) {
                    line = line.substr(pos + gPlyElementString.size());

                    // find num vertices
                    pos = line.find(gPlyVertexString);
                    if (pos != string::npos) {
                        numVertices = atoi(line.substr(pos + gPlyVertexString.size()).c_str());
                        continue;
                    }

                    // find num faces/triangles
                    pos = line.find(gPlyFaceString);
                    if (pos != string::npos) {
                        numFaces = atoi(line.substr(pos + gPlyFaceString.size()).c_str());
                        continue;
                    }
                }
            }

            vcount = numVertices;
            for (int i = 1; i <= numVertices; i++) {
                getline(inFile, line);
                for (int j = 0; j < 3; j++) {
                    size_t pos = line.find(' ');
                    currCood = vlist[i][j] = atof(line.substr(0, pos).c_str());
                    if (vlist[i][j] == 0.0f) {
                        currCood = vlist[i][j] = numeric_limits<double>::epsilon();
                    }
                    line = line.substr(pos + 1);

                    if (firstVertex)
                        lmin[j] = lmax[j] = currCood;
                    else {
                        if (lmin[j] > currCood)
                            lmin[j] = currCood;
                        if (lmax[j] < currCood)
                            lmax[j] = currCood;
                    }
                }
                firstVertex = 0;
            }

            tcount = numFaces;
            for (int i = 1; i <= numFaces; i++) {
                getline(inFile, line);
                line = line.substr(2);
                for (int j = 0; j < 3; j++) {
                    size_t pos = line.find(' ');
                    int intendedIndex = atoi(line.substr(0, pos).c_str()) + 1;
                    tlist[i][j] = intendedIndex;
                    line = line.substr(pos + 1);
                }
            }
        }
        else {
            cout << "Invalid file format given, quitting..." << endl;
            exit(1);
        }
    }

    initializeVertexColors();

    {
        ScopedTimer timer("vertex to Triangles");
        for (int i = 1; i <= tcount; i++) {
            int* vertices = tlist[i];
            vertexToTriangles[vertices[0]].emplace_back(i);
            vertexToTriangles[vertices[1]].emplace_back(i);
            vertexToTriangles[vertices[2]].emplace_back(i);
        }
    }

    computeFNext();

    //printfnList();

    // Lab 2 Optional: Computing number of components
    computeNumComponents();
    cout << "Number of components: " << numComponents << endl;

    computeTriangleNormals();

    computeVertexNormals();

    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    
    computeStat();
}

void myObjType::computeVertexNormals()
{
    ScopedTimer timer("Vertex normal computation");

    for (int i = 1; i <= vcount; i++) {
        vector<int> tris = vertexToTriangles[i];
        vec3 vNormal(0.0f, 0.0f, 0.0f);
        for (int tri : tris) {
            vNormal.x += nlist[tri][0];
            vNormal.y += nlist[tri][1];
            vNormal.z += nlist[tri][2];
        }
        vNormal /= tris.size();
        vNormal.normalize();
        vNormal.copyToArray(vnlist[i]);
    }
}

void myObjType::computeTriangleNormals()
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

void myObjType::initializeVertexColors()
{
    for (int i = 1; i <= tcount; i++) {
        for (int j = 1; j <= 3; j++) {
            colorlist[i][j] = 0.5f;
        }
        colorlist[i][3] = 1.0f;
    }
}

void myObjType::computeFNext()
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

            vector<int> orgVertexTris = vertexToTriangles[origin];
            for (int tri : orgVertexTris) {
                if (tri == i) continue;
                int* ovtVerts = tlist[tri];
                for (int v = 0; v < 3; v++) {
                    if (ovtVerts[v] == destination) {
                        fnextTriIdx = tri;
                        break;
                    }
                }
                if (fnextTriIdx != 0) break;
            }

            if (fnextTriIdx == 0) {
                // no triangle found for this fnext
                fnlist[i][k] = ot;

                setVertexColor(origin, 1.0f, 0.0f, 0.0f);
                setVertexColor(destination, 1.0f, 0.0f, 0.0f);

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

void myObjType::computeNumComponents()
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

void myObjType::computeStat()
{
    ScopedTimer timer("Stats computation");
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

void myObjType::printVertexList()
{
    for (int i = 1; i <= vcount; i++) {
        double* tri = vlist[i];

        std::cout << "Triangle " << i << ": " << tri[0] << "|" << tri[1] << "|" << tri[2];
        std::cout << endl;
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

void myObjType::printTriList()
{
    for (int i = 1; i <= tcount; i++) {
        int* tri = tlist[i];

        std::cout << "Triangle " << i << "(";
        printOrTri(makeOrTri(i, 0));
        std::cout << ")";
    }
}

void myObjType::printOrTri(OrTri ot)
{
    std::cout << org(ot) << "|" << dest(ot) << "|" << last(ot);
}

void myObjType::performRemeshing(int numIterations)
{
    for (int itNum = 0; itNum < numIterations; itNum++) {

    }


}

void myObjType::splitAllLongEdges(float threshold)
{

}

void myObjType::mergeAllShortEdges(float threshold)
{

}

void myObjType::setVertexColor(int vIdx, float r, float g, float b)
{
    colorlist[vIdx][0] = r;
    colorlist[vIdx][1] = g;
    colorlist[vIdx][2] = b;
    colorlist[vIdx][3] = 1.0f;
}
