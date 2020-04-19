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
#include <assert.h>
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

template<class T>
ostream& operator<<(ostream& os, const unordered_set<T> uos) {
    os << "[";
    for (auto it = uos.begin(); it != uos.end(); ++it) {
        os << *it << ", ";
    }
    os << "]";

    return os;
}

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
        if (!isValidTriangle(i)) continue;
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

void myObjType::goToPreviousMeshVersion()
{
    if (meshVersions.size() == 0) {
        cout << "Versioned meshes not ready yet! Please do remeshing first with the key 'G'!" << endl;
        return;
    }


    int nextVersion = (meshVersions.size() + currentMeshVersionIndex - 1) % (meshVersions.size());
    cout << "--------- Going to mesh version " << nextVersion << " --------------------" << endl;
    currentMeshVersionIndex = nextVersion;
    restoreFromSavedObject(meshVersions[nextVersion]);
}

void myObjType::goToNextMeshVersion()
{
    if (meshVersions.size() == 0) {
        cout << "Versioned meshes not ready yet! Please do remeshing first with the key 'G'!" << endl;
        return;
    }
    int nextVersion = (currentMeshVersionIndex + 1) % (meshVersions.size());
    cout << "--------- Going to mesh version " << nextVersion << " --------------------" << endl;
    currentMeshVersionIndex = nextVersion;
    restoreFromSavedObject(meshVersions[nextVersion]);
}

void myObjType::restoreFromSavedObject(const savedObject& obj)
{
    std::copy(&obj.vlist[0][0], &obj.vlist[0][0] + MAXV * 3, &vlist[0][0]);
    std::copy(&obj.tlist[0][0], &obj.tlist[0][0] + MAXT * 3, &tlist[0][0]);
    unassignedTris = unordered_set<int>(obj.unassignedTris);
    unassignedVerts = unordered_set<int>(obj.unassignedVerts);
    tcount = obj.tcount;
    vcount = obj.vcount;

    initializeMesh();

    cout << "########## Mesh version " << currentMeshVersionIndex << "'s stats ##########" << endl;
    printStats();
}

void myObjType::initializeMesh()
{
    initializeVertexColors();

    {
        ScopedTimer timer("vertex to Triangles map, and vertex degrees list");

        // reset vertexDegreeList
        std::fill_n(vertexDegreeList, vcount + 1, 0);
        std::fill_n(vertexDegreeList + vcount + 1, MAXV - (vcount + 1), -1); // non-existant vertices 

        for (int i = 1; i <= tcount; i++) {
            int* vertices = tlist[i];
            vertexToTriangles[vertices[0]].emplace(i);
            vertexToTriangles[vertices[1]].emplace(i);
            vertexToTriangles[vertices[2]].emplace(i);

            // todo: move these to be computed when we wanna do remeshing.
            ++vertexDegreeList[vertices[0]];
            ++vertexDegreeList[vertices[1]];
            ++vertexDegreeList[vertices[2]];

            vertexLinks[vertices[0]].emplace(vertices[1]);
            vertexLinks[vertices[0]].emplace(vertices[2]);
            vertexLinks[vertices[1]].emplace(vertices[0]);
            vertexLinks[vertices[1]].emplace(vertices[2]);
            vertexLinks[vertices[2]].emplace(vertices[1]);
            vertexLinks[vertices[2]].emplace(vertices[2]);

            edgeSet.emplace(vertices[0], vertices[1]);
            edgeSet.emplace(vertices[1], vertices[2]);
            edgeSet.emplace(vertices[2], vertices[0]);

            edgeLinks[{vertices[0], vertices[1]}].emplace(vertices[2]);
            edgeLinks[{vertices[1], vertices[2]}].emplace(vertices[0]);
            edgeLinks[{vertices[2], vertices[0]}].emplace(vertices[1]);
        }
    }

    computeFNext();

    // Lab 2 Optional: Computing number of components
    computeNumComponents();

    computeTriangleNormals();

    computeVertexNormals();

    computeStat();
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
        }
        else if (ext == ".ply") {
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

    initializeMesh();
}

void myObjType::computeVertexNormals()
{
    ScopedTimer timer("Vertex normal computation");

    for (int i = 1; i <= vcount; i++) {
        computeVertexNormal(i);
    }
}

void myObjType::computeTriangleNormals()
{
    ScopedTimer timer("Normal computation");

    for (int i = 1; i <= tcount; i++) {
        computeTriangleNormal(i);
    }
}

void myObjType::initializeVertexColors()
{
    for (int i = 1; i <= tcount; i++) {
        colorlist[i][0] = 0.5f;
        colorlist[i][1] = 0.5f;
        colorlist[i][2] = 0.5f;
        colorlist[i][3] = 1.0f;
    }
}

void myObjType::computeFNext()
{
    ScopedTimer timer("fnext Population");

    std::fill(&fnlist[0][0], &fnlist[0][0] + MAXT * NUM_FNEXT, 0);
    std::fill(&isBoundaryVertexList[0], &isBoundaryVertexList[0] + sizeof(isBoundaryVertexList), false);

    // Lab 2 (main): Populating fnext list
    for (int i = 1; i <= tcount; i++) {
        if (!isValidTriangle(i)) continue;
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

            unordered_set<int> orgVertexTris = vertexToTriangles[origin];
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
                isBoundaryVertexList[origin] = true;
                isBoundaryVertexList[destination] = true;

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
    int unvisitedTCount = tcount - unassignedTris.size();
    bool visitedTris[MAXT] = { 0 };

    while (unvisitedTCount > 0) {
        int startingTri = 0;
        for (int i = 1; i <= tcount; i++) {
            if (!isValidTriangle(i)) continue;
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
    minAngle = 180;
    maxAngle = 0;

    for (int i = 1; i <= tcount; i++) {
        if (!isValidTriangle(i)) continue;
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
        if (minBucketIndex >= 0 && minBucketIndex < 18) {
            ++(statMinAngle[minBucketIndex]);
        }

        if (tMaxAngle > maxAngle) {
            maxAngle = tMaxAngle;
        }
        int maxBucketIndex = tMaxAngle / 10.0f;
        if (maxBucketIndex >= 0 && maxBucketIndex < 18) {
            ++(statMaxAngle[maxBucketIndex]);
        }
    }
}

void myObjType::printStats()
{
    cout << "------------ START STATS -----------------" << endl;

    cout << "Tris: " << tcount - unassignedVerts.size() << " // Verts: " << vcount - unassignedTris.size() << endl;

    cout << "Min. angle = " << minAngle << endl;
    cout << "Max. angle = " << maxAngle << endl;

    cout << "Statistics for Maximum Angles" << endl;
    for (int i = 0; i < 18; i++)
        cout << statMaxAngle[i] << " ";
    cout << endl;
    cout << "Statistics for Minimum Angles" << endl;
    for (int i = 0; i < 18; i++)
        cout << statMinAngle[i] << " ";
    cout << endl;

    cout << "Number of components: " << numComponents << endl;

    cout << "-------------- END STATS ------------------" << endl;
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

bool myObjType::isValidTriangle(int index)
{
    if (unassignedTris.count(index) == 1) return false;
    int* vertices = tlist[index];
    return !(vertices[0] == 0 || vertices[1] == 0 || vertices[2] == 0);
}

int myObjType::AddTriangle(int vertices[3])
{
    int index;
    if (unassignedTris.size() > 0) {
        index = *unassignedTris.begin();
        unassignedTris.erase(unassignedTris.begin());
    }
    else {
        ++tcount;
        index = tcount;
    }

    vertexToTriangles[vertices[0]].emplace(index);
    vertexToTriangles[vertices[1]].emplace(index);
    vertexToTriangles[vertices[2]].emplace(index);

    tlist[index][0] = vertices[0];
    tlist[index][1] = vertices[1];
    tlist[index][2] = vertices[2];

    ++vertexDegreeList[vertices[0]];
    ++vertexDegreeList[vertices[1]];
    ++vertexDegreeList[vertices[2]];

    vertexLinks[vertices[0]].emplace(vertices[1]);
    vertexLinks[vertices[0]].emplace(vertices[2]);
    vertexLinks[vertices[1]].emplace(vertices[0]);
    vertexLinks[vertices[1]].emplace(vertices[2]);
    vertexLinks[vertices[2]].emplace(vertices[1]);
    vertexLinks[vertices[2]].emplace(vertices[2]);

    // Edges
    edgeSet.emplace(vertices[0], vertices[1]);
    edgeSet.emplace(vertices[1], vertices[2]);
    edgeSet.emplace(vertices[2], vertices[0]);

    edgeLinks[{vertices[0], vertices[1]}].emplace(vertices[2]);
    edgeLinks[{vertices[1], vertices[2]}].emplace(vertices[0]);
    edgeLinks[{vertices[2], vertices[0]}].emplace(vertices[1]);

    //cout << "Added triangle ";
    //printOrTri(makeOrTri(index, 0));
    //cout << " at index " << index << endl;

    setVertexColor(index, 0.5f, 0.5f, 0.5f);

    return index;
}

bool myObjType::AddTriangleAtIndex(int index, int vertices[3])
{
    auto it = unassignedTris.find(index);
    if (it == unassignedTris.end()) {
        cout << "Trying to add a triangle at an invalid index " << index << endl;
        return false;
    }

    unassignedTris.erase(it);

    vertexToTriangles[vertices[0]].emplace(index);
    vertexToTriangles[vertices[1]].emplace(index);
    vertexToTriangles[vertices[2]].emplace(index);

    tlist[index][0] = vertices[0];
    tlist[index][1] = vertices[1];
    tlist[index][2] = vertices[2];

    ++vertexDegreeList[vertices[0]];
    ++vertexDegreeList[vertices[1]];
    ++vertexDegreeList[vertices[2]];

    vertexLinks[vertices[0]].emplace(vertices[1]);
    vertexLinks[vertices[0]].emplace(vertices[2]);
    vertexLinks[vertices[1]].emplace(vertices[0]);
    vertexLinks[vertices[1]].emplace(vertices[2]);
    vertexLinks[vertices[2]].emplace(vertices[1]);
    vertexLinks[vertices[2]].emplace(vertices[2]);

    // Edges
    edgeSet.emplace(vertices[0], vertices[1]);
    edgeSet.emplace(vertices[1], vertices[2]);
    edgeSet.emplace(vertices[2], vertices[0]);

    edgeLinks[{vertices[0], vertices[1]}].emplace(vertices[2]);
    edgeLinks[{vertices[1], vertices[2]}].emplace(vertices[0]);
    edgeLinks[{vertices[2], vertices[0]}].emplace(vertices[1]);

    //cout << "Added triangle ";
    //printOrTri(makeOrTri(index, 0));
    //cout << " at index " << index << endl;

    setVertexColor(index, 0.5f, 0.5f, 0.5f);

    return true;
}

void myObjType::RemoveTriangleAtIndex(int index)
{
    if (!isValidTriangle(index)) {
        cout << "Tried removing an invalid triangle at index " << index << endl;
        return;
    }

    unassignedTris.emplace(index);

    int* vertices = tlist[index];

    vertexToTriangles[vertices[0]].erase(index);
    vertexToTriangles[vertices[1]].erase(index);
    vertexToTriangles[vertices[2]].erase(index);

    --vertexDegreeList[vertices[0]];
    --vertexDegreeList[vertices[1]];
    --vertexDegreeList[vertices[2]];

    vertexLinks[vertices[0]].erase(vertices[1]);
    vertexLinks[vertices[0]].erase(vertices[2]);
    vertexLinks[vertices[1]].erase(vertices[0]);
    vertexLinks[vertices[1]].erase(vertices[2]);
    vertexLinks[vertices[2]].erase(vertices[1]);
    vertexLinks[vertices[2]].erase(vertices[2]);

    // Edges
    edgeSet.erase({ vertices[0], vertices[1] });
    edgeSet.erase({ vertices[1], vertices[2] });
    edgeSet.erase({ vertices[2], vertices[0] });

    edgeLinks[{vertices[0], vertices[1]}].erase(vertices[2]);
    edgeLinks[{vertices[1], vertices[2]}].erase(vertices[0]);
    edgeLinks[{vertices[2], vertices[0]}].erase(vertices[1]);

    //cout << "Removed triangle ";
    //printOrTri(makeOrTri(index, 0));
    //cout << " at index " << index << endl;

    tlist[index][0] = 0;
    tlist[index][1] = 0;
    tlist[index][2] = 0;
}

bool myObjType::isValidVertex(int index)
{
    return !(index > vcount || unassignedVerts.count(index) > 0);
}

int myObjType::addVertex(vec3 position)
{
    int index;
    if (unassignedVerts.size() > 0) {
        index = *unassignedVerts.begin();
        unassignedVerts.erase(unassignedVerts.begin());
    }
    else {
        ++vcount;
        index = vcount;
    }

    position.copyToArray(vlist[index]);

    return index;
}

void myObjType::setVertexPosition(int index, vec3 position)
{
    if (!isValidVertex(index)) {
        //cout << "Attempted to set position of vertex that does not exist: " << index << endl;
        return;
    }

    position.copyToArray(vlist[index]);
}

void myObjType::removeVertexOnly(int index)
{
    if (!isValidVertex(index)) {
        //cout << "Attempted to remove vertex that does not exist: " << index << endl;
        return;
    }

    unassignedVerts.emplace(index);
}

bool myObjType::IsEdgeContractable(Edge edge)
{
    //cout << "--------- Checking if " << edge.ToString() << " is contractable... ---------" << endl;

    int v1 = edge.v1;
    int v2 = edge.v2;

    auto edgeLinkIt = edgeLinks.find(edge);
    if (edgeLinkIt == edgeLinks.end()) return false;
    unordered_set<int> linkEdge = edgeLinkIt->second;

    //cout << "Edge's Link: " << linkEdge << endl;

    auto v1LinkIt = vertexLinks.find(v1);
    if (v1LinkIt == vertexLinks.end()) return false;
    unordered_set<int> linkV1 = v1LinkIt->second;

    //cout << "Vertex " << v1 << "'s Link: " << linkV1 << endl;

    auto v2LinkIt = vertexLinks.find(v2);
    if (v2LinkIt == vertexLinks.end()) return false;
    unordered_set<int> linkV2 = v2LinkIt->second;

    //cout << "Vertex " << v2 << "'s Link: " << linkV2 << endl;

    unordered_set<int> intersect;
    for (int lk : linkV1)
    {
        if (linkV2.count(lk) != 0) {
            intersect.emplace(lk);
        }
    }

    //cout << "Intersect's Link: " << intersect << endl;

    // check if intersect == edge's link
    if (linkEdge.size() != intersect.size()) return false;

    for (int lk : linkEdge) {
        if (intersect.count(lk) == 0) return false;
    }

    //cout << "cxb" << endl;

    return true;
}

void myObjType::CollapseEdge(Edge edge)
{
    //cout << "--------- Contracting " << edge.ToString() << " ---------" << endl;
    int retainVertex;
    int removeVertex;

    // Keep and move the vertex that has more triangles
    if (vertexToTriangles[edge.v1].size() >= vertexToTriangles[edge.v2].size()) {
        retainVertex = edge.v1;
        removeVertex = edge.v2;
    }
    else {
        retainVertex = edge.v2;
        removeVertex = edge.v1;
    }

    unordered_set<int> toRemove;
    for (int v : vertexToTriangles[removeVertex]) {
        if (vertexToTriangles[retainVertex].find(v) != vertexToTriangles[retainVertex].end()) {
            toRemove.emplace(v);
        }
    }

    // Remove the common triangles first (that will be contracted away)
    for (int v : toRemove) {
        RemoveTriangleAtIndex(v);
    }

    // move v1 to in between v1 and v2.
    vec3 locationDiff = { vlist[removeVertex][0] - vlist[retainVertex][0],
                            vlist[removeVertex][1] - vlist[retainVertex][1],
                            vlist[removeVertex][2] - vlist[retainVertex][2] };
    locationDiff *= 0.5f; // move halfway


    // set new position
    vlist[retainVertex][0] += locationDiff.x;
    vlist[retainVertex][1] += locationDiff.y;
    vlist[retainVertex][2] += locationDiff.z;

    // Now we need to transfer removeVertex's triangles over to retainVertex
    unordered_set<int> toTransfer = vertexToTriangles[removeVertex];

    for (int triIndex : toTransfer) {
        int* vertices = tlist[triIndex];
        int savedTri[3] = { 0 };
        for (int i = 0; i < 3; i++) {
            // changing removeVertex in the triangle to retainVertex
            if (vertices[i] == removeVertex) {
                savedTri[i] = retainVertex;
            }
            else {
                savedTri[i] = vertices[i];
            }
        }

        RemoveTriangleAtIndex(triIndex);

        AddTriangleAtIndex(triIndex, savedTri);
    }

    edgeSet.erase(edge);

    removeVertexOnly(removeVertex);
}

void myObjType::SplitEdge(Edge edge)
{
    int v1 = edge.v1;
    int v2 = edge.v2;

    // Creating the new vertex
    // move v1 to in between v1 and v2.
    vec3 locationDiff = { vlist[v2][0] - vlist[v1][0],
                            vlist[v2][1] - vlist[v1][1],
                            vlist[v2][2] - vlist[v1][2] };
    locationDiff *= 0.5f; // move halfway

    double vPos[3] = { vlist[v1][0] + locationDiff.x,
                        vlist[v1][1] + locationDiff.y,
                        vlist[v1][2] + locationDiff.z };

    int v3 = addVertex(vPos);

    unordered_set<int> toSplit;
    // finding the common triangles to split. (and remove before adding the new ones)
    for (int v : vertexToTriangles[v1]) {
        if (vertexToTriangles[v2].find(v) != vertexToTriangles[v2].end()) {
            toSplit.emplace(v);
        }
    }

    // splitting a triangle into two.
    for (int t : toSplit) {
        // replacing v1 with v3
        int triV1Replace[3] = { tlist[t][0], tlist[t][1], tlist[t][2] };
        for (int i = 0; i < 3; i++) {
            if (triV1Replace[i] == v1) {
                triV1Replace[i] = v3;
                break;
            }
        }

        // replacing v2 with v3
        int triV2Replace[3] = { tlist[t][0], tlist[t][1], tlist[t][2] };
        for (int i = 0; i < 3; i++) {
            if (triV2Replace[i] == v2) {
                triV2Replace[i] = v3;
                break;
            }
        }

        RemoveTriangleAtIndex(t);
        AddTriangle(triV1Replace);
        AddTriangle(triV2Replace);
    }
}

bool myObjType::ShouldFlipEdge(const Edge& edge)
{
    auto it = edgeLinks.find(edge);
    if (it == edgeLinks.end()) {
        return false;
    }

    if (it->second.size() < 2) {
        return false;
    }

    int a1 = edge.v1;
    int a2 = edge.v2;
    int b1 = *it->second.begin();
    int b2 = *(++it->second.begin());

    bool idealDegreea1 = isBoundaryVertexList[a1] ? 4 : 6;
    bool idealDegreea2 = isBoundaryVertexList[a2] ? 4 : 6;
    bool idealDegreeb1 = isBoundaryVertexList[b1] ? 4 : 6;
    bool idealDegreeb2 = isBoundaryVertexList[b2] ? 4 : 6;

    int preflipDeviation = abs(vertexDegreeList[a1] - idealDegreea1) +
        abs(vertexDegreeList[a2] - idealDegreea2) +
        abs(vertexDegreeList[b1] - idealDegreeb1) +
        abs(vertexDegreeList[b2] - idealDegreeb2);

    int postflipDeviation = abs(vertexDegreeList[a1] - 1 - idealDegreea1) +
        abs(vertexDegreeList[a2] - 1 - idealDegreea2) +
        abs(vertexDegreeList[b1] + 1 - idealDegreeb1) +
        abs(vertexDegreeList[b2] + 1 - idealDegreeb2);

    return postflipDeviation < preflipDeviation;
}

void myObjType::FlipEdge(const Edge& edge)
{
    auto it = edgeLinks.find(edge);
    if (it == edgeLinks.end()) {
        return;
    }

    if (it->second.size() < 2) {
        return;
    }

    int a1 = edge.v1;
    int a2 = edge.v2;
    int b1 = *it->second.begin();
    int b2 = *(++it->second.begin());

    // finding the tris to remove
    unordered_set<int> toRemove;

    // finding the common triangles to split. (and remove before adding the new ones)
    for (int v : vertexToTriangles[a1]) {
        if (vertexToTriangles[a2].find(v) != vertexToTriangles[a2].end()) {
            toRemove.emplace(v);
        }
    }

    if (toRemove.size() != 2) {
        return;
    }

    int triReplace1[3];
    int triReplace2[3];

    // Checking triangle's normal in order to construct to correct triangle
    vec3 a2Pos = vec3(vlist[a2][0], vlist[a2][1], vlist[a2][2]);
    vec3 b1Pos = vec3(vlist[b1][0], vlist[b1][1], vlist[b1][2]);
    vec3 a1Pos = vec3(vlist[a1][0], vlist[a1][1], vlist[a1][2]);

    vec3 v01 = a2Pos - b1Pos;
    vec3 v02 = a1Pos - b1Pos;
    vec3 crossProd = cross(v01, v02);

    vec3 currNormal = nlist[*toRemove.begin()];

    if (dot(crossProd, currNormal) >= 0) { // pointing in same direction
        triReplace1[0] = a1; triReplace1[1] = b1; triReplace1[2] = b2;
        triReplace2[0] = a2; triReplace2[1] = b2; triReplace2[2] = b1;
    }
    else {
        triReplace1[0] = a1; triReplace1[1] = b2; triReplace1[2] = b1;
        triReplace2[0] = a2; triReplace2[1] = b1; triReplace2[2] = b2;
    }

    for (int t : toRemove) {
        RemoveTriangleAtIndex(t);
    }

    computeTriangleNormal(AddTriangle(triReplace1));
    computeTriangleNormal(AddTriangle(triReplace2));

    computeVertexNormal(a1);
    computeVertexNormal(a2);
    computeVertexNormal(b1);
    computeVertexNormal(b2);
}

vec3 myObjType::computeVertexCentroid(int vertexIndex)
{
    if (isBoundaryVertexList[vertexIndex]) return vlist[vertexIndex]; // don't move boundary vertices

    const unordered_set<int>& neighbours = vertexLinks[vertexIndex];

    vec3 currentPos = vlist[vertexIndex];

    vec3 centroid;
    for (int vert : neighbours) {
        centroid.x += vlist[vert][0];
        centroid.y += vlist[vert][1];
        centroid.z += vlist[vert][2];
    }

    centroid /= (double)neighbours.size();

    return centroid;
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

void myObjType::computeTriangleNormal(int triIndex)
{
    int* vertices = tlist[triIndex];
    vec3 v0 = vec3(vlist[vertices[0]][0], vlist[vertices[0]][1], vlist[vertices[0]][2]);
    vec3 v1 = vec3(vlist[vertices[1]][0], vlist[vertices[1]][1], vlist[vertices[1]][2]);
    vec3 v2 = vec3(vlist[vertices[2]][0], vlist[vertices[2]][1], vlist[vertices[2]][2]);

    vec3 v01 = v1 - v0;
    vec3 v02 = v2 - v0;
    vec3 crossProd = cross(v01, v02);
    crossProd.normalize();

    crossProd.copyToArray(nlist[triIndex]);
}

void myObjType::computeVertexNormal(int vertexIndex)
{
    unordered_set<int> tris = vertexToTriangles[vertexIndex];
    vec3 vNormal(0.0f, 0.0f, 0.0f);
    for (int tri : tris) {
        vNormal.x += nlist[tri][0];
        vNormal.y += nlist[tri][1];
        vNormal.z += nlist[tri][2];
    }
    vNormal /= tris.size();
    vNormal.normalize();
    vNormal.copyToArray(vnlist[vertexIndex]);
}

void myObjType::performRemeshing(int numIterations)
{
    ScopedTimer timer("Area-Equalizing Remeshing routine");
    cout << "Starting to do Remeshing with " << numIterations << " iterations." << endl;

    double runningLengthSum = 0.0f;
    for (const Edge& edge : edgeSet) {
        double len = magnitude(vec3(vlist[edge.v1]) - vec3(vlist[edge.v2]));
        runningLengthSum += len;
    }

    double meanEdgeLength = runningLengthSum / edgeSet.size();

    double splitThreshold = 4.0f * meanEdgeLength / 3.0f;
    double collapseThreshold = 4.0f * meanEdgeLength / 5.0f;

    for (int itNum = 0; itNum < numIterations; itNum++) {
        meshVersions.emplace_back(*this);

        unordered_map<Edge, double> edgeLengths;
        for (const Edge& edge : edgeSet) {
            double len = magnitude(vec3(vlist[edge.v1]) - vec3(vlist[edge.v2]));
            edgeLengths.emplace(edge, len);
        }

        unordered_set<Edge> toSplit;
        unordered_set<Edge> toCollapse;

        for (const auto& kv : edgeLengths)
        {
            if (kv.second > splitThreshold) {
                toSplit.emplace(kv.first);
            }
            else if (kv.second < collapseThreshold) {
                toCollapse.emplace(kv.first);
            }
        }

        // Step 1: Split edges that are too long
        for (const Edge& edge : toSplit) {
            SplitEdge(edge);
        }

        // Step 2: Collapse edges that are too small.
        for (const Edge& edge : toCollapse) {
            if (IsEdgeContractable(edge)) {
                CollapseEdge(edge);
            }
        }

        // Step 3: Flip edges that are good to flip
        unordered_set<Edge> edgesCopy = unordered_set<Edge>(edgeSet);
        for (const Edge& edge : edgesCopy) {
            if (ShouldFlipEdge(edge)) {
                FlipEdge(edge);
            }
        }

        // Step 4: Vertex Averaging. Moving vertex towards centroid tangentially
        for (int i = 1; i <= vcount; i++) {
            if (isValidVertex(i)) {
                vec3 centroidPos = computeVertexCentroid(i);
                centroidPos.copyToArray(centroids[i]);
            }
        }

        // Step 4.2: Apply the computed vertex positions
        for (int i = 1; i <= vcount; i++) {
            if (isValidVertex(i)) {
                vlist[i][0] = centroids[i][0];
                vlist[i][1] = centroids[i][1];
                vlist[i][2] = centroids[i][2];
            }
        }

        computeTriangleNormals();
        computeVertexNormals();
        computeStat();
        printStats();
    }

    meshVersions.emplace_back(*this);
    currentMeshVersionIndex = meshVersions.size() - 1;

}

void myObjType::setVertexColor(int vIdx, float r, float g, float b)
{
    colorlist[vIdx][0] = r;
    colorlist[vIdx][1] = g;
    colorlist[vIdx][2] = b;
    colorlist[vIdx][3] = 1.0f;
}

savedObject::savedObject(myObjType &obj)
{
    std::copy(&obj.vlist[0][0], &obj.vlist[0][0] + MAXV * 3, &vlist[0][0]);
    std::copy(&obj.tlist[0][0], &obj.tlist[0][0] + MAXT * 3, &tlist[0][0]);
    unassignedTris = unordered_set<int>(obj.unassignedTris);
    unassignedVerts = unordered_set<int>(obj.unassignedVerts);
    tcount = obj.tcount;
    vcount = obj.vcount;
}
