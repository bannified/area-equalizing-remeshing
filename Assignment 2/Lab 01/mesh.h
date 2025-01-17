#pragma once

#include <unordered_map>
#include <unordered_set>
#include <set>
#include "vec3.h"
#include "Edge.h"

// maximum number of vertices and triangles
#define MAXV 1000000
#define MAXT 1000000
#define NUM_FNEXT 6

typedef int OrTri;
typedef int tIdx;

inline OrTri makeOrTri(tIdx t, int version) { return (t << 3) | version; };
inline tIdx idx(OrTri ot) { return ot >> 3; };
inline int ver(OrTri ot) { return ot & 0b111; };
inline OrTri enext(OrTri ot) {
    int v = ver(ot);  return makeOrTri(idx(ot),
        v < 3 ? (v + 1) % 3 : 3 + ((v - 1) % 3));
};
inline OrTri sym(OrTri ot) { int v = ver(ot); return v < 3 ? ot + 3 : ot - 3; };

using namespace std;

class myObjType;

class savedObject {

public:
    int tcount = 0;
    int vcount = 0;
    double vlist[MAXV][3];   // vertices list
    int tlist[MAXT][3];      // triangle list

    /* For Addition/Removal of Tris and Vertices */
    unordered_set<int> unassignedTris;
    unordered_set<int> unassignedVerts;

    savedObject(myObjType &obj);
};

class myObjType {
public:
    int vcount = 0;
    int tcount = 0;

    double vlist[MAXV][3];   // vertices list
    int tlist[MAXT][3];      // triangle list
    int fnlist[MAXT][6];     // fnext list
    double nlist[MAXT][3];   // storing triangle normals
    double vnlist[MAXT][3];  // Vertex Normal vectors list
    bool isBoundaryVertexList[MAXV]; // isBoundary
    float colorlist[MAXT][4];   // colors list
    int triComponentNumber[MAXT];      // triangle list

    double lmax[3];          // the maximum coordinates of x,y,z
    double lmin[3];          // the minimum coordinates of x,y,z

    int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
    int statMaxAngle[18];

    float minAngle;
    float maxAngle;

    double centroids[MAXV][3];

    int numComponents;

    /* For Addition/Removal of Tris and Vertices */
    unordered_set<int> unassignedTris;
    unordered_set<int> unassignedVerts;

    /* todo: For remeshing */
    unordered_set<Edge> edgeSet;    // Set of all edges
    int vertexDegreeList[MAXV];     // keeps track of the degree of every vertex. 
                                    // if degree is -1, then vertex does not exist.
    unordered_map<int, unordered_set<int>> vertexLinks; // Lk(v), where v is a vertex
    unordered_map<Edge, unordered_set<int>> edgeLinks; // Lk(e), where e is an edge.

    std::unordered_map<int, std::unordered_set<int> > vertexToTriangles;

    int currentVertex = 1;

    vector<savedObject> meshVersions;

    int currentMeshVersionIndex = 0;
    void goToPreviousMeshVersion();
    void goToNextMeshVersion();

    void restoreFromSavedObject(const savedObject& obj);

    void initializeMesh();

public:
    bool smooth = false;

    myObjType() { vcount = 0; tcount = 0; };
    void readFile(char* filename);  // assumming file contains a manifold

    void writeFile(char* filename);
    void draw();

    int org(OrTri ot);
    int dest(OrTri ot);
    int last(OrTri ot);

    // Triangles manipulation. 
    bool isValidTriangle(int index);
    int AddTriangle(int vertices[3]); // returns the index of the triangle added
    bool AddTriangleAtIndex(int index, int vertices[3]); // returns true if successful
    void RemoveTriangleAtIndex(int index);

    // Vertex Manipulation.
    bool isValidVertex(int index);
    int addVertex(vec3 position); // returns the index of the vertex added
    void setVertexPosition(int index, vec3 position);
    void removeVertexOnly(int index); // plainly just removes the vertex

    /**
     * Edge Collapsing/Contracting
     * Should probably use IsEdgeCollapsible to check if an edge is collapsible before using CollapseEdge.
     */
    bool IsEdgeCollapsible(Edge edge);
    void CollapseEdge(Edge edge);

    /**
     * Splits an edge such that it becomes two edges.
     * Triangles that use this edge also get split accordingly.
     */
    void SplitEdge(Edge edge);

    /**
     * Edge Flipping
     * Should probably use ShouldFlip to check if an edge is collapsible before using CollapseEdge.
     */
    bool ShouldFlipEdge(const Edge& edge);
    void FlipEdge(const Edge& edge);

public:
    /**
     * Perform Area-Equalizing Remeshing.
     * Based on the paper: �A Remeshing Approach to Multiresolution Modeling� by Botsch & Kobbelt
     * 4-step approach:
     * 1. Split long edges into two shorter edges [Edge Split]
     * 2. Combine short edges into one long edge [Edge Collapse]
     * 3. For every vertex, minimize its degree�s deviation from 6 (or 4 for boundary) [Edge Flip]
     * 4. Tangentially move every vertex towards its centroid (of all its neighbours) [Laplacian vertex smoothing]
     *
     * This function generates a mesh version for every iteration that you can switch in between to examine.
     * 
     * @param numIterations The number of iterations of remeshing to run
     */
    void performRemeshing(int numIterations);

private:
    void initializeVertexColors();

    void computeVertexNormals();
    void computeVertexNormal(int vertexIndex);

    void computeTriangleNormal(int triIndex);
    void computeTriangleNormals();
    void computeFNext();
    void computeNumComponents();
    void computeStat();

    /**
     * Simple Laplacian Smoothing. ref: https://en.wikipedia.org/wiki/Laplacian_smoothing
     *
     * @param vertexIndex index of the vertex to smooth out
     */
    vec3 computeVertexCentroid(int vertexIndex);

    /* Debug functions */
    void printStats();
    void printVertexList();
    void printfnList();
    void printTriList();
    void printOrTri(OrTri ot);

    void setVertexColor(int vIdx, float r, float g, float b);
};

enum class PropertyType : unsigned char {
    POSITION_X,
    POSITION_Y,
    POSITION_Z,
    NORMAL_X,
    NORMAL_Y,
    NORMAL_Z,
    TEX_S,
    TEX_T
};

