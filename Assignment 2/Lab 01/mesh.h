#pragma once

#include <unordered_map>
#include <unordered_set>
#include "vec3.h"

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

class myObjType {

    int vcount = 0;
    int tcount = 0;

    double vlist[MAXV][3];   // vertices list
    int tlist[MAXT][3];      // triangle list
    int fnlist[MAXT][6];     // fnext list
    double nlist[MAXT][3];   // storing triangle normals
    double vnlist[MAXT][3];  // Vertex Normal vectors list
    float colorlist[MAXT][4];   // colors list
    int triComponentNumber[MAXT];      // triangle list

    double lmax[3];          // the maximum coordinates of x,y,z
    double lmin[3];          // the minimum coordinates of x,y,z

    int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
    int statMaxAngle[18];

    int numComponents;

    /* For Addition/Removal of Tris and Vertices */
    unordered_set<int> unassignedTris;
    unordered_set<int> unassignedVerts;

    /* todo: For remeshing */
    int vertexDegreeList[MAXV];     // keeps track of the degree of every vertex. 
                                    // if degree is -1, then vertex does not exist.
    
    std::unordered_map<int, std::vector<int> > vertexToTriangles;

public:
    bool smooth = false;

    myObjType() { vcount = 0; tcount = 0; };
    void readFile(char* filename);  // assumming file contains a manifold

    void computeVertexNormals();

    void computeTriangleNormals();

    void initializeVertexColors();
    void computeFNext();
    void computeNumComponents();

    void writeFile(char* filename);
    void draw();
    void computeStat();

    int org(OrTri ot);
    int dest(OrTri ot);
    int last(OrTri ot);

    //// Triangles
    //int AddTriangle(int vertices[3]); // returns the index of the triangle added
    //bool AddTriangleAtIndex(int index, int vertices[3]); // returns true if successful
    //void RemoveTriangleAtIndex(int index);

    //// Vertices
    //void AddVertexAtIndex(int index, vec3 position);
    //void RemoveVertex(int index);

private:
    void printVertexList();
    void printfnList();
    void printTriList();

    void printOrTri(OrTri ot);

    /**
     * Perform Area-Equalizing Remeshing.
     * Based on the paper: “A Remeshing Approach to Multiresolution Modeling” by Botsch & Kobbelt
     * 4-step approach:
     * 1. Split long edges into two shorter edges
     * 2. Combine short edges into one long edge
     * 3. For every vertex, minimize its degree’s deviation from 6 (or 4 for boundary)
     * 4. Tangentially move every vertex towards its centroid (of all its neighbours)
     *
     * @param numIterations The number of iterations of remeshing to run
     */
    void performRemeshing(int numIterations);

    /**
     * Step 1 of the Remeshing algorithm.
     *
     * @param threshold the threshold length of an edge. Edges above this length will be split.
     */
    void splitAllLongEdges(float threshold);

    /**
      * Step 2 of the Remeshing algorithm.
      *
      * @param threshold the threshold length of an edge. Edges under this length will be merged.
      */
    void mergeAllShortEdges(float threshold);

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

