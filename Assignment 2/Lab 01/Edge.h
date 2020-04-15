#pragma once

#include <xhash>
#include <vector>

class Edge {

public:
    int v1, v2;

    /**
     * The only allowed edge constructor, since an edge can ONLY be made out of two vertices.
     */
    Edge(int v1, int v2) : v1(v1), v2(v2) { };

    Edge() = delete;

};

template <>
struct std::hash<Edge> {
    size_t operator()(const Edge &x) const {
        return std::hash<int>{}(x.v1) ^ std::hash<int>{}(x.v2);
    }
};

template <>
struct std::equal_to<Edge> {
    bool operator()(const Edge &x, const Edge &y) const {
        return (x.v1 == y.v1 && x.v2 == y.v2) || (x.v2 == y.v1 && x.v1 == y.v2);
    }
};

template <>
struct std::hash<std::vector<Edge>> {
    size_t operator()(const std::vector<Edge> &vect) const {
        size_t accumulator = 0;
        for (auto &x : vect) accumulator = (accumulator * 16777213 + std::hash<Edge>()(x)) % 2147483647;
        return accumulator;
    }
};