#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <chrono>

using namespace std;

#define PI 3.1415926535

struct Point {
    double latitude;
    double longitude;
};

struct Edge {
    Point start;
    Point end;
    string name;
    string length;
};

unordered_map<long long int, Point> readNodes(string& path) {
    unordered_map<long long int, Point> nodes;
    ifstream file(path);
    string line, word;
    while (getline(file, line)) {
        stringstream ss(line);
        long long int id;
        Point point;
        getline(ss, word, ',');
        id = stoll(word);
        getline(ss, word, ',');
        point.latitude = stod(word);
        getline(ss, word, ',');
        point.longitude = stod(word);
        nodes[id] = point;
    }
    return nodes;
}

vector<Edge> readEdges(string& path, unordered_map<long long int, Point>& nodes) {
    vector<Edge> edges;
    ifstream file(path);
    string line, word;
    while (getline(file, line)) {
        stringstream ss(line);
        Edge edge;
        long long int startId, endId;
        getline(ss, word, ',');
        startId = stoll(word);
        getline(ss, word, ',');
        endId = stoll(word);
        getline(ss, word, ',');
        edge.name = word;
        getline(ss, word, ',');
        edge.length = word;
        edge.start = nodes.at(startId);
        edge.end = nodes.at(endId);
        edges.push_back(edge);
    }
    return edges;
}

void printEdges(vector<Edge>& edges) {
    cout << edges.size() << endl;
    for (auto& edge : edges) {
        cout << "Arista: (" << edge.start.latitude << ", " << edge.start.longitude << ") - ("<< edge.end.latitude << ", " << edge.end.longitude << "), Calle: " << edge.name << endl;
    }
}

Point Centroid(vector<Point>& points) {
    Point centroid = { 0.0, 0.0 };
    for (auto& point : points) {
        centroid.latitude += point.latitude;
        centroid.longitude += point.longitude;
    }
    centroid.latitude /= points.size();
    centroid.longitude /= points.size();
    return centroid;
}

vector<Point> subtractCentroid(vector<Point>& points, Point& centroid) {
    vector<Point> centerPoints;
    centerPoints.reserve(points.size());
    for (auto& point : points) {
        centerPoints.push_back({ point.latitude - centroid.latitude, point.longitude - centroid.longitude });
    }
    return centerPoints;
}

vector<vector<double>> covMatrix(vector<Point>& points) {
    vector<vector<double>> M(2, vector<double>(2, 0.0));
    for (auto& point : points) {
        M[0][0] += point.latitude * point.latitude;
        M[0][1] += point.latitude * point.longitude;
        M[1][0] += point.longitude * point.latitude;
        M[1][1] += point.longitude * point.longitude;
    }
    M[0][0] /= points.size();
    M[0][1] /= points.size();
    M[1][0] /= points.size();
    M[1][1] /= points.size();
    return M;
}

Point powerIteration(vector<Point>& centeredPoints, int iterations) {
    Point b = { 1.0, 0.0 };
    vector<vector<double>> cMatrix = covMatrix(centeredPoints);
    for (int iter = 0; iter < iterations; iter++) {
        Point Ab;
        Ab.latitude = cMatrix[0][0] * b.latitude + cMatrix[0][1] * b.longitude;
        Ab.longitude = cMatrix[1][0] * b.latitude + cMatrix[1][1] * b.longitude;
        double norm = sqrt(Ab.latitude * Ab.latitude + Ab.longitude * Ab.longitude);
        b.latitude = Ab.latitude / norm;
        b.longitude = Ab.longitude / norm;
    }
    return b;
}

Point PCA(vector<Point>& points, int iterations = 10) {
    Point centroid = Centroid(points);
    vector<Point> centeredPoints = subtractCentroid(points, centroid);
    return powerIteration(centeredPoints, iterations);
}

struct BallTreeNode {
    Point center;
    double radius;
    vector<Edge> edges;
    BallTreeNode* left;
    BallTreeNode* right;

    BallTreeNode() : left(nullptr), right(nullptr) {}
};

class BallTree {
private:
    BallTreeNode* root;

    BallTreeNode* buildTree(vector<Edge>& edges, int depth) {
        if (edges.size() <= 1) return nullptr;
        BallTreeNode* node = new BallTreeNode();
        vector<Point> points;
        for (auto& edge : edges) {
            points.push_back(edge.start);
            points.push_back(edge.end);
        }
        Point centroid = Centroid(points);
        node->center = centroid;
        node->radius = 0;
        for (auto& point : points) {
            double dist = geodesicDistance(point, node->center);
            if (dist > node->radius) {
                node->radius = dist;
            }
        }
        Point pc = PCA(points);
        vector<pair<double, Edge>> projections;
        for (auto& edge : edges) {
            double projStart = (edge.start.latitude - centroid.latitude) * pc.latitude + (edge.start.longitude - centroid.longitude) * pc.longitude;
            double projEnd = (edge.end.latitude - centroid.latitude) * pc.latitude + (edge.end.longitude - centroid.longitude) * pc.longitude;
            double projection = (projStart + projEnd) / 2.0;
            projections.emplace_back(projection, edge);
        }
        sort(projections.begin(), projections.end(), [](auto& a, auto& b) { return a.first < b.first; });
        size_t mid = projections.size() / 2;
        vector<Edge> leftEdges, rightEdges;
        for (size_t i = 0; i < mid; ++i) {
            leftEdges.push_back(projections[i].second);
        }
        for (size_t i = mid; i < projections.size(); ++i) {
            rightEdges.push_back(projections[i].second);
        }
        if (leftEdges.size() == edges.size() || rightEdges.size() == edges.size()) {
            node->edges = edges;
            node->left = nullptr;
            node->right = nullptr;
        }
        else {
            node->edges = edges;
            node->left = buildTree(leftEdges, depth + 1);
            node->right = buildTree(rightEdges, depth + 1);
        }
        return node;
    }

    double EdgeDistance(Point& p, Point& v, Point& w) {
        double l2 = pow(geodesicDistance(v, w), 2);
        if (l2 == 0.0) return geodesicDistance(p, v);
        double t = ((p.latitude - v.latitude) * (w.latitude - v.latitude) + (p.longitude - v.longitude) * (w.longitude - v.longitude)) / l2;
        t = max(0.0, min(1.0, t));
        Point projection = { v.latitude + t * (w.latitude - v.latitude), v.longitude + t * (w.longitude - v.longitude) };
        return geodesicDistance(p, projection);
    }

    void nearestEdgeHelper(BallTreeNode* node, Point& query, Edge& bestEdge, double& bestDist) {
        if (!node) return;

        double distToCenter = geodesicDistance(query, node->center);
        if (distToCenter - node->radius > bestDist) {
            return;
        }
        for (auto& edge : node->edges) {
            double dist = EdgeDistance(query, edge.start, edge.end);
            if (dist < bestDist) {
                bestDist = dist;
                bestEdge = edge;
            }
        }
        nearestEdgeHelper(node->left, query, bestEdge, bestDist);
        nearestEdgeHelper(node->right, query, bestEdge, bestDist);
    }

public:
    BallTree(vector<Edge>& edges) {
        root = buildTree(edges, 0);
    }

    Edge nearestEdge(Point& query, double &bestDist) {
        Edge bestEdge;
        bestDist = numeric_limits<double>::max();
        nearestEdgeHelper(root, query, bestEdge, bestDist);
        return bestEdge;
    }

    double geodesicDistance(Point& a, Point& b) {
        double R = 6371e3;
        double phi1 = a.latitude * PI / 180.0;
        double phi2 = b.latitude * PI / 180.0;
        double deltaPhi = (b.latitude - a.latitude) * PI / 180.0;
        double deltaLambda = (b.longitude - a.longitude) * PI / 180.0;
        double A = sin(deltaPhi / 2) * sin(deltaPhi / 2) + cos(phi1) * cos(phi2) *sin(deltaLambda / 2) * sin(deltaLambda / 2);
        double C = 2 * atan2(sqrt(A), sqrt(1 - A));
        return R * C;
    }
};

int main() {
    cout << fixed << setprecision(10);
    string nodespath = "nodes.csv";
    string edgespath = "edges.csv";
    unordered_map<long long int, Point> nodes = readNodes(nodespath);
    vector<Edge> edges = readEdges(edgespath, nodes);
    //printEdges(edges);
    BallTree tree(edges);
    float x, y;
    cout << "Ball* Tree" << endl;
    cout << "Ingrese un punto query:" << endl;
    cout << "Ingrese latitud: ";
    cin >> x;
    cout << "Ingrese longitud: ";
    cin >> y;
    Point query = { x, y };
    double dist;
    auto start = chrono::system_clock::now();
    Edge nearest = tree.nearestEdge(query, dist);
    cout << "Arista mas cercana: ";
    cout << "(" << nearest.start.latitude << ", " << nearest.start.longitude << ") - (" << nearest.end.latitude << ", " << nearest.end.longitude << ")" << endl;
    cout << "Calle: " << nearest.name << endl;
    cout << "Longitud de la calle: " << nearest.length << "m" << endl;
    cout << "Distancia: " << dist << "m" << endl;
    auto end = chrono::system_clock::now();
    chrono::duration<float, milli> duration = end - start;
    cout << endl << "Tiempo: " << duration.count() << "ms " << endl;
}
