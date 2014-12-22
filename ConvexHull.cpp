// Implementations of various convex hull algorithms 
// using the C++ Standard Library.
// For clarity, the implementations do not check for
// duplicate or collinear points.

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

struct point {
	float x;
	float y;

	point(float xIn, float yIn) : x(xIn), y(yIn) { } 
};

// The z-value of the cross product of segments 
// (a, b) and (a, c). Positive means c is ccw
// from (a, b), negative cw. Zero means its collinear.
float ccw(const point& a, const point& b, const point& c) {
	return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// Returns true if a is lexicographically before b.
bool isLeftOf(const point& a, const point& b) {
	return (a.x < b.x || (a.x == b.x && a.y < b.y));
}

// Used to sort points in ccw order about a pivot.
struct ccwSorter {
	const point& pivot;

	ccwSorter(const point& inPivot) : pivot(inPivot) { }

	bool operator()(const point& a, const point& b) {
		return ccw(pivot, a, b) < 0;
	}
};

// The length of segment (a, b).
float len(const point& a, const point& b) {
	return sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
}

// The unsigned distance of p from segment (a, b).
float dist(const point& a, const point& b, const point& p) {
	return fabs((b.x - a.x) * (a.y - p.y) - (b.y - a.y) * (a.x - p.x)) / len(a, b);
}

// Returns the index of the farthest point from segment (a, b).
size_t getFarthest(const point& a, const point& b, const vector<point>& v) {
	size_t idxMax = 0;
	float distMax = dist(a, b, v[idxMax]);

	for (size_t i = 1; i < v.size(); ++i) {
		float distCurr = dist(a, b, v[i]);
		if (distCurr > distMax) {
			idxMax = i;
			distMax = distCurr;
		}
	}

	return idxMax;
}


// The gift-wrapping algorithm for convex hull.
// https://en.wikipedia.org/wiki/Gift_wrapping_algorithm
vector<point> giftWrapping(const vector<point>& v) {
	// Start with the leftmost point.
	size_t startIdx = min_element(v.begin(), v.end(), isLeftOf) - v.begin();
	size_t hIdx = startIdx;

	vector<point> hull;
	do {
		// Add our current point to the hull.
		hull.push_back(v[hIdx]);

		// Find the index of the input point that is farthest
		// ccw from our last hull point.
		ccwSorter isCcw(v[hIdx]);
		size_t endIdx = 0;
		for (size_t i = 1; i < v.size(); ++i) {
			if ((endIdx == hIdx) || isCcw(v[endIdx], v[i])) {
				endIdx = i;
			}
		}

		// Make that our new 
		hIdx = endIdx;
	} while (hIdx != startIdx);

	return hull;
}


// The Graham scan algorithm for convex hull.
// https://en.wikipedia.org/wiki/Graham_scan
vector<point> GrahamScan(vector<point> v) {
	// Put our leftmost point at index 0
	swap(v[0], *min_element(v.begin(), v.end(), isLeftOf));

	// Sort the rest of the points in counter-clockwise order
	// from our leftmost point.
	sort(v.begin() + 1, v.end(), ccwSorter(v[0]));
	
	// Add our first three points to the hull.
	vector<point> hull;
	auto it = v.begin();
	hull.push_back(*it++);
	hull.push_back(*it++);
	hull.push_back(*it++);
	
	while (it != v.end()) {
		// Pop off any points that make a convex angle with *it
		while (ccw(*(hull.rbegin() + 1), *(hull.rbegin()), *it) >= 0) {
			hull.pop_back();
		}
		hull.push_back(*it++);
	}

	return hull;
}


// The monotone chain algorithm for convex hull.
vector<point> monotoneChain(vector<point> v) {
	// Sort our points in lexicographic order.
	sort(v.begin(), v.end(), isLeftOf);
	
	// Find the lower half of the convex hull.
	vector<point> lower;
	for (auto it = v.begin(); it != v.end(); ++it) {
		// Pop off any points that make a convex angle with *it
		while (lower.size() >= 2 && ccw(*(lower.rbegin() + 1), *(lower.rbegin()), *it) >= 0) {
			lower.pop_back();
		}
		lower.push_back(*it);
	}
		
	// Find the upper half of the convex hull.
	vector<point> upper;
	for (auto it = v.rbegin(); it != v.rend(); ++it) {
		// Pop off any points that make a convex angle with *it
		while (upper.size() >= 2 && ccw(*(upper.rbegin() + 1), *(upper.rbegin()), *it) >= 0) {
			upper.pop_back();
		}
		upper.push_back(*it);
	}

	vector<point> hull;
	hull.insert(hull.end(), lower.begin(), lower.end());
	// Both hulls include both endpoints, so leave them out when we 
	// append the upper hull.
	hull.insert(hull.end(), upper.begin() + 1, upper.end() - 1);
	return hull;
}


// Recursive call of the quickhull algorithm.
void quickHull(const vector<point>& v, const point& a, const point& b, 
			   vector<point>& hull) {
	if (v.empty()) {
		return;
	}

	point f = v[getFarthest(a, b, v)];

	// Collect points to the left of segment (a, f)
	vector<point> left;
	for (auto p : v) {
		if (ccw(a, f, p) > 0) {
			left.push_back(p);
		}
	}
	quickHull(left, a, f, hull);
	
	// Add f to the hull
	hull.push_back(f);

	// Collect points to the left of segment (f, b)
	vector<point> right;
	for (auto p : v) {
		if (ccw(f, b, p) > 0) {
			right.push_back(p);
		}
	}
	quickHull(right, f, b, hull);
}

// QuickHull algorithm. 
// https://en.wikipedia.org/wiki/QuickHull
vector<point> quickHull(const vector<point>& v) {
	vector<point> hull;
	
	// Start with the leftmost and rightmost points.
	point a = *min_element(v.begin(), v.end(), isLeftOf);
	point b = *max_element(v.begin(), v.end(), isLeftOf);

	// Split the points on either side of segment (a, b)
	vector<point> left, right;
	for (auto p : v) {
		ccw(a, b, p) > 0 ? left.push_back(p) : right.push_back(p);
	}
	
	// Be careful to add points to the hull
	// in the correct order. Add our leftmost point.
	hull.push_back(a);

	// Add hull points from the left (top)
	quickHull(left, a, b, hull);

	// Add our rightmost point
	hull.push_back(b);

	// Add hull points from the right (bottom)
	quickHull(right, b, a, hull);

	return hull;
}

vector<point> getPoints() {
	vector<point> v;
	
	const float lo = -100.0;
	const float hi = 100.0;

	for (int i = 0; i < 100; ++i) {
		float x = lo + 
			static_cast<float>(
				rand()) / static_cast<float>(RAND_MAX / (hi - lo));

		float y = lo + 
			static_cast<float>(
				rand()) / static_cast<float>(RAND_MAX / (hi - lo));

		v.push_back(point(x, y));
	}

	return v;
}

void print(const vector<point>& v) {
	for (auto p : v) {
		cout << p.x << ", " << p.y << endl;
	}
}

int main() { 	
	vector<point> v = getPoints();
	
	vector<point> h = quickHull(v);
	cout << "quickHull point count: " << h.size() << endl;
	print(h);

 	h = giftWrapping(v);
	cout << endl << "giftWrapping point count: " << h.size() << endl;
	print(h);

	h = monotoneChain(v);
	cout << endl << "monotoneChain point count: " << h.size() << endl;
	print(h);
	
	h = GrahamScan(v);
	cout << endl << "GrahamScan point count: " << h.size() << endl;
	print(h);

	return 0;
}
