#include <iostream>
#include <algorithm>
#include <stack>
#include <cmath>

#include "delaunay.hpp"


namespace dl = delaunay;

namespace delaunay
{
  template <typename T>
  Point2D<T> operator-(Point2D<T> const &lhs, Point2D<T> const &rhs)
  {
    return Point2D<T>{lhs.x - rhs.x, lhs.y - rhs.y};
  }

  template <typename T>
  Point2D<T> operator+(Point2D<T> const &lhs, Point2D<T> const &rhs)
  {
    return Point2D<T>{lhs.x + rhs.x, lhs.y + rhs.y};
  }

  template <typename T>
  T crossprod(Point2D<T> const &lhs, Point2D<T> const &rhs)
  {
    return lhs.x * rhs.y - lhs.y * rhs.x;
  }
}

template <typename T>
dl::Triangle<T>::Triangle(dl::Point2D<T> *p1,
                          dl::Point2D<T> *p2,
                          dl::Point2D<T> *p3)
  : corners{p1, p2, p3},
    daughters{nullptr, nullptr, nullptr},
    neighbors{nullptr, nullptr, nullptr}
{
  // I'm going to skip the normal nullptr check because this is not user-facing

  // Need to ensure that the points are in general position and ccw order
  T crossp = crossprod(*corners[1] - *corners[0], *corners[2] - *corners[1]);
  if (crossp == 0) {
    throw std::runtime_error("Cannot form a triangle; points not in general position");
  }
  if (crossp < 0) {
    std::swap(corners[0], corners[1]);
  }
}

template <typename T>
void dl::Triangle<T>::print() const
{
  std::cout << "Point 1: " << corners[0]->x << " " << corners[0]->y << "\n";
  std::cout << "Point 2: " << corners[1]->x << " " << corners[1]->y << "\n";
  std::cout << "Point 3: " << corners[2]->x << " " << corners[2]->y << "\n";
}

template <typename T>
int dl::Triangle<T>::containsPoint(dl::Point2D<T> const &p) const
{
  T crossp[3] = {crossprod(*corners[1] - *corners[0], p - *corners[0]),
                 crossprod(*corners[2] - *corners[1], p - *corners[1]),
                 crossprod(*corners[0] - *corners[2], p - *corners[2])};
  if (crossp[0] < 0 || crossp[1] < 0 || crossp[2] < 0) {
    return 0;  // outside
  } else if (crossp[0] == 0) {
    return -1;  // on the edge between corners 0 and 1
  } else if (crossp[1] == 0) {
    return -2;  // on the edge between corners 1 and 2
  } else if (crossp[2] == 0) {
    return -3;  // on the edge between corners 2 and 0
  } else {
    return 1;  // inside
  }
}

template <typename T>
dl::Triangulation<T>::Triangulation(dl::Point2D<T> const *points, size_t nPoints)
{
  // Error checking
  if (nPoints < 3) {
    throw std::runtime_error("Cannot create a Triangulation out of fewer than three points.");
  }
  if (!points) {
    throw std::runtime_error("Null pointer passed to Triangulation.");
  }

  // Keep a copy of the input set of points
  m_points.assign(points, points + nPoints);

  // Figure out the bounding box of the input points
  typename std::vector< Point2D<T> >::const_iterator minx, maxx, miny, maxy;
  std::tie(minx, maxx) = std::minmax_element(m_points.begin(), m_points.end(),
                                             [](dl::Point2D<T> const &a,
                                                dl::Point2D<T> const &b) {return a.x < b.x;});
  std::tie(miny, maxy) = std::minmax_element(m_points.begin(), m_points.end(),
                                             [](dl::Point2D<T> const &a,
                                                dl::Point2D<T> const &b) {return a.y < b.y;});

  // Make a root triangle much bigger than the set of points
  dl::Point2D<T> ctr({0.5 * (minx->x + maxx->x), 0.5 * (miny->y + maxy->y)});
  T radius = 10000 * std::max(maxx->x - minx->x, maxy->y - miny->y);
  m_virtumals[0] = ctr + dl::Point2D<T>{0, radius};
  m_virtumals[1] = ctr + dl::Point2D<T>{-0.5 * radius, -sqrt(3)/2 * radius};
  m_virtumals[2] = ctr + dl::Point2D<T>{0.5 * radius, -sqrt(3)/2 * radius};
  m_triangles.push_back(dl::Triangle<T>(&m_virtumals[0], &m_virtumals[1], &m_virtumals[2]));

  // Now add points to the triangulation
  for (auto &p: m_points) {
    dl::Triangle<T> *mother = findTriangle(p);
    m_triangles.push_back(dl::Triangle<T>(&p, mother->corners[0], mother->corners[1]));
    m_triangles.push_back(dl::Triangle<T>(&p, mother->corners[1], mother->corners[2]));
    m_triangles.push_back(dl::Triangle<T>(&p, mother->corners[2], mother->corners[0]));
    dl::Triangle<T> *daughters[3] = {&m_triangles[m_triangles.size() - 3],
                                     &m_triangles[m_triangles.size() - 2],
                                     &m_triangles[m_triangles.size() - 1]};
    std::copy(daughters, daughters + 3, mother->daughters);
    daughters[0]->neighbors[0] = mother->neighbors[2];
    daughters[1]->neighbors[0] = mother->neighbors[0];
    daughters[2]->neighbors[0] = mother->neighbors[1];
    daughters[0]->neighbors[1] = daughters[1];
    daughters[1]->neighbors[1] = daughters[2];
    daughters[2]->neighbors[1] = daughters[0];
    daughters[0]->neighbors[2] = daughters[2];
    daughters[1]->neighbors[2] = daughters[0];
    daughters[2]->neighbors[2] = daughters[1];
  }
}

template <typename T>
dl::Triangle<T> *dl::Triangulation<T>::findTriangle(dl::Point2D<T> const &point)
{
  std::stack<dl::Triangle<T> *> tristack;
  tristack.push(&m_triangles[0]);
  while (!tristack.empty()) {
    dl::Triangle<T> *tri = tristack.top();
    tristack.pop();
    if (tri->containsPoint(point) != 0) {
      bool isLeaf = true;
      for (size_t i=0; i<3; ++i) {
        if (tri->daughters[i]) {
          tristack.push(tri->daughters[i]);
          isLeaf = false;
        }
      }
      if (isLeaf) {
        return tri;
      }
    }
  }

  // No enclosing "active" (i.e., leaf node) triangle was found
  return nullptr;
}

int main(void)
{
  std::vector< dl::Point2D<double> > p{{1, 2}, {2, -1}, {-2, 1}};
  dl::Triangle<double> triangle(&p[0], &p[1], &p[2]);
  dl::Triangulation<double>(&p[0], p.size());

  triangle.print();

  std::cout << triangle.containsPoint({0, 1}) << "\n";

  return 0;
}
