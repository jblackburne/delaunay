#include <iostream>
#include <algorithm>
#include <cmath>

#include "delaunay.hpp"


namespace delaunay
{
  template <typename T>
  delaunay::Point2D<T> operator-(delaunay::Point2D<T> const &lhs, delaunay::Point2D<T> const &rhs)
  {
    return delaunay::Point2D<T>{lhs.x - rhs.x, lhs.y - rhs.y};
  }

  template <typename T>
  T crossprod(delaunay::Point2D<T> const &lhs, delaunay::Point2D<T> const &rhs)
  {
    return lhs.x * rhs.y - lhs.y * rhs.x;
  }
}

template <typename T>
delaunay::Triangle<T>::Triangle(delaunay::Point2D<T> *p1,
				delaunay::Point2D<T> *p2,
				delaunay::Point2D<T> *p3)
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
void delaunay::Triangle<T>::print() const
{
  std::cout << "Point 1: " << corners[0]->x << " " << corners[0]->y << "\n";
  std::cout << "Point 2: " << corners[1]->x << " " << corners[1]->y << "\n";
  std::cout << "Point 3: " << corners[2]->x << " " << corners[2]->y << "\n";
}

template <typename T>
int delaunay::Triangle<T>::containsPoint(delaunay::Point2D<T> const &p) const
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
delaunay::Triangulation<T>::Triangulation(delaunay::Point2D<T> const *points, size_t nPoints)
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
  T minx, maxx, miny, maxy;
  std::tie(minx, maxx) = std::minmax_element(m_points.begin(), m_points.end(),
					     [](delaunay::Point2D<T> const &a,
						delaunay::Point2D<T> const &b) {return a.x < b.x;});
  std::tie(miny, maxy) = std::minmax_element(m_points.begin(), m_points.end(),
					     [](delaunay::Point2D<T> const &a,
						delaunay::Point2D<T> const &b) {return a.y < b.y;});

  // Make a root triangle much bigger than the set of points
  delaunay::Point2D<T> ctr({0.5 * (minx + maxx), 0.5 * (miny + maxy)});
  T radius = 10000 * std::max(maxx - minx, maxy - miny);
  m_virtumals = {ctr + delaunay::Point2D<T>{0, radius},
		 ctr + delaunay::Point2D<T>{-0.5 * radius, -sqrt(3)/2 * radius},
		 ctr + delaunay::Point2D<T>{0.5 * radius, -sqrt(3)/2 * radius}};
  m_triangles.push_back(delaunay::Triangle<T>(&m_virtumals[0], &m_virtumals[1], &m_virtumals[2]));

  // Now add points to the triangulation
  for (auto const &p: m_points) {
    delaunay::Triangle<T> *mother = findTriangle(p);
    m_triangles.emplace_back(*p, *mother->corners[0], *mother->corners[1]);
    m_triangles.emplace_back(*p, *mother->corners[1], *mother->corners[2]);
    m_triangles.emplace_back(*p, *mother->corners[2], *mother->corners[0]);
    delaunay::Triangle<T> *daughters[3] = {&m_triangles[m_triangles.size() - 3],
					   &m_triangles[m_triangles.size() - 2],
					   &m_triangles[m_triangles.size() - 1]};
    std::copy(daughters, daughters + 3, mother->daughters);
    daughters[0].neighbors[0] = mother->neighbors[2];
    daughters[1].neighbors[0] = mother->neighbors[0];
    daughters[2].neighbors[0] = mother->neighbors[1];
    daughters[0].neighbors[1] = daughters[1];
    daughters[1].neighbors[1] = daughters[2];
    daughters[2].neighbors[1] = daughters[0];
    daughters[0].neighbors[2] = daughters[2];
    daughters[1].neighbors[2] = daughters[0];
    daughters[2].neighbors[2] = daughters[1];
  }
}

int main(void)
{
  delaunay::Point2D<double> p1{1, 2};
  delaunay::Point2D<double> p2{2, -1};
  delaunay::Point2D<double> p3{-2, 1};
  delaunay::Triangle<double> triangle(&p1, &p2, &p3);

  triangle.print();

  std::cout << triangle.containsPoint({0, 1}) << "\n";

  return 0;
}
