#include <iostream>
#include <algorithm>
#include <stack>
#include <cmath>
#include <cassert>

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
  m_points.resize(nPoints + 3);  // extra three to hold the root triangle's corners
  std::copy(points, points + nPoints, &m_points[2]);

  // Figure out the bounding box of the input points
  Point2D<T> const *minx, *maxx, *miny, *maxy;
  std::tie(minx, maxx) = std::minmax_element(points, points + nPoints,
                                             [](dl::Point2D<T> const &a,
                                                dl::Point2D<T> const &b) {return a.x < b.x;});
  std::tie(miny, maxy) = std::minmax_element(points, points + nPoints,
                                             [](dl::Point2D<T> const &a,
                                                dl::Point2D<T> const &b) {return a.y < b.y;});

  // Make a root triangle much bigger than the set of points
  dl::Point2D<T> ctr({0.5 * (minx->x + maxx->x), 0.5 * (miny->y + maxy->y)});
  T radius = 10000 * std::max(maxx->x - minx->x, maxy->y - miny->y);
  m_points[0] = ctr + dl::Point2D<T>{0, radius};
  m_points[1] = ctr + dl::Point2D<T>{-0.5 * radius, -sqrt(3)/2 * radius};
  m_points[2] = ctr + dl::Point2D<T>{0.5 * radius, -sqrt(3)/2 * radius};
  m_corners.push_back({0, 1, 2});
  m_daughters.push_back({-1, -1, -1});
  m_neighbors.push_back({-1, -1, -1});

  // Now add points to the triangulation
  for (int iPoint = 3; iPoint < m_points.size(); ++iPoint) {
    // Find the triangle that encloses this new point
    int iMother = findTriangle(m_points[iPoint]);
    assert(iMother >= 0);  // This is guaranteed if we constructed the root triangle correctly

    // Construct three new triangles by connecting the new point to the corners of the mother
    m_daughters[iMother][0] = static_cast<int>(m_daughters.size());
    m_daughters[iMother][1] = static_cast<int>(m_daughters.size() + 1);
    m_daughters[iMother][2] = static_cast<int>(m_daughters.size() + 2);

    m_corners.push_back({iPoint, m_corners[iMother][0], m_corners[iMother][1]});
    m_corners.push_back({iPoint, m_corners[iMother][1], m_corners[iMother][2]});
    m_corners.push_back({iPoint, m_corners[iMother][2], m_corners[iMother][0]});

    m_daughters.push_back({-1, -1, -1});
    m_daughters.push_back({-1, -1, -1});
    m_daughters.push_back({-1, -1, -1});

    m_neighbors.push_back({m_neighbors[iMother][2], m_daughters[iMother][1], m_daughters[iMother][2]});
    m_neighbors.push_back({m_neighbors[iMother][0], m_daughters[iMother][2], m_daughters[iMother][0]});
    m_neighbors.push_back({m_neighbors[iMother][1], m_daughters[iMother][0], m_daughters[iMother][1]});
  }
}

template <typename T>
int dl::Triangulation<T>::findTriangle(dl::Point2D<T> const &point)
{
  std::stack<int> tristack;
  tristack.push(0);
  while (!tristack.empty()) {
    int iTri = tristack.top();
    tristack.pop();
    if (containsPoint(iTri, point) != 0) {
      bool isLeaf = true;
      for (size_t i=0; i<3; ++i) {
        if (m_daughters[iTri][i] > -1) {
          tristack.push(m_daughters[iTri][i]);
          isLeaf = false;
        }
      }
      if (isLeaf && !isDegenerate(iTri)) {
        return iTri;
      }
    }
  }

  // No enclosing "active" (i.e., leaf node) triangle was found
  return -1;
}

template <typename T>
bool dl::Triangulation<T>::isDegenerate(int iTri) const
{
  dl::Point2D<T> const *corners[3] = {&m_points[m_corners[iTri][0]],
				      &m_points[m_corners[iTri][1]],
				      &m_points[m_corners[iTri][2]]};
  T crossp = crossprod(*corners[1] - *corners[0], *corners[2] - *corners[1]);

  return (crossp <= 0);
}

template <typename T>
int dl::Triangulation<T>::containsPoint(int iTri, dl::Point2D<T> const &p) const
{
  dl::Point2D<T> const *corners[3] = {&m_points[m_corners[iTri][0]],
				      &m_points[m_corners[iTri][1]],
				      &m_points[m_corners[iTri][2]]};
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

int main(void)
{
  std::vector< dl::Point2D<double> > p{{1, 2}, {2, -1}, {-2, 1}};
  dl::Triangle<double> triangle(&p[0], &p[1], &p[2]);
  dl::Triangulation<double>(&p[0], p.size());

  triangle.print();

  return 0;
}
