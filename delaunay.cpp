#include <iostream>
#include <algorithm>
#include <tuple>
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
dl::Triangulation<T>::Triangulation(dl::Point2D<T> const *points, size_t nPoints)
{
  // Error checking
  if (nPoints < 3) {
    throw std::runtime_error("Cannot create a Triangulation out of fewer than three points.");
  }
  if (!points) {
    throw std::runtime_error("Null pointer passed to Triangulation.");
  }

  // Figure out the bounding box of the input points
  // Throw an exception if the points are not in general position
  Point2D<T> const *minx, *maxx, *miny, *maxy;
  std::tie(minx, maxx) = std::minmax_element(points, points + nPoints,
                                             [](dl::Point2D<T> const &a,
                                                dl::Point2D<T> const &b) {return a.x < b.x;});
  std::tie(miny, maxy) = std::minmax_element(points, points + nPoints,
                                             [](dl::Point2D<T> const &a,
                                                dl::Point2D<T> const &b) {return a.y < b.y;});
  if (maxx->x == minx->x || maxy->y == miny->y) {
    throw std::runtime_error("Points are not in general position.");
  }

  // Keep a copy of the input set of points
  m_points.resize(nPoints + 3);  // extra three to hold the root triangle's corners
  std::copy(points, points + nPoints, &m_points[3]);

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

    m_corners.push_back({iPoint, m_corners[iMother][1], m_corners[iMother][2]});
    m_corners.push_back({iPoint, m_corners[iMother][2], m_corners[iMother][0]});
    m_corners.push_back({iPoint, m_corners[iMother][0], m_corners[iMother][1]});

    m_daughters.push_back({-1, -1, -1});
    m_daughters.push_back({-1, -1, -1});
    m_daughters.push_back({-1, -1, -1});

    m_neighbors.push_back({m_neighbors[iMother][0], m_daughters[iMother][1], m_daughters[iMother][2]});
    m_neighbors.push_back({m_neighbors[iMother][1], m_daughters[iMother][2], m_daughters[iMother][0]});
    m_neighbors.push_back({m_neighbors[iMother][2], m_daughters[iMother][0], m_daughters[iMother][1]});

    // Make sure the neighbors know about the new triangles
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        if (m_neighbors[m_neighbors[iMother][i]][j] == iMother) {
          m_neighbors[m_neighbors[iMother][i]][j] = m_daughters[iMother][i];
        }
      }
    }
  }
}

template <typename T>
int dl::Triangulation<T>::findTriangle(dl::Point2D<T> const &point) const
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

  // No enclosing non-degenerate leaf node triangle was found
  return -1;
}

template <typename T>
bool dl::Triangulation<T>::isLeaf(int iTri) const
{
  return (m_daughters[iTri][0] == -1 &&
          m_daughters[iTri][1] == -1 &&
          m_daughters[iTri][2] == -1);
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

template <typename T>
bool dl::Triangulation<T>::needsFlipped(int iMe, int jThem) const
{
  const int iThem = m_neighbors[iMe][jThem];

  // Never flip non-leaf triangles
  if (!isLeaf(iMe) || !isLeaf(iThem)) {
    return false;
  }

  // Do not flip if either triangle touches the root triangle
  if (m_corners[iMe][0] < 3 ||
      m_corners[iMe][1] < 3 ||
      m_corners[iMe][2] < 3 ||
      m_corners[iThem][0] < 3 ||
      m_corners[iThem][1] < 3 ||
      m_corners[iThem][2] < 3) {
    return false;
  }

  // If either triangle is degenerate, flip for sure
  if (isDegenerate(iMe) || isDegenerate(iThem)) {
    return true;
  }

  // If no other contingencies were triggered, use the normal flipping logic
  // Flip if the fourth point falls outside the circumcircle of the other three points
  dl::Point2D<T> ad = m_corners[iThem][0] - m_corners[iMe][jThem];
  dl::Point2D<T> bd = m_corners[iThem][1] - m_corners[iMe][jThem];
  dl::Point2D<T> cd = m_corners[iThem][2] - m_corners[iMe][jThem];
  // matrix = [[ad.x, ad.y, ad.x**2 + ad.y**2],
  //           [bd.x, bd.y, bd.x**2 + bd.y**2],
  //           [cd.x, cd.y, cd.x**2 + cd.y**2]]
  T det = (+(ad.x*ad.x + ad.y*ad.y) * (bd.x * cd.y - cd.x * bd.y)
           -(bd.x*bd.x + bd.y*bd.y) * (ad.x * cd.y - cd.x * ad.y)
           +(cd.x*cd.x + cd.y*cd.y) * (ad.x * bd.y - bd.x * ad.y));

  return det > 0;
}

template <typename T>
void dl::Triangulation<T>::flip(int iMe, int jThem)
{
  // NOTE: "i" indices are into the triangle vectors
  //   and "j" indices are into the triplet arrays
  // Identify which of the neighbor's neighbors I am
  int iThem = m_neighbors[iMe][jThem];
  int jMe = -1;
  for (int j=0; j<3; ++j) {
    if (m_neighbors[iThem][j] == iMe) {
      jMe = j;
    }
  }
  assert(jMe >= 0);  // Will be true unless I have a bug in my programming

  // Create two new triangles, daughters to both "me" and "them"
  m_daughters[iMe][(jThem + 2) % 3] = static_cast<int>(m_daughters.size());
  m_daughters[iThem][(jMe + 1) % 3] = static_cast<int>(m_daughters.size());
  m_daughters[iMe][(jThem + 1) % 3] = static_cast<int>(m_daughters.size() + 1);
  m_daughters[iThem][(jMe + 2) % 3] = static_cast<int>(m_daughters.size() + 1);

  m_corners.push_back({m_corners[iMe][jThem], m_corners[iMe][(jThem + 1) % 3], m_corners[iThem][jMe]});
  m_corners.push_back({m_corners[iThem][jMe], m_corners[iThem][(jMe + 1) % 3], m_corners[iMe][jThem]});

  m_daughters.push_back({-1, -1, -1});
  m_daughters.push_back({-1, -1, -1});

  m_neighbors.push_back({m_neighbors[iThem][(iMe + 1) % 3],
                         m_daughters[iMe][(jThem + 2) % 3],
                         m_neighbors[iMe][(jThem + 2) % 3]});
  m_neighbors.push_back({m_neighbors[iMe][(jThem + 1) % 3],
                         m_daughters[iMe][(jThem + 1) % 3],
                         m_neighbors[iThem][(iMe + 2) % 3]});

  // Let the other neighbors know about the two new daughters
  for (int jNeigh=jThem+1; jNeigh<jThem+3; ++jNeigh) {
    int iNeigh = m_neighbors[iMe][jNeigh % 3];
    for (int j=0; j<3; ++j) {
      if (m_neighbors[iNeigh][j] == iMe) {
        m_neighbors[iNeigh][j] = m_daughters[iMe][jNeigh % 3];
      }
    }
  }
  for (int jNeigh=jMe+1; jNeigh<jMe+3; ++jNeigh) {
    int iNeigh = m_neighbors[iThem][jNeigh % 3];
    for (int j=0; j<3; ++j) {
      if (m_neighbors[iNeigh][j] == iThem) {
        m_neighbors[iNeigh][j] = m_daughters[iThem][jNeigh % 3];
      }
    }
  }
}


int main(void)
{
  std::vector< dl::Point2D<double> > p{{1, 2}, {2, -1}, {-2, 1}};
  dl::Triangulation<double>(&p[0], p.size());

  return 0;
}
