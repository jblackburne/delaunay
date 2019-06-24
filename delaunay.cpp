#include <iostream>
#include <algorithm>

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
delaunay::Triangle<T>::Triangle(delaunay::Point2D<T> const &p1,
				delaunay::Point2D<T> const &p2,
				delaunay::Point2D<T> const &p3)
  : m_corners{p1, p2, p3},
    m_daughters{nullptr, nullptr, nullptr},
    m_neighbors{nullptr, nullptr, nullptr}
{
  // Need to ensure that the points are in general position and ccw order
  T crossp = crossprod(m_corners[1] - m_corners[0], m_corners[2] - m_corners[1]);
  if (crossp == 0) {
    throw std::runtime_error("Cannot form a triangle; points not in general position");
  }
  if (crossp < 0) {
    std::swap(m_corners[0], m_corners[1]);
  }  
}

template <typename T>
void delaunay::Triangle<T>::print() const
{
  std::cout << "Point 1: " << m_corners[0].x << " " << m_corners[0].y << "\n";
  std::cout << "Point 2: " << m_corners[1].x << " " << m_corners[1].y << "\n";
  std::cout << "Point 3: " << m_corners[2].x << " " << m_corners[2].y << "\n";
}

template <typename T>
int delaunay::Triangle<T>::containsPoint(delaunay::Point2D<T> const &p) const
{
  T crossp[3] = {crossprod(m_corners[1] - m_corners[0], p - m_corners[0]),
		 crossprod(m_corners[2] - m_corners[1], p - m_corners[1]),
		 crossprod(m_corners[0] - m_corners[2], p - m_corners[2])};
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
  delaunay::Point2D<double> p1{1, 2};
  delaunay::Point2D<double> p2{2, -1};
  delaunay::Point2D<double> p3{-2, 1};
  delaunay::Triangle<double> triangle(p1, p2, p3);

  triangle.print();

  std::cout << triangle.containsPoint({0, 1}) << "\n";
  
  return 0;
}
