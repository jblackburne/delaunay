#ifndef __DELAUNAY_H__
#define __DELAUNAY_H__

#include <vector>


namespace delaunay
{
  template <typename T>
  struct Point2D
  {
    T x;
    T y;
  };

  template <typename T>
  struct Triangle
  {
    Triangle(Point2D<T> *p1,
	     Point2D<T> *p2,
	     Point2D<T> *p3);

    void print() const;

    int containsPoint(Point2D<T> const &p) const;

    Point2D<T> *corners[3];
    Triangle<T> *daughters[3];
    Triangle<T> *neighbors[3];
  };

  template <typename T>
  class Triangulation
  {
  public:
    Triangulation(Point2D<T> const *points, size_t nPoints);

  private:
    Triangle<T> *findTriangle(Point2D<T> const &point);

    std::vector< Point2D<T> > m_points;
    std::vector< Triangle<T> > m_triangles;
    delaunay::Point2D<T> m_virtumals[3];
  };
}

#endif  // __DELAUNAY_H__
