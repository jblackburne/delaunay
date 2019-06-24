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
  class Triangle
  {
  public:
    Triangle(Point2D<T> const &p1,
	     Point2D<T> const &p2,
	     Point2D<T> const &p3);

    void print() const;
    
    int containsPoint(Point2D<T> const &p) const;
    
  private:
    Point2D<T> m_corners[3];

    Triangle<T> *m_daughters[3];
    Triangle<T> *m_neighbors[3];
  };

  template <typename T>
  class Delaunay
  {
  private:
    std::vector< Point2D<T> > m_points;
    std::vector< Triangle<T> > m_triangles;
  };
}

#endif  // __DELAUNAY_H__
