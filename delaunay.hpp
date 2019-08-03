#ifndef __DELAUNAY_H__
#define __DELAUNAY_H__

#include <vector>
#include <array>


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
    int findTriangle(Point2D<T> const &point) const;

    // Triangle operations
    bool isLeaf(int iTri) const;
    bool isDegenerate(int iTri) const;
    int containsPoint(int iTri, Point2D<T> const &p) const;

    // Flip-related operations
    bool needsFlipped(int iMe, int jThem) const;
    void flip(int iMe, int jThem);

    // We retain a copy of the input points
    // (plus three extras for the "root" triangle)
    std::vector< Point2D<T> > m_points;

    // Arrays to hold the "triangles"
    // A value of -1 indicates no data
    std::vector< std::array<int, 3> > m_corners;  // indices into m_points
    std::vector< std::array<int, 3> > m_daughters;  // indices into these triangle vectors
    std::vector< std::array<int, 3> > m_neighbors;  // indices into these triangle vectors
  };
}

#endif  // __DELAUNAY_H__
