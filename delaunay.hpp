#ifndef __DELAUNAY_H__
#define __DELAUNAY_H__

#include <vector>
#include <array>
#include <stack>


namespace delaunay
{
  template <typename T>
  struct Point2D
  {
    T x;
    T y;
  };

  template <typename T>
  class Triangulation
  {
  public:
    Triangulation(Point2D<T> const *points, size_t nPoints);

    void print() const;

  private:
    int findTriangle(Point2D<T> const &point) const;

    // Triangle operations
    bool isLeaf(int iTri) const;
    bool isDegenerate(int iTri) const;
    int containsPoint(int iTri, Point2D<T> const &p) const;

    // Flip-related operations
    bool needsFlipped(int iMe, int jThem) const;
    void flip(int iMe, int jThem, std::stack< std::pair<int, int> > &flipStack);

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
