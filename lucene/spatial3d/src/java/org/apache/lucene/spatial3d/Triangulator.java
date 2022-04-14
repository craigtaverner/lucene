package org.apache.lucene.spatial3d;

import org.apache.lucene.geo.GeoUtils;
import org.apache.lucene.geo.XYEncodingUtils;
import org.apache.lucene.spatial3d.geom.*;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

import static org.apache.lucene.geo.GeoUtils.WindingOrder;

public class Triangulator {

  /**
   * state of the triangulation split - avoids recursion
   */
  private enum State {
    INIT,
    CURE,
    SPLIT
  }

  // No Instance:
  private Triangulator() {
  }

  public static List<Triangle> triangulate(final GeoPolygon polygon, boolean checkSelfIntersections) {
    if (polygon == null) return null;
    // TODO: deal with holes
    return triangulate(polygon.getOuterShellPoints(), checkSelfIntersections);
  }

  public static List<Triangle> triangulate(final List<GeoPoint> points, boolean checkSelfIntersections) {
    // Attempt to establish a doubly-linked list of the provided shell points (should be CCW, but this will correct);
    // then filter instances of intersections.
    Node outerNode = createDoublyLinkedList(points);
    // If an outer node hasn't been detected, the shape is malformed. (must comply with OGC SFA specification)
    if (outerNode == null) {
      throw new IllegalArgumentException("Malformed shape detected in Tessellator!");
    }
    if (outerNode == outerNode.next || outerNode == outerNode.next.next) {
      throw new IllegalArgumentException("at least three non-collinear points required");
    }

    // Determine if the specified list of points contains holes
/*
    if (polygon.numHoles() > 0) {
      // Eliminate the hole triangulation.
      outerNode = eliminateHoles(polygon, outerNode);
    }
*/

    if (checkSelfIntersections) {
      checkIntersection(outerNode);
    }
    // Calculate the triangulation using the doubly LinkedList.
    List<Triangle> result = earcutLinkedList(outerNode, new ArrayList<>());
    if (result.size() == 0) {
      throw new IllegalArgumentException("Unable to Triangulate shape. Possible malformed shape detected.");
    }

    return result;
  }

  // TODO: Calculate this within the GeoPolygon and using xyz if possible
  // Right now this is just a copy of the code in the original Polygon constructor
  public static GeoUtils.WindingOrder calculateWindingOrder(List<GeoPoint> points) {
    GeoPoint firstPoint = points.get(0);
    GeoPoint lastPoint = points.get(points.size() - 1);
    double minLat = firstPoint.getLatitude();
    double maxLat = minLat;
    double minLon = firstPoint.getLongitude();
    double maxLon = minLon;
    double windingSum = 0d;
    double lastlat = lastPoint.getLatitude();
    double lastlon = lastPoint.getLongitude();
    double prevlat = firstPoint.getLatitude();
    double prevlon = firstPoint.getLongitude();
    for (GeoPoint point : points) {
      double lat = point.getLatitude();
      double lon = point.getLongitude();
      if (point != firstPoint) {
        minLat = Math.min(lat, minLat);
        maxLat = Math.max(lat, maxLat);
        minLon = Math.min(lon, minLon);
        maxLon = Math.max(lon, maxLon);
        // compute signed area
        windingSum += (prevlon - lastlon) * (lat - lastlat) - (prevlat - lastlat) * (lon - lastlon);
      }
      prevlat = lat;
      prevlon = lon;
    }
    return (windingSum >= 0) ? GeoUtils.WindingOrder.CCW : GeoUtils.WindingOrder.CW;
  }

  /**
   * Creates a circular doubly linked list using polygon points. The order is governed by the
   * specified winding order
   */
  private static Node createDoublyLinkedList(final List<GeoPoint> points) {
    WindingOrder polyWindingOrder = calculateWindingOrder(points);
    if (polyWindingOrder != WindingOrder.CCW) {
      throw new IllegalArgumentException("Cannot triangulate outer polygon not in counter-clockwise winding order");
    }
    int startIndex = 0;
    Node lastNode = null;
    // Link points into the circular doubly-linked list in the default CCW winding order
    for (GeoPoint point : points) {
      lastNode = insertNode(point, startIndex++, lastNode);
    }
    // if first and last node are the same then remove the end node and set lastNode to the start
    if (lastNode != null && isVertexEquals(lastNode, lastNode.next)) {
      removeNode(lastNode, true);
      lastNode = lastNode.next;
    }

    // Return the last node in the Doubly-Linked List
    return filterPoints(lastNode, null);
  }

  /**
   * Eliminate colinear/duplicate points from the doubly linked list
   */
  private static Node filterPoints(final Node start, Node end) {
    if (start == null) {
      return start;
    }

    if (end == null) {
      end = start;
    }

    Node node = start;
    Node nextNode;
    Node prevNode;
    boolean continueIteration;

    do {
      continueIteration = false;
      nextNode = node.next;
      prevNode = node.previous;
      // we can filter points when:
      // 1. they are the same
      // 2.- each one starts and ends in each other
      // 3.- they are co-linear and both edges have the same value in .isNextEdgeFromPolygon
      if (isVertexEquals(node, nextNode)
          || isVertexEquals(prevNode, nextNode)
          || (prevNode.isNextEdgeFromPolygon == node.isNextEdgeFromPolygon
          && colinear(prevNode.point, node.point, nextNode.point))) {
        // Remove the node
        removeNode(node, prevNode.isNextEdgeFromPolygon);
        node = end = prevNode;

        if (node == nextNode) {
          break;
        }
        continueIteration = true;
      } else {
        node = nextNode;
      }
    } while (continueIteration || node != end);
    return end;
  }

  /**
   * Main ear slicing loop which triangulates the vertices of a polygon, provided as a doubly-linked
   * list. *
   */
  private static List<Triangle> earcutLinkedList(Node currEar, final List<Triangle> triangles) {
    State state = State.INIT;
    earcut:
    while (currEar != null && currEar.previous != currEar.next) {

      Node stop = currEar;
      Node prevNode;
      Node nextNode;

      // Iteratively slice ears
      do {
        prevNode = currEar.previous;
        nextNode = currEar.next;
        // Determine whether the current triangle must be cut off.
        final boolean isConvex = isRightHandRule(prevNode.point, currEar.point, nextNode.point);
        if (isConvex && isEar(currEar)) {
          // Compute if edges belong to the polygon
          boolean abFromPolygon = prevNode.isNextEdgeFromPolygon;
          boolean bcFromPolygon = currEar.isNextEdgeFromPolygon;
          boolean caFromPolygon = isEdgeFromPolygon(prevNode, nextNode);
          // Return the triangulated data
          triangles.add(
              new Triangle(
                  prevNode, abFromPolygon, currEar, bcFromPolygon, nextNode, caFromPolygon));
          // Remove the ear node.
          removeNode(currEar, caFromPolygon);

          // Skipping to the next node leaves fewer slither triangles.
          currEar = nextNode.next;
          stop = nextNode.next;
          continue;
        }
        currEar = nextNode;
        // If the whole polygon has been iterated over and no more ears can be found.
        if (currEar == stop) {
          switch (state) {
            case INIT:
              // try filtering points and slicing again
              currEar = filterPoints(currEar, null);
              state = State.CURE;
              continue earcut;
            case CURE:
              // if this didn't work, try curing all small self-intersections locally
              currEar = cureLocalIntersections(currEar, triangles);
              state = State.SPLIT;
              continue earcut;
            case SPLIT:
              // as a last resort, try splitting the remaining polygon into two
              if (!splitEarcut(currEar, triangles)) {
                // we could not process all points. Tessellation failed
                throw new IllegalArgumentException("Unable to Tessellate shape. Possible malformed shape detected.");
              }
              break;
          }
          break;
        }
      } while (currEar.previous != currEar.next);
      break;
    }
    // Return the calculated tessellation
    return triangles;
  }

  /**
   * Determines whether a polygon node forms a valid ear with adjacent nodes. *
   */
  private static boolean isEar(final Node ear) {
    // make sure there aren't other points inside the potential ear
    Node node = ear.next.next;
    while (node != ear.previous) {
      if (pointInEar(node.point, ear.previous.point, ear.point, ear.next.point)
          && !isRightHandRule(node.previous.point, node.point, node.next.point)) {
        return false;
      }
      node = node.next;
    }
    return true;
  }

  /**
   * Iterate through all polygon nodes and remove small local self-intersections *
   */
  private static Node cureLocalIntersections(Node startNode, final List<Triangle> triangles) {
    Node node = startNode;
    Node nextNode;
    do {
      nextNode = node.next;
      Node a = node.previous;
      Node b = nextNode.next;

      // a self-intersection where edge (v[i-1],v[i]) intersects (v[i+1],v[i+2])
      if (!isVertexEquals(a, b)
          && edgeNeverIntersectsPolygon(a, a.point, b.point)
          && linesIntersect(a.point, node.point, nextNode.point, b.point)
          && isLocallyInside(a, b)
          && isLocallyInside(b, a)) {
        // compute edges from polygon
        boolean abFromPolygon =
            (a.next == node)
                ? a.isNextEdgeFromPolygon
                : isEdgeFromPolygon(a, node);
        boolean bcFromPolygon =
            (node.next == b)
                ? node.isNextEdgeFromPolygon
                : isEdgeFromPolygon(node, b);
        boolean caFromPolygon =
            (b.next == a) ? b.isNextEdgeFromPolygon : isEdgeFromPolygon(a, b);
        triangles.add(new Triangle(a, abFromPolygon, node, bcFromPolygon, b, caFromPolygon));
        // Return the triangulated vertices to the tessellation
        triangles.add(new Triangle(a, abFromPolygon, node, bcFromPolygon, b, caFromPolygon));

        // remove two nodes involved
        removeNode(node, caFromPolygon);
        removeNode(node.next, caFromPolygon);
        node = startNode = b;
      }
      node = node.next;
    } while (node != startNode);

    return node;
  }

  /**
   * Attempt to split a polygon and independently triangulate each side. Return true if the polygon
   * was split.
   */
  private static boolean splitEarcut(final Node start, final List<Triangle> triangles) {
    // Search for a valid diagonal that divides the polygon into two.
    Node searchNode = start;
    Node nextNode;
    do {
      nextNode = searchNode.next;
      Node diagonal = nextNode.next;
      while (diagonal != searchNode.previous) {
        if (searchNode.idx != diagonal.idx && isValidDiagonal(searchNode, diagonal)) {
          // Split the polygon into two at the point of the diagonal
          Node splitNode = splitPolygon(searchNode, diagonal, isEdgeFromPolygon(searchNode, diagonal));
          // Filter the resulting polygon.
          searchNode = filterPoints(searchNode, searchNode.next);
          splitNode = filterPoints(splitNode, splitNode.next);
          // Attempt to earcut both of the resulting polygons
          earcutLinkedList(searchNode, triangles);
          earcutLinkedList(splitNode, triangles);
          // Finish the iterative search
          return true;
        }
        diagonal = diagonal.next;
      }
      searchNode = searchNode.next;
    } while (searchNode != start);
    return false;
  }

  /**
   * Computes if edge defined by a and b overlaps with a polygon edge
   */
  private static void checkIntersection(Node node) {
    for (Node a = node.next; a != node.previous; a = a.next) {
      for (Node b = a.next; b != a.previous; b = b.next) {
        checkIntersectionPoint(a, b);
      }
    }
  }

  private static Vector minValues(Vector... data) {
    double minX = Double.MAX_VALUE;
    double minY = Double.MAX_VALUE;
    double minZ = Double.MAX_VALUE;
    for (Vector val : data) {
      minX = Math.min(minX, val.x);
      minY = Math.min(minX, val.y);
      minZ = Math.min(minX, val.z);
    }
    return new Vector(minX, minY, minZ);
  }

  private static Vector maxValues(Vector... data) {
    double maxX = Double.MIN_VALUE;
    double maxY = Double.MIN_VALUE;
    double maxZ = Double.MIN_VALUE;
    for (Vector val : data) {
      maxX = Math.max(maxX, val.x);
      maxY = Math.max(maxX, val.y);
      maxZ = Math.max(maxX, val.z);
    }
    return new Vector(maxX, maxY, maxZ);
  }

  private static boolean rangeDisjoint(GeoPoint a1, GeoPoint a2, GeoPoint b1, GeoPoint b2) {
    Vector minA = minValues(a1, a2);
    Vector maxA = maxValues(a1, a2);
    Vector minB = minValues(b1, b2);
    Vector maxB = maxValues(b1, b2);
    return (minA.x >= maxB.x) ||
        (minA.y >= maxB.y) ||
        (minA.z >= maxB.z) ||
        (minB.x >= maxA.x) ||
        (minB.y >= maxA.y) ||
        (minB.z >= maxA.z);
  }

  private static void checkIntersectionPoint(final Node a, final Node b) {
    if (a == b) {
      return;
    }

    if (rangeDisjoint(a.point, a.next.point, b.point, b.next.point)) {
      return;
    }

    if (lineCrossesLine(a.point, a.next.point, b.point, b.next.point)) {
      String lineA = "line[" + a.point + " to " + a.next.point + "]";
      String lineB = "line[" + b.point + " to " + b.next.point + "]";
      throw new IllegalArgumentException("Polygon self-intersection between " + lineA + " and " + lineB);
    }
    if (a.isNextEdgeFromPolygon
        && b.isNextEdgeFromPolygon
        && lineOverlapLine(a.point, a.next.point, b.point, b.next.point)) {
      throw new IllegalArgumentException("Polygon ring self-intersection at " + a.point);
    }
  }

  /**
   * Computes if edge defined by a and b overlaps with a polygon edge
   */
  private static boolean isEdgeFromPolygon(final Node a, final Node b) {
    Node next = a;
    do {
      if (isPointInLine(next, next.next, a) && isPointInLine(next, next.next, b)) {
        return next.isNextEdgeFromPolygon;
      }
      if (isPointInLine(next, next.previous, a) && isPointInLine(next, next.previous, b)) {
        return next.previous.isNextEdgeFromPolygon;
      }
      next = next.next;
    } while (next != a);
    return false;
  }

  private static boolean isPointInLine(final Node a, final Node b, final Node point) {
    return isPointInLine(a, b, point.getX(), point.getY());
  }

  private static boolean isPointInLine(
      final Node a, final Node b, final double lon, final double lat) {
    final double dxc = lon - a.getX();
    final double dyc = lat - a.getY();

    final double dxl = b.getX() - a.getX();
    final double dyl = b.getY() - a.getY();

    if (dxc * dyl - dyc * dxl == 0) {
      if (Math.abs(dxl) >= Math.abs(dyl)) {
        return dxl > 0 ? a.getX() <= lon && lon <= b.getX() : b.getX() <= lon && lon <= a.getX();
      } else {
        return dyl > 0 ? a.getY() <= lat && lat <= b.getY() : b.getY() <= lat && lat <= a.getY();
      }
    }
    return false;
  }

  /**
   * Links two polygon vertices using a bridge.
   */
  private static Node splitPolygon(final Node a, final Node b, boolean edgeFromPolygon) {
    final Node a2 = new Node(a);
    final Node b2 = new Node(b);
    final Node an = a.next;
    final Node bp = b.previous;

    a.next = b;
    a.isNextEdgeFromPolygon = edgeFromPolygon;
    a.nextZ = b;
    b.previous = a;
    b.previousZ = a;
    a2.next = an;
    a2.nextZ = an;
    an.previous = a2;
    an.previousZ = a2;
    b2.next = a2;
    b2.isNextEdgeFromPolygon = edgeFromPolygon;
    b2.nextZ = a2;
    a2.previous = b2;
    a2.previousZ = b2;
    bp.next = b2;
    bp.nextZ = b2;

    return b2;
  }

  /**
   * Determines whether a diagonal between two polygon nodes lies within a polygon interior. (This
   * determines the validity of the ray.)
   */
  private static boolean isValidDiagonal(final Node a, final Node b) {
    if (isVertexEquals(a, b)) {
      // If points are equal then use it if they are valid polygons
      return isCCWPolygon(a, b);
    }
    return a.next.idx != b.idx
        && a.previous.idx != b.idx
        && edgeNeverIntersectsPolygon(a, a.point, b.point)
        && isLocallyInside(a, b)
        && isLocallyInside(b, a)
        && middleInsert(a, a.getX(), a.getY(), b.getX(), b.getY())
        // make sure we don't introduce collinear lines
        && !colinear(a.previous.point, a.point, b.point)
        && !colinear(a.point, b.point, b.next.point)
        && !colinear(a.next.point, a.point, b.point)
        && !colinear(a.point, b.point, b.previous.point);
  }

  /**
   * Determine whether the polygon defined between node start and node end is CCW
   */
  private static boolean isCCWPolygon(final Node start, final Node end) {
    // TODO: This does not yet work in 3D
    throw new RuntimeException("Invalid call to isCCWPolygon");
  }

  private static boolean isLocallyInside(final Node a, final Node b) {
    if (colinear(a.previous.point, a.point, a.next.point)) {
      // parallel
      return false;
    } else if (isRightHandRule(a.previous.point, a.point, a.next.point)) {
      // if a is cw
      return !isRightHandRule(a.point, b.point, a.next.point)
          && !isRightHandRule(a.point, a.previous.point, b.point);
    } else {
      // ccw
      return isRightHandRule(a.point, b.point, a.previous.point)
          || isRightHandRule(a.point, a.next.point, b.point);
    }
  }

  /**
   * Determine whether the middle point of a polygon diagonal is contained within the polygon
   */
  private static boolean middleInsert(
      final Node start, final double x0, final double y0, final double x1, final double y1) {
    Node node = start;
    Node nextNode;
    boolean lIsInside = false;
    final double lDx = (x0 + x1) / 2.0f;
    final double lDy = (y0 + y1) / 2.0f;
    do {
      nextNode = node.next;
      if (node.getY() > lDy != nextNode.getY() > lDy
          && lDx
          < (nextNode.getX() - node.getX())
          * (lDy - node.getY())
          / (nextNode.getY() - node.getY())
          + node.getX()) {
        lIsInside = !lIsInside;
      }
      node = node.next;
    } while (node != start);
    return lIsInside;
  }

  /**
   * Determines if the diagonal of a polygon is intersecting with any polygon elements.
   */
  private static boolean edgeNeverIntersectsPolygon(
      final Node start, final GeoPoint a, final GeoPoint b) {
    Node node = start;
    Node nextNode;
    do {
      nextNode = node.next;
      if (!isVertexEquals(node, a.x, a.y, a.z) && !isVertexEquals(node, b.x, b.y, b.z)) {
        if (linesIntersect(node.point, nextNode.point, a, b)) {
          return false;
        }
      }
      node = nextNode;
    } while (node != start);

    return true;
  }

  /**
   * Returns a positive value if points a, b, and c are arranged in counter-clockwise order,
   * negative value if clockwise, zero if collinear.
   */
  // see the "Orient2D" method described here:
  // http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
  // https://www.cs.cmu.edu/~quake/robust.html
  // Note that this one does not yet have the floating point tricks to be exact!
  public static int orient(GeoPoint a, GeoPoint b, GeoPoint c) {
    if (colinear(a, b, c)) {
      return 0;
    } else if (isRightHandRule(a, b, c)) {
      return 1;
    } else {
      return -1;
    }
  }

  /**
   * uses orient method to compute whether two line segments cross
   */
  public static boolean lineCrossesLine(
      GeoPoint a1,
      GeoPoint b1,
      GeoPoint a2,
      GeoPoint b2) {
    return orient(a2, b2, a1) * orient(a2, b2, b1) < 0
        && orient(a1, b1, a2) * orient(a1, b1, b2) < 0;
  }

  /**
   * uses orient method to compute whether two line overlap each other
   */
  public static boolean lineOverlapLine(
      GeoPoint a1,
      GeoPoint b1,
      GeoPoint a2,
      GeoPoint b2) {
    return orient(a2, b2, a1) == 0
        && orient(a2, b2, b1) == 0
        && orient(a1, b1, a2) == 0
        && orient(a1, b1, b2) == 0;
  }


  /**
   * Determines whether two line segments intersect. *
   */
  public static boolean linesIntersect(
      final GeoPoint a0,
      final GeoPoint a1,
      final GeoPoint b0,
      final GeoPoint b1) {
    // TODO: What is the difference between this and the lineCrossesLine method above?
    return isRightHandRule(a0, a1, b0) != isRightHandRule(a0, a1, b1)
        && isRightHandRule(b0, b1, a0) != isRightHandRule(b0, b1, a1);
  }

  /**
   * Creates a node and optionally links it with a previous node in a circular doubly-linked list
   */
  private static Node insertNode(
      final GeoPoint point,
      int index,
      final Node lastNode) {
    final Node node = new Node(point, index);
    if (lastNode == null) {
      node.previous = node;
      node.previousZ = node;
      node.next = node;
      node.nextZ = node;
    } else {
      node.next = lastNode.next;
      node.nextZ = lastNode.next;
      node.previous = lastNode;
      node.previousZ = lastNode;
      lastNode.next.previous = node;
      lastNode.nextZ.previousZ = node;
      lastNode.next = node;
      lastNode.nextZ = node;
    }
    return node;
  }

  /**
   * Removes a node from the doubly linked list
   */
  private static void removeNode(Node node, boolean edgeFromPolygon) {
    node.next.previous = node.previous;
    node.previous.next = node.next;
    node.previous.isNextEdgeFromPolygon = edgeFromPolygon;

    if (node.previousZ != null) {
      node.previousZ.nextZ = node.nextZ;
    }
    if (node.nextZ != null) {
      node.nextZ.previousZ = node.previousZ;
    }
  }

  /**
   * Determines if two point vertices are equal.
   */
  private static boolean isVertexEquals(final Node a, final Node b) {
    return isVertexEquals(a, b.point.x, b.point.y, b.point.z);
  }

  /**
   * Determines if two point vertices are equal.
   */
  private static boolean isVertexEquals(final Node a, final double x, final double y, final double z) {
    return a.point.x == x && a.point.y == y && a.point.z == z;
  }

  /**
   * Determine if three points are co-linear
   */
  public static boolean colinear(final GeoPoint a, final GeoPoint b, final GeoPoint c) {
    return Plane.arePointsCoplanar(a, b, c);
  }

  /**
   * Does the order a, b, c conform to the right-handed rule (curves counter-clockwise when seen from above)
   */
  public static boolean isRightHandRule(final GeoPoint a, final GeoPoint b, final GeoPoint c) {
    SidedPlane ac = new SidedPlane(a, c);
    return !ac.strictlyWithin(b);
  }

  /**
   * Compute whether point is in a candidate ear
   */
  private static boolean pointInEar(
      final GeoPoint p,
      final GeoPoint a,
      final GeoPoint b,
      final GeoPoint c) {
    //    (c.x - p.x) * (a.y - p.y) - (a.x - p.x) * (c.y - p.y)   >= 0
    //    (a.y - p.y) * (c.x - p.x) - (a.x - p.x) * (c.y - p.y)   >= 0
    // -( (p.y - a.y) * (c.x - p.x) - (p.x - a.x) * (c.y - p.y) ) >= 0
    //    (p.y - a.y) * (c.x - p.x) - (p.x - a.x) * (c.y - p.y)   <  0
    //
    //    area(a, p, c) < 0
    // && area(b, p, a) < 0
    // && area(c, p, b) < 0
    //
    //    isRightHandRule(a, p, c)
    // && isRightHandRule(b, p, a)
    // && isRightHandRule(c, p, b)
    //
    return isRightHandRule(a, p, c)
        && isRightHandRule(b, p, a)
        && isRightHandRule(c, p, b);
  }

  /**
   * Circular Doubly-linked list used for polygon coordinates
   */
  protected static class Node {
    // node index in the linked list
    private final int idx;
    // reference to the polygon for lat/lon values;
    private final GeoPoint point;
    // encoded x value
    private final int x;
    // encoded y value
    private final int y;
    // encoded y value
    private final int z;

    // previous node
    private Node previous;
    // next node
    private Node next;
    // previous z node
    private Node previousZ;
    // next z node
    private Node nextZ;
    // if the edge from this node to the next node is part of the polygon edges
    private boolean isNextEdgeFromPolygon;

    protected Node(final GeoPoint point, final int index) {
      this.idx = index;
      this.point = point;
      // casting to float is safe as original values for non-geo are represented as floats
      this.x = XYEncodingUtils.encode((float) point.x);
      this.y = XYEncodingUtils.encode((float) point.y);
      this.z = XYEncodingUtils.encode((float) point.z);
      this.previous = null;
      this.next = null;
      this.previousZ = null;
      this.nextZ = null;
      this.isNextEdgeFromPolygon = true;
    }

    /**
     * simple deep copy constructor
     */
    protected Node(Node other) {
      this.idx = other.idx;
      this.point = other.point;
      this.x = other.x;
      this.y = other.y;
      this.z = other.z;
      this.previous = other.previous;
      this.next = other.next;
      this.previousZ = other.previousZ;
      this.nextZ = other.nextZ;
      this.isNextEdgeFromPolygon = other.isNextEdgeFromPolygon;
    }

    /**
     * get the x value
     */
    public final double getX() {
      return point.x;
    }

    /**
     * get the y value
     */
    public final double getY() {
      return point.y;
    }

    /**
     * get the z value
     */
    public final double getZ() {
      return point.z;
    }

    @Override
    public String toString() {
      StringBuilder builder = new StringBuilder();
      if (this.previous == null) builder.append("||-");
      else builder.append(this.previous.idx).append(" <- ");
      builder.append(this.idx);
      if (this.next == null) builder.append(" -||");
      else builder.append(" -> ").append(this.next.idx);
      return builder.toString();
    }
  }

  /**
   * Triangle in the tessellated mesh
   */
  public static class Triangle {
    Node[] vertex;
    boolean[] edgeFromPolygon;

    protected Triangle(
        Node a,
        boolean isABfromPolygon,
        Node b,
        boolean isBCfromPolygon,
        Node c,
        boolean isCAfromPolygon) {
      this.vertex = new Node[]{a, b, c};
      this.edgeFromPolygon = new boolean[]{isABfromPolygon, isBCfromPolygon, isCAfromPolygon};
    }

    /**
     * Get if edge is shared with the polygon for the given edge.
     * This is required by lucene for converting triangles to fields.
     */
    public boolean isEdgefromPolygon(int startVertex) {
      return edgeFromPolygon[startVertex];
    }
    /**
     * pretty print the triangle vertices
     */
    @Override
    public String toString() {
      BiFunction<Node, Boolean, String> vertexStr = (p, e) -> "" + p.x + ", " + p.y + ", " + p.z + "[" + e + "]";
      return vertexStr.apply(vertex[0], edgeFromPolygon[0])
          + vertexStr.apply(vertex[1], edgeFromPolygon[1])
          + vertexStr.apply(vertex[2], edgeFromPolygon[2]);
    }
  }
}
