/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.lucene.spatial3d;

import org.apache.lucene.geo.GeoUtils;
import org.apache.lucene.geo.Polygon;
import org.apache.lucene.spatial3d.geom.*;
import org.apache.lucene.tests.util.LuceneTestCase;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import static java.lang.Math.toRadians;

/**
 * Test case for the GeoPolygon {@link Triangulator} class
 */
public class TestTriangulator extends LuceneTestCase {

  public void testRightHandRuleLatLon() {
    GeoPoint a = new GeoPoint(PlanetModel.SPHERE, 0, 0);
    GeoPoint b = new GeoPoint(PlanetModel.SPHERE, 0, 1);
    GeoPoint c = new GeoPoint(PlanetModel.SPHERE, 1, 1);
    assertTrue("Expected CCW", Triangulator.isRightHandRule(a, b, c));
    assertFalse("Expected CW", Triangulator.isRightHandRule(a, c, b));
  }

  public void testRightHandRuleLatLonNorthPole() {
    GeoPoint a = new GeoPoint(PlanetModel.SPHERE, toRadians(89), toRadians(0));
    GeoPoint b = new GeoPoint(PlanetModel.SPHERE, toRadians(89), toRadians(120));
    GeoPoint c = new GeoPoint(PlanetModel.SPHERE, toRadians(89), toRadians(-120));
    assertTrue("Expected CCW", Triangulator.isRightHandRule(a, b, c));
    assertFalse("Expected CW", Triangulator.isRightHandRule(a, c, b));
  }

  public void testRightHandRuleLatLonSouthPole() {
    GeoPoint a = new GeoPoint(PlanetModel.SPHERE, toRadians(-89), toRadians(0));
    GeoPoint b = new GeoPoint(PlanetModel.SPHERE, toRadians(-89), toRadians(120));
    GeoPoint c = new GeoPoint(PlanetModel.SPHERE, toRadians(-89), toRadians(-120));
    assertFalse("Expected CW", Triangulator.isRightHandRule(a, b, c));
    assertTrue("Expected CCW", Triangulator.isRightHandRule(a, c, b));
  }

  public void testRightHandRuleXYZ() {
    GeoPoint x = new GeoPoint(1, 1, 0, 0);
    GeoPoint y = new GeoPoint(1, 0, 1, 0);
    GeoPoint z = new GeoPoint(1, 0, 0, 1);
    // Normal right-hand ordering
    assertTrue("Expected CCW", Triangulator.isRightHandRule(x, y, z));
    assertTrue("Expected CCW", Triangulator.isRightHandRule(y, z, x));
    assertTrue("Expected CCW", Triangulator.isRightHandRule(z, x, y));
    // Opposite ordering
    assertFalse("Expected CW", Triangulator.isRightHandRule(y, x, z));
    assertFalse("Expected CW", Triangulator.isRightHandRule(x, z, y));
    assertFalse("Expected CW", Triangulator.isRightHandRule(z, y, x));
  }

  public void testTriangleWindingOrderClockwise() {
    double[] latitudes = new double[]{0.0, 2.0, 0.0, 0.0};
    double[] longitudes = new double[]{-1.0, 0.0, 1.0, -1.0};
    Polygon polygon = new Polygon(latitudes, longitudes);
    assertEquals(GeoUtils.WindingOrder.CCW, polygon.getWindingOrder());
    ArrayList<GeoPoint> points = Geo3dTestUtil.latLonArrayToGeoPoints(latitudes, longitudes);
    GeoUtils.WindingOrder direction = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CW, direction);
  }

  public void testTriangleWindingOrderCounterClockwise() {
    double[] latitudes = new double[]{0.0, 2.0, 0.0, 0.0};
    double[] longitudes = new double[]{1.0, 0.0, -1.0, 1.0};
    Polygon polygon = new Polygon(latitudes, longitudes);
    assertEquals(GeoUtils.WindingOrder.CW, polygon.getWindingOrder());
    ArrayList<GeoPoint> points = Geo3dTestUtil.latLonArrayToGeoPoints(latitudes, longitudes);
    GeoUtils.WindingOrder direction = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CCW, direction);
  }

  public void testSquareWindingOrderClockwise() {
    double[] lats = new double[]{0.0, 2.0, 2.0, 0.0, 0.0};
    double[] lons = new double[]{-1.0, -1.0, 1.0, 1.0, -1.0};
    Polygon polygon = new Polygon(lats, lons);
    assertEquals(GeoUtils.WindingOrder.CCW, polygon.getWindingOrder());
    ArrayList<GeoPoint> points = Geo3dTestUtil.latLonArrayToGeoPoints(lats, lons);
    GeoUtils.WindingOrder direction = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CW, direction);
  }

  public void testSquareWindingOrderCounterClockwise() {
    double[] lats = new double[]{0.0, 2.0, 2.0, 0.0, 0.0};
    double[] lons = new double[]{1.0, 1.0, -1.0, -1.0, 1.0};
    Polygon polygon = new Polygon(lats, lons);
    assertEquals(GeoUtils.WindingOrder.CW, polygon.getWindingOrder());
    ArrayList<GeoPoint> points = Geo3dTestUtil.latLonArrayToGeoPoints(lats, lons);
    GeoUtils.WindingOrder direction = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CCW, direction);
  }

  public void testHighLatitudeConcave() {
    int baseLatitude = 80;
    double colinearLatitude = calculateColinearLatitude(baseLatitude, -45, 45);
    for (int offset = 0; offset <= 10; offset++) {
      int highLatitude = baseLatitude + offset;
      GeoPoint a = new GeoPoint(PlanetModel.SPHERE, toRadians(baseLatitude), toRadians(45));
      GeoPoint b = new GeoPoint(PlanetModel.SPHERE, toRadians(highLatitude), toRadians(0));
      GeoPoint c = new GeoPoint(PlanetModel.SPHERE, toRadians(baseLatitude), toRadians(-45));
      boolean expected = highLatitude > colinearLatitude;
      String orientation = expected ? "convex" : "concave";
      boolean convex = Triangulator.isRightHandRule(a, b, c);
      String message = "Should be " + orientation + " when latitude changes from " + baseLatitude + " to " + highLatitude + " and back";
      assertEquals(message, expected, convex);
    }
  }

  private double calculateColinearLatitude(double baseLatitude, double minLongitude, double maxLongitude) {
    GeoPoint a = new GeoPoint(PlanetModel.SPHERE, toRadians(baseLatitude), toRadians(minLongitude));
    GeoPoint b = new GeoPoint(PlanetModel.SPHERE, toRadians(baseLatitude), toRadians(maxLongitude));
    double colinearLatitudeFromA = Math.toDegrees(Math.atan(a.z / a.x));
    double colinearLatitudeFromB = Math.toDegrees(Math.atan(b.z / b.x));
    assertEquals("Should have same results", colinearLatitudeFromA, colinearLatitudeFromB, 0.00001);
    GeoPoint x = new GeoPoint(PlanetModel.SPHERE, toRadians(colinearLatitudeFromA), toRadians((maxLongitude + minLongitude) / 2.0));
    assertTrue("Points should be colinear at " + colinearLatitudeFromA + " latitude", Triangulator.colinear(a, x, b));
    return colinearLatitudeFromA;
  }

  public void testTriangulateTriangle() {
    ArrayList<GeoPoint> points = Geo3dTestUtil.latLonArrayToGeoPoints(
        new double[]{0.0, 2.0, 0.0},
        new double[]{1.0, 0.0, -1.0});
    GeoUtils.WindingOrder direction = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CCW, direction);
    GeoPolygon poly = GeoPolygonFactory.makeGeoPolygon(PlanetModel.SPHERE, new GeoPolygonFactory.PolygonDescription(points));
    assertEquals(1, Objects.requireNonNull(Triangulator.triangulate(poly, random().nextBoolean())).size());
  }

  public void testTriangulateSquare() {
    ArrayList<GeoPoint> points = Geo3dTestUtil.latLonArrayToGeoPoints(
        new double[]{-1.0, -1.0, 1.0, 1.0},
        new double[]{-1.0, 1.0, 1.0, -1.0});
    GeoUtils.WindingOrder direction = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CCW, direction);
    GeoPolygon poly = GeoPolygonFactory.makeGeoPolygon(PlanetModel.SPHERE, new GeoPolygonFactory.PolygonDescription(points));
    assertEquals(2, Objects.requireNonNull(Triangulator.triangulate(poly, random().nextBoolean())).size());
  }

  public void testSmallTessellation() {
    ArrayList<GeoPoint> points = Geo3dTestUtil.createRegularGeoPolygonPoints(0.0, 0.0, 100000, 80);
    GeoUtils.WindingOrder order = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CCW, order);
    List<Triangulator.Triangle> triangles = Triangulator.triangulate(points, random().nextBoolean());
    assertTrue(Objects.requireNonNull(triangles).size() > 10);
  }

  public void testSmallTessellationNorthPole() {
    ArrayList<GeoPoint> points = Geo3dTestUtil.createPolarGeoPolygonPoints(true, 100000, 80);
    GeoUtils.WindingOrder order = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CCW, order);
    List<Triangulator.Triangle> triangles = Triangulator.triangulate(points, random().nextBoolean());
    assertTrue(Objects.requireNonNull(triangles).size() > 10);
  }

  public void testSmallTessellationSouthPole() {
    ArrayList<GeoPoint> points = Geo3dTestUtil.createPolarGeoPolygonPoints(false, 100000, 80);
    GeoUtils.WindingOrder order = Triangulator.calculateWindingOrder(points);
    assertEquals(GeoUtils.WindingOrder.CCW, order);
    List<Triangulator.Triangle> triangles = Triangulator.triangulate(points, random().nextBoolean());
    assertTrue(Objects.requireNonNull(triangles).size() > 10);
  }

  public void testSimpleTessellation() {
    GeoPolygon poly = Geo3dTestUtil.createRegularGeoPolygon(0.0, 0.0, 100000, 10000);
    GeoUtils.WindingOrder order = Triangulator.calculateWindingOrder(poly.getOuterShellPoints());
    assertEquals(GeoUtils.WindingOrder.CCW, order);
    List<Triangulator.Triangle> triangles = Triangulator.triangulate(poly, random().nextBoolean());
    assertTrue(Objects.requireNonNull(triangles).size() > 10);
  }

  public void testSmallTessellationWithHole() {
    ArrayList<GeoPoint> inner = Geo3dTestUtil.latLonArrayToGeoPoints(
        new double[]{-1.0, -1.0, 0.5, 1.0, 1.0, 0.5, -1.0},
        new double[]{1.0, -1.0, -0.5, -1.0, 1.0, 0.5, 1.0});
    ArrayList<GeoPoint> inner2 = Geo3dTestUtil.latLonArrayToGeoPoints(
        new double[]{-1.0, -1.0, 0.5, 1.0, 1.0, 0.5, -1.0},
        new double[]{-2.0, -4.0, -3.5, -4.0, -2.0, -2.5, -2.0});
    GeoPolygon poly = Geo3dTestUtil.createRegularGeoPolygonWithHoles(0.0, 0.0, 100000, 80, inner, inner2);
    GeoUtils.WindingOrder order = Triangulator.calculateWindingOrder(poly.getOuterShellPoints());
    assertEquals(GeoUtils.WindingOrder.CCW, order);
    List<Triangulator.Triangle> triangles = Triangulator.triangulate(poly, random().nextBoolean());
    assertTrue(Objects.requireNonNull(triangles).size() > 10);
  }

  public void testSimpleTessellationWithHole() {
    ArrayList<GeoPoint> inner = Geo3dTestUtil.latLonArrayToGeoPoints(
        new double[]{-1.0, -1.0, 0.5, 1.0, 1.0, 0.5, -1.0},
        new double[]{1.0, -1.0, -0.5, -1.0, 1.0, 0.5, 1.0});
    ArrayList<GeoPoint> inner2 = Geo3dTestUtil.latLonArrayToGeoPoints(
        new double[]{-1.0, -1.0, 0.5, 1.0, 1.0, 0.5, -1.0},
        new double[]{-2.0, -4.0, -3.5, -4.0, -2.0, -2.5, -2.0});
    GeoPolygon poly = Geo3dTestUtil.createRegularGeoPolygonWithHoles(0.0, 0.0, 100000, 100, inner, inner2);
    GeoUtils.WindingOrder order = Triangulator.calculateWindingOrder(poly.getOuterShellPoints());
    assertEquals(GeoUtils.WindingOrder.CCW, order);
    List<Triangulator.Triangle> triangles = Triangulator.triangulate(poly, random().nextBoolean());
    assertTrue(Objects.requireNonNull(triangles).size() > 10);
  }
}
