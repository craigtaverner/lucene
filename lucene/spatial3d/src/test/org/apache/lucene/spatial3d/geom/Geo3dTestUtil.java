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
package org.apache.lucene.spatial3d.geom;

import org.apache.lucene.geo.GeoUtils;
import org.apache.lucene.spatial3d.Triangulator;
import org.apache.lucene.util.SloppyMath;

import java.util.ArrayList;

import static java.lang.Math.toRadians;
import static org.junit.Assert.assertEquals;

public class Geo3dTestUtil {

  /**
   * This is a copy of GeoTestUtils.createRegularPolygon, modified for GeoPolygon.
   * Makes an n-gon, centered at the provided lat/lon, and each vertex approximately distanceMeters
   * away from the center.
   *
   * <p>Do not invoke me across the dateline or a pole!!
   */
  public static ArrayList<GeoPoint> createRegularGeoPolygonPoints(
      double centerLat, double centerLon, double radiusMeters, int gons) {

    // Points go counterclockwise, so
    ArrayList<GeoPoint> points = new ArrayList<>();
    for (int i = 0; i < gons; i++) {
      double angle = i * (360.0 / gons);
      // System.out.println("  angle " + angle);
      double x = Math.cos(toRadians(angle));
      double y = Math.sin(toRadians(angle));
      double factor = 2.0;
      double step = 1.0;
      int last = 0;

      // Iterate out along one spoke until we hone in on the point that's nearly exactly
      // radiusMeters from the center:
      while (true) {

        // TODO: we could in fact cross a pole?  Just do what surpriseMePolygon does?
        double lat = centerLat + y * factor;
        GeoUtils.checkLatitude(lat);
        double lon = centerLon + x * factor;
        GeoUtils.checkLongitude(lon);
        double distanceMeters = SloppyMath.haversinMeters(centerLat, centerLon, lat, lon);

        if (Math.abs(distanceMeters - radiusMeters) < 0.1) {
          // Within 10 cm: close enough!
          points.add(new GeoPoint(PlanetModel.SPHERE, toRadians(lat), toRadians(lon)));
          break;
        }

        if (distanceMeters > radiusMeters) {
          // too big
          factor -= step;
          if (last == 1) {
            step /= 2.0;
          }
          last = -1;
        } else if (distanceMeters < radiusMeters) {
          // too small
          factor += step;
          if (last == -1) {
            step /= 2.0;
          }
          last = 1;
        }
      }
    }
    return points;
  }

  /**
   * Create a regular polygon around a pole
   */
  public static ArrayList<GeoPoint> createPolarGeoPolygonPoints(boolean north, double radiusMeters, int gons) {
    double offsetLat = radiusMeters/40000000 * 360.0;
    double latitude = north ? 90 - offsetLat : -90 + offsetLat;
    double offsetLon = north ? 360.0 / gons : -360.0 / gons;

    // Points go counterclockwise, so
    ArrayList<GeoPoint> points = new ArrayList<>();
    for (int i = 0; i < gons; i++) {
      double longitude = i * offsetLon;
      if (longitude < -180) {
        longitude += 360.0;
      }
      if (longitude > 180) {
        longitude -= 360.0;
      }
      points.add(new GeoPoint(PlanetModel.SPHERE, toRadians(latitude), toRadians(longitude)));
    }
    return points;
  }

  public static ArrayList<GeoPoint> latLonArrayToGeoPoints(double[] lats, double[] lons) {
    if (lats.length != lons.length) {
      throw new IllegalArgumentException("Lat/Lon arrays must have equal length: " + lats.length + " != " + lons.length);
    }
    ArrayList<GeoPoint> points = new ArrayList<>();
    for (int i = 0; i < lats.length; i++) {
      points.add(new GeoPoint(PlanetModel.SPHERE, toRadians(lats[i]), toRadians(lons[i])));
    }
    return points;
  }

  public static GeoPolygon createRegularGeoPolygon(double centerLat, double centerLon, double radiusMeters, int gons) {
    ArrayList<GeoPoint> points = createRegularGeoPolygonPoints(centerLat, centerLon, radiusMeters, gons);
    assertEquals("Expect outer polygon to have CCW winding order", GeoUtils.WindingOrder.CCW, Triangulator.calculateWindingOrder(points));
    GeoPolygonFactory.PolygonDescription pd = new GeoPolygonFactory.PolygonDescription(points);
    return GeoPolygonFactory.makeGeoPolygon(PlanetModel.SPHERE, pd);
  }

  @SafeVarargs
  public static GeoPolygon createRegularGeoPolygonWithHoles(
      double centerLat, double centerLon, double radiusMeters, int gons, ArrayList<GeoPoint>... holes_points) {
    ArrayList<GeoPoint> points = createRegularGeoPolygonPoints(centerLat, centerLon, radiusMeters, gons);
    assertEquals("Expect outer polygon to have CCW winding order", GeoUtils.WindingOrder.CCW, Triangulator.calculateWindingOrder(points));

    ArrayList<GeoPolygonFactory.PolygonDescription> holes = new ArrayList<>();
    for (ArrayList<GeoPoint> hole_points : holes_points) {
      GeoPolygonFactory.PolygonDescription holeDescription = new GeoPolygonFactory.PolygonDescription(hole_points);
      holes.add(holeDescription);
    }
    GeoPolygonFactory.PolygonDescription pd = new GeoPolygonFactory.PolygonDescription(points, holes);
    return GeoPolygonFactory.makeGeoPolygon(PlanetModel.SPHERE, pd);
  }

  public static GeoPolygon createPolarGeoPolygon(boolean north, double radiusMeters, int gons) {
    ArrayList<GeoPoint> points = createPolarGeoPolygonPoints(north, radiusMeters, gons);
    assertEquals("Expect outer polygon to have CCW winding order", GeoUtils.WindingOrder.CCW, Triangulator.calculateWindingOrder(points));
    GeoPolygonFactory.PolygonDescription pd = new GeoPolygonFactory.PolygonDescription(points);
    return GeoPolygonFactory.makeGeoPolygon(PlanetModel.SPHERE, pd);
  }

}
