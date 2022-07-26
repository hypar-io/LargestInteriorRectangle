# Largest Interior Rectangle — Hypar.Elements Version

Largest Interior Rectangle implementation in C# for use with Hypar Elements.
This code is taken nearly verbatim from [https://github.com/Evryway/lir](https://github.com/Evryway/lir) — The Unity dependencies have been replaced with Hypar.Elements. Many many thanks to @Evryway for his work, and for making it available via the MIT license! 

### A word of warning
This port is _not_ extensively tested — I pretty much just mapped from one API to another, so it should be used with caution. Feel free to file an issue if you find any problems. It is intended only for 2D use — the original code used Vector2s, so don't expect Vector3s with Z coordinates to respect/preserve that information.

## REQUIREMENTS

The Hypar.Elements NuGet package is required to use this library.
## USAGE

### Usage Example
```cs
using Elements.Geometry;
using Elements.LIR;
Polygon pgon = new Polygon
(
    (9,-7),
    (12,-6),
    (8,3),
    (10,6),
    (12,7),
    (1,9),
    (-8,7),
    (-6,6),
    (-4,6),
    (-6,2),
    (-6,0),
    (-7,-5),
    (-2,-7),
    (1,-3),
    (5,-7),
    (8,-4)
);
LargestInteriorRectangle.CalculateLargestInteriorRectangle(pgon, out var best);
var rect = best.Polygon;
var mc = new ModelCurve(pgon, BuiltInMaterials.XAxis);
var mc2 = new ModelCurve(rect, BuiltInMaterials.YAxis);

var model = new Model();
model.AddElement(mc);
model.AddElement(mc2);
return model;
```
(See also the .NET Interactive Notebook called `UsageExample.dib` to explore interactively.)



## LargestInteriorRectangle.cs

The LargestInteriorRectangle class contains various static methods for processing Vector3 arrays.
The primary methods are:

#### CalculateInteriorCells(Vector3[] vs, out double[] xs, out double[] ys, out int[,] cells)

This method takes a list of Vector3s representing a simple polygon (ordered counter-clockwise, concave, non-intersecting) and
generates three output arrays : the xs (ordered, filtered unique x coordinates for all vertices), the ys
(ordered, filtered, unique y coordinates for all vertices) and the 2D cells array representing each cell of the
projection of the xs and ys arrays into a 2D rectangle array. Each cell is marked as exterior (0) or interior
(1) based on it's position relative to the polygon.

#### CalculateLargestInteriorRectangle(double[] xs, double ys[], int[,] cells, out Bound2D best)

this method takes the xs, ys and cells arrays created by CalculateInteriorCells and calculates
the axis-aligned Largest Interior Rectangle. The bound of this rectangle is output in best.

Along with these methods, there are a variety of other geometry processing methods:

#### CalculateConcavePolygon (Vector3[] vs)

Sorts a point set into an ordered array (counter-clockwise) representing a concave polygon.

#### CalculateConvexHull(Vector3[] vs, out Vector3[] hull_vs)

Calculates a convex hull, given a point set.

#### CalculateConvexPolygonArea(Vector3[] vs)

calculates the area of a convex polygon, given a CCW ordered point array representing that convex polygon.

#### CalculatePrimaryAxis(Vector3[] vs, out Vector3 axis, out double eigenvalue)

Given a point set in vs, calculates the primary axis, and the corresponding eigenvalue.
The primary axis is the line through the point set representing the largest variance.

#### CalculateEigenValues(Double2x2 mat, out double v1, out double v2)

Given a 2x2 matrix, calculates the eigenvalues.

#### CalculateEigenVector(Double2x2 A, double eigenvalue)

Given a 2x2 matrix and an eigenvalue, calculates an eigenvector.

## Bound2D.cs

This class represents a 2D bound. The bound is not required to be axis aligned. The bound is
constructed from a centre point, a 2D size and the primary axis (the first element of the 2D size).

## Extensions.cs

This class contains some Vector3 extension methods to simplify the code (Rotate for simple
2D rotations, and F3() to make logging more concise).

## LargestInteriorRectangleTests.cs

This class contains a selection of tests across the LargestInteriorRectangle methods.

## BACKGROUND READING

please see https://www.evryway.com/largest-interior/ for background details.
