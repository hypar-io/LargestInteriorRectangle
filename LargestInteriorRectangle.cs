using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Elements.Geometry;

// static class to calculate the internal largest area rectangle of a simple polygon.
//
// reference : 

// largest interior rectangle:
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.47.3370&rep=rep1&type=pdf
// https://www.sciencedirect.com/science/article/pii/S0925772115000759
// https://mathoverflow.net/questions/105837/get-largest-inscribed-rectangle-of-a-concave-polygon
// https://www.sciencedirect.com/science/article/pii/0925772195000410
// https://mathoverflow.net/questions/105164/covering-a-polygon-with-rectangles

//
// PCA :
// https://builtin.com/data-science/step-step-explanation-principal-component-analysis
// https://medium.com/analytics-vidhya/eigenvectors-and-eigenvalues-and-there-use-in-principal-component-analysis-machine-learning-1f97fdbdb303
//
// Convex Hull :
// http://geomalgorithms.com/a10-_hull-1.html
//
// Bounding rectangle:
// https://www.tvhoang.com/articles/2018/12/rotating-calipers
// 
// method:
// https://journals.ut.ac.ir/article_71280_2a21de484e568a9e396458a5930ca06a.pdf (x)
// https://www.evryway.com/interior-rectangle/

namespace Elements.LIR
{

    public static class LargestInteriorRectangle
    {

        // I think this is obviated by the transformation from vector2 => vector3
        public static Vector3[] Vector3XZToVector3(IEnumerable<Vector3> vec3s_xz)
        {
            var c = vec3s_xz.Count();
            var outv2 = new Vector3[c];
            int i = 0;
            foreach (var v3 in vec3s_xz)
            {
                var v2 = new Vector3(v3.X, v3.Z);
                outv2[i] = v2;
                i++;
            }
            return outv2;
        }

        public static Vector3[] Vector3ToVector3XZ(IEnumerable<Vector3> vec2s, double y = 0)
        {
            var c = vec2s.Count();
            var outv3 = new Vector3[c];
            int i = 0;
            foreach (var v2 in vec2s)
            {
                var v3 = new Vector3(v2.X, y, v2.Y);
                outv3[i] = v3;
                i++;
            }
            return outv3;
        }

        // Calculate covariance matrix.
        // https://math.stackexchange.com/q/711886
        public static Double2x2 CovarianceMatrix(Vector3[] vec2s)
        {
            var n = vec2s.Length;
            double xt = 0.0f, yt = 0.0f;
            for (int i = 0; i < n; i++)
            {
                var v = vec2s[i];
                xt += v.X;
                yt += v.Y;
            }
            // x and y average values
            double xmean = xt / n;
            double ymean = yt / n;

            //Console.WriteLine($"mean : {xmean} {ymean}");

            double xvt = 0.0f, yvt = 0.0f, xycvt = 0.0f;
            for (int i = 0; i < n; i++)
            {
                var v = vec2s[i];
                var xd = v.X - xmean;
                var yd = v.Y - ymean;

                xvt += xd * xd;
                yvt += yd * yd;
                xycvt += xd * yd;
            }

            // we are assuming we have ALL the samples - so we divide by N, rather than N-1
            // if we're doing a subset sample, divide by N-1.
            // this is population covariance vs sample covariance.
            // (see https://www.visiondummy.com/2014/03/divide-variance-n-1/ )
            var divis = n;      // n-1
            double xvar = xvt / divis;
            double yvar = yvt / divis;
            double xycov = xycvt / divis;

            //Console.WriteLine($"A : [ ({xvar} {xycov}) , ({xycov}, {yvar}) ]");
            return new Double2x2(xvar, xycov, xycov, yvar);
        }

        // take quadratic equation in form ax^2 + bx + c = 0
        // return roots (values of x) in x1 and x2
        // return true if solvable, false otherwise.
        // equation is https://en.wikipedia.org/wiki/Quadratic_formula
        public static bool SolveQuadratic(double a, double b, double c, out double x1, out double x2)
        {
            x1 = 0;
            x2 = 0;

            // can't solve if a is zero.
            if (a == 0) return false;
            var part = (b * b) - (4 * a * c);
            // can't solve if part is negative (sqrt of negative number)
            if (part < 0) return false;
            var psqrt = Math.Sqrt(part);
            x1 = (-b + psqrt) / (2 * a);
            x2 = (-b - psqrt) / (2 * a);
            return true;
        }


        // v = lambda
        // det(A - vI) = 0
        // see https://www.youtube.com/watch?v=tXlMbAxbUI4
        // Khan Academy does it backwards - det(vI -A) = 0
        // these are functionally identical. I'm going with the former.
        // m =  a,b     vI =   v, 0
        //      c,d            0, v
        // det |a - v, b    |    = 0
        //     |c    , d - v|
        // det of 2x2 is ad - bc
        // so (a-v)*(d-v) - (bc) = 0
        // so (ad -dv -av + vv - bc = 0
        // characteristic polynomial:
        // vv -(a+d)v + ad - bc = 0
        //
        // the video also suggests the characteristic polynomial is
        // vv - (trace(A))v + det(A) = 0
        // which is exactly as I've worked out above, given Trace(A) is a+d.

        public static bool CalculateEigenValues(Double2x2 mat, out double v1, out double v2)
        {
            double a = mat.c0.x;
            double b = mat.c1.x;
            double c = mat.c0.y;
            double d = mat.c1.y;

            // note - coming in as a covariance matrix, b and c are identical.
            // this solves the general (not necessarily covariance matrix) case.
            double qa = 1;
            double qb = -(a + d);
            double qc = (a * d) - (b * c);
            var ok = SolveQuadratic(qa, qb, qc, out v1, out v2);
            return ok;
        }


        // Calculate EigenVector
        // W is an eigenvector. (W.x, W.y) - with components x,y.
        // A * W = vW
        // or,
        // (A - vI)W = 0
        // the EigenSpace is the line which contains all eigenvectors for a specific
        // eigenvalue of the matrix A.
        //
        // with a matrix A defined as | a , b |
        //                            | c , d |
        //
        // then we know that
        // | a , b | |x|  =   |vx|
        // | c , d | |y|      |vy|
        //
        // giving two equations:
        // ax + by = vx
        // cx + dy = vy
        // or, 
        // (a-v)x + by = 0
        // cx + (d-v)y = 0
        // hence,
        // -(a-v)x = by
        // or , (v-a)x = by
        // so (v-a)/b * x = y
        // and x = b/(v-a) * y
        // if y is 1, then x is b / (v-a)
        // cx = -(d-v)y
        // or, cx = (v-d)y
        // so x = ((v-d)/c) y
        // if y is 1, then x = (v-d)/c
        //
        // there are infinite solutions to this, so we pick a "nice" one ...
        // e.g. we set y to 1.
        // so, x = b / (v-a) and y = 1.
        // then - normalize!
        //
        // now, BOTH equations *should* give the same value.
        // but, numerical errors creep in.
        // because we're using the Covariance matrix, the Variance values
        // approximate the EigenValues as the covariance drops to 0, meaning
        // ONE of the two equations tends to BAD, and the other should be preferred.
        // see https://www.mathsisfun.com/algebra/eigenvalue.html

        public static Vector3 CalculateEigenVector(Double2x2 A, double eigenvalue)
        {
            double v = eigenvalue;
            double a = A.c0.x;
            double b = A.c1.x;
            double c = A.c0.y;
            double d = A.c1.y;

            // decide which equation to use, based on mag of variance (XX, YY)
            // against the eigenvector - we don't want a near-zero value here!
            var aq = Math.Abs(a - v);
            var dq = Math.Abs(d - v);
            var numer = (aq > dq) ? b : (v - d);
            var denom = (aq > dq) ? (v - a) : c;
            var dok = Math.Abs(denom) > Vector3.EPSILON;
            var x = dok ? numer / denom : 1.0f;
            var y = dok ? 1 : 0;

            // as eigenvectors can point either way down the eigenspace, let's
            // go for ALWAYS POSITIVE X.
            if (x < 0) { x = -x; y = -y; }
            return new Vector3(x, y).Unitized();
        }


        public static bool CalculatePrimaryAxis(Vector3[] vs, out Vector3 axis, out double eigenvalue)
        {
            // default to x-axis as primary axis.
            axis = Vector3.XAxis;
            eigenvalue = 1.0f;

            var A = CovarianceMatrix(vs);
            var ok = CalculateEigenValues(A, out double v1, out double v2);
            if (!ok) return false;

            //Console.WriteLine($"EigenValues: {v1} {v2}");
            eigenvalue = Math.Max(Math.Abs(v1), Math.Abs(v2));
            // use largest magnitude eigenvalue to calculate an eigenvector.
            axis = CalculateEigenVector(A, eigenvalue);
            //Console.WriteLine($"Vector: ({axis.x}, {axis.y}) for EigenValue {eigenvalue}");

            return true;
        }

        // Calculate convex polygon area.
        // Uses shoelace algo.
        // https://erkaman.github.io/posts/area_convex_polygon.html
        public static double CalculateConvexPolygonArea(Vector3[] vs)
        {
            var c = vs.Length;
            double s = 0.0f;
            for (int i = 0; i < c; ++i)
            {
                var v = vs[i];
                var v2 = vs[(i + 1) % c];
                s += (v.X * v2.Y) - (v2.X * v.Y);
            }
            //Console.WriteLine($"polygon area : {s*0.5f}");
            return 0.5f * s;
        }

        public static Vector3[] CalculateConcavePolygon(Vector3[] vs)
        {
            // step 1: find bottom-right point in vs.
            // minimum y (and optionally maximum x, if two or more min-y points are colinear with x-axis)
            var c = vs.Length;


            Vector3 p0 = vs[0];
            int idx = 0;
            List<Vector3> vs_sorted = new List<Vector3> { p0 };
            for (int i = 1; i < c; i++)
            {
                var pi = vs[i];
                if ((pi.Y < p0.Y) || (pi.Y == p0.Y && pi.X > p0.X)) { p0 = pi; idx = i; }
                vs_sorted.Add(pi);
            }

            vs_sorted.RemoveAt(idx);
            // at this point, vs_sorted does NOT contain p0 - we're going to put it at the front shortly.

            // sort on radial of vector p0->px.
            vs_sorted.Sort((p1, p2) =>
            {
                // 2d cross product. (a,b) , (c,d) -> ad - cb 
                // "on the left" (looking from above, anticlockwise)
                // gives a positive value for d.

                var p1x = p1.X - p0.X;
                var p2x = p2.X - p0.X;
                var p1y = p1.Y - p0.Y;
                var p2y = p2.Y - p0.Y;


                var d = (p1x * p2y) - (p2x * p1y);
                // leftmost (positive value) sorts lowest.
                var r = -d.CompareTo(0);
                // in the colinear case (p0 -> p1 -> p2 is a straight line, or p0 -> p2 -> p1 is a straight line)
                if (r == 0)
                {
                    var d1 = (p1x * p1x) + (p1y * p1y);
                    var d2 = (p2x * p2x) + (p2y * p2y);
                    // closest (smallest) sorts lowest.
                    r = d1.CompareTo(d2);
                }
                return r;
            });

            vs_sorted.Insert(0, p0);
            return vs_sorted.ToArray();

        }


        // Calculate the convex hull, given a point set.
        // Uses Graham Scan.
        public static bool CalculateConvexHull(Vector3[] vs, out Vector3[] hull_vs)
        {
            if (vs.Length < 1) { hull_vs = new Vector3[0]; return false; }


            var vs_sorted = CalculateConcavePolygon(vs);
            var c = vs_sorted.Length;

            //for (int i = 0; i < s_sorted.Count; i++) Console.WriteLine(vs_sorted[i]);


            // we're ready to construct the convex hull.

            Stack<Vector3> hull = new Stack<Vector3>();
            var v0 = vs_sorted[0];
            hull.Push(v0);

            // ensure we do NOT begin with a duplicate of our start point.
            int second = 1;
            while (second < c)
            {
                var v = vs_sorted[second];
                second++;
                if ((v - v0).LengthSquared() > 0)
                {
                    hull.Push(v);
                    break;
                }
            }
            if (hull.Count < 2 || c < 3)
            {
                hull_vs = hull.ToArray();
                return false;
            }

            // c+1 gives us a duplicate start point at the end - this should ensure
            // we don't get any colinearity on the final edge.
            // we will remove the duplicate later.
            for (int i = second; i < c + 1; i++)
            {
                var t = vs_sorted[(i % c)];

                while (hull.Count > 1)
                {
                    var s = hull.Pop();

                    // we do NOT want duplicate vertices.
                    // push the new one (same as the old one)
                    // and move on to the next vert.
                    if ((s - t) == (0, 0, 0))
                    {
                        break;
                    }

                    var r = hull.Peek();

                    var d = ((s.X - r.X) * (t.Y - r.Y)) - ((t.X - r.X) * (s.Y - r.Y));

                    // we must be POSITIVE (on the left).
                    // negative (on the right) means we want to continue back down the stack.
                    if (d < 0) continue;

                    // we may be colinear - in which case, we want to RETAIN the point (s or t) that
                    // is FURTHEST from r.
                    // this SHOULD always be t (rather than s) as long as the initial CCW sort sorts
                    // closest-first when colinear.
                    // we may also be a duplicate end/start degenerate case.
                    if (d == 0)
                    {
                        var sm = (s - r).LengthSquared();
                        var tm = (t - r).LengthSquared();
                        if (sm > tm)
                        {
                            t = s;                 // use s instead of t, as it's further away.
                        }
                        if (hull.Count < 2) break;  // cannot safely pop any more from the hull.
                        continue;
                    }

                    hull.Push(s);
                    break;
                }
                hull.Push(t);
            }

            // take off duplicate start point.
            // in the degenerate case (all colinear) this would reduce the hull to a single point,
            // so catch that.
            if (hull.Count > 2)
            {
                hull.Pop();
            }

            hull_vs = hull.Reverse().ToArray();
            //Console.WriteLine(hull_vs.Length);
            return hull_vs.Length > 2;
        }


        // Calculate Smallest Enclosing Rectangle. Not AABB.
        // takes a convex polygon.
        // REQUIRED - no duplicate points, no colinearity.
        // REQUIRED - points are ordered counter-clockwise, +X is right, +Y is up.
        // uses rotating calipers algorithm.
        public static bool CalculateSmallestEnclosingRectangle(Vector3[] vs, out Bound2D bound)
        {
            var vs_area = CalculateConvexPolygonArea(vs);
            if (vs_area <= 0)
            {
                // if we have (for example) a colinear set of points, do we want to make a stab
                // at a bound that "covers" the points? (e.g. the axis describes the colinear set,
                // size of axis_b is zero, centre is midpoint of the span, etc)?
                // for now - no.

                bound = new Bound2D(new Vector3(0, 0, 0), Vector3.XAxis, new Vector3(0, 0, 0));
                return false;
            }

            // calculate default (AABB) bounding box.
            // unlikely to be the best, but a good starting point.
            double minx = double.MaxValue;
            double maxx = double.MinValue;
            double miny = double.MaxValue;
            double maxy = double.MinValue;

            // track the vertices that touch our initial 4 edges.
            // e0 is "bottom" edge (+X)
            // e1 is "right" edge (+Y)
            // e2 is "top" edge (-X)
            // e3 is "left" edge (-Y)

            int ie0 = -1;
            int ie1 = -1;
            int ie2 = -1;
            int ie3 = -1;

            int c = vs.Length;
            for (int i = 0; i < c; i++)
            {
                var v = vs[i];
                if (v.X < minx) { minx = v.X; ie3 = i; }
                if (v.X > maxx) { maxx = v.X; ie1 = i; }
                if (v.Y < miny) { miny = v.Y; ie0 = i; }
                if (v.Y > maxy) { maxy = v.Y; ie2 = i; }
            }

            // we now have the AABB bounds, and the indices of the vertices that
            // touch each of those edges.
            var xd = maxx - minx;
            var yd = maxy - miny;

            var axis = Vector3.XAxis;
            var size = new Vector3(xd, yd);
            var centre = new Vector3(minx + (xd * 0.5f), miny + (yd * 0.5f));
            var bound_working = new Bound2D(centre, axis, size);
            bound = bound_working;

            // next, begin iterating.
            // I have four points (indices ie0-ie3) which touch the four
            // edges.

            double atot = 0.0f;
            int passes = 0;

            //ensure we don't spin around too much.
            while (atot < 90.0f && passes < (c * 4))
            {
                //Console.WriteLine($"pass {passes} : {ie0} {ie1} {ie2} {ie3}, angle {atot}");

                // generate the "following" edge vector (from the point) for each of these vertices.
                var fe0 = vs[(ie0 + 1) % c] - vs[ie0];
                var fe1 = vs[(ie1 + 1) % c] - vs[ie1];
                var fe2 = vs[(ie2 + 1) % c] - vs[ie2];
                var fe3 = vs[(ie3 + 1) % c] - vs[ie3];

                // now it's possible that using a normalized / dot with each of the respective axes
                // would be faster. HOWEVER, I need to use the angles later, so for now,
                // calculating the angles is the way to go.


                // check the size of the angle for each of these.
                var a0 = fe0.AngleTo(bound_working.axis_a);
                var a1 = fe1.AngleTo(bound_working.axis_b);
                var a2 = fe2.AngleTo(bound_working.axis_a.Negate());
                var a3 = fe3.AngleTo(bound_working.axis_b.Negate());

                // calculate minimum angle.
                var amin = Math.Min(a0, Math.Min(a1, Math.Min(a2, a3)));

                //Console.WriteLine($"{a0} {a1} {a2} {a3} : min {amin}");

                // advance one of the edge-touching indices.

                if (a0 == amin) ie0 = (ie0 + 1) % c;
                else if (a1 == amin) ie1 = (ie1 + 1) % c;
                else if (a2 == amin) ie2 = (ie2 + 1) % c;
                else ie3 = (ie3 + 1) % c;

                // could also re-calculate axis directly, rather than rotating - if we want less potential numerical error.
                // more complicated for ie1 thru ie3 as we're calculating axis_b, -axis_a and -axis_b respectively,
                // which need to be rotated to calculate axis_a into axis.
                /*
                if      (a0 == amin) { var ie0n = (ie0 + 1) % c; var edge = vs[ie0n] - vs[ie0]; axis = edge.normalized; ie0 = ie0n; }
                else if (a1 == amin) { var ie1n = (ie1 + 1) % c; var edge = vs[ie1n] - vs[ie1]; axis = new Vector3(edge.y, -edge.x).normalized; ie1 = ie1n; }
                else if (a2 == amin) { var ie2n = (ie2 + 1) % c; var edge = vs[ie2n] - vs[ie2]; axis = -edge.normalized; ie2 = ie2n; }
                else                 { var ie3n = (ie3 + 1) % c; var edge = vs[ie3n] - vs[ie3]; axis = new Vector3(-edge.y, edge.x).normalized; ie3 = ie3n; }
                 */

                if (amin > 0)
                {
                    axis = axis.Rotate(amin);
                    var axis_b = new Vector3(-axis.Y, axis.X);

                    // recalculate the bounds.
                    // ie0 touches axis_a (currently in axis) - on the "bottom" of the box.
                    // ie2 touches -axis_a - on the "top" of the box.
                    // ie1 touches axis_b - on the "right" of the box.
                    // ie3 touches -axis_b - on the "left" of the box.

                    var v0 = vs[ie0];
                    var v1 = vs[ie1];
                    var v2 = vs[ie2];
                    var v3 = vs[ie3];

                    var e02 = v2 - v0;        // vector from 0 to 2 
                    var e31 = v1 - v3;        // vector from 3 to 1 (note order!)

                    // project e02 and e31 onto the respective axes, to find the length of
                    // the edge.
                    xd = e31.Dot(axis);//Vector3.Dot(e31, axis);
                    yd = e02.Dot(axis_b); //Vector3.Dot(e02, axis_b);
                    size = new Vector3(xd, yd);

                    // calculate some corners, and average for the centre.
                    var q01 = v1 - v0;
                    var q03 = v3 - v0;
                    var c01 = (q01.Dot(axis) * axis) + v0; // (Vector3.Dot(q01, axis) * axis) + v0;
                    var c30 = (q03.Dot(axis) * axis) + v0;// (Vector3.Dot(q03, axis) * axis) + v0;

                    var q21 = v1 - v2;
                    var q23 = v3 - v2;
                    var c12 = (q21.Dot(axis.Negate()) * axis.Negate()) + v2; //(Vector3.Dot(q21, axis.Negate()) * axis.Negate()) + v2;
                    var c23 = (q23.Dot(axis.Negate()) * axis.Negate()) + v2; //(Vector3.Dot(q23, axis.Negate()) * axis.Negate()) + v2;


                    centre = (c01 + c12 + c23 + c30) * 0.25f;
                    bound_working = new Bound2D(centre, axis, size);

                    //Console.WriteLine($"new area : {bound_working.area}");

                    if (bound_working.area < bound.area)
                    {
                        //Console.WriteLine("better match.");
                        bound = bound_working;
                    }
                }
                else
                {
                    // angle is zero - a previous rotation has caused more than one edge
                    // to touch our bounding box. don't rotate, simply move that point
                    // up, along the second edge.
                    //Console.WriteLine("skipping a point ...");
                }

                // and on ...
                passes++;
                atot += amin;
            }

            // ensure bound is set such that longest axis is axis_a
            bound.AlignMajor();

            return true;
        }


        // Convert from Smallest Enclosing Rectangle to Centred Axis Aligned Bounding Box.
        // takes in a set of points (vs), and their Smallest Enclosing Rectangle bound.
        // rotates and centres the point set into vs_orient, and provides the new centred, axis-aligned
        // bound in caabb.
        //
        // REQUIRED - all points in vs fit inside ser.
        // REQUIRED - vs contains points such that a convex hull constructed from vs is minimally contained by ser.

        public static bool ConvertFromSERToCAABB(Vector3[] vs, Bound2D ser, out Vector3[] vs_orient, out Bound2D caabb)
        {

            // no-op output.
            vs_orient = vs;
            caabb = ser;

            // sanity checks, etc.
            if (vs.Length < 3) return false;
            if (ser.area <= 0) return false;
            // TODO - check ser contains vs
            // TODO - check ser minimally contains convex hull bounding vs

            var c = vs.Length;
            var centre = ser.centre;
            var angle = ser.angle;
            vs_orient = new Vector3[c];
            for (int i = 0; i < c; i++)
            {
                vs_orient[i] = (vs[i] - centre).Rotate(-angle);
            }
            caabb = new Bound2D(new Vector3(0, 0, 0), Vector3.XAxis, ser.size);

            return true;
        }


        // calculate axis-aligned cell array.
        // REQUIRED - ordered polygon points array in vs (does not need to be convex, DOES need to be ordered
        // as a simple polygon)

        public static bool CalculateInteriorCells(Vector3[] vs, out double[] xs, out double[] ys, out int[,] cells)
        {

            var vc = vs.Length;


            // optimisation : this is "pick distinct values, sort smallest to largest" - can definitely be
            // done faster for a couple of ms improvement.
            //xs = vs.Select(v => v.x).OrderBy(x => x).Distinct().ToArray();
            //ys = vs.Select(v => v.y).OrderBy(y => y).Distinct().ToArray();

            // this variant clamps with an epsilon, to ensure the math doesn't get too wacky.

            var xsl = new List<double>(vs.Length);
            var ysl = new List<double>(vs.Length);
            for (int i = 0; i < vs.Length; i++)
            {
                var v = vs[i];
                xsl.Add(v.X);
                ysl.Add(v.Y);
            }
            xsl.Sort();
            ysl.Sort();

            // de-dupe - including epsilon.

            double xmin = xsl[0];
            double xmax = xsl[xsl.Count - 1];
            double ymin = ysl[0];
            double ymax = ysl[ysl.Count - 1];
            double mmin = Math.Min(xmin, ymin);
            double mmax = Math.Max(xmax, ymax);

            var xsd = new List<double>(vs.Length) { xsl[0] };
            var ysd = new List<double>(vs.Length) { ysl[0] };

            double epsilon = (mmax - mmin) / (1024 * 1024);        // 1 millionth of the span.
            for (int i = 0; i < vs.Length - 1; i++)
            {
                if (xsl[i + 1] - xsl[i] > epsilon) xsd.Add(xsl[i + 1]);
                if (ysl[i + 1] - ysl[i] > epsilon) ysd.Add(ysl[i + 1]);
            }
            xs = xsd.ToArray();
            ys = ysd.ToArray();

            //Console.WriteLine("Xs :"); for (int i = 0; i < xs.Length; i++) Console.WriteLine($"\t{i}\t{xs[i]}");
            //Console.WriteLine("Ys :"); for (int i = 0; i < ys.Length; i++) Console.WriteLine($"\t{i}\t{ys[i]}");

            // if we want to, here we can extend xs / ys, by adding additional points (midpoints, or reproject
            // as per the paper)

            var xc = xs.Length - 1;
            var yc = ys.Length - 1;

            // cells array indicates if the cell at index X,Y (from xs[X] to xs[X+1] and ys[Y] to ys[Y+1])
            // is INTERIOR to the polygon (1) or EXTERIOR. (0 or -1)

            // 1D arrays perform faster than 2D arrays - if performance is critical, use a 1D array
            // and index into it, e.g. (y*xc) + x;
            cells = new int[xc, yc];

            // now, iterate the polygon edges. find the x span (and their associated indices in xs) and
            // the y span (and their associated indices in ys)

            var v0 = vs[0];
            var six = -1;
            var siy = -1;
            var eix = 0; while (xs[eix] < v0.X && eix < xs.Length - 1) eix++;
            var eiy = 0; while (ys[eiy] < v0.Y && eiy < ys.Length - 1) eiy++;

            for (int i = 0; i < vc; i++)
            {

                var s = vs[i];
                var e = vs[(i + 1) % vc];
                var edge = e - s;

                // get the indices for the start - it should be the end vertex of the previous edge.
                six = eix;
                siy = eiy;

                // eix = xs.IndexOf(e.x);
                // eiy = ys.IndexOf(e.y);
                // We can use the orientation of the edge to determine the direction of the search,
                // and the end vertex x,y indices should be near in the arrays to the start vertex x,y indices.
                // depending on edge length, etc - no guarantees!
                // could possibly binary search, but if the edge lengths are short, a linear scan should be
                // fairly fast anyway.
                int tx = edge.X >= 0 ? 1 : -1;             // -1 or 1
                int ty = edge.Y >= 0 ? 1 : -1;             // -1 or 1
                if (tx > 0) { while (xs[eix] < e.X && eix < xs.Length - 1) eix++; } else { while (xs[eix] > e.X && eix > 0) eix -= 1; }
                if (ty > 0) { while (ys[eiy] < e.Y && eiy < ys.Length - 1) eiy++; } else { while (ys[eiy] > e.Y && eiy > 0) eiy -= 1; }

                // we now have a span.
                var span_x_start = Math.Min(six, eix);
                var span_y_start = Math.Min(siy, eiy);

                var span_x_end = Math.Max(six, eix);
                var span_y_end = Math.Max(siy, eiy);

                //Console.WriteLine($"edge {i} : from ({s.x}, {s.y}) to ({e.x} {e.y}) - span x: {span_x_start}-{span_x_end}, y: {span_y_start}-{span_y_end}");

                // we care about the edge direction (and hence, the normal direction).

                // when testing a rectangle, we want a relevant corner to check is "inside" the edge.
                // if we pick the right one, we can be sure the other three are also inside - and any that 
                // fail mean the rectangle is partially or totally outside the edge.

                // given the winding order of the polygon is COUNTERCLOCKWISE, the "interior" side of the edge
                // is a 90 degree CCW rotation. (orthogonal to the edge, facing into the polygon)
                // I'm calling this "into" - note, this is NOT normalized! (we are only ever doing dot sign checks)
                var into = new Vector3(-edge.Y, edge.X);
                int rx = into.X >= 0 ? 0 : 1;
                int ry = into.Y >= 0 ? 0 : 1;

                // for -x,-y edges, the into direction is +X, -Y - use the TL cell vertex to test. (0,1)
                // for +x,-y edges, the into direction is +X, +Y - use the BL cell vertex to test. (0,0)
                // for -x,+y edges, the into direction is -X, -Y - use the TR cell vertex to test. (1,1)
                // for +x,+y edges, the into direction is -X, +Y - use the BR cell vertex to test. (1,0)

                // I'm not sure it's impossible for two edges to span the same cell - it may be, in which case
                // the "if interiors[p,q] < 0 continue" tests can be removed, below.

                // this covers all of the cells that the edge crosses.
                // if the edge is vertical or horizontal, the span width in that axis should be 0.
                if (span_x_end - span_x_start == 0)
                {
                    // vertical edge. pick the cells on the interior side.
                    var p = span_x_start - rx;
                    for (int q = span_y_start; q < span_y_end; q++)
                    {
                        if (cells[p, q] < 0) continue;
                        cells[p, q] = 1;
                    }
                    continue;
                }
                else if (span_y_end - span_y_start == 0)
                {
                    // horizontal edge. pick the cells on the interior side.
                    var q = span_y_start - ry;
                    for (int p = span_x_start; p < span_x_end; p++)
                    {
                        if (cells[p, q] < 0) continue;
                        cells[p, q] = 1;
                    }
                    continue;
                }
                for (int q = span_y_start; q < span_y_end; q++)
                {
                    for (int p = span_x_start; p < span_x_end; p++)
                    {
                        // if we've already marked this as exterior, then skip it.
                        // it's possible to be interior to another edge, but still exterior
                        // to this one, so continue the check in that case.
                        if (cells[p, q] < 0) continue;
                        // based on the edge direction, pick the correct corner to test against.
                        var v = new Vector3(xs[p + rx], ys[q + ry]);
                        var sv = v - s;
                        var d = sv.Dot(into);

                        // mark the cell either exterior (-1) or interior (1)
                        cells[p, q] = d < 0 ? -1 : 1;
                    }
                }
            }


            // outside region sweep.
            // some cells may not be spanned by edges.
            // start on the outside of the region, and mark everything as exterior, until
            // we come across a cell that has been explicitly marked.

            for (int q = 0; q < yc; q++)
            {
                // from the left edge.
                for (int x = 0; x < xc; x++)
                {
                    if (cells[x, q] != 0) break;
                    cells[x, q] = -1;
                }
                // from the right edge.
                for (int x = xc - 1; x >= 0; x--)
                {
                    if (cells[x, q] != 0) break;
                    cells[x, q] = -1;
                }
            }

            for (int p = 0; p < xc; p++)
            {
                // from the bottom edge.
                for (int y = 0; y < yc; y++)
                {
                    if (cells[p, y] != 0) break;
                    cells[p, y] = -1;
                }
                // from the top edge.
                for (int y = yc - 1; y >= 0; y--)
                {
                    if (cells[p, y] != 0) break;
                    cells[p, y] = -1;
                }
            }


            // sweep for interior (untested) cells.
            // in fact, mark everything.
            for (int j = 0; j < yc; j++)
            {
                for (int i = 0; i < xc; i++)
                {
                    // anything that was -1 goes to 0.
                    // anything that was 0 or 1 goes to 1.
                    // this ensures any un-tested cells are classed as interior.
                    cells[i, j] = cells[i, j] < 0 ? 0 : 1;
                }
            }

            return true;
        }


        // take the outputs from the preceeding function, and calculate the rectangular region that
        // has the largest area.

        static List<int> clir_hvec = new List<int>();
        static List<int> clir_vvec = new List<int>();
        static List<Int2> clir_spans = new List<Int2>();

        public static bool CalculateLargestInteriorRectangleUsingSpans(double[] xs, double[] ys, int[,] cells, out Bound2D best)
        {
            // cell lengths. interiors[x,y] should match [axc,ayc];
            int axc = xs.Length - 1;
            int ayc = ys.Length - 1;

            double best_area = 0.0f;
            Int2 best_origin = new Int2(-1, -1);
            Int2 best_span = new Int2(-1, -1);


            var adjacency_horizontal = new int[axc, ayc];
            var adjacency_vertical = new int[axc, ayc];

            // calculate horizontal adjacency, row by row
            for (int y = 0; y < ayc; y++)
            {
                int span = 0;
                for (int x = axc - 1; x >= 0; x--)
                {
                    if (cells[x, y] > 0) span++; else span = 0;
                    adjacency_horizontal[x, y] = span;
                }
            }

            // calculate vertical adjacency, column by column.
            for (int x = 0; x < axc; x++)
            {
                int span = 0;
                for (int y = ayc - 1; y >= 0; y--)
                {
                    if (cells[x, y] > 0) span++; else span = 0;
                    adjacency_vertical[x, y] = span;
                }
            }


            for (int y = 0; y < ayc; y++)
            {
                for (int x = 0; x < axc; x++)
                {
                    var iv = cells[x, y];
                    if (iv != 1) continue;

                    // generate H vector - this is horizontal adjacency for each step up.
                    clir_hvec.Clear();
                    // look at horizontal adjacency.
                    // step up from our initial cell, and look right.

                    var h = adjacency_horizontal[x, y];
                    clir_hvec.Add(h);
                    for (int q = y + 1; q < ayc; q++)
                    {
                        if (cells[x, q] != 1) break;
                        // each row can only be as large as the previous - a rectangle cannot push
                        // further out than a lower row.
                        h = Math.Min(adjacency_horizontal[x, q], h);
                        clir_hvec.Add(h);
                    }

                    // generate V vector. This is vertical adjacency for each step right.
                    clir_vvec.Clear();
                    // look at vertical adjacency.
                    // step right from our initial cell, and look up.

                    var v = adjacency_vertical[x, y];
                    clir_vvec.Add(v);
                    for (int p = x + 1; p < axc; p++)
                    {
                        if (cells[p, y] != 1) break;
                        // each column can only be as large as the previous - a rectangle cannot push
                        // further up than a previous column.
                        v = Math.Min(adjacency_vertical[p, y], v);
                        clir_vvec.Add(v);
                    }


                    // log the vectors.
                    // var hstr = string.Join(", ", hvec.Select(h => h.ToString()));
                    // var vstr = string.Join(", ", vvec.Select(h => h.ToString()));
                    // Console.WriteLine($"node ({x}, {y}) : H = ({hstr}), V = ({vstr})");

                    clir_spans.Clear();

                    // generate the set of valid spans.
                    Int2 span_last = new Int2(-1, -1);
                    for (int i = 0; i < clir_hvec.Count; i++)
                    {
                        int p = clir_hvec[i];
                        int q = clir_vvec[p - 1];
                        Int2 span = new Int2(p, q);
                        if (span.x != span_last.x && span.y != span_last.y)
                        {
                            clir_spans.Add(span);
                            span_last = span;
                        }
                    }


                    //Console.WriteLine($"SPANS FOR {x},{y} : {clir_spans.Count}");
                    //for (int i = 0; i < clir_spans.Count; i++)  Console.WriteLine($"\t{clir_spans[i]}");

                    // for each span, calculate the area.

                    for (int i = 0; i < clir_spans.Count; i++)
                    {
                        var span = clir_spans[i];
                        var xstart = xs[x];
                        var xend = xs[x + span.x];
                        var ystart = ys[y];
                        var yend = ys[y + span.y];
                        var xsize = xend - xstart;
                        var ysize = yend - ystart;
                        var area = xsize * ysize;
                        if (area > best_area)
                        {
                            best_area = area;
                            best_span = span;
                            best_origin = new Int2(x, y);
                        }
                    }
                }
            }

            if (best_area > 0)
            {
                //Console.WriteLine($"best area : {best_area} {best_origin} {best_span}");

                var xstart = xs[best_origin.x];
                var xend = xs[best_origin.x + best_span.x];
                var ystart = ys[best_origin.y];
                var yend = ys[best_origin.y + best_span.y];

                Vector3 centre = new Vector3((xend + xstart) * 0.5f, (yend + ystart) * 0.5f);
                Vector3 size = new Vector3((xend - xstart), (yend - ystart));

                //Console.WriteLine($"X : {xstart} {xend}, Y : {ystart} {yend}");
                Console.WriteLine($"Area : {best_area} Centre: {centre.F3()} size : {size.F3()}");

                best = new Bound2D(centre, Vector3.XAxis, size);
                return true;
            }

            best = new Bound2D(new Vector3(0, 0, 0), Vector3.XAxis, new Vector3(0, 0, 0));
            return false;
        }


        // take the outputs from the preceeding function, and calculate the rectangular region that
        // has the largest area.



        public static bool CalculateLargestInteriorRectangle(double[] xs, double[] ys, int[,] cells, out Bound2D best)
        {
            // cell lengths. interiors[x,y] should match [axc,ayc];
            int axc = xs.Length - 1;
            int ayc = ys.Length - 1;

            double best_area = 0.0f;
            Vector3 best_origin = new Vector3(1, 1, 1) * -1;
            Vector3 best_span = new Vector3(1, 1, 1) * -1;

            List<double> hspans = new List<double>();
            List<double> vspans = new List<double>();

            var lengths_horizontal = new double[axc, ayc];
            var lengths_vertical = new double[axc, ayc];


            for (int y = 0; y < ayc; y++)
            {
                double span = 0;
                for (int x = axc - 1; x >= 0; x--)
                {
                    span = (cells[x, y] <= 0) ? 0 : span + xs[x + 1] - xs[x];
                    lengths_horizontal[x, y] = span;
                }
            }

            for (int x = 0; x < axc; x++)
            {
                double span = 0;
                for (int y = ayc - 1; y >= 0; y--)
                {
                    span = (cells[x, y] <= 0) ? 0 : span + ys[y + 1] - ys[y];
                    lengths_vertical[x, y] = span;
                }
            }

            /* IF REQUIRED FOR TESTING, CALCULATE ADJACENCY.
            var adjacency_horizontal = new int[axc, ayc];
            var adjacency_vertical = new int[axc, ayc];
            for (int y = 0; y < ayc; y++)
            {
                int adj = 0;
                for (int x = axc - 1; x >= 0; x--)
                {
                    adj = (cells[x, y] <= 0) ? 0 : adj + 1;
                    adjacency_horizontal[x, y] = adj;
                }
            }

            for (int x = 0; x < axc; x++)
            {
                int adj = 0;
                for (int y = ayc - 1; y >= 0; y--)
                {
                    adj = (cells[x, y] <= 0) ? 0 : adj + 1;
                    adjacency_vertical[x, y] = adj;
                }
            }
            */





            for (int y = 0; y < ayc; y++)
            {
                for (int x = 0; x < axc; x++)
                {
                    var iv = cells[x, y];
                    if (iv == 0) continue;

                    var h = lengths_horizontal[x, y];
                    var v = lengths_vertical[x, y];

                    // if the best POSSIBLE area (which may not be valid!)
                    // is smaller than the best area, then we don't need to run any further tests.
                    if (h * v < best_area) continue;

                    // generate H vector - this is horizontal spans for each step up.
                    hspans.Clear();
                    // look at horizontal spans.
                    // step up from our initial cell, and look right.

                    hspans.Add(h);
                    for (int q = y + 1; q < ayc; q++)
                    {
                        if (cells[x, q] == 0) break;
                        var h2 = lengths_horizontal[x, q];
                        if (h2 >= h) continue;
                        h = h2;
                        hspans.Add(h);
                    }

                    // generate V vector. This is vertical spans for each step right.
                    vspans.Clear();
                    // look at vertical spans.
                    // step right from our initial cell, and look up.

                    vspans.Add(v);
                    for (int p = x + 1; p < axc; p++)
                    {
                        if (cells[p, y] == 0) break;
                        var v2 = lengths_vertical[p, y];
                        if (v2 >= v) continue;
                        v = v2;

                        vspans.Add(v);
                    }



                    // log the vectors.
                    //var hstr = string.Join(", ", hspans.Select(ht => ht.ToString()));
                    //var vstr = string.Join(", ", vspans.Select(vt => vt.ToString()));
                    //Console.WriteLine($"node ({x}, {y}) : H = ({hstr}), V = ({vstr})");

                    if (vspans.Count != hspans.Count)
                    {
                        Console.WriteLine($"span counts don't match.  {hspans.Count} {vspans.Count}");

                        /*

                        // run the adjacency code instead.
                        // THIS IS TESTING CODE. DO NOT USE THIS IN PRODUCTION - RATHER
                        // FIX THE ISSUE SO IT'S NOT REQUIRED.

                        // generate H vector - this is horizontal adjacency for each step up.
                        clir_hvec.Clear();
                        // look at horizontal adjacency.
                        // step up from our initial cell, and look right.

                        var sh = adjacency_horizontal[x, y];
                        clir_hvec.Add(sh);
                        for (int q = y + 1; q < ayc; q++)
                        {
                            if (cells[x, q] != 1) break;
                            // each row can only be as large as the previous - a rectangle cannot push
                            // further out than a lower row.
                            sh = Math.Min(adjacency_horizontal[x, q], sh);
                            clir_hvec.Add(sh);
                        }

                        // generate V vector. This is vertical adjacency for each step right.
                        clir_vvec.Clear();
                        // look at vertical adjacency.
                        // step right from our initial cell, and look up.

                        var sv = adjacency_vertical[x, y];
                        clir_vvec.Add(sv);
                        for (int p = x + 1; p < axc; p++)
                        {
                            if (cells[p, y] != 1) break;
                            // each column can only be as large as the previous - a rectangle cannot push
                            // further up than a previous column.
                            sv = Math.Min(adjacency_vertical[p, y], sv);
                            clir_vvec.Add(sv);
                        }

                        clir_spans.Clear();

                        // generate the set of valid spans.
                        Int2 span_last = new Int2(-1, -1);
                        for (int i = 0; i < clir_hvec.Count; i++)
                        {
                            int p = clir_hvec[i];
                            int q = clir_vvec[p - 1];
                            Int2 span = new Int2(p, q);
                            if (span.x != span_last.x && span.y != span_last.y)
                            {
                                clir_spans.Add(span);
                                span_last = span;
                            }
                        }

                        for (int i = 0; i < clir_spans.Count; i++)
                        {
                            var span = clir_spans[i];
                            var xstart = xs[x];
                            var xend = xs[x + span.x];
                            var ystart = ys[y];
                            var yend = ys[y + span.y];
                            var xsize = xend - xstart;
                            var ysize = yend - ystart;
                            var area = xsize * ysize;
                        }

                        */

                        continue;
                    }

                    // reverse the v spans list - this lets us trivially combine the correct
                    // spans for each rectangle combination with the same list index.
                    vspans.Reverse();

                    for (int i = 0; i < hspans.Count; i++)
                    {
                        double hl = hspans[i];
                        double vl = vspans[i];
                        double area = hl * vl;
                        if (area > best_area)
                        {
                            best_area = area;
                            best_origin = new Vector3(xs[x], ys[y]);
                            best_span = new Vector3(hl, vl);
                        }
                    }

                }
            }

            if (best_area > 0)
            {
                Vector3 centre = best_origin + (best_span * 0.5f);
                Vector3 size = best_span;
                Console.WriteLine($"best area : {best_area} {centre.F3()} {size.F3()}");

                //Console.WriteLine($"X : {xstart} {xend}, Y : {ystart} {yend}");
                //Console.WriteLine($"Centre: {centre.F3()} size : {size.F3()}");
                best = new Bound2D(centre, Vector3.XAxis, size);
                return true;
            }

            best = new Bound2D(new Vector3(0, 0, 0), Vector3.XAxis, new Vector3(0, 0, 0));
            return false;
        }

        /// <summary>
        /// Compute the largest interior XY-oriented rectangle in the specified polygon
        /// </summary>
        /// <param name="p">The 2D polygon to test</param>
        /// <param name="best">The best result, as a Bound2D.</param>
        /// <returns></returns>
        public static bool CalculateLargestInteriorRectangle(Polygon p, out Bound2D best)
        {
            best = null;
            var interiorCellsSucceeded = CalculateInteriorCells(p.Vertices.ToArray<Vector3>(), out var xs, out var ys, out int[,] cells);
            return interiorCellsSucceeded && CalculateLargestInteriorRectangle(xs, ys, cells, out best);
        }

        /// <summary>
        /// Compute the largest interior rectangle in the specified polygon,
        /// testing many different angles (in a 90-degree sweep) to achieve the best result.
        /// </summary>
        /// <param name="p">The 2D polygon to test</param>
        /// <param name="angleStep">The angle step size in degrees to use. The smaller the angle step, the longer this may take.</param>
        /// <param name="best">The best result, as a Bound2D.</param>
        public static bool CalculateLargestInteriorRectangleWithAngleSweep(Polygon p, double angleStep, out Bound2D best)
        {
            return CalculateLargestInteriorRectangleWithAngleSweep(p.Vertices.ToArray<Vector3>(), angleStep, out best, out _);
        }

        // this version does ONE lot of allocations, and tests across angles too.
        // scans from -45 to +45 degrees in angle_step increments.
        // optionally (see use_fan) splits the angle search up to give wider
        // coverage of angle options earlier, so quick reject gets more opportunity
        // to work.

        public static bool CalculateLargestInteriorRectangleWithAngleSweep(Vector3[] vs_src, double angle_step, out Bound2D best, out double best_angle)
        {
            best = new Bound2D(new Vector3(0, 0, 0), Vector3.XAxis, new Vector3(0, 0, 0));
            best_angle = 0;

            int vc = vs_src.Length;
            if (vc < 4)
            {
                return false;
            }
            Vector3[] vs = new Vector3[vc];

            double[] xs = new double[vc];
            double[] ys = new double[vc];

            double[] xds = new double[vc];
            double[] yds = new double[vc];

            List<double> xsl = new List<double>(vc);
            List<double> ysl = new List<double>(vc);

            double[] hspans = new double[vc];
            double[] vspans = new double[vc];

            var lengths_horizontal = new double[vc, vc];
            var lengths_vertical = new double[vc, vc];

            var cells = new int[vc, vc];

            double best_area = 0;
            best_angle = 0;
            Vector3 best_origin = new Vector3(1, 1, 1) * -1;
            Vector3 best_span = new Vector3(1, 1, 1) * -1;

            double ang_min = -45;
            double ang_max = 45;
            double ang_range = ang_max - ang_min;

            List<double> steps_lin = new List<double>();
            for (double f = ang_min; f < ang_max; f += angle_step) steps_lin.Add(f);

            var steps = steps_lin;

            bool use_fan = true;
            if (use_fan)
            {

                // divide and conquer the angle steps - should hopefully
                // hone in on a large value faster than just a linear sweep,
                // which should give faster quick reject.

                // ok - that's all of our steps in an array.
                // next, start pulling them out in layers (POT)
                // so we get 1, 2, 4, 8, 16, etc.

                int step_count = steps_lin.Count;
                List<int> steps_at_layer = new List<int>();
                int pull = 1;
                int steps_left = step_count;
                while (steps_left > 0)
                {
                    steps_at_layer.Add(pull);
                    steps_left -= pull;
                    pull = Math.Min(steps_left, pull << 1);
                }


                List<double> steps_fan = new List<double>();
                foreach (var lay in steps_at_layer)
                {
                    // pull out lay items from steps_lin.
                    double range_per_item = ang_range / (double)lay;
                    double mid_per_item = range_per_item * 0.5f;
                    double next = ang_min;
                    for (int i = 0; i < lay; i++)
                    {
                        double check = next + mid_per_item;

                        int idx = steps_lin.BinarySearch(check);
                        if (idx < 0) idx = Clamp((~idx) + 1, 0, steps_lin.Count - 1);
                        var val = steps_lin[idx];
                        steps_fan.Add(val);
                        steps_lin.RemoveAt(idx);
                        next += range_per_item;
                    }
                }
                steps = steps_fan;
            }


            // work out the number of layers, where we have 1 right at the top.

            foreach (var angle in steps)
            //double angle = 0f;
            {
                //Console.WriteLine(angle);

                xsl.Clear();
                ysl.Clear();

                for (int i = 0; i < vc; i++)
                {
                    var v = vs_src[i].Rotate(angle);
                    vs[i] = v;
                    xsl.Add(v.X);
                    ysl.Add(v.Y);
                }

                xsl.Sort();
                ysl.Sort();

                double xmin = xsl[0];
                double xmax = xsl[xsl.Count - 1];
                double ymin = ysl[0];
                double ymax = ysl[ysl.Count - 1];
                double mmin = Math.Min(xmin, ymin);
                double mmax = Math.Max(xmax, ymax);
                double epsilon = (mmax - mmin) / (1024 * 1024);        // 1 millionth of the span.

                var xc = 1;
                var yc = 1;
                xs[0] = xsl[0];
                ys[0] = ysl[0];

                for (int i = 1; i < vs.Length; i++)
                {
                    if (xsl[i] - xsl[i - 1] > epsilon) xs[xc++] = xsl[i];
                    if (ysl[i] - ysl[i - 1] > epsilon) ys[yc++] = ysl[i];
                }

                for (int i = 0; i < xc - 1; i++) xds[i] = xs[i + 1] - xs[i];
                for (int i = 0; i < yc - 1; i++) yds[i] = ys[i + 1] - ys[i];

                xc--;
                yc--;


                System.Array.Clear(cells, 0, cells.Length);

                var v0 = vs[0];

                var eix = 0; while (xs[eix] < v0.X && eix < xc) eix++;
                var eiy = 0; while (ys[eiy] < v0.Y && eiy < yc) eiy++;
                // in this case, it may be slightly faster - but not worth the algo complexity.
                //var eix = System.Array.BinarySearch(xs, v0.x); if (eix < 0) eix = Math.Clamp(0,xc+1,(~eix) - 1);
                //var eiy = System.Array.BinarySearch(ys, v0.y); if (eix < 0) eiy = Math.Clamp(0,yc+1,(~eiy) - 1);


                for (int i = 0; i < vc; i++)
                {
                    var s = vs[i];
                    var e = vs[(i + 1) % vc];
                    var edge = e - s;

                    var into = new Vector3(-edge.Y, edge.X);
                    int rx = into.X >= 0 ? 0 : 1;
                    int ry = into.Y >= 0 ? 0 : 1;

                    var six = eix;
                    var siy = eiy;


                    int tx = edge.X >= 0 ? 1 : -1;             // -1 or 1
                    int ty = edge.Y >= 0 ? 1 : -1;             // -1 or 1

                    if (tx > 0) { while (xs[eix] < e.X && eix <= xc) eix++; } else { while (xs[eix] > e.X && eix > 0) eix--; }
                    if (ty > 0) { while (ys[eiy] < e.Y && eiy <= yc) eiy++; } else { while (ys[eiy] > e.Y && eiy > 0) eiy--; }

                    // tried it. timed it. much slower. probably because of locality in the whiles, above.
                    //int wx = edge.x >= 0 ? 0 : 1;
                    //int wy = edge.y >= 0 ? 0 : 1;
                    //var eixb = System.Array.BinarySearch(xs, e.x); if (eixb < 0) { eixb = Math.Clamp((~eixb) - wx, 0, xc + 1); }
                    //var eiyb = System.Array.BinarySearch(ys, e.y); if (eiyb < 0) { eiyb = Math.Clamp((~eiyb) - wy, 0, yc + 1); }



                    // we now have a span.
                    var span_x_start = Math.Min(six, eix);
                    var span_y_start = Math.Min(siy, eiy);

                    var span_x_end = Math.Max(six, eix);
                    var span_y_end = Math.Max(siy, eiy);



                    if (span_x_end - span_x_start == 0)
                    {
                        var p = span_x_start - rx;
                        for (int q = span_y_start; q < span_y_end; q++)
                        {
                            if (cells[p, q] < 0) continue;
                            cells[p, q] = 1;
                        }
                        continue;
                    }
                    else if (span_y_end - span_y_start == 0)
                    {
                        var q = span_y_start - ry;
                        for (int p = span_x_start; p < span_x_end; p++)
                        {
                            if (cells[p, q] < 0) continue;
                            cells[p, q] = 1;
                        }
                        continue;
                    }
                    for (int q = span_y_start; q < span_y_end; q++)
                    {
                        for (int p = span_x_start; p < span_x_end; p++)
                        {
                            if (cells[p, q] < 0) continue;
                            var v = new Vector3(xs[p + rx], ys[q + ry]);
                            var sv = v - s;
                            var d = sv.Dot(into);
                            cells[p, q] = d < 0 ? -1 : 1;
                        }
                    }
                }

                // outside region sweep.
                // some cells may not be spanned by edges.
                // start on the outside of the region, and mark everything as exterior, until
                // we come across a cell that has been explicitly marked.

                for (int q = 0; q < yc; q++)
                {
                    for (int x = 0; x < xc; x++)
                    {
                        if (cells[x, q] != 0) break;
                        cells[x, q] = -1;
                    }
                    for (int x = xc - 1; x >= 0; x--)
                    {
                        if (cells[x, q] != 0) break;
                        cells[x, q] = -1;
                    }
                }


                for (int p = 0; p < xc; p++)
                {
                    for (int y = 0; y < yc; y++)
                    {
                        if (cells[p, y] != 0) break;
                        cells[p, y] = -1;
                    }
                    for (int y = yc - 1; y >= 0; y--)
                    {
                        if (cells[p, y] != 0) break;
                        cells[p, y] = -1;
                    }
                }

                // sweep for interior (untested) cells.
                // in fact, mark everything.
                for (int j = 0; j < yc; j++)
                {
                    for (int i = 0; i < xc; i++)
                    {
                        // anything that was -1 goes to 0.
                        // anything that was 0 or 1 goes to 1.
                        // this ensures any un-tested cells are classed as interior.
                        cells[i, j] = cells[i, j] < 0 ? 0 : 1;
                    }
                }

                for (int y = 0; y < yc; y++)
                {
                    double span = 0;
                    for (int x = xc - 1; x >= 0; x--)
                    {
                        span = (cells[x, y] <= 0) ? 0 : span + xds[x];
                        lengths_horizontal[x, y] = span;
                    }
                }

                for (int x = 0; x < xc; x++)
                {
                    double span = 0;
                    for (int y = yc - 1; y >= 0; y--)
                    {
                        span = (cells[x, y] <= 0) ? 0 : span + yds[y];
                        lengths_vertical[x, y] = span;
                    }
                }

                for (int y = 0; y < yc; y++)
                {
                    for (int x = 0; x < xc; x++)
                    {
                        var iv = cells[x, y];
                        if (iv == 0) continue;


                        var h = lengths_horizontal[x, y];
                        var v = lengths_vertical[x, y];
                        if (h * v < best_area) continue;

                        int hsc = 1;
                        hspans[0] = h;
                        for (int q = y + 1; q < yc; q++)
                        {
                            if (cells[x, q] == 0) break;
                            var h2 = lengths_horizontal[x, q];
                            if (h2 >= h) continue;
                            h = h2;
                            hspans[hsc++] = h;
                        }

                        int vsc = 1;
                        vspans[0] = v;
                        for (int p = x + 1; p < xc; p++)
                        {
                            if (cells[p, y] == 0) break;
                            var v2 = lengths_vertical[p, y];
                            if (v2 >= v) continue;
                            v = v2;
                            vspans[vsc++] = v;
                        }

                        if (hsc != vsc)
                        {
                            Console.WriteLine($"span counts don't match.  {hsc} {vsc}");
                            continue;
                        }

                        for (int i = 0; i < hsc; i++)
                        {
                            double hl = hspans[i];
                            double vl = vspans[hsc - (i + 1)];
                            double area = hl * vl;
                            if (area > best_area)
                            {
                                best_area = area;
                                best_origin = new Vector3(xs[x], ys[y]);
                                best_span = new Vector3(hl, vl);
                                best_angle = angle;
                            }
                        }
                    }
                }
            }

            //-------------------------------
            if (best_area <= 0)
            {
                return false;

            }

            var bo_r = best_origin.Rotate(-best_angle);
            var bs_r = best_span.Rotate(-best_angle);

            Vector3 centre = bo_r + (bs_r * 0.5f);
            Vector3 size = best_span;
            Vector3 axis = Vector3.XAxis.Rotate(-best_angle);

            best = new Bound2D(centre, axis, size);

            return true;
        }

        private static double Clamp(double value, double min, double max)
        {
            if (value < min) return min;
            if (value > max) return max;
            return value;
        }

        private static int Clamp(int value, int min, int max)
        {
            if (value < min) return min;
            if (value > max) return max;
            return value;
        }
    }


}