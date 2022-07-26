using System.Collections;
using System.Collections.Generic;
using Elements.Geometry;

// 2D bounds.
// NOT XY axis aligned.

namespace Evryway
{

    public class Bound2D
    {
        public Vector3 centre { get; private set; }
        public Vector3 axis_a { get; private set; }     // first axis (unit length)
        public Vector3 axis_b { get; private set; }     // second axis (unit length)
        public double length_a { get; private set; }     // first axis length
        public double length_b { get; private set; }     // second axis length (may be larger)
        public Vector3 size { get; private set; }       // (length_a,length_b)
        public Vector3 extents { get => size * 0.5f; }    // half-length, e.g. from centre.

        public Vector3[] corners { get; private set; }  // four corners. in order: BL, BR, TR, TL (assuming axis a is X axis and axis b is Y axis)
                                                        // (- axis_a * extents.X - axis_b * extents.y), 
                                                        // (+ axis_a * extents.x - axis_b * extents.y), 
                                                        // (+ axis_a * extents.x + axis_b * extents.y), 
                                                        // (- axis_a * extents.x + axis_b * extents.y), 
        public double area { get; private set; }
        public bool major_axis_is_a { get => length_a >= length_b; }
        public double angle { get; private set; }

        public Vector3 bl { get => corners[0]; }
        public Vector3 br { get => corners[1]; }
        public Vector3 tr { get => corners[2]; }
        public Vector3 tl { get => corners[3]; }

        public Bound2D(Vector3 centre, Vector3 axis, Vector3 size)
        {
            this.centre = centre;
            this.axis_a = axis.Unitized();
            this.axis_b = new Vector3(-axis_a.Y, axis_a.X);
            this.size = size;
            this.length_a = size.X;
            this.length_b = size.Y;
            Cache();

        }

        // make the bound aligned such that the major axis (longest length) is axis_a, and
        // axis_a.x is positive.

        public void AlignMajor()
        {
            bool dirty = false;
            if (!major_axis_is_a)
            {
                // rotate the axes.
                var axis_t = axis_a;
                axis_a = axis_b;
                axis_b = axis_t.Negate();
                // swap lengths.
                length_a = size.Y;
                length_b = size.X;
                size = new Vector3(length_a, length_b);
                dirty = true;
            }

            if (axis_a.X <= 0)
            {
                axis_a = axis_a.Negate();
                axis_b = axis_b.Negate();
                dirty = true;
            }

            if (dirty)
            {
                Cache();
            }

        }

        void Cache()
        {
            area = length_a * length_b;
            angle = Vector3.XAxis.PlaneAngleTo(axis_a);
            corners = new Vector3[]
            {
                centre - (axis_a * extents.X) - (axis_b * extents.Y),
                centre + (axis_a * extents.X) - (axis_b * extents.Y),
                centre + (axis_a * extents.X) + (axis_b * extents.Y),
                centre - (axis_a * extents.X) + (axis_b * extents.Y),
            };
        }

        public string info { get => $"area: {area}, axes: ({axis_a.X}, {axis_a.Y}) ; ({axis_b.X}, {axis_b.Y}), size: ({length_a}, {length_b}) angle: {angle}"; }

    }
}
