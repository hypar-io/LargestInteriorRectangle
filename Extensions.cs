using System;
using Elements.Geometry;

namespace Evryway
{

    // utility functions, extending Vector3.

    public static class Vector3Extensions
    {
        private const double Deg2Rad = Math.PI / 180.0;
        // rotate a Vector3. positive is Counter-clockwise.
        public static Vector3 Rotate(this Vector3 v, double angle_in_degrees)
        {
            double rad = angle_in_degrees * Deg2Rad;
            double cosa = Math.Cos(rad);
            double sina = Math.Sin(rad);
            var rx = (v.X * cosa) - (v.Y * sina);
            var ry = (v.X * sina) + (v.Y * cosa);
            return new Vector3(rx, ry);
        }

        // rotate a Vector3 array. positive is Counter-clockwise.
        public static void Rotate(this Vector3[] array, double angle_in_degrees)
        {
            double rad = angle_in_degrees * Deg2Rad;
            double cosa = Math.Cos(rad);
            double sina = Math.Sin(rad);

            for (int i = 0; i < array.Length; i++)
            {
                var v = array[i];
                var rx = (v.X * cosa) - (v.Y * sina);
                var ry = (v.X * sina) + (v.Y * cosa);
                array[i] = new Vector3(rx, ry);
            }
        }



        public static string F(this Vector3 source, int p) => source.ToString(); //$"F{p}");
        public static string F3(this Vector3 source) => source.ToString(); //("F3");
    }
}
