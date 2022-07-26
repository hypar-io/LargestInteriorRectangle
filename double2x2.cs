namespace Evryway
{

    public class Double2
    {
        public double x { get; set; }
        public double y { get; set; }

        public Double2(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
    }
    public class Double2x2
    {


        public double m00 { get; set; }
        public double m01 { get; set; }
        public double m10 { get; set; }
        public double m11 { get; set; }
        public Double2x2(double m00, double m01, double m10, double m11)
        {
            this.m00 = m00;
            this.m01 = m01;
            this.m10 = m10;
            this.m11 = m11;
        }

        public Double2 c0 => new Double2(m00, m01);

        public Double2 c1 => new Double2(m10, m11);
    }

}