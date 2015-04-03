/* 
 * See the "LICENSE" file for the full license governing this code. 
 */
package didroe.tnc;

/**
 * General maths functions for arrays
 *
 * @author Did
 */
interface ArrayMath {

    /**
     * Ensure each element of x is within the given bounds.
     *
     * @param x The array to be modified
     * @param lowerBounds An array of lower bounds for each x element
     * @param upperBounds An array of upper bounds for each x element
     */
    public static void clip(double[] x, double[] lowerBounds, double[] upperBounds) {
        for (int i = 0; i < x.length; i++) {
            if (x[i] < lowerBounds[i]) {
                x[i] = lowerBounds[i];
            } else if (x[i] > upperBounds[i]) {
                x[i] = upperBounds[i];
            }
        }
    }

    /**
     * Dot product of vectors x and y. Assumes that x and y are of the same
     * length.
     */
    public static double dotProduct(double[] x, double[] y) {
        double result = 0.0;
        for (int i = 0; i < x.length; i++) {
            result += x[i] * y[i];
        }
        return result;
    }

    /**
     * Copy the source array elements into destination array. Assumes the source
     * and destination are of the same length.
     */
    public static void copy(double[] src, double[] dest) {
        System.arraycopy(src, 0, dest, 0, src.length);
    }

    /**
     * Updates each element of y subject to y += a*x for each corresponding x
     * value. Assumes x and y have the same length.
     */
    public static void axPlusY(double a, double[] x, double[] y) {
        for (int i = 0; i < x.length; i++) {
            y[i] += a * x[i];
        }
    }

    /**
     * Updates each element of y subject to y += x for each corresponding x
     * value. Assumes x and y have the same length.
     */
    public static void xPlusY(double[] x, double[] y) {
        for (int i = 0; i < x.length; i++) {
            y[i] += x[i];
        }
    }

    /**
     * The euclidian norm of x. ie. sqrt(x1^2 + x2^2 + ...)
     */
    public static double euclidianNorm(double[] x) {
        // To avoid overflow, this algorithm scales each element and undoes the
        // scaling at the end. The scale is adjusted incrementally as elements 
        // are processed
        double sumOfSquares = 1.0;
        double scale = 0.0;

        for (int i = 0; i < x.length; i++) {
            double element = x[i];
            if (element != 0.0) {
                double absElement = Math.abs(element);
                if (scale < absElement) {
                    // Scale needs updating
                    double ratio = scale / absElement;
                    scale = absElement;
                    sumOfSquares = (sumOfSquares * ratio * ratio) + 1.0;
                } else {
                    // Current scale is ok
                    double ratio = absElement / scale;
                    sumOfSquares += ratio * ratio;
                }
            }
        }

        return scale * Math.sqrt(sumOfSquares);
    }

    /**
     * Negate all elements of x
     */
    public static void negate(double[] x) {
        for (int i = 0; i < x.length; i++) {
            x[i] = -x[i];
        }
    }

}
