/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package didroe.tnc;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Test the Tnc class.
 *
 * The tests taken from SciPy are originally sourced from: Prof. K.
 * Schittkowski's test examples for constrained non-linear programming.
 * http://www.ai7.uni-bayreuth.de/tpnp08.htm
 *
 * @author Did
 */
public class TncTest {

    /**
     * Assert that evaluating the function with the expected and actual 
     * parameters gives the same value.
     */
    private void assertExpectedResult(TncFunction function, TncResult result, 
            double[] expectedX) {
        double[] dummy = new double[expectedX.length];
        assertEquals(function.evaluate(expectedX, dummy), 
                function.evaluate(result.parameters(), dummy), 1e-8);
    }
    
    /**
     * @return Test setup for SciPy tests 1 and 2
     */
    private TncFunction scipy1And2Function() {
        return (double[] x, double[] gradient) -> {
            double a = 100.0;

            gradient[1] = 2.0 * a * (x[1] - Math.pow(x[0], 2));
            gradient[0] = -2.0 * (x[0] * (gradient[1] - 1.0) + 1.0);

            return (a * Math.pow((x[1] - Math.pow(x[0], 2)), 2))
                    + Math.pow(1.0 - x[0], 2);
        };
    }

    @Test
    public void testScipy1() {
        TncFunction function = scipy1And2Function();
        TncResult result = new Tnc()
                .maxFunctionEvaluations(200)
                .function(function)
                .initialGuess(new double[]{ -2.0, 1.0 })
                .lowerBounds(new double[]{ Double.NEGATIVE_INFINITY, -1.5 })
                .minimize();
        assertExpectedResult(function, result, new double[]{ 1.0, 1.0 });
    }

    @Test
    public void testScipy2() {
        TncFunction function = scipy1And2Function();
        TncResult result = new Tnc()
                .maxFunctionEvaluations(200)
                .function(function)
                .initialGuess(new double[]{ -2.0, 1.0 })
                .lowerBounds(new double[]{ Double.NEGATIVE_INFINITY, 1.5 })
                .minimize();
        assertExpectedResult(function, result, new double[]{ -1.2210262419616387, 1.5 });
    }

    @Test
    public void testScipy3() {
        TncFunction function = (double[] x, double[] gradient) -> {
            gradient[0] = -2.0 * (x[1] - x[0]) * 1.0e-5;
            gradient[1] = 1.0 - gradient[0];

            return x[1] + (Math.pow(x[1] - x[0], 2) * 1.0e-5);
        };

        TncResult result = new Tnc()
                .maxFunctionEvaluations(200)
                .function(function)
                .initialGuess(new double[]{ 10.0, 1.0 })
                .lowerBounds(new double[]{ Double.NEGATIVE_INFINITY, 0.0 })
                .minimize();
        assertExpectedResult(function, result, new double[]{ 0.0, 0.0 });
    }

    @Test
    public void testScipy4() {
        TncFunction function = (double[] x, double[] gradient) -> {
            gradient[0] = Math.pow(x[0] + 1.0, 2);
            gradient[1] = 1.0;

            return (Math.pow(x[0] + 1.0, 3) / 3.0) + x[1];
        };

        TncResult result = new Tnc()
                .maxFunctionEvaluations(200)
                .function(function)
                .initialGuess(new double[]{ 1.125, 0.125 })
                .lowerBounds(new double[]{ 1.0, 0.0 })
                .minimize();
        assertExpectedResult(function, result, new double[]{ 1.0, 0.0 });
    }

    @Test
    public void testScipy5() {
        TncFunction function = (double[] x, double[] gradient) -> {
            double v1 = Math.cos(x[0] + x[1]);
            double v2 = 2.0 * (x[0] - x[1]);
            gradient[0] = v1 + v2 - 1.5;
            gradient[1] = v1 - v2 + 2.5;

            return Math.sin(x[0] + x[1]) + Math.pow(x[0] - x[1], 2)
                    - (1.5 * x[0]) + (2.5 * x[1]) + 1.0;
        };

        TncResult result = new Tnc()
                .maxFunctionEvaluations(200)
                .function(function)
                .initialGuess(new double[]{ 0.0, 0.0 })
                .lowerBounds(new double[]{ -1.5, -3.0 })
                .upperBounds(new double[]{ 4.0, 3.0 })
                .minimize();
        assertExpectedResult(function, result, 
                new double[]{ -0.54719755119659763, -1.5471975511965976 });
    }

    @Test
    public void testScipy38() {
        TncFunction function = (double[] x, double[] gradient) -> {
            gradient[0] = ((-400.0 * x[0] * (x[1] - Math.pow(x[0], 2))) 
                    - (2.0 * (1.0 - x[0]))) * 1.0e-5;
            gradient[1] = ((200.0 * (x[1] - Math.pow(x[0], 2)))
                    + (20.2 * (x[1] - 1.0)) + (19.8 * (x[3] - 1.0))) * 1.0e-5;
            gradient[2] = ((-360.0 * x[2] * (x[3] - Math.pow(x[2], 2)))
                    - (2.0 * (1.0 - x[2]))) * 1.0e-5;
            gradient[3] = ((180.0 * (x[3] - Math.pow(x[2], 2)))
                    + (20.2 * (x[3] - 1.0)) + (19.8 * (x[1] - 1.0))) * 1.0e-5;

            return ((100.0 * Math.pow(x[1] - Math.pow(x[0], 2), 2)) + Math.pow(1.0 - x[0], 2)
                    + (90.0 * Math.pow(x[3] - Math.pow(x[2], 2), 2)) + Math.pow(1.0 - x[2], 2)
                    + (10.1 * (Math.pow(x[1] - 1.0, 2) + Math.pow(x[3] - 1.0, 2)))
                    + (19.8 * (x[1] - 1.0) * (x[3] - 1.0))) * 1.0e-5;
        };

        TncResult result = new Tnc()
                .maxFunctionEvaluations(200)
                .function(function)
                .initialGuess(new double[]{ -3.0, -1.0, -3.0, -1.0 })
                .lowerBounds(new double[]{ -10.0, -10.0, -10.0, -10.0 })
                .upperBounds(new double[]{ 10.0, 10.0, 10.0, 10.0 })
                .minimize();
        assertExpectedResult(function, result, new double[]{ 1.0, 1.0, 1.0, 1.0 });
    }

    @Test
    public void testScipy45() {
        TncFunction function = (double[] x, double[] gradient) -> {
            gradient[0] = -x[1] * x[2] * x[3] * x[4] / 120.0;
            gradient[1] = -x[0] * x[2] * x[3] * x[4] / 120.0;
            gradient[2] = -x[0] * x[1] * x[3] * x[4] / 120.0;
            gradient[3] = -x[0] * x[1] * x[2] * x[4] / 120.0;
            gradient[4] = -x[0] * x[1] * x[2] * x[3] / 120.0;

            return 2.0 - (x[0] * x[1] * x[2] * x[3] * x[4] / 120.0);
        };

        TncResult result = new Tnc()
                .maxFunctionEvaluations(200)
                .function(function)
                .initialGuess(new double[]{ 2.0, 2.0, 2.0, 2.0, 2.0 })
                .lowerBounds(new double[]{ 0.0, 0.0, 0.0, 0.0, 0.0 })
                .upperBounds(new double[]{ 1.0, 2.0, 3.0, 4.0, 5.0 })
                .minimize();
        assertExpectedResult(function, result, new double[]{ 1.0, 2.0, 3.0, 4.0, 5.0 });
    }
}
