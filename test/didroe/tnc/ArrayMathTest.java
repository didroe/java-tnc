/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package didroe.tnc;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Tests for ArrayMath class
 * @author Did
 */
public class ArrayMathTest {
    
    @Test
    public void testClip() {
        double[] lowerBounds = new double[] { 1.0, 2.0 };
        double[] upperBounds = new double[] { 8.0, 9.0 };

        double[] x1 = new double[] { 0.0, 0.0 };
        ArrayMath.clip(x1, lowerBounds, upperBounds);        
        assertArrayEquals(lowerBounds, x1, 0);        

        double[] x2 = new double[] { 10.0, 10.0 };
        ArrayMath.clip(x2, lowerBounds, upperBounds);        
        assertArrayEquals(upperBounds, x2, 0);        

        double[] x3 = new double[] { 5.0, 6.0 };
        ArrayMath.clip(x3, lowerBounds, upperBounds);        
        assertArrayEquals(new double[] { 5.0, 6.0 }, x3, 0);
    }

    @Test
    public void testDotProduct() {
        double[] x = new double[] { 1.0, 2.0 };
        double[] y = new double[] { 3.0, 4.0 };
        assertEquals(11.0, ArrayMath.dotProduct(x, y), 0);        
    }
    
    @Test
    public void testCopy() {
        double[] x = new double[] { 1.0, 2.0 };
        double[] y = new double[x.length];
        ArrayMath.copy(x, y);
        assertArrayEquals(new double[] { 1.0, 2.0 }, y, 0);
    }

    @Test
    public void testAxPlusY() {
        double a = 5.0;
        double[] x = new double[] { 6.0, 8.0 };
        double[] y = new double[] { 1.0, 2.0 };
        ArrayMath.axPlusY(a, x, y);
        assertArrayEquals(new double[] { 31.0, 42.0 }, y, 0);    
    }

    @Test
    public void testXplusY() {
        double[] x = new double[] { 6.0, 8.0 };
        double[] y = new double[] { 1.0, 2.0 };
        ArrayMath.xPlusY(x, y);
        assertArrayEquals(new double[] { 7.0, 10.0 }, y, 0);    
    }
    
    @Test
    public void testEuclideanNorm() {
        double[] x = new double[] { 4.0, 3.0, 0.0 };
        assertEquals(5.0, ArrayMath.euclidianNorm(x), 0);    
    }
    
    @Test
    public void testNegate() {
        double[] x = new double[] { 2.0, 4.0 };
        ArrayMath.negate(x);
        assertArrayEquals(new double[] { -2.0, -4.0 }, x, 0);    
    }
    
}
