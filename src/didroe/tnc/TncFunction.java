/* 
 * See the "LICENSE" file for the full license governing this code. 
 */
package didroe.tnc;

/**
 * The function to be minimized by the TNC algorithm.
 */
public interface TncFunction {

    /**
     * Evaluate the function and its gradient for the current parameter values.
     * @param x The current parameter values. Do not change this array.
     * @param gradient Array to store the result of evaluating the function's 
     * gradient for the current parameter values.
     * @return The function's value for the current parameters.
     */
    public double evaluate(double[] x, double[] gradient);
    
}
