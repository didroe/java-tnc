/*
 * See the "LICENSE" file for the full license governing this code.
 */
package didroe.tnc;

/**
 * There was an error during the minimisation.
 * 
 * @author Did
 */
public class MinimizationError extends RuntimeException {

    MinimizationError(String message) {
        super(message);
    }
    
}
