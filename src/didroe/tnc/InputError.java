/* 
 * See the "LICENSE" file for the full license governing this code. 
 */
package didroe.tnc;

/**
 * There was an error with one of the inputs to the minimisation algorithm
 * 
 * @author Did
 */
public class InputError extends RuntimeException {

    InputError(String message) {
        super(message);
    }
    
}
