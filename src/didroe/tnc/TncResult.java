/* 
 * See the "LICENSE" file for the full license governing this code. 
 */
package didroe.tnc;

/**
 *
 * @author Did
 */
public class TncResult {

    private final double[] parameters;
    private final double value;
    private final double[] gradient;
    private final int numFuncEvaluations;
    private final int numIterations;
    private final CompletionReason completionReason;
    
    TncResult(double[] parameters, double value, double[] gradient, int numFuncEvaluations, 
            int numIterations, CompletionReason completionReason) {
        this.parameters = parameters;
        this.value = value;
        this.gradient = gradient;
        this.numFuncEvaluations = numFuncEvaluations;
        this.numIterations = numIterations;
        this.completionReason = completionReason;
    }
    
    public double[] parameters() {
        return parameters;
    }
    
    public double value() {
        return value;
    }
    
    public double[] gradient() {
        return gradient;
    }
    
    public int numFuncEvaluations() {
        return numFuncEvaluations;
    }
    
    public int numIterations() {
        return numIterations;
    }
    
    public CompletionReason completionReason() {
        return completionReason;
    }
    
}
