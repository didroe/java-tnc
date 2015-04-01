/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
    
    TncResult(double[] parameters, double value, double[] gradient, int numFuncEvaluations, int numIterations) {
        this.parameters = parameters;
        this.value = value;
        this.gradient = gradient;
        this.numFuncEvaluations = numFuncEvaluations;
        this.numIterations = numIterations;
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
    
}
