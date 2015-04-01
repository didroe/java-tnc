/*
 * ORIGINAL SCIPY LICENSE HEADER:
 * Copyright (c) 2002-2005, Jean-Sebastien Roy (js@jeannot.org)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is an implementation of TNBC, a truncated newton minimization
 * package originally developed by Stephen G. Nash in Fortran.
 *
 * The original source code can be found at :
 * http://iris.gmu.edu/~snash/nash/software/software.html
 *
 * Copyright for the original TNBC fortran routines:
 *
 *   TRUNCATED-NEWTON METHOD:  SUBROUTINES
 *     WRITTEN BY:  STEPHEN G. NASH
 *           SCHOOL OF INFORMATION TECHNOLOGY & ENGINEERING
 *           GEORGE MASON UNIVERSITY
 *           FAIRFAX, VA 22030
 */
package didroe.tnc;

/**
 * Builder for TNC inputs
 * @author Did
 */
public class Tnc {
    
    private TncFunction function;
    private TncImpl.TncCallback callback;
    private double[] initialGuess;
    private double[] lowerBounds;
    private double[] upperBounds;
    private double[] scale;
    private double[] offset;
    private int maxConjugateGradientIterations = -1;
    private double eta = -1;
    private Integer maxFunctionEvaluations;
    private double maxLinearSearchStep = 0;
    private double accuracy = 0;
    private double minFunctionValueEstimate = 0;
    private double valuePrecisionGoal = -1;
    private double parameterPrecisionGoal = -1;
    private double gradientPrecisionGoal = -1;
    private double rescale = -1;
    
    /**
     * @param function The function to minimize
     */
    public Tnc function(TncFunction function) {
        this.function = function;
        return this;
    }
    
    /**
     * @param callback A callback to be invoked at each iteration
     */
    public Tnc callback(TncImpl.TncCallback callback) {
        this.callback = callback;
        return this;
    }

    /**
     * @param initialGuess Initial estimate of the minimum
     */
    public Tnc initialGuess(double[] initialGuess) {
        this.initialGuess = initialGuess;
        return this;
    }
    
    /**
     * @param lowerBounds Lower bounds for the function parameters. Set an entry
     * to Double.NEGATIVE_INFINITY for no bound on that parameter.
     */
    public Tnc lowerBounds(double[] lowerBounds) {
        this.lowerBounds = lowerBounds;
        return this;
    }
    
    /**
     * @param upperBounds Upper bounds for the function parameters. Set an entry
     * to Double.POSITIVE_INFINITY for no bound on that parameter.
     */
    public Tnc upperBounds(double[] upperBounds) {
        this.upperBounds = upperBounds;
        return this;
    }
    
    /**
     * @param scale Scaling factors to apply to each variable. If null, the
     * factors are up-low for interval bounded variables and 1+|x| for the others.
     */
    public Tnc scale(double[] scale) {
        this.scale = scale;
        return this;
    }
    
    /**
     * @param offset Value to subtract from each variable. If null, the offsets 
     * are (up+low)/2 for interval bounded variables and x for the others.
     */
    public Tnc offset(double[] offset) {
        this.offset = offset;
        return this;
    }

    /**
     * @param maxConjugateGradientIterations The maximum number of hessian*vector 
     * evaluations per main iteration. If maxCGit == 0, the direction chosen is
     * -gradient. If maxCGit &lt; 0, maxCGit is set to max(1, min(50,n/2)).  
     * Defaults to -1.
     */
    public Tnc maxConjugateGradientIterations(int maxConjugateGradientIterations) {
        this.maxConjugateGradientIterations = maxConjugateGradientIterations;
        return this;
    }
    
    /**
     * @param eta Severity of the line search. if &lt; 0 or &gt; 1, set to 0.25.
     * Defaults to -1.
     */
    public Tnc eta(double eta) {
       this.eta = eta;
       return this;
    }
    
    /**
     * @param maxFunctionEvaluations Maximum number of function evaluations. 
     * Defaults to max(100, 10*n) where n is the number of parameters.
     */
    public Tnc maxFunctionEvaluations(int maxFunctionEvaluations) {
        this.maxFunctionEvaluations = maxFunctionEvaluations;
        return this;
    }
    
    /**
     * @param maxLinearSearchStep Maximum step for the line search. May be 
     * increased during minimization. If too small, it will be set to 10.0. 
     * Defaults to 0.
     */
    public Tnc maxLinearSearchStep(double maxLinearSearchStep) {
        this.maxLinearSearchStep = maxLinearSearchStep;
        return this;
    }
    
    /**
     * @param accuracy Relative precision for finite difference calculations.
     * If &lt;= machine_precision, set to sqrt(machine_precision).
     * Defaults to 0.
     */
    public Tnc accuracy(double accuracy) {
        this.accuracy = accuracy;
        return this;
    }
    
    /**
     * @param minFunctionValueEstimate Minimum function value estimate. Defaults to 0.
     */
    public Tnc minFunctionValueEstimate(double minFunctionValueEstimate) {
        this.minFunctionValueEstimate = minFunctionValueEstimate;
        return this;
    }
    
    /**
     * @param valuePrecisionGoal Precision goal for the value of the function in
     * the stopping criterion. If &lt; 0.0, set to 0.0. Defaults to -1.
     */
    public Tnc valuePrecisionGoal(double valuePrecisionGoal) {
        this.valuePrecisionGoal = valuePrecisionGoal;
        return this;
    }
    
    /**
     * @param parameterPrecisionGoal Precision goal for the value of parameters 
     * in the stopping criterion (after applying scaling factors). If &lt; 0.0, 
     * set to sqrt(machine_precision). Defaults to -1.
     */
    public Tnc parameterPrecisionGoal(double parameterPrecisionGoal) {
        this.parameterPrecisionGoal = parameterPrecisionGoal;
        return this;
    }

    /**
     * @param gradientPrecisionGoal Precision goal for the value of the 
     * projected gradient in the stopping criterion (after applying scaling 
     * factors). If &lt; 0.0, set to 1e-2 * sqrt(accuracy). Setting it to 0.0 
     * is not recommended. Defaults to -1.
     */
    public Tnc gradientPrecisionGoal(double gradientPrecisionGoal) {
        this.gradientPrecisionGoal = gradientPrecisionGoal;
        return this;
    }
    
    /**
     * @param rescale Scaling factor (in log10) used to trigger function value
     * rescaling. If 0, rescale at each iteration. If a large value, never 
     * rescale. If &lt; 0, rescale is set to 1.3. Defaults to -1.
     */
    public Tnc rescale(double rescale) {
        this.rescale = rescale;
        return this;
    }
        
    private double[] copyArray(double[] src) {
        if (src == null) {
            return null;
        }
        
        double[] dest = new double[src.length];
        System.arraycopy(src, 0, dest, 0, src.length);
        return dest;
    }
    
    private void validate() {
        if (function == null) {
            throw new IllegalArgumentException("Function argument must be specified");
        }
        
        if (initialGuess == null) {
            throw new IllegalArgumentException("Initial guess must be specified");
        }
        
        // Check lower bound array length matches
        if (lowerBounds != null && lowerBounds.length != initialGuess.length) {
            throw new IllegalArgumentException("Lower bounds length does not match initial guess");
        }

        // Check upper bound array length matches
        if (upperBounds != null && upperBounds.length != initialGuess.length) {
            throw new IllegalArgumentException("Upper bounds length does not match initial guess");
        }

        // Check scale array length matches
        if (scale != null && scale.length != initialGuess.length) {
            throw new IllegalArgumentException("Scale length does not match initial guess");
        }

        // Check offset array length matches
        if (offset != null && offset.length != initialGuess.length) {
            throw new IllegalArgumentException("Offset length does not match initial guess");
        }
    }
    
    public TncResult minimize() {
        validate();
        
        double[] x = copyArray(initialGuess);
        
        int effectiveMaxFunctionEvaluations = (maxFunctionEvaluations != null ? 
                maxFunctionEvaluations : Math.max(100, 10 * initialGuess.length));
        
        TncRef ref = new TncRef();
        TncImpl t = new TncImpl();
        double[] g = new double[initialGuess.length];
        t.tnc(initialGuess.length, x, g, function, copyArray(lowerBounds), 
                copyArray(upperBounds), copyArray(scale), copyArray(offset), 
                TncImpl.tnc_message.allMessages(), maxConjugateGradientIterations, 
                effectiveMaxFunctionEvaluations, eta, maxLinearSearchStep, 
                accuracy, minFunctionValueEstimate, valuePrecisionGoal, 
                parameterPrecisionGoal, gradientPrecisionGoal, rescale, 
                callback, ref);
        return new TncResult(x, ref.f, g, ref.nfeval, ref.niter);
    }
    
}
