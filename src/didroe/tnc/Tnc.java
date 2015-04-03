/* 
 * See the "LICENSE" file for the full license governing this code. 
 */
package didroe.tnc;

import java.util.Arrays;

/**
 * Builder for TNC inputs and invocation of the minimisation routine.
 *
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
     * @param lowerBounds Lower bounds for the function parameters. Set an entry to
     * Double.NEGATIVE_INFINITY for no bound on that parameter.
     */
    public Tnc lowerBounds(double[] lowerBounds) {
        this.lowerBounds = lowerBounds;
        return this;
    }

    /**
     * @param upperBounds Upper bounds for the function parameters. Set an entry to
     * Double.POSITIVE_INFINITY for no bound on that parameter.
     */
    public Tnc upperBounds(double[] upperBounds) {
        this.upperBounds = upperBounds;
        return this;
    }

    /**
     * @param scale Scaling factors to apply to each variable. If null, the factors are up-low for
     * interval bounded variables and 1+|x| for the others.
     */
    public Tnc scale(double[] scale) {
        this.scale = scale;
        return this;
    }

    /**
     * @param offset Value to subtract from each variable. If null, the offsets are (up+low)/2 for
     * interval bounded variables and x for the others.
     */
    public Tnc offset(double[] offset) {
        this.offset = offset;
        return this;
    }

    /**
     * @param maxConjugateGradientIterations The maximum number of hessian*vector evaluations per
     * main iteration. If maxCGit == 0, the direction chosen is -gradient. If maxCGit &lt; 0,
     * maxCGit is set to max(1, min(50,n/2)). Defaults to -1.
     */
    public Tnc maxConjugateGradientIterations(int maxConjugateGradientIterations) {
        this.maxConjugateGradientIterations = maxConjugateGradientIterations;
        return this;
    }

    /**
     * @param eta Severity of the line search. if &lt; 0 or &gt; 1, set to 0.25. Defaults to -1.
     */
    public Tnc eta(double eta) {
        this.eta = eta;
        return this;
    }

    /**
     * @param maxFunctionEvaluations Maximum number of function evaluations. Defaults to max(100,
     * 10*n) where n is the number of parameters.
     */
    public Tnc maxFunctionEvaluations(int maxFunctionEvaluations) {
        this.maxFunctionEvaluations = maxFunctionEvaluations;
        return this;
    }

    /**
     * @param maxLinearSearchStep Maximum step for the line search. May be increased during
     * minimization. If too small, it will be set to 10.0. Defaults to 0.
     */
    public Tnc maxLinearSearchStep(double maxLinearSearchStep) {
        this.maxLinearSearchStep = maxLinearSearchStep;
        return this;
    }

    /**
     * @param accuracy Relative precision for finite difference calculations. If &lt;=
     * machine_precision, set to sqrt(machine_precision). Defaults to 0.
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
     * @param valuePrecisionGoal Precision goal for the value of the function in the stopping
     * criterion. If &lt; 0.0, set to 0.0. Defaults to -1.
     */
    public Tnc valuePrecisionGoal(double valuePrecisionGoal) {
        this.valuePrecisionGoal = valuePrecisionGoal;
        return this;
    }

    /**
     * @param parameterPrecisionGoal Precision goal for the value of parameters in the stopping
     * criterion (after applying scaling factors). If &lt; 0.0, set to sqrt(machine_precision).
     * Defaults to -1.
     */
    public Tnc parameterPrecisionGoal(double parameterPrecisionGoal) {
        this.parameterPrecisionGoal = parameterPrecisionGoal;
        return this;
    }

    /**
     * @param gradientPrecisionGoal Precision goal for the value of the projected gradient in the
     * stopping criterion (after applying scaling factors). If &lt; 0.0, set to 1e-2 *
     * sqrt(accuracy). Setting it to 0.0 is not recommended. Defaults to -1.
     */
    public Tnc gradientPrecisionGoal(double gradientPrecisionGoal) {
        this.gradientPrecisionGoal = gradientPrecisionGoal;
        return this;
    }

    /**
     * @param rescale Scaling factor (in log10) used to trigger function value rescaling. If 0,
     * rescale at each iteration. If a large value, never rescale. If &lt; 0, rescale is set to 1.3.
     * Defaults to -1.
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

    private void validateInputs() throws InputError {
        if (function == null) {
            throw new InputError("Function argument must be specified");
        }

        if (initialGuess == null) {
            throw new InputError("Initial guess must be specified");
        }

        if (initialGuess.length == 0) {
            throw new InputError("There must be at least one parameter");
        }

        if (lowerBounds != null && lowerBounds.length != initialGuess.length) {
            throw new InputError("Lower bounds length does not match initial guess");
        }

        if (upperBounds != null && upperBounds.length != initialGuess.length) {
            throw new InputError("Upper bounds length does not match initial guess");
        }

        if (scale != null && scale.length != initialGuess.length) {
            throw new InputError("Scale length does not match initial guess");
        }

        if (offset != null && offset.length != initialGuess.length) {
            throw new InputError("Offset length does not match initial guess");
        }

    }

    public TncResult minimize() throws InputError {
        validateInputs();

        double[] effectiveLowerBounds = copyArray(lowerBounds);
        if (effectiveLowerBounds == null) {
            effectiveLowerBounds = new double[initialGuess.length];
            Arrays.fill(effectiveLowerBounds, Double.NEGATIVE_INFINITY);
        }

        double[] effectiveUpperBounds = copyArray(upperBounds);
        if (effectiveUpperBounds == null) {
            effectiveUpperBounds = new double[initialGuess.length];
            Arrays.fill(effectiveUpperBounds, Double.POSITIVE_INFINITY);
        }

        for (int i = 0; i < initialGuess.length; i++) {
            if (effectiveLowerBounds[i] > effectiveUpperBounds[i]) {
                throw new InputError("One or more lower bounds are greater than the upper bounds");
            }
        }

        int effectiveMaxFunctionEvaluations = (maxFunctionEvaluations != null
                ? maxFunctionEvaluations : Math.max(100, 10 * initialGuess.length));

        if (effectiveMaxFunctionEvaluations < 1) {
            throw new InputError("Max function evaluations must be >= 1");
        }

        int numConstantParams = 0;
        for (int i = 0; i < initialGuess.length; i++) {
            if ((effectiveLowerBounds[i] == effectiveUpperBounds[i])
                    || (scale != null && scale[i] == 0.0)) {
                numConstantParams++;
            }
        }

        if (numConstantParams == initialGuess.length) {
            throw new InputError("Lower bounds are equal to upper bounds or scale is zero for all parameters");
        }

        // Starting parameters are the initial guess constrained to the bounds
        double[] x = copyArray(initialGuess);
        ArrayMath.clip(x, effectiveLowerBounds, effectiveUpperBounds);

        TncRef ref = new TncRef();
        TncImpl t = new TncImpl();
        double[] g = new double[initialGuess.length];
        t.tnc(initialGuess.length, x, g, function, effectiveLowerBounds,
                effectiveUpperBounds, copyArray(scale), copyArray(offset),
                maxConjugateGradientIterations, effectiveMaxFunctionEvaluations, eta, 
                maxLinearSearchStep, accuracy, minFunctionValueEstimate, valuePrecisionGoal,
                parameterPrecisionGoal, gradientPrecisionGoal, rescale, callback, ref);
        return new TncResult(x, ref.f, g, t.functionEvaluator.getNumEvaluations(), t.numIterations);
    }

}
