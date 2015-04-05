/* 
 * See the "LICENSE" file for the full license governing this code. 
 */
package didroe.tnc;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Did
 */
public class TncImpl {

    private final Logger logger = LoggerFactory.getLogger(TncImpl.class);

    class FunctionEvaluator {

        private int numEvaluations = 0;

        private final TncFunction function;

        private final int maxEvaluations;

        public FunctionEvaluator(TncFunction function, int maxEvaluations) {
            this.function = function;
            this.maxEvaluations = maxEvaluations;
        }

        public double evaluate(double[] x, double[] g) {
            if (numEvaluations >= maxEvaluations) {
                throw new MinimizationError("Maximum number of function evaluations reached");
            }

            numEvaluations++;
            return function.evaluate(x, g);
        }

        public int getNumEvaluations() {
            return numEvaluations;
        }
    }

    /**
     * A callback function accepting x as input parameter along with the state pointer.
     */
    public interface TncCallback {

        public void tnc_callback(double[] x);
    }

    public static final double DBL_EPSILON = Math.ulp(1.0);

    FunctionEvaluator functionEvaluator;

    /*
     * tnc : minimize a function with variables subject to bounds, using
     *       gradient information.
     *
     * n         : number of variables (must be >= 0)
     * x         : on input, initial estimate ; on output, the solution
     * f         : on output, the function value at the solution
     * g         : on output, the gradient value at the solution
     *             g should be an allocated vector of size n or NULL,
     *             in which case the gradient value is not returned.
     * function  : the function to minimize (see tnc_function)
     * state     : used by function (see tnc_function)
     * low, up   : the bounds
     *             set low[i] to Double.NEGATIVE_INFINITY to remove the lower bound
     *             set up[i] to Double.POSITIVE_INFINITY to remove the upper bound
     *             if low == NULL, the lower bounds are removed.
     *             if up == NULL, the upper bounds are removed.
     * scale     : scaling factors to apply to each variable
     *             if NULL, the factors are up-low for interval bounded variables
     *             and 1+|x] for the others.
     * offset    : constant to substract to each variable
     *             if NULL, the constant are (up+low)/2 for interval bounded
     *             variables and x for the others.
     * messages  : see the tnc_message enum
     * maxCGit   : max. number of hessian*vector evaluation per main iteration
     *             if maxCGit == 0, the direction chosen is -gradient
     *             if maxCGit < 0, maxCGit is set to max(1,min(50,n/2))
     * maxnfeval : max. number of function evaluation
     * eta       : severity of the line search. if < 0 or > 1, set to 0.25
     * stepmx    : maximum step for the line search. may be increased during call
     *             if too small, will be set to 10.0
     * accuracy  : relative precision for finite difference calculations
     *             if <= machine_precision, set to sqrt(machine_precision)
     * fmin      : minimum function value estimate
     * ftol      : precision goal for the value of f in the stoping criterion
     *             if ftol < 0.0, ftol is set to accuracy
     * xtol      : precision goal for the value of x in the stopping criterion
     *             (after applying x scaling factors)
     *             if xtol < 0.0, xtol is set to sqrt(machine_precision)
     * pgtol     : precision goal for the value of the projected gradient in the
     *             stopping criterion (after applying x scaling factors)
     *             if pgtol < 0.0, pgtol is set to 1e-2 * sqrt(accuracy)
     *             setting it to 0.0 is not recommended
     * rescale   : f scaling factor (in log10) used to trigger f value rescaling
     *             if 0, rescale at each iteration
     *             if a big value, never rescale
     *             if < 0, rescale is set to 1.3
     * nfeval    : on output, the number of function evaluations.
     *             ignored if nfeval==NULL.
     *
     * The tnc function returns a code defined in the tnc_rc enum.
     * On output, x, f and g may be very slightly out of sync because of scaling.
     *
     * This routine solves the optimization problem
     *
     *   minimize   f(x)
     *     x
     *   subject to   low <= x <= up
     *
     * where x is a vector of n real variables. The method used is
     * a truncated-newton algorithm (see "newton-type minimization via
     * the lanczos method" by s.g. nash (siam j. numer. anal. 21 (1984),
     * pp. 770-778).  this algorithm finds a local minimum of f(x). It does
     * not assume that the function f is convex (and so cannot guarantee a
     * global solution), but does assume that the function is bounded below.
     * it can solve problems having any number of variables, but it is
     * especially useful when the number of variables (n) is large.    
     */
    int numIterations = 0;

    public CompletionReason tnc(int n, double[] x, double[] g, TncFunction function,
            double[] low, double[] up, double[] scale,
            double[] offset, int maxCGit, int maxnfeval,
            double eta, double stepmx, double accuracy, double fmin,
            double ftol, double xtol, double pgtol, double rescale,
            TncCallback callback, TncRef ref) {
        /* Scaling parameters */
        double[] xscale = new double[n];
        double[] xoffset = new double[n];

        for (int i = 0; i < n; i++) {
            if (scale != null) {
                xscale[i] = Math.abs(scale[i]);
                if (xscale[i] == 0.0) {
                    xoffset[i] = low[i] = up[i] = x[i];
                }
            } else if (low[i] != Double.NEGATIVE_INFINITY && up[i] != Double.POSITIVE_INFINITY) {
                xscale[i] = up[i] - low[i];
                xoffset[i] = (up[i] + low[i]) * 0.5;
            } else {
                xscale[i] = 1.0 + Math.abs(x[i]);
                xoffset[i] = x[i];
            }
            if (offset != null) {
                xoffset[i] = offset[i];
            }
        }

        /* Default values for parameters */
        double rteps = Math.sqrt(DBL_EPSILON);

        if (stepmx < rteps * 10.0) {
            stepmx = 1.0e1;
        }
        if (eta < 0.0 || eta >= 1.0) {
            eta = 0.25;
        }
        if (rescale < 0) {
            rescale = 1.3;
        }
        if (maxCGit < 0) {          /* maxCGit == 0 is valid */

            maxCGit = n / 2;
            if (maxCGit < 1) {
                maxCGit = 1;
            } else if (maxCGit > 50) {
                maxCGit = 50;
            }
        }
        if (maxCGit > n) {
            maxCGit = n;
        }
        if (accuracy <= DBL_EPSILON) {
            accuracy = rteps;
        }
        if (ftol < 0.0) {
            ftol = accuracy;
        }
        if (pgtol < 0.0) {
            pgtol = 1e-2 * Math.sqrt(accuracy);
        }
        if (xtol < 0.0) {
            xtol = rteps;
        }

        /* Optimisation */
        tnc_minimizeRef tmpMinRef = new tnc_minimizeRef();
        CompletionReason rc = tnc_minimize(n, x, g, function, xscale, xoffset, low, up, maxCGit, 
                maxnfeval, eta, stepmx, accuracy, fmin, ftol, xtol, pgtol, rescale, callback, tmpMinRef);
        ref.f = tmpMinRef.f;

        logger.info(rc.getMessage());

        return rc;
    }

    /**
     * Unscale x
     */
    public void unscalex(double[] x, double[] xscale, double xoffset[]) {
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i] * xscale[i] + xoffset[i];
        }
    }

    /**
     * Scale x
     */
    public void scalex(double[] x, double[] xscale, double[] xoffset) {
        for (int i = 0; i < x.length; i++) {
            if (xscale[i] > 0.0) {
                x[i] = (x[i] - xoffset[i]) / xscale[i];
            }
        }
    }

    /**
     * Scale g
     */
    public void scaleg(double[] g, double[] xscale, double fscale) {
        for (int i = 0; i < g.length; i++) {
            g[i] *= xscale[i] * fscale;
        }
    }

    /**
     * Calculate the pivot vector
     */
    public void setConstraints(int n, double[] x, int[] pivot, double[] xscale,
            double[] xoffset, double[] low, double[] up) {
        for (int i = 0; i < n; i++) {
            /* tolerances should be better ajusted */
            if (xscale[i] == 0.0) {
                pivot[i] = 2;
            } else if (low[i] != Double.NEGATIVE_INFINITY
                    && (x[i] * xscale[i] + xoffset[i] - low[i]
                    <= DBL_EPSILON * 10.0 * (Math.abs(low[i]) + 1.0))) {
                pivot[i] = -1;
            } else if (up[i] != Double.POSITIVE_INFINITY
                    && (x[i] * xscale[i] + xoffset[i] - up[i]
                    >= DBL_EPSILON * 10.0 * (Math.abs(up[i]) + 1.0))) {
                pivot[i] = 1;
            } else {
                pivot[i] = 0;
            }
        }
    }

    public class tnc_minimizeRef {

        public double f;
    }

    /**
     * This routine is a bounds-constrained truncated-newton method. the truncated-newton method is
     * preconditioned by a limited-memory quasi-newton method (this preconditioning strategy is
     * developed in this routine) with a further diagonal scaling (see routine diagonalscaling).
     */
    public CompletionReason tnc_minimize(int n, double[] x, double[] gfull,
            TncFunction function, double[] xscale, double[] xoffset,
            double[] low, double[] up,
            int maxCGit, int maxnfeval, double eta, double stepmx,
            double accuracy, double fmin, double ftol, double xtol,
            double pgtol, double rescale, TncCallback callback, tnc_minimizeRef ref) {
        functionEvaluator = new FunctionEvaluator(function, maxnfeval);

        double fscale = 1.0;
        double alpha = 0.0;         /* Default unused value */

        /* Allocate temporary vectors */
        double[] oldg = new double[n];
        double[] g = new double[n];
        double[] temp = new double[n];
        double[] diagb = new double[n];
        double[] searchDirection = new double[n];
        double[] sk = new double[n];
        double[] yk = new double[n];
        double[] sr = new double[n];
        double[] yr = new double[n];
        int[] pivot = new int[n];

        /* Initialize variables */
        double difnew = 0.0;
        double epsred = 0.05;
        boolean upd1 = true;
        int icycle = n - 1;
        boolean newcon = true;

        /* Uneeded initialisations */
        boolean lreset = false;
        double yrsr = 0.0;
        double yksk = 0.0;

        /* Initial function evaluation */
        ref.f = functionEvaluator.evaluate(x, gfull);

        /* Initial scaling */
        scalex(x, xscale, xoffset);
        ref.f *= fscale;

        /* initial pivot calculation */
        setConstraints(n, x, pivot, xscale, xoffset, low, up);

        ArrayMath.copy(gfull, g);
        scaleg(g, xscale, fscale);

        /* Test the lagrange multipliers to see if they are non-negative. */
        for (int i = 0; i < n; i++) {
            if (-pivot[i] * g[i] < 0.0) {
                pivot[i] = 0;
            }
        }

        project(g, pivot);

        /* Set initial values to other parameters */
        double gnorm = ArrayMath.euclidianNorm(g);

        double fLastConstraint = ref.f;       /* Value at last constraint */

        double fLastReset = ref.f;            /* Value at last reset */

        if (logger.isTraceEnabled()) {
            logger.trace("  NIT   NF   F                       GTG");
            printCurrentIteration(n, ref.f / fscale, gfull,
                    numIterations, functionEvaluator.getNumEvaluations(), pivot);
        }

        /* Set the diagonal of the approximate hessian to unity. */
        for (int i = 0; i < n; i++) {
            diagb[i] = 1.0;
        }

        CompletionReason rc;
        /* Start of main iterative loop */
        while (true) {
            /* Local minimum test */
            if (ArrayMath.euclidianNorm(g) <= pgtol * fscale) {
                /* |PG| == 0.0 => local minimum */
                ArrayMath.copy(gfull, g);
                project(g, pivot);
                if (logger.isDebugEnabled()) {
                    logger.debug("|pg| = {} -> local minimum", ArrayMath.euclidianNorm(g) / fscale);
                }
                rc = CompletionReason.TNC_LOCALMINIMUM;
                break;
            }

            /* Rescale function if necessary */
            double newscale = ArrayMath.euclidianNorm(g);
            if ((newscale > DBL_EPSILON) && (Math.abs(Math.log10(newscale)) > rescale)) {
                newscale = 1.0 / newscale;

                ref.f *= newscale;
                fscale *= newscale;
                gnorm *= newscale;
                fLastConstraint *= newscale;
                fLastReset *= newscale;
                difnew *= newscale;

                for (int i = 0; i < n; i++) {
                    g[i] *= newscale;
                }
                for (int i = 0; i < n; i++) {
                    diagb[i] = 1.0;
                }

                upd1 = true;
                icycle = n - 1;
                newcon = true;

                logger.debug("fscale = {}", fscale);
            }

            ArrayMath.copy(x, temp);
            project(temp, pivot);
            double xnorm = ArrayMath.euclidianNorm(temp);
            int oldnfeval = functionEvaluator.getNumEvaluations();

            /* Compute the new search direction */
            tnc_directionRef tmpDirRef = new tnc_directionRef();
            tmpDirRef.diagb = diagb;
            tmpDirRef.pivot = pivot;
            tmpDirRef.sk = sk;
            tmpDirRef.sr = sr;
            tmpDirRef.x = x;
            tmpDirRef.yk = yk;
            tmpDirRef.yr = yr;
            tmpDirRef.zsol = searchDirection;
            tnc_direction(g, n, maxCGit, maxnfeval, upd1, yksk, yrsr, lreset, xscale,
                    xoffset, fscale, accuracy, gnorm, xnorm, low, up, tmpDirRef);
            diagb = tmpDirRef.diagb;
            pivot = tmpDirRef.pivot;
            sk = tmpDirRef.sk;
            sr = tmpDirRef.sr;
            x = tmpDirRef.x;
            yk = tmpDirRef.yk;
            yr = tmpDirRef.yr;
            searchDirection = tmpDirRef.zsol;

            if (!newcon) {
                if (!lreset) {
                    /* Compute the accumulated step and its corresponding gradient
                     difference. */
                    ArrayMath.xPlusY(sk, sr);
                    ArrayMath.xPlusY(yk, yr);
                    icycle++;
                } else {
                    /* Initialize the sum of all the changes */
                    ArrayMath.copy(sk, sr);
                    ArrayMath.copy(yk, yr);
                    fLastReset = ref.f;
                    icycle = 1;
                }
            }

            ArrayMath.copy(g, oldg);
            double oldf = ref.f;
            double oldgtp = ArrayMath.dotProduct(searchDirection, g);

            /* Maximum unconstrained step length */
            double ustpmax = stepmx / (ArrayMath.euclidianNorm(searchDirection) + DBL_EPSILON);

            /* Maximum constrained step length */
            double spe = stepMax(ustpmax, n, x, searchDirection, pivot, low, up, xscale, xoffset);

            if (spe > 0.0) {
                /* Set the initial step length */
                alpha = initialStep(ref.f, fmin / fscale, oldgtp, spe);

                /* Perform the linear search */
                LineSearchResult lsResult = lineSearch(low, up, xscale, xoffset, fscale, pivot,
                        eta, ftol, spe, searchDirection, x, ref.f, alpha, gfull);
                alpha = lsResult.alpha();
                ref.f = lsResult.f();
                ArrayMath.copy(lsResult.g(), gfull);
                ArrayMath.copy(lsResult.x(), x);


                /* If we went up to the maximum unconstrained step, increase it */
                if (alpha >= 0.9 * ustpmax) {
                    stepmx *= 1e2;
                    logger.debug("stepmx = {}", stepmx);
                }

                /* If we went up to the maximum constrained step,
                 a new constraint was encountered */
                if (alpha - spe >= -DBL_EPSILON * 10.0) {
                    newcon = true;
                } else {
                    newcon = false;
                }
            } else {
                /* Maximum constrained step == 0.0 => new constraint */
                newcon = true;
            }

            if (newcon) {
                if (!addConstraint(n, x, searchDirection, pivot, low, up, xscale, xoffset)) {
                    if (functionEvaluator.getNumEvaluations() == oldnfeval) {
                        throw new MinimizationError("Unable to progress");
                    }
                }
                fLastConstraint = ref.f;
            }

            numIterations++;

            /* Invoke the callback function */
            if (callback != null) {
                unscalex(x, xscale, xoffset);
                callback.tnc_callback(x);
                scalex(x, xscale, xoffset);
            }

            /* Set up parameters used in convergence and resetting tests */
            double difold = difnew;
            difnew = oldf - ref.f;

            /* If this is the first iteration of a new cycle, compute the
             percentage reduction factor for the resetting test */
            if (icycle == 1) {
                if (difnew > difold * 2.0) {
                    epsred += epsred;
                }
                if (difnew < difold * 0.5) {
                    epsred *= 0.5;
                }
            }

            ArrayMath.copy(gfull, g);
            scaleg(g, xscale, fscale);

            ArrayMath.copy(g, temp);
            project(temp, pivot);
            gnorm = ArrayMath.euclidianNorm(temp);

            /* Reset pivot */
            boolean remcon = removeConstraint(oldgtp, gnorm, pgtol * fscale, ref.f,
                    fLastConstraint, g, pivot, n);

            /* If a constraint is removed */
            if (remcon) {
                /* Recalculate gnorm and reset fLastConstraint */
                ArrayMath.copy(g, temp);
                project(temp, pivot);
                gnorm = ArrayMath.euclidianNorm(temp);
                fLastConstraint = ref.f;
            }

            if (!remcon && !newcon) {
                /* No constraint removed & no new constraint : tests for convergence */
                if (Math.abs(difnew) <= ftol * fscale) {
                    if (logger.isDebugEnabled()) {
                        logger.debug("|fn-fn-1] = {} -> convergence", Math.abs(difnew) / fscale);
                    }
                    rc = CompletionReason.TNC_FCONVERGED;
                    break;
                }
                if (alpha * ArrayMath.euclidianNorm(searchDirection) <= xtol) {
                    if (logger.isDebugEnabled()) {
                        logger.debug("|xn-xn-1] = %g -> convergence",
                                alpha * ArrayMath.euclidianNorm(searchDirection));
                    }
                    rc = CompletionReason.TNC_XCONVERGED;
                    break;
                }
            }

            project(g, pivot);

            if (logger.isTraceEnabled()) {
                printCurrentIteration(n, ref.f / fscale, gfull,
                        numIterations, functionEvaluator.getNumEvaluations(), pivot);
            }

            /* Compute the change in the iterates and the corresponding change in the
             gradients */
            if (!newcon) {
                for (int i = 0; i < n; i++) {
                    yk[i] = g[i] - oldg[i];
                    sk[i] = alpha * searchDirection[i];
                }

                /* Set up parameters used in updating the preconditioning strategy */
                yksk = ArrayMath.dotProduct(yk, sk);

                if (icycle == (n - 1) || difnew < epsred * (fLastReset - ref.f)) {
                    lreset = true;
                } else {
                    yrsr = ArrayMath.dotProduct(yr, sr);
                    if (yrsr <= 0.0) {
                        lreset = true;
                    } else {
                        lreset = false;
                    }
                }
                upd1 = false;
            }
        }

        if (logger.isTraceEnabled()) {
            printCurrentIteration(n, ref.f / fscale, gfull,
                    numIterations, functionEvaluator.getNumEvaluations(), pivot);
        }

        /* Unscaling */
        unscalex(x, xscale, xoffset);
        ArrayMath.clip(x, low, up);
        ref.f /= fscale;

        return rc;
    }

    /**
     * Print the results of the current iteration
     */
    public void printCurrentIteration(int n, double f, double[] g, int niter,
            int nfeval, int[] pivot) {
        double gtg = 0.0;
        for (int i = 0; i < n; i++) {
            if (pivot[i] == 0) {
                gtg += g[i] * g[i];
            }
        }

        logger.trace(String.format(" %4d %4d %22.15E  %15.8E", niter, nfeval, f, gtg));
    }

    /**
     * Set x[i] = 0.0 if direction i is currently constrained
     */
    public void project(double[] x, int[] pivot) {
        for (int i = 0; i < x.length; i++) {
            if (pivot[i] != 0) {
                x[i] = 0.0;
            }
        }
    }

    /**
     * Set x[i] = 0.0 if direction i is constant
     */
    public void projectConstants(int n, double[] x, double[] xscale) {
        for (int i = 0; i < n; i++) {
            if (xscale[i] == 0.0) {
                x[i] = 0.0;
            }
        }
    }

    /**
     * Compute the maximum allowable step length
     */
    public double stepMax(double step, int n, double[] x, double[] dir,
            int[] pivot, double[] low, double[] up, double[] xscale, double[] xoffset) {
        /* Constrained maximum step */
        for (int i = 0; i < n; i++) {
            if ((pivot[i] == 0) && (dir[i] != 0.0)) {
                if (dir[i] < 0.0) {
                    double t = (low[i] - xoffset[i]) / xscale[i] - x[i];
                    if (t > step * dir[i]) {
                        step = t / dir[i];
                    }
                } else {
                    double t = (up[i] - xoffset[i]) / xscale[i] - x[i];
                    if (t < step * dir[i]) {
                        step = t / dir[i];
                    }
                }
            }
        }

        return step;
    }

    /**
     * Update the constraint vector pivot if a new constraint is encountered
     */
    public boolean addConstraint(int n, double[] x, double[] p, int[] pivot,
            double[] low, double[] up, double[] xscale, double[] xoffset) {
        boolean newcon = false;
        for (int i = 0; i < n; i++) {
            if ((pivot[i] == 0) && (p[i] != 0.0)) {
                if (p[i] < 0.0 && low[i] != Double.NEGATIVE_INFINITY) {
                    double tol = DBL_EPSILON * 10.0 * (Math.abs(low[i]) + 1.0);
                    if (x[i] * xscale[i] + xoffset[i] - low[i] <= tol) {
                        pivot[i] = -1;
                        x[i] = (low[i] - xoffset[i]) / xscale[i];
                        newcon = true;
                    }
                } else if (up[i] != Double.POSITIVE_INFINITY) {
                    double tol = DBL_EPSILON * 10.0 * (Math.abs(up[i]) + 1.0);
                    if (up[i] - (x[i] * xscale[i] + xoffset[i]) <= tol) {
                        pivot[i] = 1;
                        x[i] = (up[i] - xoffset[i]) / xscale[i];
                        newcon = true;
                    }
                }
            }
        }
        return newcon;
    }

    /**
     * Check if a constraint is no more active
     */
    public boolean removeConstraint(double gtpnew, double gnorm, double pgtolfs,
            double f, double fLastConstraint, double[] g, int[] pivot, int n) {
        if (((fLastConstraint - f) <= (gtpnew * -0.5)) && (gnorm > pgtolfs)) {
            return false;
        }

        int imax = -1;
        double cmax = 0.0;

        for (int i = 0; i < n; i++) {
            if (pivot[i] == 2) {
                continue;
            }
            double t = -pivot[i] * g[i];
            if (t < cmax) {
                cmax = t;
                imax = i;
            }
        }

        if (imax != -1) {
            pivot[imax] = 0;
            return true;
        } else {
            return false;
        }

        /*
         * For details, see gill, murray, and wright (1981, p. 308) and
         * fletcher (1981, p. 116). The multiplier tests (here, testing
         * the sign of the components of the gradient) may still need to
         * modified to incorporate tolerances for zero.
         */
    }

    public class tnc_directionRef {

        public double[] zsol;
        public double[] diagb;
        public double[] x;
        public double[] sk;
        public double[] yk;
        public double[] sr;
        public double[] yr;
        public int[] pivot;
    }

    /**
     * This routine performs a preconditioned conjugate-gradient iteration in order to solve the
     * newton equations for a search direction for a truncated-newton algorithm. When the value of
     * the quadratic model is sufficiently reduced, the iteration is terminated.
     */
    public void tnc_direction(double[] g, int n, int maxCGit, int maxnfeval,
            boolean upd1, double yksk, double yrsr, boolean lreset,
            double[] xscale, double[] xoffset,
            double fscale, double accuracy, double gnorm, double xnorm,
            double[] low, double[] up, tnc_directionRef ref) {
        /* No CG it. => dir = -grad */
        if (maxCGit == 0) {
            ArrayMath.copy(g, ref.zsol);
            ArrayMath.negate(ref.zsol);
            project(ref.zsol, ref.pivot);
            return;
        }

        /* General initialization */
        double rhsnrm = gnorm;
        double tol = 1e-12;
        double qold = 0.0;
        double rzold = 0.0;                /* Uneeded */

        double[] r = new double[n]; /* Residual */

        double[] v = new double[n];
        double[] zk = new double[n];
        double[] emat = new double[n];   /* Diagonal preconditoning matrix */

        double[] gv = new double[n];       /* hessian times v */

        /* Initialization for preconditioned conjugate-gradient algorithm */
        initPreconditioner(ref.diagb, emat, n, lreset, yksk, yrsr, ref.sk, ref.yk, ref.sr,
                ref.yr, upd1);

        for (int i = 0; i < n; i++) {
            r[i] = -g[i];
            v[i] = 0.0;
            ref.zsol[i] = 0.0;          /* Computed search direction */

        }

        /* Main iteration */
        for (int k = 0; k < maxCGit; k++) {
            /* CG iteration to solve system of equations */
            project(r, ref.pivot);
            msolve(r, zk, n, ref.sk, ref.yk, ref.diagb, ref.sr, ref.yr, upd1, yksk, yrsr,
                    lreset);
            project(zk, ref.pivot);
            double rz = ArrayMath.dotProduct(r, zk);

            if (rz / rhsnrm < tol) {
                /* Truncate algorithm in case of an emergency */
                if (k == 0) {
                    ArrayMath.copy(g, ref.zsol);
                    ArrayMath.negate(ref.zsol);
                    project(ref.zsol, ref.pivot);
                }
                break;
            }
            double beta;
            if (k == 0) {
                beta = 0.0;
            } else {
                beta = rz / rzold;
            }

            for (int i = 0; i < n; i++) {
                v[i] = zk[i] + beta * v[i];
            }

            project(v, ref.pivot);
            hessianTimesVector(v, gv, n, ref.x, g, xscale, xoffset, fscale,
                    accuracy, xnorm, low, up);
            project(gv, ref.pivot);

            double vgv = ArrayMath.dotProduct(v, gv);
            if (vgv / rhsnrm < tol) {
                /* Truncate algorithm in case of an emergency */
                if (k == 0) {
                    msolve(g, ref.zsol, n, ref.sk, ref.yk, ref.diagb, ref.sr, ref.yr, upd1, yksk,
                            yrsr, lreset);
                    ArrayMath.negate(ref.zsol);
                    project(ref.zsol, ref.pivot);
                }
                break;
            }
            diagonalScaling(n, emat, v, gv, r);

            /* Compute linear step length */
            double alpha = rz / vgv;

            /* Compute current solution and related vectors */
            ArrayMath.axPlusY(alpha, v, ref.zsol);
            ArrayMath.axPlusY(-alpha, gv, r);

            /* Test for convergence */
            double gtp = ArrayMath.dotProduct(ref.zsol, g);
            double pr = ArrayMath.dotProduct(r, ref.zsol);
            double qnew = (gtp + pr) * 0.5;
            double qtest = (k + 1) * (1.0 - qold / qnew);
            if (qtest <= 0.5) {
                break;
            }

            /* Perform cautionary test */
            if (gtp > 0.0) {
                /* Truncate algorithm in case of an emergency */
                ArrayMath.axPlusY(-alpha, v, ref.zsol);
                break;
            }

            qold = qnew;
            rzold = rz;
        }

        /* Terminate algorithm */
        /* Store (or restore) diagonal preconditioning */
        ArrayMath.copy(emat, ref.diagb);
    }

    /**
     * Update the preconditioning matrix based on a diagonal version of the bfgs quasi-newton
     * update.
     */
    public void diagonalScaling(int n, double[] e, double[] v, double[] gv,
            double[] r) {
        double vr = 1.0 / ArrayMath.dotProduct(v, r);
        double vgv = 1.0 / ArrayMath.dotProduct(v, gv);
        for (int i = 0; i < n; i++) {
            e[i] += -r[i] * r[i] * vr + gv[i] * gv[i] * vgv;
            if (e[i] <= 1e-6) {
                e[i] = 1.0;
            }
        }
    }

    /**
     * Returns the length of the initial step to be taken along the vector p in the next linear
     * search.
     */
    public double initialStep(double fnew, double fmin, double gtp, double smax) {
        double d = Math.abs(fnew - fmin);
        double alpha = 1.0;
        if (d * 2.0 <= -(gtp) && d >= DBL_EPSILON) {
            alpha = d * -2.0 / gtp;
        }
        if (alpha >= smax) {
            alpha = smax;
        }

        return alpha;
    }

    /**
     * Hessian vector product through finite differences
     */
    public void hessianTimesVector(double[] v, double[] gv, int n, double[] x,
            double[] g, double[] xscale, double[] xoffset,
            double fscale, double accuracy, double xnorm, double[] low,
            double[] up) {
        double[] xv = new double[n];

        double delta = accuracy * (xnorm + 1.0);
        for (int i = 0; i < n; i++) {
            xv[i] = x[i] + delta * v[i];
        }

        unscalex(xv, xscale, xoffset);
        ArrayMath.clip(xv, low, up);
        functionEvaluator.evaluate(xv, gv);

        scaleg(gv, xscale, fscale);

        double dinv = 1.0 / delta;
        for (int i = 0; i < n; i++) {
            gv[i] = (gv[i] - g[i]) * dinv;
        }

        projectConstants(n, gv, xscale);
    }

    /**
     * This routine acts as a preconditioning step for the linear conjugate-gradient routine. It is
     * also the method of computing the search direction from the gradient for the non-linear
     * conjugate-gradient code. It represents a two-step self-scaled bfgs formula.
     */
    public void msolve(double[] g, double[] y, int n,
            double[] sk, double[] yk, double[] diagb, double[] sr,
            double[] yr, boolean upd1, double yksk, double yrsr,
            boolean lreset) {
        if (upd1) {
            for (int i = 0; i < n; i++) {
                y[i] = g[i] / diagb[i];
            }
            return;
        }

        double gsk = ArrayMath.dotProduct(g, sk);
        double[] hg = new double[n];
        double[] hyr = new double[n];
        double[] hyk = new double[n];

        /* Compute gh and hy where h is the inverse of the diagonals */
        if (lreset) {
            for (int i = 0; i < n; i++) {
                double rdiagb = 1.0 / diagb[i];
                hg[i] = g[i] * rdiagb;
                hyk[i] = yk[i] * rdiagb;
            }
            double ykhyk = ArrayMath.dotProduct(yk, hyk);
            double ghyk = ArrayMath.dotProduct(g, hyk);
            ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
        } else {
            for (int i = 0; i < n; i++) {
                double rdiagb = 1.0 / diagb[i];
                hg[i] = g[i] * rdiagb;
                hyk[i] = yk[i] * rdiagb;
                hyr[i] = yr[i] * rdiagb;
            }
            double gsr = ArrayMath.dotProduct(g, sr);
            double ghyr = ArrayMath.dotProduct(g, hyr);
            double yrhyr = ArrayMath.dotProduct(yr, hyr);
            ssbfgs(n, 1.0, sr, hg, hyr, yrsr, yrhyr, gsr, ghyr, hg);
            double yksr = ArrayMath.dotProduct(yk, sr);
            double ykhyr = ArrayMath.dotProduct(yk, hyr);
            ssbfgs(n, 1.0, sr, hyk, hyr, yrsr, yrhyr, yksr, ykhyr, hyk);
            double ykhyk = ArrayMath.dotProduct(hyk, yk);
            double ghyk = ArrayMath.dotProduct(hyk, g);
            ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
        }
    }

    /**
     * Self-scaled BFGS
     */
    public void ssbfgs(int n, double gamma, double[] sj, double[] hjv,
            double[] hjyj, double yjsj, double yjhyj, double vsj, double vhyj,
            double[] hjp1v) {
        double beta, delta;
        if (yjsj == 0.0) {
            delta = 0.0;
            beta = 0.0;
        } else {
            delta = (gamma * yjhyj / yjsj + 1.0) * vsj / yjsj
                    - gamma * vhyj / yjsj;
            beta = -gamma * vsj / yjsj;
        }

        for (int i = 0; i < n; i++) {
            hjp1v[i] = gamma * hjv[i] + delta * sj[i] + beta * hjyj[i];
        }
    }

    /**
     * Initialize the preconditioner
     */
    public void initPreconditioner(double[] diagb, double[] emat, int n,
            boolean lreset, double yksk, double yrsr, double[] sk, double[] yk,
            double[] sr, double[] yr, boolean upd1) {
        if (upd1) {
            ArrayMath.copy(diagb, emat);
            return;
        }

        double[] bsk = new double[n];

        if (lreset) {
            for (int i = 0; i < n; i++) {
                bsk[i] = diagb[i] * sk[i];
            }
            double sds = ArrayMath.dotProduct(sk, bsk);
            if (yksk == 0.0) {
                yksk = 1.0;
            }
            if (sds == 0.0) {
                sds = 1.0;
            }
            for (int i = 0; i < n; i++) {
                double td = diagb[i];
                emat[i] = td - td * td * sk[i] * sk[i] / sds
                        + yk[i] * yk[i] / yksk;
            }
        } else {
            for (int i = 0; i < n; i++) {
                bsk[i] = diagb[i] * sr[i];
            }
            double sds = ArrayMath.dotProduct(sr, bsk);
            double srds = ArrayMath.dotProduct(sk, bsk);
            double yrsk = ArrayMath.dotProduct(yr, sk);
            if (yrsr == 0.0) {
                yrsr = 1.0;
            }
            if (sds == 0.0) {
                sds = 1.0;
            }
            for (int i = 0; i < n; i++) {
                double td = diagb[i];
                bsk[i] = td * sk[i] - bsk[i] * srds / sds + yr[i] * yrsk / yrsr;
                emat[i] = td - td * td * sr[i] * sr[i] / sds + yr[i] * yr[i] / yrsr;
            }
            sds = ArrayMath.dotProduct(sk, bsk);
            if (yksk == 0.0) {
                yksk = 1.0;
            }
            if (sds == 0.0) {
                sds = 1.0;
            }
            for (int i = 0; i < n; i++) {
                emat[i] -= bsk[i] * bsk[i] / sds + yk[i] * yk[i] / yksk;
            }
        }
    }

    class LineSearchResult {
        private final double alpha;
        private final double[] x;
        private final double f;
        private final double[] g;
        
        public LineSearchResult(double alpha, double[] x, double f, double[] g) {
            this.alpha = alpha;
            this.x = x;
            this.f = f;
            this.g = g;
        }
        
        public double alpha() {
            return alpha;
        }
        
        public double[] x() {
            return x;
        }
        
        public double f() {
            return f;
        }

        public double[] g() {
            return g;
        }
    }

    
    /**
     * Line search algorithm of Gill and Murray.
     * This solves the problem:
     *      minimise g(alpha), subject to alpha > 0
     *      where g(alpha) = f(x + alpha.p) where p is the current search 
     * direction of the outer minimisation problem and x is the current 
     * solution of said problem.
     * The search direction must be a descent direction of sufficient descent
     * and must be gradient related.
     */
    public LineSearchResult lineSearch(double[] low,
            double[] up, double[] xscale, double[] xoffset, double fscale,
            int[] pivot, double eta, double ftol, double xBound, double[] searchDirection,
            double[] x, double f, double alpha, double gfull[]) {
        final int MAX_ITERATIONS = 64;

        double[] temp = new double[x.length];
        double[] tempgfull = new double[x.length];
        double[] newgfull = new double[x.length];

        ArrayMath.copy(gfull, temp);
        scaleg(temp, xscale, fscale);
        double initGu = ArrayMath.dotProduct(temp, searchDirection);

        ArrayMath.copy(x, temp);
        project(temp, pivot);
        double xnorm = ArrayMath.euclidianNorm(temp);

        double rteps = Math.sqrt(DBL_EPSILON);
        double pe = ArrayMath.euclidianNorm(searchDirection) + DBL_EPSILON;

        // Absolute and relative tolerances for the linear search
        double reltol = rteps * (xnorm + 1.0) / pe;
        double abstol = -DBL_EPSILON * (1.0 + Math.abs(f)) / (initGu - DBL_EPSILON);

        // Smallest allowable spacing between points in the linear search
        double tnytol = DBL_EPSILON * (xnorm + 1.0) / pe;

        double rtsmll = DBL_EPSILON;
        double big = 1.0 / (DBL_EPSILON * DBL_EPSILON);
        int numIterations = 0;

        GetPointCubic getPointCubic = new GetPointCubic(reltol, abstol, tnytol, ftol, eta, 1e-4, xBound, alpha, f, initGu);

        boolean requiresFurtherEvaluation = true;
        while (requiresFurtherEvaluation) {
            if (numIterations == MAX_ITERATIONS) {
                throw new MinimizationError("Line search failed. Max iterations reached");
            }
            numIterations++;

            double ualpha = getPointCubic.xmin() + getPointCubic.resultStep();
            for (int i = 0; i < x.length; i++) {
                temp[i] = x[i] + ualpha * searchDirection[i];
            }

            // Function evaluation
            unscalex(temp, xscale, xoffset);
            ArrayMath.clip(temp, low, up);
            double fu = functionEvaluator.evaluate(temp, tempgfull) * fscale;

            ArrayMath.copy(tempgfull, temp);
            scaleg(temp, xscale, fscale);
            double gu = ArrayMath.dotProduct(temp, searchDirection);

            requiresFurtherEvaluation = getPointCubic.iterate(big, rtsmll, fu, gu);

            // New best point
            if (getPointCubic.xmin() == ualpha) {
                ArrayMath.copy(tempgfull, newgfull);
            }
        } 

        double[] newX = new double[x.length];
        ArrayMath.copy(x, newX);
        ArrayMath.axPlusY(getPointCubic.xmin(), searchDirection, newX);
        return new LineSearchResult(getPointCubic.xmin(), newX, getPointCubic.fmin(), newgfull);
    }

}
