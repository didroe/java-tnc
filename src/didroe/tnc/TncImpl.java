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

import java.text.Normalizer;
import java.util.EnumSet;
import java.util.Iterator;

/**
 *
 * @author Did
 */
public class TncImpl {

    /**
     * Verbosity level
     */
    public enum tnc_message {

        /**
         * No messages
         */
        TNC_MSG_NONE,
        /**
         * One line per iteration
         */
        TNC_MSG_ITER,
        /**
         * Informational messages
         */
        TNC_MSG_INFO,
        /**
         * Exit reasons
         */
        TNC_MSG_EXIT;

        /**
         * @return Set of all messages
         */
        public static EnumSet<tnc_message> allMessages() {
            EnumSet<tnc_message> result = EnumSet.of(TNC_MSG_ITER, TNC_MSG_INFO,
                    TNC_MSG_EXIT);
            return result;
        }
    }

    /**
     * Possible return values for tnc
     */
    public enum tnc_rc {

        /**
         * Invalid parameters (n&lt;0)
         */
        TNC_EINVAL("Invalid parameters (n<0)"),
        /**
         * Infeasible (low bound &gt; up bound)
         */
        TNC_INFEASIBLE("Infeasible (low bound > up bound)"),
        /**
         * Local minima reach (|pg| ~= 0)
         */
        TNC_LOCALMINIMUM("Local minima reach (|pg| ~= 0)"),
        /**
         * Converged (|f_n-f_(n-1)| ~= 0)
         */
        TNC_FCONVERGED("Converged (|f_n-f_(n-1)| ~= 0)"),
        /**
         * Converged (|x_n-x_(n-1)| ~= 0)
         */
        TNC_XCONVERGED("Converged (|x_n-x_(n-1)| ~= 0)"),
        /**
         * Max. number of function evaluations reach
         */
        TNC_MAXFUN("Maximum number of function evaluations reached"),
        /**
         * Linear search failed
         */
        TNC_LSFAIL("Linear search failed"),
        /**
         * All lower bounds are equal to the upper bounds
         */
        TNC_CONSTANT("All lower bounds are equal to the upper bounds"),
        /**
         * Unable to progress
         */
        TNC_NOPROGRESS("Unable to progress"),
        /**
         * User requested end of minization
         */
        TNC_USERABORT("User requested end of minimization");

        private String message;

        public String getMessage() {
            return message;
        }

        private tnc_rc(String message) {
            this.message = message;
        }
    }

    /**
     * A callback function accepting x as input parameter along with the state
     * pointer.
     */
    public interface TncCallback {

        public void tnc_callback(double[] x);
    }

    /**
     * getptc return codes
     */
    public enum getptc_rc {

        /**
         * Suitable point found
         */
        GETPTC_OK,
        /**
         * Function evaluation required
         */
        GETPTC_EVAL,
        /**
         * Bad input values
         */
        GETPTC_EINVAL,
        /**
         * No suitable point found
         */
        GETPTC_FAIL;
    }

    /**
     * linearSearch return codes
     */
    public enum ls_rc {

        /**
         * Suitable point found
         */
        LS_OK,
        /**
         * Max. number of function evaluations reach
         */
        LS_MAXFUN,
        /**
         * No suitable point found
         */
        LS_FAIL,
        /**
         * User requested end of minimization
         */
        LS_USERABORT,
        /**
         * Memory allocation failed
         */
        LS_ENOMEM;
    }

    public static final double DBL_EPSILON = Math.ulp(1.0);

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

    public tnc_rc tnc(int n, double[] x, double[] g, TncFunction function,
            double[] low, double[] up, double[] scale,
            double[] offset, EnumSet<tnc_message> messages, int maxCGit, int maxnfeval,
            double eta, double stepmx, double accuracy, double fmin,
            double ftol, double xtol, double pgtol, double rescale,
            TncCallback callback, TncRef ref) {
        /*int rc, frc, i, nc, nfeval_local, free_low = false,
         free_up = false, free_g = false;
         double *xscale = NULL, fscale, rteps, *xoffset = NULL;*/

        ref.nfeval = 0;
        tnc_rc rc;

        cleanup:
        {
            /* Check for errors in the input parameters */
            if (n == 0) {
                rc = tnc_rc.TNC_CONSTANT;
                break cleanup;
            }

            if (n < 0) {
                rc = tnc_rc.TNC_EINVAL;
                break cleanup;
            }

            /* Check bounds arrays */
            if (low == null) {
                low = new double[n];
                for (int i = 0; i < n; i++) {
                    low[i] = Double.NEGATIVE_INFINITY;
                }
            }

            if (up == null) {
                up = new double[n];
                for (int i = 0; i < n; i++) {
                    up[i] = Double.POSITIVE_INFINITY;
                }
            }

            /* Coherency check */
            for (int i = 0; i < n; i++) {
                if (low[i] > up[i]) {
                    rc = tnc_rc.TNC_INFEASIBLE;
                    break cleanup;
                }
            }

            /* Coerce x into bounds */
            ArrayMath.clip(x, low, up);

            if (maxnfeval < 1) {
                rc = tnc_rc.TNC_MAXFUN;
                break cleanup;
            }

            /* Allocate g if necessary */
            if (g == null) {
                g = new double[n];
            }

            /* Initial function evaluation */
            ref.f = function.evaluate(x, g);
            ref.nfeval++;

            /* Constant problem ? */
            int nc = 0;
            for (int i = 0; i < n; i++) {
                if ((low[i] == up[i]) || (scale != null && scale[i] == 0.0)) {
                    nc++;
                }
            }

            if (nc == n) {
                rc = tnc_rc.TNC_CONSTANT;
                break cleanup;
            }

            /* Scaling parameters */
            double[] xscale = new double[n];
            double[] xoffset = new double[n];
            double fscale = 1.0;

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
            tmpMinRef.f = ref.f;
            tmpMinRef.fscale = fscale;
            tmpMinRef.nfeval = ref.nfeval;
            tmpMinRef.niter = ref.niter;
            rc = tnc_minimize(n, x, g, function, xscale, xoffset, low, up,
                    messages, maxCGit, maxnfeval, eta, stepmx,
                    accuracy, fmin, ftol, xtol, pgtol, rescale,
                    callback, tmpMinRef);
            ref.f = tmpMinRef.f;
            fscale = tmpMinRef.fscale;
            ref.nfeval = tmpMinRef.nfeval;
            ref.niter = tmpMinRef.niter;
        }
        if (messages.contains(tnc_message.TNC_MSG_EXIT)) {
            System.err.format("tnc: %s\n", rc.getMessage());
        }

        return rc;
    }


    /**
     * Unscale x
     */
    public void unscalex(int n, double[] x, double[] xscale, double xoffset[]) {
        for (int i = 0; i < n; i++) {
            x[i] = x[i] * xscale[i] + xoffset[i];
        }
    }

    /**
     * Scale x
     */
    public void scalex(int n, double[] x, double[] xscale, double[] xoffset) {
        for (int i = 0; i < n; i++) {
            if (xscale[i] > 0.0) {
                x[i] = (x[i] - xoffset[i]) / xscale[i];
            }
        }
    }

    /**
     * Scale g
     */
    public void scaleg(int n, double[] g, double[] xscale, double fscale) {
        for (int i = 0; i < n; i++) {
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
        public double fscale;
        public int nfeval;
        public int niter;
    }

    /**
     * This routine is a bounds-constrained truncated-newton method. the
     * truncated-newton method is preconditioned by a limited-memory
     * quasi-newton method (this preconditioning strategy is developed in this
     * routine) with a further diagonal scaling (see routine diagonalscaling).
     */
    public tnc_rc tnc_minimize(int n, double[] x, double[] gfull,
            TncFunction function, double[] xscale, double[] xoffset,
            double[] low, double[] up, EnumSet<tnc_message> messages,
            int maxCGit, int maxnfeval, double eta, double stepmx,
            double accuracy, double fmin, double ftol, double xtol,
            double pgtol, double rescale, TncCallback callback, tnc_minimizeRef ref) {
        /*double fLastReset, difnew, epsred, oldgtp, difold, oldf, xnorm, newscale,
         gnorm, ustpmax, fLastConstraint, spe, yrsr, yksk,
         *temp = NULL
         , *sk = NULL
         , *yk = NULL
         , *diagb = NULL
         , *sr = NULL
         ,
         *yr = NULL
         , *oldg = NULL
         , *pk = NULL
         , *g = NULL;*/
        double alpha = 0.0;         /* Default unused value */

        /*int i, icycle, oldnfeval, *pivot = NULL, frc;
         logical lreset, newcon, upd1, remcon;
         tnc_rc rc = TNC_ENOMEM;     /* Default error */
        ref.niter = 0;

        /* Allocate temporary vectors */
        double[] oldg = new double[n];
        double[] g = new double[n];
        double[] temp = new double[n];
        double[] diagb = new double[n];
        double[] pk = new double[n];
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

        /* Initial scaling */
        scalex(n, x, xscale, xoffset);
        ref.f *= ref.fscale;

        /* initial pivot calculation */
        setConstraints(n, x, pivot, xscale, xoffset, low, up);

        ArrayMath.copy(gfull, g);
        scaleg(n, g, xscale, ref.fscale);

        /* Test the lagrange multipliers to see if they are non-negative. */
        for (int i = 0; i < n; i++) {
            if (-pivot[i] * g[i] < 0.0) {
                pivot[i] = 0;
            }
        }

        project(n, g, pivot);

        /* Set initial values to other parameters */
        double gnorm = ArrayMath.euclidianNorm(g);

        double fLastConstraint = ref.f;       /* Value at last constraint */

        double fLastReset = ref.f;            /* Value at last reset */

        if (messages.contains(tnc_message.TNC_MSG_ITER)) {
            System.err.format("  NIT   NF   F                       GTG\n");
            printCurrentIteration(n, ref.f / ref.fscale, gfull,
                    ref.niter, ref.nfeval, pivot);
        }

        /* Set the diagonal of the approximate hessian to unity. */
        for (int i = 0; i < n; i++) {
            diagb[i] = 1.0;
        }

        tnc_rc rc;
        /* Start of main iterative loop */
        while (true) {
            /* Local minimum test */
            if (ArrayMath.euclidianNorm(g) <= pgtol * ref.fscale) {
                /* |PG| == 0.0 => local minimum */
                ArrayMath.copy(gfull, g);
                project(n, g, pivot);
                if (messages.contains(tnc_message.TNC_MSG_INFO)) {
                    System.err.format("tnc: |pg| = %g -> local minimum\n",
                            ArrayMath.euclidianNorm(g) / ref.fscale);
                }
                rc = tnc_rc.TNC_LOCALMINIMUM;
                break;
            }

            /* Terminate if more than maxnfeval evaluations have been made */
            if (ref.nfeval >= maxnfeval) {
                rc = tnc_rc.TNC_MAXFUN;
                break;
            }

            /* Rescale function if necessary */
            double newscale = ArrayMath.euclidianNorm(g);
            if ((newscale > DBL_EPSILON) && (Math.abs(Math.log10(newscale)) > rescale)) {
                newscale = 1.0 / newscale;

                ref.f *= newscale;
                ref.fscale *= newscale;
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

                if (messages.contains(tnc_message.TNC_MSG_INFO)) {
                    System.err.format("tnc: fscale = %g\n", ref.fscale);
                }
            }

            ArrayMath.copy(x, temp);
            project(n, temp, pivot);
            double xnorm = ArrayMath.euclidianNorm(temp);
            int oldnfeval = ref.nfeval;

            /* Compute the new search direction */
            tnc_directionRef tmpDirRef = new tnc_directionRef();
            tmpDirRef.diagb = diagb;
            tmpDirRef.nfeval = ref.nfeval;
            tmpDirRef.pivot = pivot;
            tmpDirRef.sk = sk;
            tmpDirRef.sr = sr;
            tmpDirRef.x = x;
            tmpDirRef.yk = yk;
            tmpDirRef.yr = yr;
            tmpDirRef.zsol = pk;
            tnc_rc frc = tnc_direction(g, n, maxCGit, maxnfeval, upd1, yksk, yrsr,
                    lreset, function, xscale, xoffset, ref.fscale, accuracy, gnorm,
                    xnorm, low, up, tmpDirRef);
            diagb = tmpDirRef.diagb;
            ref.nfeval = tmpDirRef.nfeval;
            pivot = tmpDirRef.pivot;
            sk = tmpDirRef.sk;
            sr = tmpDirRef.sr;
            x = tmpDirRef.x;
            yk = tmpDirRef.yk;
            yr = tmpDirRef.yr;
            pk = tmpDirRef.zsol;

            // TODO: This is not right
            if (frc != tnc_rc.TNC_LOCALMINIMUM) {
                rc = tnc_rc.TNC_USERABORT;
                break;
            }

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
            double oldgtp = ArrayMath.dotProduct(pk, g);

            /* Maximum unconstrained step length */
            double ustpmax = stepmx / (ArrayMath.euclidianNorm(pk) + DBL_EPSILON);

            /* Maximum constrained step length */
            double spe = stepMax(ustpmax, n, x, pk, pivot, low, up, xscale, xoffset);

            if (spe > 0.0) {
                /* Set the initial step length */
                alpha = initialStep(ref.f, fmin / ref.fscale, oldgtp, spe);

                /* Perform the linear search */
                linearSearchRef tmpLinSRef = new linearSearchRef();
                tmpLinSRef.alpha = alpha;
                tmpLinSRef.f = ref.f;
                tmpLinSRef.nfeval = ref.nfeval;
                ls_rc lsrc = linearSearch(n, function, low, up,
                        xscale, xoffset, ref.fscale, pivot,
                        eta, ftol, spe, pk, x, gfull,
                        maxnfeval, tmpLinSRef);
                alpha = tmpLinSRef.alpha;
                ref.f = tmpLinSRef.f;
                ref.nfeval = tmpLinSRef.nfeval;

                if (lsrc == ls_rc.LS_USERABORT) {
                    rc = tnc_rc.TNC_USERABORT;
                    break;
                }

                if (lsrc == ls_rc.LS_FAIL) {
                    rc = tnc_rc.TNC_LSFAIL;
                    break;
                }

                /* If we went up to the maximum unconstrained step, increase it */
                if (alpha >= 0.9 * ustpmax) {
                    stepmx *= 1e2;
                    if (messages.contains(tnc_message.TNC_MSG_INFO)) {
                        System.err.format("tnc: stepmx = %g\n", stepmx);
                    }
                }

                /* If we went up to the maximum constrained step,
                 a new constraint was encountered */
                if (alpha - spe >= -DBL_EPSILON * 10.0) {
                    newcon = true;
                } else {
                    /* Break if the linear search has failed to find a lower point */
                    if (lsrc != ls_rc.LS_OK) {
                        if (lsrc == ls_rc.LS_MAXFUN) {
                            rc = tnc_rc.TNC_MAXFUN;
                        } else {
                            rc = tnc_rc.TNC_LSFAIL;
                        }
                        break;
                    }
                    newcon = false;
                }
            } else {
                /* Maximum constrained step == 0.0 => new constraint */
                newcon = true;
            }

            if (newcon) {
                if (!addConstraint(n, x, pk, pivot, low, up, xscale, xoffset)) {
                    if (ref.nfeval == oldnfeval) {
                        rc = tnc_rc.TNC_NOPROGRESS;
                        break;
                    }
                }
                fLastConstraint = ref.f;
            }

            ref.niter++;

            /* Invoke the callback function */
            if (callback != null) {
                unscalex(n, x, xscale, xoffset);
                callback.tnc_callback(x);
                scalex(n, x, xscale, xoffset);
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
            scaleg(n, g, xscale, ref.fscale);

            ArrayMath.copy(g, temp);
            project(n, temp, pivot);
            gnorm = ArrayMath.euclidianNorm(temp);

            /* Reset pivot */
            boolean remcon = removeConstraint(oldgtp, gnorm, pgtol * ref.fscale, ref.f,
                    fLastConstraint, g, pivot, n);

            /* If a constraint is removed */
            if (remcon) {
                /* Recalculate gnorm and reset fLastConstraint */
                ArrayMath.copy(g, temp);
                project(n, temp, pivot);
                gnorm = ArrayMath.euclidianNorm(temp);
                fLastConstraint = ref.f;
            }

            if (!remcon && !newcon) {
                /* No constraint removed & no new constraint : tests for convergence */
                if (Math.abs(difnew) <= ftol * ref.fscale) {
                    if (messages.contains(tnc_message.TNC_MSG_INFO)) {
                        System.err.format("tnc: |fn-fn-1] = %g -> convergence\n",
                                Math.abs(difnew) / ref.fscale);
                    }
                    rc = tnc_rc.TNC_FCONVERGED;
                    break;
                }
                if (alpha * ArrayMath.euclidianNorm(pk) <= xtol) {
                    if (messages.contains(tnc_message.TNC_MSG_INFO)) {
                        System.err.format("tnc: |xn-xn-1] = %g -> convergence\n",
                                alpha * ArrayMath.euclidianNorm(pk));
                    }
                    rc = tnc_rc.TNC_XCONVERGED;
                    break;
                }
            }

            project(n, g, pivot);

            if (messages.contains(tnc_message.TNC_MSG_ITER)) {
                printCurrentIteration(n, ref.f / ref.fscale, gfull,
                        ref.niter, ref.nfeval, pivot);
            }

            /* Compute the change in the iterates and the corresponding change in the
             gradients */
            if (!newcon) {
                for (int i = 0; i < n; i++) {
                    yk[i] = g[i] - oldg[i];
                    sk[i] = alpha * pk[i];
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

        if (messages.contains(tnc_message.TNC_MSG_ITER)) {
            printCurrentIteration(n, ref.f / ref.fscale, gfull,
                    ref.niter, ref.nfeval, pivot);
        }

        /* Unscaling */
        unscalex(n, x, xscale, xoffset);
        ArrayMath.clip(x, low, up);
        ref.f /= ref.fscale;

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

        System.err.format(" %4d %4d %22.15E  %15.8E\n", niter, nfeval, f, gtg);
    }

    /**
     * Set x[i] = 0.0 if direction i is currently constrained
     */
    public void project(int n, double[] x, int[] pivot) {
        for (int i = 0; i < n; i++) {
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
        public int nfeval;
        public double[] sk;
        public double[] yk;
        public double[] sr;
        public double[] yr;
        public int[] pivot;
    }

    /**
     * This routine performs a preconditioned conjugate-gradient iteration in
     * order to solve the newton equations for a search direction for a
     * truncated-newton algorithm. When the value of the quadratic model is
     * sufficiently reduced, the iteration is terminated.
     */
    public tnc_rc tnc_direction(double[] g, int n, int maxCGit, int maxnfeval,
            boolean upd1, double yksk, double yrsr, boolean lreset,
            TncFunction function, double[] xscale, double[] xoffset,
            double fscale, double accuracy, double gnorm, double xnorm,
            double[] low, double[] up, tnc_directionRef ref) {
        /* No CG it. => dir = -grad */
        if (maxCGit == 0) {
            ArrayMath.copy(g, ref.zsol);
            ArrayMath.negate(ref.zsol);
            project(n, ref.zsol, ref.pivot);
            // TODO : Not right
            return tnc_rc.TNC_LOCALMINIMUM;
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
            project(n, r, ref.pivot);
            msolve(r, zk, n, ref.sk, ref.yk, ref.diagb, ref.sr, ref.yr, upd1, yksk, yrsr,
                    lreset);
            project(n, zk, ref.pivot);
            double rz = ArrayMath.dotProduct(r, zk);

            if ((rz / rhsnrm < tol) || (ref.nfeval >= (maxnfeval - 1))) {
                /* Truncate algorithm in case of an emergency
                 or too many function evaluations */
                if (k == 0) {
                    ArrayMath.copy(g, ref.zsol);
                    ArrayMath.negate(ref.zsol);
                    project(n, ref.zsol, ref.pivot);
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

            project(n, v, ref.pivot);
            tnc_rc frc = hessianTimesVector(v, gv, n, ref.x, g, function,
                    xscale, xoffset, fscale, accuracy, xnorm,
                    low, up);
            ref.nfeval++;
            // TODO: Not right
            if (frc != tnc_rc.TNC_LOCALMINIMUM) {
                return frc;
            }
            project(n, gv, ref.pivot);

            double vgv = ArrayMath.dotProduct(v, gv);
            if (vgv / rhsnrm < tol) {
                /* Truncate algorithm in case of an emergency */
                if (k == 0) {
                    msolve(g, ref.zsol, n, ref.sk, ref.yk, ref.diagb, ref.sr, ref.yr, upd1, yksk,
                            yrsr, lreset);
                    ArrayMath.negate(ref.zsol);
                    project(n, ref.zsol, ref.pivot);
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

        // TODO : Not right
        return tnc_rc.TNC_LOCALMINIMUM;
    }

    /**
     * Update the preconditioning matrix based on a diagonal version of the bfgs
     * quasi-newton update.
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
     * Returns the length of the initial step to be taken along the vector p in
     * the next linear search.
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
    public tnc_rc hessianTimesVector(double[] v, double[] gv, int n, double[] x,
            double[] g, TncFunction function, double[] xscale, double[] xoffset,
            double fscale, double accuracy, double xnorm, double[] low,
            double[] up) {
        double[] xv = new double[n];

        double delta = accuracy * (xnorm + 1.0);
        for (int i = 0; i < n; i++) {
            xv[i] = x[i] + delta * v[i];
        }

        unscalex(n, xv, xscale, xoffset);
        ArrayMath.clip(xv, low, up);
        double f = function.evaluate(xv, gv);

        scaleg(n, gv, xscale, fscale);

        double dinv = 1.0 / delta;
        for (int i = 0; i < n; i++) {
            gv[i] = (gv[i] - g[i]) * dinv;
        }

        projectConstants(n, gv, xscale);

        return tnc_rc.TNC_LOCALMINIMUM;
    }

    /**
     * This routine acts as a preconditioning step for the linear
     * conjugate-gradient routine. It is also the method of computing the search
     * direction from the gradient for the non-linear conjugate-gradient code.
     * It represents a two-step self-scaled bfgs formula.
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

    public class linearSearchRef {

        public double f;
        public double alpha;
        public int nfeval;
    }

    /**
     * Line search algorithm of gill and murray
     */
    public ls_rc linearSearch(int n, TncFunction function, double[] low,
            double[] up, double[] xscale, double[] xoffset, double fscale,
            int[] pivot, double eta, double ftol, double xbnd, double[] p,
            double[] x, double gfull[], int maxnfeval, linearSearchRef ref) {
        int maxlsit = 64;

        double[] temp = new double[n];
        double[] tempgfull = new double[n];
        double[] newgfull = new double[n];

        ArrayMath.copy(gfull, temp);
        scaleg(n, temp, xscale, fscale);
        double gu = ArrayMath.dotProduct(temp, p);

        ArrayMath.copy(x, temp);
        project(n, temp, pivot);
        double xnorm = ArrayMath.euclidianNorm(temp);

        /* Compute the absolute and relative tolerances for the linear search */
        double rteps = Math.sqrt(DBL_EPSILON);
        double pe = ArrayMath.euclidianNorm(p) + DBL_EPSILON;
        double reltol = rteps * (xnorm + 1.0) / pe;
        double abstol = -DBL_EPSILON * (1.0 + Math.abs(ref.f)) / (gu - DBL_EPSILON);

        /* Compute the smallest allowable spacing between points in the linear
         search */
        double tnytol = DBL_EPSILON * (xnorm + 1.0) / pe;

        double rtsmll = DBL_EPSILON;
        double big = 1.0 / (DBL_EPSILON * DBL_EPSILON);
        int itcnt = 0;

        /* Set the estimated relative precision in f(x). */
        double fpresn = ftol;

        double u = ref.alpha;
        double fu = ref.f;
        double fmin = ref.f;
        double rmu = 1e-4;

        /* Setup */
        getptcInitRef tmpPInitRef = new getptcInitRef();
        tmpPInitRef.abstol = abstol;
        tmpPInitRef.fmin = fmin;
        tmpPInitRef.fu = fu;
        tmpPInitRef.gu = gu;
        tmpPInitRef.reltol = reltol;
        tmpPInitRef.u = u;
        tmpPInitRef.xmin = ref.alpha;
        getptc_rc itest = getptcInit(tnytol, eta, rmu, xbnd, tmpPInitRef);
        double a = tmpPInitRef.a;
        abstol = tmpPInitRef.abstol;
        double b = tmpPInitRef.b;
        double b1 = tmpPInitRef.b1;
        boolean braktd = tmpPInitRef.braktd;
        double e = tmpPInitRef.e;
        double factor = tmpPInitRef.factor;
        fmin = tmpPInitRef.fmin;
        fu = tmpPInitRef.fu;
        double fw = tmpPInitRef.fw;
        double gmin = tmpPInitRef.gmin;
        double gtest1 = tmpPInitRef.gtest1;
        double gtest2 = tmpPInitRef.gtest2;
        gu = tmpPInitRef.gu;
        double gw = tmpPInitRef.gw;
        double oldf = tmpPInitRef.oldf;
        reltol = tmpPInitRef.reltol;
        double scxbnd = tmpPInitRef.scxbnd;
        double step = tmpPInitRef.step;
        double tol = tmpPInitRef.tol;
        u = tmpPInitRef.u;
        ref.alpha = tmpPInitRef.xmin;
        double xw = tmpPInitRef.xw;
        
        /* If itest == GETPTC_EVAL, the algorithm requires the function value to be
         calculated */
        while (itest == getptc_rc.GETPTC_EVAL) {
            /* Test for too many iterations or too many function evals */
            if ((++itcnt > maxlsit) || (ref.nfeval >= maxnfeval)) {
                break;
            }

            double ualpha = ref.alpha + u;
            for (int i = 0; i < n; i++) {
                temp[i] = x[i] + ualpha * p[i];
            }

            /* Function evaluation */
            unscalex(n, temp, xscale, xoffset);
            ArrayMath.clip(temp, low, up);

            fu = function.evaluate(temp, tempgfull);
            ref.nfeval++;

            fu *= fscale;

            ArrayMath.copy(tempgfull, temp);
            scaleg(n, temp, xscale, fscale);
            gu = ArrayMath.dotProduct(temp, p);

            getptcIterRef tmpPIterRef = new getptcIterRef();
            tmpPIterRef.a = a;
            tmpPIterRef.abstol = abstol;
            tmpPIterRef.b = b;
            tmpPIterRef.b1 = b1;
            tmpPIterRef.braktd = braktd;
            tmpPIterRef.e = e;
            tmpPIterRef.factor = factor;
            tmpPIterRef.fmin = fmin;
            tmpPIterRef.fu = fu;
            tmpPIterRef.fw = fw;
            tmpPIterRef.gmin = gmin;
            tmpPIterRef.gtest1 = gtest1;
            tmpPIterRef.gtest2 = gtest2;
            tmpPIterRef.gu = gu;
            tmpPIterRef.gw = gw;
            tmpPIterRef.oldf = oldf;
            tmpPIterRef.reltol = reltol;
            tmpPIterRef.scxbnd = scxbnd;
            tmpPIterRef.step = step;
            tmpPIterRef.tol = tol;
            tmpPIterRef.u = u;
            tmpPIterRef.xmin = ref.alpha;
            tmpPIterRef.xw = xw;
            itest = getptcIter(big, rtsmll, tnytol, fpresn, xbnd, tmpPIterRef);
            a = tmpPIterRef.a;
            abstol = tmpPIterRef.abstol;
            b = tmpPIterRef.b;
            b1 = tmpPIterRef.b1;
            braktd = tmpPIterRef.braktd;
            e = tmpPIterRef.e;
            factor = tmpPIterRef.factor;
            fmin = tmpPIterRef.fmin;
            fu = tmpPIterRef.fu;
            fw = tmpPIterRef.fw;
            gmin = tmpPIterRef.gmin;
            gtest1 = tmpPIterRef.gtest1;
            gtest2 = tmpPIterRef.gtest2;
            gu = tmpPIterRef.gu;
            gw = tmpPIterRef.gw;
            oldf = tmpPIterRef.oldf;
            reltol = tmpPIterRef.reltol;
            scxbnd = tmpPIterRef.scxbnd;
            step = tmpPIterRef.step;
            tol = tmpPIterRef.tol;
            u = tmpPIterRef.u;
            ref.alpha = tmpPIterRef.xmin;
            xw = tmpPIterRef.xw;

            /* New best point ? */
            if (ref.alpha == ualpha) {
                ArrayMath.copy(tempgfull, newgfull);
            }
        }

        if (itest == getptc_rc.GETPTC_OK) {
            /* A successful search has been made */
            ref.f = fmin;
            ArrayMath.axPlusY(ref.alpha, p, x);
            ArrayMath.copy(newgfull, gfull);
            return ls_rc.LS_OK;
        } else if (itcnt > maxlsit) {
            /* Too many iterations ? */
            return ls_rc.LS_FAIL;
        } else if (itest != getptc_rc.GETPTC_EVAL) {
            /* If itest=GETPTC_FAIL or GETPTC_EINVAL a lower point could not be found */
            return ls_rc.LS_FAIL;
        } else {
            /* Too many function evaluations */
            return ls_rc.LS_MAXFUN;
        }
    }

    public class getptcInitRef {

        public double reltol;
        public double abstol;
        public double u;
        public double fu;
        public double gu;
        public double xmin;
        public double fmin;
        public double gmin;
        public double xw;
        public double fw;
        public double gw;
        public double a;
        public double b;
        public double oldf;
        public double b1;
        public double scxbnd;
        public double e;
        public double step;
        public double factor;
        public boolean braktd;
        public double gtest1;
        public double gtest2;
        public double tol;
    }

    /**
     * getptc, an algorithm for finding a steplength, called repeatedly by
     * routines which require a step length to be computed using cubic
     * interpolation. The parameters contain information about the interval in
     * which a lower point is to be found and from this getptc computes a point
     * at which the function can be evaluated by the calling program.
     */
    public getptc_rc getptcInit(double tnytol, double eta, double rmu,
            double xbnd, getptcInitRef ref) {
        /* Check input parameters */
        if (ref.u <= 0.0 || xbnd <= tnytol || ref.gu > 0.0) {
            return getptc_rc.GETPTC_EINVAL;
        }
        if (xbnd < ref.abstol) {
            ref.abstol = xbnd;
        }
        ref.tol = ref.abstol;

        /* a and b define the interval of uncertainty, x and xw are points */
        /* with lowest and second lowest function values so far obtained. */
        /* initialize a,smin,xw at origin and corresponding values of */
        /* function and projection of the gradient along direction of search */
        /* at values for latest estimate at minimum. */
        ref.a = 0.0;
        ref.xw = 0.0;
        ref.xmin = 0.0;
        ref.oldf = ref.fu;
        ref.fmin = ref.fu;
        ref.fw = ref.fu;
        ref.gw = ref.gu;
        ref.gmin = ref.gu;
        ref.step = ref.u;
        ref.factor = 5.0;

        /* The minimum has not yet been bracketed. */
        ref.braktd = false;

        /* Set up xbnd as a bound on the step to be taken. (xbnd is not computed */
        /* explicitly but scxbnd is its scaled value.) Set the upper bound */
        /* on the interval of uncertainty initially to xbnd + tol(xbnd). */
        ref.scxbnd = xbnd;
        ref.b = ref.scxbnd + ref.reltol * Math.abs(ref.scxbnd) + ref.abstol;
        ref.e = ref.b + ref.b;
        ref.b1 = ref.b;

        /* Compute the constants required for the two convergence criteria. */
        ref.gtest1 = -rmu * ref.gu;
        ref.gtest2 = -eta * ref.gu;

        /* If the step is too large, replace by the scaled bound (so as to */
        /* compute the new point on the boundary). */
        if (ref.step >= ref.scxbnd) {
            ref.step = ref.scxbnd;
            /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
            ref.scxbnd -= (ref.reltol * Math.abs(xbnd) + ref.abstol) / (1.0 + ref.reltol);
        }
        ref.u = ref.step;
        if (Math.abs(ref.step) < ref.tol && ref.step < 0.0) {
            ref.u = -ref.tol;
        }
        if (Math.abs(ref.step) < ref.tol && ref.step >= 0.0) {
            ref.u = ref.tol;
        }
        return getptc_rc.GETPTC_EVAL;
    }

    public class getptcIterRef {

        public double reltol;
        public double abstol;
        public double u;
        public double fu;
        public double gu;
        public double xmin;
        public double fmin;
        public double gmin;
        public double xw;
        public double fw;
        public double gw;
        public double a;
        public double b;
        public double oldf;
        public double b1;
        public double scxbnd;
        public double e;
        public double step;
        public double factor;
        public boolean braktd;
        public double gtest1;
        public double gtest2;
        public double tol;
    }

    public getptc_rc getptcIter(double big, double rtsmll, double tnytol,
            double fpresn, double xbnd, getptcIterRef ref) {
        ConvergenceCheck:
        {
            /* Update a,b,xw, and xmin */
            if (ref.fu <= ref.fmin) {
                /* If function value not increased, new point becomes next */
                /* origin and other points are scaled accordingly. */
                double chordu = ref.oldf - (ref.xmin + ref.u) * ref.gtest1;
                if (ref.fu > chordu) {
                    /* The new function value does not satisfy the sufficient decrease */
                    /* criterion. prepare to move the upper bound to this point and */
                    /* force the interpolation scheme to either bisect the interval of */
                    /* uncertainty or take the linear interpolation step which estimates */
                    /* the root of f(alpha)=chord(alpha). */

                    double chordm = ref.oldf - ref.xmin * ref.gtest1;
                    ref.gu = -ref.gmin;
                    double denom = chordm - ref.fmin;
                    if (Math.abs(denom) < 1e-15) {
                        denom = 1e-15;
                        if (chordm - ref.fmin < 0.0) {
                            denom = -denom;
                        }
                    }
                    if (ref.xmin != 0.0) {
                        ref.gu = ref.gmin * (chordu - ref.fu) / denom;
                    }
                    ref.fu = 0.5 * ref.u * (ref.gmin + ref.gu) + ref.fmin;
                    if (ref.fu < ref.fmin) {
                        ref.fu = ref.fmin;
                    }
                } else {
                    ref.fw = ref.fmin;
                    ref.fmin = ref.fu;
                    ref.gw = ref.gmin;
                    ref.gmin = ref.gu;
                    ref.xmin += ref.u;
                    ref.a -= ref.u;
                    ref.b -= ref.u;
                    ref.xw = -ref.u;
                    ref.scxbnd -= ref.u;
                    if (ref.gu <= 0.0) {
                        ref.a = 0.0;
                    } else {
                        ref.b = 0.0;
                        ref.braktd = true;
                    }
                    ref.tol = Math.abs(ref.xmin) * ref.reltol + ref.abstol;
                    break ConvergenceCheck;
                }
            }

            /* If function value increased, origin remains unchanged */
            /* but new point may now qualify as w. */
            if (ref.u < 0.0) {
                ref.a = ref.u;
            } else {
                ref.b = ref.u;
                ref.braktd = true;
            }
            ref.xw = ref.u;
            ref.fw = ref.fu;
            ref.gw = ref.gu;
        }

        double twotol = ref.tol + ref.tol;
        double xmidpt = 0.5 * (ref.a + ref.b);

        /* Check termination criteria */
        boolean convrg = (Math.abs(xmidpt) <= twotol - 0.5 * (ref.b - ref.a))
                || (Math.abs(ref.gmin) <= ref.gtest2 && ref.fmin < ref.oldf
                && ((Math.abs(ref.xmin - xbnd) > ref.tol) || (!ref.braktd)));
        if (convrg) {
            if (ref.xmin != 0.0) {
                return getptc_rc.GETPTC_OK;
            }

            /*
             * If the function has not been reduced, check to see that the relative
             * change in f(x) is consistent with the estimate of the delta-
             * unimodality constant, tol. If the change in f(x) is larger than
             * expected, reduce the value of tol.
             */
            if (Math.abs(ref.oldf - ref.fw) <= fpresn) {
                return getptc_rc.GETPTC_FAIL;
            }
            ref.tol = 0.1 * ref.tol;
            if (ref.tol < tnytol) {
                return getptc_rc.GETPTC_FAIL;
            }
            ref.reltol = 0.1 * ref.reltol;
            ref.abstol = 0.1 * ref.abstol;
            twotol = 0.1 * twotol;
        }

        /* Continue with the computation of a trial step length */
        double r = 0.0;
        double q = 0.0;
        double s = 0.0;
        MinimumFound:
        {
            if (Math.abs(ref.e) > ref.tol) {
                /* Fit cubic through xmin and xw */
                r = 3.0 * (ref.fmin - ref.fw) / ref.xw + ref.gmin + ref.gw;
                double absr = Math.abs(r);
                q = absr;
                if (ref.gw != 0.0 && ref.gmin != 0.0) {
                    /* Compute the square root of (r*r - gmin*gw) in a way
                     which avoids underflow and overflow. */
                    double abgw = Math.abs(ref.gw);
                    double abgmin = Math.abs(ref.gmin);
                    s = Math.sqrt(abgmin) * Math.sqrt(abgw);
                    if (ref.gw / abgw * ref.gmin > 0.0) {
                        if (r >= s || r <= -s) {
                            /* Compute the square root of r*r - s*s */
                            q = Math.sqrt(Math.abs(r + s)) * Math.sqrt(Math.abs(r - s));
                        } else {
                            r = 0.0;
                            q = 0.0;
                            break MinimumFound;
                        }
                    } else {
                        /* Compute the square root of r*r + s*s. */
                        double sumsq = 1.0;
                        double p = 0.0;
                        double scale;
                        if (absr >= s) {
                            /* There is a possibility of underflow. */
                            if (absr > rtsmll) {
                                p = absr * rtsmll;
                            }
                            if (s >= p) {
                                double value = s / absr;
                                sumsq = 1.0 + value * value;
                            }
                            scale = absr;
                        } else {
                            /* There is a possibility of overflow. */
                            if (s > rtsmll) {
                                p = s * rtsmll;
                            }
                            if (absr >= p) {
                                double value = absr / s;
                                sumsq = 1.0 + value * value;
                            }
                            scale = s;
                        }
                        sumsq = Math.sqrt(sumsq);
                        q = big;
                        if (scale < big / sumsq) {
                            q = scale * sumsq;
                        }
                    }
                }

                /* Compute the minimum of fitted cubic */
                if (ref.xw < 0.0) {
                    q = -q;
                }
                s = ref.xw * (ref.gmin - r - q);
                q = ref.gw - ref.gmin + q + q;
                if (q > 0.0) {
                    s = -s;
                }
                if (q <= 0.0) {
                    q = -q;
                }
                r = ref.e;
                if (ref.b1 != ref.step || ref.braktd) {
                    ref.e = ref.step;
                }
            }
        }

        /* Construct an artificial bound on the estimated steplength */
        double a1 = ref.a;
        ref.b1 = ref.b;
        ref.step = xmidpt;
        if ((!ref.braktd) || ((ref.a == 0.0 && ref.xw < 0.0) || (ref.b == 0.0 && ref.xw > 0.0))) {
            if (ref.braktd) {
                /* If the minimum is not bracketed by 0 and xw the step must lie
                 within (a1,b1). */
                double d1 = ref.xw;
                double d2 = ref.a;
                if (ref.a == 0.0) {
                    d2 = ref.b;
                }
                /* This line might be : */
                /* if (*a == 0.0) d2 = *e */
                ref.u = -d1 / d2;
                ref.step = 5.0 * d2 * (0.1 + 1.0 / ref.u) / 11.0;
                if (ref.u < 1.0) {
                    ref.step = 0.5 * d2 * Math.sqrt(ref.u);
                }
            } else {
                ref.step = -ref.factor * ref.xw;
                if (ref.step > ref.scxbnd) {
                    ref.step = ref.scxbnd;
                }
                if (ref.step != ref.scxbnd) {
                    ref.factor = 5.0 * ref.factor;
                }
            }
            /* If the minimum is bracketed by 0 and xw the step must lie within (a,b) */
            if (ref.step <= 0.0) {
                a1 = ref.step;
            }
            if (ref.step > 0.0) {
                ref.b1 = ref.step;
            }
        }

        /*
         * Reject the step obtained by interpolation if it lies outside the
         * required interval or it is greater than half the step obtained
         * during the last-but-one iteration.
         */
        if (Math.abs(s) <= Math.abs(0.5 * q * r) || s <= q * a1 || s >= q * ref.b1) {
            ref.e = ref.b - ref.a;
        } else {
            /* A cubic interpolation step */
            ref.step = s / q;

            /* The function must not be evaluated too close to a or b. */
            if (ref.step - ref.a < twotol || ref.b - ref.step < twotol) {
                if (xmidpt <= 0.0) {
                    ref.step = -ref.tol;
                } else {
                    ref.step = ref.tol;
                }
            }
        }

        /* If the step is too large, replace by the scaled bound (so as to */
        /* compute the new point on the boundary). */
        if (ref.step >= ref.scxbnd) {
            ref.step = ref.scxbnd;
            /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
            ref.scxbnd -= (ref.reltol * Math.abs(xbnd) + ref.abstol) / (1.0 + ref.reltol);
        }
        ref.u = ref.step;
        if (Math.abs(ref.step) < ref.tol && ref.step < 0.0) {
            ref.u = -ref.tol;
        }
        if (Math.abs(ref.step) < ref.tol && ref.step >= 0.0) {
            ref.u = ref.tol;
        }
        return getptc_rc.GETPTC_EVAL;
    }


}
