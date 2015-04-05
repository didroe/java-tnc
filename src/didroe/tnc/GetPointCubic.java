/*
 * See the "LICENSE" file for the full license governing this code.
 */
package didroe.tnc;

/**
 * An algorithm for finding a step-length, called repeatedly by routines which require a step length
 * to be computed using cubic interpolation. The parameters contain information about the interval
 * in which a lower point is to be found and from this the algorithm computes a point at which the
 * function can be evaluated by the calling program.
 *
 * @author Did
 */
class GetPointCubic {

    /**
     * The original function value.
     */
    private final double oldFunctionValue;

    /**
     * Bound on the step.
     */
    private final double stepBound;

    /**
     * Smallest allowable spacing (tolerance) between points.
     */
    private final double smallestSpacing;
    
    /**
     * Estimated relative precision in f(x)
     */
    private final double functionPrecision;
    
    private final double sufficientDecreaseCondition;
    private final double sufficientDescentCondition;

    private double relativeTolerance;
    private double absoluteTolerance;
    private double resultStep;
    private double xmin;
    private double fmin;
    private double gmin;
    private double xw;
    private double fw;
    private double gw;

    /**
     * Lower bound of interval containing the minimiser.
     */
    private double intervalLowerBound;
    /**
     * Upper bound of interval containing the minimiser.
     */
    private double intervalUpperBound;
    /**
     * Is the minimum bracketed by the interval? Starts off as false.
     */
    private boolean minIsBracketed = false;

    private double b1;
    private double e;
    private double step;
    private double factor;
    private double tolerance;

    /**
     * A bound on the step to be taken.
     * FIXME: What does "scaled" mean?
     */
    private double scaledStepBound;

    /**
     * FIXME: Better name
     */
    private void updateStepValues() {
        // If the step is too large, replace by the scaled bound (so as to compute
        // the new point on the boundary). 
        if (step >= scaledStepBound) {
            step = scaledStepBound;
            // Move the step bound to the left so that scaledStepBound + tol(stepBound) = stepBound 
            scaledStepBound -= (relativeTolerance * Math.abs(stepBound) + absoluteTolerance) 
                    / (1.0 + relativeTolerance);
        }

        if (Math.abs(step) < tolerance) {
            if (step < 0.0) {
                resultStep = -tolerance;
            } else {
                resultStep = tolerance;
            }
        } else {
            resultStep = step;
        }
    }

    /**
     * 
     * @param relativeTolerance
     * @param absoluteTolerance
     * @param smallestSpacing
     * @param eta
     * @param rmu
     * @param stepBound
     * @param initialStep
     * @param fu The current function value
     * @param gu Dot product of the search direction and gradient (angle of search?)
     */
    public GetPointCubic(double relativeTolerance, double absoluteTolerance, double smallestSpacing, 
            double functionPrecision, double eta, double rmu, double stepBound, 
            double initialStep, double fu, double gu) {
        if (initialStep <= 0.0) {
            throw new IllegalArgumentException("Initial step must be > 0");
        }
        
        if (stepBound <= smallestSpacing) {
            throw new IllegalArgumentException("Step bound must be > smallest allowable spacing");
        }
        
        if (gu > 0.0) {
            throw new IllegalArgumentException("Search direction must be a descent direction");
        }

        oldFunctionValue = fu;
        this.relativeTolerance = relativeTolerance;
        this.stepBound = stepBound;
        this.smallestSpacing = smallestSpacing;
        this.functionPrecision = functionPrecision;
        
        if (stepBound < absoluteTolerance) {
            this.absoluteTolerance = this.stepBound;
        } else {
            this.absoluteTolerance = absoluteTolerance;
        }
        tolerance = this.absoluteTolerance;

        /* x and xw are points */
        /* with lowest and second lowest function values so far obtained. */
        /* initialize a,smin,xw at origin and corresponding values of */
        /* function and projection of the gradient along direction of search */
        /* at values for latest estimate at minimum. */
        intervalLowerBound = 0.0;
        xw = 0.0;
        xmin = 0.0;
        fmin = fu;
        fw = fu;
        gw = gu;
        gmin = gu;
        step = initialStep;
        factor = 5.0;

        scaledStepBound = this.stepBound;
        intervalUpperBound = scaledStepBound + (this.relativeTolerance * Math.abs(scaledStepBound)) + this.absoluteTolerance;
        e = intervalUpperBound + intervalUpperBound;
        b1 = intervalUpperBound;

        // Compute the constants required for the two convergence criteria
        sufficientDecreaseCondition = -rmu * gu;
        sufficientDescentCondition = -eta * gu;

        updateStepValues();
    }

    /**
     * Some stuff that happens before the convergence check. 
     * FIXME: Better name.
     */
    private void beforeConvergenceCheck(double fu, double gu) {
        /* Update a,b,xw, and xmin */
        if (fu <= fmin) {
            /* If function value not increased, new point becomes next */
            /* origin and other points are scaled accordingly. */
            double chordu = oldFunctionValue - (xmin + resultStep) * sufficientDecreaseCondition;
            if (fu > chordu) {
                /* The new function value does not satisfy the sufficient decrease */
                /* criterion. prepare to move the upper bound to this point and */
                /* force the interpolation scheme to either bisect the interval of */
                /* uncertainty or take the linear interpolation step which estimates */
                /* the root of f(alpha)=chord(alpha). */

                double chordm = oldFunctionValue - (xmin * sufficientDecreaseCondition);
                gu = -gmin;
                double denom = chordm - fmin;
                if (Math.abs(denom) < 1e-15) {
                    denom = 1e-15;
                    if (chordm - fmin < 0.0) {
                        denom = -denom;
                    }
                }
                if (xmin != 0.0) {
                    gu = gmin * (chordu - fu) / denom;
                }
                fu = 0.5 * resultStep * (gmin + gu) + fmin;
                if (fu < fmin) {
                    fu = fmin;
                }
            } else {
                fw = fmin;
                fmin = fu;
                gw = gmin;
                gmin = gu;
                xmin += resultStep;
                intervalLowerBound -= resultStep;
                intervalUpperBound -= resultStep;
                xw = -resultStep;
                scaledStepBound -= resultStep;
                if (gu <= 0.0) {
                    intervalLowerBound = 0.0;
                } else {
                    intervalUpperBound = 0.0;
                    minIsBracketed = true;
                }
                tolerance = Math.abs(xmin) * relativeTolerance + absoluteTolerance;
                return;
            }
        }

        /* If function value increased, origin remains unchanged */
        /* but new point may now qualify as w. */
        if (resultStep < 0.0) {
            intervalLowerBound = resultStep;
        } else {
            intervalUpperBound = resultStep;
            minIsBracketed = true;
        }
        xw = resultStep;
        fw = fu;
        gw = gu;
    }

    public boolean iterate(double big, double rtsmll, 
            double fu, double gu) {
        beforeConvergenceCheck(fu, gu);

        /* From here is ConvergenceCheck */
        
        double twiceTolerance = tolerance + tolerance;
        double xmidpt = 0.5 * (intervalLowerBound + intervalUpperBound);

        /* Check termination criteria */
        boolean maybeConverged = (Math.abs(xmidpt) <= twiceTolerance - 0.5 * (intervalUpperBound - intervalLowerBound))
                || (Math.abs(gmin) <= sufficientDescentCondition && fmin < oldFunctionValue
                && ((Math.abs(xmin - stepBound) > tolerance) || (!minIsBracketed)));
        if (maybeConverged) {
            if (xmin != 0.0) {
                return false;
            }

            /*
             * If the function has not been reduced, check to see that the relative
             * change in f(x) is consistent with the estimate of the delta-
             * unimodality constant, tol. If the change in f(x) is larger than
             * expected, reduce the value of tol.
             */
            if (Math.abs(oldFunctionValue - fw) <= functionPrecision) {
                throw new MinimizationError("GetPointubic failed. |oldf - fw| <= fpresn");
            }
            tolerance = 0.1 * tolerance;
            if (tolerance < smallestSpacing) {
                throw new MinimizationError("GetPointubic failed. tol < tnytol");
            }
            relativeTolerance = 0.1 * relativeTolerance;
            absoluteTolerance = 0.1 * absoluteTolerance;
            twiceTolerance = 0.1 * twiceTolerance;
        }

        /* Continue with the computation of a trial step length */
        double r = 0.0;
        double q = 0.0;
        double s = 0.0;
        MinimumFound:
        {
            if (Math.abs(e) > tolerance) {
                /* Fit cubic through xmin and xw */
                r = 3.0 * (fmin - fw) / xw + gmin + gw;
                double absr = Math.abs(r);
                q = absr;
                if (gw != 0.0 && gmin != 0.0) {
                    /* Compute the square root of (r*r - gmin*gw) in a way
                     which avoids underflow and overflow. */
                    double abgw = Math.abs(gw);
                    double abgmin = Math.abs(gmin);
                    s = Math.sqrt(abgmin) * Math.sqrt(abgw);
                    if (gw / abgw * gmin > 0.0) {
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
                if (xw < 0.0) {
                    q = -q;
                }
                s = xw * (gmin - r - q);
                q = gw - gmin + q + q;
                if (q > 0.0) {
                    s = -s;
                }
                if (q <= 0.0) {
                    q = -q;
                }
                r = e;
                if (b1 != step || minIsBracketed) {
                    e = step;
                }
            }
        }

        /* Construct an artificial bound on the estimated steplength */
        double a1 = intervalLowerBound;
        b1 = intervalUpperBound;
        step = xmidpt;
        if ((!minIsBracketed) || ((intervalLowerBound == 0.0 && xw < 0.0) || (intervalUpperBound == 0.0 && xw > 0.0))) {
            if (minIsBracketed) {
                /* If the minimum is not bracketed by 0 and xw the step must lie
                 within (a1,b1). */
                double d1 = xw;
                double d2 = intervalLowerBound;
                if (intervalLowerBound == 0.0) {
                    d2 = intervalUpperBound;
                }
                /* This line might be : */
                /* if (*a == 0.0) d2 = *e */
                resultStep = -d1 / d2;
                step = 5.0 * d2 * (0.1 + 1.0 / resultStep) / 11.0;
                if (resultStep < 1.0) {
                    step = 0.5 * d2 * Math.sqrt(resultStep);
                }
            } else {
                step = -factor * xw;
                if (step > scaledStepBound) {
                    step = scaledStepBound;
                }
                if (step != scaledStepBound) {
                    factor = 5.0 * factor;
                }
            }
            /* If the minimum is bracketed by 0 and xw the step must lie within (a,b) */
            if (step <= 0.0) {
                a1 = step;
            }
            if (step > 0.0) {
                b1 = step;
            }
        }

        /*
         * Reject the step obtained by interpolation if it lies outside the
         * required interval or it is greater than half the step obtained
         * during the last-but-one iteration.
         */
        if (Math.abs(s) <= Math.abs(0.5 * q * r) || s <= q * a1 || s >= q * b1) {
            e = intervalUpperBound - intervalLowerBound;
        } else {
            /* A cubic interpolation step */
            step = s / q;

            /* The function must not be evaluated too close to a or b. */
            if (step - intervalLowerBound < twiceTolerance || intervalUpperBound - step < twiceTolerance) {
                if (xmidpt <= 0.0) {
                    step = -tolerance;
                } else {
                    step = tolerance;
                }
            }
        }

        updateStepValues();

        return true;
    }

    public double xmin() {
        return xmin;
    }

    public double fmin() {
        return fmin;
    }

    public double resultStep() {
        return resultStep;
    }

}
