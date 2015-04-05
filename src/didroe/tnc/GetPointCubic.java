/*
 * See the "LICENSE" file for the full license governing this code.
 */
package didroe.tnc;

/**
 * An algorithm for finding a step-length, called repeatedly by routines which require a
 * step length to be computed using cubic interpolation. The parameters contain information
 * about the interval in which a lower point is to be found and from this the algorithm 
 * computes a point at which the function can be evaluated by the calling program.
 * @author Did
 */
class GetPointCubic {
        private final double oldf;
        private final double gtest1;
        private final double gtest2;

        private double reltol;
        private double abstol;
        private double u;
        private double xmin;
        private double fmin;
        private double gmin;
        private double xw;
        private double fw;
        private double gw;
        private double a;
        private double b;
        private double b1;
        private double scxbnd;
        private double e;
        private double step;
        private double factor;
        private boolean braktd;
        private double tol;

    
    public GetPointCubic(double reltol, double abstol, double tnytol, double eta, double rmu, 
            double xbnd, double u, double fu, double gu) {
        if (u <= 0.0 || xbnd <= tnytol || gu > 0.0) {
            throw new IllegalArgumentException("Invalid inputs to GetPointCubic");
        }
        
        this.reltol = reltol;
        
        if (xbnd < abstol) {
            this.abstol = xbnd;
        } else {
            this.abstol = abstol;
        }
        tol = this.abstol;

        /* a and b define the interval of uncertainty, x and xw are points */
        /* with lowest and second lowest function values so far obtained. */
        /* initialize a,smin,xw at origin and corresponding values of */
        /* function and projection of the gradient along direction of search */
        /* at values for latest estimate at minimum. */
        a = 0.0;
        xw = 0.0;
        xmin = 0.0;
        oldf = fu;
        fmin = fu;
        fw = fu;
        gw = gu;
        gmin = gu;
        step = u;
        factor = 5.0;

        /* The minimum has not yet been bracketed. */
        braktd = false;

        /* Set up xbnd as a bound on the step to be taken. (xbnd is not computed */
        /* explicitly but scxbnd is its scaled value.) Set the upper bound */
        /* on the interval of uncertainty initially to xbnd + tol(xbnd). */
        scxbnd = xbnd;
        b = scxbnd + this.reltol * Math.abs(scxbnd) + this.abstol;
        e = b + b;
        b1 = b;

        /* Compute the constants required for the two convergence criteria. */
        gtest1 = -rmu * gu;
        gtest2 = -eta * gu;

        /* If the step is too large, replace by the scaled bound (so as to */
        /* compute the new point on the boundary). */
        if (step >= scxbnd) {
            step = scxbnd;
            /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
            scxbnd -= (this.reltol * Math.abs(xbnd) + this.abstol) / (1.0 + this.reltol);
        }
        
        if (Math.abs(step) < tol) {
            if (step < 0.0) {
                this.u = -tol;
            } else {
                this.u = tol;
            }
        } else {
            this.u = step;
        }
    }

    public boolean iterate(double big, double rtsmll, double tnytol,
            double fpresn, double xbnd, double fu, double gu) {
        ConvergenceCheck:
        {
            /* Update a,b,xw, and xmin */
            if (fu <= fmin) {
                /* If function value not increased, new point becomes next */
                /* origin and other points are scaled accordingly. */
                double chordu = oldf - (xmin + u) * gtest1;
                if (fu > chordu) {
                    /* The new function value does not satisfy the sufficient decrease */
                    /* criterion. prepare to move the upper bound to this point and */
                    /* force the interpolation scheme to either bisect the interval of */
                    /* uncertainty or take the linear interpolation step which estimates */
                    /* the root of f(alpha)=chord(alpha). */

                    double chordm = oldf - xmin * gtest1;
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
                    fu = 0.5 * u * (gmin + gu) + fmin;
                    if (fu < fmin) {
                        fu = fmin;
                    }
                } else {
                    fw = fmin;
                    fmin = fu;
                    gw = gmin;
                    gmin = gu;
                    xmin += u;
                    a -= u;
                    b -= u;
                    xw = -u;
                    scxbnd -= u;
                    if (gu <= 0.0) {
                        a = 0.0;
                    } else {
                        b = 0.0;
                        braktd = true;
                    }
                    tol = Math.abs(xmin) * reltol + abstol;
                    break ConvergenceCheck;
                }
            }

            /* If function value increased, origin remains unchanged */
            /* but new point may now qualify as w. */
            if (u < 0.0) {
                a = u;
            } else {
                b = u;
                braktd = true;
            }
            xw = u;
            fw = fu;
            gw = gu;
        }

        double twotol = tol + tol;
        double xmidpt = 0.5 * (a + b);

        /* Check termination criteria */
        boolean convrg = (Math.abs(xmidpt) <= twotol - 0.5 * (b - a))
                || (Math.abs(gmin) <= gtest2 && fmin < oldf
                && ((Math.abs(xmin - xbnd) > tol) || (!braktd)));
        if (convrg) {
            if (xmin != 0.0) {
                return false;
            }

            /*
             * If the function has not been reduced, check to see that the relative
             * change in f(x) is consistent with the estimate of the delta-
             * unimodality constant, tol. If the change in f(x) is larger than
             * expected, reduce the value of tol.
             */
            if (Math.abs(oldf - fw) <= fpresn) {
                throw new MinimizationError("getptciter failed. |oldf - fw| <= fpresn");
            }
            tol = 0.1 * tol;
            if (tol < tnytol) {
                throw new MinimizationError("getptciter failed. tol < tnytol");
            }
            reltol = 0.1 * reltol;
            abstol = 0.1 * abstol;
            twotol = 0.1 * twotol;
        }

        /* Continue with the computation of a trial step length */
        double r = 0.0;
        double q = 0.0;
        double s = 0.0;
        MinimumFound:
        {
            if (Math.abs(e) > tol) {
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
                if (b1 != step || braktd) {
                    e = step;
                }
            }
        }

        /* Construct an artificial bound on the estimated steplength */
        double a1 = a;
        b1 = b;
        step = xmidpt;
        if ((!braktd) || ((a == 0.0 && xw < 0.0) || (b == 0.0 && xw > 0.0))) {
            if (braktd) {
                /* If the minimum is not bracketed by 0 and xw the step must lie
                 within (a1,b1). */
                double d1 = xw;
                double d2 = a;
                if (a == 0.0) {
                    d2 = b;
                }
                /* This line might be : */
                /* if (*a == 0.0) d2 = *e */
                u = -d1 / d2;
                step = 5.0 * d2 * (0.1 + 1.0 / u) / 11.0;
                if (u < 1.0) {
                    step = 0.5 * d2 * Math.sqrt(u);
                }
            } else {
                step = -factor * xw;
                if (step > scxbnd) {
                    step = scxbnd;
                }
                if (step != scxbnd) {
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
            e = b - a;
        } else {
            /* A cubic interpolation step */
            step = s / q;

            /* The function must not be evaluated too close to a or b. */
            if (step - a < twotol || b - step < twotol) {
                if (xmidpt <= 0.0) {
                    step = -tol;
                } else {
                    step = tol;
                }
            }
        }

        /* If the step is too large, replace by the scaled bound (so as to */
        /* compute the new point on the boundary). */
        if (step >= scxbnd) {
            step = scxbnd;
            /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
            scxbnd -= (reltol * Math.abs(xbnd) + abstol) / (1.0 + reltol);
        }
        u = step;
        if (Math.abs(step) < tol && step < 0.0) {
            u = -tol;
        }
        if (Math.abs(step) < tol && step >= 0.0) {
            u = tol;
        }
        
        return true;
    }
    
    public double xmin() {
        return xmin;
    }
    
    public double fmin() {
        return fmin;
    }
    
    public double u() {
        return u;
    }
    
}
