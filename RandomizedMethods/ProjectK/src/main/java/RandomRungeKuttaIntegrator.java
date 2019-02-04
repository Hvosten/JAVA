import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
//abstract class which represents random runge kutta methods
//we cannot use rungekuttaintegrator because in each iteration we randomly choose one runge kutta scheme
//so in each iteration we have to update the coefficients
public abstract class RandomRungeKuttaIntegrator extends AbstractIntegrator {
    private double[] c;
    private double[][] a;
    private double[] b;
    private RandomRungeKuttaStepInterpolator prototype;
    private double step;
    protected static double tau=Math.random();

    protected RandomRungeKuttaIntegrator(String name, double[] c, double[][] a, double[] b, RandomRungeKuttaStepInterpolator prototype, double step) {
        super(name);
        this.c = c;
        this.a = a;
        this.b = b;
        this.prototype = prototype;
        this.step = FastMath.abs(step);
    }

    public abstract double[][] update_a();
    public abstract double[] update_b();
    public abstract double[] update_c();
    public void update_coefficients(){
        UniformRealDistribution u = new UniformRealDistribution();
        tau=u.sample();
        this.c = update_c();
        this.a = update_a();
        this.b = update_b();
    }

    public void integrate(ExpandableStatefulODE equations, double t) throws NumberIsTooSmallException, DimensionMismatchException, MaxCountExceededException, NoBracketingException {
        this.sanityChecks(equations, t);
        this.setEquations(equations);
        boolean forward = t > equations.getTime();
        double[] y0 = equations.getCompleteState();
        double[] y = (double[])y0.clone();
        int stages = this.c.length + 1;
        double[][] yDotK = new double[stages][];

        for(int i = 0; i < stages; ++i) {
            yDotK[i] = new double[y0.length];
        }

        double[] yTmp = (double[])y0.clone();
        double[] yDotTmp = new double[y0.length];
        RandomRungeKuttaStepInterpolator interpolator = (RandomRungeKuttaStepInterpolator)this.prototype.copy();
        interpolator.reinitialize(this, yTmp, yDotK, forward, equations.getPrimaryMapper(), equations.getSecondaryMappers());
        interpolator.storeTime(equations.getTime());
        this.stepStart = equations.getTime();
        this.stepSize = forward ? this.step : -this.step;
        this.initIntegration(equations.getTime(), y0, t);
        this.isLastStep = false;


        do {
            this.update_coefficients();
            interpolator.shift();
            this.computeDerivatives(this.stepStart, y, yDotK[0]);
            for (int k = 1; k < stages; ++k) {

                for (int j = 0; j < y0.length; ++j) {
                    double sum = a[k-1][0] * yDotK[0][j];
                    for (int l = 1; l < k; ++l) {
                        sum += a[k-1][l] * yDotK[l][j];
                    }
                    yTmp[j] = y[j] + stepSize * sum;
                }

                this.computeDerivatives(stepStart + c[k-1] * stepSize, yTmp, yDotK[k]);
            }

            for(int j = 0; j < y0.length; ++j) {
                double sum = this.b[0] * yDotK[0][j];

                for(int l = 1; l < stages; ++l) {
                    sum += this.b[l] * yDotK[l][j];
                }

                yTmp[j] = y[j] + this.stepSize * sum;
            }

            interpolator.storeTime(this.stepStart + this.stepSize);
            System.arraycopy(yTmp, 0, y, 0, y0.length);
            System.arraycopy(yDotK[stages - 1], 0, yDotTmp, 0, y0.length);
            this.stepStart = this.acceptStep(interpolator, y, yDotTmp, t);
            if (!this.isLastStep) {
                interpolator.storeTime(this.stepStart);
                double nextT = this.stepStart + this.stepSize;
                boolean nextIsLast = forward ? nextT >= t : nextT <= t;
                if (nextIsLast) {
                    this.stepSize = t - this.stepStart;
                }
            }
        } while(!this.isLastStep);

        equations.setTime(this.stepStart);
        equations.setCompleteState(y);
        this.stepStart = 0.0D / 0.0;
        this.stepSize = 0.0D / 0.0;
    }
}