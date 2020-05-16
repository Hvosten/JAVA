package RandomizedMethods;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.EquationsMapper;
import org.apache.commons.math3.ode.sampling.AbstractStepInterpolator;

abstract class RandomRungeKuttaStepInterpolator extends AbstractStepInterpolator {
    protected double[] previousState;
    protected double[][] yDotK;
    protected AbstractIntegrator integrator;

    protected RandomRungeKuttaStepInterpolator() {
        this.previousState = null;
        this.yDotK = (double[][])null;
        this.integrator = null;
    }

    public RandomRungeKuttaStepInterpolator(RandomRungeKuttaStepInterpolator interpolator) {
        super(interpolator);
        if (interpolator.currentState != null) {
            this.previousState = (double[])interpolator.previousState.clone();
            this.yDotK = new double[interpolator.yDotK.length][];

            for(int k = 0; k < interpolator.yDotK.length; ++k) {
                this.yDotK[k] = (double[])interpolator.yDotK[k].clone();
            }
        } else {
            this.previousState = null;
            this.yDotK = (double[][])null;
        }

        this.integrator = null;
    }

    public void reinitialize(AbstractIntegrator rkIntegrator, double[] y, double[][] yDotArray, boolean forward, EquationsMapper primaryMapper, EquationsMapper[] secondaryMappers) {
        this.reinitialize(y, forward, primaryMapper, secondaryMappers);
        this.previousState = null;
        this.yDotK = yDotArray;
        this.integrator = rkIntegrator;
    }

    public void shift() {
        this.previousState = (double[])this.currentState.clone();
        super.shift();
    }

    public void writeExternal(ObjectOutput out) throws IOException {
        this.writeBaseExternal(out);
        int n = this.currentState == null ? -1 : this.currentState.length;

        int kMax;
        for(kMax = 0; kMax < n; ++kMax) {
            out.writeDouble(this.previousState[kMax]);
        }

        kMax = this.yDotK == null ? -1 : this.yDotK.length;
        out.writeInt(kMax);

        for(int k = 0; k < kMax; ++k) {
            for(int i = 0; i < n; ++i) {
                out.writeDouble(this.yDotK[k][i]);
            }
        }

    }

    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
        double t = this.readBaseExternal(in);
        int n = this.currentState == null ? -1 : this.currentState.length;
        int kMax;
        if (n < 0) {
            this.previousState = null;
        } else {
            this.previousState = new double[n];

            for(kMax = 0; kMax < n; ++kMax) {
                this.previousState[kMax] = in.readDouble();
            }
        }

        kMax = in.readInt();
        this.yDotK = kMax < 0 ? (double[][])null : new double[kMax][];

        for(int k = 0; k < kMax; ++k) {
            this.yDotK[k] = n < 0 ? null : new double[n];

            for(int i = 0; i < n; ++i) {
                this.yDotK[k][i] = in.readDouble();
            }
        }

        this.integrator = null;
        if (this.currentState != null) {
            this.setInterpolatedTime(t);
        } else {
            this.interpolatedTime = t;
        }

    }
}