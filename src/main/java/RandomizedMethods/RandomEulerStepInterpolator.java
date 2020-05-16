package RandomizedMethods;

import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import java.util.Arrays;
import java.util.Random;

class RandomEulerStepInterpolator extends RandomRungeKuttaStepInterpolator {
    private static final long serialVersionUID = 20111120L;

    public RandomEulerStepInterpolator() {
    }

    public RandomEulerStepInterpolator(RandomEulerStepInterpolator interpolator) {
        super(interpolator);
    }

    protected StepInterpolator doCopy() {
        return new RandomEulerStepInterpolator(this);
    }

    protected void computeInterpolatedStateAndDerivatives(double theta, double oneMinusThetaH) {
        int i;
        if (this.previousState != null && theta <= 0.5D) {
            for(i = 0; i < this.interpolatedState.length; ++i) {
                this.interpolatedState[i] = this.previousState[i] + theta * this.h * this.yDotK[1][i];
            }

            System.arraycopy(this.yDotK[0], 0, this.interpolatedDerivatives, 0, this.interpolatedDerivatives.length);
        } else {
            for(i = 0; i < this.interpolatedState.length; ++i) {
                this.interpolatedState[i] = this.currentState[i] - oneMinusThetaH * this.yDotK[1][i];
            }

            System.arraycopy(this.yDotK[0], 0, this.interpolatedDerivatives, 0, this.interpolatedDerivatives.length);
        }

    }

}