import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

//simple example how to use common math library to solve numerically ode
//in this example i used stephandler to print numerical solution
public class ODE_2 implements FirstOrderDifferentialEquations {

    public int getDimension() {
        return 1;
    }

    public void computeDerivatives(double t, double[] y, double[] y_prime){
        y_prime[0] = 2*(y[0]+1)/(t+1);
    }

    public static void main(String[] args) {
        FirstOrderIntegrator euler = new EulerIntegrator(0.5);
        FirstOrderDifferentialEquations ode = new ODE_1();
        double[] y = new double[] {1.0 }; // initial state


        StepHandler stepHandler = new MyStepHandler();
        euler.addStepHandler(stepHandler);
        euler.integrate(ode,0.0,y,10.0,y);
    }

    public static class MyStepHandler implements org.apache.commons.math3.ode.sampling.StepHandler {

        public void init(double t0, double[] y0, double t) {
        }
        public void handleStep(StepInterpolator interpolator, boolean isLast) {
            double   t = interpolator.getCurrentTime();
            double[] y = interpolator.getInterpolatedState();
            System.out.println(String.format("%f\t%f",t,y[0]));
        }
    }
}