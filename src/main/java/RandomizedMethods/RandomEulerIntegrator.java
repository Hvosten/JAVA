package RandomizedMethods;


import org.apache.commons.math3.ode.nonstiff.RungeKuttaIntegrator;

public class RandomEulerIntegrator extends RandomRungeKuttaIntegrator {
    private static double[] STATIC_C = new double[]{tau};
    private static double[][] STATIC_A = new double[][]{{0.0D}};
    private static double[] STATIC_B = new double[]{0.0D,1.0D};


    public RandomEulerIntegrator(double step) {
        super("RandomEuler", STATIC_C, STATIC_A, STATIC_B, new RandomEulerStepInterpolator(), step);
    }

    public double[][] update_a() {
        return STATIC_A;
    }

    public double[] update_b() {
        return STATIC_B;
    }

    public double[] update_c() {
        STATIC_C[0]=tau;
        return STATIC_C;
    }
}