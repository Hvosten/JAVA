package RandomizedMethods;

public class RandomizedRungeKuttaIntegrator extends RandomRungeKuttaIntegrator {
    private static double[] STATIC_C = new double[]{tau};
    private static double[][] STATIC_A = new double[][]{{tau}};
    private static double[] STATIC_B = new double[]{0.0D,1.0D};


    public RandomizedRungeKuttaIntegrator(double step) {
        super("RandomRungeKutta", STATIC_C, STATIC_A, STATIC_B, new RandomEulerStepInterpolator(), step);
    }

    public double[][] update_a() {
        STATIC_A[0][0]=tau;
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