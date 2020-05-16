package RandomizedMethods;

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.util.Arrays;

//example how to fit line to points with common math in java
class MyFuncFitter{
    public static void main(String[] args) {
        WeightedObservedPoints data = new WeightedObservedPoints();
        data.add(2,4);
        data.add(3,6);
        data.add(5,10);

        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(1);
        double[] coeff = fitter.fit(data.toList());
        System.out.println(Arrays.toString(coeff));
    }
}