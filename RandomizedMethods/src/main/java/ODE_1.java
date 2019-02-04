import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.*;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.io.File;
import java.io.IOException;
import java.util.*;

//simple example how to use common math library to solve numerically ode
//in this example i used stephandler to plot trajectory
public class ODE_1 implements FirstOrderDifferentialEquations {

    public int getDimension() {
        return 1;
    }

    public void computeDerivatives(double t, double[] y, double[] y_prime) throws MaxCountExceededException, DimensionMismatchException {
       y_prime[0] = 2*(y[0]+1)/(t+1);
    }

    public static void main(String[] args) {
        //FirstOrderIntegrator euler = new EulerIntegrator(0.5);
        FirstOrderIntegrator euler = new ClassicalRungeKuttaIntegrator(1);
        //FirstOrderIntegrator euler = new RandomizedRungeKuttaIntegrator(0.2);
        FirstOrderDifferentialEquations ode = new ODE_1();
        double[] y = new double[] {1.0 }; // initial state


        MyStepHandler stepHandler = new MyStepHandler();
        euler.addStepHandler(stepHandler);
        euler.integrate(ode,0.0,y,10.0,y);

        System.out.println(String.format("The numerical solution at %f is equal %f",10.0,y[0]));
        stepHandler.show_results();
        stepHandler.plot_results(new exact_solution() {
            public Double func(Double t) {
                return 1+4*t+2*t*t;
            }
        });
    }

    public static class MyStepHandler implements org.apache.commons.math3.ode.sampling.StepHandler {
        List<Double> t = new ArrayList<Double>();
        List<Double[]> x = new ArrayList<Double[]>();



        public void init(double t0, double[] y0, double t) {
        }

        public void handleStep(StepInterpolator interpolator, boolean isLast) {
            t.add(interpolator.getCurrentTime());
            x.add(ArrayUtils.toObject(interpolator.getInterpolatedState()));
        }

        public void show_results(){
            Iterator<Double> it_t = t.iterator();
            Iterator<Double[]> it_x = x.iterator();
            while(it_t.hasNext() && it_x .hasNext())
            {
                System.out.println(String.format("%f\t%s",it_t.next(), Arrays.toString(it_x.next())));
            }
        }

        public void plot_results(exact_solution sol){
            XYSeries numerical_sol = new XYSeries("Numerical solution");
            XYSeries exact_sol = new XYSeries("Exact solution");
            ListIterator<Double> it_t = t.listIterator();
            ListIterator<Double[]> it_x = x.listIterator();
            Double pom;
            while(it_t.hasNext() && it_x .hasNext()) {
                pom=it_t.next();
                numerical_sol.add(pom, it_x.next()[0]);
                exact_sol.add(pom,sol.func(pom));
            }

            XYSeriesCollection dataset = new XYSeriesCollection();
            dataset.addSeries(numerical_sol);
            dataset.addSeries(exact_sol);

            JFreeChart xyLineChart = ChartFactory.createXYLineChart("Trajectory for h=1 (RungeKutta scheme)","","",dataset, PlotOrientation.VERTICAL,true,true,false);
            File XYChart = new File("C:\\Users\\Radek\\IdeaProjects\\ProjectK\\XYline_chart.jpg");
            try {
                ChartUtilities.saveChartAsJPEG(XYChart, xyLineChart, 640, 500);
            }
            catch(IOException e){
                System.out.println("Something went wrong!");
            }
        }
    }
}
