import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
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

//class implementing FirstOrderDifferentialEquations
//we dont have to now add stephandler every time we define ode
public abstract class MyFirstOrderDifferentialEquations implements FirstOrderDifferentialEquations {
    public abstract Double sol(double t);

    public class MyStepHandler implements org.apache.commons.math3.ode.sampling.StepHandler {
        public List<Double> t = new ArrayList<Double>();
        public List<Double[]> x = new ArrayList<Double[]>();


        public void init(double t0, double[] y0, double t) {
            this.t.add(t0);
            this.x.add(ArrayUtils.toObject(y0));
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

        public void plot_results(){
            XYSeries numerical_sol = new XYSeries("Numerical solution");
            XYSeries exact_sol = new XYSeries("Exact solution");
            ListIterator<Double> it_t = t.listIterator();
            ListIterator<Double[]> it_x = x.listIterator();
            Double pom;
            while(it_t.hasNext() && it_x .hasNext()) {
                pom=it_t.next();
                numerical_sol.add(pom, it_x.next()[0]);
                exact_sol.add(pom,sol(pom));
            }

            XYSeriesCollection dataset = new XYSeriesCollection();
            dataset.addSeries(numerical_sol);
            dataset.addSeries(exact_sol);

            JFreeChart xyLineChart = ChartFactory.createXYLineChart("","","",dataset, PlotOrientation.VERTICAL,true,true,false);
            File XYChart = new File("C:\\Users\\Radek\\IdeaProjects\\ProjectK\\Plot.jpg");
            try {
                ChartUtilities.saveChartAsJPEG(XYChart, xyLineChart, 640, 500);
            }
            catch(IOException e){
                System.out.println("Something went wrong!");
            }
        }
    }

}
