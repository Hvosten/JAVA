import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.util.LogFormat;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.TextAnchor;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;

//example of ode with jumps
public class ODEwithJumps extends MyFirstOrderDifferentialEquations {
    private double T;

    public ODEwithJumps(double t) {
        T = t;
    }

    public int getDimension() {
        return 1;
    }

    public void computeDerivatives(double t, double[] y, double[] y_prime) throws MaxCountExceededException, DimensionMismatchException {
        y_prime[0]=g(t)*y[0];
    }

    public double g(double t){
        return -1.0/10*Math.signum(T/4-t)-1.0/5*Math.signum(T/2-t)-7.0/10*Math.signum(3.0*T/4-t);
    }

    public double Heav(double t){
        return 0.5*(1+Math.signum(t));
    }

    public Double sol(double t) {
        double pom=-(13.0/10)+t-1.0/20*(-1+4*t)*Heav(1.0/4-t)*Heav(1-t)-1.0/5*(-1+2*t)*Heav(1.0/2-t)*Heav(1-t)+
                21.0/20*Heav(3.0/4-t)*Heav(1-t)-
                7.0/5*t*Heav(3.0/4-t)*Heav(1-t);
        return Math.exp(pom);
    }

    public double[] error_for_step(int step,int samples){
        double[] y = new double[] {1.0 }; // initial state

        double[] out = new double[3];
        double pom = 0.0;

        FirstOrderIntegrator euler = new RandomEulerIntegrator(Math.pow(2.0,-step));
        for (int i = 0; i < samples; ++i) {
            y[0] = 1.0;
            euler.integrate(this, 0.0, y, this.T, y);
            pom += Math.pow(y[0] - this.sol(this.T), 2);
        }
        out[0]=Math.sqrt(pom/samples);

        pom = 0.0;
        FirstOrderIntegrator runge = new RandomizedRungeKuttaIntegrator(Math.pow(2.0,-step));
        for (int i = 0; i < samples; ++i) {
            y[0] = 1.0;
            runge.integrate(this, 0.0, y, this.T, y);
            pom += Math.pow(y[0] - this.sol(this.T), 2);
        }
        out[1]=Math.sqrt(pom/samples);

        pom = 0.0;
        FirstOrderIntegrator classEuler = new EulerIntegrator(Math.pow(2.0,-step));
        for (int i = 0; i < samples; ++i) {
            y[0] = 1.0;
            classEuler.integrate(this, 0.0, y, this.T, y);
            pom += Math.pow(y[0] - this.sol(this.T), 2);
        }
        out[2]=Math.sqrt(pom/samples);
        //System.out.println(Arrays.toString(out));
        return out;
    }

    public void convergance_rate(int[] steps, int samples){

        XYSeries[] error = new XYSeries[3];

        error[0] = new XYSeries("Randomized Euler");
        error[1] = new XYSeries("Randomized Runge-Kutta");
        error[2] = new XYSeries("Classical Euler");


        final double clog = 1.0/Math.log(2.0);
        WeightedObservedPoints[] data = new WeightedObservedPoints[3];
        data[0]=new WeightedObservedPoints();
        data[1]=new WeightedObservedPoints();
        data[2]=new WeightedObservedPoints();

        for(int h:steps) {
            double[] pom = error_for_step(h,samples);
            for(int i=0;i<pom.length;++i) {
                error[i].add(h, pom[i]);
                data[i].add(-h,clog*Math.log(pom[i]));
            }

        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        for(int i=0;i<error.length;++i)
            dataset.addSeries(error[i]);

        JFreeChart xyLineChart = ChartFactory.createXYLineChart("L^2 convergence","time","error",dataset, PlotOrientation.VERTICAL,true,true,false);

        XYPlot plot = xyLineChart.getXYPlot();

        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
        domain.setRange(steps[0], steps[steps.length-1]);
        domain.setStandardTickUnits(NumberAxis.createIntegerTickUnits());


        LogAxis yAxis = new LogAxis();
        yAxis.setBase(10);
        LogFormat format = new LogFormat(yAxis.getBase(), "10", "^", true);
        yAxis.setNumberFormatOverride(format);
        yAxis.setLabel("error");
        yAxis.setRange(Math.pow(10.0,-6), Math.pow(10.0,0));
        yAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

        plot.setRangeAxis(yAxis);


        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(1);

        double[] coeff = fitter.fit(data[0].toList());
        System.out.println("Convergance for RandomEuler is "+coeff[1]);
        final Marker target = new ValueMarker(1.0E-3);
        target.setPaint(Color.red);
        target.setLabel(String.format("%f",coeff[1]));
        target.setLabelAnchor(RectangleAnchor.CENTER);
        target.setLabelTextAnchor(TextAnchor.BOTTOM_RIGHT);
        plot.addRangeMarker(target);

        coeff = fitter.fit(data[1].toList());
        System.out.println("Convergance for RandomRungeKutta is "+coeff[1]);
        final Marker target2 = new ValueMarker(1.0E-5);
        target2.setPaint(Color.red);
        target2.setLabel(String.format("%f",coeff[1]));
        target2.setLabelAnchor(RectangleAnchor.TOP_RIGHT);
        target2.setLabelTextAnchor(TextAnchor.BOTTOM_RIGHT);
        plot.addRangeMarker(target2);

        coeff = fitter.fit(data[2].toList());
        System.out.println("Convergance for ClassicalEuler is "+coeff[1]);
        final Marker target3 = new ValueMarker(1.0E-1);
        target3.setPaint(Color.red);
        target3.setLabel(String.format("%f",coeff[1]));
        target3.setLabelAnchor(RectangleAnchor.TOP_LEFT);
        target3.setLabelTextAnchor(TextAnchor.BOTTOM_LEFT);
        plot.addRangeMarker(target3);


        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
        renderer.setSeriesPaint(0, Color.RED);

        renderer.setSeriesPaint(1, Color.BLUE);

        renderer.setSeriesPaint(2, Color.MAGENTA);

        plot.setRenderer(renderer);
        plot.setBackgroundPaint(Color.white);

        plot.setDomainGridlinePaint(Color.BLACK);
        plot.setRangeGridlinePaint(Color.BLACK);





        File XYChart = new File("C:\\Users\\Radek\\IdeaProjects\\ProjectK\\ODEwithJumpsError.jpg");
        try {
            ChartUtilities.saveChartAsJPEG(XYChart, xyLineChart, 640, 500);
        }
        catch(IOException e){
            System.out.println("Something went wrong!");
        }

    }

    public void plot_trajectory(double step) {
        double[] y = new double[]{1.0}; // initial state

        //XYSeries numerical = new XYSeries("numerical");
        //XYSeries exact = new XYSeries("exact");
        FirstOrderIntegrator euler = new RandomizedRungeKuttaIntegrator(step);
        MyStepHandler stepHandler = new MyStepHandler();
        euler.addStepHandler(stepHandler);
        euler.integrate(this, 0.0, y, this.T, y);
        stepHandler.plot_results();
        /*
        Iterator<Double> it_t = stepHandler.t.iterator();
        Iterator<Double[]> it_x = stepHandler.x.iterator();
        Double pom;
        while (it_t.hasNext() && it_x.hasNext()) {
            pom = it_t.next();
            numerical.add(pom, it_x.next()[0]);
            exact.add(pom, this.sol(pom));
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(numerical);
        dataset.addSeries(exact);
        XYPlot plot = new XYPlot(
                dataset,
                new NumberAxis(),
                new NumberAxis(),
                new XYLineAndShapeRenderer(true, false));
        JFreeChart chart = new JFreeChart("", JFreeChart.DEFAULT_TITLE_FONT, plot, true);

        JFrame frame = new JFrame("Trajectory");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setContentPane(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);*/
    }

    public static void main(String[] args) {
        ODEwithJumps ode = new ODEwithJumps(1.0);
        //ode.convergance_rate(new int[]{2,3,4,5,6,7,8,9,10},1000);
        ode.plot_trajectory(0.08);
    }
}
