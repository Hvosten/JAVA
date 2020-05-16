package RandomizedMethods;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
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
import java.util.ListIterator;

//example of ode with weak singularity
public class WeakSingularityExample extends MyFirstOrderDifferentialEquations {
    private double[] gamma;
    private double T;
    private int _case;


    public WeakSingularityExample(double[] gamma, double T) {
        this.gamma = gamma;
        this.T=T;
        this._case = 0;
    }

    public int getDimension() {
        return 1;
    }


    public void computeDerivatives(double t, double[] y, double[] y_prime) throws MaxCountExceededException, DimensionMismatchException {
        y_prime[0] = Math.pow(T-t,-1.0/gamma[_case]);
    }

    public Double sol(double t){
        return (Math.pow(T,1-1/gamma[_case])- Math.pow(T-t,1-1/gamma[_case]))/(T-1/gamma[_case]);
    }

    public void increase_case(){
        ++_case;
    }

    public int get_case() {
        return _case;
    }

    public double getGamma(int index) {
        return gamma[index];
    }


    public double[] error_for_step(int step,int samples){
        double[] y = new double[] {0.0 }; // initial state
        this._case=0;
        double[] out = new double[this.gamma.length];
        while(this.get_case()<this.gamma.length){
            double pom = 0.0;
            FirstOrderIntegrator euler = new RandomEulerIntegrator(Math.pow(2.0,-step));

            for (int i = 0; i < samples; ++i) {
                y[0] = 0.0;

                euler.integrate(this, 0.0, y, this.T, y);
                pom += Math.pow(y[0] - this.sol(this.T), 2);
            }
            out[this.get_case()]=Math.sqrt(pom/samples);
            this.increase_case();
        }
        //System.out.println(step);
        //System.out.println(Arrays.toString(out));
        return out;
    }

    public void convergance_rate(int[] steps, int samples){
        XYSeries[] error = new XYSeries[5];
        final double clog = 1.0/Math.log(2.0);
        WeightedObservedPoints[] data = new WeightedObservedPoints[5];

        for(int i=0;i<error.length;++i) {
            error[i] = new XYSeries(String.format("gamma=%d", new Double(this.gamma[i]).intValue()));
            data[i] = new WeightedObservedPoints();
        }

        for(int h:steps) {
            double[] pom = error_for_step(h,samples);
            for(int i=0;i<pom.length;++i) {
                error[i].add(h, pom[i]);
                data[i].add(-h, clog * Math.log(pom[i]));
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
        yAxis.setRange(Math.pow(10.0,-4), Math.pow(10.0,1));
        yAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

        plot.setRangeAxis(yAxis);

        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(1);
        double[] coeff = fitter.fit(data[0].toList());
        System.out.println("Convergance is " + coeff[1]);
        final Marker target = new ValueMarker(1.0E-1);
        target.setPaint(Color.red);
        target.setLabel(String.format("%f",coeff[1]));
        target.setLabelAnchor(RectangleAnchor.TOP_RIGHT);
        target.setLabelTextAnchor(TextAnchor.BOTTOM_RIGHT);
        plot.addRangeMarker(target);

        coeff = fitter.fit(data[1].toList());
        System.out.println("Convergance is " + coeff[1]);
        final Marker target2 = new ValueMarker(1.0E-2);
        target2.setPaint(Color.red);
        target2.setLabel(String.format("%f",coeff[1]));
        target2.setLabelAnchor(RectangleAnchor.TOP_RIGHT);
        target2.setLabelTextAnchor(TextAnchor.BOTTOM_RIGHT);
        plot.addRangeMarker(target2);

        coeff = fitter.fit(data[2].toList());
        System.out.println("Convergance is " + coeff[1]);
        final Marker target3 = new ValueMarker(1.0E-2);
        target3.setPaint(Color.red);
        target3.setLabel(String.format("%f",coeff[1]));
        target3.setLabelAnchor(RectangleAnchor.CENTER);
        target3.setLabelTextAnchor(TextAnchor.BOTTOM_RIGHT);
        plot.addRangeMarker(target3);

        coeff = fitter.fit(data[3].toList());
        System.out.println("Convergance is " + coeff[1]);
        final Marker target4 = new ValueMarker(1.0E-3);
        target4.setPaint(Color.red);
        target4.setLabel(String.format("%f",coeff[1]));
        target4.setLabelAnchor(RectangleAnchor.TOP_RIGHT);
        target4.setLabelTextAnchor(TextAnchor.TOP_RIGHT);
        plot.addRangeMarker(target4);

        coeff = fitter.fit(data[4].toList());
        System.out.println("Convergance is " + coeff[1]);
        final Marker target5 = new ValueMarker(1.0E-3);
        target5.setPaint(Color.red);
        target5.setLabel(String.format("%f",coeff[1]));
        target5.setLabelAnchor(RectangleAnchor.CENTER);
        target5.setLabelTextAnchor(TextAnchor.BOTTOM_RIGHT);
        plot.addRangeMarker(target5);





        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer( );
        renderer.setSeriesPaint(0, Color.RED);

        renderer.setSeriesPaint(1, Color.BLUE);

        renderer.setSeriesPaint(2, Color.MAGENTA);
        renderer.setSeriesPaint(3, Color.GREEN);
        renderer.setSeriesPaint(4, Color.CYAN);

        plot.setRenderer(renderer);
        plot.setBackgroundPaint(Color.white);

        plot.setDomainGridlinePaint(Color.BLACK);
        plot.setRangeGridlinePaint(Color.BLACK);

        File XYChart = new File("C:\\Users\\Radek\\IdeaProjects\\ProjectK\\WeakSingularityError.jpg");
        try {
            ChartUtilities.saveChartAsJPEG(XYChart, xyLineChart, 640, 500);
        }
        catch(IOException e){
            System.out.println("Something went wrong!");
        }

    }

    public void plot_trajectory(double step){

        double[] y = new double[] {0.0 }; // initial state


        for(int i=0;i<gamma.length;++i){
            XYSeries numerical = new XYSeries("numerical");
            XYSeries exact = new XYSeries("exact");
            FirstOrderIntegrator euler = new RandomEulerIntegrator(step);
            MyStepHandler stepHandler = new MyStepHandler();
            euler.addStepHandler(stepHandler);
            euler.integrate(this,0.0,y,this.T,y);
            Iterator<Double> it_t = stepHandler.t.iterator();
            Iterator<Double[]> it_x = stepHandler.x.iterator();
            Double pom;
            while(it_t.hasNext() && it_x .hasNext()) {
                pom=it_t.next();
                numerical.add(pom, it_x.next()[0]);
                exact.add(pom,this.sol(pom));
            }
            XYSeriesCollection dataset = new XYSeriesCollection();
            dataset.addSeries(numerical);
            dataset.addSeries(exact);
            XYPlot plot = new XYPlot(
                    dataset,
                    new NumberAxis(),
                    new NumberAxis(),
                    new XYLineAndShapeRenderer(true, false));
            JFreeChart chart = new JFreeChart(String.format("gamma=%.0f", this.getGamma(i)), JFreeChart.DEFAULT_TITLE_FONT, plot, true);

            JFrame frame = new JFrame("Trajectory");
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setContentPane(new ChartPanel(chart));
            frame.pack();
            frame.setVisible(true);
            this.increase_case();
            y[0]=0;
        }


    }

    public static void main(String[] args) {
        WeakSingularityExample ode = new WeakSingularityExample(new double[]{2,3,5,8,10},1);
        //ode.convergance_rate(new int[]{2,3,4,5,6,7,8,9,10},1000);
        ode.plot_trajectory(0.1);
    }


}
