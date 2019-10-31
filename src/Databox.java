import java.lang.reflect.Array;

public class Databox {

    public double tau, simulation_time, simulation_time_stDev, exit_time, exit_time_stDev, det_rate, threshold_K;
    public double thickness, thickness_stDev;
    public double[] event_counters;


    public Databox(double tau, double simulation_time, double exit_time, double threshold_K, double det_rate, double thickness, double[] event_counters){
        //todo - alter the code so that thresholdK is consistently placed before det rate in the other methods
        this.tau = tau;
        this.simulation_time = simulation_time;
        this.simulation_time_stDev = 0.;
        this.exit_time = exit_time;
        this.exit_time_stDev = 0.;
        this.det_rate = det_rate;
        this.threshold_K = threshold_K;
        this.thickness = thickness;
        this.thickness_stDev = 0.;
        this.event_counters = event_counters;
    }


    public double getTau(){ return tau; }
    public double getSimulation_time(){return simulation_time;}
    public double getSimulation_time_stDev(){return simulation_time_stDev;}
    public void setSimulation_time_stDev(double simulation_time_stDev){this.simulation_time_stDev = simulation_time_stDev;}
    public double getExit_time(){return exit_time;}
    public double getExit_time_stDev(){return exit_time_stDev;}
    public void setExit_time_stDev(double exit_time_stDev){this.exit_time_stDev = exit_time_stDev;}
    public double getThreshold_K(){return threshold_K;}
    public double getDet_rate(){return det_rate;}
    public double getThickness(){ return thickness; }
    public double getThickness_stDev(){ return thickness_stDev; }
    public void setThickness_stDev(double thickness_stDev){this.thickness_stDev = thickness_stDev;}
    public double[] getEvent_counters(){ return event_counters; }


    public double[] allDataInAnArray(){
        //The order of this array correlates with the order of the counter labels. make sure they match.
        double[] non_counters;

        non_counters = new double[]{tau, simulation_time, simulation_time_stDev, exit_time, exit_time_stDev, threshold_K, det_rate, thickness, thickness_stDev};

        double[] all_vals = new double[non_counters.length + event_counters.length];

        System.arraycopy(non_counters, 0, all_vals, 0, non_counters.length);
        System.arraycopy(event_counters, 0, all_vals, non_counters.length, event_counters.length);

        return all_vals;
    }




    public static Databox add(Databox d1, Databox d2){
        double[] added_counters = new double[d1.event_counters.length];
        for(int i = 0; i < d1.event_counters.length; i++){
            added_counters[i] = d1.event_counters[i] + d2.event_counters[i];
        }

        return new Databox(d1.tau+d2.tau, d1.simulation_time+d2.simulation_time, d1.exit_time+d2.exit_time, d1.threshold_K+d2.threshold_K, d1.det_rate+d2.det_rate, d1.thickness+d2.thickness, added_counters);
    }





    public static Databox divideBy(Databox db, double divisor){
        double[] scaled_counters = new double[db.event_counters.length];
        for(int i = 0; i < db.event_counters.length; i++){
            scaled_counters[i] = db.event_counters[i]/divisor;
        }
        return new Databox(db.tau/divisor, db.simulation_time/divisor, db.exit_time/divisor, db.threshold_K/divisor,db.det_rate/divisor,  db.thickness/divisor, scaled_counters);
    }





    public static Databox averagedMeasurementsAndStDev(Databox[] databoxes){
        //averages the counters in an array of databoxes and returns the st dev of the recorded thicknesses
        Databox summedDB = databoxes[0];
        for(int i = 1; i < databoxes.length; i++){
            summedDB = Databox.add(summedDB, databoxes[i]);
        }
        Databox averagedDB = Databox.divideBy(summedDB, databoxes.length);

        double sumSq_thick = 0.;
        double avgThick = averagedDB.thickness;

        double sumSq_simTime = 0.;
        double avgSimTime = averagedDB.simulation_time;

        double sumSq_exitT = 0.;
        double avgExitT = averagedDB.exit_time;
        for(int i = 0; i < databoxes.length; i++){
            sumSq_thick += (databoxes[i].thickness - avgThick)*(databoxes[i].thickness - avgThick);
            sumSq_simTime += (databoxes[i].simulation_time - avgSimTime)*(databoxes[i].simulation_time - avgSimTime);
            sumSq_exitT += (databoxes[i].exit_time - avgExitT)*(databoxes[i].exit_time - avgExitT);
        }
        averagedDB.setThickness_stDev(Math.sqrt(sumSq_thick/(databoxes.length - 1.)));
        averagedDB.setSimulation_time_stDev(Math.sqrt(sumSq_simTime/(databoxes.length - 1.)));
        averagedDB.setExit_time_stDev(Math.sqrt(sumSq_exitT/(databoxes.length - 1.)));


        return averagedDB;
    }




}
