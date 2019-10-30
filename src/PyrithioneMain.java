//this is an updated version of the deterioration rate diagnostic program, but with the poisson reseeding bug now fixed.
public class PyrithioneMain {
    public static void main(String[] args){

        //BioSystem.varyingDeteriorationAndThreshold(0.01);

        double scale_99 = 2.71760274, sigma_99 = 0.56002833;
        BioSystem.varyingTauStep(scale_99, sigma_99);
    }
}
