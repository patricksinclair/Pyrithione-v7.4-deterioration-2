import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

public class Toolbox {


    private static int largestHeaderLength(String[] headers){
        int biggun = 0;
        for(String s : headers){
            if(s.length() > biggun) biggun = s.length();
        }
        return biggun;
    }



    static void writeDataboxArraylistToFile(String directoryName, String filename, String[] headers, ArrayList<Databox> databoxes){

        int string_length = Math.max(12, Toolbox.largestHeaderLength(headers)+3);
        headers[0] = "#"+headers[0];

        File directory = new File(directoryName);
        if(!directory.exists()) directory.mkdirs();

        File file = new File(directoryName+"/"+filename+".txt");

        try{

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            int ncols = headers.length;
            String file_header = "";
            for(int i = 0; i < headers.length-1; i++){
                String heado = headers[i]+",";
                file_header += String.format("%-"+string_length+"s", heado);
            }
            String heado = headers[headers.length-1];
            file_header += String.format("%-"+string_length+"s", heado);
            bw.write(file_header);
            //bw.newLine();


            for(int i = 0; i < databoxes.size(); i++){
                bw.newLine();

                String output = "";
                double[] outputVals = databoxes.get(i).allDataInAnArray();
                for(int v = 0; v < outputVals.length-1; v++){
                    String num_val = String.format("%.4E", outputVals[v])+",";
                    output += String.format("%-"+string_length+"s", num_val);
                }
                String num_val = String.format("%.4E", outputVals[outputVals.length-1]);
                output += String.format("%-"+string_length+"s", num_val);

                bw.write(output);

            }


            bw.close();

        }catch (IOException e){}
    }






    public static String millisToShortDHMS(long duration) {
        String res = "";
        long days  = TimeUnit.MILLISECONDS.toDays(duration);
        long hours = TimeUnit.MILLISECONDS.toHours(duration)
                - TimeUnit.DAYS.toHours(TimeUnit.MILLISECONDS.toDays(duration));
        long minutes = TimeUnit.MILLISECONDS.toMinutes(duration)
                - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(duration));
        long seconds = TimeUnit.MILLISECONDS.toSeconds(duration)
                - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(duration));
        if (days == 0) {
            res = String.format("%02d:%02d:%02d", hours, minutes, seconds);
        }
        else {
            res = String.format("%dd%02d:%02d:%02d", days, hours, minutes, seconds);
        }
        return res;
    }
}
