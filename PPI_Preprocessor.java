package ppi_preprocessor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author The Boss
 */
public class PPI_Preprocessor {
    
    private int numPPIs = 3000;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            PPI_Preprocessor pp = new PPI_Preprocessor("C:\\Users\\The Boss\\Desktop\\School\\CSI 5126\\Final Project\\9606.protein.links.v10.5.txt",
                "C:\\Users\\The Boss\\Desktop\\School\\CSI 5126\\Final Project\\9606.protein.sequences.v10.5.fa",
                "C:\\Users\\The Boss\\Desktop\\School\\CSI 5126\\Final Project\\preprocessed.fa");
        } catch (IOException ex) {
            Logger.getLogger(PPI_Preprocessor.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private final String interactionFileName, seqFileName, newFileName;
    
    public PPI_Preprocessor(String interactionFileName, String seqFileName, String newFileName) throws IOException {
        this.interactionFileName = interactionFileName;
        this.seqFileName = seqFileName;
        this.newFileName = newFileName;
        preprocessDB();
    }
    
    private void preprocessDB() throws FileNotFoundException, IOException {
        Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(newFileName), "utf-8"));
        BufferedReader br = new BufferedReader(new FileReader(interactionFileName));
        int[] randLines = genRandomInts(numPPIs, getNumLines(interactionFileName));
        int lineNum = 0, i = 0;
        try {
            StringBuilder sb ;
            String line = br.readLine();
            
            while (line != null && i < randLines.length) {
                //System.out.println(lineNum + " : " + randLines[i] + " " + (lineNum == randLines[i]));
                if (lineNum == randLines[i]) {
                    //System.out.println(line);
                    sb = new StringBuilder();
                    // write sequence 1
                    sb.append(">");
                    sb.append(line.split(" ")[0]);
                    sb.append(System.lineSeparator());
                    sb.append(getSequence(seqFileName, line.split(" ")[0]));
                    sb.append(System.lineSeparator());
                    //System.out.println(sb.toString() + "\n\n");
                    writer.write(sb.toString());
                    
                    sb = new StringBuilder();
                    // write sequence 2
                    sb.append(">");
                    sb.append(line.split(" ")[1]);
                    sb.append(System.lineSeparator());
                    sb.append(getSequence(seqFileName, line.split(" ")[1]));
                    sb.append(System.lineSeparator());
                    //System.out.println(sb.toString() + "\n\n");
                    writer.write(sb.toString());
                    
                    i++;
                }
                lineNum ++;
                line = br.readLine();
            }
            
        } finally {
            br.close();
            writer.close();
        }
    
    }
    
    private int getNumLines(String fileName) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        int lines = 0;
        String line;
        while ( (line = br.readLine()) != null) {
            lines ++;
        }
        System.out.println("lines = " + lines);
        return lines;
    }
    
    private String getSequence(String fileName, String seqName) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        StringBuilder sb = new StringBuilder();
        //System.out.println("Seq: " + seqName);
        try {
            String line = br.readLine();
            while (line != null && !line.equals(">"+seqName)) {
                line = br.readLine();
            }
            
            line = br.readLine();
            while (line != null && !line.startsWith(">")) {
                //System.out.println(line);
                sb.append(line);
                line = br.readLine();
            }
            //System.out.println(sb.toString() + "\n");
            
        } finally {
            br.close();
        }
        
        //System.out.println(line);
        return sb.toString();
    }
    
    private int[] genRandomInts(int numValues, int max) {
        int[] randNums = new int[numValues];
        Random r = new Random();
        for (int i = 0; i < numValues; i++) {
            int randInt = r.nextInt(max);
            while (arrayContains(randNums, randInt)) {
                randInt = r.nextInt(max);
            }
            randNums[i] = randInt;
        }
        Arrays.sort(randNums);
        return randNums;
    }
    
    public boolean arrayContains(int[] array, int i) {
        for (int e : array){
            if (e == i) {
                return true;
            }
        }
        return false;
    }
    
}
