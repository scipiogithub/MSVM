/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package multi_svm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
//import jdk.management.resource.internal.inst.FileOutputStreamRMHooks;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 *
 * @author scipio
 */
public class Multi_SVM {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
        // TODO code application logic here
        String file[] = {/*"gsa", "usps", "Smartphone"/*,"pancan", "pancan", "pancan"*/}, nclass[] = {/*"6", "15", "20", "96", "50", "100", "200"/*, "380"*/};
        int features[] = {4, 8, 16};
        //printwithclass_all(nclass, file, features, 1, 100, 10, (int) (100000 / 150.0));
        int m[] = {100000};
        double wzs[] = {0.51};
        double ci[] = {20.5};
        int pop[] = {25};
        double cooli[] = {0.9999};
        boolean parallel = true;
        //localwithclass("50", "Smartphone", 92, 20.5, 100000, 4, 0.51, 0.9999, 1);
        for (int mind = 0; mind < 1; mind++) {
            int mi = m[0], p = pop[0];
            double wz = wzs[0];
            double c = ci[0];
            double cooling = cooli[0];
            System.out.println("iter=" + mi + "\tr=" + wz + "\tc=" + c + "\tcooling=" + cooling + "\t#Foodsources=" + p);
            System.out.println("Dataset\tnclass\t#Features\tAcc\tTime(s)\tCost Func");
            for (int f = 0; f < file.length; f++) {
                double res[][] = printwithclass(nclass, file[f], features[f], c, mi, p, wz, cooling, parallel);

                for (int i = 0; i < res.length; i++) {
                    if (res[i][0] == 0) {
                        continue;
                    }
                    System.out.println(file[f] + "\t" + nclass[i] + "\t" + Eps.n + "\t" + String.format("%.2f", res[i][0]) + "\u00B1" + String.format("%.2f", res[i][1])
                            + "\t" + String.format("%.2f", res[i][2]) + "\u00B1" + String.format("%.2f", res[i][3]) + "\t"
                            + String.format("%.2f", res[i][4]) + "\u00B1" + String.format("%.2f", res[i][5]) + "\t"
                            + String.format("%.2f", res[i][6]) + "\u00B1" + String.format("%.2f", res[i][7]) + "\t"
                            + String.format("%.2f", res[i][8]) + "\u00B1" + String.format("%.2f", res[i][9]) + "\t"
                            + String.format("%.2f", res[i][10]) + "\u00B1" + String.format("%.2f", res[i][11]) + "\t"
                            + String.format("%.2f", res[i][12]) + "\u00B1" + String.format("%.2f", res[i][13]) + "\t"
                            + String.format("%.2f", res[i][14]) + "\u00B1" + String.format("%.2f", res[i][15]) + "\t"
                            + String.format("%.2f", res[i][16]) + "\u00B1" + String.format("%.2f", res[i][17]));
                }
            }
            String file2[] = {"pems",  "Synthetic_500_150000_5", "Synthetic_500_200000_7", "Synthetic_500_250000_10","mll", "cll_sub_111"};//{/*"pems", "pems","pems", "pems", "pems", "lung", "opt","mov", "sem"/*, "mfe", "cna","K3B","sector", "Synthetic_500_150000_5", "Synthetic_500_200000_7", "Synthetic_500_250000_10",*/"mll"/* "cll_sub_111"*/};
            int features2[] = {139, 150, 200, 250, 5, 10};//{/*139, 693,1387, 2774, 6934, 5, 23, 24, 38, 16, 17, 30, 3450, 15, 20, 25, 80,*/ 2};
            for (int f = 0; f < file2.length; f++) {
                double res[] = printwithoutclass(file2[f], features2[f], c, mi, p, wz, cooling, parallel);
                //System.out.println("Dataset\tnclass\t#Features\tAcc\tTime(s)\tEmp. Bee Contribution to w\t....to z\t....to b\tOnl. Bee Contribution to w\t....to z\t....to b\tCost Func");
                if (res[0] == 0) {
                    continue;
                }
                System.out.println(file2[f] + "\t\t" + Eps.n + "\t" + String.format("%.2f", res[0]) + "\u00B1" + String.format("%.2f", res[1])
                        /*+ "\t" + String.format("%.2f", res[2]) + "\u00B1" + String.format("%.2f", res[3]) + "\t"
                        + String.format("%.2f", res[4]) + "\u00B1" + String.format("%.2f", res[5]) + "\t"
                        + String.format("%.2f", res[6]) + "\u00B1" + String.format("%.2f", res[7]) + "\t"
                        + String.format("%.2f", res[8]) + "\u00B1" + String.format("%.2f", res[9]) + "\t"
                        + String.format("%.2f", res[10]) + "\u00B1" + String.format("%.2f", res[11]) + "\t"
                        + String.format("%.2f", res[12]) + "\u00B1" + String.format("%.2f", res[13]) + "\t"
                        + String.format("%.2f", res[14]) + "\u00B1" + String.format("%.2f", res[15]) + "\t"
                        + String.format("%.2f", res[16]) + "\u00B1"*/ + String.format("%.2f", res[17]));
            }
        }

    }

    public static void readExcel(File train_file, File test_file, boolean csv,
            int myy[][], int myy_test[][], double myx[][][], double myx_test[][][]) throws FileNotFoundException, IOException {
        FileInputStream fis = new FileInputStream(train_file);
        FileInputStream fis_test = test_file == null ? null : new FileInputStream(test_file);
//creating workbook instance that refers to .xls file     //obtaining bytes from the file  
//creating Workbook instance that refers to .xlsx file  
        ArrayList<ArrayList<Double>> data = new ArrayList<>(), data_test = new ArrayList<>();
        if (!csv) {
            XSSFWorkbook wb = new XSSFWorkbook(fis), wb_test = fis_test == null ? null : new XSSFWorkbook(fis_test);
            XSSFSheet sheet = wb.getSheetAt(0), sheet_test = wb_test == null ? null : wb_test.getSheetAt(0);     //creating a Sheet object to retrieve object  
            Iterator<Row> itr = sheet.iterator(), itr_test = sheet_test == null ? null : sheet_test.iterator();    //iterating over excel file  
            while (itr.hasNext()) {
                if (itr_test != null && itr_test.hasNext()) {
                    Row row = itr_test.next();
                    Iterator<Cell> cellIterator = row.cellIterator();   //iterating over each column  
                    ArrayList<Double> datarow = new ArrayList<>();
                    while (cellIterator.hasNext()) {
                        Cell cell = cellIterator.next();
                        datarow.add(cell.getNumericCellValue());
                    }
                    data_test.add(datarow);
                }
                Row row = itr.next();
                Iterator<Cell> cellIterator = row.cellIterator();   //iterating over each column  
                ArrayList<Double> datarow = new ArrayList<>();
                while (cellIterator.hasNext()) {
                    Cell cell = cellIterator.next();
                    datarow.add(cell.getNumericCellValue());
                }
                data.add(datarow);
            }
        } else if (train_file.getAbsolutePath().contains("sector")) {
            BufferedReader train_reader = new BufferedReader(new InputStreamReader(fis)),
                    test_reader = fis_test == null ? null : new BufferedReader(new InputStreamReader(fis_test));
            String train_line = train_reader.readLine(), test_line = test_reader == null ? null : test_reader.readLine();
            int train_row = 3847, features_num = 55197, test_row = 2565, row = 0;
            int y[] = new int[train_row], y_test[] = new int[test_row];
            double x[][] = new double[train_row][features_num], x_test[][] = new double[test_row][features_num];
            while (train_line != null) {
                
                String cols[] = train_line.split(",");
                for (int i = 0; i < cols.length - 1; i++) {
                    x[row][i] = Double.parseDouble(cols[i]);
                }
                y[row] = new BigDecimal(cols[features_num]).intValue() - 1;
                train_line = train_reader.readLine();
                row += 1;
            }
            row = 0;
            while (test_line != null) {
                String cols[] = test_line.split(",");
                for (int i = 0; i < cols.length - 1; i++) {
                    x_test[row][i] = Double.parseDouble(cols[i]);
                }
                y_test[row] = new BigDecimal(cols[features_num]).intValue() - 1;
                test_line = test_reader.readLine();
                row += 1;
            }
            myy[0] = y;
            myy_test[0] = y_test;
            myx[0] = x;
            myx_test[0] = x_test;
            return;
        } else {
            BufferedReader train_reader = new BufferedReader(new InputStreamReader(fis)),
                    test_reader = fis_test == null ? null : new BufferedReader(new InputStreamReader(fis_test));
            String train_line = train_reader.readLine(), test_line = test_reader == null ? null : test_reader.readLine();
            while (train_line != null) {
                ArrayList<Double> datarow = new ArrayList<>();
                String cols[] = train_line.split(",");
                for (int i = 0; i < cols.length; i++) {
                    datarow.add(Double.parseDouble(cols[i]));
                }
                data.add(datarow);
                train_line = train_reader.readLine();
                if (test_line != null) {
                    datarow = new ArrayList<>();
                    cols = train_line.split(",");
                    for (int i = 0; i < cols.length; i++) {
                        datarow.add(Double.parseDouble(cols[i]));
                    }
                    data_test.add(datarow);
                    test_line = test_reader.readLine();
                }
            }
        }
        int y[] = new int[data.size()], y_test[] = data_test.isEmpty() ? y : new int[data_test.size()];
        double x[][] = new double[data.size()][data.get(0).size() - 1], x_test[][] = data_test.isEmpty() ? x : new double[data_test.size()][data_test.get(0).size() - 1];
        for (int i = 0; i < data.size(); i++) {
            ArrayList<Double> row = data.get(i), row_test = null;
            if (i < data_test.size()) {
                row_test = data_test.get(i);
                y_test[i] = row_test.get(row_test.size() - 1).intValue() - 1;
            }
            y[i] = row.get(row.size() - 1).intValue() - 1;
            for (int j = 0; j < row.size() - 1; j++) {
                x[i][j] = row.get(j);
                if (row_test != null) {
                    x_test[i][j] = row_test.get(j);
                }
            }
        }
        myy[0] = y;
        myy_test[0] = y_test;
        myx[0] = x;
        myx_test[0] = x_test;

    }

    public static void localwithclass(String nclass, String file, int features, double C, int max_iter, int numberoffources, double wz, double cooling, int f) throws IOException, InterruptedException {
        String sample_train = file + "_train" + nclass + "_" + f + (file.equals("pancan") ? ".csv" : ".xlsx"),
                sample_test = file + "_test" + nclass + "_" + f + (file.equals("pancan") ? ".csv" : ".xlsx");
        File train_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_train),
                test_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_test),
                ls_Analiyze = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + file + "_" + f + "_ls.csv");
        FileWriter out = new FileWriter(ls_Analiyze);
        BufferedWriter fos = new BufferedWriter(out);
        if (!train_file.exists() || !test_file.exists()) {
            return;
        }
        int y[][] = new int[1][], y_test[][] = new int[1][];
        y[0] = null;
        y_test[0] = null;
        double x[][][] = new double[1][][], x_test[][][] = new double[1][][];
        x[0] = null;
        x_test[0] = null;
        readExcel(train_file, test_file, file.equals("pancan"), y, y_test, x, x_test);
        Eps.setVars(x[0], y[0], C, features);

        double tbw[] = new double[1];

        FoodSource_LS foodsources[] = new FoodSource_LS[numberoffources];
        boolean fs_left[] = new boolean[numberoffources];
        //System.out.println("pahse 1");
        long start = System.currentTimeMillis();
        for (int i = 0; i < numberoffources; i++) {
            fs_left[i] = false;
            foodsources[i] = new FoodSource_LS(foodsources, tbw, max_iter, 10, fs_left, wz, cooling);
            tbw[0] += 1.0 / foodsources[i].fitness;
        }
        //System.out.println("Phase 2");
        ExecutorService es = Executors.newCachedThreadPool();
        for (FoodSource_LS i : foodsources) {
            es.execute(i);
        }
        es.shutdown();
        es.awaitTermination(1, TimeUnit.DAYS);
        for (int i = 0; i < max_iter; i++) {
            for (int fs = 0; fs < foodsources.length; fs++) {
                fos.write((foodsources[fs].fitnesses.get(i) + ";").replace('.', ','));
            }
            fos.newLine();
        }
        fos.close();
    }

    public static double[][] printwithclass(String nclass[], String file, int features, double C, int max_iter,
            int numberoffources, double wz, double cooling, boolean parallel) throws FileNotFoundException, IOException, InterruptedException {
        double result[][] = new double[nclass.length][18];
        for (int n = 0; n < nclass.length; n++) {
            ArrayList<Double> accuracy = new ArrayList(), times = new ArrayList<>(), obj = new ArrayList<>();
            ArrayList<Integer> empw = new ArrayList<>(), empz = new ArrayList<>(), empb = new ArrayList<>(),
                    onlw = new ArrayList<>(), onlz = new ArrayList<>(), onlb = new ArrayList<>();
            for (int f = 1; f <= 10; f++) {
                String sample_train = file + "_train" + nclass[n] + "_" + f + (file.equals("pancan") ? ".csv" : ".xlsx"),
                        sample_test = file + "_test" + nclass[n] + "_" + f + (file.equals("pancan") ? ".csv" : ".xlsx");
                File train_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_train),
                        test_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_test);
                if (!train_file.exists() || !test_file.exists()) {
                    continue;
                }
                int y[][] = new int[1][], y_test[][] = new int[1][];
                y[0] = null;
                y_test[0] = null;
                double x[][][] = new double[1][][], x_test[][][] = new double[1][][];
                x[0] = null;
                x_test[0] = null;
                readExcel(train_file, test_file, file.equals("pancan"), y, y_test, x, x_test);
                Eps.setVars(x[0], y[0], C, features);

                double tbw[] = new double[1];
                //int numberoffources = 10;
                //int max_iter = 100000;
                //int limit = (int) (max_iter / 150.0);

                //System.out.println("Phase 2");
                if (parallel) {
                    FoodSource foodsources[] = new FoodSource[numberoffources];
                    boolean fs_left[] = new boolean[numberoffources];
                    //System.out.println("pahse 1");
                    long start = System.currentTimeMillis();
                    for (int i = 0; i < numberoffources; i++) {
                        fs_left[i] = false;
                        foodsources[i] = new FoodSource(foodsources, tbw, max_iter, 10, fs_left, wz, cooling);
                        tbw[0] += 1.0 / foodsources[i].fitness;
                    }
                    ExecutorService es = Executors.newCachedThreadPool();
                    for (FoodSource i : foodsources) {
                        es.execute(i);
                    }

                    es.shutdown();
                    es.awaitTermination(1, TimeUnit.DAYS);
                    start = System.currentTimeMillis() - start;
                    int i = 0;
                    //System.out.println("Conpleted");
                    foodsources[i].updatewz();
                    //foodsources[i].abandon();
                    for (int j = 1; j < numberoffources; j++) {
                        foodsources[j].updatewz();
                        //foodsources[j].abandon();
                        if (foodsources[j].fitness < foodsources[i].fitness) {
                            i = j;
                        }
                    }
                    times.add(start / 1000.0);
                    empw.add(foodsources[i].employebeeupdatew);
                    empz.add(foodsources[i].employebeeupdatez);
                    empb.add(foodsources[i].employebeeupdateb);
                    onlw.add(foodsources[i].onlookerbeeupdatew);
                    onlz.add(foodsources[i].onlookerbeeupdatez);
                    onlb.add(foodsources[i].onlookerbeeupdateb);
                    obj.add(foodsources[i].best);
                    int y_pred[] = Eps.predict(foodsources[i].bestw, foodsources[i].bestz, x_test[0]);
                    accuracy.add(Eps.accuracy_rate(y_test[0], y_pred) * 100.0);
                } else {
                    FoodSource_SQ foodsources[] = new FoodSource_SQ[numberoffources];
                    boolean fs_left[] = new boolean[numberoffources];
                    //System.out.println("pahse 1");
                    long start = System.currentTimeMillis();
                    for (int i = 0; i < numberoffources; i++) {
                        fs_left[i] = false;
                        foodsources[i] = new FoodSource_SQ(foodsources, tbw, max_iter, 10, fs_left, wz, cooling);
                        tbw[0] += 1.0 / foodsources[i].fitness;
                    }
                    for (int counter = 0; counter < max_iter; counter++) {
                        for (FoodSource_SQ i : foodsources) {
                            i.run();
                        }
                    }
                    start = System.currentTimeMillis() - start;
                    int i = 0;
                    //System.out.println("Conpleted");
                    foodsources[i].updatewz();
                    //foodsources[i].abandon();
                    for (int j = 1; j < numberoffources; j++) {
                        foodsources[j].updatewz();
                        //foodsources[j].abandon();
                        if (foodsources[j].fitness < foodsources[i].fitness) {
                            i = j;
                        }
                    }
                    times.add(start / 1000.0);
                    empw.add(foodsources[i].employebeeupdatew);
                    empz.add(foodsources[i].employebeeupdatez);
                    empb.add(foodsources[i].employebeeupdateb);
                    onlw.add(foodsources[i].onlookerbeeupdatew);
                    onlz.add(foodsources[i].onlookerbeeupdatez);
                    onlb.add(foodsources[i].onlookerbeeupdateb);
                    obj.add(foodsources[i].best);
                    int y_pred[] = Eps.predict(foodsources[i].bestw, foodsources[i].bestz, x_test[0]);
                    accuracy.add(Eps.accuracy_rate(y_test[0], y_pred) * 100.0);
                }

                //System.out.print(f + " ");
                //System.out.println(file + " " + foodsources[i].fitness + " " + Eps.cost(foodsources[i].w, foodsources[i].z));
                //System.out.println(Eps.accuracy_rate(y_test, y_pred));
                //System.out.println(start / 1000.0 + " " + foodsources[i].r);
            }
            if (accuracy.isEmpty()) {
                continue;
            }
            double mean_time = 0, std_time = 0, mean_accuracy = 0, std_accuracy = 0, mean_empw = 0, mean_empz = 0, mean_empb = 0,
                    mean_onlw = 0, mean_onlz = 0, mean_onlb = 0, std_empw = 0, std_empb = 0, std_onlw = 0, std_empz = 0, std_onlz = 0, std_onlb = 0,
                    mean_obj = 0, std_obj = 0;
            for (int i = 0; i < times.size(); i++) {
                mean_accuracy += accuracy.get(i);
                mean_time += times.get(i);
                mean_empw += empw.get(i);
                mean_empb += empb.get(i);
                mean_onlw += onlw.get(i);
                mean_onlb += onlb.get(i);
                mean_empz += empz.get(i);
                mean_onlz += onlz.get(i);
                mean_obj += obj.get(i);
            }
            mean_accuracy /= accuracy.size();
            mean_time /= times.size();
            mean_empw /= empw.size();
            mean_empb /= empb.size();
            mean_onlw /= onlw.size();
            mean_onlb /= onlb.size();
            mean_empz /= empz.size();
            mean_onlz /= onlz.size();
            mean_obj /= obj.size();
            for (int i = 0; i < times.size(); i++) {
                std_time += (mean_time - times.get(i)) * (mean_time - times.get(i));
                std_accuracy += (mean_accuracy - accuracy.get(i)) * (mean_accuracy - accuracy.get(i));
                std_empw += (mean_empw - empw.get(i)) * (mean_empw - empw.get(i));
                std_onlw += (mean_onlw - onlw.get(i)) * (mean_onlw - onlw.get(i));
                std_empz += (mean_empz - empz.get(i)) * (mean_empz - empz.get(i));
                std_onlz += (mean_onlz - onlz.get(i)) * (mean_onlz - onlz.get(i));
                std_empb += (mean_empb - empb.get(i)) * (mean_empb - empb.get(i));
                std_onlb += (mean_onlb - onlb.get(i)) * (mean_onlb - onlb.get(i));
                std_obj += (mean_obj - obj.get(i)) * (mean_obj - obj.get(i));
            }
            std_time /= times.size();
            std_accuracy /= accuracy.size();
            std_empw /= empw.size();
            std_empb /= empb.size();
            std_onlw /= onlw.size();
            std_onlb /= onlb.size();
            std_empb = Math.sqrt(std_empb);
            std_onlb = Math.sqrt(std_onlb);
            std_empw = Math.sqrt(std_empw);
            std_onlw = Math.sqrt(std_onlw);
            std_empz /= empz.size();
            std_onlz /= onlz.size();
            std_obj /= obj.size();
            std_empz = Math.sqrt(std_empz);
            std_onlz = Math.sqrt(std_onlz);
            std_accuracy = Math.sqrt(std_accuracy);
            std_time = Math.sqrt(std_time);
            std_obj = Math.sqrt(std_obj);
            result[n][0] = mean_accuracy;
            result[n][1] = std_accuracy;
            result[n][2] = mean_time;
            result[n][3] = std_time;
            result[n][10] = mean_empw;
            result[n][11] = std_empw;
            result[n][12] = mean_empz;
            result[n][13] = std_empz;
            result[n][14] = mean_empb;
            result[n][15] = std_empb;
            result[n][4] = mean_onlw;
            result[n][5] = std_onlw;
            result[n][6] = mean_onlz;
            result[n][7] = std_onlz;
            result[n][8] = mean_onlb;
            result[n][9] = std_onlb;
            result[n][16] = mean_obj;
            result[n][17] = std_obj;
        }
        return result;
    }

    public static void printwithclass_all(String nclass[], String file, int features, double C, int max_iter, int numberoffources, int limit) throws FileNotFoundException, IOException, InterruptedException {

        for (int n = 0; n < nclass.length; n++) {
            ArrayList<Double> accuracy = new ArrayList(), times = new ArrayList<>();
            for (int f = 1; f <= 4; f++) {
                String sample_train = file + "_train" + nclass[n] + "_" + f + (file.equals("pancan") ? ".csv" : ".xlsx"),
                        sample_test = file + "_test" + nclass[n] + "_" + f + (file.equals("pancan") ? ".csv" : ".xlsx");
                File train_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_train),
                        test_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_test);
                if (!train_file.exists() || !test_file.exists()) {
                    continue;
                }
                int y[][] = new int[1][], y_test[][] = new int[1][];
                y[0] = null;
                y_test[0] = null;
                double x[][][] = new double[1][][], x_test[][][] = new double[1][][];
                x[0] = null;
                x_test[0] = null;
                readExcel(train_file, test_file, file.equals("pancan"), y, y_test, x, x_test);
                Eps_all.setVars(x[0], y[0], C, features);

                double tbw[] = new double[1];
                //int numberoffources = 10;
                //int max_iter = 100000;
                //int limit = (int) (max_iter / 150.0);
                FoodSource_all foodsources[] = new FoodSource_all[numberoffources];
                boolean fs_left[] = new boolean[numberoffources];
                //System.out.println("pahse 1");
                long start = System.currentTimeMillis();
                for (int i = 0; i < numberoffources; i++) {
                    fs_left[i] = false;
                    foodsources[i] = new FoodSource_all(foodsources, tbw, max_iter, 10, fs_left, limit);
                    tbw[0] += 1.0 / foodsources[i].fitness;
                }
                //System.out.println("Phase 2");
                ExecutorService es = Executors.newCachedThreadPool();
                for (FoodSource_all i : foodsources) {
                    es.execute(i);
                }

                es.shutdown();
                es.awaitTermination(1, TimeUnit.DAYS);
                start = System.currentTimeMillis() - start;
                int i = 0;
                //System.out.println("Conpleted");
                foodsources[i].updatewz();
                //foodsources[i].abandon();
                for (int j = 1; j < numberoffources; j++) {
                    foodsources[j].updatewz();
                    //foodsources[j].abandon();
                    if (foodsources[j].fitness < foodsources[i].fitness) {
                        i = j;
                    }
                }
                times.add(start / 1000.0);
                int y_pred[] = Eps_all.predict(foodsources[i].bestw, foodsources[i].bestz, x_test[0]);
                accuracy.add(Eps_all.accuracy_rate(y_test[0], y_pred) * 100.0);
                //System.out.print(f + " ");
                //System.out.println(foodsources[i].best);
                //System.out.println(Eps_all.cost(foodsources[i].bestw, foodsources[i].bestz));
                //System.out.println(Eps.accuracy_rate(y_test, y_pred));
                //System.out.println(start / 1000.0 + " " + foodsources[i].r);
            }
            double mean_time = 0, std_time = 0, mean_accuracy = 0, std_accuracy = 0;
            for (int i = 0; i < times.size(); i++) {
                mean_accuracy += accuracy.get(i);
                mean_time += times.get(i);
            }
            mean_accuracy /= accuracy.size();
            mean_time /= times.size();
            for (int i = 0; i < times.size(); i++) {
                std_time += (mean_time - times.get(i)) * (mean_time - times.get(i));
                std_accuracy += (mean_accuracy - accuracy.get(i)) * (mean_accuracy - accuracy.get(i));
            }
            std_time /= times.size();
            std_accuracy /= accuracy.size();
            std_accuracy = Math.sqrt(std_accuracy);
            std_time = Math.sqrt(std_time);
            System.out.println("\n" + file);
            System.out.println("nclass:" + nclass[n]);
            System.out.println("#Features:" + Eps_all.n);
            System.out.println("Accuracy:" + String.format("%.2f", mean_accuracy) + "\u00B1" + String.format("%.2f", std_accuracy));
            System.out.println("Time:" + String.format("%.2f", mean_time) + "\u00B1" + String.format("%.2f", std_time));
            System.out.println("==========================================");
        }
    }

    public static double[] printwithoutclass(String file, int features, double C, int max_iter,
            int numberoffources, double wz, double cooling, boolean parallel) throws IOException, InterruptedException {
        double result[] = new double[18];
        ArrayList<Double> accuracy = new ArrayList(), times = new ArrayList<>(), obj = new ArrayList<>();
        ArrayList<Integer> empw = new ArrayList<>(), empz = new ArrayList<>(), empb = new ArrayList<>(),
                onlw = new ArrayList<>(), onlz = new ArrayList<>(), onlb = new ArrayList<>();
        HashMap<Integer, Integer> selectedfeatures = new HashMap<>();
        for (int f = 1; f <= 10; f++) {
            String sample_train = file + (file.equals("lung") || file.equals("K3B") ? "" : file.equals("opt") || file.equals("pems") || file.equals("sector") || 
                    file.equals("mll") || file.equals("cll_sub_111") || file.contains("Synthetic") ? "_train" : "_train_" + f) + (file.equals("pancan") || file.equals("pems") || 
                    file.equals("sector") || file.equals("mll") || file.equals("cll_sub_111") ||file.contains("Synthetic") ? ".csv" : ".xlsx"),
                    sample_test = file + (file.equals("opt") || file.equals("pems") || file.equals("sector") || 
                    file.contains("Synthetic") || file.equals("mll") || file.equals("cll_sub_111") ? "_test" : "_test_" + f) + (file.equals("pancan") || file.equals("pems") || file.equals("sector") || 
                    file.equals("mll") || file.equals("cll_sub_111") || file.contains("Synthetic") ? ".csv" : ".xlsx");

            File train_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_train),
                    test_file = new File("/Users/scipio/PycharmProjects/Multiclass_low_SVM/" + sample_test);
            boolean no_test = !test_file.exists();
            if (!train_file.exists() && no_test) {
                continue;
            }

            int y[][] = new int[1][], y_test[][] = new int[1][];
            y[0] = null;
            y_test[0] = null;
            double x[][][] = new double[1][][], x_test[][][] = new double[1][][];
            x[0] = null;
            x_test[0] = null;
            readExcel(train_file, no_test ? null : test_file, (file.equals("pancan") || file.equals("pems") || file.equals("sector") || 
                    file.equals("mll") || file.equals("cll_sub_111") || file.contains("Synthetic")), y, y_test, x, x_test);
            Eps.setVars(x[0], y[0], C, features);
            double tbw[] = new double[1];
            //int numberoffources = 10;
            //int max_iter = 100000;
            //int limit = (int) (max_iter / 150.0);
            if (parallel) {
                FoodSource foodsources[] = new FoodSource[numberoffources];
                boolean fs_left[] = new boolean[numberoffources];
                //System.out.println("pahse 1");
                long start = System.currentTimeMillis();
                for (int i = 0; i < numberoffources; i++) {
                    fs_left[i] = false;
                    foodsources[i] = new FoodSource(foodsources, tbw, max_iter, 10, fs_left, wz, cooling);
                    tbw[0] += 1.0 / foodsources[i].fitness;
                }
                //System.out.println("Phase 2");
                ExecutorService es = Executors.newCachedThreadPool();
                for (FoodSource i : foodsources) {
                    es.execute(i);
                }

                es.shutdown();
                es.awaitTermination(1, TimeUnit.DAYS);
                start = System.currentTimeMillis() - start;
                int i = 0;
                //System.out.println("Conpleted");
                foodsources[i].updatewz();
                //foodsources[i].abandon();
                for (int j = 1; j < numberoffources; j++) {
                    foodsources[j].updatewz();
                    //foodsources[j].abandon();
                    if (foodsources[j].fitness < foodsources[i].fitness) {
                        i = j;
                    }
                }
                times.add(start / 1000.0);
                empw.add(foodsources[i].employebeeupdatew);
                empz.add(foodsources[i].employebeeupdatez);
                empb.add(foodsources[i].employebeeupdateb);
                onlw.add(foodsources[i].onlookerbeeupdatew);
                onlz.add(foodsources[i].onlookerbeeupdatez);
                onlb.add(foodsources[i].onlookerbeeupdateb);
                obj.add(foodsources[i].best);
                int y_pred[] = Eps.predict(foodsources[i].bestw, foodsources[i].bestz, x_test[0]);
                accuracy.add(Eps.accuracy_rate(y_test[0], y_pred) * 100.0);
                for (int j = 0; j < Eps.N; j++) {
                    if (foodsources[i].bestz[j] == 1) {
                        int val = selectedfeatures.containsKey(j) ? selectedfeatures.get(j) : 0;
                        selectedfeatures.put(j, val + 1);
                    }
                }

            } else {
                FoodSource_SQ foodsources[] = new FoodSource_SQ[numberoffources];
                boolean fs_left[] = new boolean[numberoffources];
                //System.out.println("pahse 1");
                long start = System.currentTimeMillis();
                for (int i = 0; i < numberoffources; i++) {
                    fs_left[i] = false;
                    foodsources[i] = new FoodSource_SQ(foodsources, tbw, max_iter, 10, fs_left, wz, cooling);
                    tbw[0] += 1.0 / foodsources[i].fitness;
                }
                //System.out.println("Phase 2");
                for (int counter = 0; counter < max_iter; counter++) {
                    for (FoodSource_SQ i : foodsources) {
                        i.run();
                    }
                }
                start = System.currentTimeMillis() - start;
                int i = 0;
                //System.out.println("Conpleted");
                foodsources[i].updatewz();
                //foodsources[i].abandon();
                for (int j = 1; j < numberoffources; j++) {
                    foodsources[j].updatewz();
                    //foodsources[j].abandon();
                    if (foodsources[j].fitness < foodsources[i].fitness) {
                        i = j;
                    }
                }
                times.add(start / 1000.0);
                empw.add(foodsources[i].employebeeupdatew);
                empz.add(foodsources[i].employebeeupdatez);
                empb.add(foodsources[i].employebeeupdateb);
                onlw.add(foodsources[i].onlookerbeeupdatew);
                onlz.add(foodsources[i].onlookerbeeupdatez);
                onlb.add(foodsources[i].onlookerbeeupdateb);
                obj.add(foodsources[i].best);
                int y_pred[] = Eps.predict(foodsources[i].bestw, foodsources[i].bestz, x_test[0]);
                accuracy.add(Eps.accuracy_rate(y_test[0], y_pred) * 100.0);
            }
            //System.out.print(f + " ");
            //System.out.println(foodsources[i].best);
            //System.out.println(Eps.cost(foodsources[i].bestw, foodsources[i].bestz));
            //System.out.println(Eps.accuracy_rate(y_test, y_pred));
            //System.out.println(start / 1000.0 + " " + foodsources[i].r);
        }/*
        for (int j = 0; j < Eps.N; j++) {
            if (selectedfeatures.containsKey(j)) {
                System.out.print(j + "\t");
            }
        }
        System.out.println("");
        for (int j = 0; j < Eps.N; j++) {
            if (selectedfeatures.containsKey(j)) {
                System.out.print(selectedfeatures.get(j) + "\t");
            }
        }
        System.out.println("");*/
        double mean_time = 0, std_time = 0, mean_accuracy = 0, std_accuracy = 0, mean_empw = 0, mean_empz = 0, mean_empb = 0,
                mean_onlw = 0, mean_onlz = 0, mean_onlb = 0, std_empw = 0, std_empb = 0, std_onlw = 0, std_empz = 0, std_onlz = 0, std_onlb = 0,
                mean_obj = 0, std_obj = 0;
        for (int i = 0; i < times.size(); i++) {
            mean_accuracy += accuracy.get(i);
            mean_time += times.get(i);
            mean_empw += empw.get(i);
            mean_empb += empb.get(i);
            mean_onlw += onlw.get(i);
            mean_onlb += onlb.get(i);
            mean_empz += empz.get(i);
            mean_onlz += onlz.get(i);
            mean_obj += obj.get(i);
        }
        mean_accuracy /= accuracy.size();
        mean_time /= times.size();
        mean_empw /= empw.size();
        mean_empb /= empb.size();
        mean_onlw /= onlw.size();
        mean_onlb /= onlb.size();
        mean_empz /= empz.size();
        mean_onlz /= onlz.size();
        mean_obj /= obj.size();
        for (int i = 0; i < times.size(); i++) {
            std_time += (mean_time - times.get(i)) * (mean_time - times.get(i));
            std_accuracy += (mean_accuracy - accuracy.get(i)) * (mean_accuracy - accuracy.get(i));
            std_empw += (mean_empw - empw.get(i)) * (mean_empw - empw.get(i));
            std_onlw += (mean_onlw - onlw.get(i)) * (mean_onlw - onlw.get(i));
            std_empz += (mean_empz - empz.get(i)) * (mean_empz - empz.get(i));
            std_onlz += (mean_onlz - onlz.get(i)) * (mean_onlz - onlz.get(i));
            std_empb += (mean_empb - empb.get(i)) * (mean_empb - empb.get(i));
            std_onlb += (mean_onlb - onlb.get(i)) * (mean_onlb - onlb.get(i));
            std_obj += (mean_obj - obj.get(i)) * (mean_obj - obj.get(i));
        }
        std_time /= times.size();
        std_accuracy /= accuracy.size();
        std_empw /= empw.size();
        std_empb /= empb.size();
        std_onlw /= onlw.size();
        std_onlb /= onlb.size();
        std_empb = Math.sqrt(std_empb);
        std_onlb = Math.sqrt(std_onlb);
        std_empw = Math.sqrt(std_empw);
        std_onlw = Math.sqrt(std_onlw);
        std_empz /= empz.size();
        std_onlz /= onlz.size();
        std_obj /= obj.size();
        std_empz = Math.sqrt(std_empz);
        std_onlz = Math.sqrt(std_onlz);
        std_accuracy = Math.sqrt(std_accuracy);
        std_time = Math.sqrt(std_time);
        std_obj = Math.sqrt(std_obj);
        result[0] = mean_accuracy;
        result[1] = std_accuracy;
        result[2] = mean_time;
        result[3] = std_time;
        result[10] = mean_empw;
        result[11] = std_empw;
        result[12] = mean_empz;
        result[13] = std_empz;
        result[14] = mean_empb;
        result[15] = std_empb;
        result[4] = mean_onlw;
        result[5] = std_onlw;
        result[6] = mean_onlz;
        result[7] = std_onlz;
        result[8] = mean_onlb;
        result[9] = std_onlb;
        result[16] = mean_obj;
        result[17] = std_obj;
        //System.out.println("");
        return result;
    }
}
