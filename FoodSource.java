/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package multi_svm;

import java.util.ArrayList;
import java.util.Random;

/**
 *
 * @author scipio
 */
public class FoodSource implements Runnable {

    double fitness, r, tbw[], w[][], best, bestw[][], wz, cooling;
    boolean ended, lfs[];
    int max_iter, features, employebeeupdateb, employebeeupdatew, employebeeupdatez, onlookerbeeupdateb, onlookerbeeupdatew, onlookerbeeupdatez;
    FoodSource foodsources[];
    int z[], scoutcount, bestz[];
    private static Random rnd;
    Eps E[];

    public synchronized static FoodSource selectFoodSource(FoodSource fsources[], double prob, double total, FoodSource fs) {
        double totalprobabilty = 0;
        int i = -1;
        while (totalprobabilty < prob) {
            if (fsources[++i] != fs) {
                totalprobabilty += 1 / (fsources[i].best * total);
            }
        }
        return fsources[i];
    }

    public static int[] pickSample(int popsize, int nSamplesNeeded, Random r) {
        int ret[] = new int[nSamplesNeeded];
        int nPicked = 0, i = 0, nLeft = popsize;
        while (nSamplesNeeded > 0) {
            int rand = r.nextInt(nLeft);
            if (rand < nSamplesNeeded) {
                ret[nPicked++] = i;
                nSamplesNeeded--;
            }
            nLeft--;
            i++;
            /*
            if (i > popsize) {
                System.out.println("picksample");
            }*/
        }
        return ret;
    }

    public FoodSource(FoodSource fs[], double tbw[], int max_iter, double r, boolean left_foodsources[], double wz, double cooling) {
        //Thread.__init__(self)
        this.r = r;
        this.cooling = cooling;
        rnd = new Random();
        this.wz = wz;
        ended = false;
        lfs = left_foodsources;
        this.tbw = tbw;
        this.max_iter = max_iter;
        //self.lck = lck
        foodsources = fs;
        w = new double[Eps.N + 1][Eps.K];
        bestw = new double[Eps.N + 1][Eps.K];
        z = new int[Eps.N];
        bestz = new int[Eps.N];
        features = rnd.nextInt(Eps.n) + 1;
        int indices[] = pickSample(Eps.N, features, rnd);
        for (int j = 0; j <= Eps.N; j++) {
            for (int k = 0; k < Eps.K; k++) {
                w[j][k] = rnd.nextDouble() * 2 - 1;
                bestw[j][k] = w[j][k];
            }
            if (j < features) {
                z[indices[j]] = 1;
                bestz[indices[j]] = 1;
            }
        }

        calculate();
        best = fitness;

    }

    public void calculate() {
        E = new Eps[Eps.M];
        double teta = 0;
        for (int i = 0; i < Eps.M; i++) {
            E[i] = new Eps(w, i, z);
            if (E[i].val < 1) {
                teta += 1 - E[i].val;
            }
        }
        fitness = Eps.C * teta;
        for (int j = 0; j < Eps.N; j++) {
            for (int k = 0; k < Eps.K; k++) {
                fitness += z[j] * w[j][k] * w[j][k];
            }
        }
        scoutcount = 0;
    }

    public double Dw(double dw, int j, int k) {
        double teta = 0;
        for (Eps i : E) {
            teta += i.Dw(dw, j, k);
        }
        scoutcount += 1;
        return teta;
    }

    public double Dz(int dz, int j, boolean second, double tempw[]) {
        double teta = 0;
        for (Eps i : E) {
            teta += i.Dz(dz, w, j, second, tempw);
        }
        scoutcount += 1;
        return teta;
    }

    public void update(int j, int k, double dw, double df) {
        //f = self.fitness
        for (Eps i : E) {
            i.update(k);
        }
        w[j][k] += dw;
        fitness += df;
        scoutcount = 0;
    }

    public void updatez(int j, int dz, double df, double df2, int j2, double tempw[]) {
        //#f = self.fitness
        for (Eps i : E) {
            i.update_z();
        }
        z[j] += dz;
        if (tempw != null)
            for (int i = 0; i < tempw.length; i++)
                w[j][i] = tempw[i];
        if (j2 >= 0) {
            z[j2] -= dz;
        } else {
            features += dz;
        }
        fitness += df + df2;
        scoutcount = 0;
        //#return  (1 / self.fitness) - (1 / f)
    }
/*
    public void abandon(double ro, int j, int k) {
        if (ro != 0) {
            double dw = ro * w[j][k];
            double teta = Dw(dw, j, k);
            double df = Eps.C * teta;
            if (j < Eps.N) {
                df += dw * (2.0 * w[j][k] + dw);
            }
            update(j, k, dw, df);
            if (fitness < best) {
                updatetbw(this, (1 / fitness) - (1 / best), fitness, w, z, false, 0);
            }
        } else {
            int dz = -1;
            if (z[j] == 0) {
                dz = 1;
            }
            if (features + dz > 0 && features + dz <= Eps.n) {
                double teta = Dz(dz, j, false);
                double df = 0;
                for (int kk = 0; kk < Eps.K; kk++) {
                    df += dz * w[j][kk] * w[j][kk];
                }
                updatez(j, dz, df + Eps.C * teta, 0, -1, null);
                if (fitness < best) {
                    updatetbw(this, (1 / fitness) - (1 / best), fitness, w, z, false, 0);
                }

            }
        }


        
    }*/

    /*
        FoodSource neighbour = selectFoodSource(foodsources, rnd.nextDouble(), tbw[0] - (1 / best), this);
       
        features = 0;*/
 /*for (int j = 0; j <= Eps.N; j++) {
            FoodSource nf = this;
           
            if (rnd.nextBoolean()) {
                nf = neighbour;
            } 
            if (j < Eps.N) {
                z[j] = nf.bestz[j];
                features += z[j];
            }
            for (int k = 0; k < Eps.K; k++) {
                w[j][k] = nf.bestw[j][k];
            }
        }*/
//System.out.print("This:" + me + " neighbour:" + she + " ");
/*int update = rnd.nextInt((int)(Eps.n * 0.05));
        ArrayList<Integer> updates = new ArrayList<>();
        for (int u = 0; u < update; u++) {
            int j;
            do {
                j = rnd.nextInt(Eps.N + 1);
            } while (updates.contains(j));
            updates.add(j);
            if (j < Eps.N) {
                int dz = neighbour.bestz[j] - z[j];
                z[j] = neighbour.bestz[j];
                features += dz;
            }
            for (int k = 0; k < Eps.K; k++) {
                w[j][k] = neighbour.bestw[j][k];
            }
        }
        if (features < 1) {
            int j = rnd.nextInt(Eps.N);//(0, Eps.N - 1)
            if (z[j] == 1) {
                j = (j + 1) % (Eps.N - 1);
            }
            z[j] = 1;
            features = 1;
        }
        while (features > Eps.n) {
            int j = rnd.nextInt(Eps.N);
            while (z[j] == 0) {
                j = (j + 1) % (Eps.N - 1);
            }
            z[j] = 0;
            features -= 1;
        }
        calculate();*/
//#return (1 / self.fitness) - (1 / f)
    public synchronized static void updatetbw(FoodSource fs, double tbw, double b, double bw[][], int bz[], boolean checkalltbw, int step_no) {
        if (checkalltbw) {
            if (fs.lfs[step_no]) {
                return;
            }
            boolean sefisworst = true;
            for (FoodSource i : fs.foodsources) {
                if (i.ended == false && i.best > fs.best) {
                    sefisworst = false;
                    break;
                }
            }
            if (sefisworst) {
                fs.lfs[step_no] = true;
                fs.ended = true;
            }
        } else {
            fs.tbw[0] += tbw;
            fs.best = b;
            for (int j = 0; j <= Eps.N; j++) {
                if (j < Eps.N) {
                    fs.bestz[j] = bz[j];
                }
                for (int k = 0; k < Eps.K; k++) {
                    fs.bestw[j][k] = bw[j][k];
                }
            }
        }
    }

    public void updatewz() {
        fitness = best;
        for (int j = 0; j <= Eps.N; j++) {
            if (j < Eps.N) {
                z[j] = bestz[j];
            }
            for (int k = 0; k < Eps.K; k++) {
                w[j][k] = bestw[j][k];
            }
        }
    }

    public void employedBeePhase(double ro, int j, int k) {
        FoodSource nf = selectFoodSource(foodsources, rnd.nextDouble(), tbw[0] - (1 / best), this);
        if (ro != 0) {
            ro = (ro + r) / (2 * r);
            /*
            int c = 0, j2 = rnd.nextInt(Eps.N);
            while (c++ < Eps.N && nf.bestz[j2] == 0 && z[j2] == 1) {
                j2 = (j2 + 1) % Eps.N;
            }
            if (nf.bestz[j2] == 1 && z[j2] == 0) {
                    w[j2][k] += ro * (w[j2][k] - nf.bestw[j2][k]);
            }
            */
            while (j < Eps.N && z[j] == 0) {
                j = (j + 1) % Eps.N;
                /*if (c++ > Eps.N) {
                    System.out.println("employedbee");
                }*/
            }

            double dw = ro * (w[j][k] - nf.bestw[j][k]);
            double teta = Dw(dw, j, k);
            double df = Eps.C * teta;
            if (j < Eps.N) {
                df += dw * (2.0 * w[j][k] + dw);
            }

            if (df < 0) {
                if (j == Eps.N) {
                    employebeeupdateb += 1;
                } else {
                    employebeeupdatew += 1;
                }
                update(j, k, dw, df);
                if (fitness < best) {
                    updatetbw(this, (1 / fitness) - (1 / best), fitness, w, z, false, 0);
                }
            }
        } else {
            j = rnd.nextInt(Eps.N);
            int c = 0;
            while (z[j] == nf.bestz[j]) {
                j = (j + 1) % Eps.N;
                if (c++ > Eps.N) {
                    return;
                }
            }
            if (movez(j, nf.bestw[j])) {
                employebeeupdatez += 1;
            }
        }
    }

    public boolean movez(int j, double tempw[]) {
        int dz = -1;
        double myw[] = w[j];
        if (z[j] == 0) {
            dz = 1;
            if (tempw != null)
                myw = tempw;
        }

        double teta = Dz(dz, j, false, myw);
        double df = 0;
        for (int k = 0; k < Eps.K; k++) {
            df += dz * myw[k] * myw[k];
        }
        if (features + dz > 0 && features + dz <= Eps.n && df + Eps.C * teta < 0) { //  # or uniform(0, 1) <= exp(df / -r):
            updatez(j, dz, df + Eps.C * teta, 0, -1, myw);
            if (fitness < best) {
                updatetbw(this, (1 / fitness) - (1 / best), fitness, w, z, false, 0);
            }
            return true;
        }

        int j2 = rnd.nextInt(Eps.N);
        //int c = 0;
        while (z[j2] == z[j]) {
            j2 = (j2 + 1) % Eps.N;
            /*if (c++ > Eps.N) {
                System.out.println("movez");
            }*/
        }
        teta = Dz(-dz, j2, true, null);
        double df2 = Eps.C * teta;
        for (int k = 0; k < Eps.K; k++) {
            df2 += -dz * w[j2][k] * w[j2][k];
        }
        if (df + df2 < 0) {//  # or uniform(0, 1) <= exp(df / -r):
            updatez(j, dz, df, df2, j2, myw);
            if (fitness < best) {
                updatetbw(this, (1 / fitness) - (1 / best), fitness, w, z, false, 0);
            }
            return true;
        }

        return false;
    }

    public boolean onlookerBeePhase(double ro, int j, int k) {
        if (ro != 0) {
            double dw = ro * w[j][k];
            double teta = Dw(dw, j, k);
            double df = Eps.C * teta;
            if (j < Eps.N) {
                df += dw * (2.0 * w[j][k] + dw);
            }
            if (df < 0) {
                if (j == Eps.N) {
                    onlookerbeeupdateb += 1;
                } else {
                    onlookerbeeupdatew += 1;
                }
                update(j, k, dw, df);
                if (fitness < best) {
                    updatetbw(this, (1 / fitness) - (1 / best), fitness, w, z, false, 0);
                }
                return true;
            }
        } else {
            boolean res = movez(rnd.nextInt(Eps.N), null);
            if (res) {
                onlookerbeeupdatez += 1;
            }
            return res;
        }
        return false;
    }
    /*
    public void scoutBeePhase() {
        if (limit < scoutcount) {
            int j = rnd.nextInt(Eps.N + 1);
            if (j == Eps.N || rnd.nextDouble() < 0.5) {
                double ro = rnd.nextDouble() * 2 * r - r;
                while (j < Eps.N && z[j] == 0) {
                    j = (j + 1) % (Eps.N);
                    
                }
                r *= 0.99995;
                abandon(ro, j, rnd.nextInt(Eps.K));
            } else {
                abandon(0, j, -1);
            }

            //System.out.println("abandon from " + f + " to " + fitness);
            if (fitness < best) {
                updatetbw(this, (1 / fitness) - (1 / best), fitness, w, z, false, 0);
            }
            scoutcount = 0;
        }
    }
    */
    @Override
    public void run() {
        int counter = 0;
        //int step_size = (int) (max_iter / (double) foodsources.length);
        //int step_no = 0;
        while (counter < max_iter) {
            int j = rnd.nextInt(Eps.N + 1);
            double ro;
           
            if (rnd.nextDouble() < wz) {
                ro = rnd.nextDouble() * 2 * r - r;
                while (j < Eps.N && z[j] == 0) {
                    j = (j + 1) % Eps.N;
                    
                }
                r *= cooling;
            } else {
                ro = 0;
            }
            int k = rnd.nextInt(Eps.K);
            
            employedBeePhase(ro, j, k);
            j = rnd.nextInt(Eps.N + 1);
            if (ro != 0) {
                //int c = 0;
                while (j < Eps.N && z[j] == 0) {
                    j = (j + 1) % (Eps.N);
                    /*if (c++ > Eps.N) {
                        System.out.println("run2");
                    }*/
                }
            }
            //employedBeePhase(ro, j, k);
            onlookerBeePhase(ro, j, k);
            //scoutBeePhase();
            counter += 1;

            /*
            if (lfs[step_no]) {
                step_no += 1;
            } else if (counter >= (step_no + 1) * step_size) {
                updatetbw(this, 0, 0, w, z, true, step_no);
                if (ended) {
                    System.out.println(best + "\t" + counter);
                    break;
                }
            }
             */
        }
    }

}
