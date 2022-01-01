/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package multi_svm;

import java.util.ArrayList;

/**
 *
 * @author scipio
 */
public class Eps {

    public static int K, N, M, n, y[];
    public static double x[][], C;
    public int i, z[], q[], dq[];
    public double w[][], dE, de[], e[], val;
    public boolean resort;

    public static void setVars(double independentVar[][], int dependentVar[], double penalty, int nfeatures) {
        x = independentVar;
        y = dependentVar;
        //K = np.unique(y).size
        ArrayList<Integer> unique = new ArrayList<>();
        for (int i : y) {
            if (!unique.contains(i)) {
                unique.add(i);
            }
        }
        K = unique.size();
        M = x.length;
        N = x[0].length;
        C = penalty;
        n = nfeatures;
    }

    public static int[] predict(double w[][], int z[], double x_test[][]) {
        int result[] = new int[x_test.length];
        for (int i = 0; i < result.length; i++) {
            double best = 0;
            for (int j = 0; j < z.length; j++) {
                best += z[j] * w[j][0] * x_test[i][j];
            }
            best += w[z.length][0];
            int pred = 0;
            for (int k = 1; k < w[0].length; k++) {
                double calc = 0;
                for (int j = 0; j < z.length; j++) {
                    calc += z[j] * w[j][k] * x_test[i][j];
                }
                calc += w[z.length][k];
                if (calc > best) {
                    best = calc;
                    pred = k;
                }
            }
            result[i] = pred;
        }
        return result;
    }

    public static double accuracy_rate(int y_true[], int y_prediction[]) {
        double t = 0;
        for (int i = 0; i < y_true.length; i++) {
            if (y_true[i] == y_prediction[i]) {
                t += 1;
            }
        }
        return t / y_true.length;
    }

    public static double cost(double w[][], int z[]) {
        double e[] = new double[K];
        double teta = 0;
        double f = 0;
        for (int i = 0; i < M; i++) {
            int qi = -1;
            for (int k = 0; k < K; k++) {
                e[k] = 0;
                for (int j = 0; j < z.length; j++) {
                    if (i == 0) {
                        f += z[j] * w[j][k] * w[j][k];
                    }
                    e[k] += z[j] * w[j][k] * x[i][j];
                }
                e[k] += w[z.length][k];
                if (y[i] != k && (qi == -1 || e[k] > e[qi])) {
                    qi = k;
                }
            }
            double E = e[y[i]] - e[qi];

            if (E < 1) {
                teta += 1 - E;
            }
        }
        return f + C * teta;
    }

    public static int[] getQ(double e[], int yi) {
        int q[] = new int[2];
        q[0] = -1;
        q[1] = -1;
        for (int i = 0; i < e.length; i++) {
            if (i != yi) {
                if (q[0] == -1) {
                    q[0] = i;
                } else if (e[i] > e[q[0]]) {
                    q[1] = q[0];
                    q[0] = i;
                } else if (q[1] == -1 || e[i] > e[q[1]]) {
                    q[1] = i;
                }
            }
        }
        
        return q;
    }

    public static Eps[] calculate_E(double w[][], int z[]) {
        Eps E[] = new Eps[M];
        for (int i = 0; i < M; i++) {
            E[i] = new Eps(w, i, z);
        }
        return E;
    }

    public Eps(double w[][], int i, int z[]) {
        this.i = i;
        this.resort = false;
        this.dE = 0;
        this.e = new double[K];
        this.de = new double[K];

        for (int k = 0; k < K; k++) {
            e[k] = 0;
            for (int j = 0; j < z.length; j++) {
                e[k] += z[j] * w[j][k] * x[i][j];
            }
            e[k] += w[z.length][k];
        }
        q = getQ(e, y[i]);
        dq = new int[2];
        
        val = e[y[i]] - e[q[0]];
    }

    

    public double Dz(int dz, double w[][], int j, boolean second, double tempw[]) {
        dq[0] = -1;
        dq[1] = -1;
        double myw[] = (tempw != null) ? tempw : w[j];
        for (int k = 0; k < K; k++) {
            if (second) {
                de[k] += dz * myw[k] * x[i][j];
            } else {
                de[k] = dz * myw[k] * x[i][j];
            }
            if (k != y[i]) {
                if (dq[0] == -1) {
                    dq[0] = k;
                } else if (e[k] + de[k] > e[dq[0]] + de[dq[0]]) {
                    dq[1] = dq[0];
                    dq[0] = k;
                } else if (dq[1] == -1 || e[k] + de[k] > e[dq[1]] + de[dq[1]]) {
                    dq[1] = k;
                }
            }
        }

        dE = de[y[i]] + e[q[0]] - e[dq[0]] - de[dq[0]];

        if (val < 1 && val + dE >= 1) {
            return val - 1;
        }
        if (val >= 1 && val + dE < 1) {
            return 1 - val - dE;
        }
        if (val < 1 && val + dE < 1) {
            return -dE;
        }
        return 0;
    }

    

    public double Dw(double dw, int j, int k) {
        if (j < N) {
            de[k] = dw * x[i][j];
        } else {
            de[k] = dw;
        }
        resort = true;
        if (k == y[i]) {
            dE = de[k];
            resort = false;
        } else if (k == q[0] && e[k] + de[k] >= e[q[1]]) {
            dE = -de[k];
            resort = false;
        } else if (k == q[0]) {
            dE = e[q[0]] - e[q[1]];
        } else if (e[k] + de[k] > e[q[0]]) {
            dE = e[q[0]] - e[k] - de[k];
        } else {
            dE = 0;
        }

        if (val < 1 && val + dE >= 1) {
            return val - 1;
        }
        if (val >= 1 && val + dE < 1) {
            return 1 - val - dE;
        }
        if (val < 1 && val + dE < 1) {
            return -dE;
        }
        return 0;
    }

    public void update(int k) {
        e[k] += de[k];
        if (resort) {
            if (e[k] > e[q[0]]) {
                int temp = q[0];
                q[0] = k;
                q[1] = temp;
            } else if (e[k] > e[q[1]]) {
                q[1] = k;
            } else {
                q = getQ(e, y[i]);
            }
        }

        val += dE;
    }

    public void update_z() {
        for (int k = 0; k < K; k++) {
            e[k] += de[k];
        }
        val += dE;
        q[0] = dq[0];
        q[1] = dq[1];
    }

}
