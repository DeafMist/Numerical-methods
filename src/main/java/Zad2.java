import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import java.io.FileOutputStream;
import java.io.IOException;

public final class Zad2 {
    static class Pair {
        public double[][] x;
        public int counter;

        public Pair(double[][] x, int counter) {
            this.x = x;
            this.counter = counter;
        }
    }

    public static double[][] multiplyMatrixNumber(double lambda, double[][] a) {
        int n = a.length;
        int m = a[0].length;
        double[][] res = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                res[i][j] = a[i][j] * lambda;
            }
        }

        return res;
    }

    public static double[][] multiplyMatrixMatrix(double[][] a, double[][] b) {
        int n = a.length;
        int m = a[0].length;
        int l = b[0].length;
        double[][] c = new double[n][l];

        for (int i = 0; i < n; i++) {
            for (int k = 0; k < l; k++) {
                for (int j = 0; j < m; j++) {
                    c[i][k] += a[i][j] * b[j][k];
                }
            }
        }

        return c;
    }

    public static double[][] sumMatrixMatrix(double[][] a, double[][] b) {
        int n = a.length;
        int m = a[0].length;
        double[][] c = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                c[i][j] = a[i][j] + b[i][j];
            }
        }

        return c;
    }

    public static double[][] transpose(double[][] a) {
        int n = a.length;
        int m = a[0].length;

        double[][] c = new double[m][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                c[j][i] = a[i][j];
            }
        }

        return c;
    }

    public static double determinant(double[][] a) {
        int n = a.length;
        int sign = 1;
        double res = 0;

        if (n == 1) {
            return a[0][0];
        }

        if (n == 2) {
            return (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
        }

        double[][] b = new double[n - 1][n - 1];

        for (int i = 0; i < n; i++) {
            int i_cur = 0;
            int j_cur = 0;

            for (int k = 1; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    if (l == i) {
                        continue;
                    }

                    b[i_cur][j_cur] = a[k][l];
                    j_cur++;
                    if (j_cur == n - 1) {
                        j_cur = 0;
                        i_cur++;
                    }
                }
            }

            res += sign * a[0][i] * determinant(b);
            sign *= -1;
        }

        return res;
    }

    public static double[][] reverse(double[][] a) {
        int n = a.length;
        double det = determinant(a);
        double[][] res = new double[n][n];
        double[][] b = new double[n - 1][n - 1];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int i_cur = 0;
                int j_cur = 0;

                for (int k = 0; k < n; k++) {
                    if (k == i) {
                        continue;
                    }

                    for (int l = 0; l < n; l++) {
                        if (l == j) {
                            continue;
                        }

                        b[i_cur][j_cur] = a[k][l];
                        j_cur++;
                        if (j_cur == n - 1) {
                            i_cur++;
                            j_cur = 0;
                        }
                    }
                }

                res[i][j] = ((i + j) % 2 == 0 ? determinant(b) : -determinant(b));
            }
        }

        res = multiplyMatrixNumber(1.0 / det, transpose(res));

        return res;
    }

    public static double[][] copy(double[][] a) {
        int n = a.length;
        int m = a[0].length;
        double[][] res = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                res[i][j] = a[i][j];
            }
        }

        return res;
    }

    public static double matrixNorm1(double[][] a) {
        int n = a.length;
        int m = a[0].length;
        double max = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < m; i++) {
            double cur = 0;

            for (int j = 0; j < n; j++) {
                cur += Math.abs(a[j][i]);
            }

            if (cur > max) {
                max = cur;
            }
        }

        return max;
    }

    public static double matrixNorm2(double[][] a) {
        int n = a.length;
        int m = a[0].length;
        double res = 0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                res += Math.pow(a[i][j], 2);
            }
        }
        res = Math.sqrt(res);

        return res;
    }

    public static double matrixNorm3(double[][] a) {
        int n = a.length;
        int m = a[0].length;
        double max = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < n; i++) {
            double cur = 0;

            for (int j = 0; j < m; j++) {
                cur += Math.abs(a[i][j]);
            }

            if (cur > max) {
                max = cur;
            }
        }

        return max;
    }

    private static double pickOfNorm(double[][] a) {
        double norm;

        if (matrixNorm1(a) < 1) {
            norm = matrixNorm1(a);
        }
        else if (matrixNorm2(a) < 1) {
            norm = matrixNorm2(a);
        }
        else {
            norm = matrixNorm3(a);
        }

        return norm;
    }

    public static boolean isADiagonalPredominance(double[][] a) {
        int n = a.length;
        boolean flag = true;

        for (int i = 0; i < n; i++) {
            double cur = 0;

            for (int j = 0; j < n; j++) {
                if (i != j) {
                    cur += Math.abs(a[i][j]);
                }
            }

            if (Math.abs(a[i][i]) - cur <= 0) {
                flag = false;
                break;
            }
        }

        return flag;
    }

    public static Pair simpleIterationMethod(double[][] a, double[][] b, double epsi) {
        int n = a.length;

        double[][] E = new double[n][n];
        for (int i = 0; i < n; i++) {
            E[i][i] = 1;
        }

        double[][] B;
        double k = 1.0 / matrixNorm1(a);
        B = sumMatrixMatrix(E, multiplyMatrixNumber(-k, a));

        double normB = pickOfNorm(B);
        if (normB >= 1) {
            b = multiplyMatrixMatrix(transpose(a), b);
            a = multiplyMatrixMatrix(transpose(a), a);

            k = 1.0 / matrixNorm1(a);
            B = sumMatrixMatrix(E, multiplyMatrixNumber(-k, a));
            normB = pickOfNorm(B);
        }

        double[][] c = multiplyMatrixNumber(k, b);
        double[][] xi = multiplyMatrixNumber(k, b);

        int counter = 0;
        if (normB < 1) {
            double[][] prev;

            do {
                prev = xi;
                xi = sumMatrixMatrix(multiplyMatrixMatrix(B, xi), c);
                counter++;
            } while (normB * matrixNorm1(sumMatrixMatrix(xi, multiplyMatrixNumber(-1, prev))) / (1 - normB) > epsi);
        }
        else {
            while (matrixNorm1(sumMatrixMatrix(multiplyMatrixMatrix(a, xi), multiplyMatrixNumber(-1, b))) > epsi) {
                xi = sumMatrixMatrix(multiplyMatrixMatrix(B, xi), c);
                counter++;
            }
        }

        return new Pair(xi, counter);
    }

    public static Pair methodSeidel(double[][] a, double[][] b, double epsi) {
        int n = a.length;

        if (!isADiagonalPredominance(a)) {
            b = multiplyMatrixMatrix(transpose(a), b);
            a = multiplyMatrixMatrix(transpose(a), a);
        }

        double[][] c = new double[n][n];
        double[][] d = new double[n][1];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    c[i][j] = -a[i][j] / a[i][i];
                }
            }
        }

        double[][] xi = new double[n][1];
        double[][] prev = new double[n][1];

        for (int i = 0; i < n; i++) {
            d[i][0] = b[i][0] / a[i][i];

            xi[i][0] = b[i][0] / a[i][i];
            prev[i][0] = b[i][0] / a[i][i];
        }

        int counter = 0;
        while (matrixNorm1(sumMatrixMatrix(multiplyMatrixMatrix(a, xi), multiplyMatrixNumber(-1, b))) > epsi) {
            for (int i = 0; i < n; i++) {
                xi[i][0] = 0;

                for (int j = 0; j < n; j++) {
                    xi[i][0] += c[i][j] * prev[j][0];
                }

                xi[i][0] += d[i][0];
                prev[i][0] = xi[i][0];
            }
            counter++;
        }

        return new Pair(xi, counter);
    }

    public static double[][] methodHouseholderQR(double[][] a, double[][] b) {
        int n = a.length;

        double[][] Q = new double[n][n];
        for (int i = 0; i < n; i++) {
            Q[i][i] = 1;
        }

        double[][] R = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                R[i][j] = a[i][j];
            }
        }

        int n_cur = n;
        for (int i = 0; i < n - 1; i++) {
            double[][] E = new double[n_cur][n_cur];
            for (int j = 0; j < n_cur; j++) {
                E[j][j] = 1;
            }

            double[][] zi = new double[n_cur][1];
            zi[0][0] = 1;

            double[][] yi = new double[n_cur][1];
            int i_cur = 0;
            for (int j = i; j < n; j++) {
                yi[i_cur++][0] = R[j][i];
            }

            double[][] wi = multiplyMatrixNumber(
                    1.0 / matrixNorm2(sumMatrixMatrix(yi, multiplyMatrixNumber(-matrixNorm2(yi), zi))),
                    sumMatrixMatrix(yi, multiplyMatrixNumber(-matrixNorm2(yi), zi))
            );

            double[][] Qi = sumMatrixMatrix(E, multiplyMatrixMatrix(multiplyMatrixNumber(-2, wi), transpose(wi)));

            double[][] Ri = new double[n_cur][n_cur];
            i_cur = 0;
            int j_cur = 0;
            for (int j = i; j < n; j++) {
                for (int k = i; k < n; k++) {
                    Ri[i_cur][j_cur] = R[j][k];
                    j_cur++;

                    if (j_cur == n_cur) {
                        j_cur = 0;
                        i_cur++;
                    }
                }
            }

            Ri = multiplyMatrixMatrix(Qi, Ri);

            i_cur = 0;
            j_cur = 0;
            for (int j = i; j < n; j++) {
                for (int k = i; k < n; k++) {
                    R[j][k] = Ri[i_cur][j_cur];
                    j_cur++;

                    if (j_cur == n_cur) {
                        j_cur = 0;
                        i_cur++;
                    }
                }
            }

            double[][] Q0 = new double[n][n];
            for (int j = 0; j < n; j++) {
                Q0[j][j] = 1;
            }

            i_cur = 0;
            j_cur = 0;
            for (int j = i; j < n; j++) {
                for (int k = i; k < n; k++) {
                    Q0[j][k] = Qi[i_cur][j_cur];
                    j_cur++;

                    if (j_cur == n_cur) {
                        j_cur = 0;
                        i_cur++;
                    }
                }
            }

            Q = multiplyMatrixMatrix(Q0, Q);
            n_cur--;
        }

        double[][] y = multiplyMatrixMatrix(Q, b);
        double[][] x = multiplyMatrixMatrix(reverse(R), y);

        return x;
    }

    public static double[][] methodGaussLUP(double[][] a, double[][] b) {
        int n = a.length;

        double[][] P = new double[n][n];
        for (int i = 0; i < n; i++) {
            P[i][i] = 1;
        }

        double[][] M = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] = a[i][j];
            }
        }

        for (int i = 0; i < n - 1; i++) {
            int i_max = i;
            double max = Math.abs(M[i][i]);
            for (int j = i + 1; j < n; j++) {
                if (Math.abs(M[j][i]) > max) {
                    max = M[j][i];
                    i_max = j;
                }
            }

            double temp;
            for (int j = 0; j < n; j++) {
                temp = M[i][j];
                M[i][j] = M[i_max][j];
                M[i_max][j] = temp;

                temp = P[i][j];
                P[i][j] = P[i_max][j];
                P[i_max][j] = temp;
            }

            for (int j = i + 1; j < n; j++) {
                M[j][i] /= M[i][i];
                for (int k = n - 1; k > i; k--) {
                    M[j][k] -= M[i][k] * M[j][i];
                }
            }
        }

        for (int i = 0; i < n; i++) {
            M[i][i] += 1;
        }

        double[][] L = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                if (j == i) {
                    L[i][j] = 1;
                }
                else {
                    L[i][j] = M[i][j];
                }
            }
        }

        double[][] U = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                if (j == i) {
                    U[i][j] = M[i][j] - L[i][j];
                }
                else {
                    U[i][j] = M[i][j];
                }
            }
        }

        double[][] y = multiplyMatrixMatrix(reverse(L), multiplyMatrixMatrix(P, b));
        double[][] x = multiplyMatrixMatrix(reverse(U), y);

        return x;
    }

    public static void print(double[][]a) {
        int n = a.length;
        int m = a[0].length;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                System.out.print(a[i][j] + " ");
            }
            System.out.println();
        }

        System.out.println();
    }

    private static String toString(double[][] a) {
        if (a == null)
            return "null";
        int iMax = a.length - 1;
        if (iMax == -1)
            return "";

        StringBuilder b = new StringBuilder();
        for (int i = 0; ; i++) {
            b.append(a[i][0]);
            if (i == iMax)
                return b.toString();
            b.append(", ");
        }
    }

    private static void makeTitle(Row row0, int numberSheet) {
        Cell cell;
        int i = 0;

        cell = row0.createCell(i++);
        cell.setCellValue("№ теста");

        if (numberSheet == 2) {
            cell = row0.createCell(i++);
            cell.setCellValue("n");

            cell = row0.createCell(i++);
            cell.setCellValue("E");
        }

        cell = row0.createCell(i++);
        cell.setCellValue("x_");

        cell = row0.createCell(i++);
        cell.setCellValue("e");

        cell = row0.createCell(i++);
        cell.setCellValue("MPI x");

        cell = row0.createCell(i++);
        cell.setCellValue("MPI delta");

        cell = row0.createCell(i++);
        cell.setCellValue("MPI k");

        cell = row0.createCell(i++);
        cell.setCellValue("Seidel x");

        cell = row0.createCell(i++);
        cell.setCellValue("Seidel delta");

        cell = row0.createCell(i++);
        cell.setCellValue("Seidel k");

        cell = row0.createCell(i++);
        cell.setCellValue("Gauss x");

        cell = row0.createCell(i++);
        cell.setCellValue("Gauss delta");

        cell = row0.createCell(i++);
        cell.setCellValue("Householder x");

        cell = row0.createCell(i++);
        cell.setCellValue("Householder delta");
    }

    private static void makeTest1(Sheet sheet, double[][] a, double[][] b, int currentRow, int currentTest, double[][] x) {
        int k = 0;
        double epsi = 0;
        double[][] myX = new double[a.length][1];
        Pair pair;

        for (int i = 0; i < 3; i++) {
            Row row = sheet.createRow(currentRow);
            Cell cell;

            if (currentRow % 4 == 1) {
                cell = row.createCell(0);
                cell.setCellValue(currentTest);
            }

            switch (currentRow % 4) {
                case 1 -> epsi = 1e-3;
                case 2 -> epsi = 1e-4;
                case 3 -> epsi = 1e-5;
            }

            for (int cell_i = 1; cell_i <= 12; cell_i++) {
                cell = row.createCell(cell_i);

                switch (cell_i) {
                    case 3 -> {
                        pair = simpleIterationMethod(a, b, epsi);
                        myX = pair.x;
                        k = pair.counter;
                    }
                    case 6 -> {
                        pair = methodSeidel(a, b, epsi);
                        myX = pair.x;
                        k = pair.counter;
                    }
                    case 9 -> myX = methodGaussLUP(a, b);
                    case 11 -> myX = methodHouseholderQR(a, b);
                }

                double delta = matrixNorm1(sumMatrixMatrix(x, multiplyMatrixNumber(-1, myX)));

                switch(cell_i) {
                    case 1 -> cell.setCellValue(x[0][0] + ", " + x[1][0] + ", " + x[2][0]);
                    case 2 -> cell.setCellValue(epsi);
                    case 3, 6, 9, 11 -> cell.setCellValue(myX[0][0] + ", " + myX[1][0] + ", " + myX[2][0]);
                    case 4, 7, 10, 12 -> cell.setCellValue(delta);
                    case 5, 8 -> cell.setCellValue(k);
                }
            }

            currentRow++;
        }
    }

    private static void makeTest2(Sheet sheet, double[][] a1, double[][] a2, double[][] b, int currentRow, int currentN, double[][] x) {
        int k = 0;
        double epsi = 0;
        double e = 0;
        double[][] a = new double[a1.length][a1.length];
        double[][] myX = new double[a1.length][1];
        Pair pair;

        for (int i = 0; i < 4; i++) {
            Row row = sheet.createRow(currentRow);
            Cell cell;

            if (currentRow == 1) {
                cell = row.createCell(0);
                cell.setCellValue(5);
            }

            if (currentRow % 5 == 1) {
                cell = row.createCell(1);
                cell.setCellValue(currentN);
            }

            switch (currentRow % 5) {
                case 1 -> {
                    for (int j = 0; j < a1.length; j++) {
                        for (int l = 0; l < a1.length; l++) {
                            a[j][l] = a1[j][l];
                        }
                    }

                    epsi = 1e-3;
                    e = 1e-3;
                }
                case 2, 4 -> e = 1e-4;
                case 3 -> {
                    for (int j = 0; j < a1.length; j++) {
                        for (int l = 0; l < a1.length; l++) {
                            a[j][l] = a2[j][l];
                        }
                    }

                    epsi = 1e-6;
                    e = 1e-3;
                }
            }

            for (int cell_i = 2; cell_i <= 14; cell_i++) {
                cell = row.createCell(cell_i);

                switch (cell_i) {
                    case 5 -> {
                        pair = simpleIterationMethod(a, b, e);
                        myX = pair.x;
                        k = pair.counter;
                    }
                    case 8 -> {
                        pair = methodSeidel(a, b, e);
                        myX = pair.x;
                        k = pair.counter;
                    }
                    case 11 -> myX = methodGaussLUP(a, b);
                    case 13 -> myX = methodHouseholderQR(a, b);
                }

                double delta = matrixNorm1(sumMatrixMatrix(x, multiplyMatrixNumber(-1, myX)));

                switch(cell_i) {
                    case 2 -> cell.setCellValue(epsi);
                    case 3 -> cell.setCellValue(toString(x));
                    case 4 -> cell.setCellValue(e);
                    case 5, 8, 11, 13 -> cell.setCellValue(toString(myX));
                    case 6, 9, 12, 14 -> cell.setCellValue(delta);
                    case 7, 10 -> cell.setCellValue(k);
                }
            }

            currentRow++;
        }
    }

    private static void makeFirstSheet(Sheet sheet) {
        double[][] a;
        double[][] b;
        double[][] x;

        //Test 0
        a = new double[][]{
                {0, 2, 3},
                {1, 2, 4},
                {4, 5, 6}
        };
        b = new double[][] {
                {13},
                {17},
                {32}
        };
        x = new double[][] {
                {1},
                {2},
                {3}
        };
        makeTest1(sheet, a, b, 1, 0, x);

        //Test 1
        a = new double[][]{
                {12, 1, 1},
                {1, 14, 1},
                {1, 1, 16}
        };
        b = new double[][] {
                {14},
                {16},
                {18}
        };
        x = new double[][] {
                {1},
                {1},
                {1}
        };
        makeTest1(sheet, a, b, 5, 1, x);

        //Test 2
        a = new double[][]{
                {-12, 1, 1},
                {1, -14, 1},
                {1, 1, -16}
        };
        b = new double[][] {
                {-14},
                {-16},
                {-18}
        };
        x = new double[][] {
                {916. / 661},
                {882. / 661},
                {856. / 661}
        };
        makeTest1(sheet, a, b, 9, 2, x);

        //Test 3
        a = new double[][]{
                {-12, 13, 14},
                {15, -14, 11},
                {14, 15, -16}
        };
        b = new double[][] {
                {14},
                {16},
                {18}
        };
        x = new double[][] {
                {3404. / 2577},
                {2902. / 2577},
                {2800. / 2577}
        };
        makeTest1(sheet, a, b, 13, 3, x);

        //Test 4
        a = new double[][]{
                {12, 11, 11},
                {11, 14, 11},
                {11, 11, 16}
        };
        b = new double[][] {
                {14},
                {16},
                {18}
        };
        x = new double[][] {
                {-8. / 67},
                {42. / 67},
                {52. / 67}
        };
        makeTest1(sheet, a, b, 17, 4, x);

        for (int i = 0; i <= 12; i++) {
            sheet.autoSizeColumn(i);
        }
    }

    private static void makeSecondSheet(Sheet sheet) {
        double[][] a1;
        double[][] a2;
        double[][] temp1;
        double[][] temp2;
        double[][] b;
        double[][] x;
        int currentRow = 1;

        for (int n = 4; n <= 10; n++) {
            a1 = new double[n][n];
            a2 = new double[n][n];
            temp1 = new double[n][n];
            temp2 = new double[n][n];
            b = new double[n][1];
            x = new double[n][1];
            x[n - 1][0] = 1;

            for (int i = 0; i < n; i++) {
                b[i][0] = -1;
                for (int j = 0; j < n; j++) {
                    if (j == i) {
                        a1[i][j] = 1;
                        a2[i][j] = 1;

                        temp1[i][j] = 1;
                        temp2[i][j] = 1;
                    }

                    if (j > i) {
                        a1[i][j] = -1;
                        a2[i][j] = -1;

                        temp1[i][j] = -1;
                        temp2[i][j] = -1;
                    }

                    if (j < i) {
                        a1[i][j] = 0;
                        a2[i][j] = 0;

                        temp1[i][j] = 1;
                        temp2[i][j] = 1;
                    }
                }
            }

            b[n - 1][0] = 1;

            temp1 = multiplyMatrixNumber(1e-3 * n, temp1);
            a1 = sumMatrixMatrix(a1, temp1);

            temp2 = multiplyMatrixNumber(1e-6 * n, temp2);
            a2 = sumMatrixMatrix(a2, temp2);

            makeTest2(sheet, a1, a2, b, currentRow, n, x);

            currentRow += 5;
        }

        for (int i = 0; i <= 14; i++) {
            sheet.autoSizeColumn(i);
        }
    }

    public static void main(String[] args) throws IOException {
        Workbook wb = new HSSFWorkbook();
        Sheet sheet = wb.createSheet("List №1");
        Sheet sheet1 = wb.createSheet("List №2");

        makeTitle(sheet.createRow(0), 1);

        makeFirstSheet(sheet);

        makeTitle(sheet1.createRow(0), 2);

        makeSecondSheet(sheet1);

        FileOutputStream fos = new FileOutputStream("20.Б01_Проценко_СЛАУ№2.xls");

        wb.write(fos);
        fos.close();
        wb.close();
    }
}
