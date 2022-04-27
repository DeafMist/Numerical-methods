import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import java.io.FileOutputStream;
import java.io.IOException;

public final class Zad1 {
    private static final double deltaU = 5.41e-7;
    private static final double deltaV = 5.0013e-7;
    private static final double deltaPhi = 4.0523e-5;

    private static double calcSqrt(double x) {
        double pi = 1;
        double next = (pi + x / pi) / 2;

        while (Math.abs(pi - next) > deltaPhi) {
            next = pi;
            pi = (pi + x / pi) / 2;
        }

        return pi;
    }

    private static double calcCos(double x) {
        int n = 0;
        double xn = 1;
        double cur = xn;
        double factorial = 1;

        while (Math.abs(cur) > deltaU) {
            n += 2;
            factorial *= -n * (n - 1);
            cur = Math.pow(x, n) / factorial;
            xn += cur;
        }

        return xn;
    }

    private static int sgn(double x) {
        if (x > 0)
            return 1;
        else if (x < 0)
            return -1;
        else
            return 0;
    }

    private static double calcArcTan(double x) {
        int n = 0;
        double xn;
        double sign;

        if (Math.abs(x) < 1) {
            xn = x;
            sign = 1;
        }
        else {
            xn = sgn(x) * Math.PI / 2 - x;
            sign = -1;
        }

        double cur = xn;

        while (Math.abs(cur) > deltaV) {
            n += 1;
            cur = sign * Math.pow(-1, n) * Math.pow(x, sign * (2 * n + 1)) / (2 * n + 1);
            xn += cur;
        }

        return xn;
    }

    private static void makeTitle(Row row) {
        Cell cell;

        cell = row.createCell(0);
        cell.setCellValue("x");

        cell = row.createCell(1);
        cell.setCellValue("phi(x)");

        cell = row.createCell(2);
        cell.setCellValue("deltaPhi");

        cell = row.createCell(3);
        cell.setCellValue("phi_(x)");

        cell = row.createCell(4);
        cell.setCellValue("deltaPhi_");

        cell = row.createCell(5);
        cell.setCellValue("u(x)");

        cell = row.createCell(6);
        cell.setCellValue("deltaU");

        cell = row.createCell(7);
        cell.setCellValue("u_(x)");

        cell = row.createCell(8);
        cell.setCellValue("deltaU_");

        cell = row.createCell(9);
        cell.setCellValue("v(x)");

        cell = row.createCell(10);
        cell.setCellValue("deltaV");

        cell = row.createCell(11);
        cell.setCellValue("v_(x)");

        cell = row.createCell(12);
        cell.setCellValue("deltaV_");

        cell = row.createCell(13);
        cell.setCellValue("z(x)");

        cell = row.createCell(14);
        cell.setCellValue("deltaZ");

        cell = row.createCell(15);
        cell.setCellValue("z_(x)");

        cell = row.createCell(16);
        cell.setCellValue("deltaZ_");
    }

    private static void makeTest(Sheet sheet) {
        double x = 0.1;

        for (int row_i = 1; row_i < 12; row_i++, x += 0.01) {
            Row row = sheet.createRow(row_i);

            double myPhi = 2.8 * x + calcSqrt(1 + x);
            double phi = 2.8 * x + Math.sqrt(1 + x);

            double myU = calcCos(myPhi);
            double u = Math.cos(phi);

            double myV = calcArcTan(1.5 * x + 0.2);
            double v = Math.atan(1.5 * x + 0.2);

            double myZ = myU * myV;
            double z = u * v;

            for (int cell_i = 0; cell_i < 17; cell_i++) {
                Cell cell = row.createCell(cell_i);

                switch (cell_i) {
                    case 0 -> cell.setCellValue(x);
                    case 1 -> cell.setCellValue(myPhi);
                    case 2 -> cell.setCellValue(deltaPhi);
                    case 3 -> cell.setCellValue(phi);
                    case 4 -> cell.setCellValue(Math.abs(phi - myPhi));
                    case 5 -> cell.setCellValue(myU);
                    case 6 -> cell.setCellValue(deltaU);
                    case 7 -> cell.setCellValue(u);
                    case 8 -> cell.setCellValue(Math.abs(u - myU));
                    case 9 -> cell.setCellValue(myV);
                    case 10 -> cell.setCellValue(deltaV);
                    case 11 -> cell.setCellValue(v);
                    case 12 -> cell.setCellValue(Math.abs(v - myV));
                    case 13 -> cell.setCellValue(myZ);
                    case 14 -> cell.setCellValue(1e-6);
                    case 15 -> cell.setCellValue(z);
                    case 16 -> cell.setCellValue(Math.abs(z - myZ));
                }
            }
        }
    }

    public static void main(String[] args) throws IOException {
        Workbook wb = new HSSFWorkbook();
        Sheet sheet = wb.createSheet();

        makeTitle(sheet.createRow(0));

        makeTest(sheet);

        FileOutputStream fos = new FileOutputStream("20.Б01_Проценко_Погрешности№1.xls");

        wb.write(fos);
        fos.close();
        wb.close();
    }
}
