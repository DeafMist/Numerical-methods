import java.util.Arrays;
import java.util.Random;

public final class Zad3 {
    public static double calcSin(double x) {
        int n = 1;
        double xn = x;
        double cur = x;
        double factorial = 1;

        while (Math.abs(cur) > 1e-6) {
            n += 2;
            factorial *= -n * (n - 1);
            cur = Math.pow(x, n) / factorial;
            xn += cur;
        }

        return xn;
    }

    private static double func1(double x) {
        return x * x + 4 * Math.sin(x) - 2;
    }

    private static double[] func2(double x, double y) {
        return new double[]{
                Math.sin(y - 1) + x - 1.3,
                y - Math.sin(x + 1) - 0.8
        };
    }

    private static double func1Derivative(double x) {
        return 2 * x + 4 * Math.cos(x);
    }

    private static double[][] func2Jacobi(double x, double y) {
        return new double[][]{
                {1, Math.cos(y - 1)},
                {-Math.cos(x + 1), 1}
        };
    }

    private static double[][] func2JacobiDerivative(double x, double y) {
        return new double[][] {
                {0, -Math.sin(y - 1)},
                {Math.sin(x + 1), 0}
        };
    }

    private static double pickOfNextX(double a, double b) {
        return (((new Random()).nextInt(2) + 1) == 1 ? (a + b) / 2 : a - (b - a) * func1(a) / (func1(b) - func1(a)));
    }

    private static double[][] rightCol(double x, double y) {
        return new double[][]{
                {-calcSin(y - 1) - x + 1.3},
                {-y + calcSin(x + 1) + 0.8}
        };
    }

    private static double methodNewton(double a, double b, double epsi) {
        double next = (a + b) / 2;
        double prev;

        do {
            prev = next;
            next = prev - func1(prev) / func1Derivative(prev);

            if (next < a || next > b){
                next = pickOfNextX(a, b);
            }

            if (func1(a) > 0 && func1(next) > 0 || func1(a) < 0 && func1(next) < 0) {
                a = next;
            }
            else {
                b = next;
            }
        } while (Math.abs(next - prev) > epsi);

        return next;
    }

    private static double[][] methodNewton(double[][] approxX, double epsi) {
        double[][] prev = new double[2][1];
        double[][] next = approxX;
        double[][] gh;

        do {
            prev[0][0] = next[0][0];
            prev[1][0] = next[1][0];

            gh = Zad2.methodHouseholderQR(func2Jacobi(prev[0][0], prev[1][0]), rightCol(prev[0][0], prev[1][0]));

            next[0][0] = prev[0][0] + gh[0][0];
            next[1][0] = prev[1][0] + gh[1][0];
        } while (Zad2.matrixNorm1(Zad2.sumMatrixMatrix(next, Zad2.multiplyMatrixNumber(-1, prev))) > epsi);

        return next;
    }

    private static double[] seqEnum(double a, double b) {
        int N = 10;
        double xk = a;
        double xk1 = a;

        while (func1(xk) * func1(xk1) > 0) {
            xk1 = a;
            double h = (b - xk) / N;

            for (int i = 1; i < N; i++) {
                xk = xk1;
                xk1 = a + i * h;

                if (func1(xk) * func1(xk1) < 0)
                    break;
            }

            N *= 2;
        }

        return new double[] { xk, xk1 };
    }

    private static double[][] helperFunc(double l, double x, double y) {
        return new double[][]{
                {l * Math.sin(y - 1) + x - 1.3},
                {y - l * Math.sin(x + 1) - 0.8}
        };
    }

    private static double[][] approximation() {
        double N = 10.0;
        int n = 2;
        double dx = 1e-6;

        double[][] xk = {{1}, {0.7}};

        for (int i = 0; i < N; i++) {
            double[][] der = new double[n][n];
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    double[][] xh = Zad2.copy(xk);
                    xh[j][0] += dx;
                    der[k][j] = (helperFunc(i / N, xh[0][0], xh[1][0])[k][0] - helperFunc(i / N, xk[0][0], xk[1][0])[k][0]) / dx;
                }
            }

            double[][] b = Zad2.multiplyMatrixNumber(-1, helperFunc(i / N, xk[0][0], xk[1][0]));
            double[][] deltaX = Zad2.methodGaussLUP(der, b);
            xk = Zad2.sumMatrixMatrix(xk, deltaX);
        }

        return xk;
    }

    public static void main(String[] args) {
        double a1 = -100;
        double b1 = 100;

        double[] temp = seqEnum(a1, b1);
        a1 = temp[0];
        b1 = temp[1];

        double x1 = methodNewton(a1, b1, 1e-4);
        System.out.println(
                "------------First function------------" + '\n' +
                "x = " + x1 + '\n' +
                "f(x) = " + func1(x1) + '\n'
        );

        double[][] a2 = {{0.5}, {1.5}};

        double[][] x2 = methodNewton(a2, 1e-4);
        System.out.println(
                "------------System, graphical approximation------------" + '\n' +
                        "x = [" + x2[0][0] + ", " + x2[1][0] + "]" + '\n' +
                        "f(x) = [" + func2(x2[0][0], x2[1][0])[0] + ", " + func2(x2[0][0], x2[1][0])[1] + "]" + '\n'
        );

        a2 = approximation();
        x2 = methodNewton(a2, 1e-4);
        System.out.println(
                "------------System, approximation with helper function------------" + '\n' +
                        "x = [" + x2[0][0] + ", " + x2[1][0] + "]" + '\n' +
                        "f(x) = [" + func2(x2[0][0], x2[1][0])[0] + ", " + func2(x2[0][0], x2[1][0])[1] + "]" + '\n'
        );
    }
}
