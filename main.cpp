#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double eps = 0.000001;

// Fungsi f1 dan f2
double f1(double x, double y) { return x*x + x*y - 10; }
double f2(double x, double y) { return y + 3*x*y*y - 57; }

// Fungsi iterasi g1A dan g2B
double g1A(double x, double y) { return (10.0 - x*x) / y; }
double g2B(double x, double y) { return sqrt((57 - y) / (3*x)); }

// -------- Iterasi Titik Tetap Jacobi --------
void Iterasi_Jacobi(double x0, double y0) {
    cout << "\n=== Iterasi Titik Tetap - Jacobi (g1A, g2B) ===\n";
    int i = 0;
    double x1, y1, dx, dy;
    double errPrev = 1e9, err, ratio = 0;

    cout << left << setw(10) << "Iterasi"
         << setw(12) << "x"
         << setw(12) << "y"
         << setw(12) << "delta x"
         << setw(12) << "delta y"
         << setw(12) << "Error"
         << "Rasio\n";

    do {
        x1 = g1A(x0, y0);
        y1 = g2B(x0, y0);
        dx = fabs(x1 - x0);
        dy = fabs(y1 - y0);
        err = sqrt(dx*dx + dy*dy);
        ratio = err / errPrev;
        errPrev = err;
        i++;

        cout << setw(10) << i
             << setw(12) << x1
             << setw(12) << y1
             << setw(12) << dx
             << setw(12) << dy
             << setw(12) << err
             << ratio << endl;

        if (err < eps) break;
        x0 = x1; y0 = y1;
    } while (i < 50);

    cout << "\nHasil akhir Jacobi: x=" << x1 << ", y=" << y1 << endl;
    if (isnan(x1) || isnan(y1)) cout << "Iterasi divergen (NaN terdeteksi)\n";
    else if (ratio < 1) cout << "Konvergen (rasio ~" << ratio << ")\n";
    else cout << "Cenderung divergen (rasio ≥ 1)\n";
}

// -------- Iterasi Titik Tetap Seidel --------
void Iterasi_Seidel(double x0, double y0) {
    cout << "\n=== Iterasi Titik Tetap - Seidel (g1A, g2B) ===\n";
    int iter = 0;
    double x1, y1, dx, dy;
    double errPrev = 1e9, err, ratio = 0;

    cout << left << setw(10) << "Iterasi"
         << setw(12) << "x"
         << setw(12) << "y"
         << setw(12) << "delta x"
         << setw(12) << "delta y"
         << setw(12) << "Error"
         << "Rasio\n";

    do {
        x1 = g1A(x0, y0);
        y1 = g2B(x1, y0);
        dx = fabs(x1 - x0);
        dy = fabs(y1 - y0);
        err = sqrt(dx*dx + dy*dy);
        ratio = err / errPrev;
        errPrev = err;
        iter++;

        cout << setw(10) << iter
             << setw(12) << x1
             << setw(12) << y1
             << setw(12) << dx
             << setw(12) << dy
             << setw(12) << err
             << ratio << endl;

        if (err < eps) break;
        x0 = x1; y0 = y1;
    } while (iter < 50);

    cout << "\nHasil akhir Seidel: x=" << x1 << ", y=" << y1 << endl;
    if (isnan(x1) || isnan(y1)) cout << "Iterasi divergen (NaN terdeteksi)\n";
    else if (ratio < 1) cout << "Konvergen (rasio ~" << ratio << ")\n";
    else cout << "Cenderung divergen (rasio ≥ 1)\n";
}

// -------- Metode Newton-Raphson --------
void Newton_Raphson(double x0, double y0) {
    cout << "\n=== Metode Newton-Raphson ===\n";
    int iter = 0;
    double x1, y1, dx, dy;
    double errPrev = 1e9, err, ratio = 0;

    cout << left << setw(10) << "Iterasi"
         << setw(12) << "x"
         << setw(12) << "y"
         << setw(12) << "delta x"
         << setw(12) << "delta y"
         << setw(12) << "Error"
         << "Rasio\n";

    do {
        double J = (2*x0 + y0)*(1 + 6*x0*y0) - (x0)*(3*y0*y0);
        double u = f1(x0, y0);
        double v = f2(x0, y0);
        double dfdx1 = 2*x0 + y0;
        double dfdy1 = x0;
        double dfdx2 = 3*y0*y0;
        double dfdy2 = 1 + 6*x0*y0;

        x1 = x0 - (u*dfdy2 - v*dfdy1)/J;
        y1 = y0 + (u*dfdx2 - v*dfdx1)/J;

        dx = fabs(x1 - x0);
        dy = fabs(y1 - y0);
        err = sqrt(dx*dx + dy*dy);
        ratio = err / errPrev;
        errPrev = err;
        iter++;

        cout << setw(10) << iter
             << setw(12) << x1
             << setw(12) << y1
             << setw(12) << dx
             << setw(12) << dy
             << setw(12) << err
             << ratio << endl;

        if (err < eps) break;
        x0 = x1; y0 = y1;
    } while (iter < 50);

    cout << "\nHasil akhir Newton-Raphson: x=" << x1 << ", y=" << y1 << endl;
    if (ratio < 0.5) cout << "Konvergensi cepat (quadratic)\n";
    else cout << "Konvergensi lambat atau tidak stabil\n";
}

// -------- Metode Secant --------
void Secant(double x0, double y0, double x1, double y1) {
    cout << "\n=== Metode Secant ===\n";
    int iter = 0;
    double xn, yn, dx, dy;
    double errPrev = 1e9, err, ratio = 0;

    cout << left << setw(10) << "Iterasi"
         << setw(12) << "x"
         << setw(12) << "y"
         << setw(12) << "delta x"
         << setw(12) << "delta y"
         << setw(12) << "Error"
         << "Rasio\n";

    do {
        double f1x0 = f1(x0, y0), f1x1 = f1(x1, y1);
        double f2x0 = f2(x0, y0), f2x1 = f2(x1, y1);
        xn = x1 - f1x1*(x1-x0)/(f1x1 - f1x0);
        yn = y1 - f2x1*(y1-y0)/(f2x1 - f2x0);

        dx = fabs(xn - x1);
        dy = fabs(yn - y1);
        err = sqrt(dx*dx + dy*dy);
        ratio = err / errPrev;
        errPrev = err;
        iter++;

        cout << setw(10) << iter
             << setw(12) << xn
             << setw(12) << yn
             << setw(12) << dx
             << setw(12) << dy
             << setw(12) << err
             << ratio << endl;

        if (err < eps) break;
        x0 = x1; y0 = y1;
        x1 = xn; y1 = yn;
    } while (iter < 50);

    cout << "\nHasil akhir Secant: x=" << xn << ", y=" << yn << endl;
    if (ratio < 1) cout << "Konvergen moderat (lebih lambat dari Newton)\n";
    else cout << "Cenderung divergen\n";
}

// -------- Main --------
int main() {
    double x0 = 1.5, y0 = 3.5;
    Iterasi_Jacobi(x0, y0);
    Iterasi_Seidel(x0, y0);
    Newton_Raphson(x0, y0);
    Secant(1.0, 3.0, 1.5, 3.5);
    return 0;
}
