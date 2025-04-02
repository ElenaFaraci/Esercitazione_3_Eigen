#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;

double err_r(Vector2d x, Vector2d y){
	
	double err = (x - y).norm()/x.norm();
	return err;
	}

void decomposizione(Matrix2d A, Vector2d b, Vector2d x){
	
	cout << "Soluzione del sistema lineare Ax = b\n" << endl;
	cout << "Matrice A:\n" << std::setprecision(5) << std::scientific << A << "\n" << endl;
	cout << "Vettore b:\n" << b << "\n" << endl;
	cout << "Soluzione esatta:\nx = " << x.transpose() << "\n" << endl;
	
	// soluzione con fattorizzazione QR:
	
	Vector2d x_QR = A.colPivHouseholderQr().solve(b);
	
	cout << "Soluzione con decomposizione QR:\nx = " << x_QR.transpose() << endl;
	
	double err_QR = err_r(x, x_QR);
	
	cout << "Errore relativo dec QR:\nerr = " << err_QR << "\n" << endl;
	
	// soluzione con fattorizzazione PALU:
	
	Vector2d x_LU = A.fullPivLu().solve(b);
	
	cout << "Soluzione con decomposizione PALU:\nx = " << x_LU.transpose() << endl;
	
	double err_LU = err_r(x, x_LU);
	
	cout << "Errore relativo dec PALU:\nerr = " << err_LU << "\n" << endl;
	}


int main()
{
	
	// 1:
	
	cout << "- Risoluzione del PRIMO sistema:\n" << endl;
	
	Matrix2d A1{
		{5.547001962252291e-01, -3.770900990025203e-02},
		{8.320502943378437e-01, -9.992887623566787e-01}
		};
	Vector2d b1 {-5.169911863249772e-01, 1.672384680188350e-01};
	Vector2d x_esatto {-1.0e+0, -1.0e+00};
	
	decomposizione(A1, b1, x_esatto);
	
	// 2:
	
	cout << "- Risoluzione del SECONDO sistema:\n" << endl;
	
	Matrix2d A2{
		{5.547001962252291e-01, -5.540607316466765e-01},
		{8.320502943378437e-01, -8.324762492991313e-01}
		};
	Vector2d b2 {-6.394645785530173e-04, 4.259549612877223e-04};
	
	decomposizione(A2, b2, x_esatto);
	
	// 3:
	
	cout << "- Risoluzione del TERZO sistema:\n" << endl;
	
	Matrix2d A3{
		{5.547001962252291e-01, -5.547001955851905e-01},
		{8.320502943378437e-01, -8.320502947645361e-01}
		};
	Vector2d b3 {-6.400391328043042e-10, 4.266924591433963e-10};
	
	decomposizione(A3, b3, x_esatto);
	
	return 0;
}
