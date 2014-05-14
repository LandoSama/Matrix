#include <stdio.h>
#include "Matrix.hpp"

int main() {
	Matrix x(2,2), y(2,2), z(1,1), b(4,4);
	x(0,0) = 1; x(0,1) = 1; x(1,0) = 1; x(1,1) = 1;
	y(0,0) = 3; y(0,1) = 4; y(1,0) = 1; y(1,1) = 1;
	b(0,0) = 2; b(1,0) = 3; b(2,2) = 9; b(2,1) = 5; b(1,1) = 1; b(3,3) = 3;
	cout << "Matrix x\n" << x << "\nMatrix y\n" << y << "\nMatrix b\n" << b << endl;
	cout << "x*y = \n" << x*y << endl;
	x.is_square();
	b.detmnt();
	z = x*y;
	cout << "Matrix z = x*y  (note new dimensions)\n" << z << endl;
	Matrix x2(2,1);
	x2 = 2*x;
	cout << "Matrix x2 = 2*x\n" << x2 << endl;
	cout << "x + x2 = \n" << x+x2 << endl;
	return 0;
}
