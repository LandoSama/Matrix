#include <assert.h>
#include <iostream>
#include <math.h>

using namespace std;

template<typename T>
class Matrix {
public:
	Matrix(int sizeX, int sizeY);
	Matrix();
	~Matrix();
	Matrix(const Matrix<T>& m);
	Matrix& operator=(const Matrix<T>& rhs);
	Matrix operator+(const Matrix<T> & m);
	Matrix& operator+=(const Matrix<T> & m);
	Matrix	operator-(const Matrix<T>& m);
	Matrix&	operator-=(const Matrix<T>& m);
	friend ostream &operator<<
		(ostream &out, const Matrix<T> &m);
	T &operator()(T x, T y);
	bool is_square();
	T detmnt();
	friend Matrix operator*
		(const Matrix<T> & m1, const Matrix<T> & m2);
	friend Matrix operator*
		(T c, const Matrix<T> & m2);
	friend Matrix operator*
		(const Matrix<T> & m1, T c);

private:
	int dx, dy;  // dimensions, dx by dy and determinant 
	T determ;
	T **p;	// pointer to a pointer to a long integer
	void allocArrays() {
		assert(dx>0);
		assert(dy>0);
		p = new T*[dx];
		for (int i = 0; i < dx; i++)	{
			p[i] = new T[dy]; 
		}
	}
};

template <typename T>
Matrix::Matrix(int sizeX=1, int sizeY=1)
: dx(sizeX),dy(sizeY)  {
	allocArrays();
	for (int i = 0; i < dx; i++)	{
		for (int j = 0; j < dy; j++) {
			p[i][j] = 0;
		}
	}
}

template <typename T>
bool Matrix<T>::is_square(){
		if(dx == dy){
			cout << "is square" << endl;
			return true;
		}
		else{
			cout << "not square" << endl;
			return false;
		}
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& m) : dx(m.dx), dy(m.dy) {
	allocArrays();
	for (int i=0; i<dx; ++i) {
		for (int j=0; j<dy; ++j) {
			p[i][j] = m.p[i][j];
		}
	}
}

template <typename T>
Matrix<T>::~Matrix() {
	for (int i = 0; i < dx; i++) {
		delete [] p[i]; 
	}
	delete [] p;
	p = 0;
}

template <typename T>
Matrix<T> &Matrix::operator=(const Matrix<T>& m) {
	if (this == &m) { 
		// avoid self-assignment
		return *this;
	} else {
		if (dx != m.dx || dy != m.dy) {
			this->~Matrix();
			dx = m.dx; dy = m.dy;
			allocArrays();
		}
		for (int i = 0; i < dx; i++) {
			for (int j = 0; j < dy; j++) {
				p[i][j] = m.p[i][j];
			}
		}
		return *this;
	}
}

template <typename T>
Matrix<T> &Matrix::operator+=(const Matrix<T> &m) {
	assert(dx==m.dx);
	assert(dy==m.dy);
	// x+=y adds the y-entries into the x-entries
	for (int i=0; i<dx; ++i) {
		for (int j=0; j<dy; ++j) {
			p[i][j] += m.p[i][j];
		}
	}
	return *this;
}

template <typename T>
Matrix<T> Matrix::operator+(const Matrix<T> &m) {
	Matrix<T> temp(*this); //copy constructor
	return (temp += m);
}

template <typename T>
Matrix<T> &Matrix::operator-=(const Matrix<T> &m){
	assert(dx==m.dx);
	assert(dy==m.dy);
	// x+=y adds the y-entries into the x-entires
	for(int i=0; i<dx; ++i){
		for(int j=0; j<dy; ++j){
			p[i][j] -= m.p[i][j];
		}
	}
	return *this;
}

template <typename T>
Matrix<T> Matrix::operator-(const Matrix<T> &m){
	Matrix<T> temp(*this); //copy constructor
	return(temp -= m);
}

template <typename T>
ostream &operator<<
(ostream &out, const Matrix<T> &m) 
{
	for (int i = 0; i < m.dx; ++i)	{
		for (int j = 0; j < m.dy; ++j)
			out << m.p[i][j] << "  ";
		out << endl;
	}
	return out;
}

template <typename T> 
T &Matrix::operator()(T i, T j) {
	assert(i>=0 && i<dx);
	assert(j>=0 && j<dy);
	return p[i][j];
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2) {
	Matrix<T> prod(m1.dx, m2.dy);
	for (int i=0; i<prod.dx; ++i) {
		for (int j=0; j<prod.dy; ++j) {
			for (int k=0; k<m1.dy; ++k) {
				prod.p[i][j] += m1.p[i][k] * m2.p[k][j];
			}
		}
	};
	return prod;
}

template <typename T>
Matrix<T> operator*(T c, const Matrix<T>& m2) {
	Matrix<T> prod(m2);
	for (int i=0; i<prod.dx; ++i) {
		for (int j=0; j<prod.dy; ++j) {
			prod.p[i][j] = c * m2.p[i][j];
		}
	}
	return prod;
}

template <typename T>
Matrix<T> operator*(const Matrix& m2, T c) {
	return c*m2;
}

template <typename T>
void Matrix::detmnt()
{
     T det(T **b,int m);
     cout<<"det is "<<det(p,dx) << endl;
}
 
template <typename T>
T det(T **b,int m)
{
    int sum=0,x=0,y=0,i=0,j,aa,bb;
    T **c;
    c = new T*[m];
	for (int i = 0; i < m; i++)	{
		c[i] = new T[m]; 
	}
    
    if(m==2)
            return(b[0][0]*b[1][1]-b[0][1]*b[1][0]);
    else{
        for(j=0;j<m;j++){
			for(aa=0,y=0;aa<m;aa++){
				for(bb=0;bb<m;bb++){
					if((aa!=i)&&(bb!=j)){
						c[x][y++]=b[aa][bb];
					}
                if(y>0)x++;
				y=0;
             }
             sum=sum+b[i][j]*pow(-1,i+j)*det(c,m-1);
			}
		}
	}
    return sum;
}


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
