#include <assert.h>
#include <iostream>
#include <math.h>

using namespace std;

template<typename T>
class Matrix;

template<typename T>
ostream &operator << (ostream &out, const Matrix<T> &m);

template<typename T>
Matrix<T> operator*
		(const Matrix<T> & m1, const Matrix<T> & m2);

template<typename T>
Matrix<T> operator*
		(T c, const Matrix<T> & m2);

template<typename T>		
Matrix<T> operator* 
		(const Matrix<T> & m1, T c);

template <typename T>
T det(T **b, int m);

template<typename T>
class Matrix {
public:
	Matrix(int sizeX = 1, int sizeY = 1);
	Matrix();
	~Matrix();
	Matrix(const Matrix<T>& m);
	Matrix& operator=(const Matrix<T>& rhs);
	Matrix operator+(const Matrix<T> & m);
	Matrix& operator+=(const Matrix<T> & m);
	Matrix	operator-(const Matrix<T>& m);
	Matrix&	operator-=(const Matrix<T>& m);
	friend ostream &operator<< <T>
		(ostream &out, const Matrix<T> &m);
	T &operator()(T x, T y);
	bool is_square();
	void detmnt();
	friend Matrix operator* <T>
		(const Matrix<T> & m1, const Matrix<T> & m2);
	friend Matrix operator* <T>
		(T c, const Matrix<T> & m2);               	
	friend Matrix operator* <T>
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
Matrix<T>::Matrix(int sizeX, int sizeY)
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
Matrix<T> &Matrix<T>::operator=(const Matrix<T>& m) {
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
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &m) {
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
Matrix<T> Matrix<T>::operator+(const Matrix<T> &m) {
	Matrix<T> temp(*this); //copy constructor
	return (temp += m);
}

template <typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &m){
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
Matrix<T> Matrix<T>::operator-(const Matrix<T> &m){
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
T &Matrix<T>::operator()(T i, T j) {
	assert(i>=0 && i<dx);
	assert(j>=0 && j<dy);
	return p[i][j];
}

template <typename T>
Matrix<T> operator* (const Matrix<T>& m1, const Matrix<T>& m2) {
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
Matrix<T> operator* (T c, const Matrix<T>& m2) {
	Matrix<T> prod(m2);
	for (int i=0; i<prod.dx; ++i) {
		for (int j=0; j<prod.dy; ++j) {
			prod.p[i][j] = c * m2.p[i][j];
		}
	}
	return prod;
}

template <typename T>
Matrix<T> operator* (const Matrix<T>& m2, T c) {
	return c*m2;
}

template <typename T>
void Matrix<T>::detmnt()
{

     determ = det(p, dx);
     cout<<"det is "<< determ << endl;
}

template <typename T>
T det(T **b,int m)
{
	cout << "this is something" << endl;
    int sum=0,x=0,y=0,i=0,j,aa,bb;
    T **c;
    c = new T*[m];
    cout << "this is something 2" << endl;
	for (int i = 0; i < m; i++)	{
		c[i] = new T[m]; 
		
	}
    cout << "this is something 3" << endl;
    if(m==2)
            return(b[0][0]*b[1][1]-b[0][1]*b[1][0]);
    else{
        for(j=0;j<m;j++){
			cout << "this is something 4" << endl;
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
