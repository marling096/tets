#include <iostream>  
#include <fstream>
#include "stdio.h"
#include "string.h"
#include "math.h"
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#define MATRIX_SIZE 6
#define u8 unsigned char
double m_matrix[MATRIX_SIZE][MATRIX_SIZE + 1];// 系数矩阵
double solve[MATRIX_SIZE] = { 0 };  // 方程组的解对应最小二乘椭球拟合中的，a，b，c，d，e，f， 

double m_result[MATRIX_SIZE];
int N = 0;// 计算累计的采样点次数的
double X0 = 0, Y0 = 0, Z0 = 0, A = 0, B = 0, C = 0;

// 把矩阵系数全部清除为0
void ResetMatrix(void)
{
	for (u8 row = 0; row < MATRIX_SIZE; row++)
	{
		for (u8 column = 0; column < MATRIX_SIZE + 1; column++)
			m_matrix[row][column] = 0.0f;
	}
}

// 把输入的数据先生成矩阵的元素的总和
void CalcData_Input(double x, double y, double z)
{
	double V[MATRIX_SIZE + 1];
	N++;
	V[0] = y*y;
	V[1] = z*z;
	V[2] = x;
	V[3] = y;
	V[4] = z;
	V[5] = 1.0;
	V[6] = -x*x;
	// 构建系数矩阵，并进行累加
	for (u8 row = 0; row < MATRIX_SIZE; row++)
	{
		for (u8 column = 0; column < MATRIX_SIZE + 1; column++)
		{
			m_matrix[row][column] += V[row] * V[column];
		}
	}
 
}

// 化简系数矩阵，把除以N带上  
void CalcData_Input_average()
{
	for (u8 row = 0; row < MATRIX_SIZE; row++)
	for (u8 column = 0; column < MATRIX_SIZE + 1; column++)
		m_matrix[row][column] /= N;

}

// 显示出来系数矩阵和增广矩阵[A|b]
void DispMatrix(void)
{
	for (u8 row = 0; row < MATRIX_SIZE; row++)
	{
		for (u8 column = 0; column < MATRIX_SIZE + 1; column++)
		{
			if (column == MATRIX_SIZE - 1){
				column=column;
			}

		}

	}

}

// 交换两行元素位置
void Row2_swop_Row1(int row1, int row2)
{
	double tmp = 0;
	for (u8 column = 0; column < MATRIX_SIZE + 1; column++)
	{
		tmp = m_matrix[row1][column];
		m_matrix[row1][column] = m_matrix[row2][column];
		m_matrix[row2][column] = tmp;
	}
}

// 用把row行的元素乘以一个系数k
void k_muiltply_Row(double k, int row)
{
	for (u8 column = 0; column < MATRIX_SIZE + 1; column++)
		m_matrix[row][column] *= k;
}

// 用一个数乘以row1行加到row2行上去
void Row2_add_kRow1(double k, int row1, int row2)
{
	for (u8 column = 0; column < MATRIX_SIZE + 1; column++)
		m_matrix[row2][column] += k*m_matrix[row1][column];
}


// 列主元，第k次消元之前，把k行到MATRIX_SIZE的所有行进行冒泡排出最大，排序的依据是k列的元素的大小
void MoveBiggestElement_to_Top(int k)
{
	int row = 0, column = 0;

	for (row = k + 1; row < MATRIX_SIZE; row++)
	{
		if (fabs(m_matrix[k][k]) < fabs(m_matrix[row][k]))
		{
			Row2_swop_Row1(k, row);
		}
	}
}

// 高斯消元法，求行阶梯型矩阵
u8 Matrix_GaussElimination(void)
{
	double k = 0;
	for (u8 cnt = 0; cnt < MATRIX_SIZE; cnt++)
	{
		
		MoveBiggestElement_to_Top(cnt);
		if (m_matrix[cnt][cnt] == 0)
			return(1);//error
	
		for (u8 row = cnt + 1; row < MATRIX_SIZE; row++)
		{
			k = -m_matrix[row][cnt] / m_matrix[cnt][cnt];
			Row2_add_kRow1(k, cnt, row);
		}
		DispMatrix();
	}
	return 0;
}

// 求行最简型矩阵，即把对角线的元素全部化成1
void Matrix_RowSimplify(void)
{
	double k = 0;
	for (u8 row = 0; row < MATRIX_SIZE; row++)
	{
		k = 1 / m_matrix[row][row];
		k_muiltply_Row(k, row);
	}
	DispMatrix();
}

// 求非齐次线性方程组的解
void Matrix_Solve(double* solve)
{
	for (short row = MATRIX_SIZE - 1; row >= 0; row--)
	{
		solve[row] = m_matrix[row][MATRIX_SIZE];
		for (u8 column = MATRIX_SIZE - 1; column >= row + 1; column--)
			solve[row] -= m_matrix[row][column] * solve[column];
	}

}

void Ellipsoid_fitting_Process(void)
{
    double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
	//ResetMatrix();

	CalcData_Input_average();
	DispMatrix();
	if (Matrix_GaussElimination())	{
		////printf("the marix could not be solved\r\n");
	}
		
	else
	{
		Matrix_RowSimplify();
		Matrix_Solve(solve);// 求解a,b,c,d,e,f

		
		a = solve[0];
		b = solve[1];
		c = solve[2];
		d = solve[3];
		e = solve[4];
		f = solve[5];

		
		X0 = -c / 2;
		Y0 = -d / (2 * a);
		Z0 = -e / (2 * b);
		A = sqrt(X0*X0 + a*Y0*Y0 + b*Z0*Z0 - f);
		B = A / sqrt(a);
		C = A / sqrt(b);

	}

}
void mag_correct(){
	Ellipsoid_fitting_Process();	
	cout<<"X0 " <<X0<<endl;
	cout<<"Y0 " <<Y0<<endl;
	cout<<"Z0 " <<Z0<<endl;
	cout<<"A " <<A<<endl;
	cout<<"B " <<B<<endl;
	cout<<"C " <<C<<endl;

}


int main() {
    std::ifstream file("data.txt");
    std::string line;
	int cnt=0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double num[3];

        while (iss >> num[cnt]) {
			
			if(cnt==2){
				CalcData_Input(num[0],num[1],num[2]);
				cnt=0;
			}
			else{
				cnt++;
			}			
        }
    }
	mag_correct();

	file.close();
    return 0;
}

	//




