#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)  //创建矩阵
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)  //矩阵加法
{
    if (a.rows == b.rows && a.cols == b.cols){  //判断矩阵大小是否匹配
        Matrix result = create_matrix(a.rows, a.cols);
        for(int i = 0; i < a.rows; i++) {
            for(int j = 0; j < a.cols; j++) {
                result.data[i][j] = a.data[i][j] + b.data[i][j];  //矩阵加法
            }
        }
        return result;
    }else{
        printf("Error: Matrix a and b must have the same rows and cols.\n"); //报错提示
    }
    return create_matrix(0, 0);  //不可加返回零矩阵
}

Matrix sub_matrix(Matrix a, Matrix b)  //矩阵减法
{
    if (a.rows == b.rows && a.cols == b.cols) {  //判断矩阵大小是否匹配
        Matrix result = create_matrix(a.rows, a.cols);  //创建结果矩阵
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < a.cols; j++) {
                result.data[i][j] = a.data[i][j] - b.data[i][j];  //矩阵减法
            }
        }
        return result;
    }
    else{
        printf("Error: Matrix a and b must have the same rows and cols.\n");  //报错提示
    }
    return create_matrix(0, 0);
}

Matrix mul_matrix(Matrix a, Matrix b)  //矩阵乘法
{
    if(a.cols == b.rows){
        Matrix result = create_matrix(a.rows, b.cols);
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < b.cols; j++) {
                result.data[i][j] = 0;  //初始化结果为零
            }
        }
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < b.cols; j++) {
                for (int k = 0;k < a.cols; k++) {
                    result.data[i][j] += a.data[i][k] * b.data[k][j];  //矩阵乘法
                }
            }
        }
        return result;
    }else{
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");  //报错提示
    }
    return create_matrix(0, 0);
}

Matrix scale_matrix(Matrix a, double k)  //矩阵数乘
{
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            a.data[i][j] *= k;  
        }
    }
    return a;
    return create_matrix(0, 0);
}

Matrix transpose_matrix(Matrix a)  //矩阵转置
{
    Matrix result = create_matrix(a.cols, a.rows);  //创建转置矩阵
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[j][i] = a.data[i][j];  //交换对称元素
        }
    }
    return result;
    return create_matrix(0, 0);
}

double det_matrix(Matrix a)  //矩阵行列式
{
    if (a.rows == a.cols){ //判断是否为方阵
        if (a.rows == 1) {
            return a.data[0][0];  //1x1矩阵行列式为元素本身
        } else if (a.rows == 2) {
            return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];  //2x2矩阵行列式
        } else {  //大于2阶的方阵使用Laplace展开法计算行列式
            double det = 0;
            for (int i = 0; i < a.cols; i++) {
                Matrix sub_matrix = create_matrix(a.rows - 1, a.cols - 1);  //创建子矩阵
                for (int j = 1; j < a.rows; j++) {
                    for (int k = 0; k < a.cols; k++) {  
                        if (k < i) {
                            sub_matrix.data[j - 1][k] = a.data[j][k]; 
                        } else if (k > i) {
                            sub_matrix.data[j - 1][k - 1] = a.data[j][k];  //填充子矩阵
                        }
                    }
                }
                det += pow(-1, i) * a.data[0][i] * det_matrix(sub_matrix);  //递归计算行列式
            }
            return det;
        }
    }else{
        printf("Error: The matrix must be a square matrix.\n");  //报错提示
    }
    return 0;
}

Matrix inv_matrix(Matrix a)  //矩阵的逆
{
    if(a.rows == a.cols) {  //判断是否为方阵
        double det = det_matrix(a);  //计算行列式
        if (det != 0) {  //判断行列式是否为零
            Matrix result = create_matrix(a.rows, a.cols);  //创建逆矩阵
            for (int i = 0; i < a.rows; i++) {
                for (int j = 0; j < a.cols; j++) {
                    Matrix sub_matrix = create_matrix(a.rows - 1, a.cols - 1);  //创建子矩阵
                    for (int k = 0; k < a.rows; k++) {
                        for (int l = 0; l < a.cols; l++) {
                            if (k != i && l != j) {
                                sub_matrix.data[k - (k > i)][l - (l > j)] = a.data[k][l];  //填充子矩阵
                            }
                        }
                    }
                    result.data[j][i] = pow(-1, i + j) * det_matrix(sub_matrix) / det;  //计算逆矩阵元素
                }
            }
            return result;
        } else {
            printf("Error: The matrix is singular.\n");  //报错提示
        }
    } else {
        printf("Error: The matrix must be a square matrix.\n");  //报错提示
    }
    // 如果不是方阵或行列式为零，返回一个空矩阵
    return create_matrix(0, 0);
}

int rank_matrix(Matrix a)   //矩阵的秩
{
    int rank = a.rows < a.cols ? a.rows : a.cols;  // 初始秩为行数和列数的最小值
    for (int i = 0; i < rank; i++) {
        if (a.data[i][i] == 0) {
            int j;
            for (j = i + 1; j < a.rows; j++) {
                if (a.data[j][i] != 0) {  // 寻找非零行
                    for (int k = 0; k < a.cols; k++) {
                        double temp = a.data[i][k];
                        a.data[i][k] = a.data[j][k];
                        a.data[j][k] = temp;  // 将非零行交换到当前遍历行
                    }
                    break;
                }
            }
            if (j == a.rows) {  // 如果没找到非零行
                rank--;
                continue;
            }
        }
        // 将当前列下方元素消为零
        for (int j = i + 1; j < a.rows; j++) {
            double factor = a.data[j][i] / a.data[i][i];
            for (int k = i; k < a.cols; k++) {
                a.data[j][k] -= factor * a.data[i][k];
            }
        }
    }
    // 统计非零行的数量
    int non_zero_rows = 0;
    for (int i = 0; i < a.rows; i++) {
        int  is_zero_row = 1;
        for (int j = 0; j < a.cols; j++) {
            if (a.data[i][j] != 0) {
                is_zero_row = 0;
                break;
            }
        }
        if (!is_zero_row) {
            non_zero_rows++;
        }
    }
    return non_zero_rows;
}

double trace_matrix(Matrix a)  //矩阵的迹
{
    // 严格检查方阵条件
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    double trace = 0;
    for (int i = 0; i < a.rows; ++i) {
        // 添加对角线元素存在性检查
        if (a.data[i] == NULL) {
            printf("Error: Matrix row %d is NULL.\n", i);
            return 0.0;
        }
        trace += a.data[i][i];
    }

    return trace;
}

void print_matrix(Matrix a)  
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}