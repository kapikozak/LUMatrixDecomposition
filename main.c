#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void exchange(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void swapRows(double** matrix, int i, int j) {
    double* temp = matrix[i];
    matrix[i] = matrix[j];
    matrix[j] = temp;
}

double** cloneMatrix(double** matrix, int matrixSize) {
    double** clone = (double**)malloc(matrixSize * sizeof(double*));
    for (int i = 0; i < matrixSize; i++) {
        clone[i] = (double*)malloc(matrixSize * sizeof(double));
    }
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            clone[i][j] = matrix[i][j];
        }
    }
    return clone;
}

void printMatrix(double** matrix, int matrixSize) {
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            printf("%.2lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void printVector(int* P, int matrixSize) {
    for (int i = 0; i < matrixSize; i++) {
        printf("%d ", P[i]);
    }
    printf("\n");
}

int singular;

void LUDecomposition(double** matrix, int matrixSize, double*** L, double*** U, int** P) {
    double** clone = cloneMatrix(matrix, matrixSize);
    *L = (double**)malloc(matrixSize * sizeof(double*));
    *U = (double**)malloc(matrixSize * sizeof(double*));
    *P = (int*)malloc(matrixSize * sizeof(int));

    for (int i = 0; i < matrixSize; i++) {
        (*L)[i] = (double*)calloc(matrixSize, sizeof(double));
        (*U)[i] = (double*)calloc(matrixSize, sizeof(double));
        (*P)[i] = i;
    }
    singular = 0;

    for (int k = 0; k < matrixSize; k++) {
        double maxDigit = fabs(clone[k][k]);
        int maxIndex = k;

        for (int i = k + 1; i < matrixSize; i++) {
            if (fabs(clone[i][k]) > maxDigit) {
                maxDigit = fabs(clone[i][k]);
                maxIndex = i;
            }
        }

        if (maxDigit == 0.0) {
            printf("Given Matrix is singular, LU algorithm cannot be executed\n");
            singular = 1;
            for (int i = 0; i < matrixSize; i++) {
                free(clone[i]);
            }
            free(clone);
            return;
        }

        if (maxIndex != k) {
            swapRows(clone, k, maxIndex);
            exchange(&(*P)[k], &(*P)[maxIndex]);
        }

        for (int i = k; i < matrixSize; i++) {
            (*U)[k][i] = clone[k][i];
            (*L)[i][k] = clone[i][k] / (*U)[k][k];
        }

        for (int i = k + 1; i < matrixSize; i++) {
            for (int j = k + 1; j < matrixSize; j++) {
                clone[i][j] -= (*L)[i][k] * (*U)[k][j];
            }
        }
    }
    for (int i = 0; i < matrixSize; i++) {
        free(clone[i]);
    }
    free(clone);
}

int main() {
    int matrixSize = 3;
    double** matrix = (double**)malloc(matrixSize * sizeof(double*));
    for (int i = 0; i  < matrixSize; i++) {
        matrix[i] = (double*)malloc(matrixSize * sizeof(double));
    }

    int values[3][3] = {
            {2, 2, -1},
            {1, 1, 0},
            {1, 2, 0}
    };

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            matrix[i][j] = values[i][j];
        }
    }

    double** L, **U;
    int* P;

    printf("Matrix before decomposition:\n");
    printMatrix(matrix, matrixSize);

    LUDecomposition(matrix, matrixSize, &L, &U, &P);
    if (!singular) {
        printf("\nMatrix L after decomposition:\n");
        printMatrix(L, matrixSize);

        printf("\nMatrix U after decomposition:\n");
        printMatrix(U, matrixSize);

        printf("\nPermutation vector P:\n");
        printVector(P, matrixSize);
    }

    for (int i = 0; i < matrixSize; i++) {
        free(matrix[i]);
        free(L[i]);
        free(U[i]);
    }
    free(matrix);
    free(L);
    free(U);
    free(P);


    return 0;
}