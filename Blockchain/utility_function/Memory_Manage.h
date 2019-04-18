#pragma once

void* Delete_1D_Array(void *ptr);
void* Delete_2D_Array(void **ptr, int row);

int** Allocate_2D_Array_Int(int row, int col, char *msg);
int*  Allocate_1D_Array_Int(int len, char *msg);

unsigned char** Allocate_2D_Array_UChar(int row, int col, char *msg);
unsigned char*  Allocate_1D_Array_UChar(int len, char *msg);

double** Allocate_2D_Array_Double(int row, int col, char *msg);
double*  Allocate_1D_Array_Double(int len, char *msg);

