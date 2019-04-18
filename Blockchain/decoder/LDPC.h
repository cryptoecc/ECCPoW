#pragma once
class LDPC
{
private:
	int row_deg, col_deg, block_length, message_length, seed, q;
	int redundancy, lastnonzerow = 0, Iter = 0;
	double code_rate;

	unsigned char **H;
	unsigned char **H_ORI;
	unsigned char **H_SYS;
	unsigned char **H_NEW;
	unsigned char **G_SYS;

	int **row_in_col;	// row_in_col : col_deg * block_lenegth
	int **col_in_row;	// col_in_row : row_deg * redundancy						
	double *LRft, *LRpt;
	double **LRrtl, **LRqtl;

	int  *input_word;
	int *output_word;

public:

	LDPC();
	~LDPC();

	double cross_over_probability;
	int* Get_Word(int type);
	double Get_Param(int type);
	bool Set_Param(int block_length, int message_length, int col_deg, int row_deg, int field_size, double cross_over_probability);

	bool Print_Parity_Check_Matrix(char *name);
	bool Print_Sparse_Matrix(int type, char *name);
	bool Print_Word(int type);

	// Test parts.
	void Generate_Input_Word(unsigned char*ptr);
	void Generate_Code_Word(int value);
	void Generate_Error_Word(int value);

	// Construction of both parity and geneartor matrices. 
	bool Make_Sparse_Matrix_From_Parity_Check_Matrix();
	bool Make_Parity_Check_Matrix_Sys();
	bool Make_Gallager_Parity_Check_Matrix(int seed);
	void FindGoodRowCol(int pivot, int &target_row, int &target_col);	//M=>redundancy
	bool Matrix_Swap(int i, int j, int type);


	// Message_Passing decoding
	bool LDPC_Decoding();
	bool IsCodeword();
	double func_f(double x);
	double infinity_test(double x);
};


