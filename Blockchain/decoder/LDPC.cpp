#include "LDPC.h"
#include "..\utility_function\Def_List.h"
#include "..\utility_function\Memory_Manage.h"
#include "..\utility_function\sha256.h"
#include <algorithm>    
#include <vector>
#include <time.h>

LDPC::LDPC()
{
	this->row_deg = 0;  this->col_deg = 0; this->block_length = 0, this->message_length = 0, this->seed = 0, this->q = 0;
	this->redundancy = 0; this->lastnonzerow = 0; this->Iter = 0;
	this->code_rate = 0.0;

	this->H = NULL;
	this->H_ORI = NULL;
	this->H_SYS = NULL;
	this->H_NEW = NULL;
	this->G_SYS = NULL;

	this->row_in_col = NULL;	// row_in_col : col_deg * block_lenegth
	this->col_in_row = NULL;	// col_in_row : row_deg * redundancy						
	this->LRft = NULL;
	this->LRpt = NULL;
	this->LRrtl = NULL;
	this->LRqtl = NULL;

	this->input_word = NULL;
	this->output_word = NULL;
	
}
LDPC::~LDPC()
{
	Delete_2D_Array((void**)this->H,     this->redundancy);// (char*)"Failure of freeing H variable at destruction function of LDPC");
	Delete_2D_Array((void**)this->H_ORI, this->redundancy);// (char*)"Failure of freeing H_ORI variable at destruction function of LDPC");
	Delete_2D_Array((void**)this->H_SYS, this->redundancy);// (char*)"Failure of freeing H_SYS variable at destruction function of LDPC");
	Delete_2D_Array((void**)this->H_NEW, this->redundancy);// (char*)"Failure of freeing H_NEW variable at destruction function of LDPC");
	Delete_2D_Array((void**)this->G_SYS, this->message_length);// (char*)"Failure of freeing H_ORI variable at destruction function of LDPC");

	Delete_2D_Array((void**)this->col_in_row, this->row_deg);//(char*)"Failure of freeing col_in_row variable at destruction function of LDPC");
	Delete_2D_Array((void**)this->row_in_col, this->col_deg);// (char*)"Failure of freeing row_in_col variable at destruction function of LDPC");
	Delete_2D_Array((void**)this->LRqtl, this->block_length);// (char*)"Failure of freeing LRqtl variable at destruction function of LDPC");
	Delete_2D_Array((void**)this->LRrtl, this->block_length);// (char*)"Failure of freeing LRrtl variable at destruction function of LDPC");
	Delete_1D_Array((void*)this->LRpt);// (char*)"Failure of freeing LRpt variable at destruction function of LDPC");
	Delete_1D_Array((void*)this->LRft);// (char*)"Failure of freeing LRft variable at destruction function of LDPC");

	Delete_1D_Array((void*)this->input_word);// (char*)"Failure of freeing input word variable at destruction function of LDPC");
	Delete_1D_Array((void*)this->output_word);// (char*)"Failure of freeing output word variable at destruction function of LDPC");
}

bool LDPC::Make_Gallager_Parity_Check_Matrix(int seed)
{
	this->seed = seed;
	std::vector<int> col_order;
	if (!(this->block_length * this->col_deg == this->redundancy * this->row_deg))
	{
		fprintf(stderr, "\nThis input cannot be used to make Gallagher parity check matrix\n");
		fprintf(stderr, "Blocklength*Col_deg must be equal to (Blocklength-Messagelength)*Row_deg\n");
		return false;
	}

	Delete_2D_Array((void**)this->H, this->redundancy); this->H = NULL;
	this->H = Allocate_2D_Array_UChar(this->redundancy, this->block_length, (char*)"No sufficient memory for H at Make_Gallager_Parity_Check_Matrix");
	int k = this->redundancy/this->col_deg;

	for (int i = 0; i < k; i++)
		for (int j = i*this->row_deg; j < (i + 1)*this->row_deg; j++)
			this->H[i][j] = 1;

	for (int i = 1; i < col_deg; i++)
	{
		col_order.clear();
		for (int j = 0; j <this->block_length; j++)
			col_order.push_back(j);
		std::srand(seed--);
		std::random_shuffle(col_order.begin(), col_order.end());
		for (int j = 0; j <this->block_length; j++)
		{
			int index = (col_order.at(j) / this->row_deg + k * i);
			H[index][j] = 1;
		}
	}
	if (!Make_Parity_Check_Matrix_Sys())
		return false;
	return true;
}
bool LDPC::Make_Parity_Check_Matrix_Sys()
{
	lastnonzerow = this->redundancy;
	int target_row, target_col;
	unsigned char temp_int;
	unsigned char *swapmap = NULL;
	unsigned char *temp_swapmap = NULL;
	
	Delete_2D_Array((void**)this->H_ORI, this->redundancy); this->H_ORI = NULL;
	Delete_2D_Array((void**)this->H_NEW, this->redundancy); this->H_NEW = NULL;
	Delete_2D_Array((void**)this->H_SYS, this->redundancy); this->H_SYS = NULL;
	Delete_2D_Array((void**)this->G_SYS, this->message_length); this->G_SYS = NULL;

	this->H_ORI = Allocate_2D_Array_UChar(this->redundancy, this->block_length, (char*)"No sufficient memory for H_ORI at Make_Pairty_Check_Matrix_Sys Function");
	this->H_NEW = Allocate_2D_Array_UChar(this->redundancy,   this->block_length, (char*)"No sufficient memory for H_NEW at Make_Parity_Check_Matrix_Sys Function");
	this->H_SYS = Allocate_2D_Array_UChar(this->redundancy,   this->block_length, (char*)"No sufficient memory for H_SYS at Make_Parity_Check_Matrix_Sys Function");
	this->G_SYS = Allocate_2D_Array_UChar(this->block_length, this->block_length, (char*)"No sufficient matrix memory for G_SYS at Make_Generator_Matrix_Sys Function");

	swapmap = Allocate_1D_Array_UChar(this->block_length, (char*)"No sufficient memory for swapmap at Make_Parity_Check_Matrix_Sys Function");
	temp_swapmap= Allocate_1D_Array_UChar(this->block_length, (char*)"No sufficient memory for temp_swapmap at Make_Parity_Check_Matrix_Sys Function");
	
	if (!this->H_SYS || !this->H_NEW || !this->H_ORI || !this->G_SYS || !temp_swapmap || !swapmap)
		return false;

	for (int i = 0; i < this->redundancy; i++)
		memcpy(this->H_ORI[i], this->H[i], sizeof(unsigned char)*this->block_length);

	for (int i = 0; i < this->block_length; i++)
		swapmap[i] = i;

	for (int pivot = 0; pivot < this->redundancy; pivot++)
	{
		FindGoodRowCol(pivot, target_row, target_col);
		Matrix_Swap(pivot, target_row, ROW_SWAP);
		Matrix_Swap(pivot, target_col, COLUMN_SWAP);

		temp_int = swapmap[pivot];
		swapmap[pivot] = swapmap[target_col];
		swapmap[target_col] = temp_int;               // Keep tracking the column swap.

		for (int i = pivot + 1; i < this->redundancy; i++)
		{
			if (this->H[i][pivot])
				for (int j = 0; j < this->block_length; j++)
					this->H[i][j] = (this->H[i][j] + this->H[pivot][j]) % 2;// subjection
		}
	}

	for (int i = this->redundancy - 1; i >= 0; i--)
	{
		for (int k = 0; k < i; k++)
		{
			if (this->H[k][i] == 1)
				for (int j = 0; j < this->block_length; j++)
					this->H[k][j] = (this->H[k][j] + this->H[i][j]) % 2;         // backsubtitution
		}
	}

	temp_int = 0;
	for (int i = this->redundancy - 1; i >= 0; i--)
	{
		for (int j = 0; j < this->block_length; j++)
			temp_int = temp_int + this->H[i][j];

		if (temp_int > 0)
		{
			lastnonzerow = i;                           // Find the index of last nonzero row
			break;
		}
	}
	
	memcpy(temp_swapmap, swapmap, sizeof(unsigned char)*this->block_length);
	for (int i = 0; i < lastnonzerow + 1; i++)
	{
		for (int j = 0; j < this->block_length - lastnonzerow - 1; j++)
		{
			this->H_SYS[i][j] = this->H[i][j + lastnonzerow + 1];
			swapmap[j] = temp_swapmap[j + lastnonzerow + 1];
		}
		for (int j = this->block_length - lastnonzerow - 1; j < this->block_length; j++)
		{
			this->H_SYS[i][j] = this->H[i][j - this->block_length + lastnonzerow + 1];
			swapmap[j] = temp_swapmap[j - this->block_length + lastnonzerow + 1];
		}
	}
	for (int i = 0; i < this->block_length; i++)
		for (int j = 0; j < this->redundancy; j++)
			this->H_NEW[j][i] = this->H_ORI[j][swapmap[i]];

	for (int i = 0; i < this->block_length - lastnonzerow - 1; i++)
		this->G_SYS[i][i] = 1;
	
	for (int i = 0; i < lastnonzerow + 1; i++)
	{
		for (int j = 0; j < this->block_length - lastnonzerow - 1; j++)
			this->G_SYS[j][i+ this->block_length - lastnonzerow - 1] = this->H_SYS[i][j];
	}
	
	Delete_2D_Array((void**)this->H_ORI, this->redundancy);
	Delete_2D_Array((void**)this->H_SYS, this->redundancy);
	Delete_1D_Array((void*)swapmap);
	Delete_1D_Array((void*)temp_swapmap);

	this->H = H_NEW;
	H_NEW = NULL; H_ORI = NULL; H_SYS = NULL;
	return true;
}
bool LDPC::Make_Sparse_Matrix_From_Parity_Check_Matrix()
{
	// col_in_row : row_deg *redundancy
	this->col_in_row = (int**)Delete_2D_Array((void**)this->col_in_row, this->row_deg);
	this->col_in_row = Allocate_2D_Array_Int(this->row_deg, this->redundancy,   (char*)"No sufficient memory for col_in_row at Make_Sparse_Matrix_From_Parity_Check_Matrix function");
	this->row_in_col = (int**)Delete_2D_Array((void**)this->row_in_col, this->col_deg);
	this->row_in_col = Allocate_2D_Array_Int(this->col_deg, this->block_length, (char*)"No sufficient memory for row_in_col at Make_Sparse_Matrix_From_Parity_Check_Matrix function");
	if (!col_in_row || !row_in_col)
		return false;
	
	int row_index = 0, col_index = 0;
	for (int i = 0; i < this->redundancy; i++)
	{
		for (int j = 0; j < this->block_length; j++)
		{
			if (this->H[i][j])
			{
				this->col_in_row[col_index++%this->row_deg][i] = j;
				this->row_in_col[row_index++/this->block_length][j] = i;			
			}
		}
	}
	return true;
}
bool LDPC::IsCodeword()
{
	bool flag = true;
	for (int i = 0; i < this->redundancy; i++)
	{
		int temp_iscd = 0;
		for (int j = 0; j < this->block_length; j++)
			temp_iscd += (this->output_word[j] * this->H[i][j]);
		if (temp_iscd % 2)
		{
			flag = false;
			break;
		}
	}
	return flag;
}
bool LDPC::LDPC_Decoding()
{
	double temp3, temp_sign, sign, magnitude, ee = 0.1;
	int sum = 0;
	memset(this->output_word, NULL            , sizeof(int)*this->block_length);
	memcpy(this->output_word, this->input_word, sizeof(int)*this->block_length);
	this->Iter = 0;
	if (IsCodeword())
		return true;

	// LRpt => L_j^Total, LRft = L_j, LRrtl = L_i->j
	for (int i = 0; i < this->block_length; i++)
	{
		memset(this->LRqtl[i], NULL, sizeof(double)*this->redundancy);
		memset(this->LRrtl[i], NULL, sizeof(double)*this->redundancy);
		this->LRft[i] = log((1 - this->cross_over_probability) / (this->cross_over_probability))*(double)(this->input_word[i] * 2 - 1);//infinity_test
	}
	memset(this->LRpt, NULL, sizeof(double)*this->block_length);

	int i,k, l, m, ind, t, mp;
	//Bit to Check Node Messages --> LRqtl
	for (ind = 1; ind <= ITERATIONS; ind++) 
	{	
		for (t = 0; t < this->block_length; t++)
		{
			for ( m = 0; m < this->col_deg; m++) 
			{
				temp3 = 0;
				for ( mp = 0; mp < this->col_deg; mp++)
				{
					if (mp != m)
						temp3 = infinity_test(temp3 + this->LRrtl[t][this->row_in_col[mp][t]]);
				}
				this->LRqtl[t][this->row_in_col[m][t]] = infinity_test(this->LRft[t] + temp3);
			}
		}

		//Check to Bit Node Messages --> LRrtl
		for (k = 0; k < this->redundancy; k++) 
		{
			for (l = 0; l < this->row_deg; l++)
			{
				temp3 = 0.0;
				sign = 1;
				for (m = 0; m < this->row_deg; m++)
				{
					if (m != l) {
						temp3 = temp3 + func_f(fabs(this->LRqtl[this->col_in_row[m][k]][k]));
						if (this->LRqtl[this->col_in_row[m][k]][k] > 0.0)
							temp_sign = 1.0;
						else
							temp_sign = -1.0;
						sign = sign * temp_sign;
					}
				}
				magnitude = func_f(temp3);
				this->LRrtl[this->col_in_row[l][k]][k] = infinity_test(sign*magnitude);
			}
		}
		
		//Last iteration get LR (pi)
		// LRpt => L_j^Total, LRft = L_j, LRrtl = L_i->j
		for (m = 0; m < this->block_length; m++)
		{
			this->LRpt[m] = infinity_test(this->LRft[m]);
			for (k = 0; k < this->col_deg; k++)
			{
				this->LRpt[m] += this->LRrtl[m][this->row_in_col[k][m]];
				this->LRpt[m] = infinity_test(this->LRpt[m]);
			}
		}

		//Make Decision
		for (i = 0; i < this->block_length; i++)
		{
			if (LRpt[i] >= 0)
				this->output_word[i] = 1;
			else
				this->output_word[i] = 0;
		}
		if (IsCodeword())
		{
			this->Iter = ind;
			return true;
		}
	}
	return false;
}
bool LDPC::Print_Parity_Check_Matrix(char *name)
{
	FILE *fp;
	if (!name)
		fp = stdout;
	else
		fopen_s(&fp, name, "w");

	fprintf_s(fp, "\nInformation regarding the current parity check matrix\n");
	fprintf_s(fp, "Block Lengh : %d\t Message Length : %d\t Code Rate : %f\n", this->block_length, this->message_length, this->code_rate);
	fprintf_s(fp, "# of ones in rows : %d\t # of ones in columms : %d\t Seed : %d\n", this->row_deg, this->col_deg, this->seed);
	fprintf_s(fp, "\nRegular LDPC Parity Check Matrix of size %d X %d \n", this->redundancy, this->block_length);
	for (int i = 0; i < this->redundancy; i++)
	{
		for (int j = 0; j < this->block_length; j++)
			fprintf_s(fp, "%d ", this->H[i][j]);
		fprintf_s(fp, "\n");
	}

	fprintf_s(fp, "\nRegular LDPC Generator Matrix of size %d X %d \n", this->block_length - lastnonzerow - 1, this->block_length);
	for (int i = 0; i < this->block_length - lastnonzerow - 1; i++)
	{
		for (int j = 0; j < this->block_length; j++)
			fprintf_s(fp, "%d ", this->G_SYS[i][j]);
		fprintf_s(fp, "\n");
	}
	if (name)
		fclose(fp);
	return true;
}
bool LDPC::Print_Sparse_Matrix(int type, char *name)
{
	FILE *fp;
	if (name)
		fopen_s(&fp, name, "w");
	else
		fp = stdout;
	if (type == ROW_IN_COLUMN)
	{
		fprintf_s(fp,"\nInformation regarding the row_in_col_matrix obtained from the current parity check matrix\n");
		for (int i = 0; i < this->col_deg; i++)
		{
			for (int j = 0; j < this->block_length; j++)
				fprintf_s(fp, "%d\t", this->row_in_col[i][j] + 1);
			fprintf(fp, "\n");
		}
	}
	else if (type == COLUMN_IN_ROW)
	{
		fprintf_s(fp,"\nInformation regarding the col_in_row_matrix obtained from the current parity check matrix\n");
		for (int i = 0; i < this->row_deg; i++)
		{
			for (int j = 0; j < this->redundancy; j++)
				fprintf_s(fp, "%d\t", this->col_in_row[i][j] + 1);
			fprintf_s(fp, "\n");
		}
	}
	if (name)
		fclose(fp);
	return true;
}
bool LDPC::Print_Word(int type)
{
	int *ptr = NULL;
	if (type == INPUT_WORD)
	{
		ptr = this->input_word;
		if (!ptr)
		{
			printf("There is no input word\n");
			return false;
		}
		printf("Input word\n");
	}
	else if (type == OUTPUT_WORD)
	{
		ptr = this->output_word;
		if (!ptr)
		{
			printf("There is no output word\n");
			return false;
		}
		printf("Output word\n");
	}
	for (int i = 0; i < this->block_length; i++)
		printf("%d ", ptr[i]);
	printf("\n");
	return true;
}
double LDPC::Get_Param(int type)
{
	if (type == BLOCK_LENGTH)
		return this->block_length;
	else if (type == MESSAGE_LENGTH)
		return this->message_length;
	else if (type == COLUMN_DEGREE)
		return this->col_deg;
	else if (type == ROW_DEGREE)
		return this->row_deg;
	else if (type == FIELD_SIZE)
		return this->q;
	else if (type == SEED)
		return this->seed;
	else if (type == CROSS_OVER_PROB)
		return this->cross_over_probability;
	else if (type == ITER)
		return this->Iter;
	else
		return -1;
}
bool LDPC::Set_Param(int block_length, int message_length, int col_deg, int row_deg, int q, double cross_over_probability)
{
	this->block_length = block_length;
	this->message_length = message_length;
	this->row_deg = row_deg;
	this->col_deg = col_deg;
	this->code_rate =((double)this->message_length / (double)this->block_length);
	this->redundancy = this->block_length - this->message_length;
	this->cross_over_probability = cross_over_probability;
	
	Delete_1D_Array((void*)this->input_word);
	Delete_1D_Array((void*)this->output_word);
	Delete_1D_Array((void*)this->LRft);
	Delete_1D_Array((void*)this->LRpt);
	Delete_2D_Array((void**)this->LRrtl,this->block_length);
	Delete_2D_Array((void**)this->LRqtl,this->block_length);


	this->input_word = NULL; this->output_word = NULL; this->LRft = NULL; this->LRpt = NULL; this->LRrtl = NULL; this->LRqtl = NULL;
	
	this->input_word  = Allocate_1D_Array_Int(this->block_length, (char*)"No sufficient memory for input_word at Set_Param Function");
	this->output_word = Allocate_1D_Array_Int(this->block_length, (char*)"No sufficient memory for output_word at Set_Param Function");
	this->LRft = Allocate_1D_Array_Double(this->block_length, (char*)"No sufficient memory for LRft at LPDC_Decoder Function");
	this->LRpt = Allocate_1D_Array_Double(this->block_length, (char*)"No sufficient memory for LRpt at LPDC_Decoder Function");
	this->LRrtl = Allocate_2D_Array_Double(this->block_length, this->redundancy, (char*)"No sufficient memory for LRrtl at LDPC_Decoder Function");
	this->LRqtl = Allocate_2D_Array_Double(this->block_length, this->redundancy, (char*)"No sufficient memory for LRqtl at LDPC_Decoder Function");

	return true;
}
int* LDPC::Get_Word(int type)
{
	if (type == INPUT_WORD)
		return this->input_word;
	else if (type == OUTPUT_WORD)
		return this->output_word;
	return NULL;
}
void LDPC::Generate_Input_Word(unsigned char *ptr)
{
	for (int i = 0; i < this->block_length / 8; i++)
	{
		int decimal = (int)ptr[i];
		for (int j = 7; j >= 0; j--)
		{
			this->input_word[j + 8 * (i)] = decimal % 2;
			decimal = decimal / 2;
		}
	}
}
double LDPC::func_f(double x)
{
	if (x >= BIG_INFINITY)
		return (double)(1.0 / BIG_INFINITY);

	else if (x <= (1.0 / BIG_INFINITY))
		return (double)(BIG_INFINITY);

	else
		return (double)(log((exp(x) + 1) / (exp(x) - 1)));
}
double LDPC::infinity_test(double x)
{
	if (x >= Inf)
		return Inf;
	else if (x <= -Inf)
		return -Inf;
	else
		return x;
}
void LDPC::FindGoodRowCol(int pivot, int &target_row, int &target_col)
{
	target_row = pivot; 	target_col = pivot;
	if (this->H[pivot][pivot])
		return;

	for (int j = pivot; j < this->block_length; j++)
	{
		for (int i = pivot; i < this->redundancy; i++)
		{
			if (this->H[i][j])
			{
				target_row = i;
				target_col = j;
				return;
			}
		}
	}
}
bool LDPC::Matrix_Swap(int i, int j, int type)
{
	int tmp;
	if (type == ROW_SWAP)
	{
		for (int l = 0; l < this->block_length; l++)
		{
			tmp = this->H[i][l];
			this->H[i][l] = this->H[j][l];
			this->H[j][l] = tmp;
		}
	}
	else if (type == COLUMN_SWAP)
	{
		for (int l = 0; l < this->redundancy; l++)
		{
			tmp = this->H[l][i];
			this->H[l][i] = this->H[l][j];
			this->H[l][j] = tmp;
		}
	}
	return true;
}
void LDPC::Generate_Code_Word(int value)
{
	int decimal = value;
	int *msg = NULL;
	this->G_SYS = Allocate_2D_Array_UChar(this->block_length, this->block_length, (char*)"No sufficient matrix memory for G_SYS at Make_Generator_Matrix_Sys Function");
	msg = Allocate_1D_Array_Int(this->block_length - this->lastnonzerow - 1, (char*)"No sufficient memory for msg at Generate_Codeword Function");
	for (int j = this->message_length - 1; j >= 0; j--)
	{
		msg[j] = decimal % 2;
		decimal = decimal / 2;
	}
	for (int j = 0; j < this->block_length; j++)
	{
		int temp = 0;

		for (int i = 0; i < this->block_length - this->lastnonzerow - 1; i++)
			temp = temp + msg[i] * this->G_SYS[i][j];
		
		this->input_word[j] = temp % 2;
	}
	Delete_1D_Array((void*)msg);
}
void LDPC::Generate_Error_Word(int value)
{
	int *arr = new int[block_length];
	for (int i = 0; i < this->block_length; i++)
		arr[i] = i;

	for (int i = 0; i < value; i++)
	{
		int index = arr[rand() % this->block_length];
		if (this->input_word[index] == 1)
			this->input_word[index] = 0;
		else
			this->input_word[index] = 1;
	}
	delete arr;
}
