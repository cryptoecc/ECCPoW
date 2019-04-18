#include "utility_function/memory_manage.h"
#include "decoder/LDPC.h"
#include "utility_function/Def_List.h"
#include "block/BlockHeader.h"
#include <stdio.h>
#include <time.h>
int ex(LDPC *decoder, double crossover_prob, char *name, int error_weight, int num)
{
	int succ = 0;
	int iter = 0, Iter = 0;
	for (int i = 0; i < num; i++)
	{
		decoder->Generate_Code_Word(0);
		decoder->Generate_Error_Word(p_m);
		if (decoder->LDPC_Decoding(&Iter))
		{
			succ = succ + 1;
			iter = Iter + iter;
		}
	}
	return succ;
}