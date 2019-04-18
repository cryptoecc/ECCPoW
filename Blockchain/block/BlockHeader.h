#pragma once
#include "..\utility_function\sha256.h"
#include "..\utility_function\Def_List.h"
#include "..\utility_function\Memory_Manage.h"
#include "..\decoder\LDPC.h"


class BlockHeader
{
private:
	BYTE *tmp;
	BYTE curt_hash_value[SHA256_BLOCK_SIZE];	//   current hash value 32 bytes
	BYTE prev_hash_value[SHA256_BLOCK_SIZE];	//  previous hash value 32 bytes
	BYTE merkle_tree_value[SHA256_BLOCK_SIZE];	//    merkle tree value 32 bytes
	BYTE buff[MAX_BUFF_SIZE];      //  this variable contains a byte stream constructed using Set_Buff()
	SHA256_CTX ctx;
	unsigned int version, timestamp, bits, nonce, index, seed; // version, time, difficulty, nonce, block index, seed
public:
	BlockHeader* next;
	BlockHeader();
	BlockHeader(int version, int timestamp, int bits, int nonce, int index, char prev_hash_value[SHA256_BLOCK_SIZE]);
	
	void         Set_Buff();	// Constructing a byte stream by concatenating "version" + "timestamp" + "difficulty" + "merkle tree value"
	void         Set_Merkle_Tree_Value(char MKT[SHA256_BLOCK_SIZE]);
	void         Set_Hash_Value();
	
	void         Print_Header_Info();
	bool		 Mining(LDPC *decoder);
	
	char*        Get_Hash_Value(int type);
	unsigned int Get_Param(int type);

	virtual ~BlockHeader();
};

