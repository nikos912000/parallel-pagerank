/******************** Structs ********************/

/***** Struct for timestamps *****/
struct timeval start,end;

/***** Struct used for Threads data *****/

typedef struct
{
	int tid;
	int start, end;
} Thread; 

/***** Struct used for Nodes data *****/

typedef struct
{
	double p_t0;
	double p_t1;
	double e;
	int *From_id;
	int con_size;
	int from_size;
}Node;