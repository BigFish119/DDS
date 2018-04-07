#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <string.h>
#include<string>
using namespace std;
#define NetNodeNum 34 //78 karate.txt
#define MaxNeighborNum 10000
// edge
typedef struct Edge{
	int n_1_ID;
	int n_2_ID;
	int e_ID;
	double dist;

	int * CN_Loc;
	int CN_Loc_Length;

	int * N_1_Info;
	int N_1_Info_Length;

	int * N_1_Loc;
	int N_1_Loc_Length;

	int * N_2_Info;
	int N_2_Info_Length;
	int * N_2_Loc;
	int N_2_Loc_Length;
}EA;


typedef struct EdgeLoc{
	int Row;
	int Column;
}EdgeLoc;

//community
typedef struct Community{
	int memberNum;
	int *point_to_member;
}Community;

// (NodeTable header)
typedef struct PointerToNeighbors{
	int commId;//the node's community id
	int NodeNum;//number of neighbors of the node
	int * Pointer;//points to the neighbors array
}NodeTableHead;

// (EdgeWithAttachmentTable header)
typedef struct PointerToEdges{
	int EdgeNum;
	EA * Pointer;
}EdgeTableHead;

// (EdgesOfNodeTable header)
typedef struct PointerToEdgeLoc{
	int EdgeNum;
	EdgeLoc * Pointer;
}EdgesOfNodeTableHead;

/* global variable */
NodeTableHead NodeTable[NetNodeNum];
EdgeTableHead EdgeWithAttachmentTable[NetNodeNum]; 
EdgesOfNodeTableHead EdgesOfNodeTable[NetNodeNum]; 
Community community[NetNodeNum + 1];
int CN[MaxNeighborNum];
int DiffA[MaxNeighborNum];
int DiffB[MaxNeighborNum];
int NetEdgeNum;

void establishNodeTable(char * inputfilename, int FileLine){
	//open file
	FILE * file = fopen(inputfilename, "r");
	if (file == NULL)
	{
		printf("cannot open__1 %s\n", inputfilename);
		exit(1);
	}

	int NeighborNumber[NetNodeNum];
	for (int i = 0; i < NetNodeNum; i++)
		NeighborNumber[i] = 0;
	int node_1;//n1
	int node_2;//n2

	for (int i = 0; i < FileLine; i++)
	{
		fscanf(file, "%d %d", &node_1, &node_2);
		NeighborNumber[node_1 - 1]++;
		NeighborNumber[node_2 - 1]++;
	}

	for (int i = 0; i < NetNodeNum; i++)
	{
		NodeTable[i].commId = i + 1;
		community[i + 1].memberNum = 1;
		community[i + 1].point_to_member = (int*)malloc(sizeof(int) * 1);
		community[i + 1].point_to_member[0] = i + 1;

		NodeTable[i].NodeNum = 0;
		NodeTable[i].Pointer = (int *)malloc(NeighborNumber[i] * sizeof(int));
		if (NodeTable[i].Pointer == NULL)
		{
			printf("Memory is not enough when initiating NodeTable\n");
			exit(1);
		}
	}
	fclose(file);
	file = fopen(inputfilename, "r");
	if (file == NULL)
	{
		printf("cannot open__2 %s\n", inputfilename);
		exit(1);
	}
	for (int i = 0; i < FileLine; i++)
	{
		fscanf(file, "%d %d", &node_1, &node_2);
		NodeTable[node_1 - 1].Pointer[NodeTable[node_1 - 1].NodeNum++] = node_2;
		NodeTable[node_2 - 1].Pointer[NodeTable[node_2 - 1].NodeNum++] = node_1;
	}
	for (int i = 0; i < NetNodeNum; i++)
	{
		if (NodeTable[i].NodeNum != NeighborNumber[i])
		{
			printf("NodeTable[%d] error\n", i + 1);
			exit(1);
		}
	}
	fclose(file);
	void sortFun(int *, int);
	for (int i = 0; i < NetNodeNum; i++)
	{
		sortFun(NodeTable[i].Pointer, NodeTable[i].NodeNum);
	}
}

void establishEdgeTable()
{
	int EdgeNumber[NetNodeNum];// save edge number of each node；
	//compute node[loop_i]'s neighbors number
	for (int loop_i = 0; loop_i < NetNodeNum; ++loop_i)
	{
		EdgeNumber[loop_i] = 0;
		for (int loop_j = 0; loop_j < NodeTable[loop_i].NodeNum; ++loop_j)
		{
			if (NodeTable[loop_i].Pointer[loop_j] > loop_i + 1)
			{
				++EdgeNumber[loop_i];
			}
		}
	}
	// compute n_1_ID, n_2_ID, e_ID
	int edgeid = 0;
	for (int loop_i = 0; loop_i < NetNodeNum; ++loop_i)
	{
		EdgeWithAttachmentTable[loop_i].EdgeNum = EdgeNumber[loop_i];
		EdgeWithAttachmentTable[loop_i].Pointer = (EA *)malloc(EdgeNumber[loop_i] * sizeof(EA));
		if (EdgeWithAttachmentTable[loop_i].Pointer == NULL)
		{
			printf("memory is not enough with node%d's edges\n", loop_i + 1);
			exit(1);
		}

		int edge_loc_tmp = 0;
		for (int loop_j = 0; loop_j < NodeTable[loop_i].NodeNum; ++loop_j)
		{
			if (NodeTable[loop_i].Pointer[loop_j] > loop_i + 1)
			{
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].n_1_ID = loop_i + 1;
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].n_2_ID = NodeTable[loop_i].Pointer[loop_j];
				EdgeWithAttachmentTable[loop_i].Pointer[edge_loc_tmp].e_ID = edgeid;
				++edgeid;
				NetEdgeNum = edgeid;
				++edge_loc_tmp;
			}
		}
		if (edge_loc_tmp != EdgeWithAttachmentTable[loop_i].EdgeNum)
		{
			printf("EdgeWithAttachmentTable[%d] error(line 175)\n", loop_i + 1);
			exit(1);
		}
	}
	for (int i = 0; i < NetNodeNum; ++i)
	{
		for (int j = 0; j < EdgeWithAttachmentTable[i].EdgeNum; j++)
		{
			int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
			int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
			void findNeighbors(int, int);// function declaration
			findNeighbors(node_1, node_2);
			// compute distance
			EdgeWithAttachmentTable[i].Pointer[j].dist = 1.0 - (double)(CN[0] + 2) / (double)(CN[0] + DiffA[0] + DiffB[0] + 2);

			EdgeWithAttachmentTable[i].Pointer[j].CN_Loc = (int *)malloc(4 * CN[0] * sizeof(int));
			if (EdgeWithAttachmentTable[i].Pointer[j].CN_Loc == NULL)
			{
				printf("EdgeWithAttachmentTable[%d].Pointer is NULL\n", i); exit(1);
			}
			EdgeWithAttachmentTable[i].Pointer[j].CN_Loc_Length = CN[0];
			for (int s = 1; s <= CN[0]; s++){
				int NodeCN = CN[s];
				int Loc_Node_1_R = -1;
				int Loc_Node_1_C = -1;
				int Loc_Node_2_R = -1;
				int Loc_Node_2_C = -1;
				
				int NodeMin = (node_1 < NodeCN) ? node_1 : NodeCN;
				int NodeMax = (node_1 > NodeCN) ? node_1 : NodeCN;
				Loc_Node_1_R = NodeMin - 1;
				for (int loop = 0; loop < EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
				{
					if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
					{
						Loc_Node_1_C = loop;
						break;
					}
					else
					{
						continue;
					}
				}
				
				if (Loc_Node_1_C == -1)
				{
					printf("line 207 error\n");
					exit(1);
				}
				NodeMin = (node_2 < NodeCN) ? node_2 : NodeCN;
				NodeMax = (node_2 > NodeCN) ? node_2 : NodeCN;
				Loc_Node_2_R = NodeMin - 1;
				for (int loop = 0; loop < EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
				{
					if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
					{
						Loc_Node_2_C = loop;
						break;
					}
				}
				if (Loc_Node_2_C == -1)
				{
					printf("line 261 error\n");
					exit(1);
				}
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4 * (s - 1) + 0] = Loc_Node_1_R;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4 * (s - 1) + 1] = Loc_Node_1_C;

				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4 * (s - 1) + 2] = Loc_Node_2_R;
				EdgeWithAttachmentTable[i].Pointer[j].CN_Loc[4 * (s - 1) + 3] = Loc_Node_2_C;
			}

			EdgeWithAttachmentTable[i].Pointer[j].N_1_Info = (int *)malloc(2 * 2 * DiffA[0] * sizeof(int));
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Info_Length = DiffA[0];

			if (EdgeWithAttachmentTable[i].Pointer[j].N_1_Info == NULL)
			{
				printf("memory not enough With N_1_Info\n");
				exit(1);
			}
			int Edgenum_tmp = 0;
			for (int s = 1; s <= DiffA[0]; s++)
			{
				int Node_N_1 = DiffA[s];
				void findCN(int, int);
				findCN(Node_N_1, node_2);
				Edgenum_tmp = Edgenum_tmp + CN[0];
			}
			
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc = (int *)malloc(4 * Edgenum_tmp * sizeof(int));
			if (EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc == NULL)
			{
				printf("memory not enough With N_1_Loc(line 243)\n");
				exit(1);
			}
			EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc_Length = Edgenum_tmp;//n1的所有独家邻居与n2的公共邻居的个数

			Edgenum_tmp = 0;
			for (int s = 1; s <= DiffA[0]; s++)
			{
				int Node_N_1 = DiffA[s];
				void findCN(int, int);
				findCN(Node_N_1, node_2);

				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4 * (s - 1) + 0] = Node_N_1;
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4 * (s - 1) + 1] = CN[0];

				int NodeMin;
				int NodeMax;
				NodeMin = (node_1 < Node_N_1) ? node_1 : Node_N_1;
				NodeMax = (node_1 > Node_N_1) ? node_1 : Node_N_1;
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4 * (s - 1) + 2] = NodeMin - 1;//row
				EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4 * (s - 1) + 3] = -1;
				for (int loop = 0; loop<EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
				{
					if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
					{
						EdgeWithAttachmentTable[i].Pointer[j].N_1_Info[4 * (s - 1) + 3] = loop;
						break;
					}
				}
				if (EdgeWithAttachmentTable[NodeMin - 1].EdgeNum == -1)
				{
					printf("line 270 error\n");
					exit(1);
				}
				for (int ss = 1; ss <= CN[0]; ss++)
				{
					int NodeCN = CN[ss];
					int Loc_Node_N_1_R = -1;
					int Loc_Node_N_1_C = -1;
					int Loc_Node_2_R = -1;
					int Loc_Node_2_C = -1;

					NodeMin = (Node_N_1 < NodeCN) ? Node_N_1 : NodeCN;
					NodeMax = (Node_N_1 > NodeCN) ? Node_N_1 : NodeCN;
					Loc_Node_N_1_R = NodeMin - 1;
					for (int loop = 0; loop < EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
					{
						if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
						{
							Loc_Node_N_1_C = loop;	break;
						}
					}
					NodeMin = (node_2 < NodeCN) ? node_2 : NodeCN;
					NodeMax = (node_2 > NodeCN) ? node_2 : NodeCN;
					Loc_Node_2_R = NodeMin - 1;
					for (int loop = 0; loop<EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
					{
						if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
						{
							Loc_Node_2_C = loop;
							break;
						}
					}
					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4 * Edgenum_tmp + 0] = Loc_Node_N_1_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4 * Edgenum_tmp + 1] = Loc_Node_N_1_C;
					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4 * Edgenum_tmp + 2] = Loc_Node_2_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_1_Loc[4 * Edgenum_tmp + 3] = Loc_Node_2_C;
					Edgenum_tmp++;
				}
			}
			
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Info = (int *)malloc(2 * 2 * DiffB[0] * sizeof(int));
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Info_Length = DiffB[0];

			Edgenum_tmp = 0;
			for (int s = 1; s <= DiffB[0]; s++)
			{
				int Node_N_2 = DiffB[s];
				void findCN(int, int);
				findCN(Node_N_2, node_1);
				Edgenum_tmp = Edgenum_tmp + CN[0];
			}
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc = (int *)malloc(4 * Edgenum_tmp * sizeof(int));
			EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc_Length = Edgenum_tmp;

			Edgenum_tmp = 0;
			for (int s = 1; s <= DiffB[0]; s++){
				int Node_N_2 = DiffB[s];
				void findCN(int, int);
				findCN(Node_N_2, node_1);

				EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4 * (s - 1) + 0] = Node_N_2;
				EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4 * (s - 1) + 1] = CN[0];

				int NodeMin;
				int NodeMax;

				NodeMin = (node_2 < Node_N_2) ? node_2 : Node_N_2;
				NodeMax = (node_2 > Node_N_2) ? node_2 : Node_N_2;
				EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4 * (s - 1) + 2] = NodeMin - 1;
				for (int loop = 0; loop < EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
				{
					if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
					{
						EdgeWithAttachmentTable[i].Pointer[j].N_2_Info[4 * (s - 1) + 3] = loop;
						break;
					}
				}
				for (int ss = 1; ss <= CN[0]; ss++)
				{
					int NodeCN = CN[ss];
					int Loc_Node_N_2_R;
					int Loc_Node_N_2_C;
					int Loc_Node_1_R;
					int Loc_Node_1_C;
					NodeMin = (Node_N_2 < NodeCN) ? Node_N_2 : NodeCN;
					NodeMax = (Node_N_2 > NodeCN) ? Node_N_2 : NodeCN;
					Loc_Node_N_2_R = NodeMin - 1;
					for (int loop = 0; loop < EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
					{
						if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
						{
							Loc_Node_N_2_C = loop;
							break;
						}
					}
					NodeMin = (node_1 < NodeCN) ? node_1 : NodeCN;
					NodeMax = (node_1 > NodeCN) ? node_1 : NodeCN;
					Loc_Node_1_R = NodeMin - 1;
					for (int loop = 0; loop < EdgeWithAttachmentTable[NodeMin - 1].EdgeNum; loop++)
					{
						if (EdgeWithAttachmentTable[NodeMin - 1].Pointer[loop].n_2_ID == NodeMax)
						{
							Loc_Node_1_C = loop;
							break;
						}
					}
					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4 * Edgenum_tmp + 0] = Loc_Node_N_2_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4 * Edgenum_tmp + 1] = Loc_Node_N_2_C;
					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4 * Edgenum_tmp + 2] = Loc_Node_1_R;
					EdgeWithAttachmentTable[i].Pointer[j].N_2_Loc[4 * Edgenum_tmp + 3] = Loc_Node_1_C;
					Edgenum_tmp++;
				}
			}
		}
	}
}

void establishEdgesOfNodeTable()
{
	for (int i = 0; i < NetNodeNum; ++i)
	{
		EdgesOfNodeTable[i].EdgeNum = 0;
		EdgesOfNodeTable[i].Pointer = (EdgeLoc *)malloc(NodeTable[i].NodeNum * sizeof(EdgeLoc));
		if (EdgesOfNodeTable[i].Pointer == NULL)
		{
			printf("memory not enough with EdgesOfNodeTable(line 378)\n");
			exit(1);
		}
	}
	for (int i = 0; i < NetNodeNum; i++)
	{
		for (int j = 0; j < EdgeWithAttachmentTable[i].EdgeNum; j++)
		{
			int node_1 = EdgeWithAttachmentTable[i].Pointer[j].n_1_ID;
			int node_2 = EdgeWithAttachmentTable[i].Pointer[j].n_2_ID;
			if (node_1 != i + 1)
			{
				printf("line 386 error\n");
				exit(1);
			}
			EdgesOfNodeTable[node_1 - 1].Pointer[EdgesOfNodeTable[node_1 - 1].EdgeNum].Row = i;
			EdgesOfNodeTable[node_1 - 1].Pointer[EdgesOfNodeTable[node_1 - 1].EdgeNum].Column = j;
			EdgesOfNodeTable[node_1 - 1].EdgeNum++;

			EdgesOfNodeTable[node_2 - 1].Pointer[EdgesOfNodeTable[node_2 - 1].EdgeNum].Row = i;
			EdgesOfNodeTable[node_2 - 1].Pointer[EdgesOfNodeTable[node_2 - 1].EdgeNum].Column = j;
			EdgesOfNodeTable[node_2 - 1].EdgeNum++;
		}
	}

	for (int i = 0; i < NetNodeNum; i++)
	{
		if (EdgesOfNodeTable[i].EdgeNum != NodeTable[i].NodeNum)
		{
			printf("node_%d: Edges of this node does not math neighbors of this node(line 398 error)\n", i + 1);
			exit(1);
		}
	}
}

/* interaction */
void interaction(int FileLine)
{
	if (NetEdgeNum != FileLine)
	{
		printf("NetEdgeNum error...(line 413)\n");
		exit(1);
	}
	double * D = (double *)malloc(sizeof(double) * NetEdgeNum);
	int EdgeLocCount = 0;
	for (int s_1 = 0; s_1 < NetNodeNum; s_1++)
	{
		for (int s_2 = 0; s_2 < EdgeWithAttachmentTable[s_1].EdgeNum; s_2++)
		{
			D[EdgeLocCount++] = EdgeWithAttachmentTable[s_1].Pointer[s_2].dist;
		}
	}



	if (EdgeLocCount != NetEdgeNum)
	{
		printf("line 424 error\n");
		exit(1);
	}
	int Terminate = 1;
	int Loop = 0;
	while (Terminate)
	{
		Loop++;

		void MergenceAndDivision(double *DistanceArray, int NumOfEdges);
		MergenceAndDivision(D, NetEdgeNum);
		for (int s_1 = 0; s_1 < NetNodeNum; s_1++)
		{
			for (int s_2 = 0; s_2 < EdgeWithAttachmentTable[s_1].EdgeNum; s_2++)
			{

				EA ThisEA = EdgeWithAttachmentTable[s_1].Pointer[s_2];
				int ThisEdgeId = ThisEA.e_ID;
				if (D[ThisEdgeId] > 0.0 && D[ThisEdgeId] < 1.0)
				{
					int ThisNode_1 = ThisEA.n_1_ID;
					int ThisNode_2 = ThisEA.n_2_ID;

					double DI = 0.0;
					double CI = 0.0;
					double N_1_I = 0.0;
					double N_2_I = 0.0;
					DI = (1.0 / (double)NodeTable[ThisNode_1 - 1].NodeNum + 1.0 / (double)NodeTable[ThisNode_2 - 1].NodeNum) * sin(1 - ThisEA.dist);
					for (int s_3 = 0; s_3 < ThisEA.CN_Loc_Length; s_3++)
					{
						double Distance_CI_1;
						double Distance_CI_2;
						int * CNLocTmp = ThisEA.CN_Loc;

						Distance_CI_1 = EdgeWithAttachmentTable[CNLocTmp[s_3 * 4 + 0]].Pointer[CNLocTmp[s_3 * 4 + 1]].dist;
						Distance_CI_2 = EdgeWithAttachmentTable[CNLocTmp[s_3 * 4 + 2]].Pointer[CNLocTmp[s_3 * 4 + 3]].dist;

						CI = CI + ((double)(1 - Distance_CI_2) / (double)NodeTable[ThisNode_1 - 1].NodeNum) * sin(1 - Distance_CI_1);
						CI = CI + ((double)(1 - Distance_CI_1) / (double)NodeTable[ThisNode_2 - 1].NodeNum) * sin(1 - Distance_CI_2); //耦合函数为sin(1 - x)
					}
					int s_3 = 0;
					int s_4 = 0;
					int s_4_count;
					while (s_3 < ThisEA.N_1_Info_Length)
					{
						int Ngh_N_1 = ThisEA.N_1_Info[s_3 * 4 + 0];
						int CN_Num = ThisEA.N_1_Info[s_3 * 4 + 1];
						double Distance_N_N_1 = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]].dist;
						s_4_count = s_4;
						double lambda_numerator_NghN12 = 0.0;
						for (s_4 = s_4_count; s_4 < s_4_count + CN_Num; s_4++)
						{
							int loc_1;
							int loc_2;

							loc_1 = ThisEA.N_1_Loc[s_4 * 4 + 0];
							loc_2 = ThisEA.N_1_Loc[s_4 * 4 + 1];
							lambda_numerator_NghN12 = lambda_numerator_NghN12 + (1 - EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
							loc_1 = ThisEA.N_1_Loc[s_4 * 4 + 2];
							loc_2 = ThisEA.N_1_Loc[s_4 * 4 + 3];
							lambda_numerator_NghN12 = lambda_numerator_NghN12 + (1 - EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
						}

						double lambda_denominator_NghN12 = 0;
						for (int s_5 = 0; s_5 < NodeTable[Ngh_N_1 - 1].NodeNum; s_5++)
						{
							int min = (NodeTable[Ngh_N_1 - 1].Pointer[s_5] < Ngh_N_1) ? NodeTable[Ngh_N_1 - 1].Pointer[s_5] : Ngh_N_1;
							int max = (NodeTable[Ngh_N_1 - 1].Pointer[s_5] > Ngh_N_1) ? NodeTable[Ngh_N_1 - 1].Pointer[s_5] : Ngh_N_1;
							int loc_r = min - 1;
							int loc_c;

							int findEdge(EA *pointer, int num, int n);
							loc_c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (loc_c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN12 = lambda_denominator_NghN12 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						for (int s_6 = 0; s_6 < NodeTable[ThisNode_2 - 1].NodeNum; ++s_6)
						{
							int min = (NodeTable[ThisNode_2 - 1].Pointer[s_6] < ThisNode_2) ? NodeTable[ThisNode_2 - 1].Pointer[s_6] : ThisNode_2;
							int max = (NodeTable[ThisNode_2 - 1].Pointer[s_6] > ThisNode_2) ? NodeTable[ThisNode_2 - 1].Pointer[s_6] : ThisNode_2;
							int r = min - 1;
							int c;

							int findEdge(EA *pointer, int num, int n);
							c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN12 = lambda_denominator_NghN12 + (1 - EdgeWithAttachmentTable[r].Pointer[c].dist);
						}
						double lambda_NghN12 = (double)lambda_numerator_NghN12 / (double)lambda_denominator_NghN12;
						EA temp = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]];
						double lambda_numerator_NghN11 = 0.0;
						for (int s_7 = 0; s_7 < EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]].CN_Loc_Length; ++s_7)
						{
							int loc_r;
							int loc_c;
							loc_r = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]].CN_Loc[s_7 * 4 + 0];
							loc_c = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]].CN_Loc[s_7 * 4 + 1];
							lambda_numerator_NghN11 = lambda_numerator_NghN11 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
							loc_r = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]].CN_Loc[s_7 * 4 + 2];
							loc_c = EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]].CN_Loc[s_7 * 4 + 3];
							lambda_numerator_NghN11 = lambda_numerator_NghN11 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						lambda_numerator_NghN11 = lambda_numerator_NghN11 + (1 - EdgeWithAttachmentTable[ThisEA.N_1_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_1_Info[s_3 * 4 + 3]].dist);
						double lambda_denominator_NghN11 = 0.0;
						for (int s_8 = 0; s_8 < NodeTable[Ngh_N_1 - 1].NodeNum; s_8++)
						{
							int min = (NodeTable[Ngh_N_1 - 1].Pointer[s_8] < Ngh_N_1) ? NodeTable[Ngh_N_1 - 1].Pointer[s_8] : Ngh_N_1;
							int max = (NodeTable[Ngh_N_1 - 1].Pointer[s_8] > Ngh_N_1) ? NodeTable[Ngh_N_1 - 1].Pointer[s_8] : Ngh_N_1;
							int loc_r = min - 1;
							int loc_c;

							int findEdge(EA *pointer, int num, int n);
							loc_c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (loc_c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN11 = lambda_denominator_NghN11 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						for (int s_9 = 0; s_9 < NodeTable[ThisNode_1 - 1].NodeNum; ++s_9)
						{
							int min = (NodeTable[ThisNode_1 - 1].Pointer[s_9] < ThisNode_1) ? NodeTable[ThisNode_1 - 1].Pointer[s_9] : ThisNode_1;
							int max = (NodeTable[ThisNode_1 - 1].Pointer[s_9] > ThisNode_1) ? NodeTable[ThisNode_1 - 1].Pointer[s_9] : ThisNode_1;
							int loc_r = min - 1;
							int loc_c;
							int findEdge(EA *pointer, int num, int n);
							loc_c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (loc_c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN11 = lambda_denominator_NghN11 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						lambda_denominator_NghN11 = lambda_denominator_NghN11 - (1 - temp.dist);
						double lambda_NghN11 = lambda_numerator_NghN11 / lambda_denominator_NghN11;
						double parameter = (lambda_NghN12 > lambda_NghN11) ? 1 : -1;
						if (parameter == 1)
						{
							N_1_I = N_1_I + lambda_NghN12 / (double)NodeTable[ThisNode_1 - 1].NodeNum * sin(1 - Distance_N_N_1);
						}
						else
						{
							N_1_I = N_1_I + (lambda_NghN12 - lambda_NghN11) / (double)NodeTable[ThisNode_1 - 1].NodeNum * sin(1 - Distance_N_N_1);
						}
						s_3++;
					}
					s_3 = 0;
					s_4 = 0;
					while (s_3 < ThisEA.N_2_Info_Length)
					{
						int Ngh_N_2 = ThisEA.N_2_Info[s_3 * 4 + 0];
						int CN_Num = ThisEA.N_2_Info[s_3 * 4 + 1];
						double Distance_N_N_2 = EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]].dist;
						s_4_count = s_4;
						double lambda_numerator_NghN21 = 0;
						for (s_4 = s_4_count; s_4 < s_4_count + CN_Num; ++s_4)
						{
							int loc_1;
							int loc_2;
							loc_1 = ThisEA.N_2_Loc[s_4 * 4 + 0];
							loc_2 = ThisEA.N_2_Loc[s_4 * 4 + 1];
							lambda_numerator_NghN21 = lambda_numerator_NghN21 + (1 - EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
							loc_1 = ThisEA.N_2_Loc[s_4 * 4 + 2];
							loc_2 = ThisEA.N_2_Loc[s_4 * 4 + 3];
							lambda_numerator_NghN21 = lambda_numerator_NghN21 + (1 - EdgeWithAttachmentTable[loc_1].Pointer[loc_2].dist);
						}
						double lambda_denominator_NghN21 = 0;
						for (int s_5 = 0; s_5 < NodeTable[Ngh_N_2 - 1].NodeNum; ++s_5)
						{
							int min = (NodeTable[Ngh_N_2 - 1].Pointer[s_5] < Ngh_N_2) ? NodeTable[Ngh_N_2 - 1].Pointer[s_5] : Ngh_N_2;
							int max = (NodeTable[Ngh_N_2 - 1].Pointer[s_5] > Ngh_N_2) ? NodeTable[Ngh_N_2 - 1].Pointer[s_5] : Ngh_N_2;
							int loc_r = min - 1;
							int loc_c;

							int findEdge(EA *pointer, int num, int n);
							loc_c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (loc_c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN21 = lambda_denominator_NghN21 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						for (int p = 0; p < NodeTable[ThisNode_1 - 1].NodeNum; ++p)
						{
							int min = (NodeTable[ThisNode_1 - 1].Pointer[p] < ThisNode_1) ? NodeTable[ThisNode_1 - 1].Pointer[p] : ThisNode_1;
							int max = (NodeTable[ThisNode_1 - 1].Pointer[p] > ThisNode_1) ? NodeTable[ThisNode_1 - 1].Pointer[p] : ThisNode_1;
							int r = min - 1;
							int c;

							int findEdge(EA *pointer, int num, int n);
							c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN21 = lambda_denominator_NghN21 + (1 - EdgeWithAttachmentTable[r].Pointer[c].dist);
						}
						double lambda_NghN21 = (double)lambda_numerator_NghN21 / (double)lambda_denominator_NghN21;
						EA temp = EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]];
						double lambda_numerator_NghN22 = 0.0;
						for (int s_6 = 0; s_6 < EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]].CN_Loc_Length; ++s_6)
						{
							int loc_r;
							int loc_c;
							loc_r = EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]].CN_Loc[s_6 * 4 + 0];
							loc_c = EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]].CN_Loc[s_6 * 4 + 1];
							lambda_numerator_NghN22 = lambda_numerator_NghN22 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
							loc_r = EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]].CN_Loc[s_6 * 4 + 2];
							loc_c = EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]].CN_Loc[s_6 * 4 + 3];
							lambda_numerator_NghN22 = lambda_numerator_NghN22 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						lambda_numerator_NghN22 = lambda_numerator_NghN22 + (1 - EdgeWithAttachmentTable[ThisEA.N_2_Info[s_3 * 4 + 2]].Pointer[ThisEA.N_2_Info[s_3 * 4 + 3]].dist);
						double lambda_denominator_NghN22 = 0.0;
						for (int s_5 = 0; s_5 < NodeTable[Ngh_N_2 - 1].NodeNum; s_5++)
						{
							int min = (NodeTable[Ngh_N_2 - 1].Pointer[s_5] < Ngh_N_2) ? NodeTable[Ngh_N_2 - 1].Pointer[s_5] : Ngh_N_2;
							int max = (NodeTable[Ngh_N_2 - 1].Pointer[s_5] > Ngh_N_2) ? NodeTable[Ngh_N_2 - 1].Pointer[s_5] : Ngh_N_2;
							int loc_r = min - 1;
							int loc_c;

							int findEdge(EA *pointer, int num, int n);
							loc_c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (loc_c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN22 = lambda_denominator_NghN22 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						for (int s_5 = 0; s_5 < NodeTable[ThisNode_2 - 1].NodeNum; ++s_5)
						{
							int min = (NodeTable[ThisNode_2 - 1].Pointer[s_5] < ThisNode_2) ? NodeTable[ThisNode_2 - 1].Pointer[s_5] : ThisNode_2;
							int max = (NodeTable[ThisNode_2 - 1].Pointer[s_5] > ThisNode_2) ? NodeTable[ThisNode_2 - 1].Pointer[s_5] : ThisNode_2;
							int loc_r = min - 1;
							int loc_c;
							int findEdge(EA *pointer, int num, int n);
							loc_c = findEdge(EdgeWithAttachmentTable[min - 1].Pointer, EdgeWithAttachmentTable[min - 1].EdgeNum, max);
							if (loc_c == -1)
							{
								printf("error\n");
								exit(0);
							}
							lambda_denominator_NghN22 = lambda_denominator_NghN22 + (1 - EdgeWithAttachmentTable[loc_r].Pointer[loc_c].dist);
						}
						lambda_denominator_NghN22 = lambda_denominator_NghN22 - (1 - temp.dist);
						double lambda_NghN22 = (double)lambda_numerator_NghN22 / (double)lambda_denominator_NghN22;
						double parameter = (lambda_NghN21 > lambda_NghN22) ? 1 : -1;
						if (parameter == 1)
						{
							N_2_I = N_2_I + lambda_NghN21 / (double)NodeTable[ThisNode_2 - 1].NodeNum * sin(1 - Distance_N_N_2);
						}
						else
						{
							N_2_I = N_2_I + (lambda_NghN21 - lambda_NghN22) / (double)NodeTable[ThisNode_2 - 1].NodeNum * sin(1 - Distance_N_N_2);
						}
						s_3++;
					}
					D[ThisEdgeId] = D[ThisEdgeId] + (-(N_1_I + N_2_I) - DI - CI);
					if (D[ThisEdgeId] < 0.0)
					{
						D[ThisEdgeId] = 0.0;
					}

					if (D[ThisEdgeId] > 1.0)
					{
						D[ThisEdgeId] = 1.0;
					}
				}
			}
		}
		double sum_1 = 0;
		double sum_2 = 0;
		int EdgeCounter = 0;
		for (int s_1 = 0; s_1 < NetNodeNum; s_1++)
		{
			for (int s_2 = 0; s_2 < EdgeWithAttachmentTable[s_1].EdgeNum; s_2++)
			{
				sum_1 = sum_1 + EdgeWithAttachmentTable[s_1].Pointer[s_2].dist;
				sum_2 = sum_2 + D[EdgeCounter];
				EdgeCounter++;
			}
		}
		if (sum_1 == sum_2 || Loop > 1000)
		{
			Terminate = 0;
		}
		EdgeCounter = 0;
		for (int s_1 = 0; s_1 < NetNodeNum; s_1++)
		{
			for (int s_2 = 0; s_2 < EdgeWithAttachmentTable[s_1].EdgeNum; s_2++)
			{
				EdgeWithAttachmentTable[s_1].Pointer[s_2].dist = D[EdgeCounter];
				EdgeCounter++;
			}
		}
	}
	printf("total iteration is: %d times;\n", Loop);
}



int counter = 0;
void outputCommunities(char * outputfile)
{
	FILE * fout = fopen(outputfile, "w");
	if (fout == NULL)
	{
		printf("opening outputfile fails\n");
		exit(0);
	}

	for (int i = 1; i <= NetNodeNum; ++i)
	{
		if (community[i].memberNum > 0)
		{
			++counter;
			for (int j = 0; j < community[i].memberNum - 1; ++j)
			{
				fprintf(fout, "%d	%d\n", community[i].point_to_member[j], i);
			}
			fprintf(fout, "%d	%d\n", community[i].point_to_member[community[i].memberNum - 1], i);
		}
	}
	fclose(fout);
}

void findNeighbors(int node_i, int node_j){

	int num1 = NodeTable[node_i - 1].NodeNum;
	int * A1 = NodeTable[node_i - 1].Pointer;
	int num2 = NodeTable[node_j - 1].NodeNum;
	int * A2 = NodeTable[node_j - 1].Pointer;

	int p1_loc = 0;
	int p2_loc = 0;
	int cn_length = 0;
	int diffa_length = 0;
	int diffb_length = 0;
	int cn_loc = 1;
	int diffa_loc = 1;
	int diffb_loc = 1;
	while (p1_loc < num1 && p2_loc < num2)
	{
		if (A1[p1_loc] < A2[p2_loc])
		{
			if (A1[p1_loc] != node_j)
			{
				DiffA[diffa_loc] = A1[p1_loc];
				diffa_length++;
				diffa_loc++;
				p1_loc++;
			}
			else
			{
				p1_loc++;
			}

		}
		else if (A1[p1_loc] == A2[p2_loc])
		{
			CN[cn_loc] = A1[p1_loc];
			cn_length++;
			cn_loc++;
			p1_loc++;
			p2_loc++;
		}
		else
		{
			if (A2[p2_loc] != node_i)
			{
				DiffB[diffb_loc] = A2[p2_loc];
				diffb_length++;
				diffb_loc++;
				p2_loc++;
			}
			else
				p2_loc++;
		}
	}
	if (p1_loc == num1)
	{
		while (p2_loc < num2)
		{
			if (A2[p2_loc] != node_i)
			{
				DiffB[diffb_loc] = A2[p2_loc];
				diffb_length++;
				diffb_loc++;
				p2_loc++;
			}
			else
			{
				p2_loc++;
			}
		}
	}
	else
	{
		while (p1_loc < num1)
		{
			if (A1[p1_loc] != node_j)
			{
				DiffA[diffa_loc] = A1[p1_loc];
				diffa_length++;
				diffa_loc++;
				p1_loc++;
			}
			else
			{
				p1_loc++;
			}
		}
	}

	if (cn_loc != cn_length + 1 || diffa_loc != diffa_length + 1 || diffb_loc != diffb_length + 1)
	{
		printf("error in find common neighbors(IN FUNTION findNeighbors)\n");
		exit(1);
	}
	CN[0] = cn_length;
	DiffA[0] = diffa_length;
	DiffB[0] = diffb_length;
}
void MergenceAndDivision(double *DistanceArray, int NumOfEdges)
{
	int i = 0;
	for (int s_1 = 0; s_1 < NetNodeNum; s_1++)
	{
		for (int s_2 = 0; s_2 < EdgeWithAttachmentTable[s_1].EdgeNum; s_2++)
		{
			if (DistanceArray[i] <= 0.0 && NodeTable[EdgeWithAttachmentTable[s_1].Pointer[s_2].n_1_ID - 1].commId != NodeTable[EdgeWithAttachmentTable[s_1].Pointer[s_2].n_2_ID - 1].commId)
			{
				void converge(int n1, int n2, double *d, int n);
				converge(EdgeWithAttachmentTable[s_1].Pointer[s_2].n_1_ID, EdgeWithAttachmentTable[s_1].Pointer[s_2].n_2_ID, DistanceArray, NumOfEdges);
			}
			else
			{
				if (DistanceArray[i] >= 1)
				{
					void divide(int n1, int n2, double *d, int n);
					divide(EdgeWithAttachmentTable[s_1].Pointer[s_2].n_1_ID, EdgeWithAttachmentTable[s_1].Pointer[s_2].n_2_ID, DistanceArray, NumOfEdges);
				}
			}
			++i;
		}
	}
}
void converge(int n1, int n2, double *d, int n)
{
	int small_index = (community[NodeTable[n1 - 1].commId].memberNum < community[NodeTable[n2 - 1].commId].memberNum) ? NodeTable[n1 - 1].commId : NodeTable[n2 - 1].commId;//成员较少的社团的社团标号
	int big_index = (small_index == NodeTable[n1 - 1].commId) ? NodeTable[n2 - 1].commId : NodeTable[n1 - 1].commId;//成员个数较多的社团的社团标号
	for (int i = 0; i < community[small_index].memberNum; ++i)
	{
		for (int j = 0; j < NodeTable[community[small_index].point_to_member[i] - 1].NodeNum; ++j)
		{
			if (NodeTable[NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j] - 1].commId == big_index)
			{
				int NodeMin = (community[small_index].point_to_member[i] < NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j] ? community[small_index].point_to_member[i] : NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j]);
				int NodeMax = (community[small_index].point_to_member[i] > NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j] ? community[small_index].point_to_member[i] : NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j]);
				int findIndex(int p, int q, EA *pointer, int num);
				int index = findIndex(NodeMin, NodeMax, EdgeWithAttachmentTable[NodeMin - 1].Pointer, EdgeWithAttachmentTable[NodeMin - 1].EdgeNum);
				if (index == EdgeWithAttachmentTable[NodeMin - 1].EdgeNum)
				{
					printf("error in converge communities(边不存在)\n");
					exit(1);
				}
				d[EdgeWithAttachmentTable[NodeMin - 1].Pointer[index].e_ID] = 0.0;
			}
		}
		NodeTable[community[small_index].point_to_member[i] - 1].commId = big_index;
	}
	int *temp = community[big_index].point_to_member;
	community[big_index].point_to_member = (int *)malloc(sizeof(int) * (community[big_index].memberNum + community[small_index].memberNum));
	for (int i = 0; i < community[big_index].memberNum; ++i)
	{
		community[big_index].point_to_member[i] = temp[i];
	}
	for (int j = 0; j < community[small_index].memberNum; ++j)
	{
		community[big_index].point_to_member[j + community[big_index].memberNum] = community[small_index].point_to_member[j];
	}
	community[big_index].memberNum = community[big_index].memberNum + community[small_index].memberNum;
	community[small_index].memberNum = 0;
	free(community[small_index].point_to_member);
	community[small_index].point_to_member = NULL;
	free(temp);
}

int findIndex(int n1, int n2, EA *pointer, int num)
{
	for (int i = 0; i < num; ++i)
	{
		if ((pointer[i].n_1_ID == n1 && pointer[i].n_2_ID == n2) || (pointer[i].n_1_ID == n2 && pointer[i].n_2_ID == n1))
			return i;
	}
	return num;
}

void divide(int n1, int n2, double *d, int n)
{
	int small_index = (community[NodeTable[n1 - 1].commId].memberNum < community[NodeTable[n2 - 1].commId].memberNum) ? NodeTable[n1 - 1].commId : NodeTable[n2 - 1].commId;
	int big_index = (small_index == NodeTable[n1 - 1].commId) ? NodeTable[n2 - 1].commId : NodeTable[n1 - 1].commId;
	for (int i = 0; i < community[small_index].memberNum; ++i)
	{
		for (int j = 0; j < NodeTable[community[small_index].point_to_member[i] - 1].NodeNum; ++j)
		{
			if (NodeTable[NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j] - 1].commId == big_index)
			{
				int NodeMin = (community[small_index].point_to_member[i] < NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j] ? community[small_index].point_to_member[i] : NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j]);
				int NodeMax = (community[small_index].point_to_member[i] > NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j] ? community[small_index].point_to_member[i] : NodeTable[community[small_index].point_to_member[i] - 1].Pointer[j]);
				int findIndex(int n1, int n2, EA *pointer, int num);
				int index = findIndex(NodeMin, NodeMax, EdgeWithAttachmentTable[NodeMin - 1].Pointer, EdgeWithAttachmentTable[NodeMin - 1].EdgeNum);
				if (index == EdgeWithAttachmentTable[NodeMin - 1].EdgeNum)
				{
					printf("error in divide communities(边不存在)\n");
					exit(1);
				}
				d[EdgeWithAttachmentTable[NodeMin - 1].Pointer[index].e_ID] = 1.0;
			}
		}
	}
}

int findEdge(EA *pointer, int num, int max)
{
	for (int i = 0; i < num; ++i)
	{
		if (pointer[i].n_1_ID == max || pointer[i].n_2_ID == max)
			return i;
	}
	return -1;
}

void findCN(int node_i, int node_j){

	int num1 = NodeTable[node_i - 1].NodeNum;
	int * A1 = NodeTable[node_i - 1].Pointer;
	int num2 = NodeTable[node_j - 1].NodeNum;
	int * A2 = NodeTable[node_j - 1].Pointer;
	int p1_loc = 0;
	int p2_loc = 0;
	int * p1 = &A1[0];
	int * p2 = &A2[0];
	int cn_length = 0;
	int diffa_length = 0;
	int diffb_length = 0;
	int cn_loc = 1;
	int diffa_loc = 1;
	int diffb_loc = 1;
	while (p1_loc < num1 && p2_loc < num2){
		if (p1[p1_loc] < p2[p2_loc])
			p1_loc++;
		else if (p1[p1_loc] > p2[p2_loc])
			p2_loc++;
		else
		{
			CN[cn_loc] = p1[p1_loc];
			cn_length++;
			cn_loc++;
			p1_loc++;
			p2_loc++;
		}
	}
	if (cn_loc != cn_length + 1)
	{
		printf("error in find common neighbors(IN FUNTION findCN)\n");
		exit(1);
	}
	CN[0] = cn_length;
}

void sortFun(int * Pointer, int Num){
	int i = Num - 1;
	int swap;
	while (i > 0){
		int LastChangeIndex = 0;
		for (int j = 0; j < i; j++)
		{
			if (Pointer[j] > Pointer[j + 1])
			{
				swap = Pointer[j + 1];
				Pointer[j + 1] = Pointer[j];
				Pointer[j] = swap;
				LastChangeIndex = j + 1;
			}
		}
		i = LastChangeIndex;
	}
}

int main(int argc, char * argv[])
{
	printf("begin........................\n");
	establishNodeTable(argv[1], atoi(argv[3]));
	establishEdgeTable();
	establishEdgesOfNodeTable();
	interaction(atoi(argv[3]));
	outputCommunities(argv[2]);
	printf("end..........................\n");
	return 0;
}
