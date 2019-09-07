#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <pthread.h>
#include <string.h>
#include <string>
#include <map>
#include <queue>
#include <vector>
#include <set>
#include <functional>
#include <fstream>
#include <sys/time.h>
using namespace std;

#define POOL 24

string inFileName = "data/9*9_3";
string outFileName = "result/9*9_1_24threads";

typedef struct
{
	int id;
	string top;
	string right;
	string bottom;
	string left;
	
}Tile;


struct position
{
	int x;
	int y;
	bool operator < (const position &p) const
	{
		if(x!=p.x)
			return x<p.x;
		else
			return y<p.y;
	}
};

typedef struct 
{
	int minRow;
	int minCol;
	int maxRow;
	int maxCol;
}boundaryInfo;

struct candidate
{	
	int id;
	int x;
	int y;
	int score;
	int relativeId;
	int relativeOrient;//0:T 1:R 2:B 3:L
	bool operator < (const candidate &c) const
	{
		return c.score < score ;
	}

};

struct deadEdge
{
	string edge;
	int score;
	bool operator < (const deadEdge &e) const
	{
		return score < e.score ;
	}
};

struct squareInfo
{
	string s; //square tiles ID
	int tryTime; 
	bool d0; 
	bool d1;
	bool d2;
	bool d3;
};

map<string,squareInfo> squareMap;
vector<string> squareVector;

map<string,map<int,int> > allEdges;
map<string,vector<int> > validEdges;

int size = 81;
int width = 9;
Tile tiles[81];
map<string,int> edgeVote;
set<string> invalidEdges;
int* answer;

void* searchAnswer(void *id);

void* justTest(void *id);

int* searchPosition(int pos, int theSquareSize, vector<int> squareIds, set<int> &usedTiles, vector< set<int> > &triedPos, vector<int> &posAns,int direction,string originKey);

//pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t mutex1;

void superSplit(const string& s, vector<string>& v, const string& c)
{
    string::size_type pos1, pos2;
    size_t len = s.length();
    pos2 = s.find(c);
    pos1 = 0;
    while(std::string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2-pos1));
 
        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if(pos1 != len)
        v.push_back(s.substr(pos1));
}

int main()
{
	//clock_t start,end;
	timeval t_start, t_end;
	//start = clock();
	gettimeofday( &t_start, NULL);

	ifstream input;
	input.open(inFileName);

	for(int i=0; i<size; i++)
	{
		string s;
		input>>s;
		tiles[i].id = i;
		tiles[i].top = s[0];
		tiles[i].right = s[1];
		tiles[i].bottom = s[2];
		tiles[i].left = s[3];

		string vs = to_string(i);
		squareInfo si = {s,0,false,false,false,false};
		squareMap[s] = si;
		squareVector.push_back(vs);
	}

	input.close();

    int allEdgeCount = 0;
	for(int i=0; i<size-1; i++)
		for(int j=i+1; j<size; j++)
		{
			Tile t1 = tiles[i];
			Tile t2 = tiles[j];
			string s1,s2;
			if(t1.top == t2.bottom)
			{
				s1 = to_string(t2.id)+"T-B";
				allEdges[s1][t1.id] = -1;
				s2 = to_string(t1.id)+"B-T";
				allEdges[s2][t2.id] = -1;
				allEdgeCount+=2;
			}
			if(t1.bottom == t2.top)
			{
				s1 = to_string(t1.id)+"T-B";
				allEdges[s1][t2.id] = -1;
				s2 = to_string(t2.id)+"B-T";
				allEdges[s2][t1.id] = -1;
				allEdgeCount+=2;
			}
			if(t1.left == t2.right)
			{
				s1 = to_string(t2.id)+"L-R";
				allEdges[s1][t1.id] = -1;
				s2 = to_string(t1.id)+"R-L";
				allEdges[s2][t2.id] = -1;
				allEdgeCount+=2;
			}
			if(t1.right == t2.left)
			{
				s1 = to_string(t1.id)+"L-R";
				allEdges[s1][t2.id] = -1;
				s2 = to_string(t2.id)+"R-L";
				allEdges[s2][t1.id] = -1;
				allEdgeCount+=2;
			}
		}
	cout<<"allEdgeCount:"<<allEdgeCount<<endl;

	map<string,map<int,int> > :: iterator it1;
	it1 = allEdges.begin();
	while(it1 != allEdges.end())
	{
		string tempS = it1->first;
		map<int,int> tempMap = it1->second;
		map<int,int> :: iterator it2;
		it2 = tempMap.begin();
		while(it2 != tempMap.end())
		{
			if(it2->second != -1)
			{
				it2++;
				continue;
			}
			string::size_type findResultPos;
			findResultPos = tempS.find("L-R");
			if(findResultPos != tempS.npos) //L-R edge
			{
				int theLeftId = stoi(tempS.substr(0,findResultPos));
				int theRightId = it2->first;
				bool certiValid = false;

				//top direction
				string s1 = to_string(theLeftId)+"B-T";
				string s2 = to_string(theRightId)+"B-T";
				string tempS1;
				map<int,int> topLeftValidMap = allEdges[s1];
				map<int,int> topRightValidMap = allEdges[s2];
				map<int,int>::iterator tlvm, trvm;
				tlvm = topLeftValidMap.begin();
				trvm = topRightValidMap.begin();
				while(tlvm != topLeftValidMap.end())
				{
					int topLeftId = tlvm->first;
					while(trvm!= topRightValidMap.end())
					{
						int topRightId = trvm->first;
						if(tiles[topLeftId].right == tiles[topRightId].left)
						{
							certiValid = true;
							allEdges[tempS][theRightId] = 1;
							tempS1 = to_string(theRightId)+"R-L";
							allEdges[tempS1][theLeftId] = 1;
							tempS1 = to_string(topLeftId)+"L-R";
							allEdges[tempS1][topRightId] = 1;
							tempS1 = to_string(topRightId)+"R-L";
							allEdges[tempS1][topLeftId] = 1;
							tempS1 = to_string(topLeftId)+"T-B";
							allEdges[tempS1][theLeftId] = 1;
							tempS1 = to_string(theLeftId)+"B-T";
							allEdges[tempS1][topLeftId] = 1;
							tempS1 = to_string(topRightId)+"T-B";
							allEdges[tempS1][theRightId] = 1;
							tempS1 = to_string(theRightId)+"B-T";
							allEdges[tempS1][topRightId] = 1;

						}
						trvm++;
					}
					tlvm++;
				}

				//bottom direction
				s1 = to_string(theLeftId)+"T-B";
				s2 = to_string(theRightId)+"T-B";
				map<int,int> bottomLeftValidMap = allEdges[s1];
				map<int,int> bottomRightValidMap = allEdges[s2];
				map<int,int>::iterator blvm, brvm;
				blvm = bottomLeftValidMap.begin();
				brvm = bottomRightValidMap.begin();
				while(blvm != bottomLeftValidMap.end())
				{
					int bottomLeftId = blvm->first;
					while(brvm!= bottomRightValidMap.end())
					{
						int bottomRightId = brvm->first;
						if(tiles[bottomLeftId].right == tiles[bottomRightId].left)
						{
							certiValid = true;
							allEdges[tempS][theRightId] = 1;
							tempS1 = to_string(theRightId)+"R-L";
							allEdges[tempS1][theLeftId] = 1;
							tempS1 = to_string(bottomLeftId)+"L-R";
							allEdges[tempS1][bottomRightId] = 1;
							tempS1 = to_string(bottomRightId)+"R-L";
							allEdges[tempS1][bottomLeftId] = 1;
							tempS1 = to_string(bottomLeftId)+"B-T";
							allEdges[tempS1][theLeftId] = 1;
							tempS1 = to_string(theLeftId)+"T-B";
							allEdges[tempS1][bottomLeftId] = 1;
							tempS1 = to_string(bottomRightId)+"B-T";
							allEdges[tempS1][theRightId] = 1;
							tempS1 = to_string(theRightId)+"T-B";
							allEdges[tempS1][bottomRightId] = 1;
						}
						brvm++;
					}
					blvm++;
				}
				if(!certiValid)
				{
					allEdges[tempS][theRightId] = 0;
					tempS1 = to_string(theRightId)+"R-L";
					allEdges[tempS1][theLeftId] = 0;

				}
			}
			findResultPos = tempS.find("T-B");
			if(findResultPos != tempS.npos) //T-B edge
			{
				int theTopId = stoi(tempS.substr(0,findResultPos));
				int theBottomId = it2->first;
				bool certiValid = false;

				//right direction
				string s1 = to_string(theTopId)+"L-R";
				string s2 = to_string(theBottomId)+"L-R";
				string tempS1;
				map<int,int> rightTopValidMap = allEdges[s1];
				map<int,int> rightBottomValidMap = allEdges[s2];
				map<int,int>::iterator rtvm, rbvm;
				rtvm = rightTopValidMap.begin();
				rbvm = rightBottomValidMap.begin();
				while(rtvm != rightTopValidMap.end())
				{
					int rightTopId = rtvm->first;
					while(rbvm!= rightBottomValidMap.end())
					{
						int rightBottomId = rbvm->first;
						if(tiles[rightTopId].bottom == tiles[rightBottomId].top)
						{
							certiValid = true;
							allEdges[tempS][theBottomId] = 1;
							tempS1 = to_string(theBottomId)+"B-T";
							allEdges[tempS1][theTopId] = 1;
							tempS1 = to_string(rightTopId)+"T-B";
							allEdges[tempS1][rightBottomId] = 1;
							tempS1 = to_string(rightBottomId)+"B-T";
							allEdges[tempS1][rightTopId] = 1;
							tempS1 = to_string(theTopId)+"L-R";
							allEdges[tempS1][rightTopId] = 1;
							tempS1 = to_string(rightTopId)+"R-L";
							allEdges[tempS1][theTopId] = 1;
							tempS1 = to_string(theBottomId)+"L-R";
							allEdges[tempS1][rightBottomId] = 1;
							tempS1 = to_string(rightBottomId)+"R-L";
							allEdges[tempS1][theBottomId] = 1;
						}
						rbvm++;
					}
					rtvm++;
				}

				//left direction
				s1 = to_string(theTopId)+"R-L";
				s2 = to_string(theBottomId)+"R-L";
				map<int,int> leftTopValidMap = allEdges[s1];
				map<int,int> leftBottomValidMap = allEdges[s2];
				map<int,int>::iterator ltvm, lbvm;
				ltvm = leftTopValidMap.begin();
				lbvm = leftBottomValidMap.begin();
				while(ltvm != leftTopValidMap.end())
				{
					int leftTopId = ltvm->first;
					while(lbvm!= leftBottomValidMap.end())
					{
						int leftBottomId = lbvm->first;
						if(tiles[leftTopId].bottom == tiles[leftBottomId].top)
						{
							certiValid = true;
							allEdges[tempS][theBottomId] = 1;
							tempS1 = to_string(theBottomId)+"B-T";
							allEdges[tempS1][theTopId] = 1;
							tempS1 = to_string(leftTopId)+"T-B";
							allEdges[tempS1][leftBottomId] = 1;
							tempS1 = to_string(leftBottomId)+"B-T";
							allEdges[tempS1][leftTopId] = 1;
							tempS1 = to_string(theTopId)+"R-L";
							allEdges[tempS1][leftTopId] = 1;
							tempS1 = to_string(leftTopId)+"L-R";
							allEdges[tempS1][theTopId] = 1;
							tempS1 = to_string(theBottomId)+"R-L";
							allEdges[tempS1][leftBottomId] = 1;
							tempS1 = to_string(leftBottomId)+"L-R";
							allEdges[tempS1][theBottomId] = 1;
						}
						lbvm++;
					}
					ltvm++;
				}
				if(!certiValid)
				{
					allEdges[tempS][theBottomId] = 0;
					tempS1 = to_string(theBottomId)+"B-T";
					allEdges[tempS1][theTopId] = 0;

				}
			}
			it2++;
		}
		it1++;
	}

	int validEdgeCount = 0;
	it1 = allEdges.begin();
	while(it1 != allEdges.end())
	{
		string tempS = it1->first;
		map<int,int> tempMap = it1->second;
		map<int,int> :: iterator it2;
		it2 = tempMap.begin();
		vector<int> tempVec;
		while(it2 != tempMap.end())
		{
			if(it2->second != 0)
			{
				tempVec.push_back(it2->first);
				validEdgeCount+=1;
			}
			it2++;
		}
		validEdges[tempS] = tempVec;
		it1++;
	}

	cout<<"validEdges:"<<validEdgeCount<<endl;
	return 0;


	// map<string,int>::iterator iter;
	//     iter = edgeVote.begin();
	//     while(iter != edgeVote.end()) {
	//         cout << iter->first << " : " << iter->second << endl;
	//         iter++;
	//     }

/*
	void *status;
	
	pthread_t *threads;
	threads = (pthread_t *)malloc(sizeof(pthread_t)*POOL);
	int indexes[POOL];

	answer = (int *)malloc(sizeof(int)*size);

	for(int j=0; j<size; j++)
	{
		answer[j] = -1;
	}

	pthread_mutex_init(&mutex1,NULL);

	int tid;

	for (tid = 0; tid < POOL; tid++)
	{
		indexes[tid] = tid;
		//int* pi = &tid;
		if (pthread_create(&threads[tid], NULL, searchAnswer, (void*)&indexes[tid]))
        {
                exit(-1);
        }
	}
	
	int i;
	for (i = 0; i < POOL; i++)
	{
		if (pthread_join(threads[i], &status))
		{
			exit(-1);
		}
	}
	pthread_mutex_destroy(&mutex1);

	if(answer[0]!=-1)
	{
		cout<<"the answer is:"<<endl;
		for(int j=0; j<size; j++)
		{
			cout<<answer[j]<<" ";
		}
		cout<<endl;
	}

	// for(int i=0;i<squareVector.size();i++)
	// 	cout<<squareVector[i]<<endl;

	//end = clock();
	gettimeofday( &t_end, NULL);

	double delta_t = (t_end.tv_sec-t_start.tv_sec) + (t_end.tv_usec-t_start.tv_usec)/1000000.0;

	//cout<<"the total time is:"<<((double)(end-start))/CLOCKS_PER_SEC << endl;
	cout<<"the total time is:"<<delta_t<<endl;
	ofstream output;
	
	output.open(outFileName,ios::app);
	output << delta_t<< endl;
	output.close();

	free(threads);
	free(answer);

	return 0;
	*/
}

void* justTest(void *id)
{
	cout<<"test thread id:"<<*(int*)id<<endl;
	return NULL;
}

void* searchAnswer(void *id)
{
	long seed = time(NULL);
	srand(seed+ *(int*)id * 12345);
	seed = rand();
	while(answer[0]==-1)
	{
		string randomStart = "";
		while(randomStart=="")
		{
			srand(seed);
			int randomNumber = rand()%squareVector.size();
			seed = rand();
			string vs = squareVector[randomNumber];
			if(squareMap[vs].tryTime!=4)
			{
				randomStart = vs;
				break;
			}
		}		
		//cout<<"start with:"<<randomStart<<endl;
		vector<string> splitResult;
		set<int> usedTiles;
		vector<int> squareIds;
		superSplit(randomStart,splitResult,"-");

		for(int i=0;i<splitResult.size();i++)
		{
			//cout<<"superSplit:"<<splitResult[i]<<endl;
			usedTiles.insert(stoi(splitResult[i]));
			squareIds.push_back(stoi(splitResult[i]));
		}

		int theSqSize = sqrt(squareIds.size());

		vector< set<int> > triedPos;
		vector<int> posAns;
		for(int i=0; i<theSqSize*2+1; i++)
		{
			set<int> empty;
			triedPos.push_back(empty);
			posAns.push_back(-1);
		}

		int finish = 0;
		int POS = 0;
		if(!squareMap[randomStart].d0)
		{
			
			while(finish!=2)
			{
				int *pfinish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,0,randomStart);
				finish = *pfinish;
				delete pfinish;
				if(finish == 1)
					POS+=1;
				else if(finish == 0)
					POS-=1;
				else
					break;
			}
			
		}
		else if(!squareMap[randomStart].d1)
		{
			while(finish!=2)
			{
				int *pfinish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,1,randomStart);
				finish = *pfinish;
				delete pfinish;
				if(finish == 1)
					POS+=1;
				else if(finish == 0)
					POS-=1;
				else
					break;
			}
		}
		else if(!squareMap[randomStart].d2)
		{
			while(finish!=2)
			{
				int *pfinish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,2,randomStart);
				finish = *pfinish;
				delete pfinish;
				if(finish == 1)
					POS+=1;
				else if(finish == 0)
					POS-=1;
				else
					break;
			}
		}
		else
		{
			while(finish!=2)
			{
				int *pfinish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,3,randomStart);
				finish = *pfinish;
				delete pfinish;
				if(finish == 1)
					POS+=1;
				else if(finish == 0)
					POS-=1;
				else
					break;
			}
		}
	}
	
	return NULL;
}
//0:failed 1:success 2:finish
int* searchPosition(int pos, int theSquareSize, vector<int> squareIds, set<int> &usedTiles, vector< set<int> > &triedPos, vector<int> &posAns,int direction,string originKey)
{
	// int returnStatus;
	// int *pRe;
	//cout<<"pos:"<<pos<<endl;
	//cout<<"theSquareSize:"<<theSquareSize<<endl;
	for(int i=pos+1; i<2*theSquareSize+1; i++)
	{
		if(!triedPos[i].empty())
			triedPos[i].clear(); 
	}
	bool found = false;
	//cout<<"here"<<endl;
	if(direction == 0)
	{
		if(pos == 0)
		{
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].bottom == tiles[squareIds[pos]].top)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					//cout<<"put "<<i<<" at pos "<<pos<<endl;
					break;
				}
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int leftId = posAns[pos-1];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].bottom == tiles[squareIds[pos]].top && tiles[i].left == tiles[leftId].right)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}

		}
		else if(pos == theSquareSize)
		{
			int leftId = posAns[pos-1];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].left == tiles[leftId].right)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int topId = posAns[pos-1];
			int relativeId = (pos-theSquareSize-1)*theSquareSize+theSquareSize-1;
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].top == tiles[topId].bottom && tiles[i].left == tiles[squareIds[relativeId]].right)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}


		if(found)
		{
			
			if(pos<2*theSquareSize)
			{
				//returnStatus = 1;
				int *pRe = new int(1);
				return pRe;
				//searchPosition(pos+1,theSquareSize,squareIds,usedTiles,triedPos,posAns,0,originKey);
			}
			else if(pos == 2*theSquareSize)
			{
				vector<int> newSquare;
				int k;
				for(k=0;k<theSquareSize+1;k++)
					newSquare.push_back(posAns[k]);
				for(int i=0;i<theSquareSize;i++)
				{
					for(int j=0;j<theSquareSize;j++)
					{
						newSquare.push_back(squareIds[i*theSquareSize+j]);
					}
					newSquare.push_back(posAns[k++]);
				}
				string key = "";
				for(int i=0; i<newSquare.size()-1; i++)
					key = key + to_string(newSquare[i]) + "-";
				key+=to_string(newSquare[newSquare.size()-1]);

				if(squareMap.find(key)==squareMap.end())
				{
					pthread_mutex_lock(&mutex1);
					squareInfo tempInfo = {key,0,false,false,false,false};
					squareMap[key] = tempInfo;
					squareVector.push_back(key);
					pthread_mutex_unlock(&mutex1);
				}

				

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex1);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex1);
					//returnStatus = 2;
					int *pRe = new int(2);
					return pRe;
					// return NULL;
				}

				found = false;
				pos+=1;
			}
			
		}
		if(!found)
		{
			if(pos!=0)
			{
				int preId = posAns[pos-1];
				usedTiles.erase(preId);
				posAns[pos-1]=-1;
				//returnStatus = 0;
				int *pRe = new int(0);
				return pRe;
				// searchPosition(pos-1,theSquareSize,squareIds,usedTiles,triedPos,posAns,0,originKey);
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d0 = true;
				squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				int *pRe = new int(2);
				return pRe;
			}
			
		}
	}
	else if(direction == 1)
	{
		if(pos == 0)
		{
			int relativeId = squareIds[pos*theSquareSize+theSquareSize-1];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].left == tiles[relativeId].right)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					//cout<<"put "<<i<<" at pos "<<pos<<endl;
					break;
				}
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int topId = posAns[pos-1];
			int relativeId = squareIds[pos*theSquareSize+theSquareSize-1];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].left == tiles[relativeId].right && tiles[i].top == tiles[topId].bottom)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}

		}
		else if(pos == theSquareSize)
		{
			int topId = posAns[pos-1];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].top == tiles[topId].bottom)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int rightId = posAns[pos-1];
			int relativeId = (theSquareSize-1)*theSquareSize+theSquareSize-1-(pos-theSquareSize-1);
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].right == tiles[rightId].left && tiles[i].top == tiles[squareIds[relativeId]].bottom)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}

		if(found)
		{
			
			if(pos<2*theSquareSize)
			{
				//returnStatus = 1;
				int *pRe = new int(1);
				return pRe;
				//searchPosition(pos+1,theSquareSize,squareIds,usedTiles,triedPos,posAns,1,originKey);
			}
			else if(pos == 2*theSquareSize)
			{
				vector<int> newSquare;
				int k=0;
				for(int i=0;i<theSquareSize;i++)
				{
					for(int j=0;j<theSquareSize;j++)
					{
						newSquare.push_back(squareIds[i*theSquareSize+j]);
					}
					newSquare.push_back(posAns[k++]);
				}
				for(int i=posAns.size()-1; i>=k; i--)
					newSquare.push_back(posAns[i]);
				string key = "";
				for(int i=0; i<newSquare.size()-1; i++)
					key = key + to_string(newSquare[i]) + "-";
				key+=to_string(newSquare[newSquare.size()-1]);

				if(squareMap.find(key)==squareMap.end())
				{
					pthread_mutex_lock(&mutex1);
					squareInfo tempInfo = {key,0,false,false,false,false};
					squareMap[key] = tempInfo;
					squareVector.push_back(key);
					pthread_mutex_unlock(&mutex1);
				}

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex1);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex1);
					//return NULL;
					//returnStatus = 2;
					int *pRe = new int(2);
					return pRe;
				}

				found = false;
				pos+=1;
			}
			
		}
		if(!found)
		{
			if(pos!=0)
			{
				int preId = posAns[pos-1];
				usedTiles.erase(preId);
				posAns[pos-1]=-1;
				//returnStatus = 0;
				int *pRe = new int(0);
				return pRe;
				//searchPosition(pos-1,theSquareSize,squareIds,usedTiles,triedPos,posAns,1,originKey);
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d1 = true;
				squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				int *pRe = new int(2);
				return pRe;
			}
			
		}
	}
	else if(direction == 2) //|_>
	{
		if(pos == 0)
		{
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].right == tiles[pos].left)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					//cout<<"put "<<i<<" at pos "<<pos<<endl;
					break;
				}
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int topId = posAns[pos-1];
			int relativeId = squareIds[pos*theSquareSize];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].right == tiles[relativeId].left && tiles[i].top == tiles[topId].bottom)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}

		}
		else if(pos == theSquareSize)
		{
			int topId = posAns[pos-1];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].top == tiles[topId].bottom)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int leftId = posAns[pos-1];
			int relativeId = (theSquareSize-1)*theSquareSize+pos-theSquareSize-1;
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].left == tiles[leftId].right && tiles[i].top == tiles[squareIds[relativeId]].bottom)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}
		if(found)
		{
			
			if(pos<2*theSquareSize)
			{
				//returnStatus = 1;
				int *pRe = new int(1);
				return pRe;
				//searchPosition(pos+1,theSquareSize,squareIds,usedTiles,triedPos,posAns,2,originKey);
			}
			else if(pos == 2*theSquareSize)
			{
				vector<int> newSquare;
				int k=0;
				for(int i=0;i<theSquareSize;i++)
				{
					newSquare.push_back(posAns[k++]);
					for(int j=0;j<theSquareSize;j++)
					{
						newSquare.push_back(squareIds[i*theSquareSize+j]);
					}
					
				}
				while(k<posAns.size())
					newSquare.push_back(posAns[k++]);
				string key = "";
				for(int i=0; i<newSquare.size()-1; i++)
					key = key + to_string(newSquare[i]) + "-";
				key+=to_string(newSquare[newSquare.size()-1]);

				if(squareMap.find(key)==squareMap.end())
				{
					pthread_mutex_lock(&mutex1);
					squareInfo tempInfo = {key,0,false,false,false,false};
					squareMap[key] = tempInfo;
					squareVector.push_back(key);
					pthread_mutex_unlock(&mutex1);
				}

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex1);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex1);
					//return NULL;
					//returnStatus = 2;
					int *pRe = new int(2);
					return pRe;
				}

				found = false;
				pos+=1;
			}
			
		}
		if(!found)
		{
			if(pos!=0)
			{
				int preId = posAns[pos-1];
				usedTiles.erase(preId);
				posAns[pos-1]=-1;
				//searchPosition(pos-1,theSquareSize,squareIds,usedTiles,triedPos,posAns,2,originKey);
				//returnStatus = 0;
				int *pRe = new int(0);
				return pRe;
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d2 = true;
				squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				int *pRe = new int(2);
				return pRe;
			}
			
		}

	}
	else if(direction == 3) //|<-
	{
		if(pos == 0)
		{
			int relativeId = squareIds[theSquareSize-1-pos];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].bottom == tiles[relativeId].top)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					//cout<<"put "<<i<<" at pos "<<pos<<endl;
					break;
				}
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int rightId = posAns[pos-1];
			int relativeId = squareIds[theSquareSize-1-pos];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].bottom == tiles[relativeId].top && tiles[i].right == tiles[rightId].left)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}

		}
		else if(pos == theSquareSize)
		{
			int rightId = posAns[pos-1];
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].right == tiles[rightId].left)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int topId = posAns[pos-1];
			int relativeId = (pos-theSquareSize-1)*theSquareSize;
			for(int i=0; i<size; i++)
			{
				
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				if(tiles[i].top == tiles[topId].bottom && tiles[i].right == tiles[squareIds[relativeId]].left)
				{
					found = true;
					posAns[pos] = i;
					usedTiles.insert(i);
					triedPos[pos].insert(i);
					break;
				}
			}
		}

		if(found)
		{
			
			if(pos<2*theSquareSize)
			{
				//returnStatus = 1;
				int *pRe = new int(1);
				return pRe;
				//searchPosition(pos+1,theSquareSize,squareIds,usedTiles,triedPos,posAns,3,originKey);
			}
			else if(pos == 2*theSquareSize)
			{
				vector<int> newSquare;
				for(int i=theSquareSize;i>=0;i--)
					newSquare.push_back(posAns[i]);
				int k=theSquareSize+1;
				for(int i=0;i<theSquareSize;i++)
				{
					newSquare.push_back(posAns[k++]);
					for(int j=0;j<theSquareSize;j++)
					{
						newSquare.push_back(squareIds[i*theSquareSize+j]);
					}
					
				}
				
				string key = "";
				for(int i=0; i<newSquare.size()-1; i++)
					key = key + to_string(newSquare[i]) + "-";
				key+=to_string(newSquare[newSquare.size()-1]);

				if(squareMap.find(key)==squareMap.end())
				{
					pthread_mutex_lock(&mutex1);
					squareInfo tempInfo = {key,0,false,false,false,false};
					squareMap[key] = tempInfo;
					squareVector.push_back(key);
					pthread_mutex_unlock(&mutex1);
				}

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex1);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex1);
					//return NULL;
					//returnStatus = 2;
					int *pRe = new int(2);
					return pRe;
				}

				found = false;
				pos+=1;
			}
			
		}
		if(!found)
		{
			if(pos!=0)
			{
				int preId = posAns[pos-1];
				usedTiles.erase(preId);
				posAns[pos-1]=-1;
				//searchPosition(pos-1,theSquareSize,squareIds,usedTiles,triedPos,posAns,3,originKey);
				//returnStatus = 0;
				int *pRe = new int(0);
				return pRe;
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d3 = true;
				squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				int *pRe = new int(2);
				return pRe;
			}
			
		}

	}
}

