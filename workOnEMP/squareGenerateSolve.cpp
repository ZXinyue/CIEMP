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
using namespace std;

#define POOL 1


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
	string s;
	int tryTime;
	bool d0;
	bool d1;
	bool d2;
	bool d3;
};

map<string,squareInfo> squareMap;
vector<string> squareVector;

int size = 25;
int width = 5;
Tile tiles[25];
map<string,int> edgeVote;
set<string> invalidEdges;
int* answer;

void* searchAnswer(void *id);

void* justTest(void *id);

void* searchPosition(int pos, int theSquareSize, vector<int> squareIds, set<int> &usedTiles, vector< set<int> > &triedPos, vector<int> &posAns,int direction,string originKey);

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
	clock_t start,end;
	start = clock();

	for(int i=0; i<size; i++)
	{
		string s;
		cin>>s;
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


	// for(int i=0; i<size-1; i++)
	// 	for(int j=i+1; j<size; j++)
	// 	{
	// 		Tile t1 = tiles[i];
	// 		Tile t2 = tiles[j];
	// 		string tempEdge;
	// 		if(t1.top == t2.bottom)
	// 		{
	// 			tempEdge = to_string(t2.id)+"T-B"+to_string(t1.id);
	// 			edgeVote[tempEdge] = 1;
	// 		}
	// 		if(t1.bottom == t2.top)
	// 		{
	// 			tempEdge = to_string(t1.id)+"T-B"+to_string(t2.id);
	// 			edgeVote[tempEdge] = 1;
	// 		}
	// 		if(t1.left == t2.right)
	// 		{
	// 			tempEdge = to_string(t2.id)+"L-R"+to_string(t1.id);
	// 			edgeVote[tempEdge] = 1;
	// 		}
	// 		if(t1.right == t2.left)
	// 		{
	// 			tempEdge = to_string(t1.id)+"L-R"+to_string(t2.id);
	// 			edgeVote[tempEdge] = 1;
	// 		}
	// 	}

	// map<string,int>::iterator iter;
	//     iter = edgeVote.begin();
	//     while(iter != edgeVote.end()) {
	//         cout << iter->first << " : " << iter->second << endl;
	//         iter++;
	//     }

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
	}

	for(int i=0;i<squareVector.size();i++)
		cout<<squareVector[i]<<endl;

	end = clock();

	cout<<"the total time is:"<<((double)(end-start))/CLOCKS_PER_SEC << endl;
	
	return 0;
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
	string randomStart = "";
	while(randomStart=="")
	{
		int randomNumber = rand()%squareVector.size();
		string vs = squareVector[randomNumber];
		if(squareMap[vs].tryTime!=1)
		{
			randomStart = vs;
			break;
		}
	}
	seed = rand();

	cout<<"start with:"<<randomStart<<endl;
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

	if(!squareMap[randomStart].d0)
	{
		searchPosition(0,theSqSize,squareIds,usedTiles,triedPos,posAns,0,randomStart);
	}
	// else if(!squareMap[randomStart].d1)
	// {
	// 	searchPosition(0,theSqSize,squareIds,usedTiles,triedPos,posAns,1,randomStart);
	// }
	// else if(!squareMap[randomStart].d2)
	// {
	// 	searchPosition(0,theSqSize,squareIds,usedTiles,triedPos,posAns,2,randomStart);
	// }
	// else
	// {
	// 	searchPosition(0,theSqSize,squareIds,usedTiles,triedPos,posAns,3,randomStart);
	// }
	
	return NULL;
}

void* searchPosition(int pos, int theSquareSize, vector<int> squareIds, set<int> &usedTiles, vector< set<int> > &triedPos, vector<int> &posAns,int direction,string originKey)
{
	cout<<"pos:"<<pos<<endl;
	cout<<"theSquareSize:"<<theSquareSize<<endl;
	for(int i=pos+1; i<2*theSquareSize+1; i++)
	{
		if(!triedPos[i].empty())
			triedPos[i].clear(); 
	}
	bool found = false;
	cout<<"here"<<endl;
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
					cout<<"put "<<i<<" at pos "<<pos<<endl;
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
				searchPosition(pos+1,theSquareSize,squareIds,usedTiles,triedPos,posAns,0,originKey);
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

				pthread_mutex_lock(&mutex1);
				squareInfo tempInfo = {key,0,false,false,false,false};
				squareMap[key] = tempInfo;
				squareVector.push_back(key);
				pthread_mutex_unlock(&mutex1);

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
				searchPosition(pos-1,theSquareSize,squareIds,usedTiles,triedPos,posAns,0,originKey);
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d0 = true;
				squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				return NULL;
			}
			
		}
	}
}

