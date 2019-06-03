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

int size = 9;
int width = 3;
Tile tiles[9];
map<string,int> edgeVote;
int* answer;

void* searchAnswer(void *id);

void* justTest(void *id);

void* putTile(map<int,position> &tempAnswer, position p, int tileId, set<int> &usedTiles, set<position> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary);

void* updateCandidates(map<int,position> &tempAnswer,  position p, int tileId ,set<int> &usedTiles, set<position> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary);


//pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t mutex1;

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
	}


	for(int i=0; i<size-1; i++)
		for(int j=i+1; j<size; j++)
		{
			Tile t1 = tiles[i];
			Tile t2 = tiles[j];
			string tempEdge;
			if(t1.top == t2.bottom)
			{
				tempEdge = to_string(t2.id)+"T-B"+to_string(t1.id);
				edgeVote[tempEdge] = 1;
			}
			if(t1.bottom == t2.top)
			{
				tempEdge = to_string(t1.id)+"T-B"+to_string(t2.id);
				edgeVote[tempEdge] = 1;
			}
			if(t1.left == t2.right)
			{
				tempEdge = to_string(t2.id)+"L-R"+to_string(t1.id);
				edgeVote[tempEdge] = 1;
			}
			if(t1.right == t2.left)
			{
				tempEdge = to_string(t1.id)+"L-R"+to_string(t2.id);
				edgeVote[tempEdge] = 1;
			}
		}


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
	map<int,position> tempAnswer;
	set<int> usedTiles;
	set<position> usedPos;
	priority_queue<candidate> candidates;
	boundaryInfo boundary = {0,0,0,0};

	long seed = time(NULL);
	srand(seed+ *(int*)id * 12345);
	int randomStart = rand()%size;

	position p = {0,0};

	putTile(tempAnswer,p,randomStart,usedTiles,usedPos,candidates,boundary);

	while(!candidates.empty())
	{
		candidate winner = candidates.top();
		candidates.pop();
		cout<<"winner id: "<<winner.id<<" x:"<<winner.x <<" y:"<<winner.y<<endl;
		position winnerPos;
		winnerPos.x = winner.x;
		winnerPos.y = winner.y;
		int winnerId = winner.id;
		if(usedPos.count(winnerPos) == 1)
			continue;
		if(usedTiles.count(winnerId) == 1)
			continue;
		putTile(tempAnswer,winnerPos,winnerId,usedTiles,usedPos,candidates,boundary);
		int relativeId = winner.relativeId;
		int orient = winner.relativeOrient;
		string key = "";
		switch(orient)
		{
			case 0: key = to_string(relativeId)+"T-B"+to_string(winnerId);break;
			case 1: key = to_string(winnerId)+"L-R"+to_string(relativeId);break;
			case 2: key = to_string(winnerId)+"T-B"+to_string(relativeId);break;
			case 3: key = to_string(relativeId)+"L-R"+to_string(winnerId);break;
		};
		usedPos.insert(winnerPos);
		usedTiles.insert(winnerId);
		pthread_mutex_lock(&mutex1);
		edgeVote[key]+=1;
		pthread_mutex_unlock(&mutex1);
	}


	if(tempAnswer.size() == size)
	{
		int* result;
		result = (int *)malloc(sizeof(int)*size);
		map<int,position>::iterator iter;
		iter = tempAnswer.begin();
		while(iter != tempAnswer.end())
		{
			int id = iter->first;
			position pos = iter->second;
			result[(pos.x-boundary.minCol)*width + (pos.y-boundary.minRow)] = id;
			iter++;
		}
		pthread_mutex_lock(&mutex1);
		answer = result;
		pthread_mutex_unlock(&mutex1);
	}
	else
	{
		cout<<"thread "<<*(int*)id<<" faided,the result size:"<<tempAnswer.size()<<endl;
	}
	
	return NULL;
}

void* putTile(map<int,position> &tempAnswer, position p, int tileId, set<int> &usedTiles, set<position> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary)
{
	tempAnswer[tileId]  = p;
	usedTiles.insert(tileId);
	usedPos.insert(p);
	updateCandidates(tempAnswer, p, tileId, usedTiles, usedPos, candidates, boundary);
	return NULL;

}

void* updateCandidates(map<int,position> &tempAnswer,  position p, int tileId ,set<int> &usedTiles, set<position> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary)
{
	if(tempAnswer.size()!=size)
	{
		int x = p.x;
		int y = p.y;
		//top
		position top = {x-1,y};
		if(!usedPos.count(top))
		{
			if((abs(min(boundary.minCol,top.x))+abs(max(boundary.maxCol,top.x)) < width) && (abs(min(boundary.minRow,top.y))+abs(max(boundary.maxRow,top.y)) < width) )
			{
				string shead = "T-B"+to_string(tileId);
				map<string,int>::iterator iter;
			    iter = edgeVote.begin();
			    while(iter != edgeVote.end()) 
			    {
			    	string::size_type findResultPos;
			    	findResultPos = iter->first.find(shead);
			        if(findResultPos + shead.length() == iter->first.length())
			        {
			        	int sscore = iter->second;
			        	int candidateId = stoi(iter->first.substr(0,findResultPos));
			        	candidate can;
			        	can.id = candidateId;
			        	can.x = top.x;
			        	can.y = top.y;
			        	can.score = sscore;
			        	can.relativeId = tileId;
			        	can.relativeOrient = 2; //bottom
			        	candidates.push(can);
			        	boundary.minCol = min(top.x,boundary.minCol);
			        	boundary.minRow = min(top.y,boundary.minRow);
			        	boundary.maxCol = max(top.x,boundary.maxCol);
			        	boundary.maxRow = max(top.y,boundary.maxRow);
			        }
			        iter++;
			    }
				
			}
		}
		//right
		position right = {x,y+1};
		if(!usedPos.count(right))
		{
			if((abs(min(boundary.minCol,right.x))+abs(max(boundary.maxCol,right.x)) < width) && (abs(min(boundary.minRow,right.y))+abs(max(boundary.maxRow,right.y)) < width) )
			{
				string shead = to_string(tileId)+"L-R";
				map<string,int>::iterator iter;
			    iter = edgeVote.begin();
			    while(iter != edgeVote.end()) 
			    {
			    	string::size_type findResultPos;
			    	findResultPos = iter->first.find(shead);
			        if(findResultPos == 0)
			        {
			        	int sscore = iter->second;
			        	int candidateId = stoi(iter->first.substr(shead.length()));
			        	candidate can;
			        	can.id = candidateId;
			        	can.x = right.x;
			        	can.y = right.y;
			        	can.score = sscore;
			        	can.relativeId = tileId;
			        	can.relativeOrient = 3; //left
			        	candidates.push(can);
			        	boundary.minCol = min(right.x,boundary.minCol);
			        	boundary.minRow = min(right.y,boundary.minRow);
			        	boundary.maxCol = max(right.x,boundary.maxCol);
			        	boundary.maxRow = max(right.y,boundary.maxRow);

			        }
			        iter++;
			    }
				
			}
		}
		//bottom
		position bottom = {x+1,y};
		if(!usedPos.count(bottom))
		{
			if((abs(min(boundary.minCol,bottom.x))+abs(max(boundary.maxCol,bottom.x)) < width) && (abs(min(boundary.minRow,bottom.y))+abs(max(boundary.maxRow,bottom.y)) < width) )
			{
				string shead = to_string(tileId)+"T-B";
				map<string,int>::iterator iter;
			    iter = edgeVote.begin();
			    while(iter != edgeVote.end()) 
			    {
			    	string::size_type findResultPos;
			    	findResultPos = iter->first.find(shead);
			        if(findResultPos == 0)
			        {
			        	int sscore = iter->second;
			        	int candidateId = stoi(iter->first.substr(shead.length()));
			        	candidate can;
			        	can.id = candidateId;
			        	can.x = right.x;
			        	can.y = right.y;
			        	can.score = sscore;
			        	can.relativeId = tileId;
			        	can.relativeOrient = 0; //top
			        	candidates.push(can);
			        	boundary.minCol = min(bottom.x,boundary.minCol);
			        	boundary.minRow = min(bottom.y,boundary.minRow);
			        	boundary.maxCol = max(bottom.x,boundary.maxCol);
			        	boundary.maxRow = max(bottom.y,boundary.maxRow);
			        }
			         iter++;
			    }	
			}
		}
		//left
		position left = {x,y-1};
		if(!usedPos.count(left))
		{
			if((abs(min(boundary.minCol,left.x))+abs(max(boundary.maxCol,left.x)) < width) && (abs(min(boundary.minRow,left.y))+abs(max(boundary.maxRow,left.y)) < width) )
			{
				string shead = "L-R"+to_string(tileId);
				map<string,int>::iterator iter;
			    iter = edgeVote.begin();
			    while(iter != edgeVote.end()) 
			    {
			    	string::size_type findResultPos;
			    	findResultPos = iter->first.find(shead);
			        if(findResultPos + shead.length() == iter->first.length())
			        {
			        	int sscore = iter->second;
			        	int candidateId = stoi(iter->first.substr(0,findResultPos));
			        	candidate can;
			        	can.id = candidateId;
			        	can.x = top.x;
			        	can.y = top.y;
			        	can.score = sscore;
			        	can.relativeId = tileId;
			        	can.relativeOrient = 1; //right
			        	candidates.push(can);
			        	boundary.minCol = min(left.x,boundary.minCol);
			        	boundary.minRow = min(left.y,boundary.minRow);
			        	boundary.maxCol = max(left.x,boundary.maxCol);
			        	boundary.maxRow = max(left.y,boundary.maxRow);
			        }
			        iter++;
			    }	
			}
		}
	}
	return NULL;
}

