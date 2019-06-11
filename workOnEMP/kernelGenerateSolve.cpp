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

#define POOL 40


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

int size = 25;
int width = 5;
Tile tiles[25];
map<string,int> edgeVote;
map<string,int> edgeDead;
int* answer;

void* searchAnswer(void *id);

void* justTest(void *id);

void* putTile(map<int,position> &tempAnswer, position p, int tileId, set<int> &usedTiles, map<position,int> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary);

void* updateCandidates(map<int,position> &tempAnswer,  position p, int tileId ,set<int> &usedTiles, map<position,int> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary);

int* calculateValid(int candidateId, position pos, map<int,position> &tempAnswer, map<position,int> &usedPos);

bool* validAndEdge(int candidateId, position pos, map<int,position> &tempAnswer, map<position,int> &usedPos, string &t, string &r, string &b, string &l);

void* updateBoundary(map<int,position> &tempAnswer, boundaryInfo &boundary);
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
				edgeDead[tempEdge] = 0;
			}
			if(t1.bottom == t2.top)
			{
				tempEdge = to_string(t1.id)+"T-B"+to_string(t2.id);
				edgeVote[tempEdge] = 1;
				edgeDead[tempEdge] = 0;
			}
			if(t1.left == t2.right)
			{
				tempEdge = to_string(t2.id)+"L-R"+to_string(t1.id);
				edgeVote[tempEdge] = 1;
				edgeDead[tempEdge] = 0;
			}
			if(t1.right == t2.left)
			{
				tempEdge = to_string(t1.id)+"L-R"+to_string(t2.id);
				edgeVote[tempEdge] = 1;
				edgeDead[tempEdge] = 0;
			}
		}

	map<string,int>::iterator iter;
	    iter = edgeVote.begin();
	    while(iter != edgeVote.end()) {
	        cout << iter->first << " : " << iter->second << endl;
	        iter++;
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
	map<position,int> usedPos;
	priority_queue<candidate> candidates;
	priority_queue<candidate> substitutes;

	boundaryInfo boundary = {0,0,0,0};

	long seed = time(NULL);
	srand(seed+ *(int*)id * 12345);
	int randomStart = rand()%size;
	seed = rand();
	//cout<<"randomStart: "<<randomStart<<endl;

	position p = {0,0};

	putTile(tempAnswer,p,randomStart,usedTiles,usedPos,candidates,boundary);

	while(tempAnswer.size()!=size && answer[0]==-1)
	{
		while(!candidates.empty())
		{
			candidate winner = candidates.top();
			candidates.pop();
			//cout<<"winner id: "<<winner.id<<" x:"<<winner.x <<" y:"<<winner.y<<" relativeId:"<<winner.relativeId<<" orient:"<<winner.relativeOrient<<endl;
			position winnerPos;
			winnerPos.x = winner.x;
			winnerPos.y = winner.y;
			int winnerId = winner.id;

			if((max(boundary.maxCol,winnerPos.x)-min(boundary.minCol,winnerPos.x) >= width) || (max(boundary.maxRow,winnerPos.y)-min(boundary.minRow,winnerPos.y)>= width))
			{
				//cout<<"reason0"<<endl;
				continue;
			}

			if(usedPos.find(winnerPos) != usedPos.end())
			{
				//cout<<"reason1"<<endl;
				continue;
			}

			// int *pscore;
			// pscore = calculateValid(winnerId, winnerPos, tempAnswer, usedPos);
			bool *valid;
			string t="-1",r="-1",b="-1",l="-1";
			valid = validAndEdge(winnerId,winnerPos,tempAnswer,usedPos,t,r,b,l);
			if(!*valid)
			{
				//cout<<"reason2"<<endl;
				continue;
			}

			if(usedTiles.count(winnerId) == 1)
			{
				//cout<<"reason3"<<endl;
				substitutes.push(winner);
				continue;
			}

			putTile(tempAnswer,winnerPos,winnerId,usedTiles,usedPos,candidates,boundary);
			//cout<<"put "<<winnerId<<" at "<<"( "<<winnerPos.x<<","<<winnerPos.y<<" )"<<endl;
			// int relativeId = winner.relativeId;
			// int orient = winner.relativeOrient;
			// string key = "";
			// switch(orient)
			// {
			// 	case 0: key = to_string(relativeId)+"T-B"+to_string(winnerId);break;
			// 	case 1: key = to_string(winnerId)+"L-R"+to_string(relativeId);break;
			// 	case 2: key = to_string(winnerId)+"T-B"+to_string(relativeId);break;
			// 	case 3: key = to_string(relativeId)+"L-R"+to_string(winnerId);break;
			// };
			usedPos.insert(make_pair(winnerPos,winnerId));
			usedTiles.insert(winnerId);
			pthread_mutex_lock(&mutex1);
			//edgeVote[key]+=1;
			if(t!="-1")
				edgeVote[t]+=1;
			if(r!="-1")
				edgeVote[r]+=1;
			if(b!="-1")
				edgeVote[b]+=1;
			if(l!="-1")
				edgeVote[l]+=1;
			pthread_mutex_unlock(&mutex1);
		}

		//失败的大挪移
		// while(candidates.empty() && !substitutes.empty())
		// {
		// 	candidate substitute = substitutes.top();
		// 	substitutes.pop();
		// 	position targetPos = {substitute.x,substitute.y};
		// 	if(usedPos.find(targetPos)!=usedPos.end())
		// 		continue;
		// 	int subId = substitute.id;
		// 	position subOriPos = tempAnswer[subId];
		// 	usedTiles.erase(subId);
		// 	usedPos.erase(subOriPos);
		// 	tempAnswer.erase(subId);
		// 	updateBoundary(tempAnswer,boundary);
		// 	candidates.push(substitute);
		// 	break;
		// }
		// if(substitutes.empty())
		// {
		// 	break;
		// }
		//break;




		//随机移开n块
		/*
		if(candidates.empty() && tempAnswer.size()!=size)
		{
			//cout<<"random remove!"<<endl;
			//cout<<"before remove size:"<<tempAnswer.size();
			// map<position,int>::iterator iter2;
			// iter2 = usedPos.begin();
			// while(iter2!=usedPos.end())
			// {
			// 	cout<<"pos: "<<iter2->first.x<<","<<iter2->first.y<<" id: "<<iter2->second<<endl;
			// 	iter2++;
			// }
			srand(seed);
			int randomRemoveSize = rand()%(tempAnswer.size()-1)+1;
			//option：限制随机移动块数的上限
			// int theSize = min(5,int(tempAnswer.size()-1));
			// int randomRemoveSize = rand()%theSize+1;
			seed = rand();
			//cout<<"randomRemoveSize: "<<randomRemoveSize<<endl;
			//option：全部移除
			//randomRemoveSize = tempAnswer.size()-1;
			while(randomRemoveSize)
			{
				int randomId = -1;
				while(1)
				{
					srand(seed);
					randomId = rand()%size;
					seed = rand();
					if(usedTiles.count(randomId)==1)
						break;
				}
				//cout<<"remove randomId:"<<randomId<<endl;
				if(randomId != -1)
				{
					position oriPos = tempAnswer[randomId];
					usedTiles.erase(randomId);
					usedPos.erase(oriPos);
					tempAnswer.erase(randomId);
				}
				randomRemoveSize--;
			}
			map<int,position>::iterator it;
			it = tempAnswer.begin();
			map<int,position> forDeleteTempAnswer = tempAnswer;
			while(it != tempAnswer.end() && usedTiles.size()>1)
			{
				position p = it->second;
				int tileId = it->first;
				position t={p.x-1,p.y}, r={p.x,p.y+1}, b={p.x+1,p.y}, l={p.x,p.y-1};
				if(usedPos.find(t)==usedPos.end() && usedPos.find(r)==usedPos.end() &&usedPos.find(b)==usedPos.end() &&usedPos.find(l)==usedPos.end())
				{
					//cout<<"remove lonely:"<<tileId<<endl;
					position oriPos = tempAnswer[tileId];
					usedTiles.erase(tileId);
					usedPos.erase(oriPos);
					forDeleteTempAnswer.erase(tileId);
				}
				it++;
			}
			tempAnswer = forDeleteTempAnswer;

			//cout<<"after remove size:"<<tempAnswer.size()<<endl;

			updateBoundary(tempAnswer,boundary);
			
			it = tempAnswer.begin();
			while(it != tempAnswer.end())
			{
				position p = it->second;
				int tileId = it->first;
				updateCandidates(tempAnswer, p, tileId, usedTiles, usedPos, candidates, boundary);
				it++;
			}
		}
		*/
		





		//死局次数
		
		if(candidates.empty() && tempAnswer.size()!=size)
		{
			//cout<<"before remove size:"<<tempAnswer.size();
			// map<position,int>::iterator iter2;
			// iter2 = usedPos.begin();
			// while(iter2!=usedPos.end())
			// {
			// 	cout<<"pos: "<<iter2->first.x<<","<<iter2->first.y<<" id: "<<iter2->second<<endl;
			// 	iter2++;
			// }


			map<position,int>::iterator iter;
			iter = usedPos.begin();
			vector<string> edgeArray;
			while(iter!=usedPos.end())
			{
				position tempPos = iter->first;
				//right
				position right = {tempPos.x,tempPos.y+1};
				if(usedPos.find(right)!=usedPos.end())
				{
					int rightId = usedPos[right];
					string rightDeadEdge = to_string(iter->second)+"L-R"+to_string(rightId);
					edgeArray.push_back(rightDeadEdge);
				}
				//bottom
				position bottom = {tempPos.x+1,tempPos.y};
				if(usedPos.find(bottom)!=usedPos.end())
				{
					int bottomId = usedPos[bottom];
					string bottomDeadEdge = to_string(iter->second)+"T-B"+to_string(bottomId);
					edgeArray.push_back(bottomDeadEdge);
				}
				iter++;
			}
			//cout<<"edgeArray size:"<<edgeArray.size()<<endl;

			priority_queue<deadEdge> deadEdges;
			

			pthread_mutex_lock(&mutex1);
			for(int i=0;i<edgeArray.size();i++)
			{
				string e = edgeArray[i];
				edgeDead[e]+=1;
				deadEdge de;
				de.edge = e;
				de.score = edgeDead[e];
				deadEdges.push(de);
			}
			pthread_mutex_unlock(&mutex1);

			//cout<<"deadEdgescount: "<<deadEdges.size()<<endl;
			
			int randomRemoveESize;
			if(deadEdges.size()==1)
				randomRemoveESize = 1;
			else
			{
				srand(seed);
				randomRemoveESize = rand()%(deadEdges.size()-1)+1;
				//randomRemoveESize = rand()%(deadEdges.size()*2/3)+1;
				seed = rand();
			}
			
			//cout<<"randomRemoveESize: "<<randomRemoveESize<<endl;
			while(randomRemoveESize && tempAnswer.size()>1)
			{
				deadEdge d_ = deadEdges.top();
				string removeEdge = d_.edge;
				//cout<<"the edge:"<<removeEdge<<endl;
				deadEdges.pop();
				string::size_type findResultPos;
				findResultPos = removeEdge.find("L-R");
				if(findResultPos != removeEdge.npos)
				{
					int tempIds[2];
					tempIds[0] = stoi(removeEdge.substr(0,findResultPos));
					tempIds[1] = stoi(removeEdge.substr(findResultPos+3));
					if(usedTiles.count(tempIds[0])==1 && usedTiles.count(tempIds[1])==1)
					{
						srand(seed);
						int theId = rand()%2;
						seed = rand();
						theId = tempIds[theId];
						//cout<<"remove:"<<theId<<endl;
						
						position oriPos = tempAnswer[theId];
						usedTiles.erase(theId);
						usedPos.erase(oriPos);
						tempAnswer.erase(theId);
					}
					

				}
				else
				{
					findResultPos = removeEdge.find("T-B");
					int tempIds[2];
					tempIds[0] = stoi(removeEdge.substr(0,findResultPos));
					tempIds[1] = stoi(removeEdge.substr(findResultPos+3));
					if(usedTiles.count(tempIds[0])==1 && usedTiles.count(tempIds[1])==1)
					{
						srand(seed);
						int theId = rand()%2;
						seed = rand();
						theId = tempIds[theId];
						//cout<<"remove:"<<theId<<endl;

						position oriPos = tempAnswer[theId];
						usedTiles.erase(theId);
						usedPos.erase(oriPos);
						tempAnswer.erase(theId);

					}
					
				}

				randomRemoveESize--;

			}
			//cout<<"size before jiancha:"<<tempAnswer.size()<<endl;
			map<int,position>::iterator it;
			it = tempAnswer.begin();
			map<int,position> forDeleteTempAnswer = tempAnswer;
			while(it != tempAnswer.end() && usedTiles.size()>1)
			{
				position p = it->second;
				int tileId = it->first;
				//cout<<"daijianchadeid:"<<tileId<<endl;
				position t={p.x-1,p.y}, r={p.x,p.y+1}, b={p.x+1,p.y}, l={p.x,p.y-1};
				if(usedPos.find(t)==usedPos.end() && usedPos.find(r)==usedPos.end() &&usedPos.find(b)==usedPos.end() &&usedPos.find(l)==usedPos.end())
				{
					//cout<<"remove lonely:"<<tileId<<endl;
					position oriPos = tempAnswer[tileId];
					usedTiles.erase(tileId);
					usedPos.erase(oriPos);
					forDeleteTempAnswer.erase(tileId);
				}
				it++;
			}
			tempAnswer = forDeleteTempAnswer;
			//cout<<"after remove size:"<<tempAnswer.size()<<endl;
			updateBoundary(tempAnswer,boundary);
			
			it = tempAnswer.begin();
			while(it != tempAnswer.end())
			{
				position p = it->second;
				int tileId = it->first;
				updateCandidates(tempAnswer, p, tileId, usedTiles, usedPos, candidates, boundary);
				it++;
			}
		}


		
		
	}
	


	if(tempAnswer.size() == size && answer[0]==-1)
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
		// cout<<"thread "<<*(int*)id<<" failed,the result size:"<<tempAnswer.size()<<endl;
		// map<int,position>::iterator iter;
		// iter = tempAnswer.begin();
		// while(iter!=tempAnswer.end())
		// {
		// 	cout<<"id:"<<iter->first<<" pos:"<<iter->second.x<<","<<iter->second.y<<endl;
		// 	iter++;
		// }
		// map<position,int>::iterator iter2;
		// iter2 = usedPos.begin();
		// while(iter2!=usedPos.end())
		// {
		// 	cout<<"pos: "<<iter2->first.x<<","<<iter2->first.y<<" id: "<<iter2->second<<endl;
		// 	iter2++;
		// }

	}
	
	return NULL;
}

void* putTile(map<int,position> &tempAnswer, position p, int tileId, set<int> &usedTiles, map<position,int> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary)
{
	tempAnswer[tileId]  = p;
	usedTiles.insert(tileId);
	usedPos.insert(make_pair(p,tileId));
	updateBoundary(tempAnswer,boundary);
	updateCandidates(tempAnswer, p, tileId, usedTiles, usedPos, candidates, boundary);
	return NULL;

}

void* updateCandidates(map<int,position> &tempAnswer,  position p, int tileId ,set<int> &usedTiles, map<position,int> &usedPos, priority_queue<candidate> &candidates, boundaryInfo &boundary)
{
	if(tempAnswer.size()!=size)
	{
		int x = p.x;
		int y = p.y;
		//top
		position top = {x-1,y};
		if(usedPos.find(top)== usedPos.end())
		{
			if((max(boundary.maxCol,top.x)-min(boundary.minCol,top.x) < width) && (max(boundary.maxRow,top.y)-min(boundary.minRow,top.y) < width) )
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
			        	int candidateId = stoi(iter->first.substr(0,findResultPos));
			        	int *pscore;
			        	pscore = calculateValid(candidateId, top, tempAnswer, usedPos);
			        	if(*pscore)
			        	{
			        		candidate can;
			        		can.id = candidateId;
			        		can.x = top.x;
			        		can.y = top.y;
			        		can.score = *pscore;
			        		//can.score = iter->second;
			        		can.relativeId = tileId;
			        		can.relativeOrient = 2; //bottom
			        		candidates.push(can);
			        		// boundary.minCol = min(top.x,boundary.minCol);
			        		// boundary.minRow = min(top.y,boundary.minRow);
			        		// boundary.maxCol = max(top.x,boundary.maxCol);
			        		// boundary.maxRow = max(top.y,boundary.maxRow);
			        	}
			        	
			        	
			        }
			        iter++;
			    }
				
			}
		}
		//right
		position right = {x,y+1};
		if(usedPos.find(right)==usedPos.end())
		{
			if((max(boundary.maxCol,right.x)-min(boundary.minCol,right.x) < width) && (max(boundary.maxRow,right.y)-min(boundary.minRow,right.y) < width) )
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
			        	int *pscore;
			        	pscore = calculateValid(candidateId, right, tempAnswer, usedPos);
			        	if(*pscore)
			        	{
			        		can.x = right.x;
			        		can.y = right.y;
			        		can.score = *pscore;
			        		//can.score = iter->second;
			        		can.relativeId = tileId;
			        		can.relativeOrient = 3; //left
			        		candidates.push(can);
			        		// boundary.minCol = min(right.x,boundary.minCol);
			        		// boundary.minRow = min(right.y,boundary.minRow);
			        		// boundary.maxCol = max(right.x,boundary.maxCol);
			        		// boundary.maxRow = max(right.y,boundary.maxRow);
			        	}

			        }
			        iter++;
			    }
				
			}
		}
		//bottom
		position bottom = {x+1,y};
		if(usedPos.find(bottom)==usedPos.end())
		{
			if((max(boundary.maxCol,bottom.x)-min(boundary.minCol,bottom.x) < width) && (max(boundary.maxRow,bottom.y)-min(boundary.minRow,bottom.y) < width) )
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
			        	int *pscore;
			        	pscore = calculateValid(candidateId, bottom, tempAnswer, usedPos);
			        	if(*pscore)
			        	{
			        		candidate can;
			        		can.id = candidateId;
			        		can.x = bottom.x;
			        		can.y = bottom.y;
			        		can.score = *pscore;
			        		//can.score = iter->second;
			        		can.relativeId = tileId;
			        		can.relativeOrient = 0; //top
			        		candidates.push(can);
			        		// boundary.minCol = min(bottom.x,boundary.minCol);
			        		// boundary.minRow = min(bottom.y,boundary.minRow);
			        		// boundary.maxCol = max(bottom.x,boundary.maxCol);
			        		// boundary.maxRow = max(bottom.y,boundary.maxRow);
			        	}
			        	
			        }
			         iter++;
			    }	
			}
		}
		//left
		position left = {x,y-1};
		if(usedPos.find(left)==usedPos.end())
		{
			if((max(boundary.maxCol,left.x)-min(boundary.minCol,left.x) < width) && (max(boundary.maxRow,left.y)-min(boundary.minRow,left.y) < width) )
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
			        	int *pscore;
			        	pscore = calculateValid(candidateId, left, tempAnswer, usedPos);
			        	if(*pscore)
			        	{
			        		candidate can;
			        		can.id = candidateId;
			        		can.x = left.x;
			        		can.y = left.y;
			        		can.score = *pscore;
			        		//can.score = iter->second;
			        		can.relativeId = tileId;
			        		can.relativeOrient = 1; //right
			        		candidates.push(can);
			        		// boundary.minCol = min(left.x,boundary.minCol);
			        		// boundary.minRow = min(left.y,boundary.minRow);
			        		// boundary.maxCol = max(left.x,boundary.maxCol);
			        		// boundary.maxRow = max(left.y,boundary.maxRow);
			        	}
			        	
			        }
			        iter++;
			    }	
			}
		}
	}
	return NULL;
}

int* calculateValid(int candidateId, position pos, map<int,position> &tempAnswer, map<position,int> &usedPos)
{
	int score = 0;
	position top = {pos.x-1,pos.y};
	position right = {pos.x,pos.y+1};
	position bottom = {pos.x+1,pos.y};
	position left = {pos.x,pos.y-1};
	int topId = -1;
	int rightId = -1;
	int bottomId = -1;
	int leftId = -1;
	bool valid = true;
	if(usedPos.find(top)!=usedPos.end())
	{
		topId = usedPos[top];
		if(tiles[topId].bottom != tiles[candidateId].top)
		{
			valid = false;
		}
	}
	if(usedPos.find(right)!=usedPos.end())
	{
		rightId = usedPos[right];
		if(tiles[rightId].left != tiles[candidateId].right)
		{
			valid = false;
		}
	}
	if(usedPos.find(bottom)!=usedPos.end())
	{
		bottomId = usedPos[bottom];
		if(tiles[bottomId].top != tiles[candidateId].bottom)
		{
			valid = false;
		}
	}
	if(usedPos.find(left)!=usedPos.end())
	{
		leftId = usedPos[left];
		if(tiles[leftId].right != tiles[candidateId].left)
		{
			valid = false;
		}
	}
	if(!valid)
	{
		int *p = &score;
		return p;
	}
	string s;
	if(topId != -1)
	{
		s = to_string(topId)+"T-B"+to_string(candidateId);
		score+=edgeVote[s];
	}
	if(rightId != -1)
	{
		s = to_string(candidateId)+"L-R"+to_string(rightId);
		score+=edgeVote[s];
	}
	if(bottomId != -1)
	{
		s = to_string(candidateId)+"T-B"+to_string(bottomId);
		score+=edgeVote[s];
	}
	if(leftId != -1)
	{
		s = to_string(leftId)+"L-R"+to_string(candidateId);
		score+=edgeVote[s];
	}
	int *p = &score;
	return p;
}

void* updateBoundary(map<int,position> &tempAnswer, boundaryInfo &boundary)
{
	map<int,position>::iterator iter;
	iter = tempAnswer.begin();
	int minCol=iter->second.x;
	int minRow=iter->second.y;
	int maxCol=iter->second.x;
	int maxRow=iter->second.y;
	while(iter!=tempAnswer.end())
	{
		int x = iter->second.x;
		int y = iter->second.y;
		minCol = min(x,minCol);
		minRow = min(y,minRow);
		maxCol = max(x,maxCol);
		maxRow = max(y,maxRow);
		iter++;
	}
	boundary.minCol = minCol;
	boundary.minRow = minRow;
	boundary.maxCol = maxCol;
	boundary.maxRow = maxRow;
	return NULL;
}

bool* validAndEdge(int candidateId, position pos, map<int,position> &tempAnswer, map<position,int> &usedPos, string &t, string &r, string &b, string &l)
{
	position top = {pos.x-1,pos.y};
	position right = {pos.x,pos.y+1};
	position bottom = {pos.x+1,pos.y};
	position left = {pos.x,pos.y-1};
	int topId = -1;
	int rightId = -1;
	int bottomId = -1;
	int leftId = -1;
	bool valid = true;
	if(usedPos.find(top)!=usedPos.end())
	{
		topId = usedPos[top];
		if(tiles[topId].bottom != tiles[candidateId].top)
		{
			valid = false;
		}
	}
	if(usedPos.find(right)!=usedPos.end())
	{
		rightId = usedPos[right];
		if(tiles[rightId].left != tiles[candidateId].right)
		{
			valid = false;
		}
	}
	if(usedPos.find(bottom)!=usedPos.end())
	{
		bottomId = usedPos[bottom];
		if(tiles[bottomId].top != tiles[candidateId].bottom)
		{
			valid = false;
		}
	}
	if(usedPos.find(left)!=usedPos.end())
	{
		leftId = usedPos[left];
		if(tiles[leftId].right != tiles[candidateId].left)
		{
			valid = false;
		}
	}
	if(!valid)
	{
		bool *p = &valid;
		return p;
	}
	if(topId != -1)
	{
		t = to_string(topId)+"T-B"+to_string(candidateId);
	}
	if(rightId != -1)
	{
		r = to_string(candidateId)+"L-R"+to_string(rightId);
	}
	if(bottomId != -1)
	{
		b = to_string(candidateId)+"T-B"+to_string(bottomId);
	}
	if(leftId != -1)
	{
		l = to_string(leftId)+"L-R"+to_string(candidateId);
	}
	bool *p = &valid;
	return p;
}

