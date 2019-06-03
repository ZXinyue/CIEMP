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
#define RETURN_NUM 5


typedef struct
{
	int id;
	string top;
	string right;
	string bottom;
	string left;
	
}Tile;

struct Pos
{
	friend bool operator < (Pos n1, Pos n2)
	{
	    return n1.priority < n2.priority;
	}
	int priority;
	int x;
	int y;

};



int size = 36;
int width = 6;
Tile tiles[36];
map<string,int> edgeVote;
int* answer;

void* searchAnswer(void *id);
bool* BFS(int *tempAnswer, bool *flag, int *returnTime, bool *used, int *pPos, long *pSeed);
void* backTrack(int *px, int *py, int *tempAnswer,int *returnTime, bool *used);
void* justTest(void *id);
void* findCandidates(int *px, int *py, int *tempAnswer, map<int,double> &candidates, int &leftTileId, int &rightTileId, int &topTileId, int &bottomTileId);
void* backTrack2(int *px, int *py, int *tempAnswer,bool *used, long *pSeed);


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

	cout<<"the answer is:"<<endl;
	for(int j=0; j<size; j++)
	{
		cout<<answer[j]<<" ";
	}

	end = clock();

	cout<<"the total time is:"<<((double)(end-start))/CLOCKS_PER_SEC << endl;
	
	return 0;
}

void* justTest(void *id)
{
	cout<<"test thread id:"<<*(int*)id<<endl;
}

void* searchAnswer(void *id)
{
	//cout<<"hi searchAnswer!~id:"<<*(int*)id<<endl;
	bool *flag;
	flag = (bool *)malloc(sizeof(bool)*size);
	for(int i=0; i<size; i++)
		flag[i]=false;

	int *tempAnswer;
	tempAnswer = (int *)malloc(sizeof(int)*size);
	for(int i=0; i<size; i++)
		tempAnswer[i]=-1;

	int *returnTime;
	returnTime = (int *)malloc(sizeof(int)*size);
	for(int i=0; i<size; i++)
		returnTime[i]=0;

	bool *used;
	used = (bool *)malloc(sizeof(bool)*size);
	for(int i=0; i<size; i++)
		used[i]=false;

	long seed = time(NULL);
	srand(seed+ *(int*)id * 12345);
	int randomStart = rand()%size;
	int randomPosition = rand()%size;
	seed = rand();
	//cout<<"the thread"<<*(int*)id<<"randomStart is:  "<<randomStart<<endl;
	tempAnswer[randomPosition] = randomStart;
	used[randomStart] = true;
	
	

	while(1)
	{
		if(answer[0]!=-1)
			return NULL;

		// if(tempAnswer[0] == -1)
		// {
		// 	Seed = time(NULL);
		// 	srand(Seed+ *(int*)id * 12345);
		// 	randomStart = rand()%size;
		// 	tempAnswer[0] = randomStart;
		// 	used[randomStart] = true;
		// 	//cout<<"new randomStart "<< *(int*)id <<" is!--"<<randomStart<<endl;
		// }

		int BFSStart = -1;
		vector<int> validArray;
		for(int i=0; i<size; i++)
		{
			if(tempAnswer[i]!=-1)
			{
				validArray.push_back(i);
			}
		}


		if(validArray.empty())
		{
			srand(seed);
			BFSStart = rand()%size;
			randomStart = rand()%size;
			tempAnswer[BFSStart] = randomStart;
			used[randomStart] = true;
			seed = rand();
		}
		else
		{
			srand(seed);
			int randomIndex = rand()%validArray.size();
			BFSStart = validArray[randomIndex];
			seed = rand();
		}
		int *pStart = &BFSStart;
		long *pSeed = &seed;


		bool *returnStatus = BFS(tempAnswer, flag, returnTime, used, pStart, pSeed);
		// cout<<"returnStatus is:"<<*returnStatus<<endl;
		if(*returnStatus)
			break;
	}

	
	pthread_mutex_lock(&mutex1);
	answer = tempAnswer;
	pthread_mutex_unlock(&mutex1);

	return NULL;

}

void* backTrack(int *px, int *py, int *tempAnswer,int *returnTime, bool *used)
{
	//cout<<"hi backTrack"<<endl;
	int x = *px;
	int y = *py;


	if(x-1>=0 && tempAnswer[(x-1)*width+y]!=-1) //top
	{
		
		used[tempAnswer[(x-1)*width+y]] = false;
		tempAnswer[(x-1)*width+y] = -1;

		int x_ = x-1;
		int *px_ = &x_;
		int *py_ = py;
		// int leftTileId = -1;
		// int rightTileId = -1;
		// int topTileId = -1;
		// int bottomTileId = -1;
		// map<int,double> topCandidates;
		// findCandidates(px,py,tempAnswer,topCandidates,leftTileId,rightTileId,topTileId,bottomTileId);
		

		returnTime[(x-1)*width+y] += 1;
		
		if(returnTime[(x-1)*width+y]>=RETURN_NUM)
		{
			//cout<<"aaaa"<<endl;
			
			returnTime[(x-1)*width+y] = 0;
			backTrack(px_, py_, tempAnswer,returnTime, used);
		}
	}
	if(y-1>=0 && tempAnswer[x*width+y-1]!=-1) //left
	{
		used[tempAnswer[x*width+y-1]] = false;
		returnTime[x*width+y-1] += 1;
		tempAnswer[x*width+y-1] = -1;
		if(returnTime[x*width+y-1]>=RETURN_NUM)
		{
			//cout<<"aaaa"<<endl;
			int y_ = y-1;
			int *py_ = &y_;
			int *px_ = px;
			returnTime[x*width+y-1] = 0;
			backTrack(px_, py_, tempAnswer,returnTime, used);
		}
	}
	if(x+1<width && tempAnswer[(x+1)*width+y]!=-1) //bottom
	{
		used[tempAnswer[(x+1)*width+y]] = false;
		returnTime[(x+1)*width+y] += 1;
		tempAnswer[(x+1)*width+y] = -1;
		if(returnTime[(x+1)*width+y]>=RETURN_NUM)
		{
			//cout<<"aaaa"<<endl;
			int x_ = x+1;
			int *px_ = &x_;
			int *py_ = py;
			returnTime[(x+1)*width+y] = 0;
			backTrack(px_, py_, tempAnswer,returnTime, used);
		}
	}
	if(y+1<width && tempAnswer[x*width+y+1]!=-1) //right
	{
		used[tempAnswer[x*width+y+1]] = false;
		returnTime[x*width+y+1] += 1;
		tempAnswer[x*width+y+1] = -1;
		if(returnTime[x*width+y+1]>=RETURN_NUM)
		{
			//cout<<"aaaa"<<endl;
			int y_ = y+1;
			int *py_ = &y_;
			int *px_ = px;
			returnTime[x*width+y+1] = 0;
			backTrack(px_, py_, tempAnswer,returnTime, used);
		}

	}
	return NULL;	
}

void* backTrack2(int *px, int *py, int *tempAnswer,bool *used, long *pSeed)
{
	map<int,double> topCandidates;
	map<int,double> rightCandidates;
	map<int,double> bottomCandidates;
	map<int,double> leftCandidates;
	int x = *px;
	int y = *py;
	long seed = *pSeed;
	int direIds[4];
	for(int i=0;i<4;i++)
		direIds[i] = -1;

	int x_ = x-1;
	int *px_ = &x_;
	int y_ = y-1;
	int *py_ = &y_;
	int _x = x+1;
	int *p_x = &_x;
	int _y = y+1;
	int *p_y = &_y;

	int aroundId[16];
	for(int i=0; i<16; i++)
		aroundId[i]=-1;
	//self
	if(tempAnswer[x*width+y] != -1)
	{
		int theId = tempAnswer[x*width+y];
		used[theId] = false;
		tempAnswer[x*width+y] = -1;
	}
	//top
	if(x_>=0 && tempAnswer[x_*width+y]!=-1)
	{
		direIds[0] = tempAnswer[x_*width+y];
		findCandidates(px_,py,tempAnswer,topCandidates,aroundId[3],aroundId[1],aroundId[0],aroundId[2]);
		//cout<<"topCandidatesSize!!"<<topCandidates.size()<<endl;
	}
	//right
	if(_y<width && tempAnswer[x*width+_y] != -1)
	{
		direIds[1] = tempAnswer[x*width+_y];
		findCandidates(px,p_y,tempAnswer,rightCandidates,aroundId[7],aroundId[5],aroundId[4],aroundId[6]);
	}
	//bottom
	if(_x<width && tempAnswer[_x*width+y] != -1)
	{
		direIds[2] = tempAnswer[_x*width+y];
		findCandidates(p_x,py,tempAnswer,bottomCandidates,aroundId[11],aroundId[9],aroundId[8],aroundId[10]);
	}
	//left
	if(y_>=0 && tempAnswer[x*width+y_] != -1)
	{
		direIds[3] = tempAnswer[x*width+y_];
		findCandidates(px,py_,tempAnswer,leftCandidates,aroundId[15],aroundId[13],aroundId[12],aroundId[14]);
	}
	vector<int> direCan;
	set<int> directions;
	for(int i=0;i<4;i++)
	{
		if(direIds[i] != -1)
		{
			direCan.push_back(i);
		}
	}
	if(direCan.size() == 0)
	{
		return NULL;
	}
	// cout<<"x:"<<x_<<" y:"<<y<<" upTile:"<<direIds[0]<<endl;
	// if(x_-1>=0)
	// 	cout<<"up here:"<<tempAnswer[(x_-1)*width+y]<<endl;
	// if(y-1>=0)
	// 	cout<<"left here:"<<tempAnswer[x_*width+y-1]<<endl;
	// if(y+1<width)
	// 	cout<<"right here:"<<tempAnswer[x_*width+y+1]<<endl;
	// cout<<"------------"<<endl;
	// cout<<aroundId[3]<<" "<<aroundId[1]<<" "<<aroundId[0]<<" "<<aroundId[2]<<endl;
	// cout<<"------------"<<endl;

	// cout<<"cansize"<<topCandidates.size()<<endl;


	srand(seed);
	int backTrackDireNum = rand()%direCan.size()+1;
	seed = rand();

	while(directions.size()!= backTrackDireNum)
	{
		srand(seed);
		int tempIndex = rand()%direCan.size();
		directions.insert(direCan[tempIndex]);
		seed = rand();
	}

	//backtrack top
	if(directions.count(0) == 1)
	{
		used[direIds[0]] = false;
		tempAnswer[x_*width+y] = -1;
		if(topCandidates.size() <= 1)
		{
			// pSeed = &seed;
			// backTrack2(px_,py,tempAnswer,used,pSeed);
		}
		else
		{
			map<int,double>::iterator it;
			it = topCandidates.begin();
			int minscore = -1;
			int winnerId = -1;
			while(it!=topCandidates.end())
			{
				if(it->first == direIds[0])
				{
					it++;
					continue;
				}
				if(minscore == -1)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				else if(it->second<minscore)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				it++;
			}
			//cout<<"x:"<<x_<<" y:"<<y<<"..........................winnerId:"<<winnerId<<endl;
			if(!used[winnerId])
			{
				used[winnerId] = true;
				tempAnswer[x_*width+y] = winnerId;
			}
			else
			{
				for(int i=0;i<size;i++)
				{
					if(tempAnswer[i]==winnerId)
					{
						tempAnswer[i]=-1;
						used[winnerId]=false;
						break;
					}
				}
				used[winnerId] = true;
				tempAnswer[x_*width+y] = winnerId;
			}
			pthread_mutex_lock(&mutex1);
			if(aroundId[3] != -1)
			{
				string ls = to_string(aroundId[3])+"L-R"+to_string(winnerId);
				edgeVote[ls] += 1;
			}
			if(aroundId[1] != -1)
			{
				string rs = to_string(winnerId)+"L-R"+to_string(aroundId[1]);
				edgeVote[rs] += 1;
			}
			if(aroundId[0] != -1)
			{
				string ts = to_string(aroundId[0])+"T-B"+to_string(winnerId);
				edgeVote[ts] += 1;
			}
			if(aroundId[2] != -1)
			{
				string bs = to_string(winnerId)+"T-B"+to_string(aroundId[2]);
				edgeVote[bs] += 1;
			}
			pthread_mutex_unlock(&mutex1);
		}
	}

	//backtrap right
	if(directions.count(1) == 1)
	{
		used[direIds[1]] = false;
		tempAnswer[x*width+_y] = -1;
		if(rightCandidates.size() <= 1)
		{
			// pSeed = &seed;
			// backTrack2(px,p_y,tempAnswer,used,pSeed);
		}
		else
		{
			map<int,double>::iterator it;
			it = rightCandidates.begin();
			int minscore = -1;
			int winnerId = -1;
			while(it!=rightCandidates.end())
			{
				if(it->first == direIds[1])
				{
					it++;
					continue;
				}
				if(minscore == -1)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				else if(it->second<minscore)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				it++;
			}
			if(!used[winnerId])
			{
				used[winnerId] = true;
				tempAnswer[x*width+_y] = winnerId;
			}
			else
			{
				for(int i=0;i<size;i++)
				{
					if(tempAnswer[i]==winnerId)
					{
						tempAnswer[i]=-1;
						used[winnerId]=false;
						break;
					}
				}
				used[winnerId] = true;
				tempAnswer[x*width+_y] = winnerId;
			}
			pthread_mutex_lock(&mutex1);
			if(aroundId[7] != -1)
			{
				string ls = to_string(aroundId[7])+"L-R"+to_string(winnerId);
				edgeVote[ls] += 1;
			}
			if(aroundId[5] != -1)
			{
				string rs = to_string(winnerId)+"L-R"+to_string(aroundId[5]);
				edgeVote[rs] += 1;
			}
			if(aroundId[4] != -1)
			{
				string ts = to_string(aroundId[4])+"T-B"+to_string(winnerId);
				edgeVote[ts] += 1;
			}
			if(aroundId[6] != -1)
			{
				string bs = to_string(winnerId)+"T-B"+to_string(aroundId[6]);
				edgeVote[bs] += 1;
			}
			pthread_mutex_unlock(&mutex1);
		}
	}
	
	//backtrap bottom
	if(directions.count(2) == 1)
	{
		used[direIds[2]] = false;
		tempAnswer[_x*width+y] = -1;
		if(bottomCandidates.size() <= 1)
		{
			// pSeed = &seed;
			// backTrack2(p_x,py,tempAnswer,used,pSeed);
		}
		else
		{
			map<int,double>::iterator it;
			it = bottomCandidates.begin();
			int minscore = -1;
			int winnerId = -1;
			while(it!=bottomCandidates.end())
			{
				if(it->first == direIds[2])
				{
					it++;
					continue;
				}
				if(minscore == -1)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				else if(it->second<minscore)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				it++;
			}
			if(!used[winnerId])
			{
				used[winnerId] = true;
				tempAnswer[_x*width+y] = winnerId;
			}
			else
			{
				for(int i=0;i<size;i++)
				{
					if(tempAnswer[i]==winnerId)
					{
						tempAnswer[i]=-1;
						used[winnerId]=false;
						break;
					}
				}
				used[winnerId] = true;
				tempAnswer[_x*width+y] = winnerId;
			}
			pthread_mutex_lock(&mutex1);
			if(aroundId[11] != -1)
			{
				string ls = to_string(aroundId[11])+"L-R"+to_string(winnerId);
				edgeVote[ls] += 1;
			}
			if(aroundId[9] != -1)
			{
				string rs = to_string(winnerId)+"L-R"+to_string(aroundId[9]);
				edgeVote[rs] += 1;
			}
			if(aroundId[8] != -1)
			{
				string ts = to_string(aroundId[8])+"T-B"+to_string(winnerId);
				edgeVote[ts] += 1;
			}
			if(aroundId[10] != -1)
			{
				string bs = to_string(winnerId)+"T-B"+to_string(aroundId[10]);
				edgeVote[bs] += 1;
			}
			pthread_mutex_unlock(&mutex1);
		}
	}

	//backtrap left
	if(directions.count(3) == 1)
	{
		used[direIds[3]] = false;
		tempAnswer[x*width+y_] = -1;
		if(leftCandidates.size() <= 1)
		{
			// pSeed = &seed;
			// backTrack2(px,py_,tempAnswer,used,pSeed);
		}
		else
		{
			map<int,double>::iterator it;
			it = leftCandidates.begin();
			int minscore = -1;
			int winnerId = -1;
			while(it!=leftCandidates.end())
			{
				if(it->first == direIds[3])
				{
					it++;
					continue;
				}
				if(minscore == -1)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				else if(it->second<minscore)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				it++;
			}
			if(!used[winnerId])
			{
				used[winnerId] = true;
				tempAnswer[x*width+y_] = winnerId;
			}
			else
			{
				for(int i=0;i<size;i++)
				{
					if(tempAnswer[i]==winnerId)
					{
						tempAnswer[i]=-1;
						used[winnerId]=false;
						break;
					}
				}
				used[winnerId] = true;
				tempAnswer[x*width+y_] = winnerId;
			}
			pthread_mutex_lock(&mutex1);
			if(aroundId[15] != -1)
			{
				string ls = to_string(aroundId[15])+"L-R"+to_string(winnerId);
				edgeVote[ls] += 1;
			}
			if(aroundId[13] != -1)
			{
				string rs = to_string(winnerId)+"L-R"+to_string(aroundId[13]);
				edgeVote[rs] += 1;
			}
			if(aroundId[12] != -1)
			{
				string ts = to_string(aroundId[12])+"T-B"+to_string(winnerId);
				edgeVote[ts] += 1;
			}
			if(aroundId[14] != -1)
			{
				string bs = to_string(winnerId)+"T-B"+to_string(aroundId[14]);
				edgeVote[bs] += 1;
			}
			pthread_mutex_unlock(&mutex1);
		}
	}
	return NULL;
}


bool* BFS(int *tempAnswer, bool *flag, int *returnTime, bool *used, int *pPos, long *pSeed)
{
	//cout<<"hi BFS"<<endl;
	bool RS = true;
	int pos = *pPos;
	int startx = pos/width;
	int starty = pos%width;
	//priority_queue<Pos> q;
	queue<int> xq;
	queue<int> yq;
	xq.push(startx);
	yq.push(starty);

	while(!xq.empty()&&!yq.empty())
	{
		//cout<<"hi loop in queue"<<endl;
		int x = xq.front();
		int y = yq.front();
		xq.pop();
		yq.pop();

		//cout<<"x"<<x<<"y"<<y<<endl;

		if(flag[x*width+y])
			continue;

		flag[x*width+y] = true;


		if(y+1<width && !flag[x*width+y+1])
		{
			xq.push(x);
			yq.push(y+1);
		}
		if(x+1<width && !flag[(x+1)*width+y])
		{
			xq.push(x+1);
			yq.push(y);
		}
		if(y-1>=0 && !flag[x*width+y-1])
		{
			xq.push(x);
			yq.push(y-1);
		}
		if(x-1>=0 && !flag[(x-1)*width+y])
		{
			xq.push(x-1);
			yq.push(y);
		}
		if(tempAnswer[x*width+y]!=-1)
			continue;
		
		// else
		// {
		// 	srand(*pSeed);
		// 	*pSeed = rand();
		// 	//cout<<"same seed111??"<<seed<<endl;
		// 	// cout<<"same rand??"<<rand()<<endl;
		// 	int randomSelect = rand()%size;
		// 	while(used[randomSelect])
		// 		randomSelect = rand()%size;
		// 	tempAnswer[x*width+y]=randomSelect;
		// 	//cout<<"x:"<<x<<" y:"<<y<<" randomSelect is"<<randomSelect<<endl;
		// 	used[randomSelect] = true;
		// 	continue;
		// }

		map<int,double> candidates;
		int* px = &x;
		int* py = &y;
		int leftTileId = -1;
		int rightTileId = -1;
		int topTileId = -1;
		int bottomTileId = -1;
		findCandidates(px,py,tempAnswer,candidates,leftTileId,rightTileId,topTileId,bottomTileId);

	    if(!candidates.empty())
	    {
	    	double sum = 0;
	    	map<int,double>::iterator it;

	    	//概率选取
	    	/*
	    	it = candidates.begin();
	    	while(it != candidates.end())
	    	{
	    		it->second = 1/it->second;
	    		sum+=it->second;
	    		it++;
	    	}

	    	it = candidates.begin();
	    	while(it != candidates.end())
	    	{
	    		it->second = it->second/sum;
	    		it++;
	    	}

	    	srand(seed); 
	    	seed = rand();
	    	
	    	double rnd = rand() / double(RAND_MAX);

	    	int winnerId;

	    	it = candidates.begin(); 
	    	while(it != candidates.end())
	    	{
	    		rnd -= it->second;
	    		if(rnd < 0)
	    		{
	    			winnerId = it->first;
	    			break;
	    		}
	    		it++;
	    	}
			*/

			//选择最小
			
			it = candidates.begin();
			int winnerId = it->first;
			int minscore = it->second;
			while(it != candidates.end())
			{
				//cout<<"candidates"<<it->first<<" "<<it->second<<endl;
				if(it->second<minscore)
				{
					minscore = it->second;
					winnerId = it->first;
				}
				it++;
			}
			if(used[winnerId])
			{
				for(int i=0;i<size;i++)
				{
					if(tempAnswer[i]==winnerId)
					{
						tempAnswer[i]=-1;
						used[winnerId]=false;
						flag[i] = false;
						xq.push(i/width);
						yq.push(i%width);
						break;
					}
				}
			}
	    	
	    	tempAnswer[x*width+y]=winnerId;
	    	// for(int i=0;i<size;i++)
	    	// 	cout<<tempAnswer[i]<<" ";
	    	// cout<<endl;
	    	// cout<<"x:"<<x<<"y:"<<y<<endl;
	    	
	    	used[winnerId] = true;
	    	//cout<<"winnerId:"<<winnerId<<endl;
	    	pthread_mutex_lock(&mutex1);
	    	if(leftTileId != -1)
	    	{
	    		string ls = to_string(leftTileId)+"L-R"+to_string(winnerId);
	    		edgeVote[ls] += 1;
	    	}
	    	if(rightTileId != -1)
	    	{
	    		string rs = to_string(winnerId)+"L-R"+to_string(rightTileId);
	    		edgeVote[rs] += 1;
	    	}
	    	if(topTileId != -1)
	    	{
	    		string ts = to_string(topTileId)+"T-B"+to_string(winnerId);
	    		edgeVote[ts] += 1;
	    	}
	    	if(bottomTileId != -1)
	    	{
	    		string bs = to_string(winnerId)+"T-B"+to_string(bottomTileId);
	    		edgeVote[bs] += 1;
	    	}
	    	pthread_mutex_unlock(&mutex1);
	    	
	    }
	    else
	    {
	    	// cout<<"failed tempAnswer is:"<<endl;
	    	// for(int i=0; i<size; i++)
	    	// 	cout<<tempAnswer[i]<<" ";
	    	// cout<<"edgeVote this time:"<<endl;
	    	// map<string,int>::iterator iter;
	    	//     iter = edgeVote.begin();
	    	//     while(iter != edgeVote.end()) {
	    	//         cout << iter->first << " : " << iter->second << endl;
	    	//         iter++;
	    	//     }

	    	for(int i=0; i<size; i++)
	    		flag[i]=false;

	    	// long seed = time(NULL);
	    	// //srand(seed+(int)id*12345);
	    	// int randomStart = rand()%size;
	    	// tempAnswer[0] = randomStart;
	    	// used[randomStart] = true;

	    	//backTrack
	    	int *px = &x;
	    	int *py = &y;
	    	//backTrack(px,py,tempAnswer,returnTime,used);
	    	backTrack2(px,py,tempAnswer,used,pSeed);

	    	RS = false;
	    	bool* rs = &RS;
	    	return rs;

	    }
	}
	RS = true;
	bool* rs = &RS;
	return rs;
}

void* findCandidates(int *px, int *py, int *tempAnswer, map<int,double> &candidates, int &leftTileId, int &rightTileId, int &topTileId, int &bottomTileId)
{
	int x = *px;
	int y = *py;

	Tile leftTile;
	Tile rightTile;
	Tile topTile;
	Tile bottomTile;
	//left 
	if(y-1>=0 && tempAnswer[x*width+y-1]!=-1)
		leftTileId = tempAnswer[x*width+y-1];
	//right 
	if(y+1<width && tempAnswer[x*width+y+1]!=-1)
		rightTileId = tempAnswer[x*width+y+1];
	//top 
	if(x-1>=0 && tempAnswer[(x-1)*width+y]!=-1)
		topTileId = tempAnswer[(x-1)*width+y];
	//bottom 
	if(x+1<width && tempAnswer[(x+1)*width+y]!=-1)
		bottomTileId = tempAnswer[(x+1)*width+y];
	string shead;
	string rs,ls,ts,bs;
	
	// cout<<"x::::"<<x<<" y::::"<<y<<endl;
	// cout<<"l: "<<leftTileId<<" r: "<<rightTileId<< " t: "<<topTileId<< " p: "<<bottomTileId<<endl;
	if(leftTileId != -1)
	{
		shead = to_string(leftTileId)+"L-R";
		leftTile = tiles[leftTileId];
		map<string,int>::iterator iter;
	    iter = edgeVote.begin();
	    while(iter != edgeVote.end()) 
	    {
	    	string::size_type position;
	    	position = iter->first.find(shead);
	        if(position == 0)
	        {
	        	double sscore = 0;
	        	bool valid = true;
	        	int candidateId = stoi(iter->first.substr(shead.length()));
	        	// if(used[candidateId])
	        	// {
	        	// 	iter++;
	        	// 	continue;
	        	// }
	        	Tile canTile = tiles[candidateId];
	        	if(rightTileId!=-1)
	        	{
	        		rightTile = tiles[rightTileId];
	        		if(rightTile.left != canTile.right)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			rs = to_string(candidateId)+"L-R"+to_string(rightTileId);
	        			sscore+=edgeVote[rs];
	        		}

	        	}
	        	if(topTileId!=-1)
	        	{
	        		topTile = tiles[topTileId];
	        		if(topTile.bottom != canTile.top)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			ts = to_string(topTileId)+"T-B"+to_string(candidateId);
	        			sscore+=edgeVote[ts];
	        		}
	        	}
	        	if(bottomTileId!=-1)
	        	{
	        		bottomTile = tiles[bottomTileId];
	        		if(bottomTile.top != canTile.bottom)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			bs = to_string(candidateId)+"T-B"+to_string(bottomTileId);
	        			sscore+=edgeVote[bs];

	        		}
	        	}
	        	if(valid)
	        	{
	        		sscore += iter->second;
	        		candidates.insert(make_pair(candidateId,sscore));
	        	}

	        }
	        iter++;
	    }
	    
	}
	else if(topTileId != -1)
	{
		shead = to_string(topTileId)+"T-B";
		topTile = tiles[topTileId];
		map<string,int>::iterator iter;
	    iter = edgeVote.begin();
	    while(iter != edgeVote.end()) 
	    {
	    	string::size_type position;
	    	position = iter->first.find(shead);
	        if(position == 0)
	        {
	        	double sscore = 0;
	        	bool valid = true;
	        	int candidateId = stoi(iter->first.substr(shead.length()));
	        	// if(used[candidateId])
	        	// {
	        	// 	iter++;
	        	// 	continue;
	        	// }
	        	Tile canTile = tiles[candidateId];
	        	if(rightTileId!=-1)
	        	{
	        		rightTile = tiles[rightTileId];
	        		if(rightTile.left != canTile.right)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			rs = to_string(candidateId)+"L-R"+to_string(rightTileId);
	        			sscore+=edgeVote[rs];
	        		}

	        	}
	        	if(leftTileId!=-1)
	        	{
	        		leftTile = tiles[leftTileId];
	        		if(leftTile.right != canTile.left)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			ls = to_string(leftTileId)+"L-R"+to_string(candidateId);
	        			sscore+=edgeVote[ls];
	        		}
	        	}
	        	if(bottomTileId!=-1)
	        	{
	        		bottomTile = tiles[bottomTileId];
	        		if(bottomTile.top != canTile.bottom)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			bs = to_string(candidateId)+"T-B"+to_string(bottomTileId);
	        			sscore+=edgeVote[bs];

	        		}
	        	}
	        	if(valid)
	        	{
	        		sscore += iter->second;
	        		candidates.insert(make_pair(candidateId,sscore));
	        	}

	        }
	        iter++;
	    }
	}
	else if(rightTileId != -1)
	{
		shead = "L-R"+to_string(rightTileId);
		rightTile = tiles[rightTileId];
		map<string,int>::iterator iter;
	    iter = edgeVote.begin();
	    while(iter != edgeVote.end()) 
	    {
	    	string::size_type position;
	    	position = iter->first.find(shead);
	        if(position + shead.length() == iter->first.length())
	        {
	        	double sscore = 0;
	        	bool valid = true;
	        	int candidateId = stoi(iter->first.substr(0,position));
	        	// if(used[candidateId])
	        	// {
	        	// 	iter++;
	        	// 	continue;
	        	// }
	        	Tile canTile = tiles[candidateId];
	        	if(topTileId!=-1)
	        	{
	        		topTile = tiles[topTileId];
	        		if(topTile.bottom != canTile.top)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			ts = to_string(topTileId)+"T-B"+to_string(candidateId);
	        			sscore+=edgeVote[ts];
	        		}

	        	}
	        	if(leftTileId!=-1)
	        	{
	        		leftTile = tiles[leftTileId];
	        		if(leftTile.right != canTile.left)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			ls = to_string(leftTileId)+"L-R"+to_string(candidateId);
	        			sscore+=edgeVote[ls];
	        		}
	        	}
	        	if(bottomTileId!=-1)
	        	{
	        		bottomTile = tiles[bottomTileId];
	        		if(bottomTile.top != canTile.bottom)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			bs = to_string(candidateId)+"T-B"+to_string(bottomTileId);
	        			sscore+=edgeVote[bs];

	        		}
	        	}
	        	if(valid)
	        	{
	        		sscore += iter->second;
	        		candidates.insert(make_pair(candidateId,sscore));
	        	}

	        }
	        iter++;
	    }
	}
	else if(bottomTileId != -1)
	{
		shead = "T-B"+to_string(bottomTileId);
		bottomTile = tiles[bottomTileId];
		map<string,int>::iterator iter;
	    iter = edgeVote.begin();
	    while(iter != edgeVote.end()) 
	    {
	    	string::size_type position;
	    	position = iter->first.find(shead);
	        if(position + shead.length() == iter->first.length())
	        {
	        	double sscore = 0;
	        	bool valid = true;
	        	int candidateId = stoi(iter->first.substr(0,position));
	        	// if(used[candidateId])
	        	// {
	        	// 	iter++;
	        	// 	continue;
	        	// }
	        	Tile canTile = tiles[candidateId];
	        	if(topTileId!=-1)
	        	{
	        		topTile = tiles[topTileId];
	        		if(topTile.bottom != canTile.top)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			ts = to_string(topTileId)+"T-B"+to_string(candidateId);
	        			sscore+=edgeVote[ts];
	        		}

	        	}
	        	if(leftTileId!=-1)
	        	{
	        		leftTile = tiles[leftTileId];
	        		if(leftTile.right != canTile.left)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			ls = to_string(leftTileId)+"L-R"+to_string(candidateId);
	        			sscore+=edgeVote[ls];
	        		}
	        	}
	        	if(rightTileId!=-1)
	        	{
	        		rightTile = tiles[rightTileId];
	        		if(rightTile.left != canTile.right)
	        		{
	        			valid = false;
	        		}
	        		else
	        		{
	        			rs = to_string(candidateId)+"L-R"+to_string(rightTileId);
	        			sscore+=edgeVote[rs];

	        		}
	        	}
	        	if(valid)
	        	{
	        		sscore += iter->second;
	        		//cout<<"x:"<<x<<" y:"<<y<<" candidateId here: "<<candidateId<<endl;
	        		candidates.insert(make_pair(candidateId,sscore));
	        	}

	        }
	        iter++;
	    }
	}
	return NULL;
}
