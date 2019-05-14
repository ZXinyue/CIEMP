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
using namespace std;

#define POOL 10
#define RETURN_NUM 2


typedef struct
{
	int id;
	string top;
	string right;
	string bottom;
	string left;
	
}Tile;

int size = 25;
int width = 5;
Tile tiles[25];
map<string,int> edgeVote;
int* answer;

void* searchAnswer(void *id);
bool* BFS(int *tempAnswer, bool *flag, int *returnTime, bool *used, int pos);
void* backTrack(int *px, int *py, int *tempAnswer,int *returnTime, bool *used);
void* justTest(void *id);



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
				edgeVote[tempEdge] = 0;
			}
			if(t1.bottom == t2.top)
			{
				tempEdge = to_string(t1.id)+"T-B"+to_string(t2.id);
				edgeVote[tempEdge] = 0;
			}
			if(t1.left == t2.right)
			{
				tempEdge = to_string(t2.id)+"L-R"+to_string(t1.id);
				edgeVote[tempEdge] = 0;
			}
			if(t1.right == t2.left)
			{
				tempEdge = to_string(t1.id)+"L-R"+to_string(t2.id);
				edgeVote[tempEdge] = 0;
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
	cout<<"the thread"<<*(int*)id<<"randomStart is:  "<<randomStart<<endl;
	tempAnswer[0] = randomStart;
	used[randomStart] = true;
	
	while(1)
	{
		if(answer[0]!=-1)
			return NULL;

		if(tempAnswer[0] == -1)
		{
			seed = time(NULL);
			srand(seed+ *(int*)id * 12345);
			randomStart = rand()%size;
			tempAnswer[0] = randomStart;
			used[randomStart] = true;
			//cout<<"new randomStart "<< *(int*)id <<" is!--"<<randomStart<<endl;
		}
		bool *returnStatus = BFS(tempAnswer, flag, returnTime, used, 0);
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
		returnTime[(x-1)*width+y] += 1;
		tempAnswer[(x-1)*width+y] = -1;
		if(returnTime[(x-1)*width+y]>=RETURN_NUM)
		{
			//cout<<"aaaa"<<endl;
			int x_ = x-1;
			int *px_ = &x_;
			int *py_ = py;
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


bool* BFS(int *tempAnswer, bool *flag, int *returnTime, bool *used, int pos)
{
	//cout<<"hi BFS"<<endl;
	bool RS = true;
	bool *rs = &RS;
	int startx = pos/width;
	int starty = pos%width;
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

		cout<<"x"<<x<<"y"<<y<<endl;

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
		if(tempAnswer[x*width+y]!=-1)
			continue;
		int leftTileId = -1;
		int rightTileId = -1;
		int topTileId = -1;
		int bottomTileId = -1;

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
		cout<<"l: "<<leftTileId<<" r: "<<rightTileId<< " t: "<<topTileId<< " p: "<<bottomTileId<< endl;
		if(leftTileId!=-1)
		{
			shead = to_string(leftTileId)+"L-R";
			leftTile = tiles[leftTileId];
			multimap <int,int> candidates;
			map<string,int>::iterator iter;
		    iter = edgeVote.begin();
		    while(iter != edgeVote.end()) {
		    	string::size_type position;
		    	position = iter->first.find(shead);
		        if(position == 0)
		        {
		        	bool valid = true;
		        	int candidateId = stoi(iter->first.substr(shead.length()));
		        	if(used[candidateId])
		        	{
		        		iter++;
		        		continue;
		        	}
		        	Tile canTile = tiles[candidateId];
		        	if(rightTileId!=-1)
		        	{
		        		rightTile = tiles[rightTileId];
		        		if(rightTile.left != canTile.right)
		        		{
		        			valid = false;
		        		}

		        	}
		        	if(topTileId!=-1)
		        	{
		        		topTile = tiles[topTileId];
		        		if(topTile.bottom != canTile.top)
		        		{
		        			valid = false;
		        		}
		        	}
		        	if(bottomTileId!=-1)
		        	{
		        		bottomTile = tiles[bottomTileId];
		        		if(bottomTile.top != canTile.bottom)
		        		{
		        			valid = false;
		        		}
		        	}
		        	if(valid)
		        	{
		        		int sscore = iter->second;
		        		candidates.insert(make_pair(sscore,candidateId));
		        	}

		        }
		        iter++;
		    }
		    if(!candidates.empty())
		    {
		    	int winnerId = candidates.begin()->second;
		    	tempAnswer[x*width+y]=winnerId;
		    	used[winnerId] = true;
		    	cout<<"winnerId:"<<winnerId<<endl;
		    	int score = candidates.begin()->first+1;
		    	string tempEdge = shead+to_string(winnerId);
		    	pthread_mutex_lock(&mutex1);
		    	edgeVote[tempEdge] = score;
		    	//cout<<"add score"<<endl;
		    	pthread_mutex_unlock(&mutex1);
		    }
		    else
		    {
		    	// cout<<"failed tempAnswer is:"<<endl;
		    	// for(int i=0; i<size; i++)
		    	// 	cout<<tempAnswer[i]<<" ";
		    	// cout<<"edgevote this time:"<<endl;
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
		    	backTrack(px,py,tempAnswer,returnTime,used);




		    	*rs = false;
		    	return rs;

		    }
		}
		else if(topTileId!=-1)
		{
			shead = to_string(topTileId)+"T-B";
			topTile = tiles[topTileId];
			multimap <int,int> candidates;
			map<string,int>::iterator iter;
		    iter = edgeVote.begin();
		    while(iter != edgeVote.end()) {
		    	string::size_type position;
		    	position = iter->first.find(shead);
		    	cout<<"edge is"<<iter->first<<endl;
		        if(position == 0)
		        {
		        	bool valid = true;
		        	int candidateId = stoi(iter->first.substr(shead.length()));
		        	cout<<"candidateId"<<candidateId<<endl;
		        	if(used[candidateId])
		        	{
		        		iter++;
		        		continue;
		        	}
		        	Tile canTile = tiles[candidateId];
		        	if(rightTileId!=-1)
		        	{
		        		rightTile = tiles[rightTileId];
		        		if(rightTile.left != canTile.right)
		        		{
		        			valid = false;
		        		}

		        	}
		        	if(leftTileId!=-1)
		        	{
		        		leftTile = tiles[leftTileId];
		        		if(leftTile.right != canTile.left)
		        		{
		        			valid = false;
		        		}
		        	}
		        	if(bottomTileId!=-1)
		        	{
		        		bottomTile = tiles[bottomTileId];
		        		if(bottomTile.top != canTile.bottom)
		        		{
		        			valid = false;
		        		}
		        	}
		        	if(valid)
		        	{
		        		int sscore = iter->second;
		        		candidates.insert(make_pair(sscore,candidateId));
		        	}

		        }
		        iter++;
		    }
		    if(!candidates.empty())
		    {
		    	int winnerId = candidates.begin()->second;
		    	tempAnswer[x*width+y]=winnerId;
		    	cout<<"winnerId:"<<winnerId<<endl;
		    	used[winnerId] = true;
		    	int score = candidates.begin()->first+1;
		    	string tempEdge = shead+to_string(winnerId);
		    	pthread_mutex_lock(&mutex1);
		    	edgeVote[tempEdge] = score;
		    	//cout<<"add score"<<endl;
		    	pthread_mutex_unlock(&mutex1);
		    }
		    else
		    {
		    	// cout<<"failed tempAnswer is:"<<endl;
		    	// for(int i=0; i<size; i++)
		    	// 	cout<<tempAnswer[i]<<" ";

		    	// cout<<"edgevote this time:"<<endl;
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

		    	int *px = &x;
		    	int *py = &y;
		    	backTrack(px,py,tempAnswer,returnTime,used);


		    	*rs = false;
		    	return rs;

				
			}
		}


	}
	*rs = true;
	return rs;
}

