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
#include <unordered_map>
#include <queue>
#include <vector>
#include <set>
#include <functional>
#include <fstream>
#include <sys/time.h>
#include <algorithm>
using namespace std;

#define POOL 1000
#define TIMEOUT 7200

string inFileName = "data/9*9_c10_1";
string outFileName = "result/1021/9*9_c10_1_1000threads_5";

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
	bool picked; 
	bool d0; 
	bool d1;
	bool d2;
	bool d3;
	int squareSize;
	string topColor;
	string rightColor;
	string bottomColor;
	string leftColor;
	int way;//0:normal 1:pinpin
	bool pinPicked;
};

struct myComp  
{  
    bool operator()(const string &a,const string &b)  
    {  
        int anum = count(a.begin(),a.end(),'-');
        int bnum = count(b.begin(),b.end(),'-');
        return anum-bnum>0;  
    }  
};

struct sinSquare
{
	string key;
	set<int> sTiles;
};



unordered_map<string,squareInfo> squareMap;
vector<string> squareVector;
multiset<string,myComp> squareSet;
unordered_map<int,vector<string> > squareSizeMap;
pthread_mutex_t mutex1;
unordered_map<string,vector<sinSquare> > bigTiles;
pthread_mutex_t mutex2;

bool timeOut = false;
int* answer;
pthread_mutex_t mutex3;

int justCountPinpin;
int justCountNewPinpin;
int trytry[10]={0};
int success[10]={0};


int size = 81;
int width = 9;
Tile tiles[81];
//map<string,int> edgeVote;
set<string> invalidEdges;
unordered_map<string,vector<int> > validEdges;
unordered_map<string,unordered_map<int,int> > allEdges;




void* searchAnswer(void *id);

void* justTest(void *id);

int searchPosition(int pos, int theSquareSize, vector<int> squareIds, set<int> &usedTiles, vector< set<int> > &triedPos, vector<int> &posAns,int direction,string originKey);

void* pinPin(string randomStart);

//pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;



timeval t_start, t_end;



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
	
	//start = clock();
	gettimeofday( &t_start, NULL);

	justCountPinpin=0;
	justCountNewPinpin=0;

	for(int i=1; i<width; i++)
	{
		vector<string> emptyVec;
		squareSizeMap[i] = emptyVec;
	}

	ifstream input;
	input.open(inFileName);

	//读入数据并初始化
	for(int i=0; i<size; i++)
	{
		string s;
		input>>s;
		vector<string> colors;
		superSplit(s,colors,"-");

		tiles[i].id = i;
		tiles[i].top = colors[0];
		tiles[i].right = colors[1];
		tiles[i].bottom = colors[2];
		tiles[i].left = colors[3];

		string vs = to_string(s);
		squareInfo si = {vs,0,false,false,false,false,false,1,colors[0],colors[1],colors[2],colors[3],0,false};
		squareMap[vs] = si;
		squareVector.push_back(vs);
		squareSet.insert(vs);
		squareSizeMap[1].push_back(vs);
	}


	input.close();


	//存储所有边
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

	//遍历以找到不合法边 不合法边：无法扩展成为2*2的边
	unordered_map<string,unordered_map<int,int> > :: iterator it1;
	it1 = allEdges.begin();
	while(it1 != allEdges.end())
	{
		string tempS = it1->first;
		unordered_map<int,int> tempMap = it1->second;
		unordered_map<int,int> :: iterator it2;
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
				unordered_map<int,int> topLeftValidMap = allEdges[s1];
				unordered_map<int,int> topRightValidMap = allEdges[s2];
				unordered_map<int,int>::iterator tlvm, trvm;
				tlvm = topLeftValidMap.begin();
				
				while(tlvm != topLeftValidMap.end())
				{
					int topLeftId = tlvm->first;
					trvm = topRightValidMap.begin();
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
				unordered_map<int,int> bottomLeftValidMap = allEdges[s1];
				unordered_map<int,int> bottomRightValidMap = allEdges[s2];
				unordered_map<int,int>::iterator blvm, brvm;
				blvm = bottomLeftValidMap.begin();
				
				while(blvm != bottomLeftValidMap.end())
				{
					int bottomLeftId = blvm->first;
					brvm = bottomRightValidMap.begin();
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
				unordered_map<int,int> rightTopValidMap = allEdges[s1];
				unordered_map<int,int> rightBottomValidMap = allEdges[s2];
				unordered_map<int,int>::iterator rtvm, rbvm;
				rtvm = rightTopValidMap.begin();
				
				while(rtvm != rightTopValidMap.end())
				{
					int rightTopId = rtvm->first;
					rbvm = rightBottomValidMap.begin();
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
				unordered_map<int,int> leftTopValidMap = allEdges[s1];
				unordered_map<int,int> leftBottomValidMap = allEdges[s2];
				unordered_map<int,int>::iterator ltvm, lbvm;
				ltvm = leftTopValidMap.begin();
				
				while(ltvm != leftTopValidMap.end())
				{
					int leftTopId = ltvm->first;
					lbvm = leftBottomValidMap.begin();
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
		//cout<<"key:"<<tempS<<endl;
		unordered_map<int,int> tempMap = it1->second;
		unordered_map<int,int> :: iterator it2;
		it2 = tempMap.begin();
		vector<int> tempVec;
		while(it2 != tempMap.end())
		{
			if(it2->second != 0)
			{
				tempVec.push_back(it2->first);
				//cout<<"value:"<<it2->first<<" ";
				validEdgeCount+=1;
			}
			// else
			// 	cout<<"wrong value:"<<tempS<<":"<<it2->first<<" "<<it2->second<<endl;
			it2++;
		}
		//cout<<endl;
		//set<int> tempSet(tempVec.begin(),tempVec.end());
		validEdges[tempS] = tempVec;
		it1++;
	}

	cout<<"valid edge count:"<<validEdgeCount<<endl;


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
	pthread_mutex_init(&mutex2,NULL);
	pthread_mutex_init(&mutex3,NULL);

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
	pthread_mutex_destroy(&mutex2);
	pthread_mutex_destroy(&mutex3);

	if(answer[0]!=-1)
	{
		cout<<"the answer is:"<<endl;
		for(int j=0; j<size; j++)
		{
			cout<<answer[j]<<" ";
			//cout<<answer[j]<<":"<<tiles[answer[j]].top<<"-"<<tiles[answer[j]].right<<"-"<<tiles[answer[j]].bottom<<"-"<<tiles[answer[j]].left<<" ";
		}
		cout<<endl;
		for(int j=0; j<size ;j++)
		{
			if((j+1)%width!=0)
			{
				if(tiles[answer[j]].right!=tiles[answer[j+1]].left)
					cout<<"error l-r:"<<answer[j]<<" "<<answer[j+1]<<endl;
			}
			if(j<width*(width-1))
			{
				if(tiles[answer[j]].bottom!=tiles[answer[j+width]].top)
					cout<<"error t-b:"<<answer[j]<<" "<<answer[j+width]<<endl;
			}
		}

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
	//output << justCountPinpin <<" "<<justCountNewPinpin<<endl;
	//output << "succ 2:"<<success[2]<<" succ 3:"<<success[3]<<" succ 4:"<<success[4]<<" succ 5:"<<success[5]<<endl;
	//output << "try 2:"<<trytry[2]<<" try 3:"<<trytry[3]<<" try 4:"<<trytry[4]<<" try 5:"<<trytry[5]<<endl;

	//if(answer[0]==-1)
	{
		// map<string,squareInfo>::iterator it = squareMap.begin();
		// while(it != squareMap.end())
		// {
		// 	output<<it->second.s<<" "<<it->second.tryTime<<endl;
		// 	it++;
		// }
		/*
		sort(squareVector.begin(),squareVector.end());
		for(int i=0;i<squareVector.size();i++)
			output<<squareVector[i]<<endl;
		*/
		multiset<string,myComp>::iterator setit = squareSet.begin();
		while(setit != squareSet.end())
		{
		    output<<*setit<<endl;
		    //output<<"way"<<squareMap[*setit].way<<endl;
		    setit++;
		}
	}

	//just for test new character
	/*
	output<<"块有"<<endl;
	multiset<string,myComp>::iterator setit = squareSet.begin();
	while(setit != squareSet.end())
	{
	    output<<*setit<<endl;
	    setit++;
	}
	output<<"色边"<<endl;
	unordered_map<string,vector<sinSquare> >::iterator iitt = bigTiles.begin();
	while(iitt!=bigTiles.end())
	{
		output<<iitt->first<<endl;
		for(int i=0;i<iitt->second.size();i++)
		{
			output<<"key  "<<iitt->second[i].key<<endl;
			set<int>::iterator itt = iitt->second[i].sTiles.begin();
			while(itt!=iitt->second[i].sTiles.end())
			{
				output<<*itt<<",";
				itt++;
			}
			output<<endl;
		}
		iitt++;
	}*/


	output.close();

	free(threads);
	free(answer);

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
	seed = rand();
	while(answer[0]==-1 && !timeOut)
	{
		// pthread_mutex_lock(&mutex3);
		// cout<<"thread "<<*(int*)id<<"return loop"<<endl;
		// pthread_mutex_unlock(&mutex3);
		timeval t_temp;
		gettimeofday( &t_temp, NULL);
		double theTime = (t_temp.tv_sec-t_start.tv_sec) + (t_temp.tv_usec-t_start.tv_usec)/1000000.0;
		if(theTime > TIMEOUT)
		{
			pthread_mutex_lock(&mutex2);
			timeOut = true;
			cout<<"thread "<<*(int*)id<<"find time out!"<<endl;
			pthread_mutex_unlock(&mutex2);
			break;
		}
		string randomStart = "";
		bool ffoundFlag = false;
		string vs;
		int randomNumber;

		bool changeMode = false;
		
		// if(*(int*)id % 2 == 0) //1/2
		// {
		// 	int minTryTime = -1;
		// 	//int yu = *(int*)id % (width/2-1) +2;
		// 	int max_time = 10;
		// 	srand(seed);
		// 	int yu = rand()%(width/2-1) +2;
		// 	seed = rand();


		// 	pthread_mutex_lock(&mutex1);
		// 	int theYuSize = squareSizeMap[yu].size();
		// 	pthread_mutex_unlock(&mutex1);

			
		// 	// pthread_mutex_lock(&mutex2);
		// 	// cout<<"id is :"<<*(int*)id<<endl;
		// 	// cout<<"yu:"<<yu<<" yusize:"<<theYuSize<<endl;
		// 	// pthread_mutex_unlock(&mutex2);
			
			

		// 	if(theYuSize >= 4)
		// 	{
		// 		pthread_mutex_lock(&mutex1);
		// 		while(max_time--) 
		// 		{
		// 			srand(seed);
		// 			randomNumber = rand()%squareSizeMap[yu].size();
		// 			seed = rand();
		// 			vs = squareSizeMap[yu][randomNumber];
		// 			if(squareMap[vs].pinPicked)
		// 			{
		// 				continue;
		// 			}
		// 			if(minTryTime == -1)
		// 			{
		// 				minTryTime = squareMap[vs].tryTime;
		// 				randomStart = vs;
		// 			}
		// 			else if(squareMap[vs].tryTime < minTryTime)
		// 			{
		// 				minTryTime = squareMap[vs].tryTime;
		// 				randomStart = vs;
		// 			}
		// 		}
		// 		pthread_mutex_unlock(&mutex1);
		// 		if(randomStart!="")
		// 		{
		// 			pthread_mutex_lock(&mutex1);
		// 			squareMap[randomStart].tryTime++;
		// 			squareMap[randomStart].pinPicked = true;
		// 			pthread_mutex_unlock(&mutex1);
		// 			pinPin(randomStart);
		// 		}
		// 		else
		// 			changeMode = true;			

		// 	}
		// 	else
		// 		changeMode = true;	
		// }
		//if(*(int*)id % 2 == 1 || changeMode)
		{
			pthread_mutex_lock(&mutex1);
			while(!ffoundFlag) //随机选取一个碎块作为扩展的起始
			{
				srand(seed);
				randomNumber = rand()%squareVector.size();
				seed = rand();
				vs = squareVector[randomNumber];
				if(!squareMap[vs].picked)
				{
					randomStart = vs;
					ffoundFlag = true;
					squareMap[vs].picked = true;
				}
			}
			pthread_mutex_unlock(&mutex1);

			vector<string> splitResult;
			set<int> usedTiles;
			vector<int> squareIds;
			superSplit(randomStart,splitResult,"-");
			// if(splitResult.size()==0)
			// {
			// 	cout<<"randomStart:"<<randomStart<<endl;
			// 	cout<<"vs:"<<vs<<endl;
			// 	cout<<"randomNumber:"<<randomNumber<<endl;
			// 	cout<<"squ[i]:"<<squareVector[randomNumber]<<endl;
			// }
			
			for(int i=0;i<splitResult.size();i++)
			{
				//cout<<"superSplit:"<<splitResult[i]<<endl;
				usedTiles.insert(stoi(splitResult[i]));
				squareIds.push_back(stoi(splitResult[i]));
			}
			

			int theSqSize = sqrt(squareIds.size());

			vector< set<int> > triedPos;//存储每个位置已经尝试过的块的id
			vector<int> posAns;//存储此次扩展的当前结果
			for(int i=0; i<theSqSize*2+1; i++)
			{
				set<int> empty;
				triedPos.push_back(empty);
				posAns.push_back(-1);
			}

			int finish = 0;//0:failed 1:success 2:finish
			int POS = 0;
			if(!squareMap[randomStart].d0)
			{
				
				while(finish!=2)
				{
					finish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,0,randomStart);
					if(finish == 1)
						POS+=1;
					else if(finish == 0)
						POS-=1;
					else
						break;
				}
				
			}
			finish = 0;
			POS = 0;
			triedPos.clear();
			posAns.clear();
			for(int i=0; i<theSqSize*2+1; i++)
			{
				set<int> empty;
				triedPos.push_back(empty);
				posAns.push_back(-1);
			}
			if(!squareMap[randomStart].d1)
			{
				while(finish!=2)
				{
					finish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,1,randomStart);
					if(finish == 1)
						POS+=1;
					else if(finish == 0)
						POS-=1;
					else
						break;
				}
			}
			finish = 0;
			POS = 0;
			triedPos.clear();
			posAns.clear();
			for(int i=0; i<theSqSize*2+1; i++)
			{
				set<int> empty;
				triedPos.push_back(empty);
				posAns.push_back(-1);
			}
			if(!squareMap[randomStart].d2)
			{
				while(finish!=2)
				{
					finish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,2,randomStart);
					if(finish == 1)
						POS+=1;
					else if(finish == 0)
						POS-=1;
					else
						break;
				}
			}
			finish = 0;
			POS = 0;
			triedPos.clear();
			posAns.clear();
			for(int i=0; i<theSqSize*2+1; i++)
			{
				set<int> empty;
				triedPos.push_back(empty);
				posAns.push_back(-1);
			}
			if(!squareMap[randomStart].d3)
			{
				while(finish!=2)
				{
					finish = searchPosition(POS,theSqSize,squareIds,usedTiles,triedPos,posAns,3,randomStart);
					if(finish == 1)
						POS+=1;
					else if(finish == 0)
						POS-=1;
					else
						break;
				}
			}
			
		}

		
	}
	
	return NULL;
}
//0:failed 1:success 2:finish
int searchPosition(int pos, int theSquareSize, vector<int> squareIds, set<int> &usedTiles, vector< set<int> > &triedPos, vector<int> &posAns,int direction,string originKey)
{
	for(int i=pos+1; i<2*theSquareSize+1; i++)
	{
		if(!triedPos[i].empty())
			triedPos[i].clear(); 
	}
	bool found = false;
	if(direction == 0)  //——|
	{
		if(pos == 0)
		{
			vector<int> candiVector = validEdges[to_string(squareIds[pos])+"B-T"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int leftId = posAns[pos-1];
			vector<int> v1 = validEdges[to_string(leftId)+"L-R"];
			vector<int> v2 = validEdges[to_string(squareIds[pos])+"B-T"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}

		}
		else if(pos == theSquareSize)
		{
			int leftId = posAns[pos-1];
			vector<int> candiVector = validEdges[to_string(leftId)+"L-R"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int topId = posAns[pos-1];
			int relativeId = squareIds[(pos-theSquareSize-1)*theSquareSize+theSquareSize-1];
			vector<int> v1 = validEdges[to_string(topId)+"T-B"];
			vector<int> v2 = validEdges[to_string(relativeId)+"L-R"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}


		if(found) 
		{
			
			if(pos<2*theSquareSize) 
			{
				//returnStatus = 1;
				//int *pRe = new int(1);
				//return pRe;
				return 1;
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
				set<int> sTiles;
				for(int i=0; i<newSquare.size()-1; i++)
				{
					key = key + to_string(newSquare[i]) + "-";
					sTiles.insert(newSquare[i]);
				}
				key+=to_string(newSquare[newSquare.size()-1]);
				sTiles.insert(newSquare[newSquare.size()-1]);

				pthread_mutex_lock(&mutex1);
				bool notFound = (squareMap.find(key)==squareMap.end())?true:false;
				pthread_mutex_unlock(&mutex1);

				if(notFound)
				{

					sinSquare theSQ = {key,sTiles};
					vector<sinSquare> tempSinVec;
					tempSinVec.push_back(theSQ);
					//top color
					string topColor = "T:";
					for(int i=0; i<theSquareSize; i++)
					{
						topColor += tiles[newSquare[i]].top+"-";
					}
					topColor += tiles[newSquare[theSquareSize]].top;

					//right color
					string rightColor = "R:";
					for(int i=0; i<theSquareSize; i++)
					{
						rightColor += tiles[newSquare[i*(theSquareSize+1)+theSquareSize]].right+"-";
					}
					rightColor += tiles[newSquare[(theSquareSize+1)*(theSquareSize+1)-1]].right;

					//bottom color
					string bottomColor = "B:";
					for(int i=0; i<theSquareSize; i++)
					{
						bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+i]].bottom+"-";
					}
					bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+theSquareSize]].bottom;

					//left color
					string leftColor = "L:";
					for(int i=0; i<theSquareSize; i++)
					{
						leftColor += tiles[newSquare[i*(theSquareSize+1)]].left+"-";
					}
					leftColor += tiles[newSquare[theSquareSize*(theSquareSize+1)]].left;

					squareInfo tempInfo = {key,0,false,false,false,false,false,theSquareSize+1,topColor.substr(2),rightColor.substr(2),bottomColor.substr(2),leftColor.substr(2),0,false};



					pthread_mutex_lock(&mutex1);
					squareVector.push_back(key);
					squareSet.insert(key);
					squareSizeMap[theSquareSize+1].push_back(key);
					squareMap[key] = tempInfo;
					

					if(bigTiles.count(topColor) == 0)
						bigTiles[topColor] = tempSinVec;
					else bigTiles[topColor].push_back(theSQ);

			
					if(bigTiles.count(rightColor) == 0)
						bigTiles[rightColor] = tempSinVec;
					else bigTiles[rightColor].push_back(theSQ);
					
					
					if(bigTiles.count(bottomColor) == 0)
						bigTiles[bottomColor] = tempSinVec;
					else bigTiles[bottomColor].push_back(theSQ);
					

					if(bigTiles.count(leftColor) == 0)
						bigTiles[leftColor] = tempSinVec;
					else bigTiles[leftColor].push_back(theSQ);

					pthread_mutex_unlock(&mutex1);

				}

				

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex3);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex3);
					//returnStatus = 2;
					// int *pRe = new int(2);
					// return pRe;
					return 2;
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
				// int *pRe = new int(0);
				// return pRe;
				return 0;
				// searchPosition(pos-1,theSquareSize,squareIds,usedTiles,triedPos,posAns,0,originKey);
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d0 = true;
				//squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				// int *pRe = new int(2);
				// return pRe;
				return 2;
			}
			
		}
	}
	else if(direction == 1)
	{
		if(pos == 0)
		{
			int relativeId = squareIds[pos*theSquareSize+theSquareSize-1];
			vector<int> candiVector = validEdges[to_string(relativeId)+"L-R"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				//cout<<"put "<<i<<" at pos "<<pos<<endl;
				break;
				
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int topId = posAns[pos-1];
			int relativeId = squareIds[pos*theSquareSize+theSquareSize-1];
			vector<int> v1 = validEdges[to_string(topId)+"T-B"];
			vector<int> v2 = validEdges[to_string(relativeId)+"L-R"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集

			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}

		}
		else if(pos == theSquareSize)
		{
			int topId = posAns[pos-1];
			vector<int> candiVector = validEdges[to_string(topId)+"T-B"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int rightId = posAns[pos-1];
			int relativeId = squareIds[(theSquareSize-1)*theSquareSize+theSquareSize-1-(pos-theSquareSize-1)];
			vector<int> v1 = validEdges[to_string(rightId)+"R-L"];
			vector<int> v2 = validEdges[to_string(relativeId)+"T-B"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集

			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}

		if(found)
		{
			
			if(pos<2*theSquareSize)
			{
				//returnStatus = 1;
				// int *pRe = new int(1);
				// return pRe;
				return 1;
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
				set<int> sTiles;
				for(int i=0; i<newSquare.size()-1; i++)
				{
					key = key + to_string(newSquare[i]) + "-";
					sTiles.insert(newSquare[i]);
				}
				key+=to_string(newSquare[newSquare.size()-1]);
				sTiles.insert(newSquare[newSquare.size()-1]);

				pthread_mutex_lock(&mutex1);
				bool notFound = (squareMap.find(key)==squareMap.end())?true:false;
				pthread_mutex_unlock(&mutex1);

				if(notFound)
				{

					sinSquare theSQ = {key,sTiles};
					vector<sinSquare> tempSinVec;
					tempSinVec.push_back(theSQ);
					//top color
					string topColor = "T:";
					for(int i=0; i<theSquareSize; i++)
					{
						topColor += tiles[newSquare[i]].top+"-";
					}
					topColor += tiles[newSquare[theSquareSize]].top;

					//right color
					string rightColor = "R:";
					for(int i=0; i<theSquareSize; i++)
					{
						rightColor += tiles[newSquare[i*(theSquareSize+1)+theSquareSize]].right+"-";
					}
					rightColor += tiles[newSquare[(theSquareSize+1)*(theSquareSize+1)-1]].right;

					//bottom color
					string bottomColor = "B:";
					for(int i=0; i<theSquareSize; i++)
					{
						bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+i]].bottom+"-";
					}
					bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+theSquareSize]].bottom;

					//left color
					string leftColor = "L:";
					for(int i=0; i<theSquareSize; i++)
					{
						leftColor += tiles[newSquare[i*(theSquareSize+1)]].left+"-";
					}
					leftColor += tiles[newSquare[theSquareSize*(theSquareSize+1)]].left;

					squareInfo tempInfo = {key,0,false,false,false,false,false,theSquareSize+1,topColor.substr(2),rightColor.substr(2),bottomColor.substr(2),leftColor.substr(2),0,false};

					pthread_mutex_lock(&mutex1);
					squareVector.push_back(key);
					squareSet.insert(key);
					squareSizeMap[theSquareSize+1].push_back(key);
					squareMap[key] = tempInfo;
					
					if(bigTiles.count(topColor) == 0)
						bigTiles[topColor] = tempSinVec;
					else bigTiles[topColor].push_back(theSQ);

					
					if(bigTiles.count(rightColor) == 0)
						bigTiles[rightColor] = tempSinVec;
					else bigTiles[rightColor].push_back(theSQ);

					
					if(bigTiles.count(bottomColor) == 0)
						bigTiles[bottomColor] = tempSinVec;
					else bigTiles[bottomColor].push_back(theSQ);

					
					if(bigTiles.count(leftColor) == 0)
						bigTiles[leftColor] = tempSinVec;
					else bigTiles[leftColor].push_back(theSQ);
					pthread_mutex_unlock(&mutex1);
		
				}

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex3);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex3);
					//return NULL;
					//returnStatus = 2;
					// int *pRe = new int(2);
					// return pRe;
					return 2;
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
				// int *pRe = new int(0);
				// return pRe;
				return 0;
				//searchPosition(pos-1,theSquareSize,squareIds,usedTiles,triedPos,posAns,1,originKey);
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d1 = true;
				//squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				// int *pRe = new int(2);
				// return pRe;
				return 2;
			}
			
		}
	}
	else if(direction == 2) //|_>
	{
		if(pos == 0)
		{
			vector<int> candiVector = validEdges[to_string(squareIds[pos])+"R-L"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				//cout<<"put "<<i<<" at pos "<<pos<<endl;
				break;
				
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int topId = posAns[pos-1];
			int relativeId = squareIds[pos*theSquareSize];
			vector<int> v1 = validEdges[to_string(topId)+"T-B"];
			vector<int> v2 = validEdges[to_string(relativeId)+"R-L"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}

		}
		else if(pos == theSquareSize)
		{
			int topId = posAns[pos-1];
			vector<int> candiVector = validEdges[to_string(topId)+"T-B"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int leftId = posAns[pos-1];
			int relativeId = squareIds[(theSquareSize-1)*theSquareSize+pos-theSquareSize-1];
			vector<int> v1 = validEdges[to_string(leftId)+"L-R"];
			vector<int> v2 = validEdges[to_string(relativeId)+"T-B"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}
		if(found)
		{
			
			if(pos<2*theSquareSize)
			{
				//returnStatus = 1;
				// int *pRe = new int(1);
				// return pRe;
				return 1;
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
				set<int> sTiles;
				for(int i=0; i<newSquare.size()-1; i++)
				{
					key = key + to_string(newSquare[i]) + "-";
					sTiles.insert(newSquare[i]);
				}
				key+=to_string(newSquare[newSquare.size()-1]);
				sTiles.insert(newSquare[newSquare.size()-1]);

				pthread_mutex_lock(&mutex1);
				bool notFound = (squareMap.find(key)==squareMap.end())?true:false;
				pthread_mutex_unlock(&mutex1);

				if(notFound)
				{
					
					sinSquare theSQ = {key,sTiles};
					vector<sinSquare> tempSinVec;
					tempSinVec.push_back(theSQ);
					//top color
					string topColor = "T:";
					for(int i=0; i<theSquareSize; i++)
					{
						topColor += tiles[newSquare[i]].top+"-";
					}
					topColor += tiles[newSquare[theSquareSize]].top;

					//right color
					string rightColor = "R:";
					for(int i=0; i<theSquareSize; i++)
					{
						rightColor += tiles[newSquare[i*(theSquareSize+1)+theSquareSize]].right+"-";
					}
					rightColor += tiles[newSquare[(theSquareSize+1)*(theSquareSize+1)-1]].right;

					//bottom color
					string bottomColor = "B:";
					for(int i=0; i<theSquareSize; i++)
					{
						bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+i]].bottom+"-";
					}
					bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+theSquareSize]].bottom;

					//left color
					string leftColor = "L:";
					for(int i=0; i<theSquareSize; i++)
					{
						leftColor += tiles[newSquare[i*(theSquareSize+1)]].left+"-";
					}
					leftColor += tiles[newSquare[theSquareSize*(theSquareSize+1)]].left;

					squareInfo tempInfo = {key,0,false,false,false,false,false,theSquareSize+1,topColor.substr(2),rightColor.substr(2),bottomColor.substr(2),leftColor.substr(2),0,false};


					pthread_mutex_lock(&mutex1);
					squareVector.push_back(key);
					squareSizeMap[theSquareSize+1].push_back(key);
					squareSet.insert(key);
					squareMap[key] = tempInfo;
					
					if(bigTiles.count(topColor) == 0)
						bigTiles[topColor] = tempSinVec;
					else bigTiles[topColor].push_back(theSQ);

					
					if(bigTiles.count(rightColor) == 0)
						bigTiles[rightColor] = tempSinVec;
					else bigTiles[rightColor].push_back(theSQ);

					
					if(bigTiles.count(bottomColor) == 0)
						bigTiles[bottomColor] = tempSinVec;
					else bigTiles[bottomColor].push_back(theSQ);

					
					if(bigTiles.count(leftColor) == 0)
						bigTiles[leftColor] = tempSinVec;
					else bigTiles[leftColor].push_back(theSQ);
					pthread_mutex_unlock(&mutex1);

					
					

					
				}

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex3);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex3);
					//return NULL;
					//returnStatus = 2;
					// int *pRe = new int(2);
					// return pRe;
					return 2;
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
				// int *pRe = new int(0);
				// return pRe;
				return 0;
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d2 = true;
				//squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				// int *pRe = new int(2);
				// return pRe;
				return 2;
			}
			
		}

	}
	else if(direction == 3) //|<-
	{
		if(pos == 0)
		{
			int relativeId = squareIds[theSquareSize-1-pos];
			vector<int> candiVector = validEdges[to_string(relativeId)+"B-T"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				//cout<<"put "<<i<<" at pos "<<pos<<endl;
				break;
				
			}
		}
		else if(pos>0 && pos<theSquareSize)
		{
			int rightId = posAns[pos-1];
			int relativeId = squareIds[theSquareSize-1-pos];
			vector<int> v1 = validEdges[to_string(rightId)+"R-L"];
			vector<int> v2 = validEdges[to_string(relativeId)+"B-T"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}

		}
		else if(pos == theSquareSize)
		{
			int rightId = posAns[pos-1];
			vector<int> candiVector = validEdges[to_string(rightId)+"R-L"];
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}
		else if(pos<2*theSquareSize+1)
		{
			int topId = posAns[pos-1];
			int relativeId = squareIds[(pos-theSquareSize-1)*theSquareSize];
			vector<int> v1 = validEdges[to_string(topId)+"T-B"];
			vector<int> v2 = validEdges[to_string(relativeId)+"R-L"];
			vector<int> candiVector;
			sort(v1.begin(),v1.end());   
			sort(v2.begin(),v2.end());   
			set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(candiVector));//求交集
			for(int iter=0; iter<candiVector.size(); iter++)
			{
				int i = candiVector[iter];
				if(usedTiles.count(i) == 1 || triedPos[pos].count(i) == 1)
					continue;
				
				found = true;
				posAns[pos] = i;
				usedTiles.insert(i);
				triedPos[pos].insert(i);
				break;
				
			}
		}

		if(found)
		{
			
			if(pos<2*theSquareSize)
			{
				//returnStatus = 1;
				// int *pRe = new int(1);
				// return pRe;
				return 1;
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
				set<int> sTiles;
				for(int i=0; i<newSquare.size()-1; i++)
				{
					key = key + to_string(newSquare[i]) + "-";
					sTiles.insert(newSquare[i]);
				}
				key+=to_string(newSquare[newSquare.size()-1]);
				sTiles.insert(newSquare[newSquare.size()-1]);

				pthread_mutex_lock(&mutex1);
				bool notFound = (squareMap.find(key)==squareMap.end())?true:false;
				pthread_mutex_unlock(&mutex1);

				if(notFound)
				{
					sinSquare theSQ = {key,sTiles};
					vector<sinSquare> tempSinVec;
					tempSinVec.push_back(theSQ);
					//top color
					string topColor = "T:";
					for(int i=0; i<theSquareSize; i++)
					{
						topColor += tiles[newSquare[i]].top+"-";
					}
					topColor += tiles[newSquare[theSquareSize]].top;

					//right color
					string rightColor = "R:";
					for(int i=0; i<theSquareSize; i++)
					{
						rightColor += tiles[newSquare[i*(theSquareSize+1)+theSquareSize]].right+"-";
					}
					rightColor += tiles[newSquare[(theSquareSize+1)*(theSquareSize+1)-1]].right;

					//bottom color
					string bottomColor = "B:";
					for(int i=0; i<theSquareSize; i++)
					{
						bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+i]].bottom+"-";
					}
					bottomColor += tiles[newSquare[theSquareSize*(theSquareSize+1)+theSquareSize]].bottom;

					//left color
					string leftColor = "L:";
					for(int i=0; i<theSquareSize; i++)
					{
						leftColor += tiles[newSquare[i*(theSquareSize+1)]].left+"-";
					}
					leftColor += tiles[newSquare[theSquareSize*(theSquareSize+1)]].left;

					squareInfo tempInfo = {key,0,false,false,false,false,false,theSquareSize+1,topColor.substr(2),rightColor.substr(2),bottomColor.substr(2),leftColor.substr(2),0,false};
					
					pthread_mutex_lock(&mutex1);
					squareVector.push_back(key);
					squareSet.insert(key);
					squareSizeMap[theSquareSize+1].push_back(key);
					squareMap[key] = tempInfo;
					
					if(bigTiles.count(topColor) == 0)
						bigTiles[topColor] = tempSinVec;
					else bigTiles[topColor].push_back(theSQ);

					
					if(bigTiles.count(rightColor) == 0)
						bigTiles[rightColor] = tempSinVec;
					else bigTiles[rightColor].push_back(theSQ);

					
					if(bigTiles.count(bottomColor) == 0)
						bigTiles[bottomColor] = tempSinVec;
					else bigTiles[bottomColor].push_back(theSQ);

					
					if(bigTiles.count(leftColor) == 0)
						bigTiles[leftColor] = tempSinVec;
					else bigTiles[leftColor].push_back(theSQ);
					pthread_mutex_unlock(&mutex1);
					
				}

				if(newSquare.size()==size)
				{
					pthread_mutex_lock(&mutex3);
					for(int i=0; i<newSquare.size(); i++)
					{
						answer[i]=newSquare[i];
					}
					pthread_mutex_unlock(&mutex3);
					//return NULL;
					//returnStatus = 2;
					// int *pRe = new int(2);
					// return pRe;
					return 2;
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
				// int *pRe = new int(0);
				// return pRe;
				return 0;
			}
			else
			{
				pthread_mutex_lock(&mutex1);
				squareMap[originKey].d3 = true;
				//squareMap[originKey].tryTime+=1;
				pthread_mutex_unlock(&mutex1);
				//return NULL;
				//returnStatus = 2;
				// int *pRe = new int(2);
				// return pRe;
				return 2;
			}
			
		}

	}
}

void* pinPin(string randomStart)
{
	// pthread_mutex_lock(&mutex1);
	// cout<<"start with:"<<randomStart<<endl;
	// pthread_mutex_unlock(&mutex1);
	if(squareMap.count(randomStart)==0)
		return NULL;

	pthread_mutex_lock(&mutex1);
	squareInfo startInfo = squareMap[randomStart];
	pthread_mutex_unlock(&mutex1);

	int ppsize = startInfo.squareSize;

	pthread_mutex_lock(&mutex3);
	trytry[ppsize]+=1;
	pthread_mutex_unlock(&mutex3);


	vector<string> splitResult1;
	set<int> startSet;
	superSplit(randomStart,splitResult1,"-");
	for(int i=0; i<splitResult1.size();i++)
		startSet.insert(stoi(splitResult1[i]));

	string rightBigKey = "L:"+startInfo.rightColor;
	string bottomBigKey = "T:"+startInfo.bottomColor;
	if(bigTiles.count(rightBigKey)==1 && bigTiles.count(bottomBigKey)==1)
	{
		pthread_mutex_lock(&mutex1);
		vector<sinSquare> rightTempCan = bigTiles[rightBigKey];
		pthread_mutex_unlock(&mutex1);

		vector<sinSquare> rightCan;
		for(int i=0;i<rightTempCan.size();i++)
		{
			if(timeOut)
			{
				// pthread_mutex_lock(&mutex2);
				// cout<<" righttempcan size:"<<rightTempCan.size()<<endl;
				// pthread_mutex_unlock(&mutex2);
				return NULL;
			}
			// else
			// {
			// 	pthread_mutex_lock(&mutex2);
			// 	cout<<"loop 1"<<endl;
			// 	pthread_mutex_unlock(&mutex2);
			// }
			set<int> tempSet1 = rightTempCan[i].sTiles;
			set<int> inter_Ans1;
			set_intersection(tempSet1.begin(),tempSet1.end(),startSet.begin(),startSet.end(),inserter(inter_Ans1,inter_Ans1.begin()));
			if(inter_Ans1.size()==0)
				rightCan.push_back(rightTempCan[i]);
		}

		// pthread_mutex_lock(&mutex1);
		// cout<<"rightBigKey:"<<rightBigKey<<" rightrawcan size:"<<rightTempCan.size()<<" rightcan size:"<<rightCan.size()<<endl;
		// pthread_mutex_unlock(&mutex1);
		pthread_mutex_lock(&mutex1);
		vector<sinSquare> bottomTempCan = bigTiles[bottomBigKey];
		pthread_mutex_unlock(&mutex1);

		vector<sinSquare> bottomCan;
		for(int i=0;i<bottomTempCan.size();i++)
		{
			if(timeOut)
			{
				// pthread_mutex_lock(&mutex2);
				// cout<<" bottomtempcan size:"<<bottomTempCan.size()<<endl;
				// pthread_mutex_unlock(&mutex2);
				return NULL;
			}
			// else
			// {
			// 	pthread_mutex_lock(&mutex2);
			// 	cout<<"loop 2"<<endl;
			// 	pthread_mutex_unlock(&mutex2);
			// }
			set<int> tempSet2 = bottomTempCan[i].sTiles;
			set<int> inter_Ans2;
			set_intersection(tempSet2.begin(),tempSet2.end(),startSet.begin(),startSet.end(),inserter(inter_Ans2,inter_Ans2.begin()));
			if(inter_Ans2.size()==0)
				bottomCan.push_back(bottomTempCan[i]);
		}

		// pthread_mutex_lock(&mutex1);
		// cout<<"bottomBigKey:"<<bottomBigKey<<" bottomrawcan size:"<<bottomTempCan.size()<<" bottomcan size:"<<bottomCan.size()<<endl;
		// pthread_mutex_unlock(&mutex1);

		string finalTopKey;
		string finalLeftKey;
		vector<sinSquare> finalTop;
		vector<sinSquare> finalLeft;

		for(int i=0;i<rightCan.size();i++)
			for(int j=0; j<bottomCan.size();j++)
			{
				if(timeOut)
				{
					// pthread_mutex_lock(&mutex2);
					// cout<<" bottomcan size:"<<bottomCan.size()<<" rightcan size:"<<rightCan.size()<<endl;
					// pthread_mutex_unlock(&mutex2);
					return NULL;
				}
				// else
				// {
				// 	pthread_mutex_lock(&mutex2);
				// 	cout<<"loop 3"<<endl;
				// 	pthread_mutex_unlock(&mutex2);
				// }
				set<int> rightSet = rightCan[i].sTiles;
				set<int> bottomSet = bottomCan[j].sTiles;
				set<int> interAns;
				set_intersection(rightSet.begin(),rightSet.end(),bottomSet.begin(),bottomSet.end(),inserter(interAns,interAns.begin()));
				
				if(interAns.size()==0)
				{
					finalTopKey = "T:"+squareMap[rightCan[i].key].bottomColor;
					finalLeftKey = "L:"+squareMap[bottomCan[j].key].rightColor;
					
					if(bigTiles.count(finalTopKey)==1 && bigTiles.count(finalLeftKey)==1)
					{
						pthread_mutex_lock(&mutex1);
						finalTop = bigTiles[finalTopKey];
						finalLeft = bigTiles[finalLeftKey];
						pthread_mutex_unlock(&mutex1);

						vector<sinSquare> finalCan;
						for(int u=0;u<finalTop.size();u++)
							for(int v=0;v<finalLeft.size();v++)
							{
								if(finalTop[u].key == finalLeft[v].key)
									finalCan.push_back(finalTop[u]);
							}
						
						for(int w=0;w<finalCan.size();w++)
						{
							if(timeOut)
							{
								// pthread_mutex_lock(&mutex2);
								// cout<<" bottomcan size:"<<bottomCan.size()<<" rightcan size:"<<rightCan.size()<<" finalcan size:"<<finalCan.size()<<endl;
								// pthread_mutex_unlock(&mutex2);
								return NULL;
							}
							// else
							// {
							// 	pthread_mutex_lock(&mutex2);
							// 	cout<<"loop 4"<<endl;
							// 	pthread_mutex_unlock(&mutex2);
							// }
							
							set<int> finalSet = finalCan[w].sTiles;
							set<int> interAns1;
							set<int> interAns2;
							set<int> interAns3;
							set_intersection(finalSet.begin(),finalSet.end(),startSet.begin(),startSet.end(),inserter(interAns1,interAns1.begin()));
							set_intersection(finalSet.begin(),finalSet.end(),rightSet.begin(),rightSet.end(),inserter(interAns2,interAns2.begin()));
							set_intersection(finalSet.begin(),finalSet.end(),bottomSet.begin(),bottomSet.end(),inserter(interAns3,interAns3.begin()));
							
							if(interAns1.size()==0 && interAns2.size()==0 && interAns3.size()==0)
							{
								pthread_mutex_lock(&mutex3);
								justCountPinpin++;
								success[ppsize]+=1;
								pthread_mutex_unlock(&mutex3);

								string no2 = rightCan[i].key;
								string no3 = bottomCan[j].key;
								string no4 = finalCan[w].key;

								vector<string> splitResult2;
								superSplit(no2,splitResult2,"-");

								vector<string> splitResult3;
								superSplit(no3,splitResult3,"-");

								vector<string> splitResult4;
								superSplit(no4,splitResult4,"-");

								string theBigAnsKey = "";
								set<int> theBigSet;

								for(int x=0; x<ppsize; x++)
								{
									for(int y=0; y<ppsize; y++)
									{
										theBigAnsKey += splitResult1[x*ppsize+y]+"-";
										theBigSet.insert(stoi(splitResult1[x*ppsize+y]));
									}
									for(int z=0; z<ppsize; z++)
									{
										theBigAnsKey += splitResult2[x*ppsize+z]+"-";
										theBigSet.insert(stoi(splitResult2[x*ppsize+z]));
									}
								}
								for(int x=0; x<ppsize; x++)
								{
									for(int y=0; y<ppsize; y++)
									{
										theBigAnsKey += splitResult3[x*ppsize+y]+"-";
										theBigSet.insert(stoi(splitResult3[x*ppsize+y]));
									}
									for(int z=0; z<ppsize; z++)
									{
										theBigAnsKey += splitResult4[x*ppsize+z]+"-";
										theBigSet.insert(stoi(splitResult4[x*ppsize+z]));
									}
								}
								theBigAnsKey.pop_back();

								pthread_mutex_lock(&mutex1);
								bool notFound = (squareMap.find(theBigAnsKey)==squareMap.end())?true:false;
								pthread_mutex_unlock(&mutex1);

								if(notFound)
								{
									pthread_mutex_lock(&mutex3);
									justCountNewPinpin++;
									pthread_mutex_unlock(&mutex3);

									sinSquare theSQ = {theBigAnsKey,theBigSet};
									vector<sinSquare> tempSinVec;
									tempSinVec.push_back(theSQ);
									//top color
									string topColor = startInfo.topColor+"-";
									topColor += squareMap[no2].topColor;

									//right color
									string rightColor = squareMap[no2].rightColor+"-";
									rightColor += squareMap[no4].rightColor;

									//bottom color
									string bottomColor = squareMap[no3].bottomColor+"-";
									bottomColor += squareMap[no4].bottomColor;

									//left color
									string leftColor = startInfo.leftColor+"-";
									leftColor += squareMap[no3].leftColor;
									

									squareInfo tempInfo = {theBigAnsKey,0,false,false,false,false,false,ppsize*2,topColor,rightColor,bottomColor,leftColor,1,false};
									
									pthread_mutex_lock(&mutex1);
									squareVector.push_back(theBigAnsKey);
									squareSet.insert(theBigAnsKey);
									squareSizeMap[ppsize*2].push_back(theBigAnsKey);
									squareMap[theBigAnsKey] = tempInfo;
									
									if(bigTiles.count(topColor) == 0)
										bigTiles[topColor] = tempSinVec;
									else bigTiles[topColor].push_back(theSQ);

									
									if(bigTiles.count(rightColor) == 0)
										bigTiles[rightColor] = tempSinVec;
									else bigTiles[rightColor].push_back(theSQ);

									
									if(bigTiles.count(bottomColor) == 0)
										bigTiles[bottomColor] = tempSinVec;
									else bigTiles[bottomColor].push_back(theSQ);

									
									if(bigTiles.count(leftColor) == 0)
										bigTiles[leftColor] = tempSinVec;
									else bigTiles[leftColor].push_back(theSQ);
									pthread_mutex_unlock(&mutex1);
								}

								if(ppsize*2 == size)
								{
									vector<string> splitResult;
									superSplit(theBigAnsKey,splitResult,"-");

									pthread_mutex_lock(&mutex3);
									for(int a=0; a<splitResult.size(); a++)
									{
										answer[a]=stoi(splitResult[a]);
									}
									pthread_mutex_unlock(&mutex3);
								}

							}

								
						}
						
					}

				}
			}
			
			
	}
	
	return NULL;


}

