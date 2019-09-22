#include <iostream>
#include <stdlib.h>
using namespace std;
struct tile
{
	int top;
	int right;
	int bottom;
	int left;
};
int main()
{
	long seed = time(NULL);
	int size = 9;
	int totalColor = 10;
	tile data[size][size];
	for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		{
			tile T;
			srand(seed);
			T.right = rand()%totalColor;
			T.bottom = rand()%totalColor;
			seed = rand();
			
			if(j==0)
			{
				srand(seed);
				T.left = rand()%totalColor;
				seed = rand();
			}
			else
			{
				T.left = data[i][j-1].right;
			}
			if(i==0)
			{
				srand(seed);
				T.top = rand()%totalColor;
				seed = rand();
			}
			else
			{
				T.top = data[i-1][j].bottom;
			}

			data[i][j]=T;

		}
	int times = 50;
	while(times--)
	{
		srand(seed);
		int i1 = rand()%size;
		int j1 = rand()%size;
		int i2 = rand()%size;
		int j2 = rand()%size;
		seed = rand();
		tile temp = data[i1][j1];
		data[i1][j1] = data[i2][j2];
		data[i2][j2] = temp;
	}	
	
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++)
		{
			tile T = data[i][j];
			cout<<T.top<<"-"<<T.right<<"-"<<T.bottom<<"-"<<T.left<<" ";
		}
		cout<<endl;
	}
	return 0;
}