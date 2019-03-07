#include<iostream>
#include<stdlib.h>
#include<cstring>  
#include<fstream>  
#include<cmath>
#include<ctime>
#include<string>

 
#define INT_MAX  2147483647  //INT_MAX作为最大整数(2147483647)
#define PI 3.1415926
using namespace std;


//读取文件得到距离矩阵(全连接、对称)、各点名字对应的编号(注意：因为是在内部分配内存，故distMatrix用三指针，pointName用双指针)
bool readFile(const char* filename, int filetype, double*** distMatrix, string** pointName, int* nPoint);

//禁忌搜索算法
bool tabuSearch(double** distMatrix, string* pointName, int** bestPath,
	int nPoints, int tabuLength, int nNeightbors, int nIterations, string startEndPoint,
	int* bestIteration, double* bestDist);


int main(int argc, char** argv)
{
	//声明全局变量
	int nPoints;
	int nIterations;     //迭代次数
	int nNeighbors;      //每次迭代选取的邻居数量
	int tabuLength;      //禁忌长度
	string startEndPoint;//起始点

	int* bestPath;       //用编号表示的最佳路径
	int bestIteration;   //最佳路径对应的迭代次数
	double bestDist;     //最佳路径对应的总长度

	double** distMatrix; //距离矩阵
	string* pointName;   //存放各地点的名字

	const char* filename;     //要读取的文件名
	int filetype;        //要读取的文件类型


						 //初始化全局变量
	filename = (argc >= 2) ?(argv[1]) :"test.txt";  //命令行参数：文件名
	filetype = (argc >= 3) ? atoi(argv[2]) : 1;                    //命令行参数:文件类型(默认1)
	nIterations = (argc >= 4) ? atoi(argv[3]) : 100;             //命令行参数:迭代次数
	nNeighbors = (argc >= 5) ? atoi(argv[4]) : 4;                //命令行参数:选取的邻居数目
	tabuLength = (argc >= 6) ? atoi(argv[5]) : 10;               //命令行参数:禁忌表长度!!
	startEndPoint = (argc >= 7) ? (string)(argv[6]) : "A";       //命令行参数:起始点的名字

																 //读取文件，得到距离矩阵 和 地点的名字、序号
	readFile(filename, filetype, &distMatrix, &pointName, &nPoints);

	//执行tabuSearch算法,返回：最佳路径、最佳路径的迭代次数、最佳路径对应的距离
	clock_t start, finish;
	start = clock();
	tabuSearch(distMatrix, pointName, &bestPath,
		nPoints, tabuLength, nNeighbors, nIterations, startEndPoint,
		&bestIteration, &bestDist);
	finish = clock();


	//输出结果
	cout << "best path:" << endl;                     //输出最佳路径
	for (int i = 0; i < nPoints; i++)
		cout<<(string)pointName[bestPath[i]] << "->";
	cout << pointName[bestPath[0]] << endl;
	cout << "best iteration:" << bestIteration << endl;//最佳路径对应的迭代次数
	cout << "bestDist:" << bestDist << endl;           //最佳路径对应的距离
	cout << "run time:" << (double)(finish - start) / CLOCKS_PER_SEC << endl;

	//释放内存
	for (int i = 0; i<nPoints; i++)
		delete[] distMatrix[i];
	delete[] bestPath;
	delete[] pointName;

	cout << "completed!" << endl;
	return 0;
}


//读取文件,得到距离矩阵、各点名字对应的编号
bool readFile(const char* filename, int filetype, double*** distMatrix, string** pointName, int* nPoint)
{
	//打开文件流，fstream在gcc4.9.2下编译不通过，尚未找到原因，只好使用输入重定向 
	/*fstream fin(filename);
	
	if (!fin.is_open())
	{
		cout << "can not open the file" << filename << endl;
		return 1;//错误情况返回1
	}*/
	
	freopen(filename,"r",stdin);                        //输入重定向 

	//filetype代表文件存放的是什么数据（距离 or 位置）
	//文件存放点的位置，文件格式：第一行写点数，以下各行写点的精度、维度、地点名字
	if (filetype == 0)
	{
		int nPoints;   //地点数量
		cin >> nPoints;

		(*pointName) = new string[nPoints];        //申请空间存放各点的名字
		(*distMatrix) = new double*[nPoints];      //申请距离矩阵的存放空间
		for (int i = 0; i<nPoints; i++)
		{
			(*distMatrix)[i] = new double[nPoints];
		}

		//计算距离并存入距离矩阵
		double* longitude = new double[nPoints];
		double* lantitude = new double[nPoints];
		for (int i = 0; i < nPoints; i++)
		{
			cin >> longitude[i] >> lantitude[i] >> (*pointName)[i];
		}
		for (int i = 0; i<nPoints - 1; i++)
		{
			//distMatrix[i][i]=0; //没有必要，默认就是0
			for (int j = i + 1; j<nPoints - 1; j++)   //计算(i,j)两点的距离
			{
				double deltaLon = longitude[i] * PI / 180.0 - longitude[j] * PI / 180.0;  //弧度表示的经度差
				double arcLanI = lantitude[i] * PI / 180.0;                         //i点维度的弧度表示
				double arcLanJ = lantitude[j] * PI / 180.0;
				double deltaLan = arcLanI - arcLanJ;                              //弧度表示的维度差

				double distIJ = 6378137.0 * 2 * asin(sqrt(pow(sin(deltaLan / 2), 2) + cos(arcLanI) *cos(arcLanJ)*pow(sin(deltaLon / 2), 2)));  //距离公式
				(*distMatrix[i][j]) = distIJ;
				(*distMatrix[j][i]) = distIJ;
			}
		}
		// distMatrix[nPoints-1][nPoints-1]=0;  //此句没有必要，默认的就是0
		*nPoint = nPoints;                       //最终要返回总的点数目
		delete[] longitude;
		delete[] lantitude;
	}

	//文件直接存放点与点之间的距离,文件格式:第一行写点数,下面n行写地点名字，接下来各行写点名、点名、距离
	if (filetype == 1)
	{
		int nPoints;
		cin >> nPoints;
		(*pointName) = new string[nPoints];        //申请空间存放点名
		(*distMatrix) = new double*[nPoints];      //申请空间存放距离
		for (int i = 0; i<nPoints; i++)
		{
			(*distMatrix)[i] = new double[nPoints];
		}

		for (int i = 0; i<nPoints; i++)           //读取文件，记录点名
		{
			cin >> (*pointName)[i];
		}
		string pointName1, pointName2;
		double dist;
		while (cin >> pointName1 && cin >> pointName2 && cin >> dist)  //读取后面的距离，加入距离矩阵
		{
			int I, J;
			for (int i = 0; i<nPoints; i++)   //找到该点名字的编号
			{
				if ((*pointName)[i] == pointName1)
					I = i;
				if ((*pointName)[i] == pointName2)
					J = i;
			}
			(*distMatrix)[I][J] = dist;
			(*distMatrix)[J][I] = dist;
		}
		*nPoint = nPoints;              //返回点的数目

	}
	//fin.close();
	return 0;
}


//禁忌搜索算法
//readme1.此处以目标值为禁忌对象(即把路径存放在禁忌表中)，不存在：选中禁忌后出现一个解的目标值好于前面任何一个最佳候选解，的特赦情况
//readme2.禁忌对象的选取方法以及邻居个数的统计决定了很难存在全部被禁的情况。       综上，此处没有考虑特赦的情况
//readme3.此处将禁忌表的长度作为禁忌长度,只有当禁忌表满了之后禁忌对象的被禁次数才会整体依次减小
bool tabuSearch(double** distMatrix, string* pointName, int** bestPath,
	int nPoints, int tabuLength, int nNeightbors, int nIterations, string startEndPoint,
	int* bestIteration, double* bestDist)
{
	//一些初始化工作.....
	int startNumber;
	for (int i = 0; i<nPoints; i++)               //找到起点的编号
		if (pointName[i] == startEndPoint)
		{
			startNumber = i;
			break;
		}
	int** tabuTable = new int*[tabuLength];   //为禁忌表分配空间(禁忌表中存放禁忌长度个局部最优路径)
	for (int i = 0; i<tabuLength; i++)
		tabuTable[i] = new int[nPoints];
	(*bestPath) = new int[nPoints];           //为最优路径分配空间
	*bestIteration = -1;
	srand((unsigned int)time(0));       //设置随机数种子
	int nPath = 0;                        //记录禁忌表中实际保存的路径数量

										  /**** 1.随机生成一个初始解*****/
	int* currPath = new int[nPoints];      //为当前路径分配空间
	double currDist = 0;

	currPath[0] = startNumber;             //指定开始点
	for (int i = 1; i<nPoints; i++)           //生成随机的初始路径
	{
		currPath[i] = rand() % nPoints;
		for (int j = 0; j<i; j++) //排除重复,机智！！！
		{
			if (currPath[j] == currPath[i])
			{
				i--;
				break;
			}
		}
	}

	/*cout<<"初始路径：";
	for(int i=0;i<nPoints;i++)
	cout<<pointName[currPath[i]]<<"->";
	cout<<pointName[currPath[0]]<<endl;*/

	for (int i = 0; i<nPoints - 1; i++)         //计算当前路径的长度
	{
		currDist += distMatrix[currPath[i]][currPath[i + 1]];
	}
	currDist += distMatrix[currPath[nPoints - 1]][currPath[0]];

	memcpy(*bestPath, currPath, sizeof(int)*nPoints);  //认为是当前最优
	*bestDist = currDist;


	//开始迭代
	int* localBestPath = new int[nPoints];          //存放迭代时的局部最优路径
	double localBestDist = INT_MAX;                //存放局部最优路径长度
	int* tempNeighborPath = new int[nPoints];           //存放当前邻居的路径
	double tempNeightborDist = 0;                   //存放当前邻居的路径长度
	for (int nIter = 0; nIter<nIterations; nIter++)
	{
		/**** 2.每次迭代找nNeightbors个不在禁忌表中的邻居作为候选集，选出候选集合中局部最优解******/
		int nn = 0;  //nn记录已经找到的邻居数（可能重复，这里还未改进）
		localBestDist = INT_MAX; //每一次迭代都应该重置localBest
		while (nn<nNeightbors)
		{
			// 2.1 随机交换当前解(currPath)的两个结点得到一个邻居
			//这里只写了这一种交换方法，还可以有其他方法，如随机交换一段中间距离等
			int subscript1 = rand() % nPoints;
			while (subscript1 == 0)                 //随机选择交换点1(随机数是要交换的地点的下标)，且不能是起点
				subscript1 = rand() % nPoints;
			int subscript2 = rand() % nPoints;
			while (subscript2 == 0 || subscript2 == subscript1)   //随机选择交换点2(下标)，且不能是起点或者交换点1
				subscript2 = rand() % nPoints;
			memcpy(tempNeighborPath, currPath, sizeof(int)*nPoints);  //先将当前解赋值给邻居
			int temp = tempNeighborPath[subscript1];                  //交换两个结点的位置
			tempNeighborPath[subscript1] = tempNeighborPath[subscript2];
			tempNeighborPath[subscript2] = temp;
			tempNeightborDist = 0;    //重置tempNeighborDist

									  // 2.2如果当前邻居不在禁忌表中且路径长度更短，则更新局部最优
			bool inTabuTable = 1;
			for (int i = 0; i<tabuLength; i++)    //判断当前邻居是否在禁忌表中
			{
				bool flag = 0;
				for (int j = 0; j<nPoints; j++)
				{
					if (tempNeighborPath[j] != tabuTable[i][j])
					{
						flag = 1;
						inTabuTable = 0;
						break;
					}
				}
				if (flag == 1)
					break;
			}

			for (int i = 0; i<nPoints - 1; i++)     //计算当前邻居的路径长度
				tempNeightborDist += distMatrix[tempNeighborPath[i]][tempNeighborPath[i + 1]];
			tempNeightborDist += distMatrix[tempNeighborPath[nPoints - 1]][tempNeighborPath[0]];
			if (!inTabuTable)                //如果不在禁忌表
			{
				nn++;     //不在禁忌表时才计算邻居数!!!!
				if (tempNeightborDist<localBestDist)  //更新局部最优
				{
					memcpy(localBestPath, tempNeighborPath, sizeof(int)*nPoints);
					localBestDist = tempNeightborDist;
				}
			}
		}

		/*** 3.更新：当前迭代的局部最优解加入禁忌表,更新下一次迭代的初始解,更新之前所有情况的最优解（全局最优）***/
		// 3.1如果局部最优是全局最优，则更新全局最优
		if (localBestDist<*bestDist)
		{
			memcpy(*bestPath, localBestPath, sizeof(int)*nPoints);
			*bestDist = localBestDist;
			*bestIteration = nIter;
		}

		//3.2将该次迭代的局部最优加入禁忌表
		//注：关于禁忌长度有两种说法，一种是禁忌表的长度；另一种是指禁忌对象禁止被选取的迭代次数
		//这里取的第一种说法，因为更简单（实际上当禁忌表满了之后第一种说法就等同于第二种说法）
		if (nPath<tabuLength) //还未达到禁忌表长度，直接加入禁忌表
		{
			for (int i = 0; i<nPoints; i++)
				tabuTable[nPath][i] = localBestPath[i];
			nPath++;
		}
		else                 //禁忌表已满
		{
			for (int i = 0; i<tabuLength - 1; i++) //前移一个单元
				for (int j = 0; j<nPoints; j++)
					tabuTable[i][j] = tabuTable[i + 1][j];
			for (int i = 0; i<nPoints; i++)     //局部最优存放在禁忌表末尾
				tabuTable[tabuLength - 1][i] = localBestPath[i];
		}

		//3.3跟新下次迭代的初始解
		memcpy(currPath, localBestPath, sizeof(int)*nPoints);
	}

	//释放内存
	for (int i = 0; i<tabuLength; i++)
		delete[] tabuTable[i];
	delete[] currPath;
	delete[] tempNeighborPath;
	delete[] localBestPath;
	return 0;
}
