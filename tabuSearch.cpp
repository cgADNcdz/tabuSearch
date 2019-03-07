#include<iostream>
#include<stdlib.h>
#include<cstring>  
#include<fstream>  
#include<cmath>
#include<ctime>
#include<string>

 
#define INT_MAX  2147483647  //INT_MAX��Ϊ�������(2147483647)
#define PI 3.1415926
using namespace std;


//��ȡ�ļ��õ��������(ȫ���ӡ��Գ�)���������ֶ�Ӧ�ı��(ע�⣺��Ϊ�����ڲ������ڴ棬��distMatrix����ָ�룬pointName��˫ָ��)
bool readFile(const char* filename, int filetype, double*** distMatrix, string** pointName, int* nPoint);

//���������㷨
bool tabuSearch(double** distMatrix, string* pointName, int** bestPath,
	int nPoints, int tabuLength, int nNeightbors, int nIterations, string startEndPoint,
	int* bestIteration, double* bestDist);


int main(int argc, char** argv)
{
	//����ȫ�ֱ���
	int nPoints;
	int nIterations;     //��������
	int nNeighbors;      //ÿ�ε���ѡȡ���ھ�����
	int tabuLength;      //���ɳ���
	string startEndPoint;//��ʼ��

	int* bestPath;       //�ñ�ű�ʾ�����·��
	int bestIteration;   //���·����Ӧ�ĵ�������
	double bestDist;     //���·����Ӧ���ܳ���

	double** distMatrix; //�������
	string* pointName;   //��Ÿ��ص������

	const char* filename;     //Ҫ��ȡ���ļ���
	int filetype;        //Ҫ��ȡ���ļ�����


						 //��ʼ��ȫ�ֱ���
	filename = (argc >= 2) ?(argv[1]) :"test.txt";  //�����в������ļ���
	filetype = (argc >= 3) ? atoi(argv[2]) : 1;                    //�����в���:�ļ�����(Ĭ��1)
	nIterations = (argc >= 4) ? atoi(argv[3]) : 100;             //�����в���:��������
	nNeighbors = (argc >= 5) ? atoi(argv[4]) : 4;                //�����в���:ѡȡ���ھ���Ŀ
	tabuLength = (argc >= 6) ? atoi(argv[5]) : 10;               //�����в���:���ɱ���!!
	startEndPoint = (argc >= 7) ? (string)(argv[6]) : "A";       //�����в���:��ʼ�������

																 //��ȡ�ļ����õ�������� �� �ص�����֡����
	readFile(filename, filetype, &distMatrix, &pointName, &nPoints);

	//ִ��tabuSearch�㷨,���أ����·�������·���ĵ������������·����Ӧ�ľ���
	clock_t start, finish;
	start = clock();
	tabuSearch(distMatrix, pointName, &bestPath,
		nPoints, tabuLength, nNeighbors, nIterations, startEndPoint,
		&bestIteration, &bestDist);
	finish = clock();


	//������
	cout << "best path:" << endl;                     //������·��
	for (int i = 0; i < nPoints; i++)
		cout<<(string)pointName[bestPath[i]] << "->";
	cout << pointName[bestPath[0]] << endl;
	cout << "best iteration:" << bestIteration << endl;//���·����Ӧ�ĵ�������
	cout << "bestDist:" << bestDist << endl;           //���·����Ӧ�ľ���
	cout << "run time:" << (double)(finish - start) / CLOCKS_PER_SEC << endl;

	//�ͷ��ڴ�
	for (int i = 0; i<nPoints; i++)
		delete[] distMatrix[i];
	delete[] bestPath;
	delete[] pointName;

	cout << "completed!" << endl;
	return 0;
}


//��ȡ�ļ�,�õ�������󡢸������ֶ�Ӧ�ı��
bool readFile(const char* filename, int filetype, double*** distMatrix, string** pointName, int* nPoint)
{
	//���ļ�����fstream��gcc4.9.2�±��벻ͨ������δ�ҵ�ԭ��ֻ��ʹ�������ض��� 
	/*fstream fin(filename);
	
	if (!fin.is_open())
	{
		cout << "can not open the file" << filename << endl;
		return 1;//�����������1
	}*/
	
	freopen(filename,"r",stdin);                        //�����ض��� 

	//filetype�����ļ���ŵ���ʲô���ݣ����� or λ�ã�
	//�ļ���ŵ��λ�ã��ļ���ʽ����һ��д���������¸���д��ľ��ȡ�ά�ȡ��ص�����
	if (filetype == 0)
	{
		int nPoints;   //�ص�����
		cin >> nPoints;

		(*pointName) = new string[nPoints];        //����ռ��Ÿ��������
		(*distMatrix) = new double*[nPoints];      //����������Ĵ�ſռ�
		for (int i = 0; i<nPoints; i++)
		{
			(*distMatrix)[i] = new double[nPoints];
		}

		//������벢����������
		double* longitude = new double[nPoints];
		double* lantitude = new double[nPoints];
		for (int i = 0; i < nPoints; i++)
		{
			cin >> longitude[i] >> lantitude[i] >> (*pointName)[i];
		}
		for (int i = 0; i<nPoints - 1; i++)
		{
			//distMatrix[i][i]=0; //û�б�Ҫ��Ĭ�Ͼ���0
			for (int j = i + 1; j<nPoints - 1; j++)   //����(i,j)����ľ���
			{
				double deltaLon = longitude[i] * PI / 180.0 - longitude[j] * PI / 180.0;  //���ȱ�ʾ�ľ��Ȳ�
				double arcLanI = lantitude[i] * PI / 180.0;                         //i��ά�ȵĻ��ȱ�ʾ
				double arcLanJ = lantitude[j] * PI / 180.0;
				double deltaLan = arcLanI - arcLanJ;                              //���ȱ�ʾ��ά�Ȳ�

				double distIJ = 6378137.0 * 2 * asin(sqrt(pow(sin(deltaLan / 2), 2) + cos(arcLanI) *cos(arcLanJ)*pow(sin(deltaLon / 2), 2)));  //���빫ʽ
				(*distMatrix[i][j]) = distIJ;
				(*distMatrix[j][i]) = distIJ;
			}
		}
		// distMatrix[nPoints-1][nPoints-1]=0;  //�˾�û�б�Ҫ��Ĭ�ϵľ���0
		*nPoint = nPoints;                       //����Ҫ�����ܵĵ���Ŀ
		delete[] longitude;
		delete[] lantitude;
	}

	//�ļ�ֱ�Ӵ�ŵ����֮��ľ���,�ļ���ʽ:��һ��д����,����n��д�ص����֣�����������д����������������
	if (filetype == 1)
	{
		int nPoints;
		cin >> nPoints;
		(*pointName) = new string[nPoints];        //����ռ��ŵ���
		(*distMatrix) = new double*[nPoints];      //����ռ��ž���
		for (int i = 0; i<nPoints; i++)
		{
			(*distMatrix)[i] = new double[nPoints];
		}

		for (int i = 0; i<nPoints; i++)           //��ȡ�ļ�����¼����
		{
			cin >> (*pointName)[i];
		}
		string pointName1, pointName2;
		double dist;
		while (cin >> pointName1 && cin >> pointName2 && cin >> dist)  //��ȡ����ľ��룬����������
		{
			int I, J;
			for (int i = 0; i<nPoints; i++)   //�ҵ��õ����ֵı��
			{
				if ((*pointName)[i] == pointName1)
					I = i;
				if ((*pointName)[i] == pointName2)
					J = i;
			}
			(*distMatrix)[I][J] = dist;
			(*distMatrix)[J][I] = dist;
		}
		*nPoint = nPoints;              //���ص����Ŀ

	}
	//fin.close();
	return 0;
}


//���������㷨
//readme1.�˴���Ŀ��ֵΪ���ɶ���(����·������ڽ��ɱ���)�������ڣ�ѡ�н��ɺ����һ�����Ŀ��ֵ����ǰ���κ�һ����Ѻ�ѡ�⣬���������
//readme2.���ɶ����ѡȡ�����Լ��ھӸ�����ͳ�ƾ����˺��Ѵ���ȫ�������������       ���ϣ��˴�û�п�����������
//readme3.�˴������ɱ�ĳ�����Ϊ���ɳ���,ֻ�е����ɱ�����֮����ɶ���ı��������Ż��������μ�С
bool tabuSearch(double** distMatrix, string* pointName, int** bestPath,
	int nPoints, int tabuLength, int nNeightbors, int nIterations, string startEndPoint,
	int* bestIteration, double* bestDist)
{
	//һЩ��ʼ������.....
	int startNumber;
	for (int i = 0; i<nPoints; i++)               //�ҵ����ı��
		if (pointName[i] == startEndPoint)
		{
			startNumber = i;
			break;
		}
	int** tabuTable = new int*[tabuLength];   //Ϊ���ɱ����ռ�(���ɱ��д�Ž��ɳ��ȸ��ֲ�����·��)
	for (int i = 0; i<tabuLength; i++)
		tabuTable[i] = new int[nPoints];
	(*bestPath) = new int[nPoints];           //Ϊ����·������ռ�
	*bestIteration = -1;
	srand((unsigned int)time(0));       //�������������
	int nPath = 0;                        //��¼���ɱ���ʵ�ʱ����·������

										  /**** 1.�������һ����ʼ��*****/
	int* currPath = new int[nPoints];      //Ϊ��ǰ·������ռ�
	double currDist = 0;

	currPath[0] = startNumber;             //ָ����ʼ��
	for (int i = 1; i<nPoints; i++)           //��������ĳ�ʼ·��
	{
		currPath[i] = rand() % nPoints;
		for (int j = 0; j<i; j++) //�ų��ظ�,���ǣ�����
		{
			if (currPath[j] == currPath[i])
			{
				i--;
				break;
			}
		}
	}

	/*cout<<"��ʼ·����";
	for(int i=0;i<nPoints;i++)
	cout<<pointName[currPath[i]]<<"->";
	cout<<pointName[currPath[0]]<<endl;*/

	for (int i = 0; i<nPoints - 1; i++)         //���㵱ǰ·���ĳ���
	{
		currDist += distMatrix[currPath[i]][currPath[i + 1]];
	}
	currDist += distMatrix[currPath[nPoints - 1]][currPath[0]];

	memcpy(*bestPath, currPath, sizeof(int)*nPoints);  //��Ϊ�ǵ�ǰ����
	*bestDist = currDist;


	//��ʼ����
	int* localBestPath = new int[nPoints];          //��ŵ���ʱ�ľֲ�����·��
	double localBestDist = INT_MAX;                //��žֲ�����·������
	int* tempNeighborPath = new int[nPoints];           //��ŵ�ǰ�ھӵ�·��
	double tempNeightborDist = 0;                   //��ŵ�ǰ�ھӵ�·������
	for (int nIter = 0; nIter<nIterations; nIter++)
	{
		/**** 2.ÿ�ε�����nNeightbors�����ڽ��ɱ��е��ھ���Ϊ��ѡ����ѡ����ѡ�����оֲ����Ž�******/
		int nn = 0;  //nn��¼�Ѿ��ҵ����ھ����������ظ������ﻹδ�Ľ���
		localBestDist = INT_MAX; //ÿһ�ε�����Ӧ������localBest
		while (nn<nNeightbors)
		{
			// 2.1 ���������ǰ��(currPath)���������õ�һ���ھ�
			//����ֻд����һ�ֽ������������������������������������һ���м�����
			int subscript1 = rand() % nPoints;
			while (subscript1 == 0)                 //���ѡ�񽻻���1(�������Ҫ�����ĵص���±�)���Ҳ��������
				subscript1 = rand() % nPoints;
			int subscript2 = rand() % nPoints;
			while (subscript2 == 0 || subscript2 == subscript1)   //���ѡ�񽻻���2(�±�)���Ҳ����������߽�����1
				subscript2 = rand() % nPoints;
			memcpy(tempNeighborPath, currPath, sizeof(int)*nPoints);  //�Ƚ���ǰ�⸳ֵ���ھ�
			int temp = tempNeighborPath[subscript1];                  //������������λ��
			tempNeighborPath[subscript1] = tempNeighborPath[subscript2];
			tempNeighborPath[subscript2] = temp;
			tempNeightborDist = 0;    //����tempNeighborDist

									  // 2.2�����ǰ�ھӲ��ڽ��ɱ�����·�����ȸ��̣�����¾ֲ�����
			bool inTabuTable = 1;
			for (int i = 0; i<tabuLength; i++)    //�жϵ�ǰ�ھ��Ƿ��ڽ��ɱ���
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

			for (int i = 0; i<nPoints - 1; i++)     //���㵱ǰ�ھӵ�·������
				tempNeightborDist += distMatrix[tempNeighborPath[i]][tempNeighborPath[i + 1]];
			tempNeightborDist += distMatrix[tempNeighborPath[nPoints - 1]][tempNeighborPath[0]];
			if (!inTabuTable)                //������ڽ��ɱ�
			{
				nn++;     //���ڽ��ɱ�ʱ�ż����ھ���!!!!
				if (tempNeightborDist<localBestDist)  //���¾ֲ�����
				{
					memcpy(localBestPath, tempNeighborPath, sizeof(int)*nPoints);
					localBestDist = tempNeightborDist;
				}
			}
		}

		/*** 3.���£���ǰ�����ľֲ����Ž������ɱ�,������һ�ε����ĳ�ʼ��,����֮ǰ������������Ž⣨ȫ�����ţ�***/
		// 3.1����ֲ�������ȫ�����ţ������ȫ������
		if (localBestDist<*bestDist)
		{
			memcpy(*bestPath, localBestPath, sizeof(int)*nPoints);
			*bestDist = localBestDist;
			*bestIteration = nIter;
		}

		//3.2���ôε����ľֲ����ż�����ɱ�
		//ע�����ڽ��ɳ���������˵����һ���ǽ��ɱ�ĳ��ȣ���һ����ָ���ɶ����ֹ��ѡȡ�ĵ�������
		//����ȡ�ĵ�һ��˵������Ϊ���򵥣�ʵ���ϵ����ɱ�����֮���һ��˵���͵�ͬ�ڵڶ���˵����
		if (nPath<tabuLength) //��δ�ﵽ���ɱ��ȣ�ֱ�Ӽ�����ɱ�
		{
			for (int i = 0; i<nPoints; i++)
				tabuTable[nPath][i] = localBestPath[i];
			nPath++;
		}
		else                 //���ɱ�����
		{
			for (int i = 0; i<tabuLength - 1; i++) //ǰ��һ����Ԫ
				for (int j = 0; j<nPoints; j++)
					tabuTable[i][j] = tabuTable[i + 1][j];
			for (int i = 0; i<nPoints; i++)     //�ֲ����Ŵ���ڽ��ɱ�ĩβ
				tabuTable[tabuLength - 1][i] = localBestPath[i];
		}

		//3.3�����´ε����ĳ�ʼ��
		memcpy(currPath, localBestPath, sizeof(int)*nPoints);
	}

	//�ͷ��ڴ�
	for (int i = 0; i<tabuLength; i++)
		delete[] tabuTable[i];
	delete[] currPath;
	delete[] tempNeighborPath;
	delete[] localBestPath;
	return 0;
}
