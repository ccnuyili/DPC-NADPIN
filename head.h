#pragma warning(disable : 4786)
#include<iostream>
#include<sstream>
#include<string>
#include<cstring>
#include<fstream>
#include<vector>
#include<time.h>
#include<math.h>
#include"stdio.h"

using namespace std;

typedef vector<int> VecInt;
typedef vector<double> VecDou;
typedef vector<string> VecStr;

vector<VecStr> Results;  //�����ʸ�����������
vector<int>  frequency;  //�����ʸ�������36�������г��ֵ�Ƶ��

//���ܣ���������������ֵ
void Swap(string &str1, string &str2)
{
	string temp;
	temp = str1;
	str1 = str2; 
	str2 = temp;
}

//���ܣ��������򷽷��еĻ��ַ���
int Partition(vector<string> &nodemaptemp, int p, int r)
{
	string temp = nodemaptemp[r];
	int i = p - 1, j;
	for(j=p; j<r; j++)
	{
		if(nodemaptemp[j].compare(temp)<=0)
		{
			i++;
			Swap(nodemaptemp[i], nodemaptemp[j]);
		}
	}
	Swap(nodemaptemp[i+1], nodemaptemp[r]);
	return i+1;
}

//���ܣ��������ʽڵ���������ĸ˳���������
void Quicksort(vector<string> &nodemaptemp, int p, int r)
{
	if(p<r)
	{
		int q = Partition(nodemaptemp, p, r);
		Quicksort(nodemaptemp, p, q-1);
		Quicksort(nodemaptemp, q+1, r);
	}
}

//���ܣ���ȡ��̬�����ļ���ͳ�ƽڵ���Ŀ������Ŀ������ڵ�����
void GeneratePNodeMap(vector<string> &pnodemap,string filename)
{
    ifstream in;
	in.open(filename.c_str());
	string node1, node2; 
	vector<string> nodemaptemp;
	while(in >> node1 >> node2)
	{
		nodemaptemp.push_back(node1);
	    nodemaptemp.push_back(node2);
	}
	in.close();

	Quicksort(nodemaptemp, 0, nodemaptemp.size()-1);

	//ȥ���ظ��ڵ�
	int low = 0, i;
	for(i = 0; i<nodemaptemp.size();)
	{
		low = i;
		pnodemap.push_back(nodemaptemp[low]);
        do
		{
			i += 1;
		}while(i<nodemaptemp.size() && nodemaptemp[i].compare(nodemaptemp[low]) ==0);
	}
}

//���ܣ���ʼ���໥��������
void InitialNetwork(vector<VecInt> &pnetwork, int rownum, int colnum)
{
	for(int i=0; i<rownum; ++i)
	{
		vector<int> temp(colnum, 0);
		pnetwork.push_back(temp);
	}
}

//���ܣ����������ʽڵ㣬�ҵ����ڵ����ʽڵ������ж�Ӧ�ı��
int GetIndex(vector<string> &pnodemap, string node)
{
	int pos = -1, low = 0, mid = 0, high = pnodemap.size()-1;
	do
	{
		mid = (low + high) /2;
		if(node == pnodemap[mid])
		{
			pos = mid;
			break;
		}
		else if(node.compare(pnodemap[mid]) < 0)
			high = mid - 1;
		else 
			low = mid + 1;
	}while(low <= high);
	return pos;
}

//���ܣ���ȡ�໥�����ļ��������໥����(�Գ�)����
int ReadData(string datafile, vector<string> &pnodemap, int nodenum, vector<VecInt> &pnetwork)
{
	ifstream in; 
	in.open(datafile.c_str());
	string node1, node2;
	int nodeindex1, nodeindex2;
	int edgenum = 0 ;
	while(in >> node1 >> node2)
	{
		edgenum++;
		nodeindex1 = GetIndex(pnodemap, node1);
		nodeindex2 = GetIndex(pnodemap, node2);
		pnetwork[nodeindex1][nodeindex2] += 1;  //�Գ�����
		pnetwork[nodeindex2][nodeindex1] += 1;
	}
	in.close();
	return edgenum;
}

float MiddleCC(vector<float> clusterco)
{
	int num = clusterco.size();
	int i; 
	float max = -1.0; 
	float min = 9999.0;
	for(i=0; i<num; i++)
	{
		if(clusterco[i] > max)
			max = clusterco[i];
		if(clusterco[i] < min)
			min = clusterco[i];
	}
	float middlecc = (max + min) / 4.0;
	return middlecc;
}

//���ܣ�ͳ��ÿ�������ʽڵ�ľۼ�ϵ��
float StatisticCC(vector<VecInt> &pnetwork, vector<float> &coefficient, int nodenum)
{
	int i, j, m, n, index1, index2;
	int count;  //�ھӽڵ����
	vector<int> neighbor;  //�ھӽڵ㼯��
	float cc; //�����ʽڵ�ۼ�ϵ��
	int realnum;  //�ھ�֮��ʵ���໥���ñ���
	int maxnum;  //�ھ�֮���������໥������Ŀ
	float avecc;  //����ƽ���ݼ�ϵ��

	for(i=0; i<nodenum; ++i)
	{
		count = 0;  //�ھӽڵ��ʼ��
		for(j=0; j<nodenum; ++j)
		{
			if(pnetwork[i][j] == 1)
			{
				++count;
				neighbor.push_back(j);  //�����ھӽڵ�(�ڵ��±�)
			}
		}
		realnum = 0;
		//ͳ���ھӽڵ�֮����໥���ñ���
		for(m=0; m<neighbor.size()-1; ++m)
		{
			index1 = neighbor[m];
			for(n=m+1; n<neighbor.size(); ++n)
			{
				index2 = neighbor[n];
				if(pnetwork[index1][index2] == 1)
					++realnum;
			}
		}
		if(count == 1)  //û�й����Ľڵ�
			cc = 0;	
		else
		{
			maxnum = count *(count-1);
			cc = 2*realnum / (float)(maxnum);
		}
		coefficient.push_back(cc);
		neighbor.clear();  //������ݣ��������
	}

	float sumcc = 0;
	for(i=0; i<coefficient.size(); ++i)
		sumcc += coefficient[i];
	avecc = sumcc / nodenum;  //�����ƽ���ۼ�ϵ��
	return avecc;
}

//���ܣ����ݽڵ�ϵ���Ĺ����ʼ�˼���
void ConstructCore(vector<VecInt> &pnetwork, vector<VecInt> &initialcore, vector<float> &clusterco, float avecc)
{
	//�ۼ�ϵ������ƽ���ۼ�ϵ���Ľڵ㼰���ھӽڵ�һ����Ϊ��ʼ��
	int i, j; 
	vector<int> core;  //��ʼ���б�����ǵ����ʽڵ���±꣬��������ʡ�洢�ռ�
	for(i=0; i<clusterco.size(); ++i)
	{
		if(clusterco[i] > avecc)  //��ǰ��ľۼ�ϵ������ƽ���ݼ�ϵ��
		{
			core.push_back(i);  //��ǰ�ڵ���Ϊ��ʼ�˵ĳ�Ա
			for(j=0; j<clusterco.size(); ++j)
			{
				if(pnetwork[i][j] == 1)
					core.push_back(j);  //��ǰ�ڵ���ھӽڵ�����ʼ����
			}
			initialcore.push_back(core);  //�����ʼ�˼�����
			core.clear();  //�����ʱ�˳�Ա
		}
	}
}

//���ܣ�����չ->�����ʸ�����
void ExtendCore(vector<VecInt> &pnetwork, vector<VecInt> &initialcore, int nodenum)
{
	//ÿ���������ھ���ͳ���ļ�

	vector<int> candidate;  //��ǰ�غ�ѡ�ھӽڵ㼯��
	vector<int> flag;  //��ǰ���нڵ���
	vector<int> neighbor_in;  //��ѡ�ڵ�����ھӽڵ㼯��
	vector<int> neighbor_out;  //��ѡ�ڵ�����ھӽڵ㼯��

	int i, j, k, r, index, index1, index2;

	int count = 0;  //��¼���ٸ���ʼ�˱���չ

	int initialsize, aftersize, extendlevel;  //��ʼ�˴�С����չ���С����չ����

	//���ζ�ÿ���˽�����չ(�ɽ����ĵ����ʸ�����)
	for(i=0; i<initialcore.size(); ++i)
	{
	    initialsize = 0;
		aftersize = 0;
		extendlevel= 0;
		do	
		{   
			initialsize = initialcore[i].size();
			
			//��Ǽ��ϳ�ʼ��
			for(j=0; j<nodenum; ++j)
				flag.push_back(0);  
			
			//��ǵ�ǰ���е�ÿ���ڵ�
			for(j=0; j<initialcore[i].size(); ++j)
			{
				index = initialcore[i][j];  //initialcore[i][j] ��pnodemap��ÿ�������ʽڵ��Ӧ���±�
				flag[index] = 1;  //�ڵ��ڵ�ǰ���� 
			}
			
		    //Ѱ�ҵ�ǰ�ص��ھӽڵ�(��ѡ�ڵ㼯��)
		    for(j=0; j<initialcore[i].size(); ++j)
			{
		    	index = initialcore[i][j];
		    	for(k=0; k<nodenum; ++k)
				{
			    	if(pnetwork[index][k] == 1 && flag[k] != 1)  //�ǵ�ǰ�ص��ھӽڵ㣬�Ҹýڵ㲻�ڵ�ǰ����
					{
				    	flag[k] = 1;  //�����Ϊ�˷�ֹ�ظ����ھӽڵ�
				    	candidate.push_back(k);  //�����ѡ�ڵ�
					}
				}
			}
            
			float maxmargin = 0;  //���ܶȲ�ֵ����ʼ��
			int thebest = -1;  //��Ѻ�ѡ�ڵ��ʼ��

			if(candidate.size() == 0)  //�����ں�ѡ�ڵ�
				break;
			else  //���ں�ѡ�ڵ�
			{
				flag.clear();  //�������������صĴ���
				//��Ǽ��ϳ�ʼ��
	    	    for(j=0; j<nodenum; ++j)
					flag.push_back(0); 
				
				//��ǵ�ǰ���е�ÿ���ڵ�
				for(j=0; j<initialcore[i].size(); ++j)
				{
					index = initialcore[i][j];  //initialcore[i][j] ��pnodemap��ÿ�������ʽڵ��Ӧ���±�
					flag[index] = 1;  //�ڵ��ڵ�ǰ���� 
				}

		    	//�����ж�ÿ����ѡ�ڵ��Ƿ��ܹ����뵱ǰ��
		    	for(j=0; j<candidate.size(); ++j)
				{ 
					//�ҵ������ھӽڵ㲢��������ھӽڵ����ӽ��ܶ�

			    	index = candidate[j];  //��ѡ�ڵ��Ӧ�ĵ����ʽڵ���±�

			    	for(k=0; k<initialcore[i].size(); ++k)
					{
				    	index1 = initialcore[i][k];
				    	if(pnetwork[index][index1] == 1)
					    	neighbor_in.push_back(index1);  //��ǰ��ѡ�ڵ��ڴ��ڵ��ھӽڵ㼯��	
					}
			    	int sumdegree_in = 0;  //Ҫ��ÿ���ڵ��ʼ��
			    	for(k=0; k<neighbor_in.size(); ++k)
					{
				    	index1 = neighbor_in[k];
				    	for(r=0; r<initialcore[i].size(); ++r)
						{
					    	index2 = initialcore[i][r];
					    	if(pnetwork[index1][index2] == 1)
						    	sumdegree_in += 1;
						}
					}
				    float NIAD = sumdegree_in * 1.0 / (neighbor_in.size() * 1.0);  //�����ھӽڵ�ƽ�����ܶ� 

					//�ҵ������ھӽڵ㲢��������ھӽڵ����ӽ��ܶ�

					float NOAD = 0.0;
			    	for(k=0; k<nodenum; ++k)
					{
				    	if(pnetwork[index][k] == 1 && flag[k] == 0)
				    		neighbor_out.push_back(k);  //��ǰ��ѡ�ڵ��ڴ�����ھӽڵ㼯��
					}
			    	if(neighbor_out.size() != 0)
					{
				    	int sumdegree_out = 0;
				        for(k=0; k<neighbor_out.size(); ++k)
						{
				    		index2 = neighbor_out[k];
				    	    for(r=0; r<nodenum; ++r)
							{
								if(r != index)  //ȥ����ѡ�ڵ�
								{
									if(flag[r] == 0 && pnetwork[index2][r] == 1)
										sumdegree_out += 1;
								}
							}
						}
				        NOAD =  sumdegree_out * 1.0 / (neighbor_out.size() * 1.0);
					}

			    	float margin= NIAD - NOAD;
			    	
			    	if(margin > maxmargin)
					{
				    	maxmargin = margin;
				    	thebest = candidate[j];  //��ѡ�ڵ��е���ѵ�һ��
					}
			    	neighbor_in.clear();
			    	neighbor_out.clear();
				}
				flag.clear();  //�������������صĴ���
			}			
	    	if(thebest != -1)
			{
				initialcore[i].push_back(thebest);
				extendlevel++;
			}
			candidate.clear();
			aftersize = initialcore[i].size();
		
		}while(aftersize > initialsize && extendlevel < 1);
	} 
}

//���ܣ��������ݣ���������
void Swap(int &a, int &b)
{
	int temp;
	temp = a;
	a = b; 
	b =temp;
}

//���ܣ�ð������
void bubblesort(vector<int> &temp, int p, int r)
{
	int i, j; 
	for(i=p; i<=r-1; ++i)
	{
		for(j=i+1; j<=r; ++j)
		{
			if(temp[i] > temp[j])
				Swap(temp[i], temp[j]);
		}
	}
}

//���ܣ��������ģ�������Ĳ���
vector<int> UnionCluster(vector<int> clu_1,  vector<int> clu_2)
{
	int i;
	vector<int> temp;
	for(i=0; i<clu_1.size(); ++i)
		temp.push_back(clu_1[i]);
	for(i=0; i<clu_2.size(); ++i)
		temp.push_back(clu_2[i]);
	bubblesort(temp, 0, temp.size()-1);
	int low; 
	vector<int> unionclu;
	for(i=0; i<temp.size(); )
	{
		low = i;
		unionclu.push_back(temp[i]);
		do
		{
			++i;
		}while(i<temp.size() && temp[i] == temp[low]);  //ȥ���ظ�Ԫ��
	}
	return unionclu;
}

//���ܣ���������ģ���������ص���
double OverlapRate(vector<int> clu_1, vector<int> clu_2)
{
	if(clu_1.size() == 0 || clu_2.size() == 0)  //�������һ��Ϊ�գ�����0
		return 0.0;
	vector<int> temp = UnionCluster(clu_1, clu_2);
	return (clu_1.size() + clu_2.size() - temp.size()) * 1.0 /temp.size();
}

double InitialMatchrate(vector<VecDou> &matchrate, int rowsize, int colsize)
{
	for(int i=0; i<rowsize; i++)
	{
		vector<double> temp(colsize, 0.0);
		matchrate.push_back(temp);
	}
	return 0.0;
}

//���ܣ��ϲ������ص���Ϊ1.0��ģ��
int MergeCluster(vector<VecInt> &initialcore)
{
	//ÿ���������ھ���ͳ���ļ�

	int i, j, k;
	double mr;  //mergerate �ص���
	bool flag = true;
	int one, two;
	vector<VecDou> matchrate;  //�ص��Ⱦ���

	InitialMatchrate(matchrate, initialcore.size(), initialcore.size());

	//��Ч��ʼ��
	for(i=0; i<initialcore.size(); i++)
	{
		for(j=i; j<initialcore.size(); j++)
		{
			mr = OverlapRate(initialcore[i], initialcore[j]);
			matchrate[i][j] = mr;
			matchrate[j][i] = mr;
		}
	}
    
	//�ҳ���ǰ�ص�������ģ��
	do
	{
		flag = false;
		double maxrate = 0.0;
		one = 0;
		two = 0;
	    for(i=0; i<initialcore.size(); i++)
		{
			for(j=i+1; j<initialcore.size(); j++)
			{
				if(matchrate[i][j] > maxrate)
				{
					maxrate = matchrate[i][j];
				    one = i;
				    two = j;
				}
			}
		}

	    if(fabs(maxrate-1.0) <0.000001)
		{
	    	flag = true;
     		vector<int> temp = UnionCluster(initialcore[one], initialcore[two]);
	    	initialcore[one].clear();
	    	initialcore[two].clear();
     		initialcore[one] = temp;
    
			//����matchrate[][]
	    	for(k=0; k<initialcore.size(); k++)
			{
     			matchrate[two][k] = -1.0;
    			matchrate[k][two] = -1.0;
			}
        	for(k=0; k<initialcore.size(); k++)
			{
	     		mr = OverlapRate(initialcore[one], initialcore[k]);
	    		//if(fabs(matchrate[one][k] -(-1.0)) > 0.000001 )
				if(!initialcore[k].empty())
				{
		    		matchrate[one][k] = mr;
	    		    matchrate[k][one] = mr;
				}
			}
		} 
	}while(flag == true);

	//��д���Ľ����
	vector<VecInt> finalclu;
	for(i=0; i<initialcore.size(); ++i)
	{
		if(!initialcore[i].empty())
			finalclu.push_back(initialcore[i]);
	}
	initialcore = finalclu;
	return 0;
	return 1;
}

//Ԥ�⵰���ʸ�����������
void ResultsAnalysis(vector<VecInt> &initialcore)
{
    int i;
	int num = initialcore.size();  //�ھ򸴺���ĸ���
	int number = 0;

	//�������ƽ����С
	float aversize;
	for(i=0; i<num; ++i)
		number += initialcore[i].size();
	aversize = number / num;
	
	//���������
	int maxnum = 0;
	for(i=0; i<num; ++i)
	{
		if(initialcore[i].size() > maxnum)
			maxnum = initialcore[i].size();
	}

	//��������С
	int minnum = 9999;
	for(i=0; i<num; i++)
	{
		if(initialcore[i].size() < minnum)
			minnum = initialcore[i].size();
	}
}

//���ܣ���������
void OutputResults(vector<string> &pnodemap, vector<VecInt> &initialcore, string odatafile)
{
	//ÿ���������ھ��Ľ���ļ�
	ofstream out;
	out.open(odatafile.c_str());
cout<<odatafile<<endl;
	int i, j;
	int index;
	for(i=0; i<initialcore.size(); ++i)
	{
		for(j=0; j<initialcore[i].size(); ++j)
		{
			index = initialcore[i][j];
			out<<pnodemap[index]<<"  ";
		}
		out<<"\n";
	}
	out.close();
}

//���ܣ���ȡ��ʼ���(ȥ��)
void GetResults(vector<string> &pnodemap, vector<VecInt> &initialcore, string odatafile)
{
	MergeCluster(initialcore);  //�ϲ�����ģ��
	OutputResults(pnodemap, initialcore, odatafile);  //��������
}

//���ܣ����������ھ򵰰��ʸ�����
void MiningPC(string ppifile,string complexfile)  
{
	vector<string> pnodemap;  //�����ʽڵ�����
	vector<VecInt> pnetwork;  //�໥��������(�Գ�����)
	vector<float> clusterco;  //�ڵ�ۼ�ϵ������
	vector<VecInt> initialcore;  //��ʼ�˼��� 
	GeneratePNodeMap(pnodemap, ppifile);  //���쵰���ʽڵ�������
	int nodenum = pnodemap.size();  //����ڵ���Ŀ
	int edgenum = 0;  //�������Ŀ
	float aveclustercoeff;  //����ƽ���ۼ�ϵ��
	InitialNetwork(pnetwork, nodenum, nodenum);  //��ʼ���໥��������
	edgenum = ReadData(ppifile, pnodemap, nodenum, pnetwork);  //��ȡ�ļ��������໥��������
	aveclustercoeff = StatisticCC(pnetwork, clusterco, nodenum);  //ͳ��ÿ�������ʽڵ�ľۼ�ϵ��
	float middlecc = MiddleCC(clusterco);
    ConstructCore(pnetwork, initialcore, clusterco, middlecc);  //���ݽڵ�ϵ���Ĵ�С�����ʼ�˼���
	ExtendCore(pnetwork, initialcore, nodenum);  //����չ->�����ʸ�����
   	GetResults(pnodemap, initialcore, complexfile);  //��ȡ��ʼ���(ȥ��)

	}