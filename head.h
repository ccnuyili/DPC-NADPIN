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

vector<VecStr> Results;  //蛋白质复合物向量集
vector<int>  frequency;  //蛋白质复合物在36个子网中出现的频率

//功能：交换两个变量的值
void Swap(string &str1, string &str2)
{
	string temp;
	temp = str1;
	str1 = str2; 
	str2 = temp;
}

//功能：快速排序方法中的划分方法
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

//功能：将蛋白质节点向量按字母顺序快速排序
void Quicksort(vector<string> &nodemaptemp, int p, int r)
{
	if(p<r)
	{
		int q = Partition(nodemaptemp, p, r);
		Quicksort(nodemaptemp, p, q-1);
		Quicksort(nodemaptemp, q+1, r);
	}
}

//功能：读取静态网络文件，统计节点数目，边数目，创造节点向量
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

	//去掉重复节点
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

//功能：初始化相互作用网络
void InitialNetwork(vector<VecInt> &pnetwork, int rownum, int colnum)
{
	for(int i=0; i<rownum; ++i)
	{
		vector<int> temp(colnum, 0);
		pnetwork.push_back(temp);
	}
}

//功能：给定蛋白质节点，找到其在蛋白质节点向量中对应的编号
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

//功能：读取相互作用文件，构造相互作用(对称)网络
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
		pnetwork[nodeindex1][nodeindex2] += 1;  //对称网络
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

//功能：统计每个蛋白质节点的聚集系数
float StatisticCC(vector<VecInt> &pnetwork, vector<float> &coefficient, int nodenum)
{
	int i, j, m, n, index1, index2;
	int count;  //邻居节点个数
	vector<int> neighbor;  //邻居节点集合
	float cc; //蛋白质节点聚集系数
	int realnum;  //邻居之间实际相互作用边数
	int maxnum;  //邻居之间最大可能相互作用数目
	float avecc;  //网络平均据集系数

	for(i=0; i<nodenum; ++i)
	{
		count = 0;  //邻居节点初始化
		for(j=0; j<nodenum; ++j)
		{
			if(pnetwork[i][j] == 1)
			{
				++count;
				neighbor.push_back(j);  //保存邻居节点(节点下标)
			}
		}
		realnum = 0;
		//统计邻居节点之间的相互作用边数
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
		if(count == 1)  //没有孤立的节点
			cc = 0;	
		else
		{
			maxnum = count *(count-1);
			cc = 2*realnum / (float)(maxnum);
		}
		coefficient.push_back(cc);
		neighbor.clear();  //清空数据，否则溢出
	}

	float sumcc = 0;
	for(i=0; i<coefficient.size(); ++i)
		sumcc += coefficient[i];
	avecc = sumcc / nodenum;  //网络的平均聚集系数
	return avecc;
}

//功能：依据节点系数的构造初始核集合
void ConstructCore(vector<VecInt> &pnetwork, vector<VecInt> &initialcore, vector<float> &clusterco, float avecc)
{
	//聚集系数大于平均聚集系数的节点及其邻居节点一起作为初始核
	int i, j; 
	vector<int> core;  //初始核中保存的是蛋白质节点的下标，这样做节省存储空间
	for(i=0; i<clusterco.size(); ++i)
	{
		if(clusterco[i] > avecc)  //当前点的聚集系数大于平均据集系数
		{
			core.push_back(i);  //当前节点作为初始核的成员
			for(j=0; j<clusterco.size(); ++j)
			{
				if(pnetwork[i][j] == 1)
					core.push_back(j);  //当前节点的邻居节点加入初始核中
			}
			initialcore.push_back(core);  //加入初始核集合中
			core.clear();  //清空临时核成员
		}
	}
}

//功能：核扩展->蛋白质复合物
void ExtendCore(vector<VecInt> &pnetwork, vector<VecInt> &initialcore, int nodenum)
{
	//每个子网被挖掘后的统计文件

	vector<int> candidate;  //当前簇候选邻居节点集合
	vector<int> flag;  //当前簇中节点标记
	vector<int> neighbor_in;  //候选节点簇内邻居节点集合
	vector<int> neighbor_out;  //候选节点簇外邻居节点集合

	int i, j, k, r, index, index1, index2;

	int count = 0;  //记录多少个初始核被扩展

	int initialsize, aftersize, extendlevel;  //初始核大小，扩展后大小，扩展次数

	//依次对每个核进行扩展(可交叠的蛋白质复合物)
	for(i=0; i<initialcore.size(); ++i)
	{
	    initialsize = 0;
		aftersize = 0;
		extendlevel= 0;
		do	
		{   
			initialsize = initialcore[i].size();
			
			//标记集合初始化
			for(j=0; j<nodenum; ++j)
				flag.push_back(0);  
			
			//标记当前簇中的每个节点
			for(j=0; j<initialcore[i].size(); ++j)
			{
				index = initialcore[i][j];  //initialcore[i][j] 是pnodemap中每个蛋白质节点对应的下标
				flag[index] = 1;  //节点在当前簇中 
			}
			
		    //寻找当前簇的邻居节点(候选节点集合)
		    for(j=0; j<initialcore[i].size(); ++j)
			{
		    	index = initialcore[i][j];
		    	for(k=0; k<nodenum; ++k)
				{
			    	if(pnetwork[index][k] == 1 && flag[k] != 1)  //是当前簇的邻居节点，且该节点不在当前簇中
					{
				    	flag[k] = 1;  //标记是为了防止重复的邻居节点
				    	candidate.push_back(k);  //保存候选节点
					}
				}
			}
            
			float maxmargin = 0;  //紧密度差值最大初始化
			int thebest = -1;  //最佳候选节点初始化

			if(candidate.size() == 0)  //不存在候选节点
				break;
			else  //存在候选节点
			{
				flag.clear();  //不清空则造成严重的错误
				//标记集合初始化
	    	    for(j=0; j<nodenum; ++j)
					flag.push_back(0); 
				
				//标记当前簇中的每个节点
				for(j=0; j<initialcore[i].size(); ++j)
				{
					index = initialcore[i][j];  //initialcore[i][j] 是pnodemap中每个蛋白质节点对应的下标
					flag[index] = 1;  //节点在当前簇中 
				}

		    	//依次判断每个候选节点是否能够加入当前核
		    	for(j=0; j<candidate.size(); ++j)
				{ 
					//找到簇内邻居节点并计算簇内邻居节点连接紧密度

			    	index = candidate[j];  //候选节点对应的蛋白质节点的下标

			    	for(k=0; k<initialcore[i].size(); ++k)
					{
				    	index1 = initialcore[i][k];
				    	if(pnetwork[index][index1] == 1)
					    	neighbor_in.push_back(index1);  //当前候选节点在簇内的邻居节点集合	
					}
			    	int sumdegree_in = 0;  //要对每个节点初始化
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
				    float NIAD = sumdegree_in * 1.0 / (neighbor_in.size() * 1.0);  //团内邻居节点平均紧密度 

					//找到簇外邻居节点并计算簇外邻居节点连接紧密度

					float NOAD = 0.0;
			    	for(k=0; k<nodenum; ++k)
					{
				    	if(pnetwork[index][k] == 1 && flag[k] == 0)
				    		neighbor_out.push_back(k);  //当前候选节点在簇外的邻居节点集合
					}
			    	if(neighbor_out.size() != 0)
					{
				    	int sumdegree_out = 0;
				        for(k=0; k<neighbor_out.size(); ++k)
						{
				    		index2 = neighbor_out[k];
				    	    for(r=0; r<nodenum; ++r)
							{
								if(r != index)  //去掉候选节点
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
				    	thebest = candidate[j];  //候选节点中的最佳的一个
					}
			    	neighbor_in.clear();
			    	neighbor_out.clear();
				}
				flag.clear();  //不清空则造成严重的错误
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

//功能：交换数据，函数重载
void Swap(int &a, int &b)
{
	int temp;
	temp = a;
	a = b; 
	b =temp;
}

//功能：冒泡排序
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

//功能：求得两个模块向量的并集
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
		}while(i<temp.size() && temp[i] == temp[low]);  //去掉重复元素
	}
	return unionclu;
}

//功能：返回两个模块向量的重叠度
double OverlapRate(vector<int> clu_1, vector<int> clu_2)
{
	if(clu_1.size() == 0 || clu_2.size() == 0)  //如果其中一个为空，返回0
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

//功能：合并两个重叠度为1.0的模块
int MergeCluster(vector<VecInt> &initialcore)
{
	//每个子网被挖掘后的统计文件

	int i, j, k;
	double mr;  //mergerate 重叠度
	bool flag = true;
	int one, two;
	vector<VecDou> matchrate;  //重叠度矩阵

	InitialMatchrate(matchrate, initialcore.size(), initialcore.size());

	//有效初始化
	for(i=0; i<initialcore.size(); i++)
	{
		for(j=i; j<initialcore.size(); j++)
		{
			mr = OverlapRate(initialcore[i], initialcore[j]);
			matchrate[i][j] = mr;
			matchrate[j][i] = mr;
		}
	}
    
	//找出当前重叠度最大的模块
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
    
			//更新matchrate[][]
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

	//誊写最后的结果集
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

//预测蛋白质复合物结果分析
void ResultsAnalysis(vector<VecInt> &initialcore)
{
    int i;
	int num = initialcore.size();  //挖掘复合物的个数
	int number = 0;

	//复合物的平均大小
	float aversize;
	for(i=0; i<num; ++i)
		number += initialcore[i].size();
	aversize = number / num;
	
	//复合物最大
	int maxnum = 0;
	for(i=0; i<num; ++i)
	{
		if(initialcore[i].size() > maxnum)
			maxnum = initialcore[i].size();
	}

	//复合物最小
	int minnum = 9999;
	for(i=0; i<num; i++)
	{
		if(initialcore[i].size() < minnum)
			minnum = initialcore[i].size();
	}
}

//功能：输出结果集
void OutputResults(vector<string> &pnodemap, vector<VecInt> &initialcore, string odatafile)
{
	//每个子网被挖掘后的结果文件
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

//功能：获取初始结果(去重)
void GetResults(vector<string> &pnodemap, vector<VecInt> &initialcore, string odatafile)
{
	MergeCluster(initialcore);  //合并两个模块
	OutputResults(pnodemap, initialcore, odatafile);  //输出结果集
}

//功能：在子网上挖掘蛋白质复合物
void MiningPC(string ppifile,string complexfile)  
{
	vector<string> pnodemap;  //蛋白质节点向量
	vector<VecInt> pnetwork;  //相互作用子网(对称网络)
	vector<float> clusterco;  //节点聚集系数向量
	vector<VecInt> initialcore;  //初始核集合 
	GeneratePNodeMap(pnodemap, ppifile);  //构造蛋白质节点向量集
	int nodenum = pnodemap.size();  //网络节点数目
	int edgenum = 0;  //网络边数目
	float aveclustercoeff;  //网络平均聚集系数
	InitialNetwork(pnetwork, nodenum, nodenum);  //初始化相互作用网络
	edgenum = ReadData(ppifile, pnodemap, nodenum, pnetwork);  //读取文件，构造相互作用网络
	aveclustercoeff = StatisticCC(pnetwork, clusterco, nodenum);  //统计每个蛋白质节点的聚集系数
	float middlecc = MiddleCC(clusterco);
    ConstructCore(pnetwork, initialcore, clusterco, middlecc);  //依据节点系数的大小构造初始核集合
	ExtendCore(pnetwork, initialcore, nodenum);  //核扩展->蛋白质复合物
   	GetResults(pnodemap, initialcore, complexfile);  //获取初始结果(去重)

	}