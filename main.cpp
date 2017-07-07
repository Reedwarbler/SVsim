#include<cstdio>
#include<iostream>
#include<string>
#include<fstream>
#include<stdlib.h>
#include <time.h>

#include <ctype.h>

//for linux
#include <unistd.h> 

#include"vcf_parse.h"
#include"file_io.h"
#include"public_func.h"

using namespace std;

long loadFileIntoMem(FileIO& fin, char*& buffer, char* filename)
{
	fin.setFileName(filename);

	long len=fin.getFileLen();
	if(len<0) return -1;
	buffer=new char[(sizeof(char))*len];
	fin.readFileIntoMem(buffer,len);

	return len;
}


int gnrtDeletionHap(char* vcffilename, int cnt_indi, char* chrom_name, int cnt_sim)
{
	FileIO fin; 
	char* buffer=NULL;

	//parse vcf file 
	long len_vcf=loadFileIntoMem(fin, buffer,vcffilename); 

	if(len_vcf<0) return -1;
	VCFParse vcfpser(buffer,len_vcf,cnt_indi); 
	vcfpser.parseBrkpntGenotype();
	delete[] buffer; 
	buffer=NULL;

	//first read in the whole chromsome.
	ifstream fin_chr;
	fin_chr.open(chrom_name);
	//read the first line
	char firstline[500];
	fin_chr.getline(firstline,500);

	string ref="";
	string aline;
	while(fin_chr>>aline)
	{
		ref+=aline;
	}

	int samples=cnt_sim;// how many individuals will be simulated.
	for(int s=1;s<samples;s++)
	{//how many sampels will be generated.
		string hap1="";
		string hap2="";

		int cur_pos=0;
		
		int del_start=0;
		int del_end=0;

		int size=vcfpser.brkpnts.size();
		
		for(int i=0;i<size;i++)
		{
			del_start=vcfpser.brkpnts[i].first;
			del_end=vcfpser.brkpnts[i].second;

			int g1=vcfpser.genotypes[i][s].first;
			int g2=vcfpser.genotypes[i][s].second;
			bool bg1,bg2;
			if(g1==0 && g2==0)
			{//g==0
				bg1=false;
				bg2=false;
			}
			else if(g1==1 && g2==1)
			{//g==2
				bg1=true;
				bg2=true; 
			}
			else
			{//g==1
				srand(time(NULL));
				if(rand()%2==0) 
				{
					bg1=true;
					bg2=false;
				}
				else
				{
					bg1=false;
					bg2=true;
				}
			}

			if(bg1==true)
			{
				hap1+=ref.substr(cur_pos,del_start-cur_pos);
			}
			else
			{//not delete
				hap1+=ref.substr(cur_pos,del_end-cur_pos+1);
			}
			
			if(bg2==true)
			{
				hap2+=ref.substr(cur_pos,del_start-cur_pos);
			}
			else
			{	
				hap2+=ref.substr(cur_pos,del_end-cur_pos+1);
			}
			cur_pos=del_end+1;

		}
		hap1+=ref.substr(cur_pos,ref.length()-cur_pos);
		hap2+=ref.substr(cur_pos,ref.length()-cur_pos);

		//output haplotype
		string filename1="indi"+PubFuncs::cvtInt2Str(s+1)+"_hap1.fa";
		string filename2="indi"+PubFuncs::cvtInt2Str(s+1)+"_hap2.fa";
		ofstream fout_hap1;
		fout_hap1.open(filename1.c_str());
		ofstream fout_hap2;
		fout_hap2.open(filename2.c_str());
		
		fout_hap1<<">"<<chrom_name<<endl;
		fout_hap2<<">"<<chrom_name<<endl;
		for(int h=0;h<hap1.length();h++)
		{
			if(h!=0 && h%50==0)
				fout_hap1<<endl;
			fout_hap1<<hap1.at(h);
			
		}

		for(int h=0;h<hap2.length();h++)
		{
			if(h!=0 && h%50==0)
				fout_hap2<<endl;
			fout_hap2<<hap2.at(h);
			
		}

		fout_hap1.close();
		fout_hap2.close();
	}
}


void gnrtRandStr(int len, string& str)
{
	str="";
	srand( (unsigned)time( NULL ) );

	for(int i=0;i<len;i++)
	{
		int num = rand() % 4;
		if(num==0)
			str+="A";
		else if(num==1)
			str+="C";
		else if(num==2)
			str+="G";
		else
			str+="T";
	}
}


void gnrtGntp(int row, int col, int rg0, int rg1, int rg2)
{
	ofstream fout;
	fout.open("temp_gntp.txt");

	int total=rg0 + rg1 + rg2;
	srand( (unsigned)time( NULL ) );
	string g="";
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			int num=rand() % total;
			if(num<rg0)
				g="0/0";
			else if(num>=rg0 && num<(rg0+rg1))
				g="0/1";
			else 
				g="1/1";	

			fout<<g<<" ";
		}
		fout<<endl;
	}
	fout.close();
}


int gnrtInsertionHap(char* inst_file, char* chrom_file, int cnt_indi, int rg0, int rg1, int rg2)
{	
	//read in insertion sites
	ifstream fin_inst_sites;
	fin_inst_sites.open(inst_file);
	string chrom;
	int start, inst_len;
	string inst_name;
	vector<pair<int,int> > inst_sites;
	while(fin_inst_sites>>chrom>>start>>inst_len)
	{
		inst_sites.push_back(std::make_pair(start,inst_len));
	}
	fin_inst_sites.close();

	int row=inst_sites.size();//number of insertions
	int col=cnt_indi; //number of individuals 

	gnrtGntp(row, col, rg0, rg1, rg2); //generate genotypes

	ifstream fgntp;
	fgntp.open("temp_gntp.txt");

	string sgntp;
	int** genotypes=new int*[row];
	for(int i=0;i<row;i++)
		genotypes[i]=new int[col];

	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			fgntp>>sgntp;
			if(sgntp=="0/0")
				genotypes[i][j]=0;
			else if(sgntp=="0/1")
				genotypes[i][j]=1;
			else
				genotypes[i][j]=2;
		}
	}
	fgntp.close();

	//first read in the whole chromosome.
	ifstream fin_chr;
	fin_chr.open(chrom_file);
	//read the first line
	char firstline[500];
	fin_chr.getline(firstline,500);

	string ref="";
	string aline;
	while(fin_chr>>aline)
	{
		ref+=aline;
	}
	
	int samples=col;
	for(int s=0;s<samples;s++)
	{//how many sampels will be generated.
		string hap1="";
		string hap2="";
		int cur_pos=0;
		int inst_start=0;
		int inst_length=0;

		int size=inst_sites.size();
		
		for(int i=0;i<size;i++)
		{
			inst_start=inst_sites[i].first;
			inst_length=inst_sites[i].second;

			int g=genotypes[i][s];
			
			bool bg1,bg2;
			if(g==0)
			{//g==0
				bg1=false;
				bg2=false;
			}
			else if(g==2)
			{//g==2
				bg1=true;
				bg2=true; 
			}
			else
			{//g==1
				srand(time(NULL));
				if(rand()%2==0) 
				{
					bg1=true;
					bg2=false;
				}
				else
				{
					bg1=false;
					bg2=true;
				}
			}

			string insert_str="";
			gnrtRandStr( inst_length, insert_str);//generate new string 

			if(bg1==true)
			{
				hap1+=ref.substr(cur_pos,inst_start-cur_pos);
				hap1+=insert_str;
			}
			else
			{
				hap1+=ref.substr(cur_pos,inst_start-cur_pos);
			}
			
			
			if(bg2==true)
			{
				hap2+=ref.substr(cur_pos,inst_start-cur_pos);
				hap2+=insert_str;
			}
			else
			{
				hap2+=ref.substr(cur_pos,inst_start-cur_pos);
			}
			
			cur_pos=inst_start+1;
		}
		hap1+=ref.substr(cur_pos,ref.length()-cur_pos);
		hap2+=ref.substr(cur_pos,ref.length()-cur_pos);

		//output haplotype
		string filename1="indi"+PubFuncs::cvtInt2Str(s+1)+"_hap1.fa";
		string filename2="indi"+PubFuncs::cvtInt2Str(s+1)+"_hap2.fa";
		ofstream fout_hap1;
		fout_hap1.open(filename1.c_str());
		ofstream fout_hap2;
		fout_hap2.open(filename2.c_str());
		
		fout_hap1<<">"<<chrom_file<<endl;
		fout_hap2<<">"<<chrom_file<<endl;
		for(int h=0;h<hap1.length();h++)
		{
			if(h!=0 && h%50==0)
				fout_hap1<<endl;
			fout_hap1<<hap1.at(h);
			
		}

		for(int h=0;h<hap2.length();h++)
		{
			if(h!=0 && h%50==0)
				fout_hap2<<endl;
			fout_hap2<<hap2.at(h);
			
		}

		fout_hap1.close();
		fout_hap2.close();
	}

	//release
	for(int i=0;i<row;i++)
		delete[] genotypes[i];

	delete[] genotypes;

	return 0;
}


void usage() 
{
	cout<<"Usage: ./SVHapGnrtor [options] "<<endl;
	cout<<endl;
	cout<<"                  -D/I      For deletions/insertions"<<endl;
	cout<<"                  -h        help"<<endl;
	cout<<endl;
	cout<<"Options for Deletion;"<<endl;
	cout<<"                  -v FILE   Deletion postion and genotype file"<<endl;
	cout<<"                  -i INT    Total number of individuals"<<endl;
	cout<<"                  -c FILE   Chromosome file"<<endl;
	cout<<"                  -s INT    Number of simulated individuals"<<endl;
	cout<<endl;
	cout<<"Options for Insertion;"<<endl;
	cout<<"                  -v FILE   Insertion position and length file"<<endl;
	cout<<"                  -i INT    Number of individuals to simulate"<<endl;
	cout<<"                  -c FILE   Chromosome file"<<endl;
	cout<<"                  -z INT    % of 0/0"<<endl;
	cout<<"                  -o INT    % of 0/1"<<endl;
	cout<<"                  -t INT    % of 1/1, 0/0 add 0/1 and 1/1 should be 100."<<endl;
}


int main(int argc, char* argv[])
{
	int c;
	opterr = 0;
	bool bdel=true;
	bool flag_DI=false,flag_v=false,flag_i=false,flag_c=false,flag_s=false, flag_z=false, flag_o=false, flag_t=false;
	char* gntp_file;
	char* chrom_file;
	int cnt_indi=0;
	int cnt_sim=0;
	int pct_zero=0;
	int pct_one=0;
	int pct_two=0;

	while ((c = getopt(argc, argv, "DIv:i:c:s:z:o:t:h")) != -1)
    switch (c)
    {
	case 'D': //Call genotype of deletion 
		bdel=true;
		flag_DI=true;
        break;
    case 'I': //Call genotype of insertion 
        bdel=false;
		flag_DI=true;
        break;
    case 'v': //genotype file 
        gntp_file = optarg;
		flag_v=true;
        break;
    case 'i': //number of individuals
		cnt_indi=PubFuncs::cvtStr2Int((string)optarg);
        flag_i=true;
		break;
	case 'c': //chrom file
		chrom_file=optarg;
		flag_c=true;
        break;
	case 's'://number of simulated individuals
		cnt_sim=PubFuncs::cvtStr2Int((string)optarg);
		flag_s=true;
        break;

	case 'z'://number of simulated individuals
		pct_zero=PubFuncs::cvtStr2Int((string)optarg);
		flag_z=true;
        break;
	case 'o'://number of simulated individuals
		pct_one=PubFuncs::cvtStr2Int((string)optarg);
		flag_o=true;
        break;
	case 't'://number of simulated individuals
		pct_two=PubFuncs::cvtStr2Int((string)optarg);
		flag_t=true;
        break;		

	case 'h': //help 
		usage();
        break;
    case '?':
        if (optopt == 'c')
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
            fprintf (stderr, "Unknown option -%c'.\n", optopt);
        else
            fprintf (stderr,"Unknown option character \\x%x'.\n",optopt);
        usage();
		return 1;
    default:
		usage();
        abort ();
    }

	if(bdel==true)
	{//deletion
		if(flag_DI==true && flag_v==true && flag_i==true && flag_c==true && flag_s==true)
		{
			cout<<"Deletion Individual Haplotype Generating..."<<endl;
			gnrtDeletionHap(gntp_file, cnt_indi, chrom_file, cnt_sim);
		}
		else
		{
			cout<<"Invalidated Deletion Haplotype Generation Command."<<endl;
			usage();
		}
	}
	else
	{//insertion
		if(flag_DI==true && flag_v==true && flag_i==true && flag_c==true && flag_z==true && flag_o==true && flag_t==true)
		{
			cout<<"Insertion Individual Haplotype Generating..."<<endl;
			gnrtInsertionHap(gntp_file, chrom_file, cnt_indi, pct_zero, pct_one, pct_two);
		}
		else
		{
			cout<<"Invalidated Insertion Haplotype Generation Command."<<endl;
			usage();
		}
	}
	
	return 0;
}
