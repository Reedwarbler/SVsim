
#include"public_func.h"
#include<sstream>

/*
Description: Convert int to string.
*/
std::string PubFuncs::cvtInt2Str(int i)
{
	std::stringstream ss;
	ss<<i;
	return ss.str();
}


/*
Description:  Convert string to int.
*/
int PubFuncs::cvtStr2Int(std::string str)
{
	int num=0;
	int len = str.length();
	
	for(int i=0;i<len;i++)
	{
		num*=10;
		num+=(str[i]-'0');
	}

	return num;
}