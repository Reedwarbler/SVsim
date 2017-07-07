#include"file_io.h"
#include<cstdio>
#include <stdlib.h>

using namespace std;

FileIO::FileIO()
{
	this->filename="";
}

FileIO::FileIO(char* filename)
{
	this->filename=filename;
}

/*
Description: get file length. 
*/
unsigned long FileIO::getFileLen()
{
	FILE *pFile = NULL;
	// get the file stream
	pFile=fopen(filename,"rb");
	
	if(!pFile)
	{
		perror(filename);
		cout<<"Open file failed"<<endl; 
		return -1;
	}

	// set the file pointer to end of file
	fseek( pFile, 0, SEEK_END );
	// get the file size
	unsigned long len = ftell( pFile );
	// return the file pointer to begin of file if you want to read it
	rewind( pFile );
	// close stream and release buffer 
	fclose( pFile );

	return len;
}

/*
Description: read the whole file into memory. 
*/
int FileIO::readFileIntoMem(char*& buffer,unsigned long len)
{
	FILE* pfin=NULL;

	//1.1 open file 
	pfin=fopen(filename,"r"); 
	if(!pfin)
	{
		cout<<"Open file failed"<<endl; 
		return -1;
	}

	//1.3 read the whole file into memory  
	unsigned long rst=fread(buffer,1,len,pfin);
	if(rst!=len)
	{
		cout<<"Reading error!"<<endl;
		return -1;
	}
	fflush(pfin);
	fclose(pfin);

	return 0;
}

void FileIO::setFileName(char* filename)
{
	this->filename=filename;
}