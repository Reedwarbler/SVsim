
#include<iostream>

class FileIO
{
public:
	FileIO(char* filename);
	FileIO();

public:
	unsigned long getFileLen(); //get file length
	int readFileIntoMem(char*& buffer,unsigned long len);//read the whole file into memory 

	void setFileName(char* filename);

private:
	char* filename; //file name
	
};

