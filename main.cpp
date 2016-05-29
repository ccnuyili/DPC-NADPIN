#include"head.h"

int main(int args,char* argv[])
{  
	if (args != 3) {
		cout << "parameter error!" << endl;
		return 0;
	}
	
	istringstream istr(argv[1]);
	string file;
	istr >> file;
	
	istr.clear();
	istr.str(argv[2]);
	string complexfile;
	istr>>complexfile;
	MiningPC(file,complexfile);  
	cout<<"The program is over!\n\n";
    return 1;
}