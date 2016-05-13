#include "mic_model.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;
//m_T
int training_time = 1000;

int main(int argc, char **argv )
{
	//int trainYear = 10;
	vector< vector<int> > citations;

	ifstream inFile( "/home/bidisha/workspace/training.txt" );

	if ( inFile.fail() )
	{
		cout << "fail to open files.\n";
		exit(0);
	}


	int cc;
	vector<int> citation;
	while ( inFile >> cc )
	{
		if ( cc > 0 && cc <= training_time )
			citation.push_back( cc );
	}
	citations.push_back( citation );

	inFile.close();

	
	CMicModel model( training_time );
	ofstream outFile( "a.txt" );

	for ( size_t ii = 0; ii < citations.size(); ++ii )
	{
		model.parameter_estimation(citations[ii]);
		cout << ii+1 << "\t" << citations[ii].size() << "\t" << model.get_alpha() << "\t" << model.get_beta() << "\t" << model.get_gamma() << "\n";
		//outFile << citations[ii].size() << endl;
	}

	outFile.close();
	
	return 0;
}


