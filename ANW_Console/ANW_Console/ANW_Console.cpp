// ANW_Console.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "agg_win32_bmp.h"
#include <boost/filesystem.hpp>

using namespace std;

void ReadFolder(string path) {
	namespace fs = boost::filesystem;
	ofstream pathList;
	pathList.open("result/pathListTeach.txt");
	cout << "Using files:" << endl; 
	for (fs::directory_iterator it(path), end; it != end; ++it) {
        if (it->path().extension() == ".bmp") {
			cout << it->path().filename() << endl;
			pathList << *it << endl;
		}
	}
	pathList.close();
}

int _tmain(int argc, _TCHAR* argv[])
{
	string mode;
	char answer = 'y';
	for(;;) {
	if(answer == 'y') {
		cout << "Choose a mode (teach/recognition): ";
		cin >> mode;
		if(mode == "teach") { 
			cout << "You choose Teach mode" << endl; 
			ReadFolder("teach/");  
		}
		else if(mode == "recognition") { cout << "You choose recognition mode" << endl; }
			//Recognition(), break;
		else cout << "Bad value, choose again" << endl;
	}
	cout << "Continue (y/n>? ";
	cin >> answer;
	if(answer == 'y') continue;
	else if (answer == 'n') break;
	else cout << "Bad value, choose again" << endl;
	}
	return 0;
}

