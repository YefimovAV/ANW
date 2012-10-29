// ANW_Console.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include "agg_win32_bmp.h"
#include <boost/filesystem.hpp>

using namespace std;
const int MAX_SYMBOL_COUNT = 3;

void ReadFolder(string loadPath, string savePath) {
	namespace fs = boost::filesystem;
	ofstream pathList;
	ifstream fileNameStream;
	string fileName;
	string symbolValue;
	pathList.open(savePath);
	cout << "Used files:" << endl; 
	for (fs::directory_iterator it(loadPath), end; it != end; ++it) {
		if (it->path().extension() == ".bmp") {
			int iterator = 0;
			symbolValue = "";
			fileName = it->path().filename().string();
			while(fileName[iterator] < '0' || fileName[iterator] > '9') {
				if (fileName[iterator] == '"') continue;
				string temp = string(&fileName[iterator],0,1);
				symbolValue += temp;
				++iterator;
			}
			cout << it->path().filename() << " " << symbolValue << endl;;
			pathList << *it << endl << symbolValue << endl;
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
			ReadFolder("teach/", "result/pathListTeach.txt");  
		}
		else if(mode == "recognition") {
			cout << "You choose recognition mode" << endl; 
			ReadFolder("recognition/", "result/pathListRecognition.txt");
		}
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

