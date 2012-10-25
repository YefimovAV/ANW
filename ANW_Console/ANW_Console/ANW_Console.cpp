// ANW_Console.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "agg_win32_bmp.h"
#include <boost/filesystem.hpp>

using namespace std;

void ReadFolder() {

}

int _tmain(int argc, _TCHAR* argv[])
{
	string mode;
	for(;;) {
	cout << "Choose a mode (teach/recognition): ";
	cin >> mode;
	if(mode == "teach") { cout << "You choose Teach mode" << endl; break; }
		//Teach(), break; 
	else if(mode == "recognition") { cout << "You choose recognition mode" << endl; break; }
		//Recognition(), break;
	else cout << "Bad value, choose again" << endl;
	}
	return 0;
}

