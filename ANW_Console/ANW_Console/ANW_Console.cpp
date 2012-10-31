// ANW_Console.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "agg_win32_bmp.h"
#include "NoiseEdit.h"
#include "Skelet.h"
#include <boost/filesystem.hpp>

using namespace std;
const int MAX_SYMBOL_COUNT = 3;
const int MAX_IMAGE_SIZE = 1000 * 1000;

void k_to_b(unsigned char* k, unsigned char* b,unsigned rowLength, int wI, int hI)
{
	for(int i = 0; i < hI; i++)
	{
		for(int j = 0; j < wI; j++)
		{
			b[i*rowLength + j / 8] = 0;
		}
	}

	for(int i = 0; i < hI; i++)
	{
		for(int j = 0; j < wI; j++)
		{
			b[i*rowLength + j / 8] |= !!k[(hI - i - 1)*wI + j] << (7 - j % 8);
		}
	}
}

void b_to_k(unsigned char* b, unsigned char* k,unsigned rowLength, int wI, int hI)
{
	for(int i = 0; i < hI; i++)
	{
		for(int j = 0; j < wI; j++)
		{
			k[(hI - i - 1)*wI + j] = !!(b[i*rowLength + j / 8] & (1 << (7 - j % 8)));
		}
	}
}

void SymbolScan(int iW, int iH, unsigned char *bw, unsigned char *scan) {
	int startIndex;
	bool *visited = new bool [MAX_IMAGE_SIZE];
	for (int i = 0; i < iW * iH; ++i)
		visited[i] = 0;
	for (int i = 0; i < iW; ++i) 
		for (int j = iH - 1; j > 0; --j) 
			if (bw[i + j * iW] == 1) { startIndex = i + j * iW; break; }
	

}

void ReadFolder(string loadPath, string savePath) {
	namespace fs = boost::filesystem;
	ofstream pathList;
	string filePath;
	string symbolValue;
	pathList.open(savePath);
	cout << "Used files:" << endl; 
	for (fs::directory_iterator it(loadPath), end; it != end; ++it) {
		if (it->path().extension() == ".bmp") {
			int iterator = 0;
			symbolValue = "";
			filePath = "";
			while(it->path().filename().string()[iterator] < '0' || it->path().filename().string()[iterator] > '9') {
				if (it->path().filename().string()[iterator] == '"') continue;
				string temp = string(&it->path().filename().string()[iterator],0,1);
				symbolValue += temp;
				++iterator;
			}
			iterator = 0;
			while(iterator < it->path().string().length()) {
				if (it->path().string()[iterator] == '"') continue;
				string temp = string(&it->path().string()[iterator],0,1);
				filePath += temp;
				++iterator;
			}
			cout << it->path().filename() << " " << symbolValue << endl;
			pathList << filePath << endl << symbolValue << endl;
		}
	}
	pathList.close();
}

int ImageProcessing (string filePath, string mode) {
	agg::pixel_map pmap;	
	ifstream imageFile;
	ifstream scanFile;
	string imagePath;
	string imageName;
	string fullNameResult;
	int iterator = 0;
	string strIterator;
	imageFile.open(filePath);
	while(!imageFile.eof()) {
		ostringstream streamStrIterator;
		unsigned char *blackAndWhite = new unsigned char[MAX_IMAGE_SIZE];
		unsigned char *scanArray = new unsigned char[MAX_IMAGE_SIZE];
		getline(imageFile, imagePath);
		if(imagePath.length() > 0) {
			pmap.load_from_bmp(imagePath.c_str());
			unsigned char* imageBuffer = pmap.buf();
			int imageWidth = pmap.width();
			int imageHeight = pmap.height();
			unsigned rowLength = pmap.calc_row_len(imageWidth,1);
			b_to_k(imageBuffer, blackAndWhite, rowLength, imageWidth, imageHeight);
			blackAndWhite = NoiseEdit(blackAndWhite,imageWidth,imageHeight);
			k_to_b(blackAndWhite, imageBuffer, rowLength, imageWidth, imageHeight);
			skelet(imageBuffer, imageHeight, rowLength);
			b_to_k(imageBuffer, blackAndWhite, rowLength, imageWidth, imageHeight);
			SymbolScan(imageWidth, imageHeight, blackAndWhite, scanArray);
			getline(imageFile, imageName);
			streamStrIterator << iterator;
			strIterator = streamStrIterator.str();
			fullNameResult = "result/images/" + mode + imageName + strIterator + ".bmp";
			pmap.save_as_bmp(fullNameResult.c_str());
			++iterator;
		}
		else break;
	}
	return iterator;
}

int _tmain(int argc, _TCHAR* argv[])
{
	string mode;
	char answer = 'y';
	int imageCount;
	for(;;) {
	if(answer == 'y') {
		cout << "Choose a mode (teach/recognition): ";
		cin >> mode;
		if(mode == "teach") { 
			cout << "You choose Teach mode" << endl; 
			cout << "Read folder..." << endl;
			ReadFolder("teach/", "result/pathListTeach.txt");
			cout << "Image processing... ";
			imageCount = ImageProcessing("result/pathListTeach.txt", "teach/");
			cout << "ready" <<endl << "Processed " << imageCount << " images." << endl;
		}
		else if(mode == "recognition") {
			cout << "You choose recognition mode" << endl; 
			cout << "Read folder..." << endl;
			ReadFolder("recognition/", "result/pathListRecognition.txt");
			cout << "Image processing... ";
			imageCount = ImageProcessing("result/pathListRecognition.txt", "recognition/");
			cout << "ready" <<endl << "Processed " << imageCount << " images." << endl;
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

