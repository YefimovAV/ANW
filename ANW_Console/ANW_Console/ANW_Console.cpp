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
int FindFirstNodal(int iW, int iH, unsigned char *bw, bool *visited, int currentIndex) {
	int neighborCount;
	for (int i = 0; i < iW; ++i) 
		for (int j = iH - 1; j > 0; --j) {
			if (bw[i + j * iW] == 1) {
				neighborCount = 0;
				currentIndex = i + j * iW;
				currentIndex - iW > 0 && bw[currentIndex - iW] == 1 ? ++neighborCount : 0;
				(currentIndex - iW + 1) % iW != 0 && currentIndex - iW + 1 > 0 && bw[currentIndex - iW + 1] == 1 ? ++neighborCount : 0;
				(currentIndex + 1) % iW != 0 && currentIndex + 1 < iW * iH && bw[currentIndex + 1] == 1 ? ++neighborCount : 0;
				(currentIndex + iW + 1) % iW != 0 && currentIndex + iW + 1 < iW * iH && bw[currentIndex + iW + 1] == 1 ? ++neighborCount : 0;
				currentIndex + iW < iW * iH && bw[currentIndex + iW] == 1 ? ++neighborCount : 0;
				currentIndex % iW != 0 && currentIndex + iW - 1 < iW * iH && bw[currentIndex + iW - 1] == 1 ? ++neighborCount : 0;
				currentIndex % iW != 0 && currentIndex - 1 > 0 && bw[currentIndex - 1] == 1 ? ++neighborCount : 0;
				currentIndex % iW != 0 && currentIndex - iW - 1 > 0 && bw[currentIndex - iW - 1] == 1 ? ++neighborCount : 0;
				if (neighborCount > 2) return currentIndex;
			}		
		}
	return currentIndex;	
}

int BuildScanArray(int iW, int iH, unsigned char *bw, bool *visited, int currentIndex, int *scan, int &iterator) {
	int nextIndex;
	nextIndex = currentIndex - iW;
	visited[currentIndex] = 1;
	if (nextIndex > 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else nextIndex = currentIndex - iW + 1;
	if (nextIndex > 0 && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else nextIndex = currentIndex + 1;
	if (nextIndex < iW * iH && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else nextIndex = currentIndex + iW + 1;
	if (nextIndex < iW * iH && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else nextIndex = currentIndex + iW;
	if (nextIndex < iW * iH && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else nextIndex = currentIndex + iW - 1;
	if (nextIndex < iW * iH && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else nextIndex = currentIndex - 1;
	if (nextIndex < iW * iH && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else nextIndex = currentIndex - iW - 1;
	if (nextIndex > 0 && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		currentIndex = nextIndex;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		nextIndex = BuildScanArray(iW, iH, bw, visited, currentIndex, scan, iterator);
	}
	else return currentIndex;
}

void SymbolScan(int iW, int iH, unsigned char *bw, int *scan) {
	int startIndex;
	int firstNodal;
	int iterator = 0;
	bool *visited = new bool [MAX_IMAGE_SIZE];
	for (int i = 0; i < iW * iH; ++i)
		visited[i] = 0;
	for (int i = 0; i < iW; ++i) 
		for (int j = iH - 1; j > 0; --j)	
			if (bw[i + j * iW] == 1) { startIndex = i + j * iW; break; }
	firstNodal = FindFirstNodal(iW, iH, bw, visited, startIndex);
	BuildScanArray(iW, iH, bw, visited, firstNodal, scan, iterator);
	cout << "x: " << firstNodal % iW << " y: " << (firstNodal - firstNodal % iW) / iW << endl;
	delete[] visited;
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
	ofstream pathFile;
	string imagePath;
	string imageName;
	string fullNameResult;
	int iterator = 0;
	string strIterator;
	imageFile.open(filePath);
	if (mode == "teach/")
		pathFile.open("result/SymbolPathTeach.txt");
	if (mode == "recognition/")
		pathFile.open("result/SymbolPathTeach.txt");
	while(!imageFile.eof()) {
		ostringstream streamStrIterator;
		unsigned char *blackAndWhite = new unsigned char[MAX_IMAGE_SIZE];
		int *scanArray = new int[2 * MAX_IMAGE_SIZE];
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
			for (int i = 0; i < 2 * imageWidth * imageHeight; ++i) {
				pathFile << scanArray[i] << endl;
			}
			getline(imageFile, imageName);
			streamStrIterator << iterator;
			strIterator = streamStrIterator.str();
			fullNameResult = "result/images/" + mode + imageName + strIterator + ".bmp";
			pmap.save_as_bmp(fullNameResult.c_str());
			++iterator;
		}
		else break;
	}
	pathFile.close();
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

