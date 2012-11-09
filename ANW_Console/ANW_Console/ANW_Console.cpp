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
#include "fft.h"
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
const int MAX_SYMBOL_COUNT = 3;
const int MAX_IMAGE_SIZE = 1000 * 1000;
const int SIGNAL_COUNT = 20;
const string TEACH_FOLDER = "teach/";
const string RECOGNITION_FOLDER = "recognition/";
const string PATH_LIST_TEACH = "result/pathListTeach.txt";
const string PATH_LIST_RECOGNITION = "result/pathListRecognition.txt";
const string TEACH_MODE = "teach/";
const string RECOGNITION_MODE = "recognition/";
const string IMAGE_FOLDER = "result/images/";
const string PATHES_FOLDER = "result/pathes/";
const string FFT_FOLDER = "result/fft/";

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
		for (int j = iH - 1; j >= 0; --j) {
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
	if (nextIndex > 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	nextIndex = currentIndex - iW + 1;
	if (nextIndex > 0 && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	nextIndex = currentIndex + 1;
	if (nextIndex < iW * iH && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	nextIndex = currentIndex + iW + 1;
	if (nextIndex < iW * iH && nextIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	nextIndex = currentIndex + iW;
	if (nextIndex < iW * iH && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	nextIndex = currentIndex + iW - 1;
	if (nextIndex < iW * iH && currentIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	nextIndex = currentIndex - 1;
	if (nextIndex > 0 && currentIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	nextIndex = currentIndex - iW - 1;
	if (nextIndex > 0 && currentIndex % iW != 0 && bw[nextIndex] == 1 && visited[nextIndex] == 0) {
		visited[nextIndex] = 1;
		scan[iterator] = currentIndex % iW;
		++iterator;
		scan[iterator] = (currentIndex - currentIndex % iW) / iW;
		++iterator;
		BuildScanArray(iW, iH, bw, visited, nextIndex, scan, iterator);
	}
	return 0;
}

void SymbolScan(int iW, int iH, unsigned char *bw, int *scan, int &whitePixelCount) {
	int startIndex;
	int firstNodal;
	int counter = 0;
	bool *visited = new bool [MAX_IMAGE_SIZE];
	for (int i = 0; i < iW * iH; ++i)
		visited[i] = 0;
	for (int i = 0; i < iW * iH; ++i)
		if (bw[i] == 1) { ++counter;}
	for (int i = 0; i < iW; ++i) 
		for (int j = iH - 1; j > 0; --j)	
			if (bw[i + j * iW] == 1) { startIndex = i + j * iW; i = iW; break; }				//костыль i = iW для выхода из внешнего цикла, не goto же использовать.
	firstNodal = FindFirstNodal(iW, iH, bw, visited, startIndex);
	visited[firstNodal] = 1;
	scan[whitePixelCount] = firstNodal % iW;
	++whitePixelCount;
	scan[whitePixelCount] = (firstNodal - firstNodal % iW) / iW;
	++whitePixelCount;
	BuildScanArray(iW, iH, bw, visited, firstNodal, scan, whitePixelCount);
	delete[] visited;
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
	int whitePixelCount;
	string strIterator;
	imageFile.open(filePath);
	while(!imageFile.eof()) {
		ostringstream streamStrIterator;
		unsigned char *blackAndWhite = new unsigned char[MAX_IMAGE_SIZE];
		int *scanArray = new int[2 * MAX_IMAGE_SIZE];
		whitePixelCount = 0;
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
			SymbolScan(imageWidth, imageHeight, blackAndWhite, scanArray, whitePixelCount);
			getline(imageFile, imageName);
			streamStrIterator << iterator;
			strIterator = streamStrIterator.str();
			pathFile.open(PATHES_FOLDER + mode + imageName + strIterator + ".txt");
			for (int i = 0; i < whitePixelCount; ++i) {
				pathFile << scanArray[i] << endl;
			}
			fullNameResult = IMAGE_FOLDER + mode + imageName + strIterator + ".bmp";
			pmap.save_as_bmp(fullNameResult.c_str());
			++iterator;
			pathFile.close();
		}
		else {
			delete[] blackAndWhite;
			delete[] scanArray;
			break;
		}
	}
	return iterator;
}

vector<double> FFT(vector<double> dataDouble, int size) {
	ap::complex_1d_array complexData;
	ap::complex_1d_array fftData;
	vector<double> fftDoubleData;
	complexData.setbounds(0, size);
	fftData.setbounds(0, size);
	for(int i = 0; i < size; ++i) {
		complexData(i).x = dataDouble[2 * i];
		complexData(i).y = dataDouble[2 * i + 1];
	}
	for(int i = 1; i < size; ++i)
		for(int j = 0; j < size; ++j) {
			fftData(i).x += complexData(j).x * cos(-2 * ap::pi() * i * j / size) - complexData(j).y * sin(-2 * ap::pi() * i * j / size);	
			fftData(i).y += complexData(j).x * sin(-2 * ap::pi() * i * j / size) + complexData(j).y * cos(-2 * ap::pi() * i * j / size);
		}
	fftData(0) = 0;
	double scalingFactor = sqrt(pow(ap::abscomplex(fftData(1)), 2) + pow(ap::abscomplex(fftData(size - 1)), 2));
	for(int i = 0; i < size; ++i) {
		fftData(i) /= scalingFactor;
		fftDoubleData.push_back(fftData(i).x);
		fftDoubleData.push_back(fftData(i).y);
	}
	return fftDoubleData;
}

vector<double> ReadFile(string filePath, int &size) {
	vector<double> dataDouble;
	ifstream dataFile;
	string thisData;
	dataFile.open(filePath);
	getline(dataFile, thisData);
	while ( thisData.length() > 0) {
		dataDouble.push_back(boost::lexical_cast<double>(thisData));
		getline(dataFile, thisData);
		++size;
	}
	dataFile.close();
	return dataDouble;
}

void FormSignals(string coordinatesFolder, string fftFolder, string mode) {
	namespace fs = boost::filesystem;
	string coordinatesSymbolPath;
	string imageName;
	ofstream fftFile;
	vector<double> signals;
	int whitePixelCount;
	ap::complex_1d_array z; 
	for (fs::directory_iterator it(coordinatesFolder), end; it != end; ++it) {
		if (it->path().extension() == ".txt") {
			int iterator = 0;
			whitePixelCount = 0;
			coordinatesSymbolPath = "";
			while(it->path().string()[iterator]) {
				if (it->path().string()[iterator] == '"') continue;
				string temp = string(&it->path().string()[iterator],0,1);
				coordinatesSymbolPath += temp;
				++iterator;
			}
			signals = ReadFile(coordinatesSymbolPath, whitePixelCount);
			signals = FFT(signals, whitePixelCount / 2);
			imageName = it->path().filename().string();
			fftFile.open(fftFolder + mode + imageName);
			if (whitePixelCount >= 4 * SIGNAL_COUNT) {
				for (int i = 3; i < 3 + 2 * SIGNAL_COUNT; ++i)											// Сигналы формируются начиная с первого элемента ряда Фурье, а не с нулевого.
					fftFile << signals[i] << endl;
				for (int i = whitePixelCount - 1; i > whitePixelCount - 1 - 2 * SIGNAL_COUNT; --i)		// В обратную сторону сигналы формируются начиная с последнего элемента ряда Фурье
					fftFile << signals[i] << endl;
			}
			fftFile.close();
		}	
	}
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
			ReadFolder(TEACH_FOLDER, PATH_LIST_TEACH);
			cout << "Image processing... ";
			imageCount = ImageProcessing(PATH_LIST_TEACH, TEACH_MODE);
			cout << "ready" <<endl << "Processed " << imageCount << " images." << endl << "Handling coordinates...";
			FormSignals(PATHES_FOLDER + TEACH_FOLDER, FFT_FOLDER, TEACH_MODE);
			cout << "ready" << endl << "For training used " << 2 * SIGNAL_COUNT << " signals." <<endl;
		}
		else if(mode == "recognition") {
			cout << "You choose recognition mode" << endl; 
			cout << "Read folder..." << endl;
			ReadFolder(RECOGNITION_FOLDER, PATH_LIST_RECOGNITION);
			cout << "Image processing... ";
			imageCount = ImageProcessing(PATH_LIST_RECOGNITION, RECOGNITION_MODE);
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

