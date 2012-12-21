// ANW_Console.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <locale>
#include "agg_win32_bmp.h"
#include "NoiseEdit.h"
#include "Skelet.h"
#include "fft.h"
#include "omp.h"
#include "windows.h"
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace std;
const int SYMBOL_COUNT = 32;
const int MAX_IMAGE_SIZE = 1000 * 1000;
const int SIGNAL_COUNT = 20;
const int NEURAL_COUNT = 80;
const double LARGE_RAND = 1000000;
const double TEACH_FACTOR = 0.1;
const double THRESHOLD_RECOGNIZE_VALUE = 0.7;
const double THRESHOLD_TEACH_VALUE = 0.95;
const double TEACH_INDEX = 0.92;
const string RESULT_FOLDER = "result";
const string TEACH_FOLDER = "teach/";
const string RECOGNITION_FOLDER = "recognition/";
const string PATH_LIST_TEACH = "result/pathListTeach.txt";
const string PATH_LIST_RECOGNITION = "result/pathListRecognition.txt";
const string TEACH_MODE = "teach/";
const string RECOGNITION_MODE = "recognition/";
const string IMAGE_FOLDER = "result/images/";
const string PATHES_FOLDER = "result/pathes/";
const string FFT_FOLDER = "result/fft/";
const string MATRIX_FOLDER = "result/matrix/";

class NeuralWeb {
public:
	vector<double> _firstLayer;
	vector<double> _secondLayer;
	vector<double> _standard;
	unsigned _signalCount;
	unsigned _symbolCount;
	unsigned _neuralCount;
	double bettaParam;
	double thresholdRecognizeValue;
	double thresholdTeachValue;
	double teachFactor;
	double teachIndex;
	NeuralWeb(unsigned signalCount = SIGNAL_COUNT * 4, unsigned symbolCount = SYMBOL_COUNT, unsigned neuralCount = NEURAL_COUNT) : _signalCount(signalCount), _symbolCount(symbolCount), _neuralCount(neuralCount) {
		for(int i = 0; i < signalCount * symbolCount; ++i)
			_firstLayer.push_back(rand() / LARGE_RAND);
		for(int i = 0; i < neuralCount * symbolCount; ++i)
			_secondLayer.push_back(rand() / LARGE_RAND);
		for(int i = 0; i < symbolCount; ++i)
			_standard.push_back(NULL);
		bettaParam = 1;
		thresholdRecognizeValue = THRESHOLD_RECOGNIZE_VALUE;
		thresholdTeachValue = THRESHOLD_TEACH_VALUE;
		teachFactor = TEACH_FACTOR;
		teachIndex = TEACH_INDEX;
	}

	void SaveMatrix() {
		ofstream saveStream;
		saveStream.open(MATRIX_FOLDER + "firstLayer.txt");
		for(int i = 0; i < _firstLayer.size(); ++i)
			saveStream << _firstLayer[i] << endl;
		saveStream.close();
		saveStream.open(MATRIX_FOLDER + "secondLayer.txt");
		for(int i = 0; i < _secondLayer.size(); ++i)
			saveStream << _secondLayer[i] << endl;
		saveStream.close();
	}

	void LoadMatrix() {
		ifstream loadStream;
		loadStream.open(MATRIX_FOLDER + "firstLayer.txt");
		for(int i = 0; i < _firstLayer.size(); ++i)
			loadStream >> _firstLayer[i];
		loadStream.close();
		loadStream.open(MATRIX_FOLDER + "secondLayer.txt");
		for(int i = 0; i < _secondLayer.size(); ++i)
			loadStream >> _secondLayer[i];
		loadStream.close();
	}
	
	double SigmoidalFunction(double x) {
		return 1 / (1 + exp(bettaParam * x));	
	}

	double DiffSigmoidalFunction(double x) {
		return bettaParam * SigmoidalFunction(x) * (1 - SigmoidalFunction(x));
	}

	double& GetNeural(vector<double> &layer, int layerWidth, int i, int j) {																	//ересь какая-то, без дополнительной переменной вернуть элемент нельзя
		return layer[i + j * layerWidth];
	}

	vector<double> NeuralWebFunction(vector<double> &signalArray, vector<double> &layer, int layerWidth) {
		vector<double> result;
		double temp;
		for (int i = 0; i < layer.size() / layerWidth; ++i) {
			temp = 0;
			for (int j = 0; j < signalArray.size(); ++j) {
				temp += GetNeural(layer, layerWidth,j,i) * signalArray[j]; 
			}
			result.push_back(SigmoidalFunction(temp));
		}
		return result;
	}

	vector<double> DiffNeuralWebFunction(vector<double> &signalArray, vector<double> &layer, int layerWidth) {
		vector<double> result;
		double temp;
		for (int i = 0; i < layer.size() / layerWidth; ++i) {
			temp = 0;
			for (int j = 0; j < signalArray.size(); ++j)
				temp += GetNeural(layer, layerWidth,j,i) * signalArray[j]; 
			result.push_back(DiffSigmoidalFunction(temp));
		}
		return result;
	}

	bool Teach(int currentSymbol, vector<double> signalArray, vector<double> fftArray) {
		vector<double> temp;
		vector<double> IntermediateArray;
		vector<double> diffFirstIntermediateArray;
		vector<double> diffSecondIntermediateArray;
		int resultSymbolNumber;
		bool currentResult;
		bool commonResult = true;
		for (int i = 0; i < _symbolCount; ++i) {
			if(i == currentSymbol) _standard[i] = 1;
			else _standard[i] = 0;
		}
		for(;;) {
			signalArray = fftArray;
			Recognize(signalArray, resultSymbolNumber, IntermediateArray);
			for (int i = 0; i < _symbolCount; ++i)  
				if (fabs(signalArray[i] - _standard[i]) > 1 - thresholdTeachValue) {
					currentResult = false;
					break;
				}
				else currentResult = true;
				if (currentResult == true) return true;
				diffFirstIntermediateArray = DiffNeuralWebFunction(fftArray, _firstLayer, _signalCount);
				for (int i = 0; i < _symbolCount; ++i)
					for (int j = 0; j < _signalCount; ++j)
						GetNeural(_firstLayer, _signalCount, j, i) += teachFactor * (signalArray[i] - _standard[i]) * diffFirstIntermediateArray[i] * fftArray[j];
				//for (int j = 0; j < _signalCount; ++j)
				//	for (int i = 0; i < _neuralCount; ++i)
				//			for (int k = 0; k < _symbolCount; ++k)
				//			GetNeural(_firstLayer, _signalCount, j, i) += teachFactor * (signalArray[k] - _standard[k]) * diffSecondIntermediateArray[k] * GetNeural(_secondLayer, _neuralCount, i, k) * diffFirstIntermediateArray[i] * fftArray[j];	
		}
	}

	bool Recognize(vector<double> &signalArray, int &resultSymbolNumber, vector<double> &IntermediateArray) {
		double maxResultValue = 0;
		bool result = false;
		resultSymbolNumber = -1;
		signalArray = NeuralWebFunction(signalArray, _firstLayer, _signalCount);
		for(int i = 0; i < signalArray.size(); ++i) {
			if (maxResultValue < signalArray[i]) {
				maxResultValue = signalArray[i];
				if (maxResultValue >= thresholdRecognizeValue) {
					result = true;
					resultSymbolNumber = i;
				}
			}
		}
		return result;
	}
};

int TransformLexicalToNumericValue (string symbol) {
	int numericValue;
	if(symbol == "а") numericValue = 0; else
	if(symbol == "б") numericValue = 1; else
	if(symbol == "в") numericValue = 2; else
	if(symbol == "г") numericValue = 3; else
	if(symbol == "д") numericValue = 4; else
	if(symbol == "е") numericValue = 5; else
	if(symbol == "ж") numericValue = 6; else
	if(symbol == "з") numericValue = 7; else
	if(symbol == "и") numericValue = 8; else
	if(symbol == "й") numericValue = 9; else
	if(symbol == "к") numericValue = 10; else
	if(symbol == "л") numericValue = 11; else
	if(symbol == "м") numericValue = 12; else
	if(symbol == "н") numericValue = 13; else
	if(symbol == "о") numericValue = 14; else
	if(symbol == "п") numericValue = 15; else
	if(symbol == "р") numericValue = 16; else
	if(symbol == "с") numericValue = 17; else
	if(symbol == "т") numericValue = 18; else
	if(symbol == "у") numericValue = 19; else
	if(symbol == "ф") numericValue = 20; else
	if(symbol == "х") numericValue = 21; else
	if(symbol == "ц") numericValue = 22; else
	if(symbol == "ч") numericValue = 23; else
	if(symbol == "ш") numericValue = 24; else
	if(symbol == "щ") numericValue = 25; else
	if(symbol == "ъ") numericValue = 26; else
	if(symbol == "ы") numericValue = 27; else
	if(symbol == "ь") numericValue = 28; else
	if(symbol == "э") numericValue = 29; else
	if(symbol == "ю") numericValue = 30; else
	if(symbol == "я") numericValue = 31; else
	return -1;
	return numericValue;
}

string TransformNumericToLexicalValue (int numericValue) {
	string symbol;
	if(numericValue == 0) symbol="а"; else
	if(numericValue == 1) symbol="б"; else
	if(numericValue == 2) symbol="в"; else
	if(numericValue == 3) symbol="г"; else
	if(numericValue == 4) symbol="д"; else
	if(numericValue == 5) symbol="е"; else
	if(numericValue == 6) symbol="ж"; else
	if(numericValue == 7) symbol="з"; else
	if(numericValue == 8) symbol="и"; else
	if(numericValue == 9) symbol="й"; else
	if( numericValue == 10) symbol ="к"; else
	if( numericValue == 11) symbol ="л"; else
	if( numericValue == 12) symbol ="м"; else
	if( numericValue == 13) symbol ="н"; else
	if( numericValue == 14) symbol ="о"; else
	if( numericValue == 15) symbol ="п"; else
	if( numericValue == 16) symbol ="р"; else
	if( numericValue == 17) symbol ="с"; else
	if( numericValue == 18) symbol ="т"; else
	if( numericValue == 19) symbol ="у"; else
	if( numericValue == 20) symbol ="ф"; else
	if( numericValue == 21) symbol ="х"; else
	if( numericValue == 22) symbol ="ц"; else
	if( numericValue == 23) symbol ="ч"; else
	if( numericValue == 24) symbol ="ш"; else
	if( numericValue == 25) symbol ="щ"; else
	if( numericValue == 26) symbol ="ъ"; else
	if( numericValue == 27) symbol ="ы"; else
	if( numericValue == 28) symbol ="ь"; else
	if( numericValue == 29) symbol ="э"; else
	if( numericValue == 30) symbol ="ю"; else
	if( numericValue == 31) symbol ="я"; else
	return symbol;
}

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
		else break;
	delete[] blackAndWhite;
	delete[] scanArray;
	}
	return iterator;
}

vector<double> FFT(vector<double> dataDouble, int size) {
	ap::complex_1d_array complexData;
	ap::complex_1d_array fftData;
	vector<double> fftDoubleData;
	complexData.setbounds(0, size);
	fftData.setbounds(0, size);
	for (int i = 0; i < size; ++i) {
		complexData(i).x = dataDouble[2 * i];
		complexData(i).y = dataDouble[2 * i + 1];
	}
	for (int i = 1; i < size; ++i)
		for (int j = 0; j < size; ++j) {
			fftData(i).x += complexData(j).x * cos(2 * ap::pi() * i * j / size) + complexData(j).y * sin(2 * ap::pi() * i * j / size);	
			fftData(i).y += complexData(j).y * cos(2 * ap::pi() * i * j / size) - complexData(j).x * sin(2 * ap::pi() * i * j / size);
		}
	fftData(0) = 0;
	double scalingFactor = sqrt(pow(ap::abscomplex(fftData(1)), 2) + pow(ap::abscomplex(fftData(size - 1)), 2));
	for (int i = 0; i < size; ++i) {
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
				for (int i = whitePixelCount - 2 * SIGNAL_COUNT; i < whitePixelCount; ++i)		// В обратную сторону сигналы формируются начиная с последнего элемента ряда Фурье
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
			while (it->path().filename().string()[iterator] < '0' || it->path().filename().string()[iterator] > '9') {
				if (it->path().filename().string()[iterator] == '"') continue;
				string temp = string(&it->path().filename().string()[iterator],0,1);
				symbolValue += temp;
				++iterator;
			}
			iterator = 0;
			while (iterator < it->path().string().length()) {
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

void DeleteFiles(string fftDir, string imageDir, string pathesDir, string matrixDir, string mode) {
	namespace fs = boost::filesystem;
	typedef fs::recursive_directory_iterator rdi;
	int iterator;
	string filePath;
    rdi itBegFft(fftDir + mode), itEndFft; 
    for(; itBegFft != itEndFft; ++itBegFft) {
		string filePath;
		iterator = 0;
		if (itBegFft->path().extension() == ".bmp" || itBegFft->path().extension() == ".txt" ) {
			while (iterator < itBegFft->path().string().length()) {
					string temp = string(&itBegFft->path().string()[iterator],0,1);
					filePath += temp;
					++iterator;
				}
			std::remove(filePath.c_str());
		}
    }
	rdi itBegImage(imageDir + mode), itEndImage; 
    for(; itBegImage != itEndImage; ++itBegImage) {
		string filePath;
		iterator = 0;
		if (itBegImage->path().extension() == ".bmp" || itBegImage->path().extension() == ".txt" ) {
			while (iterator < itBegImage->path().string().length()) {
					string temp = string(&itBegImage->path().string()[iterator],0,1);
					filePath += temp;
					++iterator;
				}
			std::remove(filePath.c_str());
		}
    }
	rdi itBegPathes(pathesDir + mode), itEndPathes; 
    for(; itBegPathes != itEndImage; ++itBegPathes) {
		string filePath;
		iterator = 0;
		if (itBegPathes->path().extension() == ".bmp" || itBegPathes->path().extension() == ".txt" ) {
			while (iterator < itBegPathes->path().string().length()) {
					string temp = string(&itBegPathes->path().string()[iterator],0,1);
					filePath += temp;
					++iterator;
				}
			std::remove(filePath.c_str());
		}
    }
	if (mode == TEACH_MODE)
	{
		rdi itBegMatrix(matrixDir), itEndMatrix; 
		for(; itBegMatrix != itEndMatrix; ++itBegMatrix) {
			string filePath;
			iterator = 0;
			if (itBegMatrix->path().extension() == ".bmp" || itBegMatrix->path().extension() == ".txt" ) {
				while (iterator < itBegMatrix->path().string().length()) {
						string temp = string(&itBegMatrix->path().string()[iterator],0,1);
						filePath += temp;
						++iterator;
					}
				std::remove(filePath.c_str());
			}
		}
	}
}

bool NeuralWebTeach(string fftFolder, string pathListTeach) {
	namespace fs = boost::filesystem;
	NeuralWeb ThreeLayerPerceptron;
	int fileCount = 0;
	int recognizeCount = 0;
	int resultSymbolNumber;
	int factSymbolNumber;
	vector<double> IntermediateArray;
	string factSymbolValue;
	string currentLexicalSymbol;
	vector<double> signals;
	vector<double> fftSignals;
	string fftSymbolPath;
	for (fs::directory_iterator it(fftFolder + "teach/"), end; it != end; ++it)
			if (it->path().extension() == ".txt") 
				++fileCount;
	ThreeLayerPerceptron.LoadMatrix();
	ifstream factSymbolStream;
	for(;;) {
		factSymbolStream.open(pathListTeach);
		cout << recognizeCount << endl;
		recognizeCount = 0;
		for (fs::directory_iterator it(fftFolder + "teach/"), end; it != end; ++it)
			if (it->path().extension() == ".txt") {
				int iterator = 0;
				int tempSize = 0;
				fftSymbolPath = "";
				while(it->path().string()[iterator]) {
					if (it->path().string()[iterator] == '"') continue;
					string temp = string(&it->path().string()[iterator],0,1);
					fftSymbolPath += temp;
					++iterator;
				}
				signals = ReadFile(fftSymbolPath, tempSize);
				fftSignals = signals;
				getline(factSymbolStream, factSymbolValue);
				getline(factSymbolStream, factSymbolValue);
				if (factSymbolValue.length() > 0)
				factSymbolNumber = TransformLexicalToNumericValue(factSymbolValue);
				else continue;
				cout << "Training on the symbol '" << factSymbolValue << "'\r";
				ThreeLayerPerceptron.Recognize(signals, resultSymbolNumber, IntermediateArray);
				if (resultSymbolNumber == factSymbolNumber)
					++recognizeCount; 
				else ThreeLayerPerceptron.Teach(factSymbolNumber,signals, fftSignals);
				if (recognizeCount >= ThreeLayerPerceptron.teachIndex * fileCount) {
					ThreeLayerPerceptron.SaveMatrix();
					return true;
				}
			}
			factSymbolStream.close();
			if (recognizeCount >= ThreeLayerPerceptron.teachIndex * fileCount) {
				ThreeLayerPerceptron.SaveMatrix();
				return true;
			}
			else continue;
	}
}

bool NeuralWebRecognition(string fftFolder, string pathListTeach) {
	namespace fs = boost::filesystem;
	NeuralWeb ThreeLayerPerceptron;
	int resultSymbolNumber;
	vector<double> IntermediateArray;
	string currentLexicalSymbol;
	vector<double> signals;
	vector<double> fftSignals;
	string fftSymbolPath;
	bool result = false;
	ThreeLayerPerceptron.LoadMatrix();
	for (fs::directory_iterator it(fftFolder + RECOGNITION_MODE), end; it != end; ++it)
		if (it->path().extension() == ".txt") {
			int iterator = 0;
			int tempSize = 0;
			fftSymbolPath = "";
			while(it->path().string()[iterator]) {
				if (it->path().string()[iterator] == '"') continue;
				string temp = string(&it->path().string()[iterator],0,1);
				fftSymbolPath += temp;
				++iterator;
			}
			signals = ReadFile(fftSymbolPath, tempSize);
			fftSignals = signals;
		cout << "Result symbol: ";
		if (ThreeLayerPerceptron.Recognize(signals, resultSymbolNumber, IntermediateArray)) {
			cout << TransformNumericToLexicalValue(resultSymbolNumber) << endl;
			//cout << resultSymbolNumber << endl;
			result = true;
		}
		else {result = false; cout<< "Bad recognize" << endl; continue;}
	}
	return result;
}

int _tmain(int argc, _TCHAR* argv[])
{
	setlocale(LC_ALL, "Russian");
	string mode;
	char answer = 'y';
	int imageCount;
	double start_time, end_time;
	double timer;
	for(;;) {
	if(answer == 'y') {
		cout << "Choose a mode (teach/recognition): ";
		cin >> mode;
		if(mode == "teach") { 
			DeleteFiles(FFT_FOLDER, IMAGE_FOLDER, PATHES_FOLDER, MATRIX_FOLDER, TEACH_MODE);
			cout << "You choose Teach mode" << endl; 
			cout << "Read folder..." << endl;
			ReadFolder(TEACH_FOLDER, PATH_LIST_TEACH);
			cout << "Image processing... ";
			imageCount = ImageProcessing(PATH_LIST_TEACH, TEACH_MODE);
			cout << "ready" <<endl << "Processed " << imageCount << " images." << endl << "Handling coordinates...";
			FormSignals(PATHES_FOLDER + TEACH_FOLDER, FFT_FOLDER, TEACH_MODE);
			cout << "ready" << endl << "For training used " << 2 * SIGNAL_COUNT << " signals." <<endl;
			cout << "Starting a learning process: " << endl;
			start_time = omp_get_wtime();
			if(NeuralWebTeach(FFT_FOLDER, PATH_LIST_TEACH))
				cout << "\r		Training is complete!" << endl;
			end_time = omp_get_wtime();
			timer = end_time - start_time;
			wcout << L"Время обучения: " << timer << L" секунд.\n\n";
		}
		else if(mode == "recognition") {
			DeleteFiles(FFT_FOLDER, IMAGE_FOLDER, PATHES_FOLDER, MATRIX_FOLDER, RECOGNITION_MODE);
			cout << "You choose recognition mode" << endl; 
			cout << "Read folder..." << endl;
			ReadFolder(RECOGNITION_FOLDER, PATH_LIST_RECOGNITION);
			cout << "Image processing... ";
			imageCount = ImageProcessing(PATH_LIST_RECOGNITION, RECOGNITION_MODE);
			cout << "ready" <<endl << "Processed " << imageCount << " images." << endl << "Handling coordinates...";
			FormSignals(PATHES_FOLDER + RECOGNITION_FOLDER, FFT_FOLDER, RECOGNITION_MODE);
			cout << "ready" << endl << "For recognition used " << 2 * SIGNAL_COUNT << " signals." <<endl;
			cout << "Starting a recognition process: " << endl;
			start_time = omp_get_wtime();
			if(NeuralWebRecognition(FFT_FOLDER, PATH_LIST_RECOGNITION))
				cout << "\r		Recognition is complete!" << endl;
			end_time = omp_get_wtime();
			timer = end_time - start_time;
			wcout << L"Время распознавания: " << timer << L" секунд.\n\n";
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

