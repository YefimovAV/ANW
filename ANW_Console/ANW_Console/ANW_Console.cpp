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
const int SYMBOL_COUNT = 33;
const int MAX_IMAGE_SIZE = 1000 * 1000;
const int SIGNAL_COUNT = 20;
const int NEURAL_COUNT = 80;
const double LARGE_RAND = 1000000;
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
	vector<double> _thirdLayer;
	vector<double> _standard;
	unsigned _signalCount;
	unsigned _symbolCount;
	unsigned _neuralCount;
	double bettaParam;
	double thresholdRecognizeValue;
	double thresholdTeachValue;
	NeuralWeb(unsigned signalCount = SIGNAL_COUNT * 2, unsigned symbolCount = SYMBOL_COUNT, unsigned neuralCount = NEURAL_COUNT) : _signalCount(signalCount), _symbolCount(symbolCount), _neuralCount(neuralCount) {
		for(int i = 0; i < signalCount * neuralCount; ++i)
			_firstLayer.push_back(rand() / LARGE_RAND);
		for(int i = 0; i < neuralCount * neuralCount; ++i)
			_secondLayer.push_back(rand() / LARGE_RAND);
		for(int i = 0; i < neuralCount * symbolCount; ++i)
			_thirdLayer.push_back(rand() / LARGE_RAND);
		for(int i = 0; i < symbolCount; ++i)
			_standard.push_back(NULL);
		bettaParam = 1;
		thresholdRecognizeValue = 0.7;
		thresholdTeachValue = 0.000001;
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
		saveStream.open(MATRIX_FOLDER + "thirdLayer.txt");
		for(int i = 0; i < _thirdLayer.size(); ++i)
			saveStream << _thirdLayer[i] << endl;
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
		loadStream.open(MATRIX_FOLDER + "thirdLayer.txt");
		for(int i = 0; i < _thirdLayer.size(); ++i)
			loadStream >> _thirdLayer[i];
		loadStream.close();
	}
	
	double sigmoidalFunction(double x) {
		return 1 / (1 + exp(bettaParam * x));	
	}

	double diffSigmoidalFunction(double x) {
		return bettaParam * sigmoidalFunction(x) * (1 - sigmoidalFunction(x));
	}

	double& getNeural(vector<double> layer, int layerWidth, int i, int j) {
		return layer[i + j * layerWidth];
	}

	vector<double> NeuralWebFunction(vector<double> signalArray, vector<double> layer, int layerWidth) {
		vector<double> result;
		double temp;
		for (int i = 0; layer.size() / layerWidth; ++i) {
			temp = 0;
			for (int j = 0; signalArray.size(); ++j)
				temp += getNeural(layer, layerWidth,j,i) * signalArray[j]; 
			result.push_back(sigmoidalFunction(temp));
		}
	}

	bool teach(int currentSymbol, vector<double> signalArray) {
		vector<double> temp;
		int resultSymbolNumber;
		bool currentResult;
		bool commonResult = true;
		for (int i = 0; i < _symbolCount; ++i) {
			if(i == currentSymbol) _standard[i] = 1;
			else _standard[i] = 0;
		}
		//recognize(signalArray, resultSymbolNumber);
		//if (resultSymbolNumber == currentSymbol) currentResult = true;
		//else {
		//	currentResult = false;
		//	commonResult = currentResult;
		//}
	}

	bool recognize(vector<double> &signalArray, int &resultSymbolNumber) {
		vector<double> temp;
		double maxResultValue = 0;
		bool result = false;
		resultSymbolNumber = -1;
		temp = NeuralWebFunction(signalArray, _firstLayer, _signalCount);
		temp = NeuralWebFunction(temp, _secondLayer, _neuralCount);
		temp = NeuralWebFunction(temp, _thirdLayer, _symbolCount);
		for(int i = 0; i < temp.size(); ++i) {
			if (maxResultValue < temp[i] && maxResultValue >= thresholdRecognizeValue) {
				maxResultValue = temp[i];
				result = true;
				resultSymbolNumber = i;
			}
		}
		return result;
	}
};

int TransformLexicalToNumericValue (string symbol) {
	int numericValue;
	if(symbol == "a") numericValue = 0; else
	if(symbol == "b") numericValue = 1; else
	if(symbol == "v") numericValue = 2; else
	if(symbol == "g") numericValue = 3; else
	if(symbol == "d") numericValue = 4; else
	if(symbol == "e") numericValue = 5; else
	if(symbol == "j") numericValue = 6; else
	if(symbol == "z") numericValue = 7; else
	if(symbol == "i") numericValue = 8; else
	if(symbol == "ij") numericValue = 9; else
	if(symbol == "k") numericValue = 10; else
	if(symbol == "l") numericValue = 11; else
	if(symbol == "m") numericValue = 12; else
	if(symbol == "n") numericValue = 13; else
	if(symbol == "o") numericValue = 14; else
	if(symbol == "p") numericValue = 15; else
	if(symbol == "r") numericValue = 16; else
	if(symbol == "s") numericValue = 17; else
	if(symbol == "t") numericValue = 18; else
	if(symbol == "y") numericValue = 19; else
	if(symbol == "f") numericValue = 20; else
	if(symbol == "x") numericValue = 21; else
	if(symbol == "c") numericValue = 22; else
	if(symbol == "ch") numericValue = 23; else
	if(symbol == "sh") numericValue = 24; else
	if(symbol == "sch") numericValue = 25; else
	if(symbol == "tz") numericValue = 26; else
	if(symbol == "bI") numericValue = 27; else
	if(symbol == "mz") numericValue = 28; else
	if(symbol == "ea") numericValue = 29; else
	if(symbol == "u") numericValue = 30; else
	if(symbol == "ya") numericValue = 31; else
	return -1;
	return numericValue;
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

bool NeuralWebTeach(string fftFolder, string pathListTeach) {
	namespace fs = boost::filesystem;
	NeuralWeb ThreeLayerPerceptron;
	int fileCount = 0;
	int recognizeCount = 0;
	int resultSymbolNumber;
	int factSymbolNumber;
	string factSymbolValue;
	string currentLexicalSymbol;
	vector<double> signals;
	string fftSymbolPath;
	for (fs::directory_iterator it(fftFolder), end; it != end; ++it)
			if (it->path().extension() == ".txt") 
				++fileCount;
	for(;;) {		
		for (fs::directory_iterator it(fftFolder), end; it != end; ++it)
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
			ifstream factSymbolStream;
			factSymbolStream.open(pathListTeach);
			getline(factSymbolStream, factSymbolValue);
			getline(factSymbolStream, factSymbolValue);
			factSymbolNumber = boost::lexical_cast<double>(factSymbolValue);
			if (ThreeLayerPerceptron.recognize(signals, resultSymbolNumber))
				if (resultSymbolNumber == factSymbolNumber)
					++recognizeCount;
			//NeuralWeb.teach function here;
			}
			if (recognizeCount == fileCount)
				return true;
			else continue;
	}
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
			NeuralWeb NW;
			NW.SaveMatrix();
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

