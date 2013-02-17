// K-mer.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"


using namespace std;

// Global variables 
string InputPath  = "";
string OutputPath = "";    // By default output file looks like xxxxx_K5.txt
int K = -1;                // Length of K-mer
int L = -1;                // Level
bool RECURSIVE = false;
bool COMPRESS = false;

// Help string
string HELP(" You should use these arguments:\n-f=Input_File\n-k=Kmer_Length\n[-l=Level; default=0]\n[-o=Output_File_Base; default=InputFile]\n[-r; recursive mode, counts all k-mer from 1 to Kmer_Length]\n[-c; compress revers-complementary k-mers]");

// Decoding command line
bool CheckArguments(int argc, char* argv[]) {
	if ((argc > 7) || (argc < 3)) {
		cout<<"Incorrect number of arguments!"; 
		return false;
	}

	for (int i = 1; i < argc; i++) {
		string line(argv[i]);
		int _F = line.find("-f=") + line.find("-F=");
		int _K = line.find("-k=") + line.find("-K=");
		int _L = line.find("-l=") + line.find("-L=");
		int _O = line.find("-o=") + line.find("-O=");
		int _R = line.find("-r")  + line.find("-R");
		int _C = line.find("-c")  + line.find("-C");
		
		if (_F == -1) {
			if (InputPath != "") {
				cout<<"Redefining of Input File parameter. Incorrect arguments!";
				return false;
			}
			InputPath = line.substr(3,line.length()-3);
		}
		else if (_K == -1) {
			if (K != -1) {
				cout<<"Redefining of Kmer Length parameter. Incorrect arguments!";
				return false;
			}
			K = atoi(line.substr(3,line.length()-3).c_str());
		}
		else if (_L == -1) {
			if (L != -1) {
				cout<<"Redefining of Level parameter. Incorrect arguments!";
				return false;
			}
			L = atoi(line.substr(3,line.length()-3).c_str());
		}
		else if (_O == -1) {
			if (OutputPath != "") {
				cout<<"Redefining of Output File parameter. Incorrect arguments!";
				return false;
			}
			OutputPath = line.substr(3,line.length()-3);
		}
		else if (_R == -1) {
			RECURSIVE = true;
		}
		else if (_C == -1) {
			COMPRESS = true;
		}
		else {
			cout<<"Incorrect arguments!";
			return false;
		}
	}

	// default values if needed
	if (OutputPath == "") OutputPath = InputPath.substr(0, min(InputPath.find_first_of('_'),InputPath.find_first_of('.')));
	if (L == -1) L = 0;

	if (InputPath == "" || K == -1) {
		cout<<"Incorrect arguments!";
		return false;
	}

	return true;
}

// Convert nucleotide letter to a number
int c2i (char c) {
	switch (c) {
		case 'A':
		case 'a':
			return 0;
		case 'G':
		case 'g':
			return 1;
		case 'C':
		case 'c':
			return 2;
		case 'T':
		case 't':
			return 3;
		default:
			return -1;
	}
}

long long GetHash(string s) {
	int h;
	long long hash = 0;

	for (int i = 0; i < K; i++) {
		h = c2i(s[i]);
		if (h == -1) {
			return -i - 1;
		}
		else {
			hash = (hash<<2) + h;
		}
	}
	return hash;
}

// Convert number to nucleotide letter
char i2c(int i) {
	switch(i) {
		case 0:
			return 'A';
		case 1:
			return 'G';
		case 2:
			return 'C';
		case 3:
			return 'T';
	}
}

long long* AllocateMemory(int k) {
	// allocating memory
	int SIZE = 1<<(k*2);
	long long* tmp = new long long[SIZE];

	// filling our array with zeros
	for (int i = 0; i < SIZE; i++) {
		tmp[i] = 0;
	}
	return tmp;
}

string GetSeq (long long hash, int k) {
	string s = "";
	int Mask = 0x03;
	int tmp;

	for (int i = 0; i < k;  i++) {
		tmp = (hash>>(i*2)) & Mask;
		s = i2c(tmp) + s;
	}

	return s;
}

void CountKmers(long long* KmerNumber) {
	cout<<"Start counting...\n";
	ifstream InputFile(InputPath.c_str());
	string s("");
	string ss("");
	string s_("");
	long long hash;
	int line = 0;

	getline(InputFile,s_);
	while (getline(InputFile,s_)) {
		if (s_.find('>') + s_.find('@') != -2) {
			for (int i = 0; i <= s.length() - K; i++) {
				ss = s.substr(i, K);
				hash = GetHash(ss);
				if (hash < 0) {
					i = i- hash - 1;
				}
				else {
					KmerNumber[hash]++;
				}
			}
			s = "";
		}
		else {
			s = s + s_;
		}
		line++;
		if (line % 1000000 == 0) {cout<<line<<" lines done... "<<endl;}
	}
	
	InputFile.close();
	cout<<"Counting done\n";
}

long long GetReverseComplementary(long long num, int k) {
  int* comp_num = new int[k];

  for (int i = 0; i < k; i++) {
    comp_num[i] = 3 - num % 4;
    num = num / 4;
  }

  long long cn = 0;
  for (int i = 0; i < k; i++) {
    cn += comp_num[i] * pow(4.0, k-i-1);
  }
  return cn;
}

void WriteKmers(long long* KmerNumber, string base, int k) {
	cout<<"Start writing...\n";
	string out = base + "_K" + to_string((long long)k) +".txt";
	ofstream OutputFile(out.c_str());

	long long SIZE = 1<<(k*2);

	long long TotalNumber = 0;
	bool* IsNotRC = new bool[SIZE];
	for (long long i = 0; i < SIZE; i++) {
		TotalNumber += KmerNumber[i];
		IsNotRC[i] = true;
	}

	OutputFile<<"Sequence KmerNum\n";
	// find and sum the numbers of the revers-complementary k-mers
	if (COMPRESS) {
		for (long long i = 0; i < SIZE; i++) {
			long long cn = GetReverseComplementary(i, k);
			long long cmin, cmax;
			if (cn != i) {            // if cn == i then k-mer is self-revers-complementary; there is no need to do anything
				if (cn > i) {
					cmin = i;
					cmax = cn;
				}
				else {
					cmin = cn;
					cmax = i;
				}
				if (IsNotRC[i]) OutputFile<<GetSeq(i,k)<<" "<<(((double)(KmerNumber[cmin] + KmerNumber[cmax]))/TotalNumber * SIZE)<<endl;
				IsNotRC[cmax] = false;
			}
			else OutputFile<<GetSeq(i,k)<<" "<<(((double)KmerNumber[i])/TotalNumber * SIZE)<<endl;
		}
	}
	else {
		for (long long i = 0; i < SIZE; i++) {
			OutputFile<<GetSeq(i,k)<<" "<<(((double)KmerNumber[i])/TotalNumber * SIZE)<<endl;
		}
	}

	delete[] IsNotRC;
	OutputFile.close();
	cout<<"Writing done. Total number of nucleotides: "<<TotalNumber<<"\n";
}

void RecursiveCountKmers(long long* KmerNumber, int k) {
	cout<<"Starting recursive count: k = "<<k<<endl;
	
	long long* tmp = AllocateMemory(k);
	int SIZE = 1<<(k*2);

	int G = SIZE;
	int C = 2*G;
	int T = 3*G;

#pragma omp parallel for
	for (int i = 0; i<SIZE; i++) {
		int si = i<<2;
		tmp[i] = KmerNumber[si] + KmerNumber[si+1] + KmerNumber[si+2] + KmerNumber[si+3] + 
			KmerNumber[i] + KmerNumber[G+i] + KmerNumber[C+i] + KmerNumber[T+i];
	}

	WriteKmers(tmp, OutputPath, k);
	delete[] KmerNumber;
	if (k==1) return;
	else RecursiveCountKmers(tmp, k-1);
}



int main(int argc, char* argv[]) {
	time_t t1, t2, t3;
	t1 = time(NULL);

	if (!CheckArguments(argc, argv)) {
		cout<<HELP<<endl<<"Press any key to end this programm";
		cin.get();
		return 0;
	}
	
	// If arguments are correct
	// we can proceed to counting of Kmers
	long long* KmerNumber = AllocateMemory(K);
	CountKmers(KmerNumber);
	WriteKmers(KmerNumber, OutputPath, K);

	t2 = time(NULL);
	if (RECURSIVE) {
		RecursiveCountKmers(KmerNumber, K-1);
		t3 = time(NULL);
		cout<<"First count time: "<<t2-t1<<" seconds\n";
		cout<<"Recursive count time: "<<t3-t2<<" seconds\n";
	}
	t3 = time(NULL);

	cout<<"Total time: "<<t3-t1<<" seconds\n";
	cout<<"End\n";
	return 0;
}