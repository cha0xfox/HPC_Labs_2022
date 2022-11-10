#include <mpi.h>
#include <iostream>
#include <time.h>

using namespace std;

int ProcNum;
int ProcRank;

void merge(int* a, int l, int m, int r) {
    //sizes of left and right array
    int nl = m - l + 1;
    int nr = r - m;
    int* la = new int[nl];
    int* ra = new int[nr];

    //fill left and right sub-arrays
    for (int i = 0; i < nl; i++)
        la[i] = a[l + i];
    for (int i = 0; i < nr; i++)
        ra[i] = a[m + 1 + i];

    //merge temp arrays to real array
    int i = 0, j = 0, k = l;
    while (i < nl && j < nr) {
        if (la[i] <= ra[j]) {
            a[k++] = la[i++];
        }
        else {
            a[k++] = ra[j++];
        }
    }

    //extra element in left and right temp array
    while (i < nl)
        a[k++] = la[i++];
    while (j < nr)
        a[k++] = ra[j++];


}

void mergeSort(int* a, int l, int r) {
    int m;
    if (l < r) {
        m = l + (r - l) / 2;
        mergeSort(a, l, m);
        mergeSort(a, m + 1, r);
        merge(a, l, m, r);
    }

}


void randDataInitialization(int* a, int n) {
	srand(time(0));
	for (int i = 0; i < n; i++) {
		a[i] = rand()%100;
	}
}

// Function for computational process termination
void ProcessTermination(double* pMatrix, double* pVector, double* pResult,
    double* pProcRows, double* pProcResult) {
    if (ProcRank == 0)
        delete[] pMatrix;
    delete[] pVector;
    delete[] pResult;
    delete[] pProcRows;
    delete[] pProcResult;
}



int main(int argc, char* argv[])
{
	int arr_size;
    int pArr_size;// Size of array for every procRank
	int* arr; 
	int* pArr;// array for every procRank

    double pStart, pFinish, pDuration;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

	ProccesInitialization(arr, arr_size);

    pStart = MPI_Wtime();
    /*
    somehow distribute sublists
    */
    pFinish = MPI_Wtime();

    pDuration = pFinish - pStart;


    if (ProcRank == 0) {

        cout << "Start array" << "\n";
        for (int i = 0; i < arr_size; i++)
            cout << arr[i] << " ";

        cout << '\n' << "Finish array\n";

        for (int i = 0; i < pArr_size; i++) {
            cout << pArr[i] << " ";
        }

        cout << '\n' << pArr_size;

        cout << "\nProccess finished" << '\n';
        cout << "Time of parallel execution : " << pDuration << '\n';

    }
   
	MPI_Finalize();
    
}
