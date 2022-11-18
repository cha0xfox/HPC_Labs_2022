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

void DataDistr(int* arr, int arr_size, int*& pArr, int& pArr_size) {
    int RestElem; //Number of rows, that haven't been distributed yet

    MPI_Bcast(&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    RestElem = arr_size;
    for (int i = 0; i < ProcRank; i++) {
        RestElem -= RestElem / (ProcNum - i);
    }
    pArr_size = RestElem / (ProcNum - ProcRank);

    pArr = new int[pArr_size];

    int* send_cnt = new int[ProcNum];
    int* send_disp = new int[ProcNum];
    int cnt = arr_size / ProcNum;
    RestElem = arr_size;

    send_cnt[0] = cnt;
    send_disp[0] = 0;

    for (int i = 1; i < ProcNum; i++) {
        RestElem -= cnt;
        cnt = RestElem / (ProcNum - i);
        send_cnt[i] = cnt;
        send_disp[i] = send_disp[i - 1] + send_cnt[i - 1];
    }

    MPI_Scatterv(arr, send_cnt, send_disp, MPI_INT, pArr, send_cnt[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);

    delete[] send_cnt;
    delete[] send_disp;

}

void ProccesInitialization(int*& arr, int& arr_size) {
	if (ProcRank == 0) {
		do {
			cout << "\nEnter size of array : ";
			cin >> arr_size;
			if (arr_size < ProcNum) cout << "\n Size of array must be greater than ProcNum\n";
		} while (arr_size < ProcNum);

		arr = new int[arr_size];
		randDataInitialization(arr, arr_size);
	}  
}


void MergeCalculation(int*& pArr, int& pArr_size) {
    
    int curr_ProcNum = ProcNum;
    int step = 1;

    while (curr_ProcNum > 1) {
        curr_ProcNum = curr_ProcNum / 2 + curr_ProcNum % 2;

        if ((ProcRank - step) % (2 * step) == 0) {
            MPI_Send(&pArr_size, 1, MPI_INT, ProcRank - step, 0,  MPI_COMM_WORLD);
            MPI_Send(pArr, pArr_size, MPI_INT, ProcRank - step, 0, MPI_COMM_WORLD);
        }

        if ((ProcRank % (2 * step) == 0) && (ProcNum - ProcRank > step)) {
            MPI_Status status;
            int newSize;

            MPI_Recv(&newSize, 1, MPI_INT, ProcRank + step, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int* temp = new int[pArr_size + newSize];
            MPI_Recv(temp, newSize, MPI_INT, ProcRank + step, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            
            for (int i = 0; i < pArr_size; i++) {
                temp[i + newSize] = pArr[i];
            }

            merge(temp, 0, newSize - 1, newSize + pArr_size - 1);

            pArr = new int[newSize + pArr_size];

            for (int i = 0; i < newSize + pArr_size; i++) {
                pArr[i] = temp[i];
            }           

            pArr_size += newSize;
            
            delete[] temp;
        }

        step *= 2;
    }
    
}

void ResultReceive(int* pArr, int arr_size, int* res) {
    int* recv_cnt = new int[ProcNum];
    int* recv_disp = new int[ProcNum];
    int cnt = arr_size / ProcNum;
    int RestElem = arr_size;

    recv_cnt[0] = cnt;
    recv_disp[0] = 0;

    for (int i = 1; i < ProcNum; i++) {
        RestElem -= cnt;
        cnt = RestElem / (ProcNum - i);
        recv_cnt[i] = cnt;
        recv_disp[i] = recv_cnt[i - 1] + recv_disp[i - 1];
    }

    MPI_Allgatherv(pArr, recv_cnt[ProcRank], MPI_INT,
        res, recv_cnt, recv_disp, MPI_INT, MPI_COMM_WORLD);

    delete[] recv_cnt;
    delete[] recv_disp;
}


void DistrCalculation(int* pArr, int pArr_size) {
    mergeSort(pArr, 0, pArr_size - 1);
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
    DataDistr(arr, arr_size, pArr, pArr_size);
    DistrCalculation(pArr, pArr_size);
    MergeCalculation(pArr, pArr_size);
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