#include "mpi.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include<fstream>
#include<ctime>
#include<cstring>
#include<string>

using namespace std;


void generate_by_manual(bool **data);
void generate_by_rand(bool **data);
bool ** create_2d_bool_arr(int n_rows,int n_cols);
void destory_2d_bool_arr(bool **data,int n_rows);
void print_matrix(bool** data);
void print_mat_tofile(bool **data, int step);
void update(bool ** data);
void communicate(bool **data);


//number of rows
static int n = 20;
//number of cols
static int m = 20;
//number of iterations
static int num_iter = 4;
static bool periodic = true;
static int neigh_x[8] = {-1,-1,-1,0,0,1,1,1};
static int neigh_y[8] = {-1,0,1,-1,1,-1,0,1};

static bool sendbuf[1024];
static bool recvbuf[1024];

static int myid, p, data_row;
static int tag_num = 1;


int main(int argc, char *argv[]){
    clock_t startTime, endTime;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&p);


    bool** data = create_2d_bool_arr(n,m);
    generate_by_manual(data);
    // generate_by_rand(data);

    // print_matrix(data);

    // divide the whole domain into different strips
    // decide the exact size for each strip
    int* start_row = new int[p];
    int* end_row = new int [p];

    int row_rem = n;
    start_row[0] = 0;

    for(int i = 0; i < p; i++)
    {
        int row_on_proc = row_rem / (p - i);
        end_row[i] = start_row[i] + row_on_proc -1;

        if(i + 1 < p)
            start_row[i + 1] = end_row[i] + 1;
        row_rem -= row_on_proc;
    }

    data_row = end_row[myid] - start_row[myid] + 1;


    bool** data_strip = create_2d_bool_arr((data_row + 2), m);

    for(int i = 0; i < data_row; i++){
        for(int j = 0; j < m; j++)
            data_strip[i+1][j] = data[start_row[myid] + i][j];
    }


    delete[] start_row;
    delete[] end_row;

    startTime = clock();

     for(int k=0; k<num_iter; k++){
        communicate(data_strip);
        update(data_strip);
        print_mat_tofile(data_strip,k+1);
    }
    endTime = clock();

    double time_elasped = (double)(endTime - startTime) / CLOCKS_PER_SEC;

    double recv_time;

    MPI_Reduce(&time_elasped, &recv_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(myid == 0)
        cout << "The maximum time is: " << recv_time << endl << endl;

    // cout << "Processor " << myid << " has results: " << endl ;

    //     for(int i = 1; i < (data_row+1); i++){
    //         for(int j = 0; j < m; j++){
    //             cout<< data_strip[i][j] << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;

    
    // print_mat_tofile(data_strip,num_iter);

    destory_2d_bool_arr(data_strip, data_row+2);

    MPI_Finalize();
    return 0;

}

// find the state of neighbour and apply the game rules
void update(bool ** data)
{
    bool ** temp = create_2d_bool_arr(data_row+2,m);
    int num_alive=0;
    int temp_i,temp_j;
    for(int i=1;i<data_row+1;i++)
    {
        for(int j=0;j<m;j++)
        {
            num_alive=0;
            for(int k=0;k<8;k++)
            {
                temp_i=i+neigh_x[k];
                temp_j=j+neigh_y[k];
                if(periodic)
                {
                   temp_i=(temp_i+n) % n;
                    temp_j=(temp_j+m) % m;
                }
                if(data[temp_i][temp_j])
                    num_alive++;
            }

            if(!data[i][j] && num_alive == 3) temp[i][j] = true;
            else if(data[i][j] && (num_alive < 2 || num_alive > 3)) temp[i][j] = false;
            else temp[i][j] = data[i][j];
        }
    }

    for(int i = 1; i < data_row+1; i++){
        for(int j = 0; j < m; j++){
            data[i][j] = temp[i][j];
        }
    }


    destory_2d_bool_arr(temp,data_row+2);
}


bool ** create_2d_bool_arr(int n, int m)
{
    bool** data = new bool* [n];
    for(int i = 0; i < n; i++)
    {
        data[i] = new bool [m];
        for(int j=0;j<m;j++)
            data[i][j]=0;
    }
    return data;
}

void destory_2d_bool_arr(bool **data,int n_rows)
{
    for(int i = 0; i < n_rows; i++)
        delete[] data[i];
    delete[] data;
}


void print_matrix(bool** data){
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            cout<<data[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}


void generate_by_manual(bool **data){
    fstream fin;
    fin.open("0.txt", fstream::in);
    if(!fin) {
        cout << "failed to open test file" << endl;
        return;
    }
    while (fin.good()){
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
                fin >> data[i][j];
        }
    }
}


void generate_by_rand(bool **data){
    // srand(time(NULL));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++)
            data[i][j] = rand()%2;
    }    
}

// communicate with four neighbours around the own processor
void communicate(bool **data)
{
    MPI_Status sta;
    MPI_Request req;

    int up_id,down_id;
    up_id=myid-1;
    down_id=myid+1;
    if(up_id<0)
        up_id=(p-1);
    if(down_id>=p)
        down_id=0;

    //send to down_id recv from up_id
    //pack
    for(int j=0;j < m;j++)
        sendbuf[j]=data[data_row][j];
    MPI_Irecv(recvbuf, m, MPI_CXX_BOOL, up_id, 1, MPI_COMM_WORLD, &req);
    MPI_Send(sendbuf, m, MPI_CXX_BOOL, down_id, 1, MPI_COMM_WORLD);
    MPI_Wait(&req,&sta);
    //unpack
    if(periodic || (!periodic && (up_id==myid-1)))
    {
        for(int j=0;j<m;j++)
            data[0][j]=recvbuf[j];
    }

    //send to up_id recv from down_id
    //pack
    for(int j=0;j< m;j++)
        sendbuf[j]=data[1][j];
    MPI_Irecv(recvbuf, m, MPI_CXX_BOOL, down_id, 1, MPI_COMM_WORLD, &req);
    MPI_Send(sendbuf, m, MPI_CXX_BOOL, up_id, 1, MPI_COMM_WORLD);
    MPI_Wait(&req,&sta);
    //unpack
    if(periodic || (!periodic && (down_id==myid+1)))
    {
        for(int j = 0;j < m;j++)
            data[data_row+1][j]=recvbuf[j];
    }

}


// print data into a txt file
void print_mat_tofile(bool **data, int step)
{
    string filename;
    filename+=to_string(step);
    filename+=" step";
    filename+="(";
    filename+=to_string(myid);
    filename+=").txt";

    ofstream out(filename);
    out << data_row << "," << m <<endl;

    for (int i = 1; i < data_row+1; i++)
    {
        for (int j = 0; j < m; j++)
        {
            out<<static_cast<int>(data[i][j])<<" ";
        }
        out<<endl;
    }
    out<<endl;
    out.close();
}



