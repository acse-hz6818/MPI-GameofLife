#include<iostream>
#include<cmath>
#include<fstream>
#include<ctime>
#include<cstring>
#include<string>
#include<algorithm>
#include "mpi.h"
using namespace std;

void generate_by_rand(bool **data);
bool ** create_2d_bool_arr(int nrows,int ncols);
void destory_2d_bool_arr(bool **data,int nrows);
void print_mat_tofile(bool **data,int *coord,int step);
void update(bool ** data);



void rankToCoord(int *numCPU,int rank,int *cpuCoord);
int coordToRank(int *numCPU,int x,int y);
int decomp(int n,int numCPU,int coord);
void communicate(bool **data);
void decompCPU(int totalCPU,int n,int m,int *bestCPU);



//number of rows
static int n = 20;
//number of cols
static int m = 20;
//number of iterations
static int num_iter = 20;
static bool periodic = true;
static int neigh_x[8] = {-1,-1,-1,0,0,1,1,1};
static int neigh_y[8] = {-1,0,1,-1,1,-1,0,1};

static int myid, p;
static int local_n,local_m;
static int numCPU[2],cpuCoord[2];
static int up,down,lft,rht;
static bool sendbuf[1024];
static bool recvbuf[1024];




int main(int argc,char ** argv){
    clock_t startTime, endTime;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&p);

    numCPU[0]=p;
    numCPU[1]=1;

    decompCPU(p,n,m,numCPU);


    rankToCoord(numCPU,myid,cpuCoord);
    
    local_n=decomp(n,numCPU[0],cpuCoord[0]);
    local_m=decomp(m,numCPU[1],cpuCoord[1]);

    //neighbor rank
    up=coordToRank(numCPU,cpuCoord[0]-1,cpuCoord[1]);
    down=coordToRank(numCPU,cpuCoord[0]+1,cpuCoord[1]);
    lft=coordToRank(numCPU,cpuCoord[0],cpuCoord[1]-1);
    rht=coordToRank(numCPU,cpuCoord[0],cpuCoord[1]+1);

    //halo on the boundary
    bool** data = create_2d_bool_arr(local_n+2,local_m+2);

    //generate_by_rand(data);
    generate_by_rand(data);
    print_mat_tofile(data,cpuCoord,0);


    startTime = clock();

    for(int k=0;k<num_iter;k++){
        communicate(data);
        update(data);
       print_mat_tofile(data,cpuCoord,k+1);

    }

    endTime = clock();

    double time_elasped = (double)(endTime - startTime) / CLOCKS_PER_SEC;

    double recv_time;

    MPI_Reduce(&time_elasped, &recv_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(myid == 0)
        cout << "The maximum time is: " << recv_time << endl << endl;

    // print_mat_tofile(data,cpuCoord,num_iter);

    destory_2d_bool_arr(data,local_n+2);

    MPI_Finalize();
    return 0;
}


// Divide the total processor into the appropriate proportion
void decompCPU(int totalCPU,int n,int m,int *bestCPU)
{
    int cpux,cpuy;

    double diff=9999999;
    double lenx,leny;
    for(int i=1;i<=totalCPU;i++)
    {
        if(totalCPU % i!=0)
            continue;
        cpux=i;
        cpuy=totalCPU/cpux;

        lenx=1.0*n/cpux;
        leny=1.0*m/cpuy;

        if(abs(leny-lenx)<diff)
        {
            diff=abs(leny-lenx);
            bestCPU[0]=cpux;
            bestCPU[1]=cpuy;
        }

    }
}

// transform each "id" into coordinate
void rankToCoord(int *numCPU,int myid,int *cpuCoord)
{
    //(index / the number of col) = x_coordinate
    cpuCoord[0]=myid / numCPU[1];
    //(index % the number of col) = y_coordinate
    cpuCoord[1]=myid % numCPU[1];
}

// calculate the exact size for each grid
int decomp(int n,int numCPU,int coord)
{
    int local = n/numCPU;
    if(coord < n % numCPU)
        local++;
    return local;
}


// transform coordinate into each "id"
int coordToRank(int *numCPU,int x,int y)
{

    bool exchangeFlag=false;
    if(x<0 || y<0 || x>=numCPU[0] || y>=numCPU[1])
        exchangeFlag=true;

    if(x<0)
        x+=numCPU[0];
    if(y<0)
        y+=numCPU[1];
    if(x>=numCPU[0])
        x-=numCPU[0];
    if(y>=numCPU[1])
        y-=numCPU[1];

    if(exchangeFlag)
        return -(x*numCPU[1]+y);
    else
        return x*numCPU[1]+y;
}

// find the state of neighbour and apply the game rules
void update(bool ** data)
{
    bool ** temp = create_2d_bool_arr(local_n+2,local_m+2);
    int num_alive=0;
    int temp_i,temp_j;
    for(int i=1;i<local_n+1;i++)
    {
        for(int j=1;j<local_m+1;j++)
        {
            num_alive=0;
            for(int k=0;k<8;k++)
            {
                temp_i=i+neigh_x[k];
                temp_j=j+neigh_y[k];
                if(data[temp_i][temp_j])
                    num_alive++;
            }

            if(!data[i][j] && num_alive == 3) temp[i][j] = true;
            else if(data[i][j] && (num_alive < 2 || num_alive > 3)) temp[i][j] = false;
            else temp[i][j] = data[i][j];
        }
    }

    for(int i=1;i<local_n+1;i++)
    {
        for(int j=1;j<local_m+1;j++){
            data[i][j] = temp[i][j];
        }
    }


    destory_2d_bool_arr(temp,local_n+2);
}

bool ** create_2d_bool_arr(int nrows,int ncols)
{
    bool ** data=new bool*[nrows];
    for(int i=0;i<nrows;i++)
    {
        data[i]=new bool[ncols];    
        for(int j=0;j<ncols;j++)
            data[i][j]=0;
    }
    return data;
}

void destory_2d_bool_arr(bool **data,int nrows)
{
    for(int i=0;i<nrows;i++)
    {
        delete []data[i];
    }
    delete []data;
}

void generate_by_rand(bool **data){
    //set rand seed
    srand(myid);
    for (int i = 1; i <local_n+1; i++)
    {
        for (int j = 1; j < local_m+1; j++)
            data[i][j]=rand()%2;
    }
}

void print_matrix(bool** data)
{
    for (int i = 1; i <n+1; i++)
    {
        for (int j = 1; j < m+1; j++)
        {
            cout<<data[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void print_mat_tofile(bool **data,int *coord,int step)
{
    string filename;
    filename+=to_string(step);
    filename+=" step";
    filename+="(";
    filename+=to_string(coord[0]);
    filename+=",";
    filename+=to_string(coord[1]);
    filename+=").txt";

    ofstream out(filename);
    out<<local_n<<","<<local_m<<endl;

    for (int i = 1; i <local_n+1; i++)
    {
        for (int j = 1; j < local_m+1; j++)
        {
            out<<static_cast<int>(data[i][j])<<" ";
        }
        out<<endl;
    }
    out<<endl;
    out.close();
}

// communicate with four neighbours around the own processor
void communicate(bool **data)
{
    MPI_Status sta;
    MPI_Request req;

    //send to left recv from right
    //pack
    for(int i=0;i<local_n;i++)
        sendbuf[i]=data[i+1][1];
    MPI_Irecv(recvbuf,local_n,MPI_CXX_BOOL,abs(rht),abs(rht),MPI_COMM_WORLD,&req);
    MPI_Send(sendbuf,local_n,MPI_CXX_BOOL,abs(lft),myid,MPI_COMM_WORLD);
    MPI_Wait(&req,&sta);
    //unpack
    if(periodic || (!periodic && (rht==myid + 1) && (rht>=0)))
    {
        for(int i=0;i<local_n;i++)
            data[i+1][local_m+1]=recvbuf[i];
    }

    //send to right recv from left
    //pack
    for(int i=0;i<local_n;i++)
        sendbuf[i]=data[i+1][local_m];
    MPI_Irecv(recvbuf,local_n,MPI_CXX_BOOL,abs(lft),abs(lft),MPI_COMM_WORLD,&req);
    MPI_Send(sendbuf,local_n,MPI_CXX_BOOL,abs(rht), myid,MPI_COMM_WORLD);
    MPI_Wait(&req,&sta);
    //unpack
    if(periodic || (!periodic && (lft==myid - 1) && (lft>=0)))
    {
        for(int i=0;i<local_n;i++)
            data[i+1][0]=recvbuf[i];
    }


    //send to down recv from up
    //pack
    for(int j=0;j<local_m+2;j++)
        sendbuf[j]=data[local_n][j];
    MPI_Irecv(recvbuf,local_m+2,MPI_CXX_BOOL,abs(up),abs(up),MPI_COMM_WORLD,&req);
    MPI_Send(sendbuf,local_m+2,MPI_CXX_BOOL,abs(down),myid,MPI_COMM_WORLD);
    MPI_Wait(&req,&sta);
    //unpack
    if(periodic || (!periodic && (up==myid - numCPU[1]) && (up>=0)))
    {
        for(int j=0;j<local_m+2;j++)
            data[0][j]=recvbuf[j];
    }

    //send to up recv from down
    //pack
    for(int j=0;j<local_m+2;j++)
        sendbuf[j]=data[1][j];
    MPI_Irecv(recvbuf,local_m+2,MPI_CXX_BOOL,abs(down),abs(down),MPI_COMM_WORLD,&req);
    MPI_Send(sendbuf,local_m+2,MPI_CXX_BOOL,abs(up),myid,MPI_COMM_WORLD);
    MPI_Wait(&req,&sta);
    //unpack
    if(periodic || (!periodic && (down==myid + numCPU[1]) && (down>=0)))
    {
        for(int j=0;j<local_m+2;j++)
            data[local_n+1][j]=recvbuf[j];
    }
 
}

