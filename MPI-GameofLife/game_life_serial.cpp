#include<iostream>
#include<cmath>
#include<fstream>
#include<ctime>
#include<cstring>
using namespace std;

void generate_by_manual(bool **data);
bool ** create_2d_bool_arr(int n_rows,int n_cols);
void destory_2d_bool_arr(bool **data);
void print_matrix(bool** data);
void update(bool ** data);

//number of rows
static int n = 6;
//number of cols
static int m = 6;
//number of iterations
static int num_iter = 4;
static bool periodic = true;
static int neigh_x[8] = {-1,-1,-1,0,0,1,1,1};
static int neigh_y[8] = {-1,0,1,-1,1,-1,0,1};

int main(){
    bool ** data = create_2d_bool_arr(n,m);
    generate_by_manual(data);
    print_matrix(data);

    for(int k=0;k<num_iter;k++){
        update(data);
    }
    print_matrix(data);
    destory_2d_bool_arr(data);
    return 0;
}

// find the state of neighbour and apply the game rules
void update(bool ** data)
{
    bool ** temp = create_2d_bool_arr(n,m);
    int num_alive=0;
    int temp_i,temp_j;
    for(int i=0;i<n;i++)
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
                }else{
                    if((temp_i<0) || (temp_i>=n) || (temp_j<0) || (temp_j>=m))
                        continue;
                }

                if(data[temp_i][temp_j])
                    num_alive++;
            }

            if(!data[i][j] && num_alive == 3) temp[i][j] = true;
            else if(data[i][j] && (num_alive < 2 || num_alive > 3)) temp[i][j] = false;
            else temp[i][j] = data[i][j];
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            data[i][j] = temp[i][j];
        }
    }


    destory_2d_bool_arr(temp);
}

bool ** create_2d_bool_arr(int n, int m)
{
    bool** data = new bool* [n];
    for(int i = 0; i < n; i++)
    {
        data[i] = new bool [m];
    }
    return data;
}

void destory_2d_bool_arr(bool **data)
{
    for(int i = 0; i < m; i++)
        delete[] data[i];
    delete[] data;
}

void generate_by_manual(bool **data){
    fstream fin;
    fin.open("test2.txt", fstream::in);
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

void print_matrix(bool** data)
{
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
