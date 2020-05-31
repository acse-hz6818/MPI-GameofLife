import os
import numpy as np


def merge(path:str,n:int,m:int,step:int):
    totalmap=[]
    for i in range(n):
        matline = []
        for j in range(m):
            name=path+str(step)+" step"+"("+str(i)+","+str(j)+").txt"
            with open(name,"r") as fp:
                mat = []
                sz = fp.readline()
                sz = sz.replace('\n','')
                sz = sz.split(',')
                coordi = int(sz[0])
                coordj = int(sz[1])
                data = fp.read()
                if data[-1]=='\n':
                    data = data[:-1]
                data = data.split('\n')
                if data[-1]=='':
                    data = data[:-1]
                
                for strline in data:
                    tmp = strline.split(" ")
                    if tmp[-1]=='':
                        tmp=tmp[:-1]
                    mat.append(list(map(int,tmp)))
            matline.append(mat)
        totalmap.append(matline)
    

    for i in range(0,n):
        for j in range(0,m):
            if j==0:
                mergedline = totalmap[i][j]
            else:
                mergedline = np.append(mergedline,totalmap[i][j],axis=1)
        if i==0:
            merge = mergedline
        else:
            merge = np.append(merge,mergedline,axis=0)

    with open(path+str(step)+".txt",'w') as fp:
        for i in merge:
            for j in i:
                fp.write(str(j))
                fp.write(' ')
            fp.write("\n")
    


                


def main():
    n = 2
    m = 2
    step = 4
    BASEDIR = os.getcwd()
    BASEDIR=BASEDIR+"/"
    for i in range(step+1):
        merge(BASEDIR,n,m,i)

main()
