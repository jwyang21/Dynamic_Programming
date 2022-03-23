#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Goal: Implement (1) DP and (2) Greedy Algorithm in "The Manhattan Tourist Problem


# In[1]:


import pandas as pd
import numpy as np


# In[2]:


def dynamic(n, m, h, v): #n: number of rows, m: number of columns, h: horizontal weight matrix, v: vertical weight matrix
    S = np.zeros((n,m), dtype=int) #Make array
    #initialize first row and column
    S[0,:]=np.concatenate(([0],np.cumsum(h[0])))
    S[:,0]=np.concatenate(([0],np.cumsum(v[0])))
    for i in range(1,n):
        for j in range(1,m):
            FromUp=S[i-1][j]+v[j][i-1]
            FromLeft=S[i][j-1]+h[i][j-1]
            S[i][j]=max(FromUp, FromLeft)
    return S, S[n-1][m-1]


# In[3]:


horizontal_weight=[[3,2,4,0],[3,2,4,2],[0,7,3,4],[3,3,0,2],[1,3,2,2]]
vertical_weight=[[1,4,4,5],[0,6,4,6],[2,5,5,8],[4,2,2,5],[3,1,1,3]]


# In[4]:


print("The Manhattan Tourist Problem")


# In[5]:


print("1. Dynamic Programming")


# In[6]:


S_dp, max_dp_score = dynamic(5, 5, horizontal_weight, vertical_weight)


# In[7]:


print("1-1. Scoring matrix S calculated from Dynamic Programming:")
print(S_dp)
print("\n")
print("1-2. Max score calculated from Dynamic Programming:")
print(max_dp_score)
print("------------------------------")


# In[8]:


print("2. Greedy Algorithm")


# In[9]:


def greedy(n, m, h, v): #n: number of rows, m: number of columns, h: horizontal weight matrix, v: vertical weight matrix
    S = np.zeros((n,m), dtype=int) #Make array
    i = 0
    j = 0
    while i<n and j<m:
        if j==m-1:
            #print("Reached rightmost column. Going down.")
            S[i+1][j] = S[i][j]+v[j][i]
            i+=1
            if i==n-1: break
        elif i==n-1:
            #print("Reached bottom line. Going right.")
            S[i][j+1] = S[i][j]+h[i][j]
            j+=1
            if j==m-1: break
        else:
            #print("current (i,j):",(i,j))
            #print("S[i][j]:",S[i][j])
            ToRight = S[i][j]+h[i][j]
            #print("ToRight:",ToRight)
            ToDown = S[i][j]+v[j][i]
            #print("ToDown:",ToDown)
            if ToDown>=ToRight:
                #print("Going Down")
                S[i+1][j] = ToDown
                i+=1
                #print("-----------------")
            else:
                #print("Going Right")
                S[i][j+1] = ToRight
                j+=1            
                #print("-----------------")
    return S, S[n-1][m-1]


# In[10]:


S_greedy, max_greedy_score = greedy(5, 5, horizontal_weight, vertical_weight)


# In[11]:


print("2-1. Scoring matrix S calculated from Greedy Algorithm:")
print(S_greedy)
print("\n")
print("2-2. Max score calculated from Greedy Algorithm:")
print(max_greedy_score)
print("------------------------------")
print("\n")


# In[12]:


print("Max score between DP and Greedy Algorithm:", max(max_greedy_score, max_dp_score))
print("\n")


# In[13]:


if max_dp_score>max_greedy_score:
    print("DP yielded optimal solution.")
elif max_dp_score == max_greedy_score:
    print("Both DP and Greedy Algorithm yielded same optimal solution.")
else:
    print("Greedy Algorithm yielded optimal solution.")


# In[16]:


#shell!jupyter nbconvert --to script Manhattan_Tourist_Problem.ipynb

