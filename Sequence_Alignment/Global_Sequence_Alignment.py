#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Needleman-Wunsch Algorithm (NW)


# In[2]:


import pandas as pd
import numpy as np
import blosum as bl


# In[3]:


def NW_fill(x, y, d, sbt):
    '''
    Fill in scoring matrix S using NW algorithm
    d: gap penalty * (-1)
    sbt: substitution matrix (which is BLOSUM62 in this case)
    x: query sequence
    y: target sequence
    '''    
    S = np.zeros((len(x)+1,len(y)+1),dtype=int)
    S[0,:] = [d*i for i in range(len(y)+1)]
    S[:,0] = [d*i for i in range(len(x)+1)]
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            match = S[i-1][j-1]+sbt[str(x[i-1]+y[j-1])]
            FromUp = S[i-1][j]+d
            FromLeft = S[i][j-1]+d
            S[i][j] = max(match, FromUp, FromLeft)
    return S


# In[4]:


def NW_backtrace(x,y,S,sbt,d):
    '''
    Backtrace to get aligned sequence.
    d: gap penalty * (-1)
    sbt: substitution matrix (which is BLOSUM62 in this case)
    x: query sequence
    y: target sequence
    '''
    alignment = []
    current_x = len(x)
    current_y = len(y)
    while current_x>0 and current_y>0:
        if S[current_x][current_y]==S[current_x-1][current_y-1]+sbt[str(x[current_x-1]+y[current_y-1])]:
            alignment.append(x[current_x-1])
            current_x-=1
            current_y-=1
        elif S[current_x][current_y]==S[current_x-1][current_y]+d:
            current_x-=1
            alignment.append('-')
        else:
            current_y-=1
            alignment.append('-')
    alignment.reverse()
    return ('').join(alignment)


# In[5]:


def NW_pointer(x,y,sbt,d):
    '''
    Make pointer dataframe.
    d: gap penalty * (-1)
    sbt: substitution matrix (which is BLOSUM62 in this case)
    x: query sequence
    y: target sequence
    '''
    #d=-5
    P = np.array([['']*(len(y)+1)]*(len(x)+1),dtype=str)
    #P[0,1:] = [y[i] for i in range(len(y))]
    #P[1:,0] = [x[i] for i in range(len(x))]
    current_x = len(x)
    current_y = len(y)
    while current_x>0 and current_y>0:
        if S[current_x][current_y]==S[current_x-1][current_y-1]+sbt[str(x[current_x-1]+y[current_y-1])]:
            P[current_x][current_y] = '↖'
            current_x-=1
            current_y-=1
        elif S[current_x][current_y]==S[current_x-1][current_y]+d:
            P[current_x][current_y] = '↑'
            current_x-=1
        else:
            P[current_x][current_y] = '←'
            current_y-=1
    P_df = pd.DataFrame(data=P[1:,1:], columns=[y[i] for i in range(len(y))], index=[x[i] for i in range(len(x))])
    return P_df


# In[6]:


#import substitution matrix BLOSUM62 from web
sbt = dict(bl.BLOSUM(62)) #use BLOSUM62 as dictionary


# In[7]:


x = "ALVNKDK"
y = "ASLVNDK"
gap_penalty=5
print("Query sequence: {}".format(x))
print("Target sequence: {}".format(y))
print("Gap penalty: {}".format(gap_penalty))
print("\n")
print("Sequence alignment using Needleman-Wunsch Algorithm(Global Alignment)")
print("-------------------------------------")


# In[8]:


S = NW_fill(x,y,(-1)*gap_penalty, sbt)


# In[9]:


print("1. Scoring matrix: ")
print(S)
print("\n")


# In[10]:


seq = NW_backtrace(x,y,S,sbt,(-1)*gap_penalty)


# In[11]:


print("2. Aligned sequence: {}".format(seq))
print("\n")


# In[12]:


P_df = NW_pointer(x,y,sbt,(-1)*gap_penalty)


# In[13]:


print("3. Pointer matrix: ")
print(P_df)

