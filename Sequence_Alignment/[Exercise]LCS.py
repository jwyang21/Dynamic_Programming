#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Implementation of LCS: longest common sequence


# In[2]:


import pandas as pd
import numpy as np


# In[3]:


def lcs(x, y):
    len_x = len(x)
    len_y = len(y)
    S = np.zeros((len_x+1, len_y+1), dtype=int)
    # Initialize scoring matrix (first row and column)
    S[0,:] = np.fromfunction(lambda x,y: 0, (1,len_y+1), dtype=int)
    S[:,0] = np.fromfunction(lambda x,y: 0, (1, len_x+1), dtype=int)
    for i in range(1,len_x+1):
        for j in range(1,len_y+1):
            if x[i-1]==y[j-1]: #match
                S[i][j] = S[i-1][j-1]+1
            else: #mismatch
                S[i][j] = max(S[i-1][j],S[i][j-1])
    return S, S[len_x][len_y]           


# In[4]:


def backtrace(S, x, y):
    len_x = S.shape[0]-1
    len_y = S.shape[1]-1
    len_LCS = S[len_x][len_y]
    LCS = ['']*len_LCS
    index = len_LCS-1
    while len_x>0 and len_y>0:
        if x[len_x-1]==y[len_y-1]:
            LCS[index] = x[len_x-1]
            index-=1
            len_x-=1
            len_y-=1
        else:
            if S[len_x-1][len_y]>=S[len_x][len_y-1]:
                len_x-=1
            else:
                len_y-=1
    return ('').join(LCS)


# In[5]:


x = "ABCBDAB"
y = "BDCABA"
S, len_LCS = lcs(x,y)


# In[6]:


print("LCS: Longest Common Sequence")
print("\n")
print("1. Scoring matrix S: ")
print(S)
print("\n")
print("2. Length of LCS: ", len_LCS)
print("\n")


# In[7]:


print("3. Backtracking...")
print("\n")


# In[8]:


LCS = backtrace(S, x, y)


# In[9]:


print("Longest Common Sequence between x={} and y={} is {}.".format(x, y, LCS))


# In[13]:


get_ipython().system('jupyter nbconvert --to script [Exercise]LCS.ipynb')


# In[ ]:




