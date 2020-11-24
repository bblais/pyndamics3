#!/usr/bin/env python
# coding: utf-8

# In[1]:


# default_exp testit


# In[2]:


#hide
import numpy as np


# In[3]:


#hide
def blue():
    print("blue")


# In[4]:


#export
def red():
    print("red")
    
def cyan():
    print("light ",end="")
    blue()


# In[ ]:




