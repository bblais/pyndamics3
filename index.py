#!/usr/bin/env python
# coding: utf-8

# In[2]:


#hide
from pyndamics3.core import *
from pyndamics import Simulation


# # Pyndamics3
# 
# > An Update of the Python Numerical Dynamics library pyndamics

# This file will become your README and also the index of your documentation.

# ## Install

# `pip install pyndamics3`

# ## How to use

# Some simple examples.

# In[3]:


sim=Simulation()
sim.add("p'=a*p*(1-p/K)",1,plot=True)
sim.params(a=1,K=50)
sim.run(50)


# In[ ]:




