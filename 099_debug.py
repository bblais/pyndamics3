#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('pylab', 'inline')


# In[2]:


from pyndamics3 import Simulation
from pyndamics3.fit import fit,Parameter


# In[3]:


import lmfit


# In[4]:


t=array([7,14,21,28,35,42,49,56,63,70,77,84],float)
h=array([17.93,36.36,67.76,98.10,131,169.5,205.5,228.3,247.1,250.5,253.8,254.5])


# In[6]:


sim=Simulation()
sim.add("h'=a",1,plot=True)
sim.add_data(t=t,h=h,plot=True)
sim.params(a=1)
sim.run(0,90)


# In[7]:


results=fit(sim,
   Parameter('a',value=1,min=0),
   Parameter('initial_h',value=10,min=0))

results


# In[ ]:




