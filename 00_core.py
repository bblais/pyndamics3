#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# default_exp core


# # Simulation
# 
# > The main simulation functions used for adding, running, and visualizing ODEs

# In[1]:


#hide
from nbdev.showdoc import *


# In[2]:


#export
def say_hello(to):
    "Say hello to somebody"
    return f'Hello {to}!'


# Some examples
# 
# $$E=mc^2 \times \sqrt{\frac{1}{1-v^2/c^2}}$$

# In[3]:


say_hello("Me!")


# make sure to include the parameter, or you'll get an error.

# In[4]:


say_hello()


# In[5]:


assert say_hello("Jeremy")=="Hello Jeremy!"


# In[ ]:




