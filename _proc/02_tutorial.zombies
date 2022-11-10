#!/usr/bin/env python
# coding: utf-8

# # Introduction to MCMC on Dynamical Systems Using Zombies

# In[1]:


#| hide
#skip
get_ipython().system(' [ -e /content ] && pip install -Uqq pyndamics3 emcee # upgrade pyndamics3 on colab')


# In[2]:


get_ipython().run_line_magic('pylab', 'inline')


# In[3]:


from pyndamics3 import Simulation


# In[4]:


from pyndamics3.mcmc import *


# ## SIR Model

# In[5]:


sim=Simulation()
sim.add("S'=-β*S*I",1,plot=1)
sim.add("I'=β*S*I-ζ*I",.001,plot=1)
sim.add("R'=ζ*I",0,plot=1) 
sim.params(β=5,ζ=1)
sim.run(0,10) 


# ## SEIR Model

# In[6]:


sim=Simulation()
sim.add("S'=-β*S*I",1,plot=1)
sim.add("E'=β*S*I-ζ*E",0,plot=1)
sim.add("I'=ζ*E-α*I",.001,plot=1) 
sim.add("R'=α*I",0,plot=1)  
sim.params(α=.3,β=10,ζ=.5)
sim.run(0,10) 


# ## SZR Model from Munz et al. (2009)
# 
# Notice that no matter what the parameters are changed to, Z (zombies) always win.

# In[7]:


sim=Simulation()
sim.add("S'=Π-β*S*Z-δ*S",500,plot=1)                   #S (Susceptible)
sim.add("Z'=β*S*Z+ζ*R-α*S*Z",.002,plot=1)            #Z (Zombie)
sim.add("R'=δ*S+α*S*Z-ζ*R",1,plot=False)            #R (Removed)
sim.params(α=.005,β=.0095,ζ=.05, δ=.01,Π=0)     #parameters changed to match the Munz et al. (2009) figures
sim.run(0,30)


# ## SEZR Model based on dynamics observed in 'Night of the Living Dead'

# Movie "data" from Night of the Living Dead

# In[8]:


t=array([0,1,1.5,3,4.5,5,5.75,5.9,10])
zombies=array([1,1,3,8,10,20,28,30,40])


# In[9]:


sim=Simulation()
sim.add("S'=-β*S*Z-δ*S",178.5,plot=1)
sim.add("E'=β*S*Z-ζ*E",0,plot=False)
sim.add("Z'=ζ*E-α*S*Z",1,plot=1) 
sim.add("R'=α*S*Z+δ*S",0,plot=False) 
sim.params(α=.0342,β=.0445,ζ=4.63, δ=0.0)
sim.add_data(t=t,Z=zombies,plot=1)
sim.run(0,10)


# MCMC parameter estimation for $\alpha$ (rate of zombies being permanently removed), $\beta$ (rate of susceptibles becoming infected), $\zeta$ (the rate of infected into becoming zombies), and $\delta$ (suicide rate among susceptibles)

# In[18]:


model=MCMCModel(sim,
                α=Uniform(0,.5),
                β=Uniform(0,.5),
                ζ=Uniform(0,10),
                δ=Uniform(0,.01),
               )
                


# In[19]:


number_of_iterations=500 # use 500 or so for the figures below, but for CI timeout reasons I include only 5
model.run_mcmc(number_of_iterations,repeat=3)
model.plot_chains()


# In[20]:


sim.run(0,10)


# In[22]:


Ro=model.eval('β/α')


# In[23]:


model.plot_distributions(Ro)


# In[24]:


model.plot_many(0,13,'Z')


# In[25]:


model.triangle_plot()


# In[26]:


model.plot_distributions()


# ## SEZR Model based on dynamics observed in 'Shaun of the Dead'
# 
# Data from Shaun of the Dead

# In[5]:


t=array([0,3,5,6,8,10,22,22.2,22.5,24,25.5,26,26.5,27.5,27.75,28.5,29,29.5,31.5])
zombies=array([0,1,2,2,3,3,4,6,2,3,5,12,15,25,37,25,65,80,100])


# In[6]:


sim=Simulation()
sim.add("S'=-β*S*Z",508.2,plot=1)
sim.add("E'=β*S*Z-ζ*E",0,plot=0)
sim.add("Z'=ζ*E-α*S*Z",.000347759,plot=1)
sim.add("R'=α*S*Z",0,plot=False)
sim.params(α=2.96e-8,β=0.000808995,ζ=60)
sim.add_data(t=t,Z=zombies,plot=1)
sim.run(0,50)


# In[86]:


model=MCMCModel(sim,
                α=Uniform(0,.01),
                β=Uniform(0,.01),
                ζ=Uniform(0,100),
               )


# In[87]:


model.run_mcmc(2*number_of_iterations,repeat=3)
model.plot_chains()


# In[88]:


model.plot_distributions()


# In[90]:


model.plot_many(0,35,'Z')


# In[91]:


model.triangle_plot()


# ## With different priors

# In[5]:


t=array([0,3,5,6,8,10,22,22.2,22.5,24,25.5,26,26.5,27.5,27.75,28.5,29,29.5,31.5])
zombies=array([0,1,2,2,3,3,4,6,2,3,5,12,15,25,37,25,65,80,100])

sim=Simulation()
sim.add("S'=-β*S*Z",508.2,plot=1)
sim.add("E'=β*S*Z-ζ*E",0,plot=0)
sim.add("Z'=ζ*E-α*S*Z",.000347759,plot=1)
sim.add("R'=α*S*Z",0,plot=False)
sim.params(α=2.96e-8,β=0.000808995,ζ=60)
sim.add_data(t=t,Z=zombies,plot=1)
sim.run(0,50)

model=MCMCModel(sim,
                α=Uniform(0,.01),
                β=Uniform(0,.01),
                ζ=Normal(10,10,all_positive=True)
               )


# In[6]:


model.run_mcmc(800,repeat=2)
model.plot_chains()


# In[9]:


model.plot_many(0,35,'Z')


# In[8]:


model.plot_distributions()


# In[ ]:




