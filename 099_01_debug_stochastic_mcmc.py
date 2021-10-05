#!/usr/bin/env python
# coding: utf-8

# # Debugging Stochastic Models

# In[1]:


#hide
#skip
get_ipython().system(' [ -e /content ] && pip install -Uqq pyndamics3 emcee # upgrade pyndamics3 on colab')


# In[2]:


get_ipython().run_line_magic('pylab', 'inline')


# In[3]:


from pyndamics3 import Simulation,Stochastic_Simulation


# In[4]:


from pyndamics3.mcmc import *


# In[5]:


flut = array([0,1,2,3,4,5,6,7,8,9,10,11,12,13])
flui = array([3,8,26,76,225,298,258,233,189,128,68,29,14,4])


# In[26]:


β=1.9732213241997467
γ=0.47521873806558335

β=.5
γ=1

So=763
Io=1

stoch_sim=Stochastic_Simulation()
stoch_sim.add("-S+I",'β*S*I/N',S=So,I=Io)
stoch_sim.add("-I +R",'γ*I',R=0)
stoch_sim.add("N=S+I+R")
stoch_sim.params(β=β,γ=γ)
stoch_sim.add_data(t=flut,I=flui)
stoch_sim.run(20,Nsims=100)


# In[27]:


stoch_sim['I']


# In[28]:


I=stoch_sim.components[1]


# In[ ]:





# In[ ]:





# In[29]:


for i in range(100):    
    plot(stoch_sim.t,stoch_sim.I[i],'ro',alpha=0.05)
    
plot(flut,flui,'ko',ms=10,lw=3,)    


# In[30]:


print(stoch_sim.func_str)


# In[ ]:





# In[ ]:





# In[31]:


dynamic_sim=sim=Simulation()
sim.add("N=S+I+R")
sim.add("S'=-β*S*I/N",So)
sim.add("I'=+β*S*I/N-γ*I",Io)
sim.add("R'=+γ*I",0)
sim.params(β=β,γ=γ)
sim.add_data(t=flut,I=flui)
sim.run(20)

plot(sim.t,sim.I)
plot(flut,flui,'ko',ms=10,lw=3,)


# In[32]:


stoch_sim.I[i]


# In[33]:


stoch_model=MCMCModel(stoch_sim,β=Uniform(0,5),
               γ=Uniform(0,5))


# In[34]:


number_of_iterations=500
stoch_model.run_mcmc(number_of_iterations,repeat=3)

stoch_model.plot_chains()


# In[17]:


stoch_model.plot_distributions()


# In[23]:


stoch_sim._params['β']


# In[24]:


stoch_sim.run(20,Nsims=100)


# In[25]:


for i in range(100):    
    plot(stoch_sim.t,stoch_sim.I[i],'ro',alpha=0.05)
    
plot(flut,flui,'ko',ms=10,lw=3,)    


# In[ ]:




