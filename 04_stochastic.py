#!/usr/bin/env python
# coding: utf-8

# # Exploring Stochastic Models

# In[1]:


#hide
#skip
get_ipython().system(' [ -e /content ] && pip install -Uqq pyndamics3 emcee # upgrade pyndamics3 on colab')


# In[2]:


get_ipython().run_line_magic('pylab', 'inline')


# In[3]:


from pyndamics3 import Simulation,Stochastic_Simulation


# ## Stochastic SIR Model

# In[4]:


β=0.2
γ=0.1
So=990
Io=10

dynamic_sim=sim=Simulation()
sim.add("N=S+I+R")
sim.add("S'=-β*S*I/N",So)
sim.add("I'=+β*S*I/N-γ*I",Io)
sim.add("R'=+γ*I",0)
sim.params(β=β,γ=γ)
sim.run(200)


stoch_sim=sim=Stochastic_Simulation()
sim.add("-S+I",'β*S*I/N',S=So,I=Io)
sim.add("-I +R",'γ*I',R=0)
sim.add("N=S+I+R")
sim.params(β=β,γ=γ)
sim.run(200,Nsims=100)

for i in range(100):
    
    plot(sim.t,sim.S[i],'bo',alpha=0.05)
    plot(sim.t,sim.I[i],'ro',alpha=0.05)

plot(dynamic_sim.t,dynamic_sim.S,'c-')
plot(dynamic_sim.t,dynamic_sim.I,'m-')

print(sim.func_str)


# In[4]:


flut = array([0,1,2,3,4,5,6,7,8,9,10,11,12,13])
flui = array([3,8,26,76,225,298,258,233,189,128,68,29,14,4])


# In[5]:


from pyndamics3.mcmc import *


# In[10]:


β=1.9
γ=0.5
So=763
Io=1

dynamic_sim=sim=Simulation()
sim.add("N=S+I+R")
sim.add("S'=-β*S*I/N",So)
sim.add("I'=+β*S*I/N-γ*I",Io)
sim.add("R'=+γ*I",0)
sim.params(β=β,γ=γ)
sim.add_data(t=flut,I=flui)
sim.run(20)


# In[11]:


model=MCMCModel(sim,β=Uniform(0,5),
               γ=Uniform(0,5))


# In[12]:


number_of_iterations=100
model.run_mcmc(number_of_iterations,repeat=3)
model.plot_chains()


# In[13]:


plot(sim.t,sim.I)
plot(flut,flui,'ko',ms=10,lw=3,)


# In[19]:


sim.β,sim.γ


# In[17]:


stoch_sim=Stochastic_Simulation()
stoch_sim.add("-S+I",'β*S*I/N',S=So,I=Io)
stoch_sim.add("-I +R",'γ*I',R=0)
stoch_sim.add("N=S+I+R")
stoch_sim.params(β=1.9732213241997467,γ=1.9732213241997467)
stoch_sim.add_data(t=flut,I=flui)
stoch_sim.run(20,Nsims=100)


# In[15]:


for i in range(100):    
    plot(stoch_sim.t,stoch_sim.I[i],'ro',alpha=0.05)
    
plot(flut,flui,'ko',ms=10,lw=3,)    


# In[16]:


stoch_model=MCMCModel(stoch_sim,β=Uniform(0,5),
               γ=Uniform(0,5))


# In[ ]:




