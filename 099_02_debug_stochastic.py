#!/usr/bin/env python
# coding: utf-8

# # More Debugging Stochastic Models
# 

# In[1]:


#hide
#skip
get_ipython().system(' [ -e /content ] && pip install -Uqq pyndamics3 emcee # upgrade pyndamics3 on colab')


# In[2]:


get_ipython().run_line_magic('pylab', 'inline')


# In[3]:


from pyndamics3 import Simulation,Stochastic_Simulation


# In[7]:


So=5364500
So=100
Eo=1
Io=0
β=0.2
q=2
ρ=1/5
γ=1/7
ts=130

sim=Simulation()
sim.add("S'=-β*S*I/N",So)
sim.add("E'=+β*S*I/N-ρ*E",Eo,plot=1)
sim.add("I'=+ρ*E-γ*I",Io,plot=1)
sim.add("R'=+γ*I",0)
sim.add("N=S+E+I+R")
sim.params(β=β,γ=γ,q=q,ρ=ρ,ts=ts)
sim.run(400)


# In[19]:


So=5364500
#So=100000
Eo=1
Io=0
β=0.2
q=2
ρ=1/5
γ=1/7
ts=130

sim=Stochastic_Simulation()
sim.add("-S+E",'β*S*I/N',S=So,E=Eo,I=Io)
sim.add("-E+I",'ρ*E')
sim.add("-I+R",'γ*I',R=0)
sim.add("N=S+E+I+R")
sim.params(β=β,γ=γ,q=q,ρ=ρ,ts=ts)
sim.run(500,Nsims=100)


# In[20]:


for i in range(100):
    
    plot(sim.t,sim.E[i],'bo',alpha=0.05)
    plot(sim.t,sim.I[i],'ro',alpha=0.05)


# In[ ]:




