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


# In[10]:


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


# this raises an error

# In[11]:


β=1.9732213241997467
γ=0.47521873806558335

β=.5
γ=1

So=763
Io=1

stoch_sim_err=Stochastic_Simulation()
stoch_sim_err.add("-S+I",'β*S*I/N',S=So,I=Io)
stoch_sim_err.add("-I +R",'γ*I',R=0)
stoch_sim_err.add("N=S+I+R")
stoch_sim_err.params(β=β)
stoch_sim_err.add_data(t=flut,I=flui)
stoch_sim_err.run(20,Nsims=100)


# In[12]:


stoch_sim['I']


# In[13]:


I=stoch_sim.components[1]


# In[ ]:





# In[ ]:





# In[14]:


for i in range(100):    
    plot(stoch_sim.t,stoch_sim.I[i],'ro',alpha=0.05)
    
plot(flut,flui,'ko',ms=10,lw=3,)    


# In[15]:


print(stoch_sim.func_str)


# In[ ]:





# In[ ]:





# In[16]:


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


# In[17]:


stoch_sim.I[i]


# In[18]:


stoch_model=MCMCModel(stoch_sim,β=Uniform(0,5),
               γ=Uniform(0,5))


# In[19]:


number_of_iterations=500
stoch_model.run_mcmc(number_of_iterations,repeat=3)

stoch_model.plot_chains()


# In[15]:


stoch_model.plot_distributions()


# In[20]:


stoch_sim._params['β']


# In[21]:


stoch_sim.run(20,Nsims=100)


# In[22]:


for i in range(100):    
    plot(stoch_sim.t,stoch_sim.I[i],'ro',alpha=0.05)
    
plot(flut,flui,'ko',ms=10,lw=3,)    


# ## debug with flu data

# In[ ]:


get_ipython().run_line_magic('pylab', 'inline')


# In[5]:


from pyndamics3 import Simulation,Stochastic_Simulation


# In[6]:


flut = array([0,1,2,3,4,5,6,7,8,9,10,11,12,13])
flui = array([3,8,26,76,225,298,258,233,189,128,68,29,14,4])


# In[7]:


flui=array([0,72,112,145,194])
flut=array([1,2,3,4,5])


# In[ ]:


β=1.9732213241997467
γ=0.47521873806558335

β=.5
γ=1

So=763
Io=1

stoch_sim=Stochastic_Simulation()
stoch_sim.add("-S+I",'β*S*E/N',S=So,I=Io)
stoch_sim.add("-E+I",'ζ*S*I',E=0)
stoch_sim.add("-I +R",'γ*I',R=0)
stoch_sim.add("N=S+I+R")
stoch_sim.params(β=β,γ=γ,ζ=.1)
stoch_sim.add_data(t=flut,I=flui)
stoch_sim.run(20,Nsims=100)


# In[ ]:





# ## Debug with vampire data

# In[28]:


get_ipython().run_line_magic('pylab', 'inline')


# In[29]:


from pyndamics3 import Simulation,Stochastic_Simulation


# In[30]:


tbt=array([0,72,112,145,194])
tbv=array([1,2,3,4,5])


# In[31]:


So=100
Vo=1
Eo=0
t_max=1.1*tbt.max()
β=0.5
γ=0.5
ζ=0.5
δ=0.5

stoch_sim=sim=Stochastic_Simulation()
sim.add("-S+E",'β*S*V/N',S=So,V=Vo)
sim.add("-E+V",'γ*S*V',E=Eo)
sim.add("-E+X",'ζ*S*V',X=1)
sim.add("-V+R",'δ*S*V',R=0)
sim.add("N=S+E+V+X")
sim.params(β=0.03,γ=0.00047,ζ=ζ,δ=δ)
sim.add_data(t=tbt,V=tbv)
sim.run(t_max,Nsims=100)


# In[ ]:




