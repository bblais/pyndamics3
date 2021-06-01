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


# In[1]:


get_ipython().run_line_magic('pylab', 'inline')


# In[2]:


from pyndamics3 import Struct


# In[3]:



# copied from http://be150.caltech.edu/2018/handouts/l12_stochastic_simulation.html

import numba

@numba.jit(nopython=True)
def _sample_discrete(probs, probs_sum):
    q = np.random.rand() * probs_sum

    i = 0
    p_sum = 0.0
    while p_sum < q:
        p_sum += probs[i]
        i += 1
    return i - 1


class Stochastic_Simulation(object):
    
    def __init__(self):
        self.components=[]
        self.equations=[]
        self.initial_values={}
        self.current_values={}
        self.ν=None
        self.state_change_strings=[]
        self.rate_equations=[]
        self.quasi=[]
        self._params={}
        self._params_keys=()
        self._params_vals=()
    
    def params(self,**kwargs):
        self._params.update(kwargs)
        self._params_keys=tuple(self._params.keys())
        self._params_vals=tuple([self._params[_] for _ in self._params_keys])
        
    def add(self,component_change_equation,rate_equation=None,plot=False,quasi=None,**kwargs):
        
        
        if "=" in component_change_equation:
            self.equations.append(component_change_equation)
            return 
        
        component_change_equation=component_change_equation.replace('+',' +')
        component_change_equation=component_change_equation.replace('-',' -')
        
        parts=component_change_equation.split()
        for part in parts:
            if not (part.startswith('-') or part.startswith('+')):
                raise SyntaxError("State change strings must start with + or -: %s" % component_change_equation)
            name=part[1:]
            if name not in self.components:
                self.components.append(name)
            
        self.state_change_strings.append(component_change_equation)            
        self.rate_equations.append(rate_equation)
        self.initial_values.update(kwargs)
        self.current_values.update(kwargs)
        self.quasi.append(quasi)

    def initialize(self):
        import numba
        import numpy as np
        num_components=len(self.components)
        num_reactions=len(self.rate_equations)
        self.ν=np.zeros((num_reactions,num_components),int)
        
        for j,(state_change,rate) in enumerate(zip(self.state_change_strings,self.rate_equations)):
            parts=state_change.split()
            for part in parts:
                if not (part.startswith('-') or part.startswith('+')):
                    raise SyntaxError("State change strings must start with + or -: %s" % component_change_equation)
                name=part[1:]
                if part[0]=='-':
                    val=-1
                else:
                    val=+1
                
                i=self.components.index(name)
                self.ν[j,i]=val
                    

        for c in self.initial_values:
            if not c in self.components:
                raise ValueError("%s not in components values." % c)
                    
        for c in self.components:
            if not c in self.initial_values:
                raise ValueError("%s not in initial values." % c)
            
            
        #func_str="@numba.jit(nopython=True)\ndef _propensity_function(population, args):\n"
        func_str="@numba.jit(nopython=True)\ndef _propensity_function_abcde(population, args):\n"

        func_str+="    "
        
        if len(self.components)>1:        
            func_str+=",".join(self.components) + " = population\n"
        else:
            func_str+=self.components[0] + ", = population\n"
                        
        if self._params_keys:
            func_str+="    "
            if len(self._params_keys)>1:        
                func_str+=",".join(self._params_keys)+ " = args\n"
            else:
                func_str+=self._params_keys[0]+ ", = args\n"
            
        func_str+="    "+"\n"

        for eq in self.equations:
            func_str+="    "+eq+"\n"


        func_str+="    "+"\n"


        func_str+="    "+"val = np.array([\n"
        for a in self.rate_equations:
            func_str+="        "+a+",\n"
        func_str+="    "+"],float)\n"

        for qi,q in enumerate(self.quasi):
            if not q:
                continue
                
            func_str+="    "+f"if ({q}):\n"
            func_str+="    "+"    "+f"val[{qi}]=0\n"
        
            func_str+="    "+f"if ((A==0) or (B==0)):\n"
            func_str+="    "+"    "+f"raise ValueError()\n"
                
        
        func_str+="    "+"return val"
        
        
        self.func_str=func_str
            
        exec (func_str,globals())                      
        self.propensity_function=_propensity_function_abcde
        
    def run(self,t_max,Nsims=1,num_iterations=1001,):
        from tqdm import tqdm
        
        
        if self.ν is None:
            self.initialize()

        _propensity_function=self.propensity_function
        
        
        @numba.jit(nopython=True)
        def _ssa(update, population_0, time_points, args):
            # Initialize output
            pop_out = np.empty((len(time_points), update.shape[1]), dtype=np.int64)

            # Initialize and perform simulation
            i_time = 1
            i = 0
            t = time_points[0]
            population = population_0.copy()
            pop_out[0,:] = population
            extinction_time=-1.0
            previous_t=t
            while i < len(time_points):
                while t < time_points[i_time]:
                    # draw the event and time step
                    event, dt = _draw(population, args)

                    # Update the population
                    population_previous = population.copy()
                    population += update[event,:]

                    # Increment time
                    previous_t=t
                    t += dt


                if dt==1e500 and extinction_time<0.0:
                    extinction_time=previous_t

                # Update the index (Have to be careful about types for Numba)
                i = np.searchsorted((time_points > t).astype(np.int64), 1)

                # Update the population
                for j in np.arange(i_time, min(i, len(time_points))):
                    pop_out[j,:] = population_previous

                # Increment index
                i_time = i

            return pop_out,extinction_time
        
        @numba.jit(nopython=True)
        def _draw(population, args):
            """
            Draws a reaction and the time it took to do that reaction.

            Assumes that there is a globally scoped function
            `prop_func` that is Numba'd with nopython=True.
            """
            # Compute propensities
            props = _propensity_function(population, args)

            # Sum of propensities
            props_sum = np.sum(props)

            if props_sum==0:
                time=1e500
                rxn=0
            else:

                # Compute time
                time = np.random.exponential(1 / props_sum)

                # Draw reaction given propensities
                rxn = _sample_discrete(props, props_sum)

            return rxn, time

        
        self.all_storage=[]
        
        disable=Nsims==1
        
        population_0=np.array([self.initial_values[c] for c in self.components], dtype=int)
        time_points=np.linspace(0,t_max,num_iterations)        
        args = np.array(self._params_vals)
        n_simulations = Nsims

        # Initialize output array
        pops = np.empty((n_simulations, len(time_points), len(population_0)), dtype=int)
        extinction_time=np.empty(n_simulations,dtype=np.float64)

        # Run the calculations
        for _i in tqdm(range(n_simulations),disable=disable):
            pops[_i,:,:],extinction_time[_i] = _ssa(self.ν, 
                                        population_0, time_points, args=args)            

        self.t=time_points
        self.extinction_times=extinction_time
        D={}
        for _i,c in enumerate(self.components):
            v=pops[:,:,_i]
            if v.shape[0]==1:
                v=v.ravel()
                
            setattr(self, c,v)
            D[c]=v
        
        for eq in self.equations:
            exec(eq,D)
            name=eq.split('=')[0].strip()
            setattr(self, name,D[name])

  


# In[4]:


β=0.2
γ=0.1
So=990
Io=10

stoch_sim=sim=Stochastic_Simulation()
sim.add("+A",'A',A=1)
sim.add("-A",'A**2/10')
sim.run(10,Nsims=100)

for i in range(100):    
    plot(sim.t,sim.A[i],'b-o',alpha=0.05)


# In[11]:


β=0.2
γ=0.1
So=990
Io=10

stoch_sim=sim=Stochastic_Simulation()
sim.add("+A",'A',A=1)
sim.add("-A",'A**2/10',quasi='A==1')
sim.run(10,Nsims=100)

for i in range(100):    
    plot(sim.t,sim.A[i],'b-o',alpha=0.05)


# In[5]:


print(sim.func_str)


# In[7]:


@numba.jit(nopython=True)
def _propensity_function_abcde(population, args):
    A, = population
    
    
    val = np.array([
        A,
        A**2/10,
    ],float)
    return val


# In[ ]:




