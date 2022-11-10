#!/usr/bin/env python
# coding: utf-8

# In[1]:


#| default_exp fit


# From https://people.duke.edu/~ccc14/sta-663/CalibratingODEs.html

# In[2]:


#| export
import lmfit


# # Parameter Fitting (using lmfit)

# In[3]:


get_ipython().run_line_magic('pylab', 'inline')


# In[4]:


from pyndamics3 import Simulation


# In[5]:


#| export
import numpy as np
from lmfit import minimize, Parameters, report_fit


# In[6]:


from scipy.integrate import odeint


# In[ ]:


def f(xs, t, ps):
    """Receptor synthesis-internalization model."""
    try:
        a = ps['a'].value
        b = ps['b'].value
    except:
        a, b = ps
    x = xs

    return a - b*x

def g(t, x0, ps):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(f, x0, t, args=(ps,))
    return x

def residual(ps, ts, data):
    x0 = ps['x0'].value
    model = g(ts, x0, ps)
    return (model - data).ravel()


# In[7]:


a = 2.0
b = 0.5
true_params = [a, b]
x0 = 10.0


t = np.linspace(0, 10, 10)
data = g(t, x0, true_params)
data += np.random.normal(size=data.shape)

# set parameters incluing bounds
params = Parameters()
params.add('x0', value=float(data[0]), min=0, max=100)
params.add('a', value= 1.0, min=0, max=10)
params.add('b', value= 1.0, min=0, max=10)

# fit model and find predicted values
result = minimize(residual, params, args=(t, data), method='leastsq')
final = data + result.residual.reshape(data.shape)

# plot data and fitted curves
plot(t, data, 'o')
plot(t, final, '--', linewidth=2, c='blue');

# display fitted statistics
report_fit(result)


# In[8]:


sim=Simulation()
sim.add("x'=a-b*x",10,plot=True)
sim.add_data(t=t,x=data,plot=True)
sim.params(a=5.,b=5)
sim.run(10)


# In[8]:


#| export
def Parameter(name,**kwargs):
    from lmfit import Parameters
    params = Parameters()
    params.add(name, **kwargs)
    
    return params


# In[9]:


#| export
def residual(ps, sim):
    
    params={}
    for key in ps.keys():
        if key.startswith('initial_'):
            name=key.split('initial_')[1]
            _c=sim.get_component(name)
            _c.initial_value=ps[key].value
            
        else:
            params[key]=ps[key].value
            
    sim.params(**params)
    
    # run the sim
    sim.run_fast()
    
    # compare with data
    values=[]
    for _c in sim.components+sim.assignments:
        if not _c.data:
            continue
        t=np.array(_c.data['t']).ravel()
        y=np.array(_c.data['value']).ravel()
        y_fit=sim.interpolate(t,_c.name)

        if any(np.isnan(y_fit)):
            values.append([-np.inf])
        elif any(abs(y_fit)>1e100):
            values.append([-np.inf])
        else:
            values.append(y-y_fit)

    return np.concatenate(values).ravel()
    

def fit(sim,
       *args,method='leastsq'):
    
    found_data=False
    for _c in sim.components+sim.assignments:
        if not _c.data:
            continue
        found_data=True
        
    if not found_data:
        raise ValueError("No data given for the simulation...can't fit.")
    
    
    
    from lmfit import Parameters,minimize
    
    fitparams=Parameters()
    for arg in args:
        fitparams+=arg
    
    for key in fitparams.keys():
        if key.startswith('initial_'):
            name=key.split('initial_')[1]
            try:
                _c=sim.get_component(name)
            except IndexError:
                raise ValueError("%s is a bad initial variable because %s is not a variable in the dynamical model." % (key,name))
        else:
            if not key in sim.original_params:
                raise ValueError("%s is not a parameter in the dynamical model.  Parameters are %s" % (key,str(sim.original_params)))

    
    result = minimize(residual, fitparams, args=(sim,), method=method)    
    
    params={}
    ps=result.params
    for key in ps.keys():
        if key.startswith('initial_'):
            name=key.split('initial_')[1]
            _c=sim.get_component(name)
            _c.initial_value=ps[key].value
            
        else:
            params[key]=ps[key].value
            
    sim.params(**params)
    
    
    return result


# In[10]:


results=fit(sim,
   Parameter('a',value=1,min=0),
   Parameter('b',value=5,min=0),
   Parameter('initial_x',value=10,min=0),
           method='nelder')

report_fit(result)


# In[11]:


sim.run(10)


# In[12]:


residual(result.params,sim)


# In[13]:


result


# In[14]:


a=result.params['a']


# In[15]:


a.stderr


# ## Example with Growth

# In[16]:


t=array([7,14,21,28,35,42,49,56,63,70,77,84],float)
h=array([17.93,36.36,67.76,98.10,131,169.5,205.5,228.3,247.1,250.5,253.8,254.5])


# ### Linear Model

# In[17]:


sim=Simulation()
sim.add("h'=a",1,plot=True)
sim.add_data(t=t,h=h,plot=True)
sim.params(a=1)
sim.run(0,90)


# In[18]:


results=fit(sim,
   Parameter('a',value=1,min=0),
   Parameter('initial_h',value=10,min=0))

results


# In[19]:


sim.run(0,90)


# ### Logistic

# In[20]:


sim=Simulation()
sim.add("h'=a*h*(1-h/K)",1,plot=True)
sim.add_data(t=t,h=h,plot=True)
sim.params(a=2,K=100)
sim.run(0,90)


# In[21]:


results=fit(sim,
   Parameter('a',value=1,min=0.01,max=2),
   Parameter('K',value=1,min=0,max=400),
   Parameter('initial_h',value=10,min=0,max=50))

report_fit(results)

sim.run(0,90)


# the fit is lousy -- bad initial guesses, method possibly a problem.   Retrying with powell method.

# In[22]:


results=fit(sim,
   Parameter('a',value=1,min=0.01,max=2),
   Parameter('K',value=1,min=0,max=400),
   Parameter('initial_h',value=10,min=0,max=50),method='powell')

report_fit(results)

sim.run(0,90)


# much better!

# In[ ]:




