#!/usr/bin/env python
# coding: utf-8

# In[1]:


# default_exp chem


# In[2]:


get_ipython().run_line_magic('pylab', 'inline')


# In[ ]:





# In[3]:


#export
def ChemSimulation(eqnstr,**kwargs):
    from pyndamics3 import Simulation
    lines=eqnstr.strip().split('\n')
    
    components=set([])
    parameters=set([])
    for line in lines:
        lhs=[_.strip() for _ in line.split('--')[0].strip().split('+')]
        rhs=[_.strip() for _ in line.split('->')[-1].strip().split('+')]
        middle=line.split('--')[1].lstrip('\t ->').split('->')[0].rstrip('\t ->')

        components=components.union(set(lhs+rhs))
        parameters=parameters.union([middle])
        
        print(lhs,middle,rhs)
        
        
    components=sorted(components)
    parameters=sorted(parameters)
    diffeqs=[]
    for c in components:
        eqn="%s' = " % c

        for line in lines:
            lhs=[_.strip() for _ in line.split('--')[0].strip().split('+')]
            rhs=[_.strip() for _ in line.split('->')[-1].strip().split('+')]
            middle=line.split('--')[1].lstrip('\t ->').split('->')[0].rstrip('\t ->')


            sign='0'
            if c in lhs:
                sign='-'
                if len(rhs)>1:
                    sign+='%d*' % len(rhs)


            if c in rhs:
                sign='+'
                if len(lhs)>1:
                    sign+='%d*' % len(lhs)

            if sign=='0':
                continue

            plhs='*'.join(lhs)

            eqn+=f" {sign}{middle}*{plhs}"


        diffeqs.append(eqn)   
        
    sim=Simulation()
        
    for c,d in zip(components,diffeqs):
        if not c in kwargs:
            c0=0
        else:
            c0=kwargs[c]
        sim.add(d,c0)
        
    for p in parameters:
        if p in kwargs:
            sim.params(**{p:kwargs[p]})
        
    return sim


# In[4]:


sim=ChemSimulation("""
A   -->k1->   B
B   -->k_1->  A
A+B -->k2->   C
C   -->k_2->  A+B
""",A=1,B=1,C=0,k1=1,k_1=2,k2=3,k_2=4)

print()

print(sim.equations())


# In[5]:


for c in sim.components:
    c.plot=1


# In[6]:


sim.run(10)


# In[ ]:




