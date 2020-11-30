# Pyndamics3
> An Update of the Python Numerical Dynamics library pyndamics


This file will become your README and also the index of your documentation.

## Install

`pip install pyndamics3`

## How to use

Some simple examples.

```python
from pyndamics3 import Simulation
```

    pyndamics3  version  0.0.2


```python
sim=Simulation()
sim.add("p'=a*p*(1-p/K)",1,plot=True)
sim.params(a=1,K=50)
sim.run(50)
```


![png](docs/images/output_6_0.png)



    <Figure size 432x288 with 0 Axes>

