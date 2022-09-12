---
jupytext:
  formats: notebooks//ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
```

# Computational Mechanics Project #01 - Heat Transfer in Forensic Science

We can use our current skillset for a macabre application. We can predict the time of death based upon the current temperature and change in temperature of a corpse. 

Forensic scientists use Newton's law of cooling to determine the time elapsed since the loss of life, 

$\frac{dT}{dt} = -K(T-T_a)$,

where $T$ is the current temperature, $T_a$ is the ambient temperature, $t$ is the elapsed time in hours, and $K$ is an empirical constant. 

Suppose the temperature of the corpse is 85$^o$F at 11:00 am. Then, 2 hours later the temperature is 74$^{o}$F. 

Assume ambient temperature is a constant 65$^{o}$F.

## 1.

1. Use Python to calculate $K$ using a finite difference approximation, $\frac{dT}{dt} \approx \frac{T(t+\Delta t)-T(t)}{\Delta t}$.

+++

We can solve for $\frac{dT}{dt}$ in the equation $\frac{dT}{dt} \approx \frac{T(t+\Delta t)-T(t)}{\Delta t}$ with known variables of the temperature of the body at 11:00am and 2 hours later. 

Then, this value of $\frac{dT}{dt}$ can be plugged into the Newton's Law of Cooling equation to solve for K, this time using the bodytemperature after 2 hours $T$, as well as the ambient temperature $T_a$. 

I.e., 

$$\frac{T(t+\Delta t)-T(t)}{\Delta t} = -K(T - T_a)$$

Then solving for K,

$$K = -\frac{\frac{T(t+\Delta t)-T(t)}{\Delta t}}{T - T_a}$$


```{code-cell} ipython3
dT_dt = (74 - 85) / 2

K = - dT_dt / (74 - 65)

print("K =", K)
```

## 2.

2. Change your work from problem 1 to create a function that accepts the temperature at two times, ambient temperature, and the time elapsed to return $K$.

```{code-cell} ipython3
def K_constant(T_1, T_2, T_a, dt):
    '''This function returns the value of the constant K in Newton's Law of Cooling, 
    where T_1 is the first known body temperature, T_2 is the later known body temperature, 
    T_a is the ambient temperature, and dt is the time elapsed'''
    
    dT_dt = (T_2 - T_1) / dt
    K = - dT_dt / (T_2 - T_a)
    return K
    
```

```{code-cell} ipython3
# verifying the function, 

K_constant(85, 74, 65, 2)
```

## 3.

3. A first-order thermal system has the following analytical solution, 

    $T(t) =T_a+(T(0)-T_a)e^{-Kt}$

    where $T(0)$ is the temperature of the corpse at t=0 hours i.e. at the time of discovery and $T_a$ is a constant ambient temperature. 

    a. Show that an Euler integration converges to the analytical solution as the time step is decreased. Use the constant $K$ derived above and the initial temperature, T(0) = 85$^o$F. 

    b. What is the final temperature as t$\rightarrow\infty$?
    
    c. At what time was the corpse 98.6$^{o}$F? i.e. what was the time of death?

+++

Performing the Euler integration, let's take the equation defined above from before: 
$$\frac{T(t+\Delta t)-T(t)}{\Delta t} = -K[T(t+\Delta t) - T_a]$$


```{code-cell} ipython3
print('a.')
K = K_constant(85, 74, 65, 2)

# euler approx.
t = np.linspace(0, 10, 10)       # initial guess of 10 time steps, and plotting the solution for the first 10 hours
T = np.zeros(len(t))
dt = t[1] - t[0]
T[0] = 85
for i in range(1, len(t)):
    T[i] = T[i-1] - K*(T[i-1]-T_a)*dt

    
# analytical solution
t_an = np.linspace(0, 10, 50)
T_an = np.zeros(len(t_an))
t_0 = t_an[0]
T_0 = 85
T_a = 65
T_an = T_a + (T_0 - T_a)*np.exp(-K*(t_an - t_0))


plt.plot(t_an, T_an, label = 'Analytical')
plt.plot(t, T, 'o-', label = 'Euler Integration')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('10 Time Steps')
plt.legend()
plt.show()


##################################################################
K = K_constant(85, 74, 65, 2)

# euler approx.
t = np.linspace(0, 10, 40)       # 40 time steps, and plotting the solution for the first 10 hours
T = np.zeros(len(t))
dt = t[1] - t[0]
T[0] = 85
for i in range(1, len(t)):
    T[i] = T[i-1] - K*(T[i-1]-T_a)*dt

    
# analytical solution
t_an = np.linspace(0, 10, 50)
T_an = np.zeros(len(t_an))
t_0 = t_an[0]
T_0 = 85
T_a = 65
T_an = T_a + (T_0 - T_a)*np.exp(-K*(t_an - t_0))


plt.plot(t_an, T_an, label = 'Analytical')
plt.plot(t, T, 'o-', label = 'Euler Integration')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('40 Time Steps')
plt.legend()
plt.show()


##################################################################
K = K_constant(85, 74, 65, 2)

# euler approx.
t = np.linspace(0, 10, 100)       # 100 time steps, and plotting the solution for the first 10 hours
T = np.zeros(len(t))
dt = t[1] - t[0]
T[0] = 85
for i in range(1, len(t)):
    T[i] = T[i-1] - K*(T[i-1]-T_a)*dt

    
# analytical solution
t_an = np.linspace(0, 10, 50)
T_an = np.zeros(len(t_an))
t_0 = t_an[0]
T_0 = 85
T_a = 65
T_an = T_a + (T_0 - T_a)*np.exp(-K*(t_an - t_0))


plt.plot(t_an, T_an, label = 'Analytical')
plt.plot(t, T, 'o-', label = 'Euler Integration')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('100 Time Steps')
plt.legend()
plt.show()

print('It can be seen in the above three plots that the Euler integration converges to the analytical solution as the time step is decreased (i.e., the amount of time steps utilized in the Euler Approximation is increased).')
```

b. As time goes to $\infty$, the final body temperature is 65$^o$F  (i.e., the ambient temperature). This can be seen in the above plots, in which the temperature as time goes on begins to flatline at a constant 65$^o$F.

+++

c.

To determine the time of death (i.e., when the body is 98.6$^o$F), the analytical solution will be used. In the analytical solution, $T(0)$ = 98.6$^o$F, as this is the starting temperature of the body. Then, the analytical solution can be rearranged to solve for time t. In this situation, it is assumed that a body temperature of 98.6 occurs at time t = 0. 

$$t = -\frac{1}{K} * ln(\frac{T(t) - T_a}{T(0) - T_a})$$

The time after death at which the body temperature is 85$^o$F can be found. We know a body temp of T = 85 occurs at 11:00 o'clock. By knowing the time it takes from the time of death to reach a body temp of 85$^o$F, and also knowing that T = 85 occurs at 11:00 o'clock, the time of death can be found.

```{code-cell} ipython3
K = K_constant(85, 74, 65, 2)


# euler approx.
t = np.linspace(0, 10, 100)       # initial guess of 10 time steps, and plotting the solution for the first 10 hours
T = np.zeros(len(t))
dt = t[1] - t[0]
T[0] = 98.6
for i in range(1, len(t)):
    T[i] = T[i-1] - K*(T[i-1]-T_a)*dt

    
# analytical solution
t_an = np.linspace(0, 10, 50)
T_an = np.zeros(len(t_an))
t_0 = t_an[0]
T_0 = 98.6
T_a = 65
T_an = T_a + (T_0 - T_a)*np.exp(-K*(t_an - t_0))


plt.plot(t_an, T_an, label = 'Analytical')
plt.plot(t, T, 'o-', label = 'Euler Integration')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('Determing Time of Death')
plt.axhline(85, color = 'aqua', label = "85 degF line")
plt.legend()

```

```{code-cell} ipython3
K = K_constant(85, 74, 65, 2)

t = -(1/K) * np.log((85 - 65)/(98.6 - 65))

t_death = 11 - t

hours = int(t_death)
minutes = (t_death*60) % 60
seconds = (t_death*3600) % 60

print('The time of death occurs at', "%d:%02d:%02d" % (hours, minutes, seconds))
```

```{code-cell} ipython3

```
