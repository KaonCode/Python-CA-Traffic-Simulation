#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:38:14 2019

@author: lukesmith
"""
from random import uniform, shuffle
import matplotlib.pyplot as plt
import numpy as np
"""
import matplotlib.cm as cm
import matplotlib.animation as animation
"""

L = 100 # number of cells on road
n_iters = 100 # no. iterations
density = 0.2 # percentage of cars

vmax = 5 #maximum velocity
p = 0.3 #probability that a car will slow down


car_num = int(density * L)
ini = [0] * car_num + [-1] * (L - car_num) 
    # creates an array with cars and empty spaces
shuffle(ini) 
    # This randomises the array to the sars are spread along the road.

iterate = [ini]

for i in range(n_iters):
    prev,cur = iterate[-1],[-1] * L

    for k in range(L):
        if prev[k] > -1:
            vi = prev[k]
            d = 1
            while prev[(k + d) % L] < 0:
                d += 1

            vtemp = min(vi+1, d - 1, vmax)  #increase speed up to max. 
                                            #cars do not move further than next car
            v = max(vtemp - 1, 0) if uniform(0,1) < p else vtemp #probability p a car hits the brakes otherwise velocity is sustained
#            if (k+v) < L:
#                    cur[(k + v)] = v
            if (k+v) < L:
                cur[(k + v)] = v # allows the cars to exit the screen
            else:
                cur[(k)] = -1
      

    iterate.append(cur)


a = np.zeros(shape=(n_iters,L))
for i in range(L):
    for j in range(n_iters):
        a[j,i] = 1 if iterate[j][i] > -1 else 0
 
#plotting
test = np.array(iterate)

fig, ak = plt.subplots()
im = ak.imshow(test, cmap="gist_heat", interpolation="nearest")
plt.xlabel("Space")
plt.ylabel("Time")
plt.show()


