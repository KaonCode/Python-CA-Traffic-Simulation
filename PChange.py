#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 00:23:55 2019

@author: lukesmith
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 19:06:39 2019

@author: lukesmith
"""

from random import uniform, shuffle
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as npran
import matplotlib.animation as animation

L = 100 # number of cells on road
n_iters = 200 # no. iterations
density = 0 # percentage of cars
vmax = 4 #maximum velocity
vmaxR = 6 #maximum velocity
Pnew = 0.5
Pslow = 0.3 #probability that a car will slow down
Pc = 0.9
Pcr = 0.4
tlt = 20
changedRL=np.zeros((100,100))
changedLR=np.zeros((100,100))
jk = np.linspace(0.1,Pc,100)
for IT in range(100):
    print(IT)
    for z in range(100):
    
    
        # Initial Data Array ----------------------------------------------------------
        car_num = int(density * L)
        ini = [0] * car_num + [-1] * (L - car_num) # creates an array with cars and empty spaces
        shuffle(ini) # This randomises the array to the sars are spread along the road.
        iterateL = [ini]
        block = iterateL[0]
        T = [0] * n_iters
        tll = int(L/2.)
        iterateR = [[-1] * L]
        #Creating an array which has a value of one when the traffic light is on ------
        for i in range(n_iters):
            if n_iters/2 <= i <= (n_iters/2) + tlt:
                T[i] = 1
        
        #------------------------------------------------------------------------------
        changeLR = 0
        changeRL = 0
        
        def LaneUpdate(previous,Pnew,current,vmax,iterate,xn,xp):            
            if previous[0] == -1 and uniform(0,1) < Pnew:
                if n_iters > 1:
                    previous[0]= npran.randint(1,vmax)      
                    
            for k in range(L-1, -1 , -1):  
        #        print(k)           
                if previous[k] > -1:
                    vi = previous[k]
                    d, b = 1,0
                    
                    while (k + b) < tll :
                        b+= 1
        
                    while previous[(k + d) % L] < 0:
                        d += 1
                    
                    #find the values of velocity of the next car so that the follower
                    #can take into account the next position on the of the leader.
                    vn=0
                    if k+d < L:
                        vn = previous[(k+d)]
                        
                    if current[xn] <= previous[xp]:
                        vn = 0
                    if (k+d) >= L:
                        d = 100
                    #increase speed up to max.
                    #cars do not move further than next car (Velocity of next car considered)
                    if T[i] > 0 and k  < tll:
                        vtemp = min(vi+1, (d + vn ) - 1, vmax , (b-1))
                    else:
                        vtemp = min(vi+1, (d + vn ) - 1, vmax)
                        
                                                    
                    v = max(vtemp - 1, 0) if uniform(0,1) < Pslow else vtemp
                    #prob p that a car slows otherwise velocity maintained
                        
                    if (k+v) < L:
                        current[(k + v)] = v # allows the cars to exit the screen
                        xn = k + v
                        xp = k
                    else:
                        current[(k)] = -1
        
        
                
            
            iterate.append(current)
            return previous,current,iterate,xn,xp,
        
        
        
        
        #------------------------------------------------------------------------------       
        for i in range(n_iters):
            previousL , changedL = iterateL[-1] , [-1] * (L)
            previousR , changedR = iterateR[-1] , [-1] * (L)
            currentL = [-1] * (L)
            currentR = [-1] * (L)
            for k in range(L-1, -1, -1):
                if previousL[k] > -1:
        #left lane change
                #obtaining the distances to the car in front, the car in fron in the next lane
                #and the car behind inthe next lane.
                    Ldn, Ld2n, Ld2b =1 , 0, 0
                    while previousL[(k+Ldn)% L]<0 and k + Ldn <= L:
                        Ldn += 1
        
                    while previousR[(k+Ld2n)% L]<0 and k + Ld2n <= L:
                        Ld2n +=1  
        
                        
                    while previousR[(k-Ld2b)% L]<0 and k - Ld2b >= 0:
                        Ld2b +=1  
                  
                        
                    vL,vLn,vL2n,vL2b = 0,0,0,0
                    vL = previousL[k] #V of current car
                    
                    if k + Ldn < L:
                        vLn = previousL[k + Ldn]  #V of car in fron in same lane
                        
                    if k + Ld2n < L:            
                        vL2n = previousR[k+Ld2n] 
                        
                    if k - Ld2b > 0:                
                        vL2b = previousR[k-Ld2b]
                    
                    if vL > Ldn + vLn  and vL < Ld2n + vL2n and Ld2b + vL > vL2b and uniform(0,1) < Pc:
                        changedR[k] = previousL[k]
        
                        changeLR += 1
                        
                    else:
                        changedL[k] = previousL[k]
                        
                if previousR[k] > -1:
               
        #right lane change
                        
                    Rdn,Rdb, Rd2n, Rd2b =1, 1, 0, 0
                    while previousR[(k+Rdn)% L]<0 and k + Rdn <= L:
                        Rdn += 1  
                        
                    while previousR[(k-Rdb)% L]<0 and k - Rdb >= 0:
                        Rdb += 1  
        
                    while previousL[(k+Rd2n)% L]<0 and k + Rd2n <= L:
                        Rd2n +=1
                      
                    while previousL[(k-Rd2b)% L]<0 and k - Rd2b >= 0:
                        Rd2b +=1   
                      
                        
                    vR,vRn,vRb,vR2n,vR2b = 0,0,0,0,0
                   
                    vR = previousR[k] #V of current car
                    if k + Rdn < L:            
                        vRn = previousR[k + Rdn]  #V of car in fron in same lane
                    if k -Rdb > 0:
                        vRb = previousR[k - Rdb]  # V of the car behind in the same lane
                    if k + Rd2b > 0:
                        vR2b =previousL[k - Rd2b]  #v of the car behind in the next lane
                    if k + Rd2n < L:
                        vR2n =previousL[k + Rd2n]  #v of the car in front in the net lane
                    
                    if vR < 3  and vR < vRb + Rdb and vR2b < Rd2b + vR and vR < Rd2n + vR2n and uniform(0,1) < Pcr:
                        changedL[k] = previousR[k]
                        changeRL += 1
                        
        
                    else:
                        changedR[k] = previousR[k]
        
        
            
            xnL, xpL, xnR , xpR = 0, 0, 0, 0
            
    
            zx=int(z*10)
            previousL,currentL,iterateL,xnL,xpL = LaneUpdate(changedL,Pnew,currentL,vmax,iterateL,xnL,xpL)
            previousR,currentR,iterateR,xnR,xpR = LaneUpdate(changedR,Pnew,currentR,vmax,iterateR,xnR,xpR)
            changedRL[z][IT]=changeRL
            changedLR[z][IT]=changeLR


avchangedLR= np.zeros(100)
avchangedRL=np.zeros(100)
    
for i in range(100):
    avchangedRL[i]=np.mean(changedRL[i])
    avchangedLR[i]=np.mean(changedLR[i])


plt.plot(avchangedRL,avchangedLR,"bx")
X= range(100)
Y=X
plt.plot(X,Y,'k-')   
plt.xlabel("no. Changes R to L")    
plt.ylabel("no. Changes L to R")    
    