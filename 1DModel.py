#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:38:14 2019

@author: lukesmith
"""
from random import uniform, shuffle
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as npran
import matplotlib.animation as animation

L = 400 # number of cells on road
n_iters = 400 # no. iterations
density = 0 # percentage of cars
vmax = 5 #maximum velocity
Pslow = 0.2 #probability that a car will slow down
Pnew = 0.8
tlt = 10

# Initial Data Array ----------------------------------------------------------
car_num = int(density * L)
ini = [0] * car_num + [-1] * (L - car_num) # creates an array with cars and empty spaces
shuffle(ini) # This randomises the array to the sars are spread along the road.
iterate = [ini]
block = iterate[0]
T = [0] * n_iters
tll = int(L/2.)

lane2 = [0] * L
#Creating an array which has a value of one when the traffic light is on ------
for i in range(n_iters):
    if n_iters/2 <= i <= (n_iters/2) + tlt:
        T[i] = 1


for i in range(n_iters):
    previous , current = iterate[-1] , [-1] * (L)
    
    if previous[0] == -1 and uniform(0,1) <= Pnew:
        if n_iters > 1:
            previous[0]= npran.randint(1,vmax)
    xn, xp = 0, 0
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
            if T[i] > 0 and k < tll:
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


# converts vehicles with velocity to a 1 and empty spaces to 0s----------------
test = np.array(iterate)
a = np.zeros(shape=(n_iters,L))
for i in range(L):
    for j in range(n_iters):
        a[j,i] = 1 if iterate[j][i] > -1 else 0





#Plotting --------------------------------------------------------------------


#plot for 2x2 representation       
fig, ax = plt.subplots()
im = ax.imshow(test, cmap="gist_heat", interpolation="nearest")

plt.xlabel("Space")
plt.ylabel("Time")
X=range(L)
Y= np.ones(L) * int(n_iters/2)
Y2= np.ones(L) * (int(n_iters/2) + tlt)
plt.plot(X,Y, "r-", label="Red Light")
plt.plot(X,Y2,"g-", label="Green Light")
plt.legend()

Writer = animation.writers['ffmpeg']
writer = Writer(fps=2, metadata=dict(artist='luke'), bitrate=3600)
"""
#------------------------------------------------------------------------------
##plot for animation imshow
#fig2, ax2 = plt.subplots()
#im = plt.imshow(np.expand_dims(test[0],axis=0), cmap = "gist_heat")
#ax2.grid(which='major', axis='x', linestyle='-', color='k', linewidth=0)
#ax2.axis("off")
#
#
#i=0
#def updatefig(*args):
#    global i
#    if (i<n_iters-1):
#        i += 1
#    else:
#        i=0
#    im.set_array(np.expand_dims(test[i],axis=0))
#    return im,
#ani = animation.FuncAnimation(fig2, updatefig, blit=True)
#
##writing the animation to an mp4 file

#ani.save('im.mp4', writer=writer)
#plt.show()

#------------------------------------------------------------------------------
#plot for Triangular representation

fig3, ax3 = plt.subplots()
ax3.set_ylim(1.00001,0.99999)
ax3.axes.get_yaxis().set_visible(False)
ax3.set_xlabel("Space")

x0 = range(L)
y0 = a[0,:]

line, = ax3.plot(x0, y0, '>')
Y2=np.arange(0.5,1.6,0.1)
X2 = np.full(len(Y2),(L/2))
Light, = ax3.plot(X2,Y2, "-g")

def animate(i):
    if n_iters/2 <= i <= (n_iters/2) + tlt:
        Light.set_color("red")
    else:
        Light.set_color("green")

    
    line.set_xdata(range(L))
    line.set_ydata(a[i,:])
    return line,

anim = animation.FuncAnimation(fig3, animate, frames=n_iters, blit=True)

anim.save('Road Animation.mp4', writer=writer)
plt.show()

#------------------------------------------------------------------------------
# plot for initial grid.
#fig4 , ax4 = plt.subplots()
#
#IR = plt.imshow(np.expand_dims(block,axis=0), cmap = 'inferno')
#plt.xlabel("Space")
#plt.ylabel("")
#plt.colorbar()
#fig4.set_figheight(1)
#fig4.set_figwidth(10)
#
#plt.show()



#------------------------------------------------------------------------------


#plot for animation imshow
fig2, ax2 = plt.subplots()
im2 = ax2.imshow(test[:][:L], cmap="gist_heat", interpolation="nearest")
ax2.grid(which='major', axis='x', linestyle='-', color='k', linewidth=0)
ax2.axis("off")


i=0
def updatefig(*args):
    global i
    if (i<n_iters-(L-1)):
        i += 1
    else:
        i=0
    im2.set_array(test[:][i-1:i+(L-1)])
    return im,
ani = animation.FuncAnimation(fig2, updatefig, blit=True, frames = (n_iters-L))

#writing the animation to an mp4 file

ani.save('im.mp4', writer=writer)
plt.show()
"""


