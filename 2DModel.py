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

L = 100 # number of cells on road
n_iters = 200 # no. iterations
density = 0.0 # percentage of cars
vmax = 6 #maximum velocity
vmaxR = 6 #maximum velocity

Pslow = 0.3 #probability that a car will slow down
Pnew = 0.02
Pc = 0.9
Pcr = 0.5
tlt = 0

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
                print("change")
                print(k)
                print(i)

                changeLR += 1
                
            else:
                changedL[k] = previousL[k]
                
        if previousR[k] > -1:
       
#right lane change
                
            Rdn, Rdb, Rd2n, Rd2b =1, 1, 0, 0
            while previousR[(k+Rdn)% L]<0 and k + Rdn <= L:
                Rdn += 1  
                
            while previousR[(k-Rdb)% L]<0 and k - Rdb >= 0:
                Rdb += 1  

            while previousL[(k+Rd2n)% L]<0 and k + Rd2n <= L:
                Rd2n +=1
              
            while previousL[(k-Rd2b)% L]<0 and k - Rd2b >= 0:
                Rd2b +=1   
              
                
            vR, vRn, vRb, vR2n, vR2b = 0,0,0,0,0
           
            vR = previousR[k] #V of current car
            if k + Rdn < L:            
                vRn = previousR[k + Rdn]  #V of car in from in same lane
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
    
    
    previousL,currentL,iterateL,xnL,xpL = LaneUpdate(changedL,Pnew,currentL,vmax,iterateL,xnL,xpL)
    previousR,currentR,iterateR,xnR,xpR = LaneUpdate(changedR,Pnew,currentR,vmax,iterateR,xnR,xpR)

    
    
uniqueL , countsL = np.unique(iterateL, return_counts=True)   
uniqueR , countsR = np.unique(iterateR, return_counts=True)

VL = uniqueL + 1
VR = uniqueR + 1

VelocityL = VL*countsL
VelocityR = VR*countsR

avVR = sum(VelocityR)/sum(countsR[1:])
avVL = sum(VelocityL)/sum(countsL[1:])
         
    # here is where the update rules need to go for two lanes     
    # take the row ,iterate along to decide is lanes need switching
    # need a row for right lane and left lane.
    # make two new arrays each time for this purpose
    # the output of this becomes previous ffor the next part.


# converts vehicles with velocity to a 1 for the left lane and 0.5 for the right. and empty spaces to 0s----------------
testL = np.array(iterateL)
testR = np.array(iterateR)

a = np.zeros(shape=(n_iters,L))
for i in range(L):
    for j in range(n_iters):
        a[j,i] = 0.75 if iterateL[j][i] > -1 else 0
        
aR = np.zeros(shape=(n_iters,L))
for i in range(L):
    for j in range(n_iters):
        aR[j,i] = 1 if iterateR[j][i] > -1 else 0

#testR = np.array(iterateL)
#a = np.zeros(shape=(n_iters,L))
#for i in range(L):
#    for j in range(n_iters):
#        a[j,i] = 1 if iterateR[j][i] > -1 else 0



#Plotting --------------------------------------------------------------------


#plot for 2x2 representation--------------------------------------------------       
fig, ax = plt.subplots()
im = ax.imshow(testL, cmap="gist_heat", interpolation="nearest")

plt.xlabel("Space")
plt.ylabel("Time")
X=range(L)
Y= np.ones(L) * int(n_iters/2)
Y2= np.ones(L) * (int(n_iters/2) + tlt)
plt.plot(X,Y, "r-", label="Stop")
plt.plot(X,Y2,"g-", label="Go")
plt.title("Left lane")

plt.legend()
fig.set_figheight(10)
fig.set_figwidth(10)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=2, metadata=dict(artist='luke'), bitrate=3600)


fig2, ax2 = plt.subplots()
im2 = ax2.imshow(testR, cmap="gist_heat", interpolation="nearest")
fig2.set_figheight(10)
fig2.set_figwidth(10)

plt.xlabel("Space")
plt.ylabel("Time")
X=range(L)
Y= np.ones(L) * int(n_iters/2)
Y2= np.ones(L) * (int(n_iters/2) + tlt)
plt.title("Right lane")
plt.plot(X,Y, "r-", label="Stop")
plt.plot(X,Y2,"g-", label="Go")
plt.legend()


# Animation showing the two lanes side by side --------------------------------
fig3, ax3 = plt.subplots()
ax3.set_ylim(1.5,0.25)
ax3.axes.get_yaxis().set_visible(False)
ax3.set_xlabel("Space")
x0 = range(L)
y0 = a[0,:]
y0R = aR[0,:]

line, = ax3.plot(x0, y0, '>')
lineR, = ax3.plot(x0, y0R, '>')

Y2=np.arange(0.2 ,1.6,0.1)
X2 = np.full(len(Y2),(L/2))
Light, = ax3.plot(X2,Y2, "-g")

def animate(i):
    if n_iters/2 <= i <= (n_iters/2) + tlt:
        Light.set_color("red")
    else:
        Light.set_color("green")

    
    line.set_xdata(range(L))
    line.set_ydata(a[i,:])
    lineR.set_xdata(range(L))
    lineR.set_ydata(aR[i,:])
    return line,

anim = animation.FuncAnimation(fig3, animate, frames=n_iters, blit=True)

anim.save('Road Animation.mp4', writer=writer)
plt.show()



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
