#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 21:56:16 2017

@author: jitendrian
"""
import csv
import numpy as np
import os
from multiprocessing import Pool
import subprocess
import shutil
import pandas as pd

#flock behaviour
def exchange(arr):
        J = arr[0]
        A = arr[1]
        chi = arr[2]
        eta = arr[3]
        Am = arr[4]
        dt = arr[5]
        T = arr[6]
        L = arr[7]
        IR = arr[8]
        mx = arr[9]
        my = arr[10]
        mu = arr[11]
        f = arr[12]
        Alpha = arr[13]
        Beta = arr[14]
        name = arr[15]

	G = int(T/dt)
	N = int(L**2)
	# J,A = exchange interaction
	# chi = spin inertia
	# eta =
	# dt = grid size of time
	# T = total time
	# N = flock size
	# L = length of a square box
	# IR = interaction radius
	# G =  total grids in T time
	# f = random gaussian noise
	# mx = kx wave number
	# my = ky wave number
	# f = strength of random noise
	# Alpha = strength of cohesion
    	"*************************************************************************"
	## make arrays of size N.
	vxi = np.zeros(N)
	vyi = np.zeros(N)
	vxf = np.zeros(N)
	vyf = np.zeros(N)
	rxi = np.zeros(N)
	ryi = np.zeros(N)
	rxf = np.zeros(N)
	ryf = np.zeros(N)
	thei = np.zeros(N)
	thef = np.zeros(N)
	si = np.zeros(N)
	sf = np.zeros(N)

	# vxi = x-component of initial velocity
	# vyi = y-component of initial velocity
	# vxf = x-component of next velocity
	# vyf = y-component of next velocity
	# rxi = x-component of initial position
	# ryi = y-component of initial position
	# rxf = x-component of next position
	# ryf = y-component of next position
	# thi = initial theta
	# thf = next theta
	# si = initial spin
	# sf = spin in next instant
	#images

	"*************************************************************************"
	#distance between ath bird and bth bird on a torus
	def dist_ab(rax,ray,rbx,rby):

		if  rax <= IR:
			if L-IR<= rbx <= L:
				Dx = (rbx - rax) - L
			else:
				Dx = rbx - rax
		elif L-IR <= rax <= L:
			if 0 <= rbx <= IR:
				Dx = rbx - rax + L
			else:
				Dx = rbx -rax
		else:
			Dx = rbx - rax

		if ray <= IR:
			if L-IR <= rby <= L:
				Dy = rby - ray - L
			else:
				Dy = rby - ray
		elif L-IR <= ray <= L:
			if 0 <= rby <= IR:
				Dy = rby - ray + L
			else:
				Dy = rby -ray
		else:
			Dy = rby - ray

		D = np.sqrt(Dx**2 + Dy**2)

		return [D, Dx, Dy]
	"*************************************************************************"
        "*************************************************************************"
        #r = distance between particles
        #rH = Radius hardball
        #VH = Potential Hardball
        #rC = Radius comfort zone
        #FC = slope of comfort zone
        #C = potential intercept
        #rA = radius attraction zone
        #VA = potential attraction zone
        #VIR = potential inside IR

        def coffecient(r):
	    rH = 0.2
	    VH = -1000
	    FC = 0.5
	    C = 0.25
	    rA = 0.8
	    IR = 1.2
	    VIR = 1
	    if 0<r<rH:
		    return VH
	    elif rH<r<rA:
		    return FC*r + C
	    elif rA<r<IR:
		    return VIR
	    elif IR<r:
		    return 0
        "**************************************************************************"
	#cohesion>summation over all b's in neighbour of a, f_ab*(r_baXV_b)
	def cohesion(a):
		sC_ab = 0
                vax = np.cos(thei[a])
		vay = np.sin(thei[a])
                for i in range(N):
			if i != a:
				d = dist_ab(rxi[a],ryi[a],rxi[i],ryi[i])
				rbax = (d[1]/d[0]) 
                                rbay = (d[2]/d[0])	
				sC_ab = sC_ab + coffecient(d[0])*((vax*rbay)-(vay*rbax))
		return sC_ab
	"*************************************************************************"
	#exchange interaction>summation over all b's in neighbour of a, (J_ab*(v_aXv_b))
	def interaction(a):
		sJ_ab = 0
		for i in range(N):
			if i!=a:
				d = dist_ab(rxi[a],ryi[a],rxi[i],ryi[i])
				if d[0] <= IR:
					phi = ((d[1]*np.cos(thei[a])) + (d[2]*np.sin(thei[a])))/d[0]
					sJ_ab = sJ_ab + (J - phi*A)*np.sin(thei[i]-thei[a])
		return sJ_ab

	"*************************************************************************"
	"             initiation of velocity, spin and position                   "
	"*************************************************************************"
	q = 0
	for t in range(int(L)):
		for w in range(int(L)):
			rxi[q], ryi[q] = t, w
			q = q+1

	for t in range(N):
		thei[t] = np.pi/2 + (Am*np.sin(mu + (2*np.pi*mx*rxi[t]/L) + (2*np.pi*my*ryi[t]/L)))
		vxi[t] = np.cos(thei[t])
		vyi[t] = np.sin(thei[t])
		si[t] = Am*np.sin(mu + (2*np.pi*mx*rxi[t]/L) + (2*np.pi*my*ryi[t]/L))



	"*************************************************************************"
	"                              saving location                            "
	"*************************************************************************"
	Z = eta*np.sqrt(J/chi)
	zeta = np.sqrt((J*chi/(eta**2))*(1+((A/eta)*np.sqrt(chi/J))))
	p1_path = os.getcwd()
	folder = name
	os.mkdir(folder)
	n_path = os.path.join(p1_path,folder)
	os.chdir(n_path)
        tobe_copied1 = os.path.join(p1_path, "imageplot.sh")
        copied_to1 = os.path.join(n_path, "imageplot.sh")
        tobe_copied2 = os.path.join(p1_path, "plts.sh")
        copied_to2 = os.path.join(n_path, "plts.sh")
        shutil.copy(tobe_copied1, copied_to1)
        shutil.copy(tobe_copied2, copied_to2)
	"*************************************************************************"
	"                            details in text file                         "
	"*************************************************************************"
	text = open("parameters.txt", "w")
	text.write(str('flocks0.py'))
	text.write('Below are the parameters used\n')
	text.write("J,A = exchange interaction\n")
	text.write('J=' + str(J) + "\n")
	text.write('A=' + str(A) + "\n")
	text.write('Ac=' + str(Z) + "\n")
	text.write("chi = spin inertia\n")
	text.write('chi=' + str(chi) + "\n")
	text.write("# eta = viscosity")
	text.write('eta=' + str(eta) + "\n")
	text.write("Am = Amplitude of perturbation\n")
	text.write('Am = ' + str(Am) + '\n')
	text.write("dt = grid size of time\n")
	text.write('dt=' + str(dt) + "\n")
	text.write("# T = total time\n")
	text.write('T=' + str(T) + "\n")
	text.write("# N = flock size\n")
	text.write('N=' + str(N) + "\n")
	text.write("# L = length of a square box\n")
	text.write('L=' + str(L) + "\n")
	text.write("IR = interaction radius\n")
	text.write('IR=' + str(IR) + "\n")
	text.write("MX = L/wavelength of perturbation in X\n")
	text.write('MX=' + str(float(mx)) + "\n")
	text.write('MY = L/wavelength of perturbation in Y\n')
	text.write('MY=' + str(float(my)) + "\n")
	text.write("Zeta = crossover length\n")
	text.write('zeta=' + str(zeta) + "\n")
	text.write("f = random gaussian noise strength\n")
	text.write('f=' + str(f) + "\n")
	text.write("Alpha = strength\n")
	text.write('Alpha=' + str(Alpha) + "\n")
	text.close()

	"*************************************************************************"
	"                             time evolution                              "
	"*************************************************************************"
        order_parameter = []
	for g in range(G):
                for t in range(N):
			rndmnois = np.random.normal(0,1)
			sf[t] = si[t] + dt*(Alpha*interaction(t) + Beta*cohesion(t) - ((eta/chi)*si[t]) + f*rndmnois)
			thef[t] = thei[t] + (si[t]/chi)*dt
			vxf[t] = np.cos(thef[t])
			vyf[t] = np.sin(thef[t])
			rxf[t] = rxi[t] + vxf[t]*dt
			ryf[t] = ryi[t] + vyf[t]*dt
		"*********************************************************************"
		"                               periodicity                           "
		"*********************************************************************"
		for t in range(N):
			if rxf[t]<0:
				rxf[t] = rxf[t] + L
			if ryf[t]<0:
				ryf[t] = ryf[t] + L
			if rxf[t]>L:
				rxf[t] = rxf[t] - L
			if ryf[t]>L:
				ryf[t] = ryf[t] - L

		"*********************************************************************"
		for t in range(N):
			vxi[t] = vxf[t]
			vyi[t] = vyf[t]
			rxi[t] = rxf[t]
			ryi[t] = ryf[t]
			si[t] = sf[t]
			thei[t] = thef[t]

                "************************************************************************"
		"                                Saving data                             "
		"************************************************************************"
		dr = zip(rxi,ryi,vxi,vyi)
		nav = str(g)+'.dat'
		with open(nav, 'wb') as fp:
			a = csv.writer(fp, delimiter='\t')
			for x in dr:
				a.writerow(x)
		fp.close()
                "************************************************************************"
                "                   avg velocity, variance, avg_Spin                     "
                "************************************************************************"
                vx_aver = 0
                vy_aver = 0
                vx2_aver = 0
                vy2_aver = 0
                s_aver = 0
                s2_aver = 0
                
                for t in range(N):
                    vx_aver = vx_aver + vxi[t]
                    vy_aver = vy_aver + vyi[t]
                    vx2_aver = vx2_aver + vxi[t]**2
                    vy2_aver = vy2_aver + vyi[t]**2
                    s_aver = s_aver + si[t]
                    s2_aver = s2_aver +si[t]**2

                vx_average = (1/float(N))*vx_aver
                vy_average = (1/float(N))*vy_aver
                vx2_average = (1/float(N))*vx2_aver
                vy2_average = (1/float(N))*vy2_aver
                s_average = (1/float(N))*s_aver
                s2_average = (1/float(N))*s2_aver
                
                v_average = np.sqrt(vx_average**2 + vy_average**2)
                vx_variance = vx2_average - vx_average**2
                vy_variance = vy2_average - vy_average**2
                s_variance = s2_average - s_average**2
                
                if vx_average != 0:
                    orientation = (180/np.pi)*np.arctan(vy_average/vx_average)
                else:
                    orientation = 90
                timestep = g*dt
                order_parameter.append([timestep, v_average, orientation, vx_average, vy_average, vx_variance, vy_variance, s_average, s_variance])
                
        ord_param = np.matrix(order_parameter)
        df_op = pd.DataFrame(data=ord_param.astype(float))
        df_op.to_csv("order_parameter.dat", sep=' ', header=False, float_format='%.15f', index=False)

        "*********************************************************************************"
        subprocess.call(['bash', 'imageplot.sh'])
        "*********************************************************************************"
        subprocess.call(['bash', 'plts.sh'])
        "*********************************************************************************"
	os.chdir(p1_path)
        "*********************************************************************************"


if __name__=="__main__":
    chi = 10.0
    eta = 5.0
    Am = 1*np.pi/180
    dt = 0.05
    T = 80.0
    L = 40
    IR = 1.2
    J = 10
    F = 0
    mu = np.pi/6
    MX = 1.0
    MY = 1.0
    k = eta*np.sqrt(J/chi)
    Asy = [0.6*k, 0.7*k, 1*k, 1.1*k]
    beta = [1.0,2.0,5.0] #[2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
    Alps = [1.0,2.0,5.0] #[2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
    Arr = []
    for Beta in beta:
        for Alpha in Alps:
            for A in Asy:
                name ="Alp" + str(Alpha) + "Bet" + str(Beta) +  "nis" + str(F) + "A" + str(A)
                arr = [J,A,chi,eta,Am,dt,T,L,IR,MX,MY,mu,F,Alpha,Beta,name]
                Arr.append(arr)
    
    p = Pool()
    p.map(exchange,Arr)
    #exchange(arr)
