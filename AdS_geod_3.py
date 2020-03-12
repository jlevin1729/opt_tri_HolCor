#!/usr/bin/env python3
import numpy as np
import math as m
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import scipy.optimize as opt
import random as rand

#start with three non-overlapping boundary regions and an optimized correlation measure

CM = 6

theta_A1 = 0
theta_A2 = m.pi/3

theta_B1 = 2*m.pi/3
theta_B2 = m.pi

theta_C1 = 4*m.pi/3
theta_C2 = 5*m.pi/3

x_A1 = m.cos(theta_A1)
y_A1 = m.sin(theta_A1)
x_A2 = m.cos(theta_A2)
y_A2 = m.sin(theta_A2)

x_B1 = m.cos(theta_B1)
y_B1 = m.sin(theta_B1)
x_B2 = m.cos(theta_B2)
y_B2 = m.sin(theta_B2)

x_C1 = m.cos(theta_C1)
y_C1 = m.sin(theta_C1)
x_C2 = m.cos(theta_C2)
y_C2 = m.sin(theta_C2)


theta_A = theta_A2 - theta_A1
theta_B = theta_B2 - theta_B1
theta_C = theta_C2 - theta_C1

theta_AB = theta_B1 - theta_A2
theta_BC = theta_C1 - theta_B2
theta_CA = theta_A1 - theta_C2 + 2*m.pi

#set cut-off and define distance function

def distance(x_1, y_1, x_2, y_2):
	
	cutoff = 0.000000000000001
	max_r = 1 - cutoff

	if m.sqrt(x_1**2 + y_1**2) > max_r:	

		x_1 = x_1*(max_r/m.sqrt(x_1**2 + y_1**2))	
		y_1 = y_1*(max_r/m.sqrt(x_1**2 + y_1**2))
	
	if m.sqrt(x_2**2 + y_2**2) > max_r:

		x_2 = x_2*(max_r/m.sqrt(x_2**2 + y_2**2))	
		y_2 = y_2*(max_r/m.sqrt(x_2**2 + y_2**2))
	
	num = (x_2 - x_1)**2 + (y_2 - y_1)**2
	den = (1 - x_1**2 - y_1**2)*(1 - x_2**2 - y_2**2)
	delta = 2*num/den
	return m.acosh(1+delta)

#make a figure

fig, ax = plt.subplots()

ax.set_xlim((-1,1))
ax.set_ylim((-1,1))

#draw boundary

boundary = plt.Circle((0,0), 1, color = 'black', fill = False, clip_on = False)
ax.add_artist(boundary)

#plot boundary region endpoints

ax.plot((x_A1), (y_A1), 'o', color = 'black')
ax.plot((x_A2), (y_A2), 'o', color = 'black')
ax.plot((x_B1), (y_B1), 'o', color = 'black')
ax.plot((x_B2), (y_B2), 'o', color = 'black')
ax.plot((x_C1), (y_C1), 'o', color = 'black')
ax.plot((x_C2), (y_C2), 'o', color = 'black')

ax.set_aspect('equal')

#define function that draws the geodesic between two given points, with given color and thickness

def geod(x_1, y_1, x_2, y_2, n):

	theta_ = m.atan2(y_2,x_2) - m.atan2(y_1,x_1)
	
	if theta_ > m.pi:

		X_1 = x_2
		Y_1 = y_2
		X_2 = x_1
		Y_2 = y_1
		theta = 2*m.pi - theta_

	elif theta_ < -m.pi:

		X_1 = x_1
		Y_1 = y_1
		X_2 = x_2
		Y_2 = y_2
		theta = 2*m.pi - abs(theta_)

	elif 0 < theta_ < m.pi:
		
		X_1 = x_1
		Y_1 = y_1
		X_2 = x_2
		Y_2 = y_2
		theta = theta_

	else:

		X_1 = x_2
		Y_1 = y_2
		X_2 = x_1
		Y_2 = y_1
		theta = abs(theta_)

	if theta == m.pi:

		if n == 1:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 1, color = 'g')

		if n == 2:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 4, color = 'g')

		if n == 3:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 7, color = 'g')

		if n == 4:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 10, color = 'g')

		if n == -1:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 1, color = 'r')

		if n == -2:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 4, color = 'r')

		if n == -3:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 7, color = 'r')

		if n == -4:

			plt.plot((x_1, x_2), (y_1, y_2), linewidth = 10, color = 'r')

	else:

		

		c_x = -(Y_1*(1 + X_2**2 + Y_2**2) - Y_2*(1 + X_1**2 + Y_1**2))/(2*(m.sqrt(X_1**2 + Y_1**2)*m.sqrt(X_2**2 + Y_2**2))*m.sin(theta))
		c_y = (X_1*(1 + X_2**2 + Y_2**2) - X_2*(1 + X_1**2 + Y_1**2))/(2*(m.sqrt(X_1**2 + Y_1**2)*m.sqrt(X_2**2 + Y_2**2))*m.sin(theta))
		R = m.sqrt((c_x - x_1)**2 + (c_y - y_1)**2)

		theta_2 = m.atan2(Y_2 - c_y, X_2 - c_x)*180/m.pi
		theta_1 = m.atan2(Y_1 - c_y, X_1 - c_x)*180/m.pi

	
		if n == 1:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 1, color = 'g')

		if n == 2:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 4, color = 'g')

		if n == 3:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 7, color = 'g')

		if n == 4:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 10, color = 'g')

		if n == -1:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 1, color = 'r')

		if n == -2:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 4, color = 'r')

		if n == -3:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 7, color = 'r')

		if n == -4:

			geod = pat.Arc((c_x,c_y), 2*R, 2*R, angle = 0, theta1 = theta_2, theta2 = theta_1, linewidth = 10, color = 'r')


		ax.add_patch(geod)
	return

#draw entanglement wedge

EW = distance(x_A2, y_A2, x_B1, y_B1) + distance(x_B2, y_B2, x_C1, y_C1) + distance(x_C2, y_C2, x_A1, y_A1) 

if EW <= distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x_C2, y_C2) and EW <= distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B2, y_B2, x_C1, y_C1) + distance(x_B1, y_B1, x_C2, y_C2) and EW <= distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C2, y_C2, x_A1, y_A1) + distance(x_C1, y_C1, x_A2, y_A2) and EW <= distance(x_C1, y_C1, x_C2, y_C2) + distance(x_A2, y_A2, x_B1, y_B1) + distance(x_A1, y_A1, x_B2, y_B2):

	print('I(A:B) = ', max(0, distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x_B2, y_B2) - distance(x_A2, y_A2, x_B1, y_B1) - distance(x_A1, y_A1, x_B2, y_B2)))
	print('I(B:C) = ', max(0, distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x_C2, y_C2) - distance(x_B2, y_B2, x_C1, y_C1) - distance(x_B1, y_B1, x_C2, y_C2)))
	print('I(A:C) = ', max(0, distance(x_A1, y_A1, x_A2, y_A2) + distance(x_C1, y_C1, x_C2, y_C2) - distance(x_A2, y_A2, x_C1, y_C1) - distance(x_A1, y_A1, x_C2, y_C2)))

	c_x1 = -(y_A2*2 - y_B1*2)/(2*m.sin(theta_AB))
	c_y1 = (x_A2*2 - x_B1*2)/(2*m.sin(theta_AB))
	R1 = m.sqrt((c_x1 - x_A2)**2 + (c_y1 - y_A2)**2)

	c_x2 = -(y_B2*2 - y_C1*2)/(2*m.sin(theta_BC))
	c_y2 = (x_B2*2 - x_C1*2)/(2*m.sin(theta_BC))
	R2 = m.sqrt((c_x2 - x_B2)**2 + (c_y2 - y_B2)**2)

	c_x3 = -(y_C2*2 - y_A1*2)/(2*m.sin(theta_CA))
	c_y3 = (x_C2*2 - x_A1*2)/(2*m.sin(theta_CA))
	R3 = m.sqrt((c_x3 - x_C2)**2 + (c_y3 - y_C2)**2)

	geod1 = pat.Arc((c_x1,c_y1), 2*R1, 2*R1, angle = 0, theta1 = m.atan2(y_B1 - c_y1, x_B1 - c_x1)*180/m.pi, theta2 = m.atan2(y_A2 - c_y1, x_A2 - c_x1)*180/m.pi, linestyle = '--')
	geod2 = pat.Arc((c_x2,c_y2), 2*R2, 2*R2, angle = 0, theta1 = m.atan2(y_C1 - c_y2, x_C1 - c_x2)*180/m.pi, theta2 = m.atan2(y_B2 - c_y2, x_B2 - c_x2)*180/m.pi, linestyle = '--')
	geod3 = pat.Arc((c_x3,c_y3), 2*R3, 2*R3, angle = 0, theta1 = m.atan2(y_A1 - c_y3, x_A1 - c_x3)*180/m.pi, theta2 = m.atan2(y_C2 - c_y3, x_C2 - c_x3)*180/m.pi, linestyle = '--')
	

	ax.add_patch(geod1)
	ax.add_patch(geod2)
	ax.add_patch(geod3)

#find points where EWCS is achieved

	def f_P(x):

		x1 = c_x1 + R1*m.cos(x[0])
		y1 = c_y1 + R1*m.sin(x[0])
		x2 = c_x2 + R2*m.cos(x[1])
		y2 = c_y2 + R2*m.sin(x[1])
		x3 = c_x3 + R3*m.cos(x[2])
		y3 = c_y3 + R3*m.sin(x[2])

		S_Aa = distance(x1, y1, x3, y3)
		S_Bb = distance(x1, y1, x2, y2)
		S_Cc = distance(x2, y2, x3, y3)

		f_P = S_Aa + S_Bb + S_Cc

		return f_P

	E_P = opt.minimize(f_P, [theta_A2 + theta_AB/2 - m.pi, theta_B2 + theta_BC/2 - m.pi, theta_C2 + theta_CA/2 - m.pi])
	#, bounds = opt.Bounds([theta_A2 + theta_AB - 3*m.pi/2, theta_A2 - m.pi/2], [theta_B2 + theta_BC - 3*m.pi/2, theta_B2 - m.pi/2], [theta_C2 + theta_CA - 3*m.pi/2, theta_C2 - m.pi/2])
	
	print('E_P =', E_P.fun)
	print(E_P.x)
	
		
#define objective function to be minimized

	def f_CM(x):

		x1 = c_x1 + R1*m.cos(x[0])
		y1 = c_y1 + R1*m.sin(x[0])
		x2 = c_x2 + R2*m.cos(x[1])
		y2 = c_y2 + R2*m.sin(x[1])
		x3 = c_x3 + R3*m.cos(x[2])
		y3 = c_y3 + R3*m.sin(x[2])

		#S_A-------------------------------

		S_A = distance(x_A1, y_A1, x_A2, y_A2)

		#S_B-------------------------------

		S_B = distance(x_B1, y_B1, x_B2, y_B2)

		#S_C-------------------------------

		S_C = distance(x_C1, y_C1, x_C2, y_C2)

		#S_a-------------------------------

		S_a = distance(x_A1, y_A1, x3, y3) + distance(x_A2, y_A2, x1, y1)

		#S_b-------------------------------

		S_b = distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x2, y2)	
	
		#S_AB------------------------------	

		S_AB = min(distance(x_A1, y_A1, x_B2, y_B2) + distance(x_A2, y_A2, x_B1, y_B1), distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x_B2, y_B2))

		#S_AC------------------------------

		S_AC = min(distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x_C1, y_C1), distance(x_A1, y_A1, x_A2, y_A2) + distance(x_C1, y_C1, x_C2, y_C2))

		#S_Aa------------------------------

		S_Aa = distance(x1, y1, x3, y3)

		#S_Ab------------------------------

		S_Ab = min(distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x2, y2), distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x2, y2) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x_B1, y_B1) + distance(x_B2, y_B2, x2, y2), distance(x_A2, y_A2, x_B2, y_B2) + distance(x_A1, y_A1, x2, y2) + distance(x_B1, y_B1, x1, y1))

		#S_BC------------------------------

		S_BC = min(distance(x_B1, y_B1, x_C2, y_C2) + distance(x_B2, y_B2, x_C1, y_C1), distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x_C2, y_C2))

		#S_Ba------------------------------

		S_Ba = min(distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x3, y3) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x3, y3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_B1, y_B1, x1, y1) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_A1, y_A1, x3, y3), distance(x_A1, y_A1, x_B1, y_B1) + distance(x_B2, y_B2, x3, y3) + distance(x_A2, y_A2, x1, y1))

		#S_Bb------------------------------

		S_Bb = distance(x1, y1, x2, y2)

		#S_Ca------------------------------

		S_Ca = min(distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x3, y3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x1, y1))

		#S_Cb------------------------------

		S_Cb = min(distance(x_C1, y_C1, x_C2, y_C2) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x2, y2), distance(x_C2, y_C2, x1, y1) + distance(x_C1, y_C1, x2, y2) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_C2, y_C2, x1, y1) + distance(x_C1, y_C1, x_B1, y_B1) + distance(x_B2, y_B2, x2, y2), distance(x_C2, y_C2, x_B2, y_B2) + distance(x_C1, y_C1, x2, y2) + distance(x_B1, y_B1, x1, y1))

		#S_ab------------------------------

		S_ab = distance(x_A2, y_A2, x_B1, y_B1) + distance(x_A1, y_A1, x3, y3) + distance(x_B2, y_B2, x2, y2)

		#S_ABC-----------------------------

		S_ABC = distance(x_A2, y_A2, x_B1, y_B1) + distance(x_B2, y_B2, x_C1, y_C1) + distance(x_C2, y_C2, x_A1, y_A1)

		#S_ABa-----------------------------

		S_ABa = min(distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x3, y3), distance(x1, y1, x3, y3) + distance(x_B1, y_B1, x_B2, y_B2))

		#S_ABb-----------------------------

		S_ABb = min(distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x2, y2), distance(x1, y1, x2, y2) + distance(x_A1, y_A1, x_A2, y_A2))

		#S_ACa-----------------------------

		S_ACa = min(distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x3, y3), distance(x1, y1, x3, y3) + distance(x_C1, y_C1, x_C2, y_C2))

		#S_ACb-----------------------------

		S_ACb = min(distance(x_C2, y_C2, x_A1, y_A1) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x2, y2), distance(x_A2, y_A2, x_A1, y_A1) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2) + distance(x_B1, y_B1, x1, y1), distance(x_A1, y_A1, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2), distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2), distance(x_A1, y_A1, x2, y2) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x1, y1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_B1, y_B1), distance(x_A2, y_A2, x_A1, y_A1) + distance(x_C2, y_C2, x1, y1) + distance(x_B1, y_B1, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2))

		#S_Aab-----------------------------

		S_Aab = min(distance(x_B1, y_B1, x_B2, y_B2) + distance(x2, y2, x3, y3), distance(x_B1, y_B1, x3, y3) + distance(x_B2, y_B2, x2, y2))

		#S_BCa-----------------------------

		S_BCa = min(distance(x_B2, y_B2, x_C1, y_C1) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x1, y1), distance(x_C2, y_C2, x_C1, y_C1) + distance(x_A1, y_A1, x3, y3) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_B2, y_B2) + distance(x_C2, y_C2, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x3, y3), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1), distance(x_C1, y_C1, x_C2, y_C2) + distance(x_A1, y_A1, x3, y3) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1), distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_B2, y_B2) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_B2, y_B2) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_A1, y_A1), distance(x_C2, y_C2, x_C1, y_C1) + distance(x_B2, y_B2, x3, y3) + distance(x_A1, y_A1, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1))

		#S_BCb-----------------------------

		S_BCb = min(distance(x_C1, y_C1, x_C2, y_C2) + distance(x2, y2, x1, y1), distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x1, y1))

		#S_Bab-----------------------------

		S_Bab = min(distance(x_A1, y_A1, x_A2, y_A2) + distance(x2, y2, x3, y3), distance(x_A1, y_A1, x3, y3) + distance(x_A2, y_A2, x2, y2))

		#S_Cab-----------------------------

		S_Cab = min(distance(x_A2, y_A2, x_B1, y_B1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_C2, y_C2) + distance(x_A1, y_A1, x3, y3), distance(x_B2, y_B2, x_B1, y_B1) + distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_A2, y_A2) + distance(x_B2, y_B2, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_C1, y_C1, x2, y2), distance(x_B1, y_B1, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3), distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3), distance(x_B1, y_B1, x3, y3) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_A2, y_A2) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_A2, y_A2) + distance(x_B2, y_B2, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_C1, y_C1), distance(x_B2, y_B2, x_B1, y_B1) + distance(x_A2, y_A2, x2, y2) + distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3))

		#S_ABCa----------------------------

		S_ABCa = distance(x_B2, y_B2, x_C1, y_C1) + distance(x_B1, y_B1, x1, y1) + distance(x_C2, y_C2, x3, y3)

		#S_ABCb----------------------------

		S_ABCb = distance(x_C2, y_C2, x_A1, y_A1) + distance(x_A2, y_A2, x1, y1) + distance(x_C1, y_C1, x2, y2)

		#S_ABab----------------------------

		S_ABab = distance(x2, y2, x3, y3)

		#S_ACab----------------------------
	
		S_ACab = min(distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_B2, y_B2, x2, y2) + distance(x_B1, y_B1, x3, y3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_B2, y_B2, x2, y2) + distance(x_B1, y_B1, x_C1, y_C1) + distance(x_C2, y_C2, x3, y3), distance(x_B2, y_B2, x_C2, y_C2) + distance(x_B1, y_B1, x3, y3) + distance(x_C1, y_C1, x2, y2))

		#S_BCab----------------------------

		S_BCab = min(distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_A2, y_A2, x2, y2) + distance(x_A1, y_A1, x3, y3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_A2, y_A2, x2, y2) + distance(x_A1, y_A1, x_C1, y_C1) + distance(x_C2, y_C2, x3, y3), distance(x_A2, y_A2, x_C2, y_C2) + distance(x_A1, y_A1, x3, y3) + distance(x_C1, y_C1, x2, y2))

		#S_ABCab---------------------------

		S_ABCab = distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x3, y3)

		
		
		S_BCabc = S_A

		S_ACabc = S_B

		S_ABabc = S_C

		S_ABCbc = S_a

		S_ABCac = S_b

		S_Cabc = S_AB

		S_Babc = S_AC

		S_BCbc = S_Aa

		S_BCac = S_Ab

		S_Aabc = S_BC

		S_ACbc = S_Ba

		S_ACac = S_Bb

		S_ABbc = S_Ca 

		S_ABac = S_Cb  

		S_ABCc = S_ab  

		S_abc = S_ABC  

		S_Cbc = S_ABa  

		S_Cac = S_ABb  

		S_Bbc = S_ACa  

		S_Bac = S_ACb  

		S_BCc = S_Aab  

		S_Abc = S_BCa  

		S_Aac = S_BCb  

		S_ACc = S_Bab  

		S_ABc = S_Cab  

		S_bc = S_ABCa 

		S_ac = S_ABCb  

		S_Cc = S_ABab  
	
		S_Bc = S_ACab  

		S_Ac = S_BCab  

		S_c = S_ABCab 

		

		alpha = S_Aa + S_Bb + S_Cc

		beta = S_A + S_B + S_C 

		gamma = S_AB + S_AC + S_BC

		delta = S_ABC

		epsilon = S_ABc + S_ACb + S_BCa

		eta = S_Ab + S_Ac + S_Ba + S_Bc + S_Ca + S_Cb

		theta = S_ABa + S_ABb + S_ACa + S_ACc + S_BCb + S_BCc

		omega = S_ab + S_ac + S_bc

		rho = S_a + S_b + S_c


		if CM == 1:

			f_CM = 2*beta - gamma

		elif CM == 2:

			f_CM = 2*beta + theta - eta

		elif CM == 3:

			f_CM = gamma - 2*delta

		elif CM == 4:

			f_CM = alpha

		elif CM == 5:

			f_CM = theta - omega

		elif CM == 6:

			f_CM = eta + theta - 2*omega - 2*rho

		elif CM == 7:

			f_CM = eta + theta + beta - omega - 2*rho - epsilon

		elif CM == 8:

			f_CM = theta - 2*rho

		elif CM == 9:

			f_CM = theta + beta + gamma - epsilon - omega

		elif CM == 10:

			f_CM = theta + epsilon - rho - 2*omega

		elif CM == 11:

			f_CM = theta + 2*delta - 2*omega

		elif CM == 12:

			f_CM = 2*alpha + eta - theta - 2*rho

		elif CM == 13:

			f_CM = 2*alpha + gamma - theta

		elif CM == 14:

			f_CM = alpha + 2*delta - omega

		elif CM == 15:

			f_CM = 2*alpha + 2*beta + epsilon - eta - rho

		elif CM == 16:

			f_CM = alpha + epsilon - omega - rho

		elif CM == 17:

			f_CM = alpha + beta + gamma - epsilon

		elif CM == 18:

			f_CM = alpha + omega - 2*rho

		else:

			f_CM = alpha + beta + eta - epsilon - 2*rho


		return f_CM
		
	
#minimize the objective function f_CM, plot the minimizing points	

	#for i in range(10):	

		#A = [rand.randint(450, 550), rand.randint(450, 550), rand.randint(450, 550)]
		
		#theta_1 = m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi) + (A[0]/1000)*((m.atan2(y_A2 - c_y1, x_A2 - c_x1)%(2*m.pi) - m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi))%(2*m.pi))
		#theta_2 = m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi) + (A[1]/1000)*((m.atan2(y_B2 - c_y2, x_B2 - c_x2)%(2*m.pi) - m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi))%(2*m.pi))
		#theta_3 = m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi) + (A[2]/1000)*((m.atan2(y_C2 - c_y3, x_C2 - c_x3)%(2*m.pi) - m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi))%(2*m.pi))
		
		#ax.plot(c_x1 + R1*m.cos(theta_1), c_y1 + R1*m.sin(theta_1), 'o', color = 'black')
		#ax.plot(c_x2 + R2*m.cos(theta_2), c_y2 + R2*m.sin(theta_2), 'o', color = 'black')
		#ax.plot(c_x3 + R3*m.cos(theta_3), c_y3 + R3*m.sin(theta_3), 'o', color = 'black')

		#lb_1 = m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi)
		#ub_1 = m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi) + (m.atan2(y_A2 - c_y1, x_A2 - c_x1)%(2*m.pi) - m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi))%(2*m.pi)

		#lb_2 = m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi)
		#ub_2 = m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi) + (m.atan2(y_B2 - c_y2, x_B2 - c_x2)%(2*m.pi) - m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi))%(2*m.pi)

		#lb_3 = m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi)
		#ub_3 = m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi) + (m.atan2(y_C2 - c_y3, x_C2 - c_x3)%(2*m.pi) - m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi))%(2*m.pi)

		#f__CM = opt.minimize(f_CM, [theta_1, theta_2, theta_3])
		#, bounds = ((lb_1,ub_1), (lb_2,ub_2), (lb_3,ub_3))

		#if i == 0:

			#E_CM = f__CM

		#else:
	
			#if f__CM.fun < E_CM.fun:
			
				#E_CM = f__CM

	#lb_1 = m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi)
	#ub_1 = m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi) + (m.atan2(y_A2 - c_y1, x_A2 - c_x1)%(2*m.pi) - m.atan2(y_B1 - c_y1, x_B1 - c_x1)%(2*m.pi))%(2*m.pi)

	#lb_2 = m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi)
	#ub_2 = m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi) + (m.atan2(y_B2 - c_y2, x_B2 - c_x2)%(2*m.pi) - m.atan2(y_C1 - c_y2, x_C1 - c_x2)%(2*m.pi))%(2*m.pi)

	#lb_3 = m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi)
	#ub_3 = m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi) + (m.atan2(y_C2 - c_y3, x_C2 - c_x3)%(2*m.pi) - m.atan2(y_A1 - c_y3, x_A1 - c_x3)%(2*m.pi))%(2*m.pi)			
		
	E_CM = opt.minimize(f_CM, [theta_A2 + theta_AB/2 - m.pi, theta_B2 + theta_BC/2 - m.pi, theta_C2 + theta_CA/2 - m.pi])

	#initial guess E_P.x
	#, bounds = opt.Bounds([theta_A2 + theta_AB - 3*m.pi/2, theta_A2 - m.pi/2], [theta_B2 + theta_BC - 3*m.pi/2, theta_B2 - m.pi/2], [theta_C2 + theta_CA - 3*m.pi/2, theta_C2 - m.pi/2])
	#, bounds = ((lb_1,ub_1), (lb_2,ub_2), (lb_3,ub_3))

	print(E_CM.x)
	print('E_CM =', E_CM.fun)


	x_CM1 = c_x1 + R1*m.cos(E_CM.x[0])
	y_CM1 = c_y1 + R1*m.sin(E_CM.x[0])
	x_CM2 = c_x2 + R2*m.cos(E_CM.x[1])
	y_CM2 = c_y2 + R2*m.sin(E_CM.x[1])
	x_CM3 = c_x3 + R3*m.cos(E_CM.x[2])
	y_CM3 = c_y3 + R3*m.sin(E_CM.x[2])

	#print('S(Ca) = ', min(distance(x_A2, y_A2, x_CM1, y_CM1) + distance(x_A1, y_A1, x_CM3, y_CM3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_C1, y_C1, x_CM1, y_CM1) + distance(x_C2, y_C2, x_CM3, y_CM3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_C1, y_C1, x_CM1, y_CM1) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x_CM3, y_CM3), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x_CM3, y_CM3) + distance(x_A2, y_A2, x_CM1, y_CM1)))

	#print(distance(x_A2, y_A2, x_CM1, y_CM1) + distance(x_A1, y_A1, x_CM3, y_CM3) + distance(x_C2, y_C2, x_C1, y_C1))

	#print('Surfaces Ca = ', np.argmin((distance(x_A2, y_A2, x_CM1, y_CM1) + distance(x_A1, y_A1, x_CM3, y_CM3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_C1, y_C1, x_CM1, y_CM1) + distance(x_C2, y_C2, x_CM3, y_CM3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_C1, y_C1, x_CM1, y_CM1) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x_CM3, y_CM3), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x_CM3, y_CM3) + distance(x_A2, y_A2, x_CM1, y_CM1))))

	#print('S(Cb) = ', min(distance(x_C1, y_C1, x_C2, y_C2) + distance(x_B1, y_B1, x_CM1, y_CM1) + distance(x_B2, y_B2, x_CM2, y_CM2), distance(x_C2, y_C2, x_CM1, y_CM1) + distance(x_C1, y_C1, x_CM2, y_CM2) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_C2, y_C2, x_CM1, y_CM1) + distance(x_C1, y_C1, x_B1, y_B1) + distance(x_B2, y_B2, x_CM2, y_CM2), distance(x_C2, y_C2, x_B2, y_B2) + distance(x_C1, y_C1, x_CM2, y_CM2) + distance(x_B1, y_B1, x_CM1, y_CM1)))

	#print(distance(x_C1, y_C1, x_C2, y_C2) + distance(x_B1, y_B1, x_CM1, y_CM1) + distance(x_B2, y_B2, x_CM2, y_CM2))

	#print('Surfaces Cb = ', np.argmin((distance(x_C1, y_C1, x_C2, y_C2) + distance(x_B1, y_B1, x_CM1, y_CM1) + distance(x_B2, y_B2, x_CM2, y_CM2), distance(x_C2, y_C2, x_CM1, y_CM1) + distance(x_C1, y_C1, x_CM2, y_CM2) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_C2, y_C2, x_CM1, y_CM1) + distance(x_C1, y_C1, x_B1, y_B1) + distance(x_B2, y_B2, x_CM2, y_CM2), distance(x_C2, y_C2, x_B2, y_B2) + distance(x_C1, y_C1, x_CM2, y_CM2) + distance(x_B1, y_B1, x_CM1, y_CM1))))

	
	
	ax.plot((x_CM1), (y_CM1), 'o', color = 'black')
	ax.plot((x_CM2), (y_CM2), 'o', color = 'black')
	ax.plot((x_CM3), (y_CM3), 'o', color = 'black')

#define function which gives bulk surfaces contributing to f_CM (order = (A1-A2, A1-B1, A1-B2, A1-C1, A1-1, A1-2, A1-3, A2-B2, A2-C1, A2-C2, A2-1, A2-2, A2-3, B1-B2, B1-C1, B1-C2, B1-1, B1-2, B1-3, B2-C2, B2-1, B2-2, B2-3, C1-C2, C1-1, C1-2, C1-3, C2-1, C2-2, C2-3, 1-2, 1-3, 2-3))

	def surf_CM(x):
		
		x1 = c_x1 + R1*m.cos(x[0])
		y1 = c_y1 + R1*m.sin(x[0])
		x2 = c_x2 + R2*m.cos(x[1])
		y2 = c_y2 + R2*m.sin(x[1])
		x3 = c_x3 + R3*m.cos(x[2])
		y3 = c_y3 + R3*m.sin(x[2])


		v0 = np.array([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v1 = np.array([0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v2 = np.array([0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v3 = np.array([0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v4 = np.array([0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v5 = np.array([0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v6 = np.array([0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v7 = np.array([0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v8 = np.array([0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v9 = np.array([0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v10 = np.array([0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v11 = np.array([0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v12 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v13 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v14 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v15 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v16 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v17 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v18 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v19 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
		v20 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0])
		v21 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0])
		v22 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0])
		v23 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0])
		v24 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0])
		v25 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0])
		v26 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0])
		v27 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0])
		v28 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0])
		v29 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0])
		v30 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0])
		v31 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0])
		v32 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])

		#A--------------------------------			

		A = v0

		#B--------------------------------

		B = v13

		#C--------------------------------

		C = v23

		#a--------------------------------

		a = v6 + v10

		#b--------------------------------

		b = v16 + v21	

		#AB-------------------------------

		surf_AB = np.argmin((distance(x_A1, y_A1, x_B2, y_B2) + distance(x_A2, y_A2, x_B1, y_B1), distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x_B2, y_B2)))

		if surf_AB == 0:

			AB = v2 + v10 + v16

		else:

			AB = v0 + v13

		#AC-------------------------------

		surf_AC = np.argmin((distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x_C1, y_C1), distance(x_A1, y_A1, x_A2, y_A2) + distance(x_C1, y_C1, x_C2, y_C2)))

		if surf_AC == 0:

			AC = v6 + v29 + v8

		else:

			AC = v0 + v23

		#Aa-------------------------------

		Aa = v31

		#Ab-------------------------------

		surf_Ab = np.argmin((distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x2, y2), distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x2, y2) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x_B1, y_B1) + distance(x_B2, y_B2, x2, y2), distance(x_A2, y_A2, x_B2, y_B2) + distance(x_A1, y_A1, x2, y2) + distance(x_B1, y_B1, x1, y1)))		

		if surf_Ab == 0:

			Ab = v0 + v16 + v21

		elif surf_Ab == 1:

			Ab = v10 + v5 + v13

		elif surf_Ab == 2:

			Ab = v10 + v1 + v21

		else:

			Ab = v7 + v5 + v16

		#BC-------------------------------

		surf_BC = np.argmin((distance(x_B1, y_B1, x_C2, y_C2) + distance(x_B2, y_B2, x_C1, y_C1), distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x_C2, y_C2)))		

		if surf_BC == 0:	

			BC = v15 + v21 + v25

		else:

			BC = v13 + v23

		#Ba-------------------------------

		surf_Ba = np.argmin((distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x3, y3) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x3, y3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_B1, y_B1, x1, y1) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_A1, y_A1, x3, y3), distance(x_A1, y_A1, x_B1, y_B1) + distance(x_B2, y_B2, x3, y3) + distance(x_A2, y_A2, x1, y1)))

		if surf_Ba == 0:

			Ba = v10 + v6 + v13

		elif surf_Ba == 1:

			Ba = v16 + v22 +v0

		elif surf_Ba == 2:

			Ba = v16 + v7 + v6

		else: 

			Ba = v1 + v22 + v10

		#Bb--------------------------------

		Bb = v30

		#Ca--------------------------------

		surf_Ca = np.argmin((distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x3, y3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x1, y1)))		

		if surf_Ca == 0:

			Ca = v10 + v6 + v23

		elif surf_Ca == 1:

			Ca = v24 + v29 + v0

		elif surf_Ca == 2:

			Ca = v24 + v9 + v6

		else:

			Ca = v3 + v29 + v10

		#Cb--------------------------------

		surf_Cb = np.argmin((distance(x_C1, y_C1, x_C2, y_C2) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x2, y2), distance(x_C2, y_C2, x1, y1) + distance(x_C1, y_C1, x2, y2) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_C2, y_C2, x1, y1) + distance(x_C1, y_C1, x_B1, y_B1) + distance(x_B2, y_B2, x2, y2), distance(x_C2, y_C2, x_B2, y_B2) + distance(x_C1, y_C1, x2, y2) + distance(x_B1, y_B1, x1, y1)))		

		if surf_Cb == 0:

			Cb = v23 + v16 + v21

		elif surf_Cb == 1:

			Cb = v27 + v25 + v13	

		elif surf_Cb == 2:

			Cb = v14 + v27 + v21

		else:

			Cb = v19 + v25 + v16

		#ab--------------------------------

		ab = v6 + v21 + v10 + v16

		#ABC-------------------------------

		ABC = v10 + v16 + v21 + v25 + v6 + v29

		#ABa-------------------------------

		surf_ABa = np.argmin((distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x3, y3), distance(x1, y1, x3, y3) + distance(x_B1, y_B1, x_B2, y_B2)))		

		if surf_ABa == 0:	

			ABa = v16 + v22

		else:

			ABa = v13 + v31

		#ABb-------------------------------
		
		surf_ABb = np.argmin((distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x2, y2), distance(x1, y1, x2, y2) + distance(x_A1, y_A1, x_A2, y_A2))) 		

		if surf_ABb == 0:

			ABb = v5 + v10

		else:

			ABb = v0 + v30

		#ACa-------------------------------

		surf_ACa = np.argmin((distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x3, y3), distance(x1, y1, x3, y3) + distance(x_C1, y_C1, x_C2, y_C2))) 		

		if surf_ACa == 0:	

			ACa = v24 + v29

		else:

			ACa = v23 + v31

		#ACb-------------------------------

		surf_ACb = np.argmin((distance(x_C2, y_C2, x_A1, y_A1) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x2, y2), distance(x_A2, y_A2, x_A1, y_A1) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2) + distance(x_B1, y_B1, x1, y1), distance(x_A1, y_A1, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2), distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x1, y1) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2), distance(x_A1, y_A1, x2, y2) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1) + distance(x_C1, y_C1, x_C2, y_C2), distance(x_A1, y_A1, x_C2, y_C2) + distance(x_A2, y_A2, x1, y1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_B1, y_B1), distance(x_A2, y_A2, x_A1, y_A1) + distance(x_C2, y_C2, x1, y1) + distance(x_B1, y_B1, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2)))

		if surf_ACb == 0:

			ACb = v6 + v29 + v10 + v13 + v25

		elif surf_ACb == 1:

			ACb = v0 + v16 + v21 + v23

		elif surf_ACb == 2:

			ACb = v6 + v29 + v8 + v21 + v16

		elif surf_ACb == 3:

			ACb = v1 + v10 + v21 + v23

		elif surf_ACb == 4:

			ACb = v1 + v10 + v19 + v25

		elif surf_ACb == 5:

			ACb = v0 + v16 + v19 + v25

		elif surf_ACb == 6:

			ACb = v5 + v7 + v16 + v23

		elif surf_ACb == 7:

			ACb = v6 + v29 + v7 + v16 + v23

		elif surf_ACb == 8:

			ACb = v6 + v29 + v10 + v21 + v14

		else:

			ACb = v0 + v27 + v14 + v21

		
		#Aab-------------------------------

		surf_Aab = np.argmin((distance(x_B1, y_B1, x_B2, y_B2) + distance(x2, y2, x3, y3), distance(x_B1, y_B1, x3, y3) + distance(x_B2, y_B2, x2, y2))) 		

		if surf_Aab == 0:

			Aab = v13 + v32

		else:

			Aab = v18 + v21

		#BCa-------------------------------

		surf_BCa = np.argmin((distance(x_B2, y_B2, x_C1, y_C1) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_A2, y_A2) + distance(x_B1, y_B1, x1, y1), distance(x_C2, y_C2, x_C1, y_C1) + distance(x_A1, y_A1, x3, y3) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_B2, y_B2) + distance(x_C2, y_C2, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1) + distance(x_A1, y_A1, x3, y3), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1), distance(x_C1, y_C1, x_C2, y_C2) + distance(x_A1, y_A1, x3, y3) + distance(x_A2, y_A2, x_B2, y_B2) + distance(x_B1, y_B1, x1, y1), distance(x_C1, y_C1, x1, y1) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_B2, y_B2) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3) + distance(x_B1, y_B1, x_B2, y_B2), distance(x_C1, y_C1, x_B2, y_B2) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x1, y1) + distance(x_B1, y_B1, x_A1, y_A1), distance(x_C2, y_C2, x_C1, y_C1) + distance(x_B2, y_B2, x3, y3) + distance(x_A1, y_A1, x_B1, y_B1) + distance(x_A2, y_A2, x1, y1)))

		if surf_BCa == 0:

			BCa = v21 + v25 + v29 + v0 + v16

		elif surf_BCa == 1:

			BCa = v23 + v6 + v10 + v13

		elif surf_BCa == 2:

			BCa = v21 + v25 + v15 + v10 + v6

		elif surf_BCa == 3:

			BCa = v3 + v29 + v10 + v13

		elif surf_BCa == 4:

			BCa = v3 + v29 + v7 + v16

		elif surf_BCa == 5:

			BCa = v23 + v6 + v7 + v16

		elif surf_BCa == 6:

			BCa = v24 + v9 + v6 + v13

		elif surf_BCa == 7:

			BCa = v21 + v25 + v9 + v6 + v13

		elif surf_BCa == 8:

			BCa = v21 + v25 + v29 + v10 + v1

		else:

			BCa = v23 + v22 + v1 + v10

		

		#BCb-------------------------------

		surf_BCb = np.argmin((distance(x_C1, y_C1, x_C2, y_C2) + distance(x2, y2, x1, y1), distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x1, y1))) 		

		if surf_BCb == 0:

			BCb = v23 + v30

		else:

			BCb = v25 + v27

		#Bab-------------------------------

		surf_Bab = np.argmin((distance(x_A1, y_A1, x_A2, y_A2) + distance(x2, y2, x3, y3), distance(x_A1, y_A1, x3, y3) + distance(x_A2, y_A2, x2, y2))) 		

		if surf_Bab == 0:

			Bab = v0 + v32

		else:
	
			Bab = v6 + v11

		#Cab-------------------------------

		surf_Cab = np.argmin((distance(x_A2, y_A2, x_B1, y_B1) + distance(x_B2, y_B2, x2, y2) + distance(x_C1, y_C1, x_C2, y_C2) + distance(x_A1, y_A1, x3, y3), distance(x_B2, y_B2, x_B1, y_B1) + distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_A2, y_A2) + distance(x_B2, y_B2, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3) + distance(x_C1, y_C1, x2, y2), distance(x_B1, y_B1, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_C1, y_C1) + distance(x_B2, y_B2, x2, y2) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3), distance(x_B1, y_B1, x_B2, y_B2) + distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x_A2, y_A2) + distance(x_A1, y_A1, x3, y3), distance(x_B1, y_B1, x3, y3) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_A2, y_A2) + distance(x_B2, y_B2, x_C2, y_C2) + distance(x_C1, y_C1, x2, y2) + distance(x_A1, y_A1, x_A2, y_A2), distance(x_B1, y_B1, x_A2, y_A2) + distance(x_B2, y_B2, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A1, y_A1, x_C1, y_C1), distance(x_B2, y_B2, x_B1, y_B1) + distance(x_A2, y_A2, x2, y2) + distance(x_C1, y_C1, x_A1, y_A1) + distance(x_C2, y_C2, x3, y3)))

		if surf_Cab == 0:

			Cab = v10 + v16 + v21 + v23 + v6

		elif surf_Cab == 1:

			Cab = v13 + v25 + v29 + v0

		elif surf_Cab == 2:

			Cab = v10 + v16 + v2 + v29 + v25

		elif surf_Cab == 3:

			Cab = v14 + v21 + v29 + v0

		elif surf_Cab == 4:

			Cab = v14 + v21 + v9 + v6

		elif surf_Cab == 5:

			Cab = v13 + v25 + v9 + v6

		elif surf_Cab == 6:

			Cab = v18 + v19 + v25 + v0

		elif surf_Cab == 7:

			Cab = v10 + v16 + v19 + v25 + v0

		elif surf_Cab == 8:

			Cab = v10 + v16 + v21 + v29 + v3

		else:

			Cab = v13 + v11 + v3 + v29

		

		#ABCa------------------------------

		ABCa = v16 + v29 + v21 + v25

		#ABCb------------------------------

		ABCb = v10 + v25 + v6 + v29

		#ABab------------------------------

		ABab = v32

		#ACab------------------------------
	
		surf_ACab = np.argmin((distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_B2, y_B2, x_B1, y_B1), distance(x_B2, y_B2, x2, y2) + distance(x_B1, y_B1, x3, y3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_B2, y_B2, x2, y2) + distance(x_B1, y_B1, x_C1, y_C1) + distance(x_C2, y_C2, x3, y3), distance(x_B2, y_B2, x_C2, y_C2) + distance(x_B1, y_B1, x3, y3) + distance(x_C1, y_C1, x2, y2))) 		

		if surf_ACab == 0:

			ACab = v25 + v29 + v13

		elif surf_ACab == 1:

			ACab = v18 + v21 + v23

		elif surf_ACab == 2:

			ACab = v21 + v14 + v29

		else:

			ACab = v18 + v19 + v25

		#BCab------------------------------

		surf_BCab = np.argmin((distance(x_C1, y_C1, x2, y2) + distance(x_C2, y_C2, x3, y3) + distance(x_A2, y_A2, x_A1, y_A1), distance(x_A2, y_A2, x2, y2) + distance(x_A1, y_A1, x3, y3) + distance(x_C2, y_C2, x_C1, y_C1), distance(x_A2, y_A2, x2, y2) + distance(x_A1, y_A1, x_C1, y_C1) + distance(x_C2, y_C2, x3, y3), distance(x_A2, y_A2, x_C2, y_C2) + distance(x_A1, y_A1, x3, y3) + distance(x_C1, y_C1, x2, y2))) 		

		if surf_BCab == 0:

			BCab = v25 + v29 + v0

		elif surf_BCab == 1:

			BCab = v6 + v11 + v23

		elif surf_BCab == 2:

			BCab = v3 + v11 + v29

		else:

			BCab = v6 + v9 + v25

		#ABCab-----------------------------

		ABCab = v25 + v29

		
		
		BCabc = A

		ACabc = B

		ABabc = C

		ABCbc = a

		ABCac = b

		Cabc = AB

		Babc = AC

		BCbc = Aa

		BCac = Ab

		Aabc = BC

		ACbc = Ba

		ACac = Bb

		ABbc = Ca 

		ABac = Cb  

		ABCc = ab  

		abc = ABC  

		Cbc = ABa  

		Cac = ABb  

		Bbc = ACa  

		Bac = ACb  

		BCc = Aab  

		Abc = BCa  

		Aac = BCb  

		ACc = Bab  

		ABc = Cab  

		bc = ABCa 

		ac = ABCb  

		Cc = ABab  
	
		Bc = ACab  

		Ac = BCab  

		c = ABCab 

		

		alpha = Aa + Bb + Cc

		beta = A + B + C 

		gamma = AB + AC + BC

		delta = ABC

		epsilon = ABc + ACb + BCa

		eta = Ab + Ac + Ba + Bc + Ca + Cb

		theta = ABa + ABb + ACa + ACc + BCb + BCc

		omega = ab + ac + bc

		rho = a + b + c


		if CM == 1:

			surf_CM = 2*beta - gamma

		elif CM == 2:

			surf_CM = 2*beta + theta - eta

		elif CM == 3:

			surf_CM = gamma - 2*delta

		elif CM == 4:

			surf_CM = alpha

		elif CM == 5:

			surf_CM = theta - omega

		elif CM == 6:

			surf_CM = eta + theta - 2*omega - 2*rho

		elif CM == 7:

			surf_CM = eta + theta + beta - omega - 2*rho - epsilon

		elif CM == 8:

			surf_CM = theta - 2*rho

		elif CM == 9:

			surf_CM = theta + beta + gamma - epsilon - omega

		elif CM == 10:

			surf_CM = theta + epsilon - rho - 2*omega

		elif CM == 11:

			surf_CM = theta + 2*delta - 2*omega

		elif CM == 12:

			surf_CM = 2*alpha + eta - theta - 2*rho

		elif CM == 13:

			surf_CM = 2*alpha + gamma - theta

		elif CM == 14:

			surf_CM = alpha + 2*delta - omega

		elif CM == 15:

			surf_CM = 2*alpha + 2*beta + epsilon - eta - rho

		elif CM == 16:

			surf_CM = alpha + epsilon - omega - rho

		elif CM == 17:

			surf_CM = alpha + beta + gamma - epsilon

		elif CM == 18:

			surf_CM = alpha + omega - 2*rho

		else:

			surf_CM = alpha + beta + eta - epsilon - 2*rho

		#surf_CM = Cb

		return surf_CM

#draw bulk surfaces contributing to E_CM, the minimized f_CM

	surf_CM = surf_CM(E_CM.x)
	print(surf_CM)

	
	
	if surf_CM[0] != 0:
	
		geod(x_A1, y_A1, x_A2, y_A2, surf_CM[0])

	if surf_CM[1] != 0:
	
		geod(x_A1, y_A1, x_B1, y_B1, surf_CM[1])

	if surf_CM[2] != 0:
	
		geod(x_A1, y_A1, x_B2, y_B2, surf_CM[2])

	if surf_CM[3] != 0:
	
		geod(x_A1, y_A1, x_C1, y_C1, surf_CM[3])

	if surf_CM[4] != 0:
	
		geod(x_A1, y_A1, x_CM1, y_CM1, surf_CM[4])

	if surf_CM[5] != 0:
	
		geod(x_A1, y_A1, x_CM2, y_CM2, surf_CM[5])

	if surf_CM[6] != 0:
	
		geod(x_A1, y_A1, x_CM3, y_CM3, surf_CM[6])

	if surf_CM[7] != 0:
	
		geod(x_A2, y_A2, x_B2, y_B2, surf_CM[7])

	if surf_CM[8] != 0:
	
		geod(x_A2, y_A2, x_C1, y_C1, surf_CM[8])

	if surf_CM[9] != 0:
	
		geod(x_A2, y_A2, x_C2, y_C2, surf_CM[9])

	if surf_CM[10] != 0:
	
		geod(x_A2, y_A2, x_CM1, y_CM1, surf_CM[10])

	if surf_CM[11] != 0:
	
		geod(x_A2, y_A2, x_CM2, y_CM2, surf_CM[11])

	if surf_CM[12] != 0:
	
		geod(x_A2, y_A2, x_CM3, y_CM3, surf_CM[12])

	if surf_CM[13] != 0:
	
		geod(x_B1, y_B1, x_B2, y_B2, surf_CM[13])

	if surf_CM[14] != 0:
	
		geod(x_B1, y_B1, x_C1, y_C1, surf_CM[14])

	if surf_CM[15] != 0:
	
		geod(x_B1, y_B1, x_C2, y_C2, surf_CM[15])

	if surf_CM[16] != 0:
	
		geod(x_B1, y_B1, x_CM1, y_CM1, surf_CM[16])

	if surf_CM[17] != 0:
	
		geod(x_B1, y_B1, x_CM2, y_CM2, surf_CM[17])

	if surf_CM[18] != 0:
	
		geod(x_B1, y_B1, x_CM3, y_CM3, surf_CM[18])

	if surf_CM[19] != 0:
	
		geod(x_B2, y_B2, x_C2, y_C2, surf_CM[19])

	if surf_CM[20] != 0:
	
		geod(x_B2, y_B2, x_CM1, y_CM1, surf_CM[20])

	if surf_CM[21] != 0:
	
		geod(x_B2, y_B2, x_CM2, y_CM2, surf_CM[21])

	if surf_CM[22] != 0:
	
		geod(x_B2, y_B2, x_CM3, y_CM3, surf_CM[22])

	if surf_CM[23] != 0:
	
		geod(x_C1, y_C1, x_C2, y_C2, surf_CM[23])

	if surf_CM[24] != 0:
	
		geod(x_C1, y_C1, x_CM1, y_CM1, surf_CM[24])

	if surf_CM[25] != 0:
	
		geod(x_C1, y_C1, x_CM2, y_CM2, surf_CM[25])

	if surf_CM[26] != 0:
	
		geod(x_C1, y_C1, x_CM3, y_CM3, surf_CM[26])

	if surf_CM[27] != 0:
	
		geod(x_C2, y_C2, x_CM1, y_CM1, surf_CM[27])

	if surf_CM[28] != 0:
	
		geod(x_C2, y_C2, x_CM2, y_CM2, surf_CM[28])

	if surf_CM[29] != 0:
	
		geod(x_C2, y_C2, x_CM3, y_CM3, surf_CM[29])

	if surf_CM[30] != 0:
	
		geod(x_CM1, y_CM1, x_CM2, y_CM2, surf_CM[30])

	if surf_CM[31] != 0:
	
		geod(x_CM1, y_CM1, x_CM3, y_CM3, surf_CM[31])

	if surf_CM[32] != 0:
	
		geod(x_CM2, y_CM2, x_CM3, y_CM3, surf_CM[32])
			
	
	
	plt.show()

else:
	
	print("Disconnected entanglement wedge") 

#fig.savefig('AdS_geod.png')

plt.show()



 
