import os
import numpy as np
import matplotlib.pyplot as plt

x,y,z = [],[],[]
BTS_peak_points_path_and_name = os.path.join(os.getcwd(),'BTS_peak_points.dat')
with open(BTS_peak_points_path_and_name, "r") as f_w_peak_points:
    for line in f_w_peak_points:
        values = [float(s) for s in line.split()]
        x.append(values[0])
        y.append(values[1])
        z.append(values[2])

X, Y = np.meshgrid(x, y)
#Z = z.reshape(21, 21)
Z = z

plt.contourf(X, Y, Z, 20, cmap=plt.get_cmap('YlGn'))
#plt.contourf(X, Y, Z, 20)
#plt.pcolor(xi, yi, zi, cmap=plt.get_cmap('hot'))
plt.show()