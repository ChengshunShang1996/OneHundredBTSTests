import os
import numpy as np
import matplotlib.pyplot as plt

#levels = np.linspace(43000,45000,10)
levels = 50

x,y,z = [],[],[]
BTS_peak_points_path_and_name = os.path.join(os.getcwd(),'BTS_peak_points.dat')
with open(BTS_peak_points_path_and_name, "r") as f_w_peak_points:
    for line in f_w_peak_points:
        values = [float(s) for s in line.split()]
        x.append(values[0])
        y.append(values[1])
        z.append(values[2])

#plt.tricontour(x, y, z, levels=[46193, 46195], colors = 'blue', linewidths = 0.5)

'''
x,y,z = [],[],[]
BTS_peak_points_path_and_name = os.path.join(os.getcwd(),'BTS_peak_points_01.dat')
with open(BTS_peak_points_path_and_name, "r") as f_w_peak_points:
    for line in f_w_peak_points:
        values = [float(s) for s in line.split()]
        x.append(values[0])
        y.append(values[1])
        z.append(values[2])
'''
#plt.contourf(X, Y, Z, 20, cmap=plt.get_cmap('YlGn'))
#cs = plt.tricontour(x, y, z, levels=levels, colors = 'white', linewidths = 0.1)
plt.tricontour(x, y, z, levels=[46193, 46195], colors = 'white', linewidths = 0.5)
plt.tricontourf(x, y, z, levels=levels, cmap='coolwarm')
#cs.clabel(inline=True, fmt='%d', fontsize = 'smaller', manual=true)

plt.xlabel('Tension_limit')
plt.ylabel('Shear_limit')
plt.title('BTS strength')

plt.colorbar()
plt.show()