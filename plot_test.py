import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

x,y,z = [],[],[]
BTS_peak_points_path_and_name = os.path.join(os.getcwd(),'BTS_peak_points.dat')
with open(BTS_peak_points_path_and_name, "r") as f_w_peak_points:
    for line in f_w_peak_points:
        values = [float(s) for s in line.split()]
        x.append(values[0])
        y.append(values[1])
        z.append(values[2])

#levels = np.linspace(43000,45000,10)
levels = 50

#plt.contourf(X, Y, Z, 20, cmap=plt.get_cmap('YlGn'))
#cs = plt.tricontour(x, y, z, levels=levels, colors = 'white', linewidths = 0.1)
#plt.tricontour(x, y, z, levels=[43999, 44001], colors = 'white', linewidths = 0.5)
#plt.tricontourf(x, y, z, levels=levels, cmap='coolwarm')
plotx,ploty, = np.meshgrid(np.linspace(np.min(x),np.max(x),10),\
                           np.linspace(np.min(y),np.max(y),10))
plotz = interp.griddata((x,y),z,(plotx,ploty),method='linear')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(plotx,ploty,plotz,cstride=1,rstride=1,cmap='coolwarm')

plt.xlabel('sigma_limit')
plt.ylabel('tension_limit')
plt.title('BTS strength')

#plt.colorbar()
plt.show()