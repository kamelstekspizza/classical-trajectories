import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import sys

dirs = sys.argv[1:]

colormap = plt.cm.jet #or any other colormap
normalize = matplotlib.colors.Normalize(vmin=-0.5, vmax=0)


for dir in dirs:
    DATA = np.loadtxt(f'{dir}/output.dat')
    fig,ax = plt.subplots()
    sc = ax.scatter(DATA[:,0],DATA[:,1],c = DATA[:,2],cmap = colormap,norm= normalize,s = 1) #Make marker color depend on energy!
    cb = fig.colorbar(sc)


    
plt.show()
