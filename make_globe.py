from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import sys

lat_0 = float(sys.argv[1])
lon_0 = float(sys.argv[2])

plt.figure(figsize=(4,4))
map = Basemap(projection='ortho',lat_0=lat_0,lon_0=lon_0,resolution='l')
map.drawcoastlines(linewidth=0.25)

map.fillcontinents(color='coral')
map.drawmapboundary(fill_color='aqua')
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))

plt.show()
