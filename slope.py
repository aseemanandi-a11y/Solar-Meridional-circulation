#use xarray
#find m and yintercept
#create list of values for line
#plot separately then together
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Load the NetCDF file
file_path = '/home/aseema/Projects/CSIR2025/IO.nc'
ds = xr.open_dataset(file_path)

# Extract coordinates and data
lon = ds['lon'].values
lat = ds['lat'].values
tos = ds['tos'].values  # Shape: (time, lat, lon)

# loop to create list of all years
av=list()
x=list()
for yr in range(1964,2015):
   x.append(yr)
   
   m=list()
   t=list()

   ye=(yr-1964)*12+1368
   for i in range (0,12):
      m.append(ye+i)
      
   for k in m:
      for i in range(10,30):
         for j in range(20,41):
            c=tos[k][i][j]
            t.append(c)
            valid_mask = ~np.isnan(t) & (t != 1.e+20)
            tval = np.array(t)[valid_mask]
   a=sum(tval)/len(tval)
   av.append(a)
coefficients = np.polyfit(x, av, 1)
y_pred=list()
for i in x:
   y_pred.append(coefficients[0] * i + coefficients[1])

plt.scatter(x,av,color='red')   
plt.plot(x,y_pred)
plt.xlabel('years')
plt.ylabel('temperature in celsius')
plt.title('average temperature over BoB from 1964 to 2014')
plt.figtext(0.8, 0.98, "CanESM5 CMIP6 historical simulation", ha="center", fontsize=7)
plt.figtext(0.8, 0.96, "resolution = 1 degree", ha="center", fontsize=7)
plt.show()


