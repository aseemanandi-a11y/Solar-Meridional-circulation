#current update:  capped intensity
# The plan is to see at what cap it works with dr then just directly apply to mc

from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import glob 
import csv
import numpy as np 
from datetime import datetime


#functions
def load_map(file_path):

  with fits.open(file_path) as hdul:
      image = hdul[0].data.copy().astype(float)
      header = hdul[0].header.copy()
      valid = np.isfinite(image)

      within_cap = valid & (np.abs(image) <= 500)

      percentage = 100 * np.sum(within_cap) / np.sum(valid)

      print(f"{percentage:.2f}% of valid pixels are within ±800 G")
      print(f"{100-percentage:.2f}% of valid pixels are being capped")
      image = np.clip(image, -Bmax, Bmax)
  return image, header

def lattorow(latitude):
  return np.argmin(np.abs(latitudes - latitude)) 

def patch(img, alat):
  rows = np.where((latitudes >= alat-5) & (latitudes <= alat+5))[0]
  patch = img[rows,0:360]
  return patch,rows
  
def vget(best_shift, latitude):
    row = lattorow(latitude)
    row2 = row + best_shift
    delta_lat = latitudes[row2] - latitudes[row]
    distance = 6.957e8 * np.radians(delta_lat)
    velocity = distance / (rot_t * 86400)
    return velocity


#constants  
lat_theory = np.linspace(-50, 50, 1000)
nlat = 180
sinlat = np.linspace(-1 + 1/nlat, 1 - 1/nlat, nlat)
latitudes = np.degrees(np.arcsin(sinlat))
omega_carr = 360/27.2753
analysis_lats = np.arange(-50,51,5)

Bmax = 500


approved_dates = [ '1999-12-19T00:00:00', '2000-01-21T16:23:23', '2000-02-21T16:01:10', '2000-02-28T21:03:57', '2000-03-29T15:53:25']

finalset = []
finalcorr = []
   
#=========================== DIFFERENTIAL ROTATION STARTS HERE


#loop starts here separates every pair
for date in approved_dates:
  with open("observations.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
      if row["DATE"] == date:
        row2 = next(reader)
        image1, header1 = load_map(row["filepath"])
        image2, header2 = load_map(row2["filepath"])
        
        #fixing the time difference issue
        date2 = row2["DATE"]
        t1 = datetime.fromisoformat(date)
        t2 = datetime.fromisoformat(date2)
        print(date,date2)
        rot_t = (t2 - t1).total_seconds()/86400
               
  #lists to store the values of all the latitudes in a pair
  latshifts = []
  latcorr = []
  wmeasured = []
  wtheory = []
  
  #starting for loop to go through all latitudes in analysis_lats for a given pair
  #analysis_lats = [50]
  for a in analysis_lats:
    
    #getting the differential rotation shifts
    patch1,rows1 = patch(image1,a)
    nrows = patch1.shape[0]
    centre_row = lattorow(a)

    shifts = []
    correlations = []
    patch2,rows2 = patch(image2,a)
    for dx in range(-50, 51):

      if dx > 0:
        p1 = patch1[:, dx:]
        p2 = patch2[:, :-dx]

      elif dx < 0:
        p1 = patch1[:, :dx]
        p2 = patch2[:, -dx:]

      else:
        p1 = patch1
        p2 = patch2

      corr = np.corrcoef(p1.ravel(),p2.ravel())[0,1]
      correlations.append(corr)
      shifts.append(dx)

    #separating best shift values for pair
    best = np.argmax(correlations)
    best_shift = shifts[best] 
    delta_omega = best_shift / rot_t
    omega_measured = omega_carr - delta_omega
    best_corr = correlations[best]
    latshifts.append(best_shift)
    latcorr.append(best_corr)
    wmeasured.append(omega_measured)
  
  #compiling results of all pairs
  finalset.append(np.array(wmeasured))
  finalcorr.append(np.array(best_corr))
  
#finding the mean diff.rot values  
finalcorr = np.array(finalcorr)
finalset = np.array(finalset)
avgcorr = np.mean(finalcorr, axis = 0)
avgrot = np.mean(finalset, axis=0)
stddevrot = np.std(finalset, axis=0)

#improvement checker
mean_std = np.mean(stddevrot)
print('mean_std = ',mean_std)
print('curve',np.mean(np.abs(np.diff(avgrot, n=2))))

#the csv 
with open("avg_dr23.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Latitude", "Omega"])

    for lat, omega in zip(analysis_lats, avgrot):
        writer.writerow([lat, omega])

  
#getting theoretical values
for lat in lat_theory:
    omega_theory = (14.38 - 2.72*np.sin(np.radians(lat))**2 - 1.42*np.sin(np.radians(lat))**4 - 0.986)
    wtheory.append(omega_theory)
  

plt.figure(figsize=(8,5))

# Your measurements with error bars
plt.errorbar(
    analysis_lats,
    avgrot,
    yerr=stddevrot,
    fmt='o',
    capsize=4,
    label='Measured'
)

# Theoretical curve
plt.plot(
    lat_theory,
    wtheory,
    'r-',
    linewidth=2,
    label='Theoretical (Carrington frame)'
)

plt.xlabel("Latitude (degrees)")
plt.ylabel("Rotation rate (deg/day)")
plt.title("Solar Differential Rotation")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig("diffrot_500_rising.png", dpi=300)



#=========================== DIFFERENTIAL ROTATION ENDS HERE

#mean_std0 =  0.5992403109544265
#mean_1k =  0.5911359105405174
#mean_800 =  0.5887530937952548

#curve0 = 0.07752150749513008
#curve1k = 0.07383527039214444
#curve800 = 0.07044092129532997


#rising
#mean_std =  0.6515003687394715
#curve = 0.14013850715979648
#mean_std_800 =  0.6460445574976111
#curve_800 = 0.1286626702538038
#mean_std_500 =  0.6456377154104513
#curve_500 = 0.1306646972854955


#falling
#mean_std =  0.5529044900867282
#curve = 0.08829037389375216
#mean_std 800 =  0.5399221732063724
#curve 800 = 0.07664885214475123
#mean_std_500 =  0.5366347461854849
#curve_500 = 0.08178119772765993



#=========================== MERIDIONAL CIRCULATTION STARTS HERE


finalset = []
finalcorr = []

#pull up list of differential rotation shifts
dr_shifts = []
with open("avg_dr23.csv") as f:
   reader = csv.DictReader(f)
   i=0
   for row in reader:
      i=i+1
      dr_shifts.append(row["Omega"])    


#start loop for pairs
for dat in range(0,len(approved_dates)):
  date = approved_dates[dat]
  with open("observations.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
      if row["DATE"] == date:
        row2 = next(reader)
        image1, header1 = load_map(row["filepath"])
        image2, header2 = load_map(row2["filepath"])
        
        #fixing the time difference issue
        date2 = row2["DATE"]
        t1 = datetime.fromisoformat(date)
        t2 = datetime.fromisoformat(date2)
        print(date,date2)
        rot_t = (t2 - t1).total_seconds()/86400
        break

  latshifts = []
  latcorr = []
  latv = []
  
  #take slices and start loop on analysis lats
  for i in range(0,len(analysis_lats)):
    a  = analysis_lats[i]
    patch1,row1 = patch(image1,a)
    
    #starting correlation
    shifts = []
    correlations = []
    center_row = lattorow(a)
    nrows = patch1.shape[0]  
    
    for dy in range(-10, 11):
      shifted_center = center_row + dy
      
      if nrows%2 == 0:
        start = shifted_center - nrows//2
        end = shifted_center + nrows//2 
      else:
        start = shifted_center - nrows//2
        end = shifted_center + nrows//2 + 1
      #print(start,end)
      
      if start < 0 or end > image2.shape[0]:
        continue

      patch2 = image2[start:end, :]
      
      if patch2.shape != patch1.shape:
        print(patch1.shape,patch2.shape,dy)
      
      #shifting patch  2 
      dr = float(dr_shifts[i])
      shift_pixels = int(round((omega_carr - dr) * rot_t))
      patch2_rot = np.roll(patch2, shift=shift_pixels, axis=1)
      
      #correlating
      corr = np.corrcoef(patch1.ravel(), patch2_rot.ravel())[0, 1]

      shifts.append(dy)
      correlations.append(corr)
    
    #append to latshifts and stuff
    correlations = np.array(correlations)
    best = np.argmax(correlations)     
    best_shift = shifts[best]
    velocity = vget(best_shift,a)
    latshifts.append(best_shift)
    latcorr.append(correlations[best])
    latv.append(velocity)
    
  #append latshifts to the final set 
  finalset.append(latv)  
  finalcorr.append(latcorr)


#get avg
finalset = np.array(finalset)
avg_v = np.mean(finalset, axis=0)
finalcorr = np.array(finalcorr)
avg_corr = np.nanmean(finalcorr, axis=0)
stddev_v = np.std(finalset, axis=0)

#improvement checker
mean_std = np.mean(stddev_v)
print('mean_std = ',mean_std)
print('curve',np.mean(np.abs(np.diff(avg_v, n=2))))

#for j in range(0,len(analysis_lats)):
#  print(analysis_lats[j],avg_v[j],avg_corr[j])

#plotting
plt.figure(figsize=(6,4))

plt.errorbar(
    analysis_lats,
    avg_v,
    yerr=stddev_v,
    fmt='o',
    capsize=4,
    label='Measured'
)

plt.plot(analysis_lats, avg_v, 'o-')

plt.xlabel("Latitude (degrees)")
plt.ylabel("Meridional velocity (m/s)")
plt.title("Meridional flow")
plt.grid(True)

plt.savefig("merflow_500_rising.png")
plt.close()

#mean_std800 =  19.316966511404935
#mean_std =  19.25192526499945
#curve800 3.9892210303299325
#curve = 3.902802122115327

#rising
#mean_std =  23.140622181233717
#curve = 7.381553330141507
#mean_std_800 =  23.146211178473035
#curve_800 = 7.166525253954736
#mean_std_500 =  22.92129645664564
#curve_500 = 6.866638173916186

#falling
#mean_std =  16.250432544174956
#curve = 6.1094478034081465
#mean_std_800 =  16.31315614573021
#curve_800 = 6.0955921779946625
#mean_std_500 =  16.206978661992636
#curve_500 = 5.5764611782980555

