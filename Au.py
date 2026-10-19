#current update:  Gaussian smoothing 
# The plan is to see if it is at all useful with our current dataset then to check what value of sigma gives the best results and if that works we will try to include the 'bad' maps that were previously excluded

from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import glob 
import csv
import numpy as np 
from datetime import datetime
from scipy.ndimage import gaussian_filter


#functions
def load_map(file_path):

  with fits.open(file_path) as hdul:
      image = hdul[0].data.copy().astype(float)
      header = hdul[0].header.copy()

  # Identify valid pixels
  valid = np.isfinite(image)

  # Replace NaNs with zero temporarily
  image_zero = np.where(valid, image, 0.0)

  # Gaussian smoothing of the data
  smoothed = gaussian_filter(image_zero, sigma=3)

  # Gaussian smoothing of the validity mask
  weights = gaussian_filter(valid.astype(float), sigma=3)
  smoothed = smoothed / weights
  smoothed[weights == 0] = np.nan

  return smoothed, header

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


approved_dates = ['1996-08-09T00:00:00', '1996-09-01T00:00:00', '1996-09-23T00:00:00', '1996-10-23T00:00:00', '1996-11-24T00:00:00', '1996-12-13T00:00:00', '1997-01-13T00:00:00', '1997-02-05T00:00:00', '1997-03-12T00:00:00', '1997-04-06T00:00:00', '1997-04-30T00:00:00', '1997-05-29T00:00:00', '1997-06-25T00:00:00', '1997-07-17T00:00:00', '1997-08-18T00:00:00', '1997-09-18T00:00:00', '1997-10-09T00:00:00', '1997-11-07T00:00:00', '1997-12-01T00:00:00', '1998-01-05T00:00:00', '1998-01-28T00:00:00', '1998-02-27T00:00:00', '1998-03-25T00:00:00', '1998-05-05T00:00:00', '1998-05-19T00:00:00', '1998-07-28T00:00:00', '1998-07-18T00:00:00', '1998-08-13T00:00:00', '1998-09-02T00:00:00', '1998-10-01T00:00:00', '1998-10-28T00:00:00', '1998-11-30T00:00:00', '1998-12-24T00:00:00', '1999-02-19T00:00:00', '1999-02-19T00:00:00', '1999-03-10T00:00:00', '1999-04-08T00:00:00', '1999-05-04T00:00:00', '1999-06-08T00:00:00', '1999-06-30T00:00:00', '1999-07-23T00:00:00', '1999-08-30T00:00:00', '1999-09-22T00:00:00', '1999-11-04T00:00:00', '1999-11-17T00:00:00', '1999-12-19T00:00:00', '2000-01-21T16:23:23', '2000-02-21T16:01:10', '2000-02-28T21:03:57', '2000-03-29T15:53:25', '2000-04-26T21:53:00', '2000-05-18T23:26:42', '2000-06-15T20:10:30', '2000-07-18T21:25:24', '2000-08-18T20:34:50', '2000-09-08T19:52:47', '2000-10-05T16:54:13', '2000-11-17T18:15:40', '2000-11-30T16:54:29', '2000-12-28T18:39:08', '2001-01-22T23:31:32', '2001-02-14T20:13:54', '2001-03-19T18:23:59', '2001-04-26T00:21:18', '2001-05-14T16:17:28', '2001-06-05T23:12:15', '2001-07-09T23:47:08', '2001-08-01T18:21:59', '2001-08-29T18:59:37', '2001-09-26T17:35:17', '2001-11-02T21:25:38', '2001-11-19T19:59:40', '2001-12-12T23:58:56', '2002-01-08T21:26:04', '2002-02-14T20:13:00', '2002-03-08T20:06:32', '2002-04-12T19:41:58', '2002-05-14T21:43:40', '2002-06-04T16:36:02', '2002-06-26T20:09:16', '2002-07-19T22:30:13', '2002-08-22T00:41:35', '2002-09-18T21:11:17', '2002-10-30T11:07:51', '2002-11-08T16:50:51', '2002-12-06T20:55:17', '2003-01-07T18:42:26', '2003-01-27T23:57:00', '2003-02-25T00:06:51', '2003-03-21T16:27:41', '2003-04-17T22:39:01', '2003-05-16T23:39:50', '2003-06-13T00:38:12', '2003-07-10T23:04:11','2003-10-04T16:57:00', '2003-10-29T17:36:00', '2003-11-27T17:02:00', '2003-12-21T19:21:31', '2004-01-19T18:18:21', '2004-02-17T17:44:43', '2004-03-16T16:52:53', '2004-07-02T18:52:12', '2004-07-29T17:23:41', '2004-08-26T17:10:29', '2004-09-22T18:14:33', '2004-10-20T16:33:30', '2004-11-14T16:22:18', '2004-12-13T23:15:33', '2005-04-02T15:31:17', '2005-04-29T16:37:59', '2005-05-26T15:41:57', '2005-06-22T14:25:45', '2005-08-16T23:15:14', '2005-12-30T18:07:28', '2006-05-16T17:29:50', '2006-06-12T21:15:45', '2006-07-09T14:44:21', '2006-09-01T16:57:18', '2006-09-29T19:41:18', '2006-11-22T20:15:03', '2007-04-08T20:11:02', '2007-05-05T17:58:57', '2007-06-01T16:11:42', '2007-09-18T17:35:51', '2007-10-15T18:30:50', '2008-01-03T17:23:44', '2008-02-02T17:48:07', '2008-02-29T15:20:33', '2008-04-24T20:12:36', '2008-05-20T15:11:27', '2008-06-17T18:35:32', '2008-09-06T21:24:45', '2008-10-03T18:47:56', '2008-10-31T16:05:23', '2008-11-24T17:08:17']

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
plt.savefig("diffrot_Gauss3.png", dpi=300)


#=========================== DIFFERENTIAL ROTATION ENDS HERE

#mean_std0 =  0.5992403109544265
#mean_stdg1 =  0.5490982642845362
#mean_stdg2 =  0.5007768006555436
#mean_stdg3 =  0.4671079326311413
#mean_stdg4 =  0.4108248855756707 no plot

#curve0 = 0.07752150749513008
#curveg1 = 0.07739430764815851
#curveg2 = 0.0721453269106727
#curveg3 = 0.04162389770858227
#curveg4 = 0.03649629785270429 no plot

#rising
#mean_std_g3 =  0.5473080973830576
#curve_g3 = 0.06358186627009983
#mean_std =  0.589143150115131
#curve 0.09072543040863261
#mean_std_G1 =  0.6156558514124875
#curve_G1 0.1140606612703457
#mean_std =  0.6515003687394715
#curve = 0.14013850715979648

#falling
#mean_std =  0.5529044900867282
#curve = 0.08829037389375216
#mean_std_g1 =  0.48637704998935904
#curve_g1 = 0.10007618260070404
#mean_std_g2 =  0.4240186770332807
#curve_g2 = 0.07377794955877065
#mean_std_g3 =  0.3962808957713054
#curve_g3 = 0.05114351369975969



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

plt.savefig("merflow_Gauss3.png")
plt.close()

#=========================== MERIDIONAL CIRCULATION ENDS HERE
#for the whole cycly the std dev reduced as sigma increased but curve increased so we stop at sigma = 1
#mean_std =  19.25192526499945
#mean_stdg1 =  18.99170241575076
#curve = 3.902802122115327
#curveg1 = 4.614539184209977

#rising the best for curveseems to be at sigma = 2 here 
#mean_std_g3 =  21.316590011166028
#curve_g3 = 6.702480177394812
#mean_std_g2 =  22.50597724780133
#curve_g2 = 5.118322160480154
#mean_std_g1_rising =  22.586324701139798
#curve_g1_rising = 6.218042946412938
#mean_std =  23.140622181233717
#curve = 7.381553330141507

#falling
#mean_std =  16.250432544174956
#curve = 6.1094478034081465
#mean_std_g1 =  15.230096663541918
#curve_g1 = 6.210407158483825
#mean_std_g2 =  14.57504682910296
#curve_g2 = 7.142173398753805
#mean_std_g3 =  13.698491853685804
#curve_g3 = 6.721336066560546
