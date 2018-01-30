import reweight as rr
import numpy as np
import matplotlib.pyplot as plt

exp0 = "data/saxs_exp.dat"
sim0 = "data/saxs_calc.dat"

# initalize class, reads experimental and simulated data. 
rew = rr.Reweight([exp0],[sim0])

# optimize using different theta
thetas = range(15,150,5)
data = []

for t in thetas:
    rew.optimize(theta=t,method="MAXENT")
    # i. saxs_theta.stats.dat, with averages, chi_squared, and rmsd before and after minimization
    # ii. saxs_theta.weights.dat, with weights before/after minimization
    rew.weight_exp([exp0],[sim0],'saxs_%04.0f' % t)

    # read number of effective frames, chi2, srel from stats file
    fh = open("saxs_%04.0f.stats.dat" % t)
    for line in fh:
        if("neff" in line): neff = float(line.split()[2])
        if("Srel" in line): srel = float(line.split()[3])
        if("chi2" in line):
            chi2_b = float(line.split()[2])
            chi2_a = float(line.split()[3])
    data.append([t,chi2_a,chi2_b,srel,neff])
data = np.array(data)

# now do some plotting
f, axarr = plt.subplots(2,2,figsize=(7,6))

# plot chi squared as a function of theta
axarr[0,0].plot(data[:,0],data[:,1],c='k',label='after')
axarr[0,0].plot(data[:,0],data[:,2],c='k',ls='--',label="before")
axarr[0,0].set_xlabel("Theta")
axarr[0,0].set_ylabel("chi_sq")
axarr[0,0].legend()

# plot relative entropy as a function of theta
axarr[0,1].plot(data[:,0],data[:,3],c='k')
axarr[0,1].set_xlabel("Theta")
axarr[0,1].set_ylabel("relative entropy")

# plot fraction of effctive frames as a function of theta
axarr[1,0].plot(data[:,0],data[:,4],c='k')
axarr[1,0].set_xlabel("Theta")
axarr[1,0].set_ylabel("n eff")

# plot srel vs chi_squared
axarr[1,1].plot(data[:,3],data[:,1],c='k')
axarr[1,1].set_xlabel("relative entropy")
axarr[1,1].set_ylabel("chi_sq")
# label every 4 theta
for j in range(0,data.shape[0],4):
    axarr[1,1].text(data[j,3],data[j,1],"%.0f" % data[j,0])

plt.savefig("saxs.png",dpi=300)
