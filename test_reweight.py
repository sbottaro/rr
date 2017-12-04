import reweight as rr


###################################################################
# EXAMPLE 1: use NOE for reweighting.

# This file contains NOE m experimental averages in the following format:
# DATA=NOE PRIOR=GAUSS (first line only)
# COL1=PROTON1 COL2=PROTON2 COL3=R_MIN COL4=R_AVG COL5=RMAX
# Also accepts data in the format (to be tested)
# COL1=PROTON1 COL2=PROTON2 COL3=R_AVG COL4=SIGMA_R
exp0 = "data/noe_CCCC_exp.dat"

# This file contains the m noe distances calculated for all frames in a simulation
# This files has m+1 column (first is a label, i.e. frame number) and n rows, where n is
# the number of samples. In general, it makes sense to have n > m.
# It is VERY important that the experimental data in row i corresponds to column j
sim0 = "data/NOE.noe.calc.dat"

# initalize class, reads experimental and simulated data. 
rew = rr.Reweight([exp0],[sim0])

# do optimization. You can set theta here. 
# For testing purposes, and if n<m, you may want to try method="BER". Very slow.
#rew.optimize(theta=1.,method="MAXENT")

# Back calculate data using new/old weights.
# this command will write two files: 
# i. example1.stats.dat, with averages, chi_squared, and rmsd before and after minimization
# ii. examples1.weights.dat, with weights before/after minimization
#rew.weight_exp([exp0],[sim0],'example1_noe')

# one can also calculate other averages before/after reweighting, e.g. 3J scalar couplings
# experimental data with format
# DATA=JCOUPLINGS PRIOR=GAUSS
# COL1=LABEL COL2=AVG COL3=SIGMA
exp1 = "data/j3_CCCC.exp.dat"
sim1 = "data/J3.calc.dat"

# do weighting. If plot=True, full distributions are plotted to example1_3j.pdf
#rew.weight_exp([exp1],[sim1],'example1_3j',plot=True)

# Finally, one can calculate any other quantity. Here, the RMSD distance from some structure is calculated
# before/after reweighting.
distance="data/distance.dat"
#rew.weight([distance],'distance1',plot=True)


###################################################################
# EXAMPLE 2: use NOE+J3 for reweighting.
# initalize class, reads experimental and simulated data. 
rew2 = rr.Reweight([exp0,exp1],[sim0,sim1])
rew2.optimize(theta=1.0,method="MAXENT")


# Back calculate data using new/old weights.
# this command will write two files: 
# i. example1.stats.dat, with averages, chi_squared, and rmsd before and after minimization
# ii. examples1.weights.dat, with weights before/after minimization
# in this case, NOE, 3J and UPPER/LOWER NOE are calculated
exp2 = "data/noe_CCCC_exp_bound.dat"
# note that the experimental datafile is in the format
# DATA=NOE PRIOR=GAUSS (first line only)
# COL1=PROTON1 COL2=PROTON2 COL3=R COL4=SIGMA COL5=UPPER (for upper boundary) or LOWER (for lower boundary)
sim2 = "data/NOE_BOUND.noe.calc.dat"

rew2.weight_exp([exp0,exp1,exp2],[sim0,sim1,sim2],'example2')

###############################################################
