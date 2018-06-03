import reweight as rr


###################################################################
# EXAMPLE 1: use NOE for reweighting.

# This file contains m+1 rows. The first line is a header describing the type of data
# and the error prior. 
# # DATA=NOE PRIOR=GAUSS 
# The other m rows correspond to each NOE experimental average in the following format:
# COL1=LABEL COL2=R_AVERAGE COL3=SIGMA 
exp0 = "data/noe_CCCC_exp.dat"

# The following file contains m noe distances calculated for all frames in a simulation
# The file has m+1 column (first is a label, i.e. frame number) and n rows, where n is
# the number of samples. In general, it makes sense to have n > m.
# It is VERY important that the experimental data in row i corresponds to column j
sim0 = "data/NOE.noe.calc.dat"

# Here initalize class, reads experimental and simulated data. 
rew = rr.Reweight([exp0],[sim0])

# do optimization. Theta is a tunable parameter that can be set at optimization
rew.optimize(theta=1.,method="MAXENT")

# Back calculate data using new and old weights.
# this command will write two files: 
# i. example1.stats.dat, with averages, chi_squared, and rmsd before and after minimization
# ii. examples1.weights.dat, with weights before and after minimization
rew.weight_exp([exp0],[sim0],'example1_noe_rescale')

# It is useful to calculate the agreement with experimental averages not used in the reweighting.
# In this example, we used  3J scalar couplings to cross-validate
# The experimental data are in the format
# DATA=JCOUPLINGS PRIOR=GAUSS
# COL1=LABEL COL2=AVG COL3=SIGMA
exp1 = "data/j3_CCCC.exp.dat"
sim1 = "data/J3.calc.dat"

# do weighting. If plot=True, full distributions are plotted to example1_3j.pdf
rew.weight_exp([exp1],[sim1],'example1_3j_rescale',plot=True)

# Finally, one can calculate any other quantity. Here, the RMSD distance from some structure is calculated
# before/after reweighting.
distance="data/distance.dat"
rew.weight([distance],'distance1',plot=True)


###################################################################
# EXAMPLE 2: use NOE+J3 for reweighting.
# initalize class, reads experimental and simulated data. 
rew2 = rr.Reweight([exp0,exp1],[sim0,sim1])
rew2.optimize(theta=1.0,method="MAXENT")


##################################################################
# EXAMPLE 3: use NOE, 3J and 
# in this case, NOE, 3J and NOE boundaries are used in the minimization
exp2 = "data/noe_CCCC_exp_bound.dat"

# note that the experimental datafile exp2 is in the format
# DATA=NOE PRIOR=GAUSS (first line only)
# COL1=LABEL  COL2=R COL3=SIGMA COL4=UPPER (for upper boundary) or LOWER (for lower boundary)
sim2 = "data/NOE_BOUND.noe.calc.dat"

rew2.weight_exp([exp0,exp1,exp2],[sim0,sim1,sim2],'example2_rescale')

###############################################################
