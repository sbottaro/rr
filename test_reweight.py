import reweight_04 as rr
import sys

exp0 = "../CCCC/j3_CCCC.exp.dat"
exp1 = "../CCCC/noe_CCCC_exp.dat"
exp2 = "../CCCC/noe_CCCC_exp_bound.dat"

sim0 = "../CCCC_TIP3P_300/J3.calc.dat"
sim1 = "../CCCC_TIP3P_300/NOE.noe.calc.dat"
sim2="../CCCC_TIP3P_300/NOE_BOUND.noe.calc.dat"


# initalize class
rew = rr.Reweight([exp0],[sim0])
# do optimization
rew.optimize(theta=1.,method="MAXENT")
# back calculate averages (with experimental data)
rew.weight_exp([exp0,exp1,exp2],[sim0,sim1,sim2],'example1')


# initalize class, use exp0 and exp1 together for reweighting
rew = rr.Reweight([exp0,exp1],[sim0,sim1])

# do optimization
rew.optimize(theta=1.,method="MAXENT")

# back calculate averages (with experimental data)
rew.weight_exp([exp0,exp1,exp2],[sim0,sim1,sim2],'example2')


# initalize class, use exp0 exp1 exp2 together for reweighting
rew = rr.Reweight([exp0,exp1,exp2],[sim0,sim1,sim2])

# do optimization
rew.optimize(theta=1.,method="MAXENT")

# back calculate averages (with experimental data)
rew.weight_exp([exp0,exp1,exp2],[sim0,sim1,sim2],'example3')


