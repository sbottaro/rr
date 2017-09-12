import sys
import numpy as np
import mdtraj as md
from scipy import optimize
import argparse
import re

noe_factor = 0.4


exp_types = ["NOE","JCOUPLINGS","CS","SAXS"]
prior_types = ["GAUSS","LAPLACE"]


def parse():

    parser = argparse.ArgumentParser(description='Reweight here.')   

    parser.add_argument("-o", dest="name",help="output_name",default=None,required=True)

    parser.add_argument("--exp", dest="exp_file",help="experimental_data",default=None,required=True,nargs="*")
    parser.add_argument("--bcalc", dest="bcalc_file",help="back-calculated quantity",default=None,required=True,nargs="*")

    #parser.add_argument("-n", dest="nsamples",help="number of samples ",type=int,default=None,required=False)
    #parser.add_argument("--start", dest="start",help="start",required=False,type=int,default=-1)
    #parser.add_argument("--stride", dest="stride",help="stride",required=False,type=int,default=-1
    parser.add_argument("--power", dest="noe_power",help="power for NOE (default=6)",type=float,default=6.0,required=False)
    parser.add_argument("--theta", dest="theta",help="theta value",type=float,default=1.0,required=False)
    parser.add_argument("--blocks", dest="blocks",help="number of blocks",type=int,default=1,required=False)
    parser.add_argument("--verbose", dest="verbose",help="verbose",action='store_true',default=False)

    #parser.add_argument("--seed", dest="seed",help="seed",type=int,default=0,required=False)
    
    args = parser.parse_args()
    
    return args


def read_sim(f_sim):

    frames = []
    for i,f in enumerate(f_sim):
        fh = open(f)
        tmp = []
        for line in fh:
            if("#" not in line):
                frames.append(float(line.split()[0]))
                tmp.append([float(x) for x in line.split()[1:]])
                #data_source_new.append(data_source[i])
        tmp = np.array(tmp)
        if(i==0):
            data = np.array(tmp)
        else:
            assert tmp.shape[0] == data.shape[0]
            data = np.concatenate((data,tmp),axis=1)
        fh.close()
    return frames, np.array(data)
    
def read_exp(f_exp):

    labels = []
    exp_data = []
    bound = []
    data_source = []
    for f in f_exp:

        fh = open(f)
        for i,line in enumerate(fh):
            if(len(line.split())<3): continue
            # Read header
            if(i==0):
                assert(line.split()[0]=="#")
                exp_type = line.split("DATA=")[1].split()[0]
                if(exp_type not in exp_types):
                    print "# I don't know ", exp_type, ". Exiting"
                    sys.exit(1)
                    
                prior_type = line.split("PRIOR=")[1].split()[0]
                if(prior_type not in prior_types):
                    print "# I don't know ", prior_type, ". Exiting"
                    sys.exit(1)
                
            else:
                
                # NOES 
                if(exp_type=="NOE"):
                    data_source.append(exp_type)
                    labels.append([x for x in line.split()[:2]])
                    if(len(line.split())==4):
                        r = float(line.split()[2])
                        sigma = float(line.split()[3])
                        bound.append([None, None])

                        #exp_data.append([r,sigma**2])
                        
                    if(len(line.split())==5):
                        
                        try:
                            r_low = float(line.split()[2])
                            r = float(line.split()[3])
                            r_up = float(line.split()[4])
                            assert(r_up >= r)
                            assert(r >= r_low)
                            # sigma is taken as the minimum distance from boundary
                            sigma = min(r_up-r, r-r_low)
                            bound.append([None, None])
                        except:
                            
                            r = float(line.split()[2])
                            sigma = float(line.split()[3])
                            lim = line.split()[4]
                            if(lim=="UPPER"):
                                # only negative lagrange multipliers
                                bound.append([None, 0.0])
                                
                            else:
                                if(lim=="LOWER"):
                                    # only positive lagrange multipliers
                                    bound.append([0.0,None])
                                else:
                                    print "# invalid bound. Must be UPPER, LOWER"
                                    sys.exit(1)
                    
                    rpow = np.power(noe_factor*r,-noe_power)
                    err = (noe_power*rpow*sigma/(r))**2
                    exp_data.append([rpow,err])
                    
                    
                # chemical shift
                if(exp_type=="CS" or exp_type=="JCOUPLINGS" or exp_type == "SAXS"):
                    if("#" in line):
                        continue
                    else:
                        data_source.append(exp_type)
                        labels.append(line.split()[0])
                        val = float(line.split()[1])
                        sigma = float(line.split()[2])
                        exp_data.append([val,sigma**2])
                        bound.append([None, None])
                        
        fh.close()
    return labels,np.array(exp_data),bound,data_source


def read_labels(top,labels):
    
    # read in trajectory
    labels_idx = [None]*len(labels)
    for res in top.residues:
        for at in res.atoms:
            for j,ll in enumerate(labels):
                if(str(at)==ll): labels_idx[j] = at.index
                
    # check that all atoms are found
    assert(all(labels_idx) != None)
    labels_idx_p = [[labels_idx[j],labels_idx[j+1]] for j in range(0,len(labels_idx),2)]
    print "# pairs in trajectory correctly assigned"
    return labels_idx_p

def chi(weights,bounds):

    bcalc = np.sum((weights[:,np.newaxis]*sim_data_block),axis=0)
    diff = (bcalc-exp_avg[:,0])
    # account for boundaries
    idxs = [0.0 if(( diff[i] >= 0 and bounds[i][1] == 0) or (diff[i] <= 0 and bounds[i][0] == 0)) else 1.0 for i in range(len(bounds))]
    diff *= idxs
    chi = np.sum(((diff**2)/(exp_avg[:,1])))
    return chi,bcalc


def func_gh(w,chi):

    # ATT! assume uniform input weights
    num = w.shape[0]
    bcalc = np.sum((w[:,np.newaxis]*sim_data_block),axis=0)
    chi_half = 0.5*chi
    idxs = np.where(w>1.0E-30)
    srel = theta*np.sum(w[idxs]*np.log(num*w[idxs]))
    return chi_half+srel


def func_maxent_gauss(lambdas):

    # weights
    ww  = np.exp(-np.sum(lambdas[np.newaxis,:]*sim_data_block,axis=1))
    # normalization 
    zz = np.sum(ww)
    ww /= zz
    # new averages
    avg = np.sum((ww[:,np.newaxis]*sim_data_block), axis=0)

    # errors are rescaled by factor theta
    err = theta*exp_avg[:,1]
    # gaussian integral
    eps2 = 0.5*theta*np.sum((lambdas**2)*(err))
    # experimental value 
    sum1 = np.dot(lambdas,exp_avg[:,0])
    fun = sum1 + eps2+ np.log(zz)
    # gradient
    jac = exp_avg[:,0] + theta*lambdas*err - avg

    return  fun,jac


def func_maxent_laplace(lambdas):

    # weights
    ww  = np.exp(-np.sum(lambdas[np.newaxis,:]*sim_data_block,axis=1))
    # normalization 
    zz = np.sum(ww)
    ww /= zz
    # new averages
    avg = np.sum((ww[:,np.newaxis]*sim_data_block), axis=0)
    # errors are rescaled by factor theta
    err = theta*exp_avg[:,1]

    # integral error
    eps2 = np.sqrt(2.*err)/(1.-0.5*(err*lambdas**2))
    eps2 = np.sum(np.log(eps2))
    
    # experimental value 
    sum1 = np.dot(lambdas,exp_avg[:,0])
    fun = sum1 + eps2+ np.log(zz)
    # gradient
    lap0 = -lambdas*err
    lap1 = lap0/(1.+0.5*(lambdas*lap0))
    jac = exp_avg[:,0] - lap1 - avg

    return  fun,jac



    
###################################



#def chi1(weights):#
#
#    ss = np.sum((weights[:,np.newaxis]*sim_data_block),axis=0)
#    bcalc=np.power(ss,-1./noe_power)
#    return bcalc

def rmsd(weights):
    
    bcalc = np.sum((weights[:,np.newaxis]*sim_data_block),axis=0)
    #bcalc=np.power(ss,-1./6.)
    #return np.sqrt(np.average(((bcalc-exp_avg)**2)))
    return np.sqrt(np.average(((bcalc-exp_avg[:,0])**2)))

def rmsd1(weights):
    
    ss=  np.sum((weights[:,np.newaxis]*sim_data_block),axis=0)
    bcalc=np.power(ss,-1./noe_power)
    #return np.sqrt(np.average(((bcalc-exp_avg)**2)))
    return np.sqrt(np.average(((bcalc-np.power(exp_avg[:,0],-1./noe_power))**2)))

def srel(weights):

    # ATT! assume uniform input weights
    num = weights.shape[0]
    idxs = np.where(weights>1.0E-15)
    return  np.sum(weights[idxs]*np.log(num*weights[idxs]))


def stats(weights,bounds,source):
    

    chi2,crap = chi(weights,bounds)
    fun = func_gh(weights,chi2)
    sreli = srel(weights)
    norm = np.sum(weights)
    #stri = "# %10s # \n" % (s)
    stri = "# %10s %10s %10s %10s\n" % ("function","chi","srel","normw")
    stri += "# %10.5f %10.5f %10.5f %10.5f \n" % (fun,chi2,sreli,norm)
    stri += "# %20s %8s %8s %8s \n" % ("Label","predict","experim","sigma")
    ss = np.sum((weights[:,np.newaxis]*sim_data_block),axis=0)
    for p in range(exp_avg.shape[0]):
        #print  "%10s %10s %8.5f %8.5f\n" % (labels[p][0],labels[p][1],np.power(bcalcn[p],-1./noe_power), exp_avg[p,0])
        #print  "%10s %10s %8.5f %8.5f\n" % (labels[p][0],labels[p][1],bcalcn[p], np.power(exp_avg[p,0],-1./noe_power))
        
        if(source[p]=="NOE"):
            bcalc = np.power(ss[p],-1./noe_power)/noe_factor
            stri += "%10s-%10s %8.5f %8.5f\n" % (labels[p][0],labels[p][1],bcalc, np.power(exp_avg[p,0],-1./noe_power)/noe_factor)
        else:
            stri += "%10s %8.5f %8.5f\n" % (labels[p][0],ss[p], exp_avg[p,0])
    stri += " \n"    
    stri += " \n"    
    return stri

def write(name,weights,bounds,indeces,source,mode="a"):
    
    fho = open("stats_%s.dat" % (name),mode)
    fho.write(stats(weights,bounds,source))
    fho.close()

    # print weights
    fho = open("weights_%s.dat" % (name),mode)  
    for f in range(len(weights)):
        #print 0.5*(ee[f+1]+ee[f]), hh[f], hhw[f]
        fho.write("%10d %e \n" % (indeces[f],weights[f]))
    fho.write("\n")
    fho.write("\n")
    fho.close()

    
####################### MAIN #########################

def main():

    global exp_avg
    global labels
    global theta
    global sim_data_block
    
    # Parse options
    args = parse()
    theta = args.theta
    global noe_power
    noe_power =  args.noe_power

    # read experimental data
    labels,exp_avg,bounds,data_source  = read_exp(args.exp_file)

    # read back-calculated
    frames, sim_data = read_sim(args.bcalc_file)
    print "A1", sim_data.shape
    print "A2", len(data_source)
    # scale NOEs by some factor to avoid numerical issues and calculate the -noe_power
    sim_data = np.array([np.power(noe_factor*sim_data[:,j],-noe_power) if data_source[j]=="NOE" else sim_data[:,j] for j in range(len(data_source))]).T
    print "# ", sim_data.shape, exp_avg.shape

    assert(sim_data.shape[1] == exp_avg.shape[0])

    #print bounds

    opt={'maxiter':50000,'disp': args.verbose,'gtol':1.0e-15}
    meth = "L-BFGS-B"

    # initialize with zero lambdas and zero weights
    #lambdas=np.zeros(exp_avg.shape[0])
    lambdas=np.zeros(exp_avg.shape[0])

    len_b = sim_data.shape[0]/args.blocks
    for b in range(args.blocks):
        start = b*len_b
        stop = (b+1)*len_b
        if(b==args.blocks-1):
            stop = sim_data.shape[0]

        name = "%s" % (args.name)
        if(args.blocks>1):
            name = "%s_%d" % (args.name,b)

        sim_data_block = sim_data[start:stop]
        w_init = np.exp(-np.sum(lambdas[np.newaxis,:]*sim_data_block,axis=1))
        w_init /= np.sum(w_init)
        write(name,w_init,bounds,frames,data_source,mode="w")
        
        # write some statistics to screen 
        chi0,bcalc0 = chi(w_init,bounds)
        print "# Initial conditions"
        print "# Chi squared: %10.6f" %  chi0 
        #print "# RMSD:        %10.6f" %  rmsd(w_init)
        print "# Srel:        %10.6f" %  srel(w_init) 
        print "# Functional:  %10.6f" %  func_gh(w_init,chi0)
        
        result1 = optimize.minimize(func_maxent_gauss,lambdas,options=opt,method=meth,jac=True,bounds=bounds)
        
        w_opt = np.exp(-np.sum(result1.x[np.newaxis,:]*sim_data_block,axis=1))
        w_opt /= np.sum(w_opt)
        
        #print result1.x
        eff = np.exp(np.sum(-w_opt*np.log(w_opt)))/len(w_opt)
        # write some statistics to screen 
        chi_opt,bcalc_opt = chi(w_opt,bounds)
        
        print "#  After minimization"
        print "# Chi squared: %10.6f" %  chi_opt 
        #print "# RMSD:        %10.6f" %  rmsd(w_opt)
        print "# Srel:        %10.6f" %  srel(w_opt) 
        print "# Functional:  %10.6f" %  func_gh(w_opt,chi_opt)
        print "# Fraction of used frames %4.3f" % eff
        print "# ", result1.success, result1.message
        #print "# AAAA  %10.6f %10.6f %10.6f %10.6f %4.3f" % (chi_opt ,rmsd1(w_opt),srel(w_opt), func_gh(w_opt),eff)
        # and to output file
        write(name,w_opt,bounds,frames,data_source,mode="a")
        
        fho = open("weights_%s.dat" % name,"a")  
        fho.write("# %s %s \n" % (result1.success, result1.message))
        fho.close()





    
    
if __name__ == "__main__":
    main()

