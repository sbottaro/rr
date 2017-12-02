import numpy as np
from scipy import optimize
import argparse
import sys



class Reweight:


    def __init__(self,exp_files,sim_files,noe_power=6.0):

        self.exp_data = []
        self.sim_data = []
        self.data_source = []
        self.labels = []
        self.bounds = []
        
        # this is to avoid numerics in NOE
        self.noe_factor = 0.4
        print "###### INITIALIZATION ########"
        
        # these are known data types and error priors
        self.exp_types = ["NOE","JCOUPLINGS","CS","SAXS"]
        self.prior_types = ["GAUSS"]

        assert len(exp_files)==len(sim_files), "# Error. Number of experimental (%d) and simulation (%d) files must be equal" % (len(exp_files),len(sim_files))
        for k in range(len(exp_files)):
            nexp = self.read_exp(exp_files[k],noe_power)
            if(self.data_source[-1]=="NOE"):
                nsim,nframes = self.read_sim(sim_files[k],noe_power)
            else:
                nsim,nframes = self.read_sim(sim_files[k])
            assert nexp ==nsim, "# Number of rows in %s not equal to number of columns in %s" % (nexp,nsim)


        self.exp_data = np.array(self.exp_data)
        self.sim_data = np.array(self.sim_data)
        #self.bounds = list(np.array(self.bounds)[ii])
        #self.data_source = list(np.array(self.data_source)[ii])
        
        print "# Exp data shape ",self.exp_data.shape
        print "# Simulation data shape ",self.sim_data.shape
        
        # check that each expt point is within simulation data. Otherwise there might be something wrong?
        mins = np.min(self.sim_data,axis=0)
        maxs = np.max(self.sim_data,axis=0)
        for k in range(len(self.exp_data)):
            if(self.bounds[k][0] == None and self.bounds[k][1] == None):
                if(self.exp_data[k,0] > maxs[k]):
                    print "# Warning: expt average %s=%-10.4f is larger than maximum value in simulation %-10.4f"\
                        % (self.labels[k],self.exp_data[k,0],maxs[k])
                if(self.exp_data[k,0] < mins[k]):
                    print "# Warning: expt average %s=%-10.4f is smaller than minimum value in simulation %-10.4f"\
                        % (self.labels[k],self.exp_data[k,0],mins[k])
            else:
                if(self.data_source[k] == "NOE"):
                    # lower boundary - check that at least some points are above 
                    exp_r = np.power(self.exp_data[k,0],-1./noe_power)/self.noe_factor
                    max_r = np.power(mins[k],-1./noe_power)/self.noe_factor
                    min_r = np.power(maxs[k],-1./noe_power)/self.noe_factor
                    
                    if(self.bounds[k][0] == 0.0):
                        if( max_r < exp_r):
                            print "# Warning: maximum value in simulation %-10.4f lower than exp. lower bound %s=%-10.4f"\
                                % (max_r,self.labels[k],exp_r)
                    if(self.bounds[k][1] == 0.0):
                        if( min_r > exp_r):
                            print "# Warning: minimum value in simulation %-10.4f larger than exp. upper bound %s=%-10.4f"\
                                % (max_r,self.labels[k],exp_r)
        
                else:
                    # lower boundary - check that at least some points are above 
                    if(self.bounds[k][0] == 0.0):
                        if( maxs[k] < self.exp_data[k,0]):
                            print "# Warning: expt lower boundary %s=%-10.4f is larger than maximum value in simulation %-10.4f"\
                                % (self.labels[k],self.exp_data[k,0],maxs[k])
                    if(self.bounds[k][1] == 0.0):
                        if( mins[k] < self.exp_data[k,0]):
                            print "# Warning: expt upper boundary %s=%-10.4f is larger than minimum value in simulation %-10.4f"\
                                % (self.labels[k],self.exp_data[k,0],maxs[k])
        print "###### OK. ########"
        self.w0 = np.ones(self.sim_data.shape[0])/self.sim_data.shape[0]
        
    # read experimental data
    def read_exp(self,filename,noe_power=6.0):

        fh = open(filename)
        first = fh.readline().split()
        
        # do some checks on first line
        assert first[0] == "#", "Error. First line of exp file %s must be in the format \# DATA=XXX PRIOR=XXX" % filename
        data_type = (first[1].split("=")[1]).strip()
        prior_type = (first[2].split("=")[1]).strip()
        assert data_type in self.exp_types , "Error. DATA must be one of the following: %s " % (self.exp_types)
        assert prior_type in self.prior_types , "Error. PRIOR must be one of the following: %s " % (self.prior_types)

        ln = 0
        # read data - NOE is a special case
        if(data_type=="NOE"):
            for line in fh:
                # skip commented out data
                if("#" in line): continue

                ln += 1
                lsplit = line.split()
                self.labels.append( "%s/%s"  % (lsplit[0],lsplit[1]))
                self.data_source.append(data_type)
                
                # if lenght is 4, third and fourth columns are  r_avg and sigma
                if(len(lsplit)==4):
                    avg = float(lsplit[2])
                    sigma = float(lsplit()[3])
                    self.bounds.append([None,None])

                # if lenght is 5
                elif(len(lsplit)==5):

                    # upper/lower bound
                    if(lsplit[4] == "UPPER" or lsplit[4] == "LOWER"):
                        avg = float(lsplit[2])
                        sigma = float(lsplit[3])
                        # upper/lower is inverted because of 1/r^n dependence of NOE
                        if(lsplit[4] == "UPPER"): self.bound.append([None, 0.0])
                        else: self.bounds.append([0.0,None])
                        
                    # rmin,r,rmax
                    else:
                        r_low = float(line.split()[2])
                        avg = float(line.split()[3])
                        r_up = float(line.split()[4])
                        assert(r_up >= avg)
                        assert(avg >= r_low)
                        # sigma is taken as the minimum distance from boundary
                        sigma = min(r_up-avg, avg-r_low)
                        self.bounds.append([None, None])
                else:
                    print "# Fatal error.Something wrong in your expt. data"
                    sys.exit(1)
                # now convert distance to intensities. take the -n power and multiply by factor to avoid numerical problems
                avg_int = np.power(self.noe_factor*avg,-noe_power)
                sigma_int = (noe_power*avg_int*sigma/(avg))**2
                self.exp_data.append([avg_int,sigma_int])
                
        else:
            for line in fh:
                # skip commented out data
                if("#" in line): continue
                
                lsplit = line.split()
                ln += 1
                self.labels.append( "%s"  % (lsplit[0]))
                self.data_source.append(data_type)
                self.exp_data.append([float(lsplit[1]),float(lsplit[2])*float(lsplit[2])])
                                        
                if(len(lsplit)==4):
                    if(lsplit[4] == "LOWER"):
                        self.bounds.append([None, 0.0])
                    elif (lsplit[3] == "UPPER"):
                        self.bounds.append([0.0,None])
                    else:
                        print "# Fatal error.Something wrong in your expt. data"
                        sys.exit(1)
                else:
                    self.bounds.append([None, None])
        print "# Read %d %s experimental datapoints" % (ln,data_type)
        assert(len(self.labels) == len(self.bounds))
        assert(len(self.labels) == len(self.data_source))
        assert(len(self.labels) == len(self.exp_data))
        
        return ln

    def read_sim(self,filename,noe_power=None):
        
        data = np.array([[float(x) for x in line.split()[1:]] for line in open(filename) if("#" not in line)])
        if(noe_power != None):
            data = np.power(self.noe_factor*data,-noe_power)
        print "# Read %d simulated datapoints from %d frames" % (data.shape[1],data.shape[0])
        self.sim_data.extend(list(data))
        return data.shape[1],data.shape[0]

    def set_w0(self,filename):
        w0 = np.array([float(line.split()[0]) for line in open(filename) if("#" not in line)])
        print "# Set non-uniform initial weights from file. Sum=", np.sum(w0)
        print "# not tested yet."
        self.w0 = w0


    def chi_squared(self,ww):

        bcalc = np.sum(ww[:,np.newaxis]*self.sim_data,axis=0)
        diff = bcalc-self.exp_data[:,0]
        idxs = [0.0 if(( diff[i] >= 0 and self.bounds[i][1] == 0) or (diff[i] <= 0 and self.bounds[i][0] == 0)) else 1.0 for i in range(len(self.bounds))]
        diff *= idxs
        return np.sum(diff*diff/self.exp_data[:,1])

                
    def optimize(self,theta,method="MAXENT"):

        def func_maxent_gauss(lambdas):
            
            # weights
            ww  = self.w0*np.exp(-np.sum(lambdas[np.newaxis,:]*self.sim_data,axis=1))
            # normalization 
            zz = np.sum(ww)
            ww /= zz
            # new averages
            avg = np.sum((ww[:,np.newaxis]*self.sim_data), axis=0)
            
            # errors are rescaled by factor theta
            err = (self.theta)*(self.exp_data[:,1])
            # gaussian integral
            eps2 = 0.5*(self.theta)*np.sum((lambdas*lambdas)*err)
            # experimental value 
            sum1 = np.dot(lambdas,self.exp_data[:,0])
            fun = sum1 + eps2+ np.log(zz)
            # gradient
            jac = self.exp_data[:,0] + self.theta*lambdas*err - avg

            return  fun,jac

        
        
        def func_ber_gauss(w):
            
            bcalc = np.sum(w[:,np.newaxis]*self.sim_data,axis=0)
            diff = bcalc-self.exp_data[:,0]
            idxs = [0.0 if(( diff[i] >= 0 and self.bounds[i][1] == 0) or (diff[i] <= 0 and self.bounds[i][0] == 0)) else 1.0 for i in range(len(self.bounds))]
            diff *= idxs
            chi2_half =  0.5*np.sum(((diff**2)/(self.exp_data[:,1])))

            idxs = np.where(w>1.0E-50)
            srel = self.theta*np.sum(w[idxs]*np.log(w[idxs]/self.w0[idxs]))
            return chi2_half+srel

        def srel(w):
            idxs = np.where(w>1.0E-50)
            return np.sum(w[idxs]*np.log(w[idxs]/self.w0[idxs]))
            
        self.theta = theta
        # first, calculate initial chi squared and RMSD
        
        if(method=="MAXENT"):
            opt={'maxiter':50000,'disp': False,'gtol':1.0e-20}
            #opt={'maxiter':50000,'disp': False,'ftol':1.0e-50}
            meth = "L-BFGS-B"

            lambdas=np.zeros(self.exp_data.shape[0])

            result = optimize.minimize(func_maxent_gauss,lambdas,options=opt,method=meth,jac=True,bounds=self.bounds)
            
            w_opt = np.exp(-np.sum(result.x[np.newaxis,:]*self.sim_data,axis=1))
            w_opt /= np.sum(w_opt)
            
        if(method=="BER"):
            opt={'maxiter':2000,'disp': True,'ftol':1.0e-20}
            cons = {'type': 'eq', 'fun':lambda x: np.sum(x)-1.0}
            bounds = [(0.,1.)]*len(self.w0)  
            meth = "SLSQP"
            print "# Bayesian Ensemble Refinement. Useful for testing purposes and "
            print "# When the number of experimental data is larger than the number of samples."
            w_init = np.ones(self.sim_data.shape[0])/self.sim_data.shape[0]
            result = optimize.minimize(func_ber_gauss,w_init,constraints=cons,options=opt,method=meth,bounds=bounds)
            w_opt = result.x

        chi_sq = self.chi_squared(self.w0)/len(self.exp_data)
        print "Initial average chi squared %10.4f, srel %10.4f " % (chi_sq, srel(self.w0))

        chi_sq = self.chi_squared(w_opt)/len(self.exp_data)
        print "Final average chi squared   %10.4f, srel %10.4f " % (chi_sq, srel(w_opt))
        self.w_opt = w_opt


    def weight_exp(exp_files, sim_files):
        
'''        
        def func_maxent_laplace(lambdas):

            # weights
            ww  = self.w0*np.exp(-np.sum(lambdas[np.newaxis,:]*self.sim_data,axis=1))
            # normalization 
            zz = np.sum(ww)
            ww /= zz
            # new averages
            avg = np.sum((ww[:,np.newaxis]*self.sim_data), axis=0)
            # errors are rescaled by factor theta
            err = theta*self.exp_data[:,1]
            
            # integral error
            eps2 = np.sqrt(2.*err)/(1.-0.5*(err*lambdas**2))
            eps2 = np.sum(np.log(eps2))
    
            # experimental value 
            sum1 = np.dot(lambdas,self.exp_data[:,0])
            fun = sum1 + eps2+ np.log(zz)
            # gradient
            lap0 = -lambdas*err
            lap1 = lap0/(1.+0.5*(lambdas*lap0))
            jac = exp_avg[:,0] - lap1 - avg

            return  fun,jac
'''
