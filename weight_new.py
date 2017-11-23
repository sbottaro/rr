import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

sugars = ["H1H2","H2H3","H3H4"]

def parse():

    parser = argparse.ArgumentParser(description='Reweight here.')   

    parser.add_argument("-o", dest="name",help="output_name",default=None,required=True)
    parser.add_argument("--w", dest="w",help="weight file",default=None,required=True)
    parser.add_argument("--obs-exp", dest="obs_exp",help="observable with experimental value (exp_file1, bcalc_file1, exp_file2, bcalc_file2...)",default=None,required=False,nargs="*")
    parser.add_argument("--obs", dest="obs",help="observable without experimental value (bcalc_file1, bcalc_file2...)",default=None,required=False,nargs="*")
    parser.add_argument("--power", dest="noe_power",help="power for NOE (default=6)",type=float,default=6.0,required=False)
    parser.add_argument("--blocks", dest="blocks",help="number of blocks",type=int,default=1,required=False)

    parser.add_argument("--verbose", dest="verbose",help="verbose",action='store_true',default=False)
    parser.add_argument("--plot-avg", dest="plot",help="plot",action='store_true',default=False)
    parser.add_argument("--plot-hist", dest="plot",help="plot",action='store_true',default=False)
    
    args = parser.parse_args()
    
    return args



def read_weigths(w_file,blocks,verbose):

    name = w_file 

    w0 = []
    w1 = []
    for j in range(blocks):
        if(blocks!=1):
            name = "%s_%d_weights.dat" % (w_file,j)
        else:
            name = "%s_weights.dat" % (w_file,j)
        fh = open(name)
        ww = []
        for line in fh:
            if(len(line.split())==2):
                ww.append(float(line.split()[1]))
            if("#" in line):
                # check exit message from minimization
                msg = line.split()[1]
                if(msg!="True"):
                    print "# Minimization not converged."
                    sys.exit(1)
        ll = len(ww)
        w0t = np.array(ww[:ll/2])
        w1t = np.array(ww[ll/2:])
        w0.append(w0t)
        w1.append(w1t)
        if(verbose):
            print "# read %d/%d weights " % (w0t.shape[0],w1t.shape[0])

    return w0,w1
        
#print np.average(diffsq_total_b), np.average(diffsq_total_a)

def read_exp(file_name):
    
    labels_tmp = []
    exp_avg_tmp = []
    exp_sigma_tmp = []
    bounds_tmp = []
    is_noe = False
            
    fh = open(file_name)
    for ii,line in enumerate(fh):
        if("#" in line):
            if(ii==0):
                # check if is noe
                if("NOE" in line):
                    is_noe = True
                data_type = line.split()[1].split("=")[1]
            else:
                continue
        else:
            ls = line.split()
            # noes
            if(is_noe):
                labels_tmp.append("%s/%s" % (ls[0],ls[1]))
                if(len(ls)==5):
                    if(ls[4] =="UPPER"):
                        exp_avg_tmp.append(float(ls[2]))
                        exp_sigma_tmp.append(float(ls[3]))
                        bounds_tmp.append(+1)
                    else:
                        if(ls[4] =="LOWER"):
                            exp_avg_tmp.append(float(ls[2]))
                            exp_sigma_tmp.append(float(ls[3]))
                            bounds_tmp.append(-1)
                        else:
                            # lab1 lab2 r_min r_avg r_max
                            exp_avg_tmp.append(float(ls[3]))
                            exp_sigma_tmp.append(np.min([float(ls[3])-float(ls[2]),float(ls[4])-float(ls[3])]))
                            bounds_tmp.append(0)
                            #print exp_avg_tmp[-1], exp_sigma_tmp[-1]
                # lab1 lab2 r_avg sigma
                if(len(ls)==4):
                    exp_avg_tmp.append(float(ls[2]))
                    exp_sigma_tmp.append(float(ls[3]))
                    bounds_tmp.append(0)
                            
            # other types of data
            else:
                labels_tmp.append("%s" % (ls[0]))
                exp_avg_tmp.append(float(ls[1]))
                exp_sigma_tmp.append(float(ls[2]))
                bounds_tmp.append(0)
    fh.close()
    return [data_type,labels_tmp, exp_avg_tmp, exp_sigma_tmp,bounds_tmp]
    
def read_obs_exp(files,verbose):

    assert(len(files)%2==0)
    
    labels = []
    exp_avg = []
    exp_sigma = []
    bounds = []
    bcalc = []
    data_type = []
    
    for i,f in enumerate(files):
        
        if(i%2==0):
            
            # read experimental data first
            dd = read_exp(files[i])
            if(verbose):
                print "# read %s with %d experimental datapoints" % (files[i],len(dd[1]))
            labels.append(dd[1])
            exp_avg.append(np.array(dd[2]))
            exp_sigma.append(np.array(dd[3]))
            bounds.append(dd[4])
            data_type.append(dd[0])
            
        # read back-calculated data
        else:
            bcalc_tmp = np.array([[float(x) for x in line.split()[1:]] for line in open(files[i]) if "#" not in line])
            if(verbose):
                print "# read %s with %d/%d back-calculated data" % (files[i],bcalc_tmp.shape[0],bcalc_tmp.shape[1])
            assert(bcalc_tmp.shape[1] == len(dd[1]))
            bcalc.append(bcalc_tmp)
            
    return [data_type,labels,exp_avg,exp_sigma,bounds], bcalc

def read_obs(files):
    bcalc = []
    for f in files:
        bcalc.append(np.array([[float(x) for x in line.split()[1:]] for line in open(f) if "#" not in line]))
    return bcalc
            
def main():
    
    args = parse()
    
    w0,w1 = read_weigths(args.w,args.blocks,args.verbose)
    print "### ", args.w, "###"
    string = ""
    # read observables with experimental data
    if(args.obs!=None):
        data = read_obs(args.obs)

        bins = np.linspace(0,np.max(data[0][:,0]),50)
        ii = np.where(bins<0.75)[0][-1]
        
        start = 0
        stop = 0
        fraction_old = []
        fraction_new = []
        
        for k in range(args.blocks):
            
            stop += len(w0[k])
            ermsd = data[0][start:stop,0]

            hh1,ee1 = np.histogram(ermsd,weights=w0[k],bins=bins,normed=True)
            hh2,ee2 = np.histogram(ermsd,weights=w1[k],bins=bins,normed=True)
            #for k in range(len(hh1)):
            #    print 0.5*(ee1[k]+ee1[k+1]), hh1[k],hh2[k]
            #print " "
            delta = ee1[1]-ee1[0]
            #plt.plot(0.5*(ee1[1:]+ee1[:-1]),hh1,c='k')
            #plt.plot(0.5*(ee1[1:]+ee1[:-1]),hh2,c='r')
            fraction_old.append(100.*np.sum(hh1[:ii])*delta)
            fraction_new.append(100.*np.sum(hh2[:ii])*delta)
            
            start = stop
        print "%40s " % args.obs[0],
        print " A-form before : %6.2f+/-%6.2f  " % (np.average(fraction_old),np.std(fraction_old,ddof=1)/np.sqrt(len(fraction_old))),
        print " A-form after  : %6.2f+/-%6.2f  " % (np.average(fraction_new),np.std(fraction_new,ddof=1)/np.sqrt(len(fraction_new)))
        #plt.savefig("crap.png")
        #exit()
        plt.close()
        
    # read observables with experimental data
    if(args.obs_exp!=None):
        
        exp_data, bcalc_data = read_obs_exp(args.obs_exp,args.verbose)
        # calculate statistics and write to file
        for i in range(len(bcalc_data)):

            start = 0
            stop = 0
            avg_old = []
            avg_new = []
            ncalc = bcalc_data[i].shape[1]
            for k in range(args.blocks):
                stop += len(w0[k])

                if(exp_data[0][i] == "NOE"):
                    bcalc_tmp = [np.power(bcalc_data[i][start:stop,j],-args.noe_power) for j in range(ncalc)]
                else:
                    bcalc_tmp = [bcalc_data[i][start:stop,j] for j in range(ncalc)]
                    
                avg_old.append([np.average(bcalc_tmp[j],weights=w0[k],axis=0) for j in range(ncalc)])
                avg_new.append([np.average(bcalc_tmp[j],weights=w1[k],axis=0) for j in range(ncalc)])

                start = stop
            avg_old = np.array(avg_old)
            avg_new = np.array(avg_new)

            chi_before = []
            chi_after = []
            
            if(exp_data[0][i] == "NOE"):
                avg_old = np.power(avg_old,-1./args.noe_power)
                avg_new = np.power(avg_new,-1./args.noe_power)

            # difference
            for k in range(args.blocks):
                diff1 = avg_old[k]-exp_data[2][i] 
                diff2 = avg_new[k]-exp_data[2][i] 
            
                # set to zero when is upper or lower bound
                idx1 = [1.0 if ((diff1[j] >0 and exp_data[4][i][j] == 1) or (diff1[j] < 0 and exp_data[4][i][j] == -1) or exp_data[4][i][j] ==0 ) else 0.0  for j in range(len(exp_data[4][i]))]
                idx2 = [1.0 if ((diff2[j] >0 and exp_data[4][i][j] == 1) or (diff2[j] < 0 and exp_data[4][i][j] == -1) or exp_data[4][i][j] ==0 ) else 0.0  for j in range(len(exp_data[4][i]))]
                diff1 *= idx1
                diff2 *= idx2
                
                rmsd1 = np.sqrt(np.average(diff1**2))
                rmsd2 = np.sqrt(np.average(diff2**2))
                
                chi1 = np.average(diff1**2/exp_data[3][i]**2)
                chi2 = np.average(diff2**2/exp_data[3][i]**2)

                #for k in range(len(idx1)):
                #    print exp_data[1][i][k], avg_old[k],avg_new[k], exp_data[2][i][k]
                #print chi1, chi2, diff1.shape
                #print rmsd1, rmsd2
                # splin in backbone and sugar if Jcouplings
                if(exp_data[0][i] == "NOE"):
                    chi_before.append([chi1])
                    chi_after.append([chi2])
                    
                if(exp_data[0][i] == "JCOUPLINGS"):
                    idx_bb = [l for l,lab in enumerate(exp_data[1][i]) if lab.split("-")[1] not in sugars ]
                    idx_su = [l for l,lab in enumerate(exp_data[1][i]) if lab.split("-")[1] in sugars ]
                
                    rmsd1_bb = np.sqrt(np.average(diff1[idx_bb]**2))
                    rmsd1_su = np.sqrt(np.average(diff1[idx_su]**2))
                    
                    rmsd2_bb = np.sqrt(np.average(diff2[idx_bb]**2))
                    rmsd2_su = np.sqrt(np.average(diff2[idx_su]**2))
                    
                    chi1_bb = np.average(diff1[idx_bb]**2/exp_data[3][i][idx_bb]**2)
                    chi1_su = np.average(diff1[idx_su]**2/exp_data[3][i][idx_su]**2)
                    
                    chi2_bb = np.average(diff2[idx_bb]**2/exp_data[3][i][idx_bb]**2)
                    chi2_su = np.average(diff2[idx_su]**2/exp_data[3][i][idx_su]**2)
                    
                    chi_before.append([chi1_bb, chi1_su])
                    chi_after.append([chi2_bb, chi2_su])
                
                
                if(exp_data[0][i] == "CS"):
                    idx_h = [l for l,lab in enumerate(exp_data[1][i]) if "H" in lab ]
                    idx_c = [l for l,lab in enumerate(exp_data[1][i]) if "H" not in lab ]
                    
                    rmsd1_h = np.sqrt(np.average(diff1[idx_h]**2))
                    rmsd1_h = np.sqrt(np.average(diff1[idx_c]**2))
                    
                    rmsd2_h = np.sqrt(np.average(diff2[idx_h]**2))
                    rmsd2_c = np.sqrt(np.average(diff2[idx_c]**2))
                    
                    chi1_h = np.average(diff1[idx_h]**2/exp_data[3][i][idx_h]**2)
                    chi1_c = np.average(diff1[idx_c]**2/exp_data[3][i][idx_c]**2)
                    
                    chi2_h = np.average(diff2[idx_h]**2/exp_data[3][i][idx_h]**2)
                    chi2_c = np.average(diff2[idx_c]**2/exp_data[3][i][idx_c]**2)
                    chi_before.append([chi1_h, chi1_c])
                    chi_after.append([chi2_h, chi2_c])

            print "%40s " % args.obs_exp[2*i],
            chi_before = np.array(chi_before)
            chi_after = np.array(chi_after)
            spec = ["BB","SU"]
            for j in range(chi_before.shape[1]):
                #print np.average(chi_before[:,j]), chi_before.shape
                print " chi2 before %s : %6.2f+/-%6.2f  " % (spec[j],np.average(chi_before[:,j]),np.std(chi_before[:,j],ddof=1)/np.sqrt(chi_before.shape[0])),
                print " chi2 after %s : %6.2f+/-%6.2f  " % (spec[j],np.average(chi_after[:,j]),np.std(chi_after[:,j],ddof=1)/np.sqrt(chi_after.shape[0])),
            print ""

            string += "#%40s \n" % args.obs_exp[2*i]
            for j in range(len(exp_data[1][i])):
                string += "%20s " % exp_data[1][i][j]
                string += "%8.4f %8.4f " % ( exp_data[2][i][j],exp_data[3][i][j])
                string += "%8.4f %8.4f " % (np.average(avg_old[:,j]), np.std(avg_old[:,j],ddof=1)/np.sqrt(args.blocks))
                string += "%8.4f %8.4f \n" % ( np.average(avg_new[:,j]), np.std(avg_new[:,j],ddof=1)/np.sqrt(args.blocks))

    fhw = open(args.name + ".dat","w")
    fhw.write(string)
    fhw.close()
        #print "###############"
            #for k in range(len(exp_data[4][i])):
            #    print exp_data[1][i][k],exp_data[2][i][k],avg_old[k],diff1[k],avg_new[k],diff2[k],idx1[k]
            #print "###############"
    

    
if __name__ == "__main__":
    main()



