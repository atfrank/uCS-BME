import os, argparse, shutil, matplotlib
from BME import BME as bme
import pandas as pd
import numpy as np
from py.reweighting import  *
from py.PyRNA import *
from itertools import combinations, permutations
from collections import namedtuple
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.pyplot import hist2d

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Tahoma']

if __name__ == "__main__":        
    parser = argparse.ArgumentParser()
    parser.add_argument("-ex","--experimental", help="CSV file with experimental peaks (peaks those peaks to have names", required = True)
    parser.add_argument("-si","--simulated", help="CSV file with simulated chemical shift table", required = True)
    parser.add_argument("-n1","--name_one", help="Column name for dimenison 1 (C or N) in the experimental peak file", default = "(F1) [ppm]", required = True, type = str)
    parser.add_argument("-n2","--name_two", help="Column name for dimenison 2 (H) in the experimental peak file ", default = "(F2) [ppm]", required = True, type = str)
    parser.add_argument("-e1","--error_one", help="Estimated chemical shift prediction error for dimenison 1", required = True, type = float)
    parser.add_argument("-e2","--error_two", help="Estimated chemical shift prediction error for dimenison 2", required = True, type = float)
    parser.add_argument("-ou","--output", help="Output prefix for generated files", type = str, default = "test")
    parser.add_argument("-im","--imino_only", help="Use only imino simulated chemical shifts (GUA: N1/H1 and URA: N3/H3)", action = "store_true")
    parser.add_argument("-tm","--tmpdir", help="Location used to store auxillary files")
    parser.add_argument("-se","--separation", help="Separation character in CSV", default = ",")
    
    
    # parse command line
    a = parser.parse_args()  
        
    # read in experimental 2D peaks
    expcs_paired = read_peaks(a.experimental, a.name_one, a.name_two, a.separation )

    # read in simulated chemical shifts and convert to 2D peaks
    simcs = read_computed_cs(a.simulated, sep = a.separation)
    if a.imino_only:
        print("NOTE: Only retaining simulated chemical shifts for imino nuclei")
        pairing = {"H1":"N1", "H3":"N3"}
        simcs = pd.concat([simcs[simcs.resname.isin(['GUA']) & simcs.nucleus.isin(['N1', 'H1'])], simcs[simcs.resname.isin(['URA']) & simcs.nucleus.isin(['N3', 'H3'])]])
    else:
        pairing = {"H1'" : "C1'", "H2'" : "C2'", "H3'" : "C3'", "H4'" : "C4'", "H5'" : "C5'", "H5''" : "C5'", "H2" : "C2", "H5" : "C5", "H6" : "C6", "H8" : "C8"}
    simcs_paired = create_paired_data_simulated(simcs, csname = "simcs", debug = False, pairing = pairing)

    # set errors
    error_CN, error_H =  a.error_one, a.error_two

    # create histograms
    expcs_hist, simcs_hist = create_histogram_data(expcs_paired, simcs_paired, xrange = [expcs_paired['F1'].min(),expcs_paired['F1'].max()], xgrid_size = error_CN, yrange = [expcs_paired['F2'].min(), expcs_paired['F2'].max()], ygrid_size = error_H, x_name = "F1", y_name = "F2")

    # write files for BME
    bme_exp = "%s_exp.dat"%(a.output)
    bme_sim = "%s_sim.dat"%(a.output)

    expcs_bme, simcs_bme = write_BME_chemical_shifts_unassigned(expcs_hist, simcs_hist, output_name_exp = bme_exp, output_name_sim = bme_sim, error = 1)
    name = "ucsbme_%s_%s"%(error_CN, error_H)
    bmea = bme.Reweight(name)     
    bmea.load(bme_exp , bme_sim)

    theta, avg_phi = bmea.theta_scan([i for i in np.linspace(1, 100, 200)], tmp_dir = "%s/"%(a.tmpdir))
    optimal_weights, chi2_before, chi2_after, srel = find_weights(bme_exp , bme_sim, theta = theta, name = name)

    # collect data
    n_conformers = len(optimal_weights)
    weights = pd.DataFrame({"model": [i for i in range(1, n_conformers+1)], "w": optimal_weights, "errors": [ np.sum(np.abs(expcs_hist - simcs_hist[i])) for i in range(1, 1+n_conformers)]})
    weights = weights.round(5)
    weights.sort_values(by="w")
    weights.to_csv("%s_ucsbme_weights.csv"%(a.output))
    
    # plot of histogram errors vs weights
    weights.plot(x = "w", y = "errors", kind = "scatter")       
    plt.savefig("%s_ucsbme_histogram_errors_%s_%s.pdf"%(a.output, error_CN, error_H))

    # plot 2D histogram
    xlim = (expcs_paired['F1'].min(), expcs_paired['F1'].max())
    ylim = (expcs_paired['F2'].min(), expcs_paired['F2'].max())
    xgrid_size = error_CN
    ygrid_size = error_H
    xbins = int((xlim[1] - xlim[0])/xgrid_size)
    ybins = int((ylim[1] - ylim[0])/ygrid_size)
    best = weights.model[np.argmax(weights.w.values)]
    diff = np.sum(np.abs(expcs_hist - simcs_hist[best]))
    title = "highest weighted [%s]: Histogram error = %i"%(best, diff)
    plt.figure()  # create a new figure
    ax = plt.subplot(111)
    expcs_paired.plot(ax=ax, kind = "scatter", x="F1", y="F2", c = "red", xlim = xlim, ylim = ylim, title = title)
    simcs_paired[best].plot(ax=ax, kind = "scatter", x="F1", y="F2", c = "blue", xlim = xlim, ylim = ylim)
    ax.set_xticks(np.geomspace(xlim[0],xlim[1],xbins), minor=True )
    ax.set_yticks(np.geomspace(ylim[0],ylim[1],ybins), minor=True )
    ax.grid('on', which='minor', axis='x' )
    ax.grid('on', which='minor', axis='y' )
    ax.tick_params(which='major', length=0, axis='x')
    ax.tick_params(which='major', length=0, axis='y')   
    plt.savefig("%s_ucsbme_2D_histogram.pdf"%(a.output))
    
    # clean up
    os.rename('%s.log'%(name), '%s_ucsbme.log'%(a.output))
    os.rename('%s/crossval_%s.png'%(a.tmpdir, name), '%s_ucsbme_crossval.png'%(a.output))
    os.rename(bme_sim, '%s/bme_sim.txt'%(a.tmpdir))
    os.rename(bme_exp, '%s/bme_exp.txt'%(a.tmpdir))
    