# Analysis codes, data files and plotting scripts for the paper arXiv:2202.05981 [nucl-th]
---
This repository contains the source code and data for the paper titled "Scaling properties of exclusive vector meson production cross section from gluon saturation". The usage of several of these macros, in particular the source code, depends on a correctly built version of Sartre 1.33, available online at https://sartre.hepforge.org/

Within this repository, there are four separate main directories. These are, along with their contents...

- **./src_code**
  - This directory contains two macros which are meant to be compiled via Sartre in a similar fashion to those, such as sartreMain.cpp, in Sartre's 'examples' directory. Both of these macros are fed a *sartreRuncard.txt* script in the command line (ex: ./sartreWhitepaper.cpp sartreRuncard.txt) to set event wide parameters, such as beam energies, ion species, and decay VMs. 
  - The first of these macros, *sartreQUICKSMEAR3.cpp*, is responsible for generating several "A" and "Q2" scaled cross section/differential cross section 1-d histograms, at both t=0 and several t values. Event-by-event, particle kinematics are passed through a hand written experimental filter (accounting for EIC Handbook/EIC Yellow Report detector resolutions and acceptances). These histograms are all saved into a single .root file.
  - In a similar way, the send macro, *sartreWhitepaper.cpp* saves histograms corresponding to the cuts and binning of the diffractive peak cross section plots, seen both in the EIC whitepaper and in our report.
  
- **./t_equal_zero_plots**
  - This directory contains several subdirectories with self-explanatory titles. The *./data* folder stores *.root* files compiled from multiple executions of *./src_code/sartreQUICKSMEAR3.cpp* with different runcards. The *./new_paper_plots* folder stores plots containing Monte Carlo cross section overlays (to compare with the pseudodata). These plots can be seen in Appendix C of the paper. The *./paper_plots* folder stores plots with only pseudodata data points, no Monte Carlo. It can be quickly noted that these plots are not the final visual result within the paper (see Figures 5-9 and 14-15) yet nevertheless are created with the same histograms in *./data*. The *./plots/* folder is fairly old and can be ignored. Lastly, the plotting macro *createPlots.C* macro, executed via root, pulls in the histograms from the *./data* folder and compiles them as figures for the user.
  
- **./t_not_equal_zero_plots**
  - Similar to *./t_equal_zero_plots*, this directory simply creates figures for the t-scaling plots in Figures 10-13.

- **./YellowReport_plots/**
  - This directory's structure is similar to the above, but in this case the *./data* folder's histograms are stored in such a way to make the plots in figures 3-4. The *./data* directory contains *.root* files created by the *sartreWhitepaper.cpp* program. The *./has_truth_data* directory is identical to *./data*, with the subtle difference that histograms containing the suffix *_0* in their TObject name have no smearing whatsoever (hence they are used for Monte Carlo comparison). The directory *./newplots* is where figure 4 is generated (using the true_whitepaperplot.C script), and *./plots* is where figure 3 is generated (using the whitepaperplot.C script).
  
If you have any questions about the workings of these macros, how the data is generated, or how to compile the source code, feel free to contact me at my email (gregory.matousek@duke.edu)


