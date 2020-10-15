# cluster_fusion_precision
This repository corresponds to code corresponding to "Estimating Multiple Precision Matrices with Cluster Fusion Regularization". by Price, Molstad, and Sherwood.


This folder has 4 main subfolders including 2 files plus this readme. 

The files cluster_da.R and rda.R are implementations of our methods and a corresponding RDA method that may be used inside of each simulation.  The cluster_da.R file will be used to create an R Package upon publication of this work.  

We note most if not all of these methods were run on the HPC system at West Virginia University, thus the code base is built to interact with array jobs in that way.  We have done our best to take out the file system requirements of the HPC system, but the code is built to interact that system.  Hence we have included the .RData files and .Rout files where applicable.  We note that any .pbs file is used on the HPC system to execute the job and will define the R environment used and the .txt files provided provide the parameters that are used in the array jobs.  The RJob.o files provide outputs of each element of the array.  We note that these are quite large and do not recommend reproducing this entire simulation without access to a large HPC system.  

The folder graph_sim contains the code, results, and output files for the simulations on Gaussian Graphical models.  Inside we have two subfolders that contain results for accuracy and estimation, which are the results described in the manuscript.  While the speed and convergence folder contains the accuracy and estimation information contained in  the file.  Each subfolder contains an outputfiles folder that contain the information and execution of the HPC environment. The a or b in the file names comes from if it is replications 1 to 25 (a) or replications 26:50.  This is due to time and space constraints on the HPC system and time run over a given number of minuets.  In the timing files we add timing before a or b.  Note that setting 1,2,4 correspond to simulations in the manuscript while setting 4 is the simulation that is found in the supplemental material.  

The qda_sim folder has two folders sim1 and sim3.  The sim1 folder contains the simulation study we explore in the supplemental material, while sim3 folder contains the simulation we explore in the manuscript.  Some files were run from a terminal batch mode and others the HPC system at West Virginia University. The batch mode runs were run with "R CMD BATCH --vanilla filename.R".  The output files are contained within.  

The libras subfolder contains two files.  The first is the .R file that analyzes the libra data and the second is the .Rout containing the results from that run.  Since 5 fold cross validation is used and the amount of tuning parameters investigated, the time to run this code is considerable.  This is part of the reason we put the .Rout files in for transparency.  

The pulmonary_graph folder contains 4 files.  Two .R files, one .RDS file, and one .RData file.  The two .R files implement the competing methods and provide results. The .RDS file is the raw data used for the analysis.  The .RData file contain the best results from PCEN.  


