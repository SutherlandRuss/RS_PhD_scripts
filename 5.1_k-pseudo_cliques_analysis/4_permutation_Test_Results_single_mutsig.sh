##############################################
#!/bin/sh
#$ -S /bin/sh
#
##$ -cwd
#$ -N mutsig_perm_results
##$ -M russel.sutherland@kcl.ac.uk
##$ -m be
##$ -j y
#$ -V
#
#$ -o /home/rsutherlandbrc/networks/huppi2_noSL/mutsig_perm_results.std.out
##$ -t 1-6

##the start and end permutation numbers for each R script

#index=$(($SGE_TASK_ID - 1))

#the density of cliques I investiagted
run=0.75



#############################################################################################
############################################################################################


Rscript /home/rsutherlandbrc/networks/huppi2_noSL/DME_permutation_test_mutsig3.R  ${run}



##############################################################################################


