###########################################
#!/bin/sh
#$ -S /bin/sh
#
##$ -cwd
#$ -N huppi_network_permutation
##$ -M russel.sutherland@kcl.ac.uk
##$ -m be
#$ -j y
#$ -V
#
#$ -o /home/rsutherlandbrc/networks/std.out.txt
#$ -t 1-10

##The numeber of permutations to tun on each R script
permPerRun=1000

##the start and end permutation numbers for each R script

N2=$(($permPerRun * $SGE_TASK_ID))
N1=$(($N2 - 999))


Rscript /home/rsutherlandbrc/networks/network_permutation_cluster_script.R  ${N1} ${N2}
