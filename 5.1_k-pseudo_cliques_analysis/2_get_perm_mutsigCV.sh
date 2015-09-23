###########################################
#!/bin/sh
#$ -S /bin/sh
#
##$ -cwd
#$ -N mutsigCV
##$ -M russel.sutherland@kcl.ac.uk
##$ -m be
#$ -j y
#$ -V
#
#$ -o /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/mutsig.out.txt
#$ -t 1-10

##The numeber of permutations to tun on each R script
permPerRun=1000

##the start and end permutation numbers for each R script

N2=$(($permPerRun * $SGE_TASK_ID))
N1=$(($N2 - 999))



Rscript /home/rsutherlandbrc/networks/huppi2_noSL/getPermutedMutSigCV_cluster.R  ${N1} ${N2} ${permPerRun}


done
