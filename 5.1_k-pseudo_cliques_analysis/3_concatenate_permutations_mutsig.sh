##############################################
#!/bin/sh
#$ -S /bin/sh
#
##$ -cwd
#$ -N mutsig_concatenate_3
##$ -M russel.sutherland@kcl.ac.uk
##$ -m be
#$ -j y
#$ -V
#
#$ -o /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/concatenate.std.out.txt
#$ -t 1-6

##the start and end permutation numbers for each R script

index=$(($SGE_TASK_ID - 1))

#the density of cliques I investiagted
runs=(0.95 0.90 0.85 0.8 0.75 0.70)

#the number of k-cliques I investigated
k=(3 4 5 6 7)


#############################################################################################
############################################################################################
# the overall results

mkdir /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated

##copy the files to the concanate folder
find  /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]} -name  "permutedMutSigCV_huppi2_noSL*" -type f -exec cp {} /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated \;

#concatenating the sample diversity files together in to one file of 10000 lines. one line per permutation
cat /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/*.txt > /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/permutedMutSigCV_huppi2_noSL_a${runs[$index]}_10000

#removing the individual permutation files after concatenating them toegther 
rm -f /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/*txt

#adding a .txt extension to the concatenated samplediversity 10000 permutation file
#mv /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/permutedMutSigCV_huppi2_noSL_10000 /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/permutedMutSigCV_huppi2_noSL_a${runs[$index]}_10000.txt

#################################################################################################
#################################################################################################
### the code for each of the k-pseduocliques



####for num in $($k); do

for i in 3 4 5 6 7
do

echo "$i"

##copy the files to the concanate folder
find  /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]} -name  "permutedMutSigCV_huppi2noSL_M0_SD5_a${runs[$index]}_pseudo_${i}"* -type f -exec cp {} /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated \;

##concatenating the sample diversity files together in to one file of 10000 lines. one line per permutation
cat /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/permutedMutSigCV_huppi2noSL_M0_SD5_a${runs[$index]}_pseudo_${i}* > /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/permutedMutSigCV_huppi2noSL_a${runs[$index]}_pseudo_${i}_10000

##removing the individual permutation files after concatenating them toegther
rm -f /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/*M0_SD5*

##adding a .txt extension to the concatenated samplediversity 10000 permutation file
#mv /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSLa${runs[$index]}/concatenated/permutedMutSigCV_huppi2noSL_a${runs[$index]}_pseudo_${i}_10000 /home/rsutherlandbrc/networks/huppi2_noSL/mutsig/huppi2_noSL${runs[$index]}/concatenated/permutedMutSigCV_huppi2noSL_a${runs[$index]}_pseudo_${i}_10000.txt

done
#######################################################################################################



