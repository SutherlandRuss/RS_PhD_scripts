# RS_PhD_scripts
R scripts used thoughout my PhD split according to chapter
The .R scripts should all work correctly after running "https://github.com/SutherlandRuss/RS_PhD_scripts/my_functions.R"

5.1 contains the scripts used to run k-psuedocliques analysis. 
    This is a method that uses a protein interaction network to prioritise gene level score data. 
    The steps to analysis are as follows:
      1. Obtain the dense module enumeration:http://people.kyb.tuebingen.mpg.de/georgii/dme.html.
      2. Use DME to analyse a biological network and extract the communities that satisfy various density thresholds ranging from 0.95 to 0.75.
      3. Conduct a gene-level test. (I used the GenePattern server [http://www.broadinstitute.org/cancer/software/genepattern/] to analyse colorectal cancer data from the TCGA PanCancer dataset])
      3. Using the resultant file and the biological network (with gene names), run 1_net_perm.sh on sung grid engine. This will call the network_permutation_cluster_script.R and generate 10 000 permuted networks.
      4. Run the 2_get_perm_mutsigCV.sh script. This will access the getPermutedMutSigCV_cluster.R script. This will obtain a median gene-level score for each community discovered by the DME algorith for each of the 10 000 permuted networks.
      5. Run 4_permutation_Test_Results_single_mutsig.sh. This will compare the median gene-level test statistic for each community in the observed network, to that of the 10 000 permutations.
        A community is considered enriched with higih scoring genes if the median gene-level test statistic is within the top 5% of the permuted scores once adjusted for multiple testing across communities.
        A statistical report is generated and the significant communities are output in a network plot where the genes are colored from white through to green through to red indicating the gene level scores.
        
        
NB.
The .sh scripts are to be used with sun grid engine and may require tinkering.
The .R scripts should all work correctly after running "https://github.com/SutherlandRuss/RS_PhD_scripts/my_functions.R", but will require edits to file paths.
