# OKseqHMM: a R package for profiling OK-Seq data to study the genome-wide DNA Replication Fork Directionality (RFD)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7428883.svg)](https://doi.org/10.5281/zenodo.7428883)


Package: OKseqHMM

Version: v2.0.0

Depends: R (>= 3.1.0)

deepTools (https://deeptools.readthedocs.io/en/develop/content/installation.html)

Python 3

Imports: HMM, Rsamtools, GenomicAlignments

License: GNU General Public License v3.0

Encoding: UTF-8

LazyData: true

RoxygenNote: 7.0.2

Original designer: Yves d'AUBENTON-CARAFA (CNRS, I2BC, C.Thermes's team)

Contributors and Maintainers: [Chun-Long CHEN](mailto:chunlong.chen@curie.fr) (Institut Curie), [Yaqun LIU](mailto:liu.yaqun19@gmail.com)

## Abstract: 

This R package is served for analyzing OK-Seq data from original mapping bam file to count matrices, RFD calculation, and inititation/termination zone calling. Firstly it transforms data into RFD (replication fork directionality) profiles for a primary visualisation (e.g. with IGV), and then by using the HMM package of R (http://www.r-project.org/), it can identify accurately most of the replication initiation zones (upward transitions on RFD profile), termination zones (downward transitions on RFD profile) and also the intermediate states (flat RDF profile).

## Description:

Despite intense investigation, human replication origins and termini remain elusive. Existing data have shown strong discrepancies. Here  N.Petryk et al sequenced highly purified Okazaki fragments from two cell types and, for the first time, quantitated replication fork directionality and delineated initiation and termination zones genome-wide.

Definition of replication fork directionality:

RFD was computed for each 1 kb window as follows:

RFD = (C - W) / (C + W)

where C and W correspond to the number of reads mapped on the Crick and the Watson strand, that reveal, respectively, the proportions of rightward- and leftward- moving forks within each window. Since the total amount of replication on both strands should be constant across the genome, we normalized the difference between two strands by the total read count to account for variation in read depth due to copy number, sequence bias and so on. RFD ranges from -1 (100% leftward moving forks) to +1 (100% rightward moving forks) and 0 means equal proportions of leftward- and rightward-moving forks. Replicate experiments produced RFD profiles that strongly correlated to each other (for HeLa cells, R=0.92 (Pearson), P<10^-15 (t-test) and for GM06990 cells, R=0.93, P<10^-15). Similar correlations were observed between RFD profiles with EdU or EdC labelling.

Segmentation of RFD profiles:

A four-state HMM was used to detect within the RFD profiles the AS, DS and FS segments representing regions of predominant initiation (‘Up’ state), predominant termination (‘Down’ state) and constant RFD (‘Flat1’ and ‘Flat2’ states). 


<img align="left" src="https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig4.png" hspace="30" width="320" height="300"/>
<img align="left" src="https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig1.png" hspace="30" width="320" height="300"/>
<br/><br/><br/><br/><br/>

!["Fig.3 emission (left table) and transition (right table) probabilities used in the HMM segmentation." ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig2.png) 


In segmentation process, the RFD values were computed within 15 kb sliding windows (stepped by 1 kb across the autosomes). The HMM used the Delta RFD values between adjacent windows (that is, Delta RFDn = (RFD(n+1)  - RFD(n))/2 for window n). Windows with <30 reads on both strands were masked. The Delta RFD values were divided into five quantiles and the HMM package of R (http://www.r-project.org/) was used to perform the HMM prediction with probabilities of transition and emission indicated as below. 

Only segments reproducibly identified in both biological replicates were retained. The choice of a 15 kb sliding window is based on a compromise between spatial resolution and reproducibility of AS detection among biological replicates. The efficiency of AS was estimated as:

Delta RFD (segment) = (RFD(end) - RFD(start)) / 2

where RFD(start) and RFD(end) correspond, respectively, to the RFD values computed in 5 kb windows around left and right extremities of each segment.

The final RFD profile and the relication origin/ termination zones calling are like this: Red (blue) lines above (below) HeLa RFD profile indicate ascending (descending) HMM-detected segments, corresponding to initiation (termination) zones; magenta and cyan arrows indicate genes.


![    Fig.1 Red (blue) lines above (below) HeLa RFD profile indicate ascending (descending) HMM-detected segments (see Methods section); magenta and cyan arrows indicate genes. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig3.png) 

## R package OKseqHMM function usage: 

Please make sure that you already install and import the associated R package HMM, Rsamtools and GenomicAlignments in R or Rstudio, and also install samtools in console or terminal before executing this function.

### For the input:

Please preprare the aligned OK-Seq data, which is a Bam file with the corresponding indexed file (.bai) and the annotation coordinates which containing all chromosomes and their lengths (can be easily downloaded from UCSC, e.g. hg19.chr.sizes.txt).

Then just enter the file path linking to your bam file and the annotation file, read coverage threshold(e.g. 10), the smoothing window size (e.g. 15 as kilobase) and then also define the prefix of the output files:

(e.g. OKseqHMM(bamfile = "PATH/my.bam",thresh = 10, chrsizes = "hg19.chr.size.txt", binSize= 1000, winS=15, fileOut = "PATH/my_hmm"))

The program can automatically identify that the input bam file is paired-ends or single-end data, then seperate the bam into W and C strands, respectively, to calculate the 1kb binsize coverage. 

This R function takes the threshold as user set to remove the abnormal counts for each strand and then smooths the data into 15 kb sliding windows (with a step of 1 kb) to get the best RFD profile and the corresponding HMM calling results.

### For the output:

You will obtain a serie of output files, which are:

1-4, two bam files for the forward and reverse strand seperated from the input bam and their corresponding index files.

5-6, two RFD files ("_RFD.bedgraph", one is 1kb binsize and the smoothing one) in bedgraph format that allow you to visualize directly the replication fork direction profile in some integrative genomic viewers or browsers like IGV/IGB. You can also transform the bedgraphs into bigWig by the UCSC tool bedGraphToBigWig (http://hgdownload.soe.ucsc.edu/admin/exe/) to get binary compressed files.

7,log file ("_log.txt") that records all of the parameters you use and also the default setting information.

8,HMM result text file ("_HMM.txt") that records all of the global optimal hidden states calculated by HMM Viterbi algorithm.

9,HMM result text file ("_HMMpropa.txt") that records all of the previous state positions that caused the maximum local probability of a state by HMM posterior algorithm.

10-17, generate 8 result text files that records the genomic positions and also the other corresponding probabilities for the final identified optimal states:

"_HMMsegments_IZ.txt/bed" is for the replication initiation zone calling result.

"_HMMsegments_TZ.txt/bed" is for the replication termination zone calling result.

"_HMMsegments_highFlatZone.txt/bed" and "_HMMsegments_LowFlatZone.txt/bed" are for the two replication intermediate (flat) zones calling results.

### Test data and published results
You can find a test dataset and the corresponding results in the "templates_bam2hmm" folder.

And the OK-Seq data and results from our previous publications, as well as our re-analysis results of the OK-seq data published by others, are available at the "published_results" folder.

## R package OKseqOEM function usage: 

For a further investigation of origin expanding efficiency (i.e. deltaRFD), we provide here a second function to visualize it directly in multiple scales.

### For the input:

Please preprare the two bam files of Forward and Reverse strand generated by previous OKseqHMM and again the annotation coordinates which containing all chromosomes and their lengths. Entering the path linking to your bam files, annotation file, the prefix of the output file and define a series of window sizes as different visualisation scales that you would like to check in the genomic viewer. 

(e.g. for human cells: 

OKseqOEM(bamInF ="path_to_bam_Forward_strand", bamInR="path_to_bam_Reverse_strand", chrsizes = "hg19.chr.size.txt", fileOut = "path/name_of_my_OEM", binsize=1000, binList=c(1,10,20,50,100,250,500,1000)))

Then this function can calculate directly in different scales, the origin expanding tracks as you set.

### For the output:

You will obtain a serie of wiggle files calculated by using different sliding window sizes (always with 1kb step), which allow you to load them directly into IGV for visualisation.  To improve the speed of loading and decrease the memory used, you can also transform the wiggle into bigWig by the UCSC tool wigToBigWig (http://hgdownload.soe.ucsc.edu/admin/exe/) to get binary compressed files.


![    Fig.5 OEM profiles for HeLa. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig5.png) 

In the genomic viewer, you can set different colors (here, positive values in pink and negative values in blue) to see the RFD transitions and the pink enrichment bulks are perfectly matched with our peak calling results. It should be noted that the replication initiation zones show a highest positive transition at the small scales (most at 10kb-100kb) with a faster expanding fashion, while the termination zones expand larger (100kb-1MB) regions in a much slower manner of delta RFD transition.

## Module AveragePlot for the metaplot of average profile/heatmap: 

Please make sure that you already installed the deepTools (https://deeptools.readthedocs.io/en/develop/index.html) and the python environment (the recommanded version is python 3.6.4. The latest python version could cause some incompatible issue with deepTools.) The shell-based script "average_profile_heatmap.sh" show us the template about how to use the deepTools to generate the average profile and heatmap among the annotated genes or the regions of interest by using the computeMatrix and plotProfile/plotHeatmap function, defining the two upstream and downstream border size and intragenic body size and also the other parameters indicated in the script.

## Reference:

Petryk N., Kahli M., d'Aubenton-Carafa Y., Jaszczsyzyn Y., Shen Y., Sylvain M., Thermes C., CHEN C.L.#, and Hyrien O.# (#co-last authors). Replication landscape of the human genome. ***Nat. Commun.*** 7:10208 (2016). https://doi.org/10.1038/ncomms10208

Promonet A.\*, Padioleau I.\*, LIU Y.\* (\*co-first authors), Sanz L., Schmitz A., Skrzypczak M., Sarrazin A., Ginalski K., Chedin F., Rowicka  M., CHEN C.L.#, Lin Y.L.# and Pasero P.# (#co-last authors). Topoisomerase 1 prevents R-loop mediated replication stress at transcription termination sites. ***Nat. Commun.*** 11:3940 (2020) https://doi.org/10.1038/s41467-020-17858-2

Liu, Y., Wu, X., D'Aubenton-Carafa, Y., Thermes, C. & Chen, C.-L. OKseqHMM: a genome-wide replication fork directionality analysis toolkit. ***Nucleic Acids Research*** (2023) https://doi.org/10.1093/nar/gkac1239

Wu, X., Liu, Y., d’Aubenton-Carafa, Y. et al. Genome-wide measurement of DNA replication fork directionality and quantification of DNA replication initiation and termination with Okazaki fragment sequencing. ***Nat. Protoc.*** (2023). https://doi.org/10.1038/s41596-022-00793-5
