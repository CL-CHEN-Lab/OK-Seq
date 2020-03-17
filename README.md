# hmmPolarity: R package for profiling OK-Seq data to study the genome-wide DNA replication fork directionality

Package: hmmPolarity

Version: 0.1.0

Depends: R (>= 3.1.0)

Imports: HMM

License: GNU General Public License v3.0

Encoding: UTF-8

LazyData: true

RoxygenNote: 7.0.2

Original Authors: Yves d'AUBENTON-CARAFA <yves@daubenton.net>, Chun-long CHEN <chunlong.chen@curie.fr>

Second author and Maintainer: Yaqun LIU <yaqun.liu@curie.fr>

## Abstract: 

This R package is served for analyzing OK-Seq data from count matrices, RFD calculation to inititation/termination zone calling. Firstly it could transform data into RFD (replication fork directionality) profiles for a primary visualisation and then by using the HMM package of R (http://www.r-project.org/), it can identify accurately most of the replication initiation zones, termination zones and also the intermediate states.

## Description:

Despite intense investigation, human replication origins and termini remain elusive. Existing data have shown strong discrepancies. Here  N.Petryk et al sequenced highly purified Okazaki fragments from two cell types and, for the first time, quantitated replication fork directionality and delineated initiation and termination zones genome-wide.

Definition of replication fork directionality:

RFD was computed for each 1 kb window as follows:

RFD = (C - W) / (C + W)

where C and W correspond to the number of reads mapped on Crick and Watson strand, that reveal, respectively, the proportions of rightward- and leftward- moving forks within each window. Since the total amount of replication on both strands should be constant across the genome, we normalized the difference between two strands by the total read count to account for variation in read depth due to copy number, sequence bias and so on. RFD ranges from -1 (100% leftward moving forks) to +1 (100% rightward moving forks) and 0 means equal proportions of leftward- and rightward-moving forks. Replicate experiments produced RFD that strongly correlated to each other (for HeLa cells, R=0.92 (Pearson), P<10^15 (t-test) and for GM06990 cells, R=0.93, P<10^15). Similar correlations were observed between EdU and EdC profiles

Segmentation of RFD profiles:

A four-state HMM was used to detect within the RFD profiles the AS, DS and FS segments representing regions of predominant initiation (‘Up’ state), predominant termination (‘Down’ state) and constant RFD (‘Flat1’ and ‘Flat2’ states). 


<img align="left" src="https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig4.png" hspace="30" width="350" height="300"/>
<img align="left" src="https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig1.png" hspace="30" width="350" height="300"/>
<br/><br/><br/><br/><br/>


The RFD values were computed within 15 kb sliding windows (stepped by 1 kb across the autosomes). The HMM used the Delta RFD values between adjacent windows (that is, Delta RFDn = (RFD(n+1)  - RFD(n))/2 for window n). Windows with <30 reads on one strand were masked. The Delta RFD values were divided into five quantiles and the HMM package of R (http://www.r-project.org/) was used to perform the HMM prediction with probabilities of transition and emission indicated as below. 


!["Fig.3 emission (left table) and transition (right table) probabilities used in the HMM segmentation." ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig2.png) 


Only segments reproducibly identified in both biological replicates were retained. The choice of a 15 kb sliding window is based on a compromise between spatial resolution and reproducibility of AS detection among biological replicates. The efficiency of AS was estimated as:

Delta RFD (segment) = (RFD(end) - RFD(start)) / 2

where RFD(start) and RFD(end) correspond, respectively, to the RFD values computed in 5 kb windows around left and right extremities of each segment.

The final RFD profile and the relication origin/ termination zones calling are like this: Red (blue) lines above (below) HeLa RFD profile indicate ascending (descending) HMM-detected segments, corresponding to initiation (termination) zones; magenta and cyan arrows indicate genes.


![    Fig.1 Red (blue) lines above (below) HeLa RFD profile indicate ascending (descending) HMM-detected segments (see Methods section); magenta and cyan arrows indicate genes. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig3.png) 

## R package hmmPolarity function usage: 

Please make sure that you already install and import the associated R package HMM before.

### For the input:

Please preprare two .bedgraph files (The format is by defaut as 4 columns : chr, start, end and counts) for each strand (watson and crick) of OK-Seq data. And OK-seq data needs to be calculated into 1kb adjacent binsize, which can be easily got by using some bio-informatic tools like deeptools.

Then just easily enter the pathways linking to your 2 bedgraph files and also define the prefix of the output files:

(e.g. hmmPolarity(fileW= ".../my_watson_1kb_okseq.bedgraph", fileC=".../my_crick_1kb_okseq.bedgraph", fileOut=".../my_hmm_OK_result"))

By default, this R function takes the threshold as 30 to remove the abnormal counts for each strand and then smooths the data into 15kb windoz size to get the best RFD profile and the corresponding HMM calling results.

### For the output:

You will obtain 8 output files which are:

1, RFD file ("_RFD.wig") in wiggle format that allows you to visualize directly the replication fork direction profile in some integrative genomic viewers or browsers like IGV/IGB.

2,log file ("_log.txt") that records all of the parameters you use and also the default setting information.

3, text file ("_HMM.txt") that records all of the global optimal hidden states calculated by HMM Viterbi algorithm.

4,text file ("_HMMpropa.txt") that records all of the previous state positions that caused the maximum local probability of a state by HMM posterior algorithm.

5-8, generate 4 text files that records the genomic positions and also the other corresponding possibilities for the final identified optimal states:

"_HMMsegments_IZ.txt" is for the replication initiation zone calling result.

"_HMMsegments_TZ.txt" is for the replication termination zone calling result.

"_HMMsegments_highFlatZone.txt" and "_HMMsegments_LowFlatZone.txt" are for the two replication intermediate zones calling results.

## Reference:

Petryk, Nataliya, et al. "Replication landscape of the human genome." Nature communications 7 (2016): 10208.
