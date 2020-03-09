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

A four-state HMM was used to detect within the RFD profiles the AS, DS and FS segments representing regions of predominant initiation (‘Up’ state), predominant termination (‘Down’ state) and constant RFD (‘Flat1’ and ‘Flat2’ states; see Fig.1 for a graphic representation). 

![Fig.1 state model used in the HMM segmentation (Methods); Up, regions of predominant initiation, i.e. AS; Down, regions of predominant termination, i.e. DS; Flat1 and Flat2, FS. ](/Users/yliu2/Documents/20171010_Cnflict_Repli_Trans_Analysis_CellCycle_Yaqun/script/hmmPolarity/fig1.png) 

The RFD values were computed within 15 kb sliding windows (stepped by 1 kb across the autosomes). The HMM used the Delta RFD values between adjacent windows (that is, Delta RFDn = (RFD(n+1)  - RFD(n))/2 for window n). Windows with <30 reads on one strand were masked. The Delta RFD values were divided into five quantiles and the HMM package of R (http://www.r-project.org/) was used to perform the HMM prediction with probabilities of transition and emission indicated in the Fig. 2. Only segments reproducibly identified in both biological replicates were retained. The choice of a 15 kb sliding window is based on a compromise between spatial resolution and reproducibility of AS detection among biological replicates (Fig. 3). The efficiency of AS was estimated as:

Delta RFD (segment) = (RFD(end) - RFD(start)) / 2

where RFD(start) and RFD(end) correspond, respectively, to the RFD values computed in 5 kb windows around left and right extremities of each segment.

## Reference:

Petryk, Nataliya, et al. "Replication landscape of the human genome." Nature communications 7 (2016): 10208.
