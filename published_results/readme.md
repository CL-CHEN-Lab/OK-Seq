The raw OK-Seq data of HeLa and GM06990 cells published in (Petryk, et al. Replication landscape of the human genome. Nature communications 7 (2016): 10208.) have been deposited at NCBI Sequence Read Archive with number SRP065949: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP065949.

The profiles of results visualized in IGV is like below:

![    Fig.1 1kb binsize RFD profiles for HeLa and GM06990 and the correspondaing initiation zones (IZs) by OKseqHMM. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/igv_snapshot_okseq_Hela_GM.png) 


And also the corresponding origin expanding effficiency profiles calculated by OKseqOEM (HeLa):

![    Fig.2 1kb binsize OEM profile by OKseqOEM for HeLa. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/fig5.png) 

For more OK-seq of different cell lines (i.e. BL79, IARC385, K562, IB118, IMR90, Raji, TF1, TLSE cells), the raw data source is available in (Wu X, et al. Developmental and cancer-associated plasticity of DNA replication preferentially targets GC-poor, lowly expressed and late-replicating regions. Nucleic Acids Res. 2018.), which has been deposited at ENA database with accession number PRJEB25180: https://www.ebi.ac.uk/ena/browser/view/PRJEB25180.

![    Fig.3 1kb binsize RFD profiles for different cell lines and the correspondaing initiation zones (IZs) by OKseqHMM. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/igv_snapshot_okseq_Xia_diff_cell.png) 

OEM profile by OKseqOEM for IARC385:

![    Fig.4 1kb binsize OEM profile by OKseqOEM for IARC385. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/igv_snapshot_OEM_TF1_GFP.png) 

OEM profile by OKseqOEM for TF1-GFP:

![    Fig.5 1kb binsize OEM profile by OKseqOEM for TF1-GFP. ](https://github.com/CL-CHEN-Lab/OK-Seq/blob/master/img/igv_snapshot_IARC385_OEM_HMM.png) 

Additional processing RFD data of HeLa and GM06990 cells (e.g. read counts on W and C strands, biological replicates, etc., in bigwig format, 1 kb resolution, human genome version hg19) are available at: http://xfer.curie.fr/get/SbyS2ueTaTz/OKSeq.Data.zip.

You can also downlaod at the current folder the initiaton zones, termination zones and highRFD regions of corresponding cell types optained by the HMM process.

If the link is not accessible, please contact the maintainer for the update.
