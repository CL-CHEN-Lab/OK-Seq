
#########################################################################################################################  Please make sure you already install Deeptools in your terminal and don't use the newest version of python environment! less then Python 3.6.4 should be fine ###########################################################################################

# PLOT average profile of RFD around TSS and TTS

computeMatrix scale-regions --regionsFileName {your bed file of interested regions/genes PATH e.g.codingGenes.bed} --beforeRegionStartLength {e.g. 10000} --afterRegionStartLength {e.g. 10000} --regionBodyLength {e.g. 20000} --binSize {e.g. 1000} --scoreFileName {RFD bigwig file PATH e.g. Hela.EdC.Combined_OkaSeq.RFD.bw} --outFileName {e.g. "OUTPUT.matrix"} --missingDataAsZero --skipZeros

plotProfile --matrixFile {e.g. "OUTPUT.matrix"} --outFileName {e.g. "RFD_averageProfile.stGeneLength.png"}  --averageType mean --startLabel {e.g. start/TSS} --endLabel {e.g. end/TTS} --plotType se
############################################################################################################

# PLOT Heatmap of RFD IZ center +/-100kb
computeMatrix reference-point --regionsFileName {your IZ bed file PATH e.g. HeLa_hmm_HMMsegments_IZ.bed} --beforeRegionStartLength {e.g. 100000} --afterRegionStartLength {e.g. 100000} --binSize {e.g. 1000} --scoreFileName {RFD bigwig file PATH e.g. Hela.EdC.Combined_OkaSeq.RFD.bw} --outFileName {e.g. "OUTPUT.matrix"} --missingDataAsZero --skipZeros --referencePoint center
plotHeatmap --matrixFile {e.g. "OUTPUT.matrix"} --outFileName {e.g. "RFD_sortbyLength.png"} --whatToShow "plot, heatmap and colorbar" --samplesLabel {e.g. HeLa} --refPointLabel center --sortUsing region_length --sortRegions ascend

### PLOT Heatmap of OEM IZ center +/-100kb
computeMatrix reference-point --regionsFileName {your IZ bed file PATH e.g. HeLa_hmm_HMMsegments_IZ.bed} --beforeRegionStartLength {e.g. 100000} --afterRegionStartLength {e.g. 100000} --binSize {e.g. 1000} --scoreFileName {series of OEM bigwig file PATH e.g. 20130819CGM130726.Hela_OEM_10kb.bw 20130819CGM130726.Hela_OEM_20kb.bw 20130819CGM130726.Hela_OEM_50kb.bw 20130819CGM130726.Hela_OEM_100kb.bw 20130819CGM130726.Hela_OEM_250kb.bw 20130819CGM130726.Hela_OEM_500kb.bw  20130819CGM130726.Hela_OEM_1Mb.bw} --outFileName {e.g. "OUTPUT.matrix"} --missingDataAsZero --skipZeros --referencePoint center
plotHeatmap --matrixFile {e.g. "OUTPUT.matrix"} --outFileName {e.g. "OEM_sortbyLength.png"} --whatToShow "plot, heatmap and colorbar" --refPointLabel center --samplesLabel {e.g. "HeLa 10kb" "HeLa 20kb" "HeLa 50kb" "HeLa 100kb" "HeLa 250kb" "HeLa 500kb" "HeLa 1Mb"} --sortUsing region_length --sortRegions ascend
