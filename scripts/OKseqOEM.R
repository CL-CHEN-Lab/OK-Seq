#force R to use regular numbers instead of using the e+10- (exponential)like notation
options(scipen = 999)

OKseqOEM <- function(bamInF, bamInR, chrsizes, fileOut, binSize, binList)
{

  chromNames <-  read.table(chrsizes,header=FALSE,sep="\t",comment.char = "#",stringsAsFactors = FALSE)
  chr.sizes <- data.frame(chr=chromNames[,1],size=chromNames[,2])
  require(Rsamtools)
  paired <- testPairedEndBam(bamInF)
  for (i in c(1:nrow(chr.sizes))){
    chr.name <- chr.sizes$chr[i]
    print(chr.name)
    chr.length <- chr.sizes$size[i]
    print(chr.length)
    if (paired)
    {
      print(paste0("It's pair-end. Calculating ",binSize,"bp binsize coverage for forward strand."))
      system(paste0("samtools view -q 1 -f 0x42 -F 0x4 ",bamInF," ",chr.name," > fwd_",chr.name,".sam"))
      system(paste0("awk '$3~/^", chr.name, "$/ {print $2 \"\t\" $4}' fwd_",chr.name,".sam > fwd_",chr.name,".txt"))

      print(paste0("Calculating ",binSize,"bp binsize coverage for reverse strand."))
      system(paste0("samtools view -q 1 -f 0x42 -F 0x4 ",bamInR," ",chr.name," > rev_",chr.name,".sam"))
      system(paste0("awk '$3~/^", chr.name, "$/ {print $2 \"\t\" $4}' rev_",chr.name,".sam > rev_",chr.name,".txt"))
    }else{
      print(paste0("It's single-end. Calculating ",binSize,"bp binsize coverage for forward strand."))
      system(paste0("samtools view ",bamInF," ",chr.name," > fwd_",chr.name,".sam"))
      system(paste0("awk '$3~/^", chr.name, "$/ {print $2 \"\t\" $4}' fwd_",chr.name,".sam > fwd_",chr.name,".txt"))

      print(paste0("Calculating ",binSize,"bp binsize coverage for reverse strand."))
      system(paste0("samtools view ",bamInR," ",chr.name," > rev_",chr.name,".sam"))
      system(paste0("awk '$3~/^", chr.name, "$/ {print $2 \"\t\" $4}' rev_",chr.name,".sam > rev_",chr.name,".txt"))
    }

    fileInF <- paste0("fwd_",chr.name,".txt")
    tmp<- read.table(fileInF, header=F, comment.char="",colClasses=c("integer","integer"),fill=TRUE)
    tags <- tmp[,2]
    tags[tags<3] <- 0
    breaks <- seq(0, chr.length+binSize, by=binSize)
    h <- hist(tags, breaks=breaks, plot=FALSE)
    Temp.chr.F <- h$counts

    fileInR <- paste0("rev_",chr.name,".txt")
    tmp <- read.table(fileInR, header=F, comment.char="",colClasses=c("integer","integer"),fill=TRUE)
    tags <- tmp[,2]
    tags[tags<3] <- 0
    breaks <- seq(0, chr.length+binSize, by=binSize)
    h <- hist(tags, breaks=breaks, plot=FALSE)
    Temp.chr.R <- h$counts

    # delete intermediate files. If keep them, mask these 4 lines with #
    system(paste0("rm fwd_",chr.name,".sam"))
    system(paste0("rm rev_",chr.name,".sam"))
    system(paste0("rm ",fileInF))
    system(paste0("rm ",fileInR))

    #get the cumulative sum
    Temp.chr.F <- cumsum(Temp.chr.F)
    Temp.chr.R <- cumsum(Temp.chr.R)

    print("Calculating OEM.")
    for (n in c(1:length(binList)))
    {

      print(paste0("The smoothing window size for OEM is ",binList[n]*binSize/1000,"kb."))

      Data.chr.F <- Temp.chr.F[(binList[n]+1):length(Temp.chr.F)]-Temp.chr.F[1:(length(Temp.chr.F)-binList[n])]
      Data.chr.R <- Temp.chr.R[(binList[n]+1):length(Temp.chr.R)]-Temp.chr.R[1:(length(Temp.chr.R)-binList[n])]

      Data.chr.Smooth <- Data.chr.F/(Data.chr.F+Data.chr.R)
      Data.chr <- Data.chr.Smooth[(binList[n]+1):length(Data.chr.Smooth)]-Data.chr.Smooth[1:(length(Data.chr.Smooth)-binList[n])]

      Data.chr <- c(rep(NA,binList[n]-1),Data.chr)
      Data.chr[which(is.na(Data.chr))] <- 0
      # Data.chr<- Data.chr[1:chr.length]

      ##Save file in wig format
      if (i==1) {
        Title <- paste0("fixedStep chrom=", chr.name, " start=1 step=",binSize," span=",binSize, sep="")
        fileOutWig <- paste0(fileOut,"_OEM_",binList[n]*binSize/1000,"kb.wig")
        write.table(Data.chr, file=fileOutWig, quote = FALSE, row.names = FALSE, col.names=Title, append = FALSE)
      } else {

        Title <- paste0("fixedStep chrom=", chr.name, " start=1 step=",binSize," span=",binSize, sep="")
        fileOutWig <- paste0(fileOut,"_OEM_",binList[n]*binSize/1000,"kb.wig")
        write.table(Data.chr, file=fileOutWig, quote = FALSE, row.names = FALSE, col.names=Title, append = TRUE)

      }
    }
  }
}



















