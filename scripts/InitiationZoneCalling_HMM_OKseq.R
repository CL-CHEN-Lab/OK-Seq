####  Replication initiation zone calling
####  OK-Seq
####  2020-01-21
####  lab Chen
#################################################################################################################
#  
# for 15kb bin of HeLa cell line
# set the quantile as 
quant <- c(-1, -0.0082058939609862, -0.00141890249101162, 0.00103088286465956, 0.00800467305420799, 1)   


#     Initialize HMM 4 states, observations, start probability, emission probability, transition probability 
#============================================================================================================

hmmPolarity <- function(fileW, fileC, fileOut, binSize=1000, chrom, thresh=30, winS=15,hwinS=winS/2,
                           st=c("D", "L", "H", "U"),
                           sym=c("V", "W", "X", "Y", "Z"),
                           pstart=rep(1/4, 4),
                           pem=t(matrix(c(0.383886256, 0.255924171, 0.170616114, 0.113744076, 0.075829384,
                                          .10,.20,.40,.20,.10, 
                                          .10,.20,.40,.20,.10,
                                          0.022222222, 0.033333333, 0.066666667, 0.211111111, 0.666666667), 
                                        ncol=4)),
                           ptrans=t(matrix(c(0.9999,0.000020,0,0.000080, 
                                             0,0.999,0,0.001, 
                                             0.001,0,0.999,0, 
                                             0.000080,0,0.000020,0.9999), 
                                           ncol=4)),
                           quant=NA)

{
  
  require(HMM)
  
  # 1kb bs of Watson strand		
  ta <- read.table(file=fileW, sep="\t", header=F, as.is=T)
  w <- ta[,1]
  
  # 1kb bs of Crick strand		
  ta <- read.table(file=fileC, sep="\t", header=F, as.is=T)
  c <- ta[,1]
  
  # raw polarity for later
  polar <- c/(c+w)
  polar[c<thresh & w<thresh] <- NA
  
  # smoothing into 15kb binsize
  print(paste("window size :", winS))
  sw <- cumsum(w)
  lg <- length(w)
  from <- (-hwinS+2):(lg-hwinS+1)
  to <- from+winS-1
  from[from<1] <- 1
  to[to>lg] <- lg
  
  print("")
  print(paste("number of bins :", length(v)))
  print("")
  
  win <- matrix(c(from, to), ncol=2)
  ws <- apply(win ,1, function(x) { (sw[x[2]]-sw[x[1]])/winS } )
  
  sc <- cumsum(c)
  lg <- length(c)
  cs <- apply(win ,1, function(x) { (sc[x[2]]-sc[x[1]])/winS } )
  
  thresh = thresh/winS
  print(paste("cutoff is :",thresh))
  bias <- cs/(ws+cs) 
  bias[is.na(bias)] <- 0.5
  bias[ws<thresh & cs<thresh] <- 0.5
  
  write.table(bias, file = paste(fileOut, "_bias.txt", sep=""), append = F, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  
  # HMM from new deltas =============================
  
  # derive
  delta <- c(0,bias[-1]-bias[-length(bias)])
  delta[is.na(delta)] <- 0.5
  
  # affect symbols
  if (is.na(quant[1])) { quant <- quantile(delta, probs = seq(0, 1, 0.20)) }
  quant[1] <- -1
  quant[length(quant)] <- 1
  
  print("quantile borders :")
  print(quant)
  print("")
  dx  <- unlist(sapply(delta, function(x) { ix <- which(x>=quant); ix[length(ix)] })) 
  dx[dx>5] <- 5
  
  # write log ==================================
  
  logFile <- paste(fileOut, "_log.txt", sep="")
  write.table(data.frame(c("fileW",fileW)), file = logFile, 
              append = F, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  write.table(data.frame(c("fileC",fileC)), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  write.table(data.frame(c("fileOut",fileOut)), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  write.table(data.frame("ptrans"), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(ptrans, file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  write.table(data.frame("pem"), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(pem, file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(data.frame("pstart"), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(t(data.frame(pstart)), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(data.frame("st"), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(t(data.frame(st)), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  write.table(data.frame("sym"), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(t(data.frame(sym)), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  write.table(data.frame("quant"), file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)			
  write.table(quant, file = logFile, 
              append = T, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  
  # HMM ==============================
  
  hmm1 <- initHMM(States=st, Symbols=sym, startProbs=pstart, transProbs=ptrans, emissionProbs=pem)
  print(hmm1)	
  print("wait, viterbi...")	
  seg <- viterbi(hmm=hmm1, observation=dx)
  
  # pour affichage du profil des etats
  prof <- rep(NA, times=length(seg))
  prof[seg=="U"] <- 1
  prof[seg=="H"] <- 0.2
  prof[seg=="L"] <- -0.2
  prof[seg=="D"] <- -1
  prof[is.na(v[-length(v)])] <- -1.1
  prof[is.na(w[-length(w)])] <- -1.1
  prof[is.na(prof)] <- -.1
  
  write.table(prof, file = paste(fileOut, "_HMM.txt", sep=""), append = F, quote = FALSE, sep = "\t", col.names=F, row.names=F)
  
  # adding the probability curve =========================
  print("wait, probabilities…")
  post <- posterior(hmm1,dx)
  
  prb <- rep(NA, length(seg))
  for (i in 1:length(seg))
  {
    prb[i] <- post[seg[i],i]
  }
  prb[prb>1] <- 1	
  
  write.table(as.integer(prb*1000), file = paste(fileOut, "_HMMproba.txt", sep=""), append = F, 
              quote = FALSE, sep = "\t", col.names=F, row.names=F)
  
  
  # calculate the region coordinates by profil viterbi ========================
  
  left <- seg[-length(seg)]
  right <- seg[-1]
  ix <- which(left!=right)
  from <- c(1,ix+1)
  to <- c(ix,length(seg))
  
  # for the table, adding states and calculating lengths, slopes and associated probabilities ====================
  # calculate the slope of the polarity of the states (between the left and right part of the state))
  
  states <- meanPol <- ymax <- ymin <- cpp <- p <- inc <- napc <- corr <- rep(NA, length(from))
  states <- seg[from]
  
  for (i in 1:length(from))
  {
    pos <- from[i]:to[i]
    lgpos <- length(pos)
    
    m <- prb[pos[!is.na(prb[pos])]]
    if (length(m)>0)
    {
      p[i] <- mean(m, na.rm=T)
      cpp[i] <- 1-(sum((1-m), na.rm=T)/((1-min(m, na.rm=T))*length(m)))
    }
    meanPol[i] <- mean(polar[pos], na.rm=T)
    
    pos2 <- pos[!is.na(polar[pos])]
    lgpos2 <- length(pos2)
    napc[i] <- 100-round(lgpos2/lgpos*100)
    if (lgpos2<2) 
    { 
      res <- data.frame(NA, NA)
      inc[i] <- ymin[i] <- ymax[i] <- NA 
    } else { 
      realPos <- pos2*binSize
      res <- lm(polar[pos2] ~ realPos)
      inc[i] <- res[[1]][2]
      ymin[i] <- 	(inc[i]*realPos[1])+res[[1]][1]
      ymax[i] <- 	(inc[i]*realPos[lgpos2])+res[[1]][1]
      corr[i] <- cor(x=realPos, y=polar[pos2])
    }
  }
  
  ymin[ymin<0] <- 0
  ymax[ymax<0] <- 0
  ymin[ymin>1] <- 1
  ymax[ymax>1] <- 1
  
  # for display ==================================
  
  inc1 <- round(inc*10^8)		# as a percentage of RFD per megabase
  cpp <- round(100*cpp)
  p <- round(p*100)
  chr <- rep(chrom,length(from))
  from1 <- (from-1)*binSize+1
  to1 <- to*binSize
  lg1 <- to1-from1+1
  meanPol1 <- round(meanPol*100)
  ymin1 <- round(ymin*100)
  ymax1 <- round(ymax*100)
  corr1 <- round(corr*100)
  
  # adjustment of the border polarities by the average and recalculation of the slope ======================
  
  polL <- ymin1
  polR <- ymax1
  
  polD <- polR[-length(polR)]
  polG <- polL[-1]
  
  polD[is.na(polD)] <- polG[is.na(polD)]
  polG[is.na(polG)] <- polD[is.na(polG)]
  
  polM <- round((polG+polD)/2)
  polL[-1] <- polR[-length(polR)] <- polM
  
  # adjusted slope: 10 ^ 6 for clarity of the table display (% by megabase)
  slope_adj <- round(10^6*(polR-polL)/(to1-from1))
  
  
  # writing =============================
  dataOut <- data.frame(chr, from=from1, to=to1, state=states, length=lg1, slope=inc1, 
                        p, fcp=cpp, pol_mean=meanPol1, pol_left=ymin1, pol_right=ymax1, 
                        na=napc, cor=corr1, slope_adj=slope_adj, pol_adj_left=polL, pol_adj_right=polR)
  write.table(dataOut, file = paste(fileOut, "_HMMsegments.txt", sep=""), append = F, 
              quote = FALSE, sep = "\t", col.names=T, row.names=F)
  
}

# end of the function =========================================================================
