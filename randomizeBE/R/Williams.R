# Construct sequences of a williams design
# 
# Author: dlabes
###############################################################################

# Sheehe PR, Bross IDJ (1961)
# "Latin Squares to Balance Immediate Residual and Other Effects." 
# Biometrics, 17, 405-414.
#
# Jones B, Kenward MG (2003). 
# Design and Analysis of Cross-Over Trials. 2nd edition.
# Chapman & Hall, London.

williams <- function(ntmt=4, tmts=NULL)
{ # some error checking 
  if (is.numeric(ntmt) == FALSE || 
      any(length(ntmt) != 1, ntmt%%1 != 0, ntmt < 2) == TRUE) {
    stop("Number of treatments is not an integer >1.")
  }
  if (!is.null(tmts) & length(tmts)!=ntmt){
    stop ("Treatment codes must have ", ntmt, " entries!")
  }
  if (ntmt==2) return(c("AB","BA"))
  # Sheehe and Bross's algorithm
  # start vector by random
  sv <- sample(1:ntmt)
  # cyclic Latin square
  des   <- matrix(nrow=ntmt, ncol=ntmt)
  for (i in 1:ntmt){
   if(i==1) des[i,] <- sv[i:ntmt] else des[i,] <-c(sv[i:ntmt],sv[1:(i-1)]) 
  }
  # mirror the matrix
  mir <- des
  for (i in 1:ntmt){
    mir[i,] <- rev(des[i,])
  }
  # interlace both
  inter <- matrix(nrow=ntmt, ncol=2*ntmt)
  cnt <- 0
  for (i in 1:ntmt){
    cnt <- cnt + 1
    inter[,cnt] <- des[,i]
    cnt <- cnt + 1
    inter[,cnt] <- mir[,i]
  }
  # slice at middle
  m1 <- inter[,1:ntmt]
  m2 <- inter[,(ntmt+1):(2*ntmt)]
  
  # if ntmt is odd combine both
  if (2*trunc(ntmt/2) != ntmt) m1 <- rbind(m1, m2)
  
  # order by first column
  m1   <- m1[order(m1[,1]),]
  # coded with AB...
  m1c  <- matrix(ncol=ntmt, LETTERS[m1])
  
  seqs <- vector(mode="character", length = nrow(m1c))
  for (i in seq_along(m1c[,1])){
    seqs[i] <- paste(m1c[i,],sep="",collapse="")
  }
  # order lexically?
  seqs <- sort(seqs)
  # if tmts are given
  if (! is.null(tmts)){
    #replace the a codes by tmts
    for (i in seq_along(tmts)){
      st <- LETTERS[i]
      seqs <- gsub(st,tmts[i],seqs)
    }
  }
  
  seqs
  
}