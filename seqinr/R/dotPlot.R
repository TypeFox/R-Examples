dotPlot <- function(seq1, seq2, wsize = 1, wstep = 1, nmatch = 1, col = c("white", "black"), 
xlab = deparse(substitute(seq1)), ylab = deparse(substitute(seq2)), ...){
 #
 # Check arguments:
 #
 if(nchar(seq1[1]) > 1) stop("seq1 should be provided as a vector of single chars")
 if(nchar(seq2[1]) > 1) stop("seq2 should be provided as a vector of single chars")
 if(wsize < 1) stop("non allowed value for wsize")
 if(wstep < 1) stop("non allowed value for wstep")
 if(nmatch < 1) stop("non allowed value for nmatch")
 if(nmatch > wsize) stop("nmatch > wsize is not allowed")
 #
 # sliding window on sequences:
 #
 mkwin <- function(seq, wsize, wstep){
  sapply(seq(from = 1, to = length(seq) - wsize + 1, by = wstep), function(i) c2s(seq[i:(i + wsize - 1)]))
 }
 wseq1 <- mkwin(seq1, wsize, wstep)
 wseq2 <- mkwin(seq2, wsize, wstep)
 if( nmatch == wsize ){
   # perfect match case
   xy <- outer(wseq1, wseq2, "==")
 } else {
   # partial match case
   "%==%" <- function(x, y) colSums(sapply(x, s2c) == sapply(y, s2c)) >= nmatch
   xy <- outer(wseq1, wseq2, "%==%")
 }
 image(x = seq(from = 1, to = length(seq1), length = length(wseq1)), 
       y = seq(from = 1, to = length(seq2), length = length(wseq2)),
       z = xy, col = col, xlab = xlab, ylab = ylab, ...)
  box()
}
