
test_TKF91 <- function(){
  library(seqinr)
  data(GONNET)
  data(GONNETBF)
  fasta <- read.fasta(file.path(system.file("extdata", package="TKF"), 
                                "pair1.fasta"),
                      seqtype="AA", set.attributes=FALSE)

  ## 1D estimation: only distance
  seq1 <- fasta[[1]]
  seq2 <- fasta[[2]]
  ans <- TKF91Pair(seq1, seq2, mu=5.920655e-04,
                   substModel=GONNET, substModelBF=GONNETBF)
  checkEqualsNumeric(116.3416784006, ans["PAM"], tolerance=1e-3)
}
