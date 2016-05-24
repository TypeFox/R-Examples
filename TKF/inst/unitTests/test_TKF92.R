
test_TKF92 <- function(){
  library(seqinr)
  data(GONNET)
  data(GONNETBF)
  fasta <- read.fasta(file.path(system.file("extdata", package="TKF"), 
                                "pair1.fasta"),
                      seqtype="AA", set.attributes=FALSE)

  ## 1D estimation: only distance
  seq1 <- fasta[[1]]
  seq2 <- fasta[[2]]
  ans <- TKF92Pair(seq1, seq2, mu=0.0006137344, r=0.7016089061,
                   substModel=GONNET, substModelBF=GONNETBF)
  checkEqualsNumeric(116.6130887028, ans["PAM"], tolerance=1e-3)

}
