
test_TKF92HG <- function(){
  library(seqinr)
  data(GONNET)
  data(GONNETBF)
  fasta <- read.fasta(file.path(system.file("extdata", package="TKF"), 
                                "pair1.fasta"),
                      seqtype="AA", set.attributes=FALSE)

  ## 1D estimation: only distance
  seq1 <- fasta[[1]]
  seq2 <- fasta[[2]]
  ans <- TKF92HGPair(seq1, seq2, mu=5.920655e-04, r=0.8, Ps=1, Kf=1.2,
                   substModel=GONNET, substModelBF=GONNETBF)
  checkEqualsNumeric(119.3832517, ans["PAM"], tolerance=1e-3)

}
