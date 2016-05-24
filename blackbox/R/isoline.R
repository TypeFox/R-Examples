isoline <- function(latt2Ns2) { ## util. to add Nb levels to an existing (2Nm, g) surface plot
  seq2Nm <- seq(blackbox.getOption("FONKgLow")["twoNm"], latt2Ns2/blackbox.getOption("mincondS2"), length.out=100)
  islog2Ns2 <- latt2Ns2
  if (islogscale("latt2Ns2")) islog2Ns2 <- log(islog2Ns2) ## because tofullKrigingspace then assumes that latt2Ns2 is logscale
  seqg <- sapply(seq2Nm, ## the twoNmu value because tofullK catches (twoNmu=NA & Nratio=NA)
               function(v) {tofullKrigingspace(list(twoNmu=0, twoNm=v), fixedlist=list(latt2Ns2=islog2Ns2))["g"]}
  )
  lines(seq2Nm, seqg, type="l", lty=2)
}
