cabasicresults<-setClass("cabasicresults",
representation(
  RX="matrix", CX="matrix", Rweights="matrix", Cweights="matrix",
  Raxes="matrix", Caxes="matrix", mu="numeric",mu2="numeric",catype="character",
tauDen="numeric",Z="matrix",ZtZ="matrix",tZZ="matrix"))
