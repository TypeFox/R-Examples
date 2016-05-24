MVA.ident <- function(x,...) {
  x.class <- MVA.class(x)
  class(x) <- c(class(x),x.class)
  return(x)
}

MVA.class <- function(x,...) {
  UseMethod("MVA.class")
}

MVA.class.default <- function(x,...) {
  if (is.list(x) && "GOF" %in% names(x)) {res <- "PCoA.stats"} else
  if (is.list(x) && "stress" %in% names(x)) {res <- "nMDS.MASS"} else
  {res <- "unknown"}
  return(res)
}

MVA.class.dudi <- function(x,...) {
  if (inherits(x,"pca")) {res <- "PCA.ade4"} else
  if (inherits(x,"pco")) {res <- "PCoA.ade4"} else
  if (inherits(x,"pcaiv")) {
    if (inherits(x,"cca")) {res <- "CCA.ade4"} else
	{res <- "RDA.ade4"}
  } else
  if (inherits(x,"pcaivortho")) {res <- "RDAortho.ade4"} else
  if (inherits(x,"nsc")) {res <- "NSCOA.ade4"} else
  if (inherits(x,"dec")) {res <- "DCOA.ade4"} else
  if (inherits(x,"acm")) {res <- "MCA.ade4"} else
  if (inherits(x,"mix")) {res <- "Mix.ade4"} else
  if (inherits(x,"coa")) {res <- "COA.ade4"} else
  if (inherits(x,"coinertia")) {res <- "CIA.ade4"}
  return(res)
}

MVA.class.pca <- function(x,...) {
  if (inherits(x,"dudi")) {res <- MVA.class.dudi(x)} else
  if (inherits(x,"prcomp")) {res <- MVA.class.prcomp(x)} else
    {res <- "PCA.labdsv"}
  return(res)
}

MVA.class.prcomp <- function(x,...) {
  if (inherits(x,"pca")) {
    if (inherits(x,"spca")) {res <- "sPCA.mixOmics"} else {res <- "PCA.mixOmics"}
  } else {
    res <- "PCA.prcomp"
  }
  return(res)
}

MVA.class.princomp <- function(x,...) {"PCA.princomp"}

MVA.class.cca <- function(x,...) {
  if (inherits(x,"dudi")) {res <- MVA.class.dudi(x)} else
  if (inherits(x,"rda")) {
    if (inherits(x,"capscale")) {
	if (is.null(x$CCA)) {res <- "PCoA.vegan"} else {res <- "dbRDA.vegan"}
    } else {
	if (is.null(x$CCA)) {res <- "PCA.vegan"} else {res <- "RDA.vegan"}
    }
  } else {
	if (is.null(x$CCA)) {res <- "COA.vegan"} else {res <- "CCA.vegan"}
  }
  return(res)
}

MVA.class.ipca <- function(x,...) {"IPCA.mixOmics"}

MVA.class.sipca <- function(x,...) {"sIPCA.mixOmics"}

MVA.class.pcoa <- function(x,...) {"PCoA.ape"}

MVA.class.pco <- function(x,...) {
  if (inherits(x,"dudi")) {res <- MVA.class.dudi(x)} else {res <- "PCoA.labdsv"}
  return(res)
}

MVA.class.wcmdscale <- function(x,...) {"PCoA.vegan"}

MVA.class.dpcoa <- function(x,...) {"DPCoA.ade4"}

MVA.class.monoMDS <- function(x,...) {"nMDS.mono.vegan"}

MVA.class.metaMDS <- function(x,...) {
  if (inherits(x,"monoMDS")) {res <- MVA.class.monoMDS(x)} else {res <- "nMDS.iso.vegan"}
  return(res)
}

MVA.class.nmds <- function(x,...) {"nMDS.labdsv"}

MVA.class.lda <- function(x,...) {"LDA.MASS"}

MVA.class.discrimin <- function(x,...) {
  if (inherits(x,"coadisc")) {res <- "CDA.ade4"} else {
	if (inherits(eval(x$call$dudi,parent.frame()),"coa")) {res <- "CDA.ade4"} else {res <- "LDA.ade4"}
  }
  return(res)
}

MVA.class.plsda <- function(x,...) {"PLSDA.mixOmics"}

MVA.class.splsda <- function(x,...) {
  if (inherits(x,c("splsda1fact","splsda2fact"))) {res <- "Multilevel.sPLSDA.mixOmics"} else {res <- "sPLSDA.mixOmics"}
  return(res)
}

MVA.class.mvr <- function(x,...) {
  if (x$method %in% c("kernelpls","widekernelpls","simpls","oscorespls")) {res <- "PLSR.pls"} else
  if (x$method=="cppls") {res <- "CPPLS.pls"} else
  if (x$method=="svdpc") {res <- "PCR.pls"}
  return(res)
}

MVA.class.pls <- function(x,...) {
  if (inherits(x,"spls")) {
    res <- MVA.class.spls(x)
  } else {
    if (x$mode=="canonical") {res <- "2BPLS.mixOmics"} else {res <- "PLSR.mixOmics"}
  }
  return(res)
}

MVA.class.spls <- function(x,...) {
  if (inherits(x,"mlspls")) {
    if (x$mode=="canonical") {res <- "Multilevel.2BsPLS.mixOmics"} else {res <- "Multilevel.sPLSR.mixOmics"}
  } else {
    if (x$mode=="canonical") {res <- "2BsPLS.mixOmics"} else {res <- "sPLSR.mixOmics"}
  }
  return(res)
}

MVA.class.plsRmodel <- function(x,...) {"PLSR.plsRglm"}

MVA.class.plsRglmmodel <- function(x,...) {"PLSGLR.plsRglm"}

MVA.class.procuste <- function(x,...) {"PCIA.ade4"}

MVA.class.CCorA <- function(x,...) {"CCorA.vegan"}

MVA.class.rcc <- function(x,...) {"rCCorA.mixOmics"}

MVA.class.rgcca <- function(x,...) {
  if ("variates" %in% names(x)) {res <- "rGCCA.mixOmics"} else {res <- "rGCCA.RGCCA"}
  return(res)
}

MVA.class.sgcca <- function(x,...) {
  if ("variates" %in% names(x)) {res <- "sGCCA.mixOmics"} else {res <- "sGCCA.RGCCA"}
  return(res)
}

