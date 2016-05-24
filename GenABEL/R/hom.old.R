"hom.old" <- 
function(data,snpsubset,idsubset,weight="no") {
  if (is(data,"gwaa.data")) {
    data <- data@gtdata
  }
  if (!is(data,"snp.data")) {
    stop("data argument should have gwaa.data-class or snp.data-class")
  }
  wargs <- c("no","freq")
  if (!(match(weight,wargs,nomatch=0)>0)) {
    out <- paste("weight argument should be one of",wargs,"\n")
    stop(out)
  }

  if (!missing(snpsubset)) data <- data[,snpsubset]
  if (!missing(idsubset)) data <- data[idsubset,]

  totsnps <- data@nsnps
  totids <- data@nids
  opt <- match(weight,wargs)-1
  out <- .C("homold",as.raw(data@gtps),as.integer(data@nids),as.integer(data@nsnps),as.integer(opt),sout = double((2+opt)*data@nids), PACKAGE="GenABEL")$sout
  if (weight=="freq") {
    dim(out) <- c(data@nids,3)
    F <- (out[,2]-out[,3])/(out[,1]-out[,3])
    out <- cbind(out,F)
    out[,2] <- out[,2]/out[,1]
    colnames(out) <- c("NoMeasured","Hom","E(Hom)","F")
  } else {
    dim(out) <- c(data@nids,2)
    out[,2] <- out[,2]/out[,1]
    colnames(out) <- c("NoMeasured","Hom")
  }
  rownames(out) <- data@idnames
  out
}
