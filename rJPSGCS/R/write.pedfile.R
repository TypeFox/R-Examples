write.pedfile <-
function(pedinf, snp.data, file, transpose=FALSE, sep=" ", eol="\n", na="0") {
  # Error checking
  if (!is(snp.data, "snp.matrix")) {
    stop("argument snp.data must be a snp.matrix object")
  }
  if(length(pedinf)==1 && pedinf=="unrelated") {
    pedinf<-make.unrel.pedinf(nrow(snp.data))
  }
  if(is.data.frame(pedinf)) {
    pedinf <- as.matrix(pedinf)
  }
  if(nrow(pedinf) != nrow(snp.data)) {
    stop("pedinf and snp.data must have the same number of rows")
  }
  for(i in 5:dim(pedinf)[2]){
    if (!(all(pedinf[,i] %in% c(0,1,2)))) {
        stop("pedinf columns 5 onward must be numeric with values 0, 1 or 2")
    }
  }

  # Done error checking, now write appropriate file
  if(transpose) {
    res <- .C("write_as_transposed_pedfile", as.character(file),
       as.integer(pedinf), snp.data@.Data, as.integer(nrow(snp.data)), as.integer(ncol(snp.data)),
       as.integer(ncol(pedinf)), as.character(sep), as.character(eol),
       as.character(na), logical(1), PACKAGE="rJPSGCS")
  } else {
    res <- .C("write_as_pedfile", as.character(file),
       as.character(pedinf), snp.data@.Data, as.integer(nrow(snp.data)), as.integer(ncol(snp.data)),
       as.integer(ncol(pedinf)), as.character(sep), as.character(eol),
       as.character(na), logical(1), PACKAGE="rJPSGCS")
  }
  error <- res[[9]]
  if (error==1)
    stop("Couldn't open output file")
  else
    c(nrow(snp.data), ncol(snp.data))
}


# function that creates post-MAKEPED pedfile from SNP data on 
# unrelated individuals
make.unrel.pedinf=function(n){

    # adapted from pedfunc function in CrypticIBDcheck

    #snobj: a snp.matrix object
    #file: name of the ped file to be generated

    ones <- rep(1,n)
    zeros <- rep(0,n)
    ids <- (1:n)
    pedinf <- cbind(ids,ids,zeros,zeros,zeros,zeros,zeros,ones,ones)
    return(pedinf)
}

