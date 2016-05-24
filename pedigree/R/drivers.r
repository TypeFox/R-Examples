countOff <- function(ped)
  {
    if(!is.data.frame(ped))stop("ped should be data.frame")
    ord <- orderPed(ped)
    if(!identical(ord,1:nrow(ped)))
        stop("pedigree is not ordered.")
    ped <- ped[order(ord),]
    idNames <- ped[,1]
    id <- 1:nrow(ped)
    dam <- match(ped[,2],ped[,1],nomatch = 0)
    sire <- match(ped[,3],ped[,1],nomatch = 0)
    n <- length(id)
    nOff <- .C("countOff",ind = as.integer(id),dam = as.integer(dam),sire = as.integer(sire),
               n = as.integer(n),nOff = as.integer(rep(0,n)))$nOff
    return(nOff[ord])
  }

calcInbreeding <- function(ped)
  {
    if(!is.data.frame(ped))stop("ped should be data.frame")
    ord <- orderPed(ped)
    if(!identical(ord,1:nrow(ped)))
        stop("pedigree is not ordered.")
    ped <- ped[order(ord),]
    idNames <- ped[,1]
    id <- 1:nrow(ped)
    dam <- match(ped[,2],ped[,1],nomatch = 0)
    sire <- match(ped[,3],ped[,1],nomatch = 0)
    n <- length(id)
    f <- .C("calcInbreeding",ind = as.integer(id),dam = as.integer(dam),sire = as.integer(sire),
            n = as.integer(n),f = as.double(rep(0,n)))$f
    return(f[ord])
  }

orderPed <- function(ped)
  {
    if(!is.data.frame(ped))
        stop("ped should be data.frame")

    id <- 1:nrow(ped)
    dam <- match(ped[,2],ped[,1],nomatch = 0)
    sire <- match(ped[,3],ped[,1],nomatch = 0)
    n <- length(id)
    ord <- .C("orderPed",ind = as.integer(id),dam = as.integer(dam),sire = as.integer(sire),
              n = as.integer(n),order = as.integer(rep(0,n)))$order
    if(-1%in%ord)
        warning("Be carefull, there are loops in the pedigree, individuals involved in the loop are indicated with a -1\n")
    return(ord)
  }

countGen <- function(ped)
  {
    if(!is.data.frame(ped))
        stop("ped should be data.frame")
    ord <- orderPed(ped)
    if(!identical(ord,1:nrow(ped)))
        stop("pedigree is not ordered.")
    id <- 1:nrow(ped)
    dam <- match(ped[,2],ped[,1],nomatch = 0)
    sire <- match(ped[,3],ped[,1],nomatch = 0)
    n <- length(id)
    .C("countGen",ind = as.integer(id),dam = as.integer(dam),sire = as.integer(sire),
       n = as.integer(n),gen = as.integer(rep(0,n)))$gen
  }

makeAinv <- function(ped)
  {
    if(!is.data.frame(ped))stop("ped should be data.frame")
    ord <- orderPed(ped)
    if(!identical(ord,1:nrow(ped)))
      stop("pedigree is not ordered.")
    idNames <- ped[,1]
    id <- 1:nrow(ped)
    dam <- match(ped[,2],ped[,1],nomatch = 0)
    sire <- match(ped[,3],ped[,1],nomatch = 0)
    n <- length(id)
    res <- .C("getAinv",ind = as.integer(id),dam = as.integer(dam),sire = as.integer(sire),
              n = as.integer(n))
    return(TRUE)
  }

makeA <- function(ped,which)
  {
    if(!is.data.frame(ped))stop("ped should be data.frame")
    if(!is.logical(which))stop("which should be a logical")
    if(length(which)!=nrow(ped))stop("length which should coincide with nrow(ped)")
    which <- as.numeric(which)
    ord <- orderPed(ped)
    if(!identical(ord,1:nrow(ped)))
      stop("pedigree is not ordered.")
    idNames <- ped[,1]
    id <- 1:nrow(ped)
    dam <- match(ped[,2],ped[,1],nomatch = 0)
    sire <- match(ped[,3],ped[,1],nomatch = 0)
    n <- length(id)
    res <- .C("getA",ind = as.integer(id),dam = as.integer(dam),sire = as.integer(sire),
              n = as.integer(n),which = as.integer(which))
    return(TRUE)
  }

trimPed <- function(ped,data,ngenback = NULL)
  {
    if(!is.data.frame(ped))
        stop("ped should be data.frame")
    if(length(data) != nrow(ped))
        stop("length of data should coincide with nrow of pedigree")
    ord <- orderPed(ped)
    if(!identical(ord,1:nrow(ped)))
      stop("pedigree is not ordered.")
    id <- 1:nrow(ped)
    dam <- match(ped[,2],ped[,1],nomatch = 0)
    sire <- match(ped[,3],ped[,1],nomatch = 0)
    data <- as.integer(data>0)
    n <- length(id)
    if(is.null(ngenback))
      ngenback <- as.integer(max(countGen(ped)))
    .C("trimPed",ind = as.integer(id),dam = as.integer(dam),sire = as.integer(sire),
       data = as.integer(data),ngenback = as.integer(ngenback),n = as.integer(n))$data==1
  }


