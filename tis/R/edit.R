edit.tis <- function(name, transpose = F, ...){
  ## Create a matrix that isn't a time series
  inmat <- unclass(cbind(name))
  tsp(inmat) <- NULL

  ## Fix up some dimnames for name
  dn.t <- dimnames(name)
  dimnames(inmat)[[1]] <- ymd(ti(name))
  if((length(dn.t) != 2) || (length(dn.t[2]) == 0)){
    nstr <- deparse(substitute(name))
    if(ncol(inmat) > 1)
      dimnames(inmat)[[2]] <- paste(nstr, 1:ncol(inmat), sep = ".")
    else
      dimnames(inmat)[[2]] <- nstr
  }

  if(transpose){
    z <- t(edit(t(inmat), ...))
  }
  else {
    z <- edit(inmat, ...)
  }
  mode(z) <- mode(name)
  outser <- tis(z, start = as.numeric(dimnames(z)[[1]][1]), tif = tif(name))

  ## May drop dim attribute for outser
  o.m  <- ncol(outser)
  o.t  <- max(ncol(name), 1, na.rm = T)
  if((o.m == 1) && is.null(dim(name))) 
    dim(outser) <- NULL

  ## Fix up dimnames for outser
  dn.o <- dimnames(outser)
  if(length(dn.t) != 2)  
    dimnames(outser) <- NULL
  else
    dimnames(outser) <- list(character(0), c(dn.o[[2]], rep("", o.m))[1:o.m])
  outser
}

edit.ti <- function(name, header = character(0), ...){
  ## edit a ti
  z <- asTi(as.matrix(name))
  inmat <- matrix(paste(ymd(z), format(tifName(z)), sep = " : "), ncol = ncol(z))
  if(!is.null(rn <- rownames(z))) rownames(inmat) <- rn
  if(!is.null(cn <- colnames(z))) colnames(inmat) <- cn
  outmat <- edit(inmat, header = c(header,
                                  "ti dates shown as ymd : tifName"), ...)
  splitlist <- strsplit(outmat, split = ":")
  ymds <- as.numeric(stripBlanks(sapply(splitlist, "[", 1)))
  tifs <- stripBlanks(sapply(splitlist, "[", 2))
  ans <- ymds*0
  for(i in 1:length(ans)) ans[i] <- ti(ymds[i], tifs[i])
  ans <- asTi(ans)
  if(is.matrix(name)){
    dim(ans) <- dim(outmat)
    dimnames(ans) <- dimnames(outmat)
  }
  else
    if(length(rn <- rownames(outmat)) == length(ans))
      names(ans) <- rn
  ans
}
