.dimSums_fallback<-function (x, na.rm = FALSE, dims = NULL, sep = ".", ...) 
{
  dims <- as.integer(round((dims-3)*10 + 2))
  #fallback option if data dimension is sparse
  if (!is.magpie(x)) stop("Input is not a MAgPIE object!")
  f <- fulldim(x)
  ndim <- length(f[[1]])
  if (any(dims > ndim)) stop("Invalid dimension(s) specified")
  
  if(any(dims<3)) stop("not implemented yet")
  
  #remove names from dimensions that should be summed up
  s <- paste0("^",paste(rep("([^\\.]*)",ndim-2),collapse="\\."),"$")
  d <- setdiff(1:ndim,dims)-2
  d <- d[d>0]
  r <- paste0("\\",d,collapse="\\.")
  
  getNames(x) <- sub(s,r,getNames(x))
  getSets(x,fulldim=FALSE)[3] <- sub(s,r,getSets(x,fulldim=FALSE)[3])
  if(getSets(x,fulldim=FALSE)[3]=="") getSets(x,fulldim=FALSE)[3] <- "data1"
  un <- unique(getNames(x))
  
  out <- new.magpie(cells_and_regions=getCells(x),years=getYears(x),names=un,sets=getSets(x))
  names(dimnames(out)) <- names(dimnames(x))
  x <- as.array(x)
  for(i in un) {
    j <- which(dimnames(x)[[3]]==i)
    out[,,i] <- dimSums(x[,,j,drop=FALSE],na.rm=na.rm,dim=3,sep=sep,...)
  }
  if(ndata(out)==1) if(getNames(out)=="") getNames(out) <- NULL
    
  return(out)
}

dimSums<-function (x, na.rm = FALSE, dims = NULL, dim = 3, sep = ".", ...) 
{
  if(!is.null(dims)) {
    warning('Argument "dims" is depreceated, please use "dim" instead. See ?dimSums for more information!')
    if(is.character(dims)) {
      dim <- dims
    } else {
      dim <- dims
      dim[dim>=3] <- 3 + (dim[dim>=3]-2)/10
    }
  }
  if (is.magpie(x)){
    dim <- dimCode(dim,x)
    if(prod(fulldim(x)[[1]])!=prod(dim(x))) {
      if(any(dim>3)) {
        tmp <- .dimSums_fallback(x, na.rm = na.rm, dims = dim[dim>3], sep = sep, ...)
      } else {
        tmp <- x
      } 
      dim <- dim[dim<=3]
      if(length(dim)==0) {
        return(tmp)
      } else {
        tmp <- as.array(tmp)
      }
    } else {
      if(all(dim!=3)) {
        tmp <- unwrap(x)
        dim[dim>3] <- as.integer(round((dim[dim>3]-3)*10+2))
      } else {
        tmp <- as.array(x)
        dim <- dim[dim<=3]
      }
    }
  } else if (is.array(x)) {
    tmp<-x
  } else {
    stop("Input is neiter an array nor a MAgPIE object!")
  }
  
  if (any(dim > length(dim(tmp)))) 
    stop("Invalid dimension(s) specified")
  unchanged_dims <- which(!1:length(dim(tmp)) %in% dim)
  out <- aperm(tmp, perm = c(unchanged_dims, dim))
  out <- rowSums(out, na.rm = na.rm, dims = length(unchanged_dims), 
                 ...)
  remaining_dims <- match(1:length(dim(tmp)), unchanged_dims, 
                          nomatch = 0)
  remaining_dims <- remaining_dims[remaining_dims > 0]
  
  if (is.magpie(x)){
    spatial <- ifelse(1 %in% dim, 0, 1)
    temporal <- ifelse(2 %in% dim, 0, ifelse(1 %in% dim, 1, 2))
    out <- as.magpie(aperm(as.array(out), perm = remaining_dims),spatial=spatial,temporal=temporal)
    if (1 %in% dim && nregions(x) == 1) 
      dimnames(out)[[1]] <- getRegions(x)
    if (2 %in% dim && nyears(x) == 1) 
      dimnames(out)[[2]] <- getYears(x)
    out <- clean_magpie(out)
  } else {
    out <- aperm(as.array(out), perm = remaining_dims)
  }
  return(out)
}