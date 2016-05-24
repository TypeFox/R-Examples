tidy2magpie <- function(x,spatial=NULL,temporal=NULL) {
  # assumption: dataframe format in which only the very last
  #             column contains values!
  sep <- "."
  
  colnames(x) <- make.unique(colnames(x),sep="")
  
  if(dim(x)[1]==0) return(copy.attributes(x,new.magpie(NULL)))
  
  if(is.null(spatial))     spatial  <- colnames(x)[apply(x,2,is.spatial)]
  if(is.null(temporal))    temporal <- colnames(x)[apply(x,2,is.temporal)]                                                                                         
  if(is.numeric(spatial))  spatial  <- colnames(x)[spatial]
  if(is.numeric(temporal)) temporal <- colnames(x)[temporal]
    
  .collapsecol <- function(x,which,sep=".") {
    xname <- paste(colnames(x)[which],collapse=sep)
    args <- list()
    for(i in 1:length(which)) {
      args[[i]] <- x[,which[i]]
    }
    args["sep"] <- sep
    out <- as.data.frame(do.call(paste,args))
    colnames(out) <- xname
    return(out)
  }  
  
  if(sum(colnames(x) %in% temporal)>1) {
    t <- .collapsecol(x,which(colnames(x) %in% temporal),sep)
  } else if(sum(colnames(x) %in% temporal)==1) {
    t <- x[,which(colnames(x) %in% temporal),drop=FALSE]
  } else { 
    t <- data.frame(year=rep("NOTIME",dim(x)[1]))  
  }
  t[[1]] <- as.character(t[[1]])
  
  if(sum(colnames(x) %in% spatial)>1) {
    s <- .collapsecol(x,which(colnames(x) %in% spatial),sep) 
  } else if(sum(colnames(x) %in% spatial)==1) {
    s <- x[,which(colnames(x) %in% spatial),drop=FALSE]
  } else {
    s <- data.frame(region=rep("GLO",dim(x)[1]))
  }
  s[[1]] <- as.character(s[[1]])
  
  if(sum(!(colnames(x)[-dim(x)[2]] %in% c(temporal,spatial)))>1) {
    d <- .collapsecol(x,which(!(colnames(x)[-dim(x)[2]] %in% c(temporal,spatial))),sep)
  } else if(sum(!(colnames(x)[-dim(x)[2]] %in% c(temporal,spatial)))==1) {
    d <- x[,which(!(colnames(x)[-dim(x)[2]] %in% c(temporal,spatial))),drop=FALSE]
  } else {
    d <- data.frame(data=rep("NODATA",dim(x)[1]))
  }
  d[[1]] <- as.character(d[[1]])
  
  u_spat <- as.character(unique(s[,1]))
  u_temp <- as.character(unique(t[,1]))
  u_data <- as.character(unique(d[,1]))
  dimnames <- list(u_spat,u_temp,u_data)
  m <- array(dim=c(length(u_spat),length(u_temp),length(u_data)),dimnames=dimnames)
  coord <- as.matrix(cbind(s,t,d))
  .duplicates_check(coord)
  m[coord] <- x[,dim(x)[2]]
  if(dim(m)[2]==1) if(dimnames(m)[[2]]=="NOTIME") dimnames(m) <- list(dimnames(m)[[1]],NULL,dimnames(m)[[3]])
  if(dim(m)[3]==1) if(dimnames(m)[[3]]=="NODATA") dimnames(m) <- list(dimnames(m)[[1]],dimnames(m)[[2]],NULL)
  
  names(dimnames(m)) <- c(names(s),names(t),names(d))
  m <- as.magpie(m,spatial=1,temporal=2)
  return(copy.attributes(x,m))
}