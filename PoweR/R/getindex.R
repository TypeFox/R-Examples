getindex <- function(law.indices=NULL,stat.indices=NULL) {
  
  if(getRversion() < "3.1.0") dontCheck <- identity
  
  tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])

## We deal with the informations on the laws.
  
  ind.laws <- grep("law",tmp)
  nb.laws <- length(ind.laws)
  
  mat.laws <- as.data.frame(matrix(NA,nrow=nb.laws,ncol=7))
  colnames(mat.laws) <- c("Index","Law","Nbparams","Default1","Default2","Default3","Default4")
  
  getname <- TRUE # to retrieve (or not) the law name. It takes some time to retrieve the law name ...

  for (i in 1:nb.laws) {
    Claw.name <- tmp[ind.laws][i]
    out <- .C(dontCheck(Claw.name),xlen=0L,x=0.0,name=c("1",rep(" ",49)),as.integer(getname),params=rep(0.0,4),nbparams=0L,setseed=0L)
    nbparams <- out$nbparams
    if (nbparams == 0) params <- rep(NA,4)
    if (nbparams == 1) params <- c(out$params[1],rep(NA,3))
    if (nbparams == 2) params <- c(out$params[1:2],rep(NA,2))
    if (nbparams == 3) params <- c(out$params[1:3],NA)
    if (nbparams == 4) params <- out$params
    mat.laws[i,1] <- as.numeric(substring(Claw.name,4))
    mat.laws[i,2] <- gsub('\\','',gsub('$','',sub(' +$', '', paste(out$name,collapse="")),fixed=TRUE),fixed=TRUE)
    mat.laws[i,3] <- as.numeric(nbparams)
    mat.laws[i,4] <- as.double(params[1])
    mat.laws[i,5] <- as.double(params[2])
    mat.laws[i,6] <- as.double(params[3])
    mat.laws[i,7] <- as.double(params[4])
  }

## We deal with the informations on the stats.
  
  ind.stats <- grep("stat",tmp)
  nb.stats <- length(ind.stats)

  mat.stats <- as.data.frame(matrix(NA,nrow=nb.stats,ncol=4))
  colnames(mat.stats) <- c("Index","Stat","Alter","Nbparams")
  mat.stats[,4] <- rep(NA,nb.stats)
  
  for (i in 1:nb.stats) {
    Cstat.name <- tmp[ind.stats][i]
    out <- .C(dontCheck(Cstat.name),0.0,0L,0.0,0L,statname=rep(" ",50),1L,0.0,0L,0.0,0.0,0.0,0L,alter=0L,0L,rep(0.0,4),1L)
    mat.stats[i,1:3] <- c(as.numeric(substring(Cstat.name,5)),sub(' +$', '', paste(out$statname,collapse="")),out$alter)
    mat.stats[i,2] <- gsub('\\','',gsub('$','',sub(' +$', '', paste(mat.stats[i,2],collapse="")),
                                        fixed=TRUE),fixed=TRUE) # Remove $ and backlash signs
    mat.stats[i,3][mat.stats[i,3]==0] <- c("0,1,2")	
  }
  
  mat.stats[,4] <- getnbparstats()

  if (is.null(law.indices) && is.null(stat.indices)) {
    res <- list(mat.laws=mat.laws,mat.stats=mat.stats)
    class(res) <- "index"
    return(res)
  }
  
  if (is.null(law.indices)) {
    res <- list(mat.stats=mat.stats[stat.indices,])
    class(res) <- "index"
    return(res)
  }
  
  if (is.null(stat.indices)) {
    res <- list(mat.laws=mat.laws[law.indices,])
    class(res) <- "index"
    return(res)
  }
  
  res <- list(mat.laws=mat.laws[law.indices,],mat.stats=mat.stats[stat.indices,])
  class(res) <- "index"
  return(res)
  
}

