getnbparstats <- function(stat.indices=NULL) {

  if(getRversion() < "3.1.0") dontCheck <- identity
  
  tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
  ind.stats <- grep("stat",tmp[grep("stat",tmp)])
  if (!all(stat.indices %in% ind.stats)) stop(paste("The values",paste(stat.indices[which(!(stat.indices %in% ind.stats))],collapse=" ")," in 'stat.indices' do not correspond to any defined statistic!",collapse=""))

  if (is.null(stat.indices)) stat.indices <- ind.stats

  lst <- length(stat.indices)
  nbparstats.list <- as.vector(rep(0,lst))
  
  for (i in 1:lst) {
    Cstat.name <- paste("stat",stat.indices[i],sep="")
    nbparstats.list[i] <- (.C(dontCheck(Cstat.name),0.0,1L,0.0,1L,rep(" ",50),1L,0.0,0L,0.0,0.0,0.0,0L,0L,0L,0.0,nbparamstat=0L,PACKAGE="PoweR"))$nbparamstat
  }
  
  return(nbparstats.list)
  
}
