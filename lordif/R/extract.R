extract <-
function(ipar) {
    calib<-coef(ipar,IRTpars=TRUE)
    ni<-length(calib)-1
    calib<-calib[1:ni]
    ncat<-unlist(lapply(calib,length)) 
    maxCat<-max(ncat)
    out.ipar<-data.frame(matrix(NA,ni,maxCat))
    names(out.ipar)<-c("a",paste0("cb",1:(maxCat-1)))
    for (i in 1:ni) {
      out.ipar[i,1]<-calib[[i]][1]
      for (j in 1:(ncat[i]-1)) {
        out.ipar[i,j+1]<-calib[[i]][j+1]
      }
    }
    if (all(out.ipar[,1]<0)) {
      out.ipar<- -out.ipar
    } else if (any(out.ipar[,1]<0)) {
      cat("ERROR: The following items had negative slope parameters.\n")
	  cat(paste0(which(out.ipar[,1]<0),collapse=","))
	  cat(paste0("ERROR: The following items had negative slope parameters (",paste0(which(out.ipar[,1]<0),collapse=","),").\n"))
      cat("\n")
    }
    return(out.ipar)
  }
