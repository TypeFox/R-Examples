calcprob <-
function(ipar,theta,model="GRM"){
    if (!(model %in% c("GRM","GPCM"))) {
      warning("model must be either \"GRM\" or \"GPCM\"; will be reset to default")
      model<-"GRM"
    }
    ni<-nrow(ipar) 
    maxCat<-ncol(ipar) 
    NCAT<-apply(ipar,1,function (x) sum(!is.na(x))) 
    DISC<-ipar[,"a"] 
    CB<-ipar[,paste0("cb",1:(maxCat-1)),drop=FALSE] 
    pp<-array(NA,c(length(theta),ni,maxCat)) 
    if (model=="GPCM"){
      for (i in 1:ni) {
        pp[,i,1:(NCAT[i])]<-probgpcm(theta,DISC[i],CB[i,])
      }
    } else {
      for (i in 1:ni) {
        pp[,i,1:(NCAT[i])]<-probgrm(theta,DISC[i],CB[i,])
      }
    }
    return(pp)
  }
