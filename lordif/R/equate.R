equate <-
function(ipar.to,ipar.from,theta,model="GRM",start.AK=c(1.0,0.0),lower.AK=c(0.5,-2.0),upper.AK=c(2.0,2.0)) {
    if (!(model %in% c("GRM","GPCM"))) {
      warning("model must be either \"GRM\" or \"GPCM\", will be reset to \"GRM\"")
      model<-"GRM"
    }
    SL<-function(AK) {
      a.to<-ipar.to$a
      a.from<-ipar.from$a
      cb.to<-ipar.to[-1]
      cb.from<-ipar.from[-1]
      a.from<-a.from/AK[1]
      cb.from<-cb.from*AK[1]+AK[2]
      return(sum((tcc(a.to,cb.to,theta,model=model)-tcc(a.from,cb.from,theta,model=model))^2))
    }
    return (nlminb(start.AK,SL,lower=lower.AK,upper=upper.AK)$par)
  }
