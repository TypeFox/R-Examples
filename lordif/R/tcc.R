tcc <-
function(a,cb,theta,model="GRM") {
    if (!(model %in% c("GRM","GPCM"))) {
      warning("model must be either \"GRM\" or \"GPCM\", will be reset to default")
      model<-"GRM"
    }
    ni<-length(a)
    T<-numeric(length(theta))
    if (model=="GPCM") {
      for (i in 1:ni) {
        T<-T+probgpcm(theta,a[i],cb[i,])%*%seq(0,sum(!is.na(cb[i,])))
      }
    } else {
      for (i in 1:ni) {
        T<-T+probgrm(theta,a[i],cb[i,])%*%seq(0,sum(!is.na(cb[i,])))
      }
    }
    return(T)
  }
