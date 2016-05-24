# gfac version first index fastest
gfac2<-function(covlevels){

      nrows<-1
      ncov<-length(covlevels)
      levmult <- cumprod(covlevels)
      totlev<-levmult[ncov]
      levmult<-c(1,levmult[-ncov])

      # generate covariates

      cov<-NULL
      for (j in 1:ncov) {
          scov<-gl(covlevels[j],levmult[j],totlev*nrows)
          cov<-cbind(cov,scov)
      }
      cov
}
