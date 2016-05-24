opt_trophy<-function(M=NULL,N=1000,recage=NULL,entage=NULL,trage=NULL,stF=0,
        endF=2,intF=0.05){
    if(is.null(M)) stop("Natural mortality estimate is missing.")
    if(is.null(recage)) stop("Age at first recruitment is missing.")
    if(is.null(entage)) stop("Age at entry into the exploited stock is missing.")
    if(is.null(trage)) stop("Age when fish become trophy fish is missing.")
      F<-seq(stF,endF,intF)
      total<-(F*N*exp(-M*(entage-recage)))/(F+M)
      trophy<-(F*N*exp(-M*(entage-recage)-(F+M)*(trage-entage)))/(F+M)
      Fmax<-(-M*(trage-entage)+sqrt(M^2*(trage-entage)^2+4*M*(trage-entage)))/(2*(trage-entage))
      outpt<-list(as.data.frame(cbind(F,total,trophy)),Fmax)
      names(outpt)<-c("Catch","Fmax")
     return(outpt)
}