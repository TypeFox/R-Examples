opt_slot<-function(M=NULL,N=1000,recage=NULL,entage=NULL,trage=NULL,slage=NULL,stF=0,
        endF=2,intF=0.05){
    if(is.null(M)) stop("Natural mortality estimate is missing.")
    if(is.null(recage)) stop("Age at which fish could be recruited is missing.")
    if(is.null(entage)) stop("Age at entry into the exploited stock is missing.")
    if(is.null(trage)) stop("Age when fish become trophy fish is missing.")
     if(is.null(trage)) stop("Age for upper slot is missing.")
       F<-seq(stF,endF,intF)
       trophy<-(F*N*exp(-M*(entage-recage)-(F+M)*(slage-entage)-M*(trage-slage)))/(F+M)
       nontrophy<-(F*N*exp(-M*(entage-recage)))/(F+M)*(1-exp(-(F+M)*(slage-entage)))
       total<-trophy+nontrophy
       Fmax<-(-M*(slage-entage)+sqrt(M^2*(slage-entage)^2+4*M*(slage-entage)))/(2*(slage-entage))
       outpt<-list(as.data.frame(cbind(F,total,nontrophy,trophy)),Fmax)
       names(outpt)<-c("Catch","Fmax")
       return(outpt)
}