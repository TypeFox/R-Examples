JacksonParameters<-function(nPCs, TFBS){
    
    modelAll<-prcomp(TFBS)
    sdevAll<-modelAll$sdev
    sdevAll<-sdevAll[(nPCs+1):length(sdevAll)]
    sigma1<-sum(sdevAll)
    sigma2<-sum(sdevAll^2)
    sigma3<-sum(sdevAll^3)
    h0<-1-(2*sigma1*sigma3)/(3*sigma2*sigma2)
    x1<-1/sigma1
    x2<- -h0*(h0-1)*sigma2/(sigma1*sigma1)-1

    x3<-sigma1/sqrt(2*sigma2*h0*h0)
    y<-list(h0=h0, x1=x1, x2=x2, x3=x3)
    return(y)


}
