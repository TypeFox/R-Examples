"shannonbio" <-
function(data1){

    tot<-tapply(data1[,2],data1[,1],sum)
    tot<-tot[!is.na(tot)]
    tot<-tot/sum(tot)
    tot<-tot*(log(tot))/log(2)
    h<-sum(tot[is.finite(tot)])
    hmax<-log(1/length(tot))/log(2)
    c(H=-h,J=h/hmax)
}
