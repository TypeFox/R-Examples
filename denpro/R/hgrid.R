hgrid<-function(h1,h2,lkm,base=10)
{
step<-(h2-h1)/(lkm-1)

if (is.null(base)){
   hseq<-seq(h2,h1,-step)
}
else{
   a<-(h2-h1)/(base^(h2)-base^(h1))
   b<-h1-a*base^(h1)
   un<-seq(h2,h1,-step)
   hseq<-a*base^(un)+b
}

return(hseq)
}

