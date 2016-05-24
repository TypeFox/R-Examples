N2Exact<-function(p=rep(0.25,4))
{
    if (!is.vector(p)) 
	{
        stop("a vector is expected for 'p'")
    }
   
   if (sum(p)< 0.99|sum(p)>1.01)
   {
        stop("sum of allele frequencies must be between 0.99 and 1.01")
   }
        if (!is.numeric(p) || any(is.na(p)) || any(p < 0) || 
            any(p > 1)) {
            stop("alleles frequencies in  'p' must be numbers strictly between 0 and 1")
        }
    m<-length(p)
    alfa1<-sum(p^4)

    alfa2<-alfa3<-alfa4<-0
    for(i in 1:m){
         for (j in setdiff(1:m,i)){
         alfa2<-alfa2+p[i]*p[j]*(2*p[i]^2+3*p[i]*p[j]+2*p[j]^2)
    }}

    for(i in 1:m){
         for (j in setdiff(1:m,i)){
              for (k in setdiff(1:m,c(i,j))){
                 alfa3<-alfa3+2*p[i]*p[j]*p[k]*(p[i]+p[j]+p[k])
     }}}

     for(i in 1:m){
         for (j in setdiff(1:m,i)){
              for (k in setdiff(1:m,c(i,j))){
                  for (l in setdiff(1:m,c(i,j,k))){
                    alfa4<-alfa4+p[i]*p[j]*p[k]*p[l]
     }}}}


     alfa<-c(alfa1,alfa2,alfa3,alfa4)
     names(alfa)<-c("P(N=1)","P(N=2)","P(N=3)","P(N=4)")
     alfa
}




