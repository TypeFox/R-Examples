R_w_ma <-
function(theta,nstart,nlen){

    T=nrow(theta)
    q=ncol(theta)
    q=q-1  
 
    B=matrix(0,nlen,nlen+q)
    n=nstart
    for (j in 1:nlen)
        {time_index=matlab::mod(n-1,T)+1
         B[j,j:(j+q)]=matlab::fliplr(theta[time_index,])
         n=n+1}
     r=B%*%t(B)

     rindex<-NULL
     R<-NULL

     for (i in 1:nlen){
     for (j in 1:nlen){
     if (r[i,j]!=0) 
        {index=c(i,j)
         R=c(R,r[i,j])
         rindex=rbind(rindex,c(i,j))
     } else {
          rindex=rindex
      }
   
    }}
      
    result = list(R= R, rindex=rindex) 
    class(result) = "R_w_ma"
    result 
}

