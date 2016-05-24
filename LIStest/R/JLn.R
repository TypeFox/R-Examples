JLn <-
function(x,y){
    jklid<-jkl<-0  
    M<-length(x)   
    if(M >= 20){
      d<-x              
      D<-ceiling(0.4*M^(5/6))
      for(i in 1:M){d[i]<-abs(rank(x)[i]-rank(y)[i])} 
      v<-x[d<D]
      w<-y[d<D]
      N<-length(v)
    } else{
      v<-x
      w<-y
      N=M
    } 
    if (N>1){ jklid<-Ln(v[2:N],w[2:N])
              l<-Ln(v[1:(N-1)],w[1:(N-1)]) 
              jklid<-jklid+l
    } 
    if (N>2){
        for(i in 2:(N-1)){ 
            l<-Ln(c(v[1:(i-1)],v[(i+1):N]),c(w[1:(i-1)],w[(i+1):N]))
            jklid<-jklid+l
        }        
    }
    if(N==0){ jkl=0 }
    else if(N==1){jkl<-1}
    else{ 
        jkl<-jklid/N  }
    res<-jkl 
}
