adj2A <-
function(B,type="full"){
  if(nrow(B)!=ncol(B)){
    print("not a symmetric matrix")
    return;
  }
  A=diag(1,nrow=nrow(B),ncol=ncol(B))
  for(i in 1:(nrow(B)-1)){
    for( j in (i+1):ncol(B)){
      
      if(type=="full"){
        tmp=rep(0,nrow(B))
        tmp[c(i,j)]=1
        A=cbind(A,tmp)	
        
      }else if(B[i,j]==1){
        tmp=rep(0,nrow(B))
        tmp[c(i,j)]=1
        A=cbind(A,tmp)				
      }
      
    }
    
    
  }
  #colnames(A)=c()
  return(A)
  
}
