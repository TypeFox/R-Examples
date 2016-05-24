GCD=function(x,y,k=0){
  if (y!=0){
    k=k+1
    GCD(y,x%%y,k)
  }else{    
    result=list(k=k,g=x)
    return(result)
  }
  
}