quitarrep <-
function(matriz){
    if(is.null(dim(matriz))){
      a2=c(FALSE,matriz[-1]==matriz[-length(matriz)])
    }else{
      if(nrow(matriz)>=2){
        a2=c()
        for(k in 1:ncol(matriz)){
          a2c=FALSE
          for(j in 2:nrow(matriz)){ 
            if(k>=2){
              a2c=c(a2c,matriz[j,k]==matriz[j-1,k] & (a2[j,k-1]))
            }else{
              a2c=c(a2c,matriz[j,k]==matriz[j-1,k])
            }  
          }
          a2=cbind(a2,a2c)
        }
      }else{
        a2=rep(TRUE,ncol(matriz))
      }  
    }
    matriz[a2]=NA
    return(matriz)
}
