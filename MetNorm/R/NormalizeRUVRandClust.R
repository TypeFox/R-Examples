NormalizeRUVRandClust<-function(RUVRand,
                              maxIter, 
                              nUpdate=maxIter+1, 
                              lambdaUpdate=TRUE, 
                              p=p,...){
  output<-RuvRandIter(RUVRand=RUVRand,
                          maxIter=maxIter, 
                          wUpdate=nUpdate, 
                          lambdaUpdate=lambdaUpdate, 
                          p=p,...)
  return(structure(output, class="normdata"))   
}

                        