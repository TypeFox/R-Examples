CPposition <-
function(nbCP,CPpositionPostDist,CPsamples,segMinLength){
    if(nbCP>0) {
      totalPos=length(CPpositionPostDist)
      nbCPtoFind=nbCP
      CPtemp=c(CPsamples[1,1],totalPos)
      orderCPtemp=order(CPpositionPostDist,decreasing=TRUE)
      
      while(nbCPtoFind>0){
        while(sum(abs(CPtemp-orderCPtemp[1])<segMinLength)>0){
          orderCPtemp=orderCPtemp[-1]
          if(length(orderCPtemp)==0){
            stop(paste("Please choose a larger value for segMinLength: there is no CP position vector of estimated length such that the minimum segment length is below",segMinLength))
          }

          
        }

        CPtemp=sort(c(CPtemp,orderCPtemp[1]))

          
        orderCPtemp=orderCPtemp[-1]
        nbCPtoFind=nbCPtoFind-1
      }
      CPpos=sort(CPtemp)
      
    }else{
      CPpos=c(CPsamples[1,1],max(CPsamples[1,]))	
    }
 
 return(CPpos)
}
