calcSelcrit <-
function(allmodeldata,include,sd.target) {
  allmodeldata$selcrit <- rep(NA,nrow(allmodeldata)) 
  include[is.na(include)] <- F
  if (sum(include)>0) {
    if (is.na(sd.target)) {
      allmodeldata$selcrit[include] <- allmodeldata$BIC[include]
    } else {
      suppressWarnings({ maxbic <- max(allmodeldata$BIC[include],na.rm=T) })
        #otherwise warning if no data; which is not a problem as maxbic becomes -Inf
      if (is.na(maxbic) || maxbic>0) { #doesn't work with BIC>0; anyway BIC>0 is not common
        allmodeldata$selcrit[include] <- allmodeldata$BIC[include]
      } else {
        sdvalues <- allmodeldata$sdtrans0 #all sdtrans are equal, take the first
        mulf <- rep(1,nrow(allmodeldata))
        mulf[!is.na(sdvalues) & sdvalues>sd.target] <- sd.target / sdvalues[!is.na(sdvalues) & sdvalues>sd.target] #larger sdvalue -> smaller mulf
        allmodeldata$selcrit[include] <-  allmodeldata$BIC[include] * mulf[include] #with larger sdvalue selcrit moves towards zero (becomes bigger, as BIC is negative)
      }  
    }
  } 
  allmodeldata   
}
