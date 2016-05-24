"multiLin" <-
function (thetaClass, diffs, parvec) 
{
  cnt <- 1
  newpar <- 0
  if(is.list(diffs$ind2)) {
    for(i in 1:length(parvec)) { 
      newpar <- newpar + ( parvec[[i]] * slot(thetaClass, diffs$what2[[i]])[[diffs$ind2[[i]]]])
    }
  #} 
  #else {
  #  ## this is not working right when what2 indexes a list
  #  if(is.list(slot(thetaClass, diffs$what2))) {
  #    for(i in 1:length(diffs$ind2)) { 
  #    newpar <- newpar + ( parvec[i] * slot(thetaClass, diffs$what2[i])[diffs$ind2[i]])      
  #    }
     
    } else {
      for(i in 1:length(diffs$ind2)) {
        newpar <- newpar + ( parvec[cnt] * slot(thetaClass, 
                                                diffs$what2[i])[diffs$ind2[i]]) 
        cnt <- cnt + 1 
      }
     }
   #}
  newpar 
}
