"compoundInterest" <-
function(interest, periods=1, 
                frequency=1, net.value=FALSE){
#
  cI <- exp(periods*frequency*
      log(1+interest/frequency))
  cI[frequency==Inf] <- exp(periods*interest)
#
  if(net.value)cI-(net.value)
  cI
}

