see<-
function(true,pred) {
  if(length(true)!=length(pred)) stop("Length of 'true' must equal 'pred'")
  res<-true-pred
  sum.sq<-sum(res^2)
  see<-sqrt(sum.sq/length(true))
  return(see)
}
