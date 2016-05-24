CohensD.unpairedT <-
function(t.val,n1,n2){
  d<-t.val*sqrt((n1+n2)/(n1*n2))
  names(d)<-"d"
  return(d)
}
