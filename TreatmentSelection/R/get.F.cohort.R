get.F.cohort <-
function(marker, event, trt, rho, return.fun = FALSE){
  # rank(marker, ties.method="max")/length(marker) #older way of doing it
 
  # remember that ecdf defines a function
  if(!return.fun){
  
     ecdf(marker)(marker)
  
  }else {
  
     ecdf(marker)
  
  }

}
