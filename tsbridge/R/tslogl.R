tslogl <-
function(bug, ymean, sigma){
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  n<-bug$info$n
  beg<-bug$info$args$beg
  y<-bug$data$y[beg:n]
  #need the NA remove for missing data
  ll<-sum(dnorm(y, mean = ymean , sd = sigma, log = TRUE),na.rm=TRUE)
  ll
}
