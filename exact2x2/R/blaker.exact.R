`blaker.exact` <-
function(x, y = NULL, or = 1, alternative = "two.sided", 
    conf.int = TRUE, conf.level = 0.95, tol=0.00001){
   call <- as.list(match.call())
   OUT<-do.call("exact2x2",c(call[-1],tsmethod="blaker"),envir=parent.frame())
   OUT
}

