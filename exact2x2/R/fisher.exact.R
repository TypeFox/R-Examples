`fisher.exact` <-
function(x, y = NULL, or = 1, alternative = "two.sided",
    tsmethod="minlike", 
    conf.int = TRUE, conf.level = 0.95, tol=0.00001){
   call<-as.list(match.call())
   if (!is.null(call$tsmethod) && call$tsmethod=="blaker") warning("output is blaker's exact test")
   OUT<-do.call("exact2x2",call[-1],envir=parent.frame())
   OUT
}

