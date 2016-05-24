`mcnemar.exact` <-
function(x, y = NULL, conf.level = 0.95){
   call<-as.list(match.call())
   OUT<-do.call("exact2x2",c(call[-1], alternative="two.sided", tsmethod="central",paired=TRUE),
        envir=parent.frame())
   OUT
}
