DR_coxph <- function(time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleY=TRUE, plot=FALSE,...){
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
#library(survival)
if((scaleY & missing(time2))){time <- scale(time)}
mf2 <- mf <- match.call(expand.dots = FALSE)
m <- match(c("time", "time2", "event", "type", "origin"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("Surv")
YCsurv <- eval(mf, parent.frame())
if(plot){plot(survival::survfit(YCsurv~1))}
mf2 <- match.call(expand.dots = FALSE)
m2 <- match(c("weighted", "collapse", "origin"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$type <- typeres
mf2$object <- coxph(YCsurv~1,...)
mf2[[1L]] <- as.name("residuals")
return(eval(mf2, parent.frame()))
}
