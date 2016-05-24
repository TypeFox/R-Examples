coxpls3.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, nt=min(7,ncol(Xplan)), typeVC="none", plot=FALSE, allres=FALSE,sparse=FALSE,sparseStop=TRUE,...) {
if(scaleX){Xplan <- scale(Xplan)}
if((scaleY & missing(time2))){time <- scale(time)}
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("plsRglm"),silent=TRUE)
on.exit(try(unloadNamespace("plsRglm"),silent=TRUE))

mf <- match.call(expand.dots = FALSE)
m <- match(c("time", "time2", "event", "type", "origin"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("Surv")
YCsurv <- eval(mf, parent.frame())

mf2 <- match.call(expand.dots = FALSE)
m2 <- match(c("nt","typeVC","sparse","sparseStop"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$dataY <- time
mf2$dataX <- Xplan
mf2[[1L]] <- as.name("PLS_lm")
pls3_mod <- eval(mf2, parent.frame())
tt_pls3 <- data.frame(pls3_mod$tt)

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_pls3
mf2b[[1L]] <- as.name("coxph")
cox_pls3 <- eval(mf2b, parent.frame())
cox_pls3$call$data <- as.name("tt_pls3")

if(!allres){return(cox_pls3)}
else {return(list(tt_pls3=tt_pls3, cox_pls3=cox_pls3, pls3_mod=pls3_mod))}
}
