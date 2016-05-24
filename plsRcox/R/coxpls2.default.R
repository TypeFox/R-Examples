coxpls2.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, ncomp=min(7,ncol(Xplan)), methodpls="kernelpls", validation = "CV", plot=FALSE, allres=FALSE,...) {
if(scaleX){Xplan <- as.data.frame(scale(Xplan))} else {Xplan <- as.data.frame(Xplan)}
if((scaleY & missing(time2))){time <- scale(time)}
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("pls"),silent=TRUE)
on.exit(try(unloadNamespace("pls"),silent=TRUE),add=TRUE)

mf <- match.call(expand.dots = FALSE)
m <- match(c("time", "time2", "event", "type", "origin"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("Surv")
YCsurv <- eval(mf, parent.frame())

mf2 <- match.call(expand.dots = FALSE)
m2 <- match(c("ncomp", "validation"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$formula <- as.formula(time~.)
mf2$data <- Xplan
mf2$method <- methodpls
mf2[[1L]] <- as.name("plsr")
pls_mod <- eval(mf2, parent.frame())
tt_pls <- data.frame(scores(pls_mod)[,])

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_pls
mf2b[[1L]] <- as.name("coxph")
cox_pls <- eval(mf2b, parent.frame())
cox_pls$call$data <- as.name("tt_pls")

if(!allres){return(cox_pls)}
else {return(list(tt_pls=tt_pls, cox_pls=cox_pls, pls_mod=pls_mod))}
}
