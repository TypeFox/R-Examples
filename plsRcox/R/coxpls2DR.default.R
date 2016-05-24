coxpls2DR.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, ncomp=min(7,ncol(Xplan)), methodpls="kernelpls", validation = "CV", plot=FALSE, allres=FALSE,...) {
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

mf1 <- match.call(expand.dots = TRUE)
m1 <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf1), 0L)
mf1 <- mf1[c(1L, m1)]
mf1$formula <- as.formula(YCsurv~1)
mf1[[1L]] <- as.name("coxph")
coxDR <- eval(mf1, parent.frame())

mf2 <- match.call(expand.dots = FALSE)
m2 <- match(c("weighted", "collapse", "origin"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$type <- typeres
mf2$object <- coxDR
mf2[[1L]] <- as.name("residuals")
DR_coxph <- eval(mf2, parent.frame())

mf3 <- match.call(expand.dots = FALSE)
m3 <- match(c("ncomp", "validation"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3$formula <- as.formula(DR_coxph~.)
mf3$data <- Xplan
mf3$method <- methodpls
mf3[[1L]] <- as.name("plsr")
pls2DR_mod <- eval(mf3, parent.frame())
tt_pls2DR <- data.frame(scores(pls2DR_mod)[,])

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_pls2DR
mf2b[[1L]] <- as.name("coxph")
cox_pls2DR <- eval(mf2b, parent.frame())
cox_pls2DR$call$data <- as.name("tt_pls2DR")


if(!allres){return(cox_pls2DR)}
else {return(list(tt_pls2DR=tt_pls2DR, cox_pls2DR=cox_pls2DR, pls2DR_mod=pls2DR_mod))}
}
