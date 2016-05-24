coxpls3DR.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, nt=min(7,ncol(Xplan)),typeVC="none", plot=FALSE, allres=FALSE,sparse=FALSE,sparseStop=TRUE,...) {
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

mf1 <- match.call(expand.dots = TRUE)
m1 <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf1), 0L)
mf1 <- mf1[c(1L, m1)]
mf1$formula <- as.formula(YCsurv~1)
mf1[[1L]] <- as.name("coxph")
cox3DR <- eval(mf1, parent.frame())

mf2 <- match.call(expand.dots = FALSE)
m2 <- match(c("weighted", "collapse", "origin"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$type <- typeres
mf2$object <- cox3DR
mf2[[1L]] <- as.name("residuals")
DR3_coxph <- eval(mf2, parent.frame())

mf3 <- match.call(expand.dots = FALSE)
m3 <- match(c("nt", "typeVC","sparse","sparseStop"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3$dataY <- DR3_coxph
mf3$dataX <- Xplan
mf3[[1L]] <- as.name("PLS_lm")
pls3DR_mod <- eval(mf3, parent.frame())
tt_pls3DR <- data.frame(pls3DR_mod$tt)

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_pls3DR
mf2b[[1L]] <- as.name("coxph")
cox_pls3DR <- eval(mf2b, parent.frame())
cox_pls3DR$call$data <- as.name("tt_pls3DR")


if(!allres){return(cox_pls3DR)}
else {return(list(tt_pls3DR=tt_pls3DR, cox_pls3DR=cox_pls3DR, pls3DR_mod=pls3DR_mod))}
}
