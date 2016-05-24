coxpls.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, ncomp=min(7,ncol(Xplan)), modepls="regression", plot=FALSE, allres=FALSE,...) {
if(scaleX){Xplan <- scale(Xplan); XplanScal <- attr(Xplan,"scaled:scale"); XplanCent <- attr(Xplan,"scaled:center"); Xplan <- as.data.frame(Xplan)} else {Xplan <- as.data.frame(Xplan);XplanScal <- rep(1,ncol(Xplan)); XplanCent <- rep(0,ncol(Xplan))}
if((scaleY & missing(time2))){time <- scale(time)}
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("mixOmics"),silent=TRUE)
on.exit(try(unloadNamespace("mixOmics"),silent=TRUE),add=TRUE)



mf <- match.call(expand.dots = FALSE)
m <- match(c("time", "time2", "event", "type", "origin"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("Surv")
YCsurv <- eval(mf, parent.frame())

mf2 <- match.call(expand.dots = FALSE)
m2 <- match(c("ncomp"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$ncomp <- eval.parent(mf2$ncomp)
mf2$X <- eval.parent(Xplan)
mf2$Y <- eval.parent(time)
mf2$mode <- eval.parent(modepls)
mf2$scale.X = FALSE
mf2$scale.Y = FALSE
mf2[[1L]] <- as.name("pls.cox")
if(mf2$ncomp==0){
pls_mod <- NULL
} else {
pls_mod <- eval(mf2)
}
tt_pls <- data.frame(pls_mod$variates$X)
if(mf2$ncomp>0){
colnames(tt_pls) <- paste("dim",1:ncol(tt_pls),sep=".")
}

if(mf2$ncomp==0){
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~1)
mf2b$data <- tt_pls
mf2b[[1L]] <- as.name("coxph")
cox_pls <- eval(mf2b, parent.frame())
cox_pls$call$data <- as.name("tt_pls")
} else {
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_pls
mf2b[[1L]] <- as.name("coxph")
cox_pls <- eval(mf2b, parent.frame())
cox_pls$call$data <- as.name("tt_pls")
}

if(!allres){return(cox_pls)}
else {
CoeffCFull = matrix(NA,nrow=ncomp,ncol=ncomp)
if(mf2$ncomp>0){
for(iii in 1:ncomp)
{
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_pls[,1:iii,drop=FALSE]
mf2b[[1L]] <- as.name("coxph")
cox_pls <- eval(mf2b, parent.frame())
cox_pls$call$data <- as.name("tt_pls")
CoeffCFull[,iii] <- c(cox_pls$coefficients,rep(NA,ncomp-iii))
}
}
return(list(tt_pls=tt_pls, cox_pls=cox_pls, pls_mod=pls_mod, XplanScal=XplanScal, XplanCent=XplanCent, CoeffCFull=CoeffCFull))}
}
