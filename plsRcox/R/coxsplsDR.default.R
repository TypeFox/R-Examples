coxsplsDR.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, ncomp=min(7,ncol(Xplan)), modepls="regression", plot=FALSE, allres=FALSE, eta=.5, trace=FALSE,...) {
if(scaleX){Xplan <- scale(Xplan); XplanScal <- attr(Xplan,"scaled:scale"); XplanCent <- attr(Xplan,"scaled:center"); Xplan <- as.data.frame(Xplan)} else {Xplan <- as.data.frame(Xplan);XplanScal <- rep(1,ncol(Xplan)); XplanCent <- rep(0,ncol(Xplan))}
if((scaleY & missing(time2))){time <- scale(time)}
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("mixOmics"),silent=TRUE)
on.exit(try(unloadNamespace("mixOmics"),silent=TRUE),add=TRUE)
suppressMessages(try(attachNamespace("spls"),silent=TRUE))
on.exit(try(unloadNamespace("spls"),silent=TRUE),add=TRUE)


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
m3 <- match(c("ncomp"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3$x <- eval.parent(Xplan)
mf3$y <- eval.parent(DR_coxph)
mf3$fit <- eval.parent(modepls)
mf3$K <- eval.parent(mf3$ncomp)
mf3$scale.x = FALSE
mf3$ncomp <- NULL
mf3$eta <- eval.parent(eta)
mf3$trace <- eval.parent(trace)
mf3[[1L]] <- as.name("spls.cox")
if(mf3$K==0){
splsDR_mod <- NULL
} else {
splsDR_mod <- eval(mf3)
}
splsDR_modplsr <- splsDR_mod$plsmod
splsDR_mod$plsmod <- NULL
tt_splsDR = data.frame(splsDR_modplsr$variates$X)
if(mf3$K>0){
colnames(tt_splsDR) <- paste("dim",1:ncol(tt_splsDR),sep=".")
}

if(mf3$K==0){
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~1)
mf2b$data <- tt_splsDR
mf2b[[1L]] <- as.name("coxph")
cox_splsDR <- eval(mf2b, parent.frame())
cox_splsDR$call$data <- as.name("tt_splsDR")
} else {
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_splsDR
mf2b[[1L]] <- as.name("coxph")
cox_splsDR <- eval(mf2b, parent.frame())
cox_splsDR$call$data <- as.name("tt_splsDR")
}

if(!allres){return(cox_splsDR)}
else {
CoeffCFull = matrix(NA,nrow=ncomp,ncol=ncomp)
if(mf3$K>0){
for(iii in 1:ncomp)
{
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_splsDR[,1:iii,drop=FALSE]
mf2b[[1L]] <- as.name("coxph")
cox_splsDR <- eval(mf2b, parent.frame())
cox_splsDR$call$data <- as.name("tt_splsDR")
CoeffCFull[,iii] <- c(cox_splsDR$coefficients,rep(NA,ncomp-iii))
}
}
return(list(tt_splsDR=tt_splsDR, cox_splsDR=cox_splsDR, splsDR_mod=splsDR_mod, splsDR_modplsr=splsDR_modplsr, XplanScal=XplanScal, XplanCent=XplanCent, CoeffCFull=CoeffCFull))}
}



