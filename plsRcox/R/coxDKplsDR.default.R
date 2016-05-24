coxDKplsDR.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, ncomp=min(7,ncol(Xplan)), modepls="regression", plot=FALSE, allres=FALSE, kernel="rbfdot", hyperkernel, verbose=TRUE,...) {
if(scaleX){Xplan <- scale(Xplan); XplanScal <- attr(Xplan,"scaled:scale"); XplanCent <- attr(Xplan,"scaled:center"); Xplan <- as.data.frame(Xplan)} else {Xplan <- as.data.frame(Xplan);XplanScal <- rep(1,ncol(Xplan)); XplanCent <- rep(0,ncol(Xplan))}
if((scaleY & missing(time2))){time <- scale(time)}
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("mixOmics"),silent=TRUE)
on.exit(try(unloadNamespace("mixOmics"),silent=TRUE),add=TRUE)
try(attachNamespace("kernlab"),silent=TRUE)
on.exit(try(unloadNamespace("kernlab"),silent=TRUE),add=TRUE)


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

if(verbose){cat("Kernel : ",kernel,"\n")}
kernel2c <- get(kernel)
if(missing(hyperkernel)){if(kernel=="rbfdot"){
mf2c <- match.call(expand.dots = FALSE)
m2c <- match(NULL, names(mf2c), 0L)
mf2c <- mf2c[c(1L, m2c)]
mf2c$x <- as.matrix(Xplan)
mf2c$scaled <- FALSE
mf2c[[1L]] <- as.name("sigest")
srangeDKplsDR_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKplsDR_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKplsDR_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}
if(kernel=="laplacedot"){
mf2c <- match.call(expand.dots = FALSE)
m2c <- match(NULL, names(mf2c), 0L)
mf2c <- mf2c[c(1L, m2c)]
mf2c$x <- as.matrix(Xplan)
mf2c$scaled <- FALSE
mf2c[[1L]] <- as.name("sigest")
srangeDKplsDR_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKplsDR_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKplsDR_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}} else {formals(kernel2c) <- hyperkernel
if(verbose){if(kernel=="rbfdot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}
if(verbose){if(kernel=="laplacedot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}}
kernDKplsDR_mod <- eval(call(as.character(quote(kernel2c))))
Xplan_kernDKplsDR_mod <- kernelMatrix(kernDKplsDR_mod, as.matrix(Xplan))

mf3 <- match.call(expand.dots = FALSE)
m3 <- match(c("ncomp"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3$ncomp <- eval.parent(mf3$ncomp)
mf3$X <- eval.parent(Xplan_kernDKplsDR_mod)
mf3$Y <- eval.parent(DR_coxph)
mf3$mode<- eval.parent(modepls)
mf3$scale.X = FALSE
mf3$scale.Y = FALSE
mf3[[1L]] <- as.name("pls.cox")
if(mf3$ncomp==0){
DKplsDR_mod <- NULL
} else {
  DKplsDR_mod <- eval(mf3)
}
tt_DKplsDR <- data.frame(DKplsDR_mod$variates$X)
if(mf3$ncomp>0){
colnames(tt_DKplsDR) <- paste("dim",1:ncol(tt_DKplsDR),sep=".")
}

if(mf3$ncomp==0){
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~1)
mf2b$data <- tt_DKplsDR
mf2b[[1L]] <- as.name("coxph")                                     
cox_DKplsDR <- eval(mf2b, parent.frame())
cox_DKplsDR$call$data <- as.name("tt_DKplsDR")
} else {
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_DKplsDR
mf2b[[1L]] <- as.name("coxph")                                     
cox_DKplsDR <- eval(mf2b, parent.frame())
cox_DKplsDR$call$data <- as.name("tt_DKplsDR")
}

if(!allres){return(cox_DKplsDR)}
else {
CoeffCFull = matrix(NA,nrow=ncomp,ncol=ncomp)
if(mf3$ncomp>0){
for(iii in 1:ncomp)
{
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_DKplsDR[,1:iii,drop=FALSE]
mf2b[[1L]] <- as.name("coxph")
cox_DKplsDR <- eval(mf2b, parent.frame())
cox_DKplsDR$call$data <- as.name("tt_DKplsDR")
CoeffCFull[,iii] <- c(cox_DKplsDR$coefficients,rep(NA,ncomp-iii))
}
}
return(list(tt_DKplsDR=tt_DKplsDR, cox_DKplsDR=cox_DKplsDR, DKplsDR_mod=DKplsDR_mod, XplanScal=XplanScal, XplanCent=XplanCent, CoeffCFull=CoeffCFull, kernDKplsDR_mod=kernDKplsDR_mod))}
}
