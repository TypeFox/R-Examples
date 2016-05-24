coxDKsplsDR.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, ncomp=min(7,ncol(Xplan)), modepls="regression", plot=FALSE, allres=FALSE, eta, trace=FALSE,kernel="rbfdot",hyperkernel, verbose=TRUE,...) {
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
srangeDKsplsDR_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKsplsDR_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKsplsDR_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}
if(kernel=="laplacedot"){
mf2c <- match.call(expand.dots = FALSE)
m2c <- match(NULL, names(mf2c), 0L)
mf2c <- mf2c[c(1L, m2c)]
mf2c$x <- as.matrix(Xplan)
mf2c$scaled <- FALSE
mf2c[[1L]] <- as.name("sigest")
srangeDKsplsDR_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKsplsDR_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKsplsDR_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}} else {formals(kernel2c) <- hyperkernel
if(verbose){if(kernel=="rbfdot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}
if(verbose){if(kernel=="laplacedot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}}
kernDKsplsDR_mod <- eval(call(as.character(quote(kernel2c))))
Xplan_kernDKsplsDR_mod <- kernelMatrix(kernDKsplsDR_mod, as.matrix(Xplan))

mf3 <- match.call(expand.dots = FALSE)
m3 <- match(c("ncomp"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3$x <- eval.parent(Xplan_kernDKsplsDR_mod)
mf3$y <- eval.parent(DR_coxph)
mf3$fit <- eval.parent(modepls)
mf3$K <- eval.parent(mf3$ncomp)
mf3$scale.x = FALSE
mf3$ncomp <- NULL
mf3$eta <- eval.parent(eta)
mf3$trace <- eval.parent(trace)
mf3[[1L]] <- as.name("spls.cox")
if(mf3$K==0){
DKsplsDR_mod <- NULL
} else {
DKsplsDR_mod <- eval(mf3)
}
DKsplsDR_modplsr <- DKsplsDR_mod$plsmod
DKsplsDR_mod$plsmod <- NULL
tt_DKsplsDR <- data.frame(DKsplsDR_modplsr$variates$X)
if(mf3$K>0){
colnames(tt_DKsplsDR) <- paste("dim",1:ncol(tt_DKsplsDR),sep=".")
}

if(mf3$K==0){
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~1)
mf2b$data <- tt_DKsplsDR
mf2b[[1L]] <- as.name("coxph")
cox_DKsplsDR <- eval(mf2b, parent.frame())
cox_DKsplsDR$call$data <- as.name("tt_DKsplsDR")
} else {
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_DKsplsDR
mf2b[[1L]] <- as.name("coxph")
cox_DKsplsDR <- eval(mf2b, parent.frame())
cox_DKsplsDR$call$data <- as.name("tt_DKsplsDR")
}

if(!allres){return(cox_DKsplsDR)}
else {
CoeffCFull = matrix(NA,nrow=ncomp,ncol=ncomp)
if(mf3$K>0){
for(iii in 1:ncomp)
{
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_DKsplsDR[,1:iii,drop=FALSE]
mf2b[[1L]] <- as.name("coxph")
cox_DKsplsDR <- eval(mf2b, parent.frame())
cox_DKsplsDR$call$data <- as.name("tt_DKsplsDR")
CoeffCFull[,iii] <- c(cox_DKsplsDR$coefficients,rep(NA,ncomp-iii))
}
}
return(list(tt_DKsplsDR=tt_DKsplsDR, cox_DKsplsDR=cox_DKsplsDR, DKsplsDR_mod=DKsplsDR_mod, DKsplsDR_modplsr=DKsplsDR_modplsr, XplanScal=XplanScal, XplanCent=XplanCent, CoeffCFull=CoeffCFull, kernDKsplsDR_mod=kernDKsplsDR_mod))}
}



