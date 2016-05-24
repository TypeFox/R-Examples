DKplsRcoxmodel.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, nt=min(2,ncol(Xplan)),limQ2set=.0975, dataPredictY=Xplan, pvals.expli=FALSE,alpha.pvals.expli=.05,tol_Xi=10^(-12),weights,control, sparse=FALSE,sparseStop=TRUE, plot=FALSE, allres=FALSE, kernel="rbfdot", hyperkernel, verbose=TRUE,...) {

dataX<-Xplan
dataY<-time

try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("kernlab"),silent=TRUE)
on.exit(try(unloadNamespace("kernlab"),silent=TRUE),add=TRUE)

if(missing(weights)){NoWeights=TRUE} else {if(all(weights==rep(1,length(dataY)))){NoWeights=TRUE} else {NoWeights=FALSE}}

if ((scaleY & missing(time2))) {if(NoWeights){RepY <- scale(dataY)} else {meanY <- weighted.mean(dataY,weights); stdevY <- sqrt((length(dataY)-1)/length(dataY)*weighted.mean((dataY-meanY)^2,weights)); RepY <- (dataY-meanY)/stdevY; attr(RepY,"scaled:center") <- meanY ; attr(RepY,"scaled:scale") <- stdevY}}
else {
    RepY <- time
    attr(RepY,"scaled:center") <- 0
    attr(RepY,"scaled:scale") <- 1
}
if (scaleX) {if(NoWeights){ExpliX <- scale(dataX)} else {meanX <- apply(dataX,2,weighted.mean,weights); stdevX <- sqrt((length(time)-1)/length(time)*apply((sweep(dataX,2,meanX))^2,2,weighted.mean,weights)); ExpliX <- sweep(sweep(dataX, 2, meanX), 2 ,stdevX, "/"); attr(ExpliX,"scaled:center") <- meanX ; attr(ExpliX,"scaled:scale") <- stdevX}
}
else {
    ExpliX <- dataX
    attr(ExpliX,"scaled:center") <- rep(0,ncol(dataX))
    attr(ExpliX,"scaled:scale") <- rep(1,ncol(dataX))
}

mf <- match.call(expand.dots = FALSE)
m <- match(c("time", "time2", "event", "type", "origin"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("Surv")
mf$time <- RepY
YCsurv <- eval(mf, parent.frame())
attr(YCsurv,"scaled:center") <- attr(RepY,"scaled:center")
attr(YCsurv,"scaled:scale") <- attr(RepY,"scaled:scale")
RepY <- YCsurv

if(verbose){cat("Kernel : ",kernel,"\n")}
kernel2c <- get(kernel)
if(missing(hyperkernel)){if(kernel=="rbfdot"){
mf2c <- match.call(expand.dots = FALSE)
m2c <- match(NULL, names(mf2c), 0L)
mf2c <- mf2c[c(1L, m2c)]
mf2c$x <- ExpliX
mf2c$scaled <- FALSE
mf2c[[1L]] <- as.name("sigest")
srangeDKplsRcox_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKplsRcox_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKplsRcox_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}
if(kernel=="laplacedot"){
mf2c <- match.call(expand.dots = FALSE)
m2c <- match(NULL, names(mf2c), 0L)
mf2c <- mf2c[c(1L, m2c)]
mf2c$x <- ExpliX
mf2c$scaled <- FALSE
mf2c[[1L]] <- as.name("sigest")
srangeDKplsRcox_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKplsRcox_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKplsRcox_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}} else {formals(kernel2c) <- hyperkernel
if(verbose){if(kernel=="rbfdot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}
if(verbose){if(kernel=="laplacedot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}}
kernDKplsRcox_mod <- eval(call(as.character(quote(kernel2c))))
Xplan_kernDKplsRcox_mod <- kernelMatrix(kernDKplsRcox_mod, as.matrix(ExpliX))

mf3 <- match.call(expand.dots = FALSE)
m3 <- match(c("time","time2","event","type","origin","collapse","weighted","nt","limQ2set","dataPredictY","scaleY","pvals.expli","alpha.pvals.expli","tol_Xi","weights","subset","control","sparse","sparseStop"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3$Xplan <- as.formula(~.)
mf3$dataXplan <- Xplan_kernDKplsRcox_mod
mf3$model_frame <- FALSE
mf3$scaleX <- FALSE
mf3[[1L]] <- as.name("plsRcox")
DKplsRcox_mod <- eval(mf3, parent.frame())
tt_DKplsRcox <- data.frame(DKplsRcox_mod$tt)

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_DKplsRcox
mf2b[[1L]] <- as.name("coxph")                                     
cox_DKplsRcox <- eval(mf2b, parent.frame())
cox_DKplsRcox$call$data <- as.name("tt_DKplsRcox")


if(!allres){return(cox_DKplsRcox)}
else {return(list(tt_DKplsRcox=tt_DKplsRcox, cox_DKplsRcox=cox_DKplsRcox, DKplsRcox_mod=DKplsRcox_mod, kernDKplsRcox_mod=kernDKplsRcox_mod, XplanScal=attr(ExpliX,"scaled:scale"), XplanCent= attr(ExpliX,"scaled:center") ))}
}
