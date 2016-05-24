coxDKpls2DR.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, ncomp=min(7,ncol(Xplan)), methodpls="kernelpls", validation = "CV", plot=FALSE, allres=FALSE, kernel="rbfdot", hyperkernel, verbose=TRUE,...) {
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("pls"),silent=TRUE)
on.exit(try(unloadNamespace("pls"),silent=TRUE),add=TRUE)
try(attachNamespace("kernlab"),silent=TRUE)
on.exit(try(unloadNamespace("kernlab"),silent=TRUE),add=TRUE)


if(scaleX){Xplan <- as.data.frame(scale(Xplan))} else {Xplan <- as.data.frame(Xplan)}
if((scaleY & missing(time2))){time <- scale(time)}

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
srangeDKpls2DR_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKpls2DR_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKpls2DR_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}
if(kernel=="laplacedot"){
mf2c <- match.call(expand.dots = FALSE)
m2c <- match(NULL, names(mf2c), 0L)
mf2c <- mf2c[c(1L, m2c)]
mf2c$x <- as.matrix(Xplan)
mf2c$scaled <- FALSE
mf2c[[1L]] <- as.name("sigest")
srangeDKpls2DR_mod <- eval(mf2c, parent.frame())
hyperkernel=list(sigma = srangeDKpls2DR_mod[2])
if(verbose){cat("Estimated_sigma ",srangeDKpls2DR_mod[2],"\n")}
formals(kernel2c) <- hyperkernel
}} else {formals(kernel2c) <- hyperkernel
if(verbose){if(kernel=="rbfdot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}
if(verbose){if(kernel=="laplacedot"){cat("Used_sigma ",hyperkernel$sigma,"\n")}}}
kernDKpls2DR_mod <- eval(call(as.character(quote(kernel2c))))
Xplan_kernDKpls2DR_mod <- kernelMatrix(kernDKpls2DR_mod, as.matrix(Xplan))

mf3 <- match.call(expand.dots = FALSE)
m3 <- match(c("ncomp", "validation"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3$formula <- as.formula(DR_coxph~.)
mf3$data <- Xplan_kernDKpls2DR_mod
mf3$method<-methodpls
mf3[[1L]] <- as.name("plsr")
DKpls2DR_mod <- eval(mf3, parent.frame())
tt_DKpls2DR <- data.frame(scores(DKpls2DR_mod)[,])

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- tt_DKpls2DR
mf2b[[1L]] <- as.name("coxph")                                     
cox_DKpls2DR <- eval(mf2b, parent.frame())
cox_DKpls2DR$call$data <- as.name("tt_DKpls2DR")


if(!allres){return(cox_DKpls2DR)}
else {return(list(tt_DKpls2DR=tt_DKpls2DR, cox_DKpls2DR=cox_DKpls2DR, DKpls2DR_mod=DKpls2DR_mod, kernDKpls2DR_mod=kernDKpls2DR_mod))}
}
