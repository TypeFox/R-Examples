larsDR_coxph.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=FALSE, scaleY=TRUE, plot=FALSE, typelars="lasso", normalize = TRUE, max.steps, use.Gram = TRUE, allres=FALSE, verbose=TRUE,...){
try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
try(attachNamespace("lars"),silent=TRUE)
on.exit(try(unloadNamespace("lars"),silent=TRUE),add=TRUE)

if(!is.matrix(Xplan)){Xplan <- as.matrix(Xplan);if(verbose){cat("scales\n")}}
if(scaleX){Xplan <- scale(Xplan)}
if((scaleY & missing(time2))){time <- scale(time)}


mf <- match.call(expand.dots = FALSE)
m <- match(c("time", "time2", "event", "type", "origin"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("Surv")
YCsurv <- eval(mf, parent.frame())

if(plot){plot(survival::survfit(YCsurv~1))}

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
m3 <- match(c("normalize", "use.Gram", "max.steps"), names(mf3), 0L)
mf3 <- mf3[c(1L, m3)]
mf3[[1L]] <- as.name("lars")
mf3$x <- Xplan
mf3$y <- DR_coxph
mf3$type <- typelars
larsDR <- eval(mf3, parent.frame())
X_larsDR <- data.frame(Xplan[,as.numeric(names(table(unlist(larsDR$actions))[as.numeric(names(table(unlist(larsDR$actions))))>0]))])

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~.)
mf2b$data <- X_larsDR
mf2b[[1L]] <- as.name("coxph")
cox_larsDR <- eval(mf2b, parent.frame())
cox_larsDR$call$data <- as.name("X_larsDR")


if(!allres){return(cox_larsDR)}
else {return(list(DR_coxph=DR_coxph, larsDR=larsDR, X_larsDR=X_larsDR, cox_larsDR=cox_larsDR))}
}
