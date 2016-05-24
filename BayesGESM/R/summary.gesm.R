summary.gesm <-
function(object, ...){
chains <- object$chains     
p <- object$p
q <- object$q
ks <- object$ks
ks2 <- object$ks2
nu0 <- object$nu0


quant005 <- function(x){quantile(x,prob=0.025)}
quant095 <- function(x){quantile(x,prob=0.975)}

reem <- function(aa,b){
                ag <- aa
                ag <- sub("(", "", ag,fixed=TRUE)
                ag <- sub(")", "", ag,fixed=TRUE)
                ag <- sub(b, "", ag,fixed=TRUE)
                ag <- strsplit(ag, ",")
                ag <- ag[[1]][1]
                ag
        }


if(object$family=="Normal" || object$family=="Laplace"  || attr(object$eta,"know")==0){
cat("\n       Error distribution:", object$family)
}else{
cat("\n       Error distribution:", object$family,"(",object$eta,")")
}

cat("\n              Sample size:", length(object$y))
cat("\n Size of posterior sample:", object$post.sam.s, "\n")

cat("\n =================   Location Submodel   =================")

if(p >0){
cat("\n =========   Parametric part     \n")
a <- round(apply(as.matrix(chains[,1:p]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,1:p]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,1:p]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,1:p]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,1:p]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$X)
printCoefmat(e)
}
if(sum(ks) >0){
cat("=========   Nonparametric part   \n")
nvar <- matrix(0,length(ks), 1)
nps <- colnames(object$nps)
for(i in 1:length(ks)){
nvar[i] <- reem(nps[i],"bsp")}

cat(" Effects    internal knots      ")
b <- (object$ks - 3)
rownames(b) <- paste("  ",nvar)
colnames(b) <- "             "
printCoefmat(b)
}

cat("\n =================  Dispersion Submodel  =================")
if(q >0){
cat("\n =========   Parametric part    \n")
if(sum(ks)==0) length(ks) <- 0 
a <- round(apply(as.matrix(chains[ ,(p+sum(ks)+length(ks)+1):(p+sum(ks)+length(ks)+q)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+1):(p+sum(ks)+length(ks)+q)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+1):(p+sum(ks)+length(ks)+q)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+1):(p+sum(ks)+length(ks)+q)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+1):(p+sum(ks)+length(ks)+q)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$Z)
printCoefmat(e)
}

if(sum(ks2) >0){
cat("=========   Nonparametric part   \n")
nvar2 <- matrix(0,length(ks2), 1)
nps2 <- colnames(object$nps2)
for(i in 1:length(ks2)){
nvar2[i] <- reem(nps2[i],"bsp")}

cat(" Effects    internal knots      ")
b <- (object$ks2 - 3)
rownames(b) <- paste(" ",nvar2)
colnames(b) <- "             "
printCoefmat(b)
}

if(object$family!="Normal" && object$family!="Laplace"  && attr(object$eta,"know")==0){
cat("\n =========   Eta parameter    \n")
if(sum(ks)==0) length(ks) <- 0
if(sum(ks2)==0) length(ks2) <- 0
a <-  round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+1):(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+length(nu0))]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+1):(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+length(nu0))]), 2, median), digits=4)
b <-  round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+1):(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+length(nu0))]), 2, sd), digits=4)
c <-  round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+1):(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+length(nu0))]), 2, quant005), digits=4)
d <-  round(apply(as.matrix(chains[,(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+1):(p+sum(ks)+length(ks)+q+sum(ks2)+length(ks2)+length(nu0))]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- cbind(t(paste(" eta",1:length(nu0),"   ")))
printCoefmat(e)
}


if(sum(ks)>0 || sum(ks2)>0){
cat("\n Graphs of the nonparametric effects are provided by \n")
cat(" using the function 'bsp.graph.gesm'    \n")}

cat("\n ============   Model Selection Criteria   ============")
cat("\n DIC=" , object$DIC,   "    LMPL=" , object$LMPL, "\n")

}
