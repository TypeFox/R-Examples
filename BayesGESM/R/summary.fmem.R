summary.fmem <-
function(object, ...){
chains <- object$chains     
p <- object$p
q <- object$q
ks <- object$ks
nu0 <- object$nu0
homo <- object$homo


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


cat("\n ===============    Parametric part   ===============  ")
if(p >0){
cat("\n =======   Covariates measured without error     \n")
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

cat("\n =======   Covariates measured with error     \n")
a <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$M)
printCoefmat(e)

if(sum(ks) >0){
cat("\n ===============   Nonparametric part ===============  \n")

nvar <- matrix(0,length(ks),1)
for(i in 1:length(ks)){
nvar[i] <- reem(colnames(object$nps)[i],"bsp")}

cat(" Effects    internal knots      ")
b <- (object$ks - 3)
rownames(b) <- nvar
colnames(b) <- "             "
printCoefmat(b)

cat("\n Graphs of the nonparametric effects are provided by \n")
cat(" using the function 'bsp.graph.fmem'    \n")
}

if(object$family!="Normal" && object$family!="Laplace"  && attr(object$eta,"know")==0){
cat("\n =======   Eta parameter    \n")
if(homo==1){
if(sum(ks)==0) length(ks) <- 0
a <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+1+length(ks)+1):(p+3*q+sum(ks)+1+length(ks)+length(nu0))]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+3*q+sum(ks)+1+length(ks)+1):(p+3*q+sum(ks)+1+length(ks)+length(nu0))]), 2, median), digits=4)
b <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+1+length(ks)+1):(p+3*q+sum(ks)+1+length(ks)+length(nu0))]), 2, sd), digits=4)
c <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+1+length(ks)+1):(p+3*q+sum(ks)+1+length(ks)+length(nu0))]), 2, quant005), digits=4)
d <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+1+length(ks)+1):(p+3*q+sum(ks)+1+length(ks)+length(nu0))]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$nu0)
printCoefmat(e)
}else{
if(sum(ks)==0) length(ks) <- 0
a <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+length(ks)+1):(p+3*q+sum(ks)+length(ks)+length(nu0))]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+3*q+sum(ks)+length(ks)+1):(p+3*q+sum(ks)+length(ks)+length(nu0))]), 2, median), digits=4)
b <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+length(ks)+1):(p+3*q+sum(ks)+length(ks)+length(nu0))]), 2, sd), digits=4)
c <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+length(ks)+1):(p+3*q+sum(ks)+length(ks)+length(nu0))]), 2, quant005), digits=4)
d <-  round(apply(as.matrix(chains[,(p+3*q+sum(ks)+length(ks)+1):(p+3*q+sum(ks)+length(ks)+length(nu0))]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$nu0)
printCoefmat(e)
}
}

if(homo ==1){
cat("\n =======   Dispersion parameter    \n")
a <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- "Sigma2_y  "
printCoefmat(e)
cat("\n  Ratio of the error variances: " , (object$omeg))
}


cat("\n\n ==========   Model Selection Criteria   ==========")
cat("\n DIC=" , object$DIC,  "   LMPL=" , object$LMPL, "\n")
}
