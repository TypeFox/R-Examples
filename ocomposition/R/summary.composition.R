summary.composition <- function(object, individual = NULL, ...){

print.r <- function(x){
cat(" ", '\n')
cat(attr(x, "ext.name"), '\n')
x <- as.mcmc(x)
ci <- HPDinterval(x)
R <- cbind(apply(x,2,mean), apply(x,2,sd), ci, batchSE(x,nrow(x)/10), heidel.diag(x)[,4])
R <- round(R, 3)
rownames(R) <- colnames(x)
colnames(R) <- c("Estimate", "SE", "0.025-HPD", "0.975-HPD","MCMC error", "Heidel-conv")
print(R)
}

if (!is.null(individual)) {
for (i in individual){ 
y <- t(object$g[,i,])
attr(y, "ext.name") <- paste("Component size model parameter gamma_", i, sep="")
R = print.r(y)
return(invisible(R))
}
}

if(is.null(individual)){
R1 = print.r(object$b)
R2 = print.r(object$mu)
return(invisible(list(R1, R2)))
}

}
