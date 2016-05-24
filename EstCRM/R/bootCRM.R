bootCRM <- function(data,max.item,min.item,max.EMCycle=500,converge=.01,type="Shojima",BFGS=TRUE,nsample=50) {

org <- EstCRMitem(data, max.item, min.item, max.EMCycle, converge,type="Shojima",BFGS=TRUE)$param

pars <- vector("list",nsample)
for(i in 1:nsample) {
	d <- data[sample(1:nrow(data),nrow(data),replace = TRUE),]
	pars[[i]] <- EstCRMitem(d, max.item, min.item, max.EMCycle, converge,type="Shojima",BFGS=TRUE)$param
}

average <- matrix(nrow=nrow(org),ncol=ncol(org))
sds <- matrix(nrow=nrow(org),ncol=ncol(org))

for(i in 1:nrow(org)) {
for(j in 1:ncol(org)) { 
	estimate <- c()
	for(k in 1:nsample) { 
		estimate[k] <- pars[[k]][i,j] 
	}
		average[i,j]=mean(estimate,na.rm=TRUE)
		sds[i,j]=sd(estimate,na.rm=TRUE)
}}

disc <- as.data.frame(matrix(nrow=nrow(org), ncol=3))
disc[,1] <- org[,1]
disc[,2] <- average[,1]
disc[,3] <- sds[,1]
colnames(disc) <- c("Original Estimate","Bootstrap Estimate","Bootstrap Std.Err.")
rownames(disc) <- colnames(data)


diff <- as.data.frame(matrix(nrow=nrow(org), ncol=3))
diff[,1] <- org[,2]
diff[,2] <- average[,2]
diff[,3] <- sds[,2]
colnames(diff) <- c("Original Estimate","Bootstrap Estimate","Bootstrap Std.Err.")
rownames(diff) <- colnames(data)

alpha <- as.data.frame(matrix(nrow=nrow(org), ncol=3))
alpha[,1] <- org[,3]
alpha[,2] <- average[,3]
alpha[,3] <- sds[,3]
colnames(alpha) <- c("Original Estimate","Bootstrap Estimate","Bootstrap Std.Err.")
rownames(alpha) <- colnames(data)

return(list(Discrimination=disc,Difficulty=diff,Alpha=alpha))
}







