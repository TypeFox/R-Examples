apval_Chen2010_diffcov <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	p <- dim(sam1)[2]
	T1 <- sam1 %*% t(sam1)
	T2 <- sam2 %*% t(sam2)
	P1 <- (sum(T1) - sum(diag(T1)))/(n1*(n1 - 1))
	P2 <- (sum(T2) - sum(diag(T2)))/(n2*(n2 - 1))
	P3 <- -2*sum(sam1 %*% t(sam2))/(n1*n2)
	T <- P1 + P2 + P3

	tr.cov1.sq <- tr.cov2.sq <- tr.cov1.cov2 <- 0
	for(j in 1:n1){
		for(k in 1:n1){
			if(j != k){
				tempmean <- (colSums(sam1) - sam1[j,] - sam1[k,])/(n1 - 2)
				P1 <- sum(sam1[j,]*(sam1[k,] - tempmean))
				P2 <- sum(sam1[k,]*(sam1[j,] - tempmean))
				tr.cov1.sq <- tr.cov1.sq + P1*P2
			}
		}
	}
	tr.cov1.sq <- tr.cov1.sq/(n1*(n1 - 1))
	for(j in 1:n2){
		for(k in 1:n2){
			if(j != k){
				tempmean <- (colSums(sam2) - sam2[j,] - sam2[k,])/(n2 - 2)
				P1 <- sum(sam2[j,]*(sam2[k,] - tempmean))
				P2 <- sum(sam2[k,]*(sam2[j,] - tempmean))
				tr.cov2.sq <- tr.cov2.sq + P1*P2
			}
		}
	}
	tr.cov2.sq <- tr.cov2.sq/(n2*(n2 - 1))
	for(j in 1:n1){
		for(k in 1:n2){
			tempmean1 <- (colSums(sam1) - sam1[j,])/(n1 - 1)
			tempmean2 <- (colSums(sam2) - sam2[k,])/(n2 - 1)
			P1 <- sum(sam1[j,]*(sam2[k,] - tempmean2))
			P2 <- sum(sam2[k,]*(sam1[j,] - tempmean1))
			tr.cov1.cov2 <- tr.cov1.cov2 + P1*P2
		}
	}
	tr.cov1.cov2 <- tr.cov1.cov2/(n1*n2)

	test.stat <- T/sqrt(2/(n1*(n1 - 1))*tr.cov1.sq + 2/(n2*(n2 - 1))*tr.cov2.sq + 4/(n1*n2)*tr.cov1.cov2)
	test.stat <- as.numeric(test.stat)
	pval <- 1 - pnorm(test.stat)
	names(pval) <- "Chen2010"
	out <- NULL
	out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
	out$cov.assumption <- "the two groups have different covariances"
	out$method <- "asymptotic distribution"
	out$pval <- pval
	return(out)
}