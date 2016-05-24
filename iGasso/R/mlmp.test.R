mlmp.test = function(g, y, weights=NULL, stat="score")
{
	DNAME = paste(deparse(substitute(g)), "and", deparse(substitute(y)))
	y = scale(y)    # dimension n x q
	n = nrow(y)
	q = ncol(y)
	A = crossprod(y)

	vg = apply(g, 2, var)    
	g = as.matrix(g[,is.finite(1/vg)], nrow=n)                        # remove non-informative SNPs
	if (is.null(weights)) weights = 1/vg[is.finite(1/vg)]     # Not the square root, different from SKAT
	else weights = weights[is.finite(1/vg)]
	g = g %*% diag(sqrt(weights))
	Sigma0 = var(g, g)
	p = ncol(g)

    B.t = solve(A, crossprod(y, g))
    H = crossprod(B.t, A %*% B.t)
	Sigma = ((n-1)*Sigma0 - H)/(n-1-q)     # the unbiased estimate of Sigma

    if (stat=="F"){
    	S = sum(diag(H))/q/sum(diag(Sigma))
    	names(S) = "F stat"
    	lmbd = eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values
    	p.value = davies(0, lambda=c(lmbd, -S*q/(n-1-q)*lmbd), h=rep(c(q, n-1-q), c(p, p)))$Qq
    }
    else if (stat=="Wald"){
    	S = sum(diag(H))/q/sum(diag(Sigma))
    	names(S) = "Wald stat"
    	lmbd = eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values
    	p.value = davies(S*q*sum(diag(Sigma)), lambda=lmbd, h=rep(q, p))$Qq
    }
    else if (stat=="score"){
    	S = sum(diag(H))/q/sum(diag(Sigma0))
    	names(S) = "Score stat"
    	lmbd = eigen(Sigma0, symmetric=TRUE, only.values=TRUE)$values
    	p.value = davies(S*q*sum(diag(Sigma0)), lambda=lmbd, h=rep(q, p))$Qq
    }

    PAR = c(p, q)
    names(PAR) = c("# loci", "# traits")
    structure(list(statistic = S, p.value = p.value, parameter = PAR, 
        method = "Multilocus multiple phenotypes test of association", data.name = DNAME), 
        class = "htest")
}        
