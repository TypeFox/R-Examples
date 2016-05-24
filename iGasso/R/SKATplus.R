SKATplus = function(y, G, X=NULL, out_type="D", tau=NULL, permutation=FALSE, B=1000){
	DNAME = paste(deparse(substitute(G)), "and", deparse(substitute(y)), "and", deparse(substitute(X)))
    n = length(y)
    G = as.matrix(G, nrow = n)      # genotype matrix
	if(is.null(X)) X = matrix(1, ncol=1, nrow=n) else X = cbind(1, X)
    k = ncol(X)                               # k=p+1

    if (out_type == "D"){
     	fit0 = glm(y ~ X, family="binomial")
    	    hpi = fitted.values(fit0)
    	    res.y = y - hpi
    	    V  = hpi*(1-hpi)
    }
    if (out_type == "C"){
       	fit0 = lm(y ~ X)
	    res.y = residuals(fit0)
	    	V = rep(summary(fit0)$sigma^2, n)
    }
    stat = sum(crossprod(G, res.y)^2)

    if (is.null(tau)) selected = (res.y <= 0) else selected = (rank(res.y) <= tau*n)
    G0 = G[selected,]  
    n0 = sum(selected)
    if (permutation){
        stat.p = rep(1,B)
        for (i in 1:B) {
          	res.y0 = sample(res.y, n0)
          	stat.p[i] = sum(crossprod(G0, res.y0 - mean(res.y0))^2)*(n-k)/(n0-k)
        }
        pvalue = mean(stat.p>=stat)
    }
    else{
        
        PPP = function(V, X){
     	    XV = crossprod(X, diag(V))
    	    diag(V) - crossprod(XV, solve(XV %*% X, XV))
        }      

     	P0 = PPP(V[selected], X[selected,])
    	    GPG0 = crossprod(G0, P0 %*% G0)
        pca0 = eigen(GPG0, symmetric=TRUE)
        pvalue = davies(stat, lambda=pca0$values*(n-k)/(n0-k))$Qq
    }

    PAR = c(n, k-1)
    names(PAR) = c("# subjects", "# SNPs")
    structure(list(statistic = stat, p.value = pvalue, parameter = PAR, 
        method = "SKAT+ gene- or pathway-based test of association", data.name = DNAME), 
        class = "htest")
}
