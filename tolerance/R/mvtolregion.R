mvtol.region <- function (x, alpha = 0.05, P = 0.99, B = 1000, M = 1000, 
	method = c("KM", "AM", "GM", "HM", "MHM", "V11", "HM.V11", "MC")) 
{
	method <- match.arg(method)
	P <- sort(P)
	alpha <- sort(alpha, decreasing = TRUE)
    n <- nrow(x)
    p <- ncol(x)
	if(method == "KM"){
    q.squared <- matrix(rchisq(p * B, df = 1), ncol = p)/n
    L <- t(sapply(1:B, function(i) eigen(rwishart(n - 1, p))$values))
    c1 <- apply((1 + q.squared)/L, 1, sum)
    c2 <- apply((1 + 2 * q.squared)/L^2, 1, sum)
    c3 <- apply((1 + 3 * q.squared)/L^3, 1, sum)
    a <- (c2^3)/(c3^2)
    T <- matrix(sapply(1:length(P), function(i) (n - 1) * (sqrt(c2/a) * 
        (qchisq(P[i], a) - a) + c1)), ncol = length(P))
    tol <- matrix(sapply(1:length(alpha), function(i) apply(T, 
        2, quantile, 1 - alpha[i])), ncol = length(alpha))
	} else if(method == "AM"){
		tol <- sapply(1:length(alpha), function(i) (p*(n-1)*qchisq(P,p,p/n))/(qchisq(alpha[i],(n-1)*p)))
		} else if(method == "GM"){
			g1 <- (p/2)*(1-((p-1)*(p-2))/(2*n))^(1/p)
			tol <- sapply(1:length(alpha), function(i) g1*(n-1)*qchisq(P,p,p/n)/qgamma(alpha[i],p*(n-p)/2))
			} else if(method == "HM"){
				tol <- sapply(1:length(alpha), function(i) (p*(n-1)*qchisq(P,p,p/n))/(qchisq(alpha[i],(n-1)*p-p*(p+1)+2)))
				} else if(method == "MHM"){
					b <- (p*(n-p-1)*(n-p-4)+4*(n-2))/(n-2)
					a <- (p*(b-2))/(n-p-2)
					tol <- sapply(1:length(alpha), function(i) (a*(n-1)*qchisq(P, p, p/n))/(p*qchisq(alpha[i],b)))
					} else if(method == "V11"){
						tol <- sapply(1:length(alpha), function(i) (n-1)*qchisq(P,p,p/n)/qchisq(alpha[i],n-p))
						} else if(method == "HM.V11"){
							e <- (4*p*(n-p-1)*(n-p)-12*(p-1)*(n-p-2))/(3*(n-2)+p*(n-p-1))
							d <- (e-2)/(n-p-2)
							tol <- sapply(1:length(alpha), function(i) d*(n-1)*qchisq(P,p,p/n)/qchisq(alpha[i],e))
							} else if(method == "MC"){
								U <- matrix(rnorm(B*p,0,1/n),nrow=B,ncol=p)
								V <- lapply(1:B, function(i) solve(rwishart(n - 1, p)))
								Y <- lapply(1:B, function(i) matrix(rnorm(M*p,0,1),nrow=p,ncol=M)) 
								RES <- lapply(1:B, function(i) Y[[i]]-U[i,])
								T <- matrix(sapply(1:length(P), function(j) sapply(1:B, function(i) quantile((n-1)*diag(t(RES[[i]])%*%V[[i]]%*%(RES[[i]])),P[j]))),ncol=length(P))
								tol <- matrix(sapply(1:length(alpha), function(i) apply(T, 2, quantile, 1 - alpha[i])), ncol = length(alpha))
								}
    colnames(tol) <- alpha
    rownames(tol) <- P
    tol
}
