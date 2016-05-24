SPIn <-
function (x, conf = 0.95, bw = 0, lb = -Inf, ub = Inf, l=NA, u=NA) 
{
    if (lb != -Inf) 
        x <- c(x, lb)
    if (ub != Inf) 
        x <- c(x, ub)
    dens <- density(x)
    if (is.na(l)){
    n.sims <- length(x)
    conf <- 1 - conf
    nn <- round(n.sims * conf)
    x <- sort(x)
    xx <- x[(n.sims - nn):n.sims] - x[1:(nn + 1)]
    m <- min(xx)
    k <- which(xx == m)[1]
    l <- x[k]
    ui <- n.sims - nn + k - 1
    u <- x[ui]
    } else{
		if (sum(x==l)==0){
			x <- c(x,l)
		}
		if (sum(x==u)==0){
			x <- c(x,u)
		}
		x <- sort(x)
		n.sims <- length(x)
    }
    library(quadprog)
    if (bw <= 0) 
        bw <- round((sqrt(n.sims)-1)/2)
    
#	x <- sort(x)
#	n.sims <- length(x)
	k <- which(x==l)[1]
	ui <- which(x==u)[1]
	l.l <- max(1,k-bw)
	l.u <- k+(k-l.l)
	u.u <- min(n.sims,ui+bw)
	u.l <- ui-(u.u-ui)
	
#	D <- matrix(nrow=n.sims,ncol=n.sims)
	n.l <- l.u-l.l+1
	n.u <- u.u-u.l+1
	D.l <- matrix(nrow=n.l,ncol=n.l)
	D.u <- matrix(nrow=n.u,ncol=n.u)
	
	p <- (l.l:l.u)/(n.sims+1)
	q <- 1 - p
	Q <- quantile(x,p)
	d.q <- rep(0,n.l)
	for (r in 1:n.l) d.q[r] <- dens$y[which.min(abs(dens$x-Q[r]))]
	Q. <- 1/d.q
	diag(D.l) <- 2 * (Q^2 + p*q*Q.^2/(n.sims+2))
	d.l <- 2*Q*l
	if(n.l > 1){
	for (j in 1:(n.l-1))
		for (m in (j+1):n.l){
			D.l[j,m] <- Q.[j]*Q.[m]*p[j]*q[m]*2/(n.sims+2)+Q[j]*Q[m]*2
			D.l[m,j] <- D.l[j,m]
		}
	}
	
	p <- (u.l:u.u)/(n.sims+1)
	q <- 1 - p
	Q <- quantile(x,p)
	d.q <- rep(0,n.u)
	for (r in 1:n.u) d.q[r] <- dens$y[which.min(abs(dens$x-Q[r]))]
	Q. <- 1/d.q
	diag(D.u) <- 2 * (Q^2 + p*q*Q.^2/(n.sims+2))
	d.u <- 2*Q*u
	if(n.u > 1){
	for (j in 1:(n.u-1))
		for (m in (j+1):n.u){
			D.u[j,m] <- Q.[j]*Q.[m]*p[j]*q[m]*2/(n.sims+2)+Q[j]*Q[m]*2
			D.u[m,j] <- D.u[j,m]
		}
	}
	
	if (k == 1){
		x1 <- l
		w.l <- 1
	} else{
	A.l <- matrix(0, nrow= l.u-l.l+3, ncol= l.u-l.l+1)
	A.l[1,] <- 1
	if (bw > 1){
	if(k > 2){
	for (j in 1:(k-l.l-1)){
		if (x[l.l+j+1] == x[l.l+j]){
			A.l[1+j,j+1] <- 1
			A.l[1+j,j+2] <- -1
		} else{
		aa <- (x[l.l+j]-x[l.l+j-1])/(x[l.l+j+1]-x[l.l+j])
		A.l[1+j,j] <- 1
		A.l[1+j,j+1] <- -(aa+1)
		A.l[1+j,j+2] <- aa
		}
	}

	for (j in 0:(l.u-k-2)){
		if (x[k+j+1] == x[k+j+2]){
			A.l[k-l.l+1+j,k-l.l+2+j] <- 1
			A.l[k-l.l+1+j,k-l.l+3+j] <- -1
		} else{
		aa <- (x[k+j]-x[k+j+1])/(x[k+j+1]-x[k+j+2])
		A.l[k-l.l+1+j,k-l.l+1+j] <- -1
		A.l[k-l.l+1+j,k-l.l+2+j] <- aa+1
		A.l[k-l.l+1+j,k-l.l+3+j] <- -aa
		}
	}
	}
	}
	if (x[k+1] == x[k]){
		aa <- (x[k]-x[k-1])/(x[k+1]-x[k]+.000001)
	} else{
		aa <- (x[k]-x[k-1])/(x[k+1]-x[k])
	}
	A.l[l.u-l.l,k-l.l+1] <- aa-1
	A.l[l.u-l.l,k-l.l] <- 1
	A.l[l.u-l.l,k-l.l+2] <- -aa
	
	A.l[l.u-l.l+1,l.u-l.l] <- 1
	A.l[l.u-l.l+1,l.u-l.l+1] <- -1
	A.l[l.u-l.l+2,1] <- 1
	A.l[l.u-l.l+3,l.u-l.l+1] <- 1
	A.l <- t(A.l)
	w.l <- solve.QP(D.l,d.l,A.l,c(1,rep(0,l.u-l.l+2)),l.u-l.l)
	w.l <- w.l$solution
	x1 <- w.l%*%x[l.l:l.u]
	}
	
	if (ui == n.sims){
		x2 <- u
		w.u <- 1
	} else{
	A.u <- matrix(0, nrow= u.u-u.l+3, ncol= u.u-u.l+1)
	A.u[1,] <- 1
	if (bw > 1){
	if (ui-u.l>1){
	for (j in 1:(ui-u.l-1)){
		if (x[u.l+j+1] == x[u.l+j]){
			A.u[1+j,j+1] <- 1
			A.u[1+j,j+2] <- -1
		} else{
		aa <- (x[u.l+j]-x[u.l+j-1])/(x[u.l+j+1]-x[u.l+j])
		A.u[1+j,j] <- 1
		A.u[1+j,j+1] <- -(aa+1)
		A.u[1+j,j+2] <- aa
		}
	}

	i <- 0
	for (j in (ui-u.l):(u.u-u.l-2)){
		if (x[ui+i+1] == x[ui+i+2]){
			A.u[1+j,j+2] <- 1
			A.u[1+j,j+3] <- -1
		} else{
		aa <- (x[ui+i]-x[ui+i+1])/(x[ui+i+1]-x[ui+i+2])
		A.u[1+j,j+1] <- -1
		A.u[1+j,j+2] <- aa+1
		A.u[1+j,j+3] <- -aa
		}
		i <- i + 1
	}
	}
	}
	if (x[ui+1] == x[ui]){
#		aa <- (x[ui]-x[ui-1])/(x[ui+1]-x[ui]+.000001)
		aa <- (x[ui]-x[ui-1])/(x[ui+2]-x[ui])
		A.u[u.u-u.l,ui-u.l] <- 1
		A.u[u.u-u.l,ui-u.l+1] <- aa-1
		A.u[u.u-u.l,ui-u.l+3] <- -aa
	} else{
		aa <- (x[ui]-x[ui-1])/(x[ui+1]-x[ui])
		A.u[u.u-u.l,ui-u.l] <- 1
		A.u[u.u-u.l,ui-u.l+1] <- aa-1
		A.u[u.u-u.l,ui-u.l+2] <- -aa
	}

	A.u[u.u-u.l+1,u.u-u.l] <- 1
	A.u[u.u-u.l+1,u.u-u.l+1] <- -1
	A.u[u.u-u.l+2,1] <- 1
	A.u[u.u-u.l+3,u.u-u.l+1] <- 1
	A.u <- t(A.u)
	w.u <- solve.QP(D.u,d.u,A.u,c(1,rep(0,u.u-u.l+2)),u.u-u.l)
	w.u <- w.u$solution
	x2 <- w.u%*%x[u.l:u.u]
	}
    hpd <- list(spin = c(x1, x2), conf = 1 - conf, x = x, w.l=w.l, w.u=w.u, l.l=l.l, l.u=l.u, u.l=u.l, u.u=u.u)
    class(hpd) <- "SPIn"
    return(hpd)
}
