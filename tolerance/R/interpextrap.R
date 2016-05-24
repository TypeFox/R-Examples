#Internal Functions

two.sided <- function (x, alpha, P){
	n <- length(x)
	x <- sort(x)
	gamma <- 1-alpha
	out <- nptol.int(1:n, alpha = alpha, P = P, side = 2,method = "HM")[,3:4]
	r <- as.numeric(out[,1])
	s <- as.numeric(out[,2])
	if(nrow(out) == 2){
		X1.L <- x[c(r[1],r[1]+1)]
		X2.L <- x[c(r[2],r[2]+1)]
		X1.U <- x[c(s[1],s[1]-1)]
		X2.U <- x[c(s[2],s[2]-1)]
		g <- c(pbinom(s[1]-r[1]-1,n,P),pbinom(s[1]-(r[1]+1)-1,n,P))
		out1.L <- as.numeric(predict(lm(X1.L~g),newdata=list(g=gamma)))
		out2.L <- as.numeric(predict(lm(X2.L~g),newdata=list(g=gamma)))
		out1.U <- as.numeric(predict(lm(X1.U~g),newdata=list(g=gamma)))
		out2.U <- as.numeric(predict(lm(X2.U~g),newdata=list(g=gamma)))
		temp <- cbind(c(out1.L,out2.L,x[r[1]],x[r[2]]),c(x[s[1]],x[s[2]],out1.U,out2.U))
		temp <- cbind(temp,apply(temp,1,diff))
		if(pbinom(s[1]-r[1]-1,n,P)>=gamma){
			ind <- which.max(temp[,3])
			temp <- temp[ind,1:2]
				if(ind==1|ind==3) ord <- 1 else ord <- 2
				} else{
					ind <- which.max(temp[,3])
					temp <- temp[ind,1:2]
					if(ind==1|ind==3) ord <- 1 else ord <- 2
					}
			temp <- matrix(temp,nrow=1)
		} else{
		X.L <- x[c(r,r+1)]
		X.U <- x[c(s,s-1)]
		g <- c(pbinom(s-r-1,n,P),pbinom(s-(r+1)-1,n,P))
		out.L <- as.numeric(predict(lm(X.L~g),newdata=list(g=gamma)))
		out.U <- as.numeric(predict(lm(X.U~g),newdata=list(g=gamma)))
		temp <- cbind(c(out.L,x[r]),c(x[s],out.U))
		temp <- cbind(temp,apply(temp,1,diff))
		if(pbinom(s[1]-r[1]-1,n,P)>=gamma){
			temp <- temp[which.min(temp[,3]),1:2]
			} else{
 				temp <- cbind(out.L,out.U) #Extrapolates both sides.
				}
		temp <- matrix(temp,nrow=1)
		}
	temp <- cbind(alpha,P,temp)
	colnames(temp) <- c("alpha","P","2-sided.lower","2-sided.upper")
	rownames(temp) <- "OS-Based"
	temp
}

interp <- function (x, alpha, P){
	n <- length(x)
	x <- sort(x)
	gamma <- 1-alpha
	out <- as.numeric(nptol.int(1:n,alpha=alpha,P=P,side=1)[3:4])
	s <- as.numeric(out[1])
	r <- as.numeric(out[2])
	###Beran-Hall
	pi.l <- (gamma-pbinom(n-s-1,n,P))/dbinom(n-s,n,P)
	pi.u <- (gamma-pbinom(r-2,n,P))/dbinom(r-1,n,P)
	Q.l <- pi.l*x[s+1]+(1-pi.l)*x[s]
	Q.u <- pi.u*x[r]+(1-pi.u)*x[r-1]
	###Hutson
	f.l <- function(u1,u,n,alpha) pbeta(u,(n+1)*u1,(n+1)*(1-u1))-1+alpha
	f.u <- function(u2,u,n,alpha) pbeta(u,(n+1)*u2,(n+1)*(1-u2))-alpha
	u1 <- uniroot(f.l,interval=c(0.00001,0.99999),u=1-P,n=n,alpha=alpha)$root
	u2 <- uniroot(f.u,interval=c(0.00001,0.99999),u=P,n=n,alpha=alpha)$root
	eps <- function(u,n) (n+1)*u-floor((n+1)*u)
	Qh.l <- (1-eps(u1,n))*x[s]+eps(u1,n)*x[s+1]
	Qh.u <- (1-eps(u2,n))*x[r-1]+eps(u2,n)*x[r]
	temp <- data.frame(rbind(c(alpha,P,Q.l,Q.u),c(alpha,P,Qh.l,Qh.u)))
	colnames(temp) <- c("alpha","P","1-sided.lower","1-sided.upper")
	rownames(temp) <- c("OS-Based","FOS-Based")
	temp
}

extrap <- function (x, alpha, P){
	n <- length(x)
	x <- sort(x)
	gamma <- 1-alpha
	out.exp <- as.numeric(nptol.int(x,alpha=alpha,P=P,side=1)[3:4])
	###Beran-Hall
	pi.b <- -(gamma-pbinom(n-1,n,P))/dbinom(n-1,n,P)
	Qexp.l <- pi.b*x[2]+(1-pi.b)*x[1]
	Qexp.u <- pi.b*x[n-1]+(1-pi.b)*x[n]
	###Hutson
	f.lb <- function(u1,u,n,alpha) pbeta(u,(n+1)*u1,(n+1)*(1-u1))-1+alpha
	f.ub <- function(u2,u,n,alpha) pbeta(u,(n+1)*u2,(n+1)*(1-u2))-alpha
	u1.b <- uniroot(f.lb,interval=c(0.00001,0.99999),u=1-P,n=n,alpha=alpha)$root
	u2.b <- uniroot(f.ub,interval=c(0.00001,0.99999),u=P,n=n,alpha=alpha)$root
	eps <- function(u,n) -((n+1)*u-floor((n+1)*u))
	Qhexp.l <- (1-eps(u1.b,n))*x[1]+eps(u1.b,n)*x[2]
	Qhexp.u <- (1-eps(u2.b,n))*x[n]+eps(u2.b,n)*x[n-1]
	temp <- data.frame(rbind(c(alpha,P,Qexp.l,Qexp.u),c(alpha,P,Qhexp.l,Qhexp.u)))
	colnames(temp) <- c("alpha","P","1-sided.lower","1-sided.upper")
	rownames(temp) <- c("OS-Based","FOS-Based")
	temp
}











