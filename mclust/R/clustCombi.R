clustCombi <- function(data, MclustOutput = NULL, ...)
{
	if (is.null(MclustOutput)) 
		{
			MclustOutput <- Mclust(data, ...)
		}
	
	combiRes <- combi(data, MclustOutput)
	
	return(combiRes)
}

combMat <- function(K,l1,l2)
{
	l=c(min(l1,l2), max(l1,l2))
	if(any(length(l1) == 0, length(l2) == 0)){
		l1 = numeric(0)
	  l2 = l[2]}
	 else {
	 	l1 = l[1]
	 	l2 = l[2]} 
	
	M <- rbind(cbind(diag(l2-1),matrix(rep(0,(K-l2+1)*(l2-1)), nrow=l2-1, ncol=K-l2+1)), cbind(matrix(rep(0,l2*(K-l2)), nrow=K-l2, ncol=l2),diag(K-l2)))
	M[l1,l2] <- 1
	return(M)
}

## Define xlog to handle x*log(x) as x=0

xlog <- function(x) 
	{
		xlog1d <- function (xi) if (xi == 0) 0 else (xi*log(xi))
		
		if (is.null(dim(x)))
			{
				return(sapply(x,xlog1d))
			}
		else
			{
				return(matrix(sapply(x,xlog1d),dim(x)))
			}
	}

combi <- function(data, MclustOutput, n = nrow(data), d = ncol(data))
{
	combiM <- list()
	combiM[[MclustOutput$G]] <- diag(MclustOutput$G)
	tau <- list()
	tau[[MclustOutput$G]] = MclustOutput$z
	classif <- list()
	classif[[MclustOutput$G]] = map(tau[[MclustOutput$G]])

	for (K in MclustOutput$G:2)
		{
			dEnt <- matrix(0,nrow=K-1, ncol=K)
			preCombiTau <- tau[[K]]
			for (l1 in 1:(K-1))
				{
					for (l2 in (l1+1):K)
						{
							postCombiTau <- t(combMat(K,l1,l2) %*% t(preCombiTau))
							dEnt[l1,l2] <- sum(xlog(postCombiTau[,l1])) - sum(xlog(preCombiTau[,l1])+xlog(preCombiTau[,l2]))
						}	
				}	
			l1=which(dEnt==max(dEnt),arr.ind=TRUE)[1]
			l2=which(dEnt==max(dEnt),arr.ind=TRUE)[2]
			
			combiM[[K-1]] <- combMat(K,l1,l2)
			tau[[K-1]] = t(combiM[[K-1]] %*% t(tau[[K]]))
			classif[[K-1]] = map(tau[[K-1]])
		}
	
	output <- list(classification = classif, combiM = combiM, combiz = tau, MclustOutput = MclustOutput)
	class(output) <- "clustCombi"
	return(output)
}

plot.clustCombi <- function(x, data = NULL, what = c("classification", "entropy"), reg = 2, ...)
{
  output <- x # Argh. Really want to use 'output'

	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar))
	
	if(any(what == "classification") && is.null(data)) 
	  warning("Data not supplied: classification plots can't be displayed")
		
	if (any(what == "classification") && !is.null(data)) 
		{
			## Sort z columns so that one of the two combined column is the last one at each step (prevents the colors and symbols to be mixed as K -> K-1)
			
			curr = 1:output$MclustOutput$G
			i <- numeric()
			j <- numeric()
			
			for (K in (output$MclustOutput$G):2)
				{
						l1 = which(!output$combiM[[K-1]] %*% rep(1,K) == 1)
						l2 = (output$combiM[[K-1]] %*% curr)[l1] - curr[l1]
						
						i <- c(curr[l1],i)
						j <- c(l2,j)

						curr <- output$combiM[[K-1]] %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(K-1-l1)))
				}
				
			j <- c(1,j)
			i <- c(0,i)
			
			permutz <- output$MclustOutput$z[,j] 
			
			permutMat <- function(j,K) 
				{
					M <- diag(K)
					M[j,j] <- 0
					M[K,K] <- 0
					M[j,K] <- 1
					M[K,j] <- 1	
					
					return(M)
				}	
			
			combiM = diag(output$MclustOutput$G)
			
			for (K in output$MclustOutput$G:1)
				{					
					combiPlot(data = data, z = permutz, combiM = combiM, ...)
					curr_title <- if (K == output$MclustOutput$G) paste( "BIC solution (",as.character(K)," clusters)") else paste( "Combined solution with ",as.character(K)," clusters")
					if (ncol(as.matrix(data)) > 2) par(oma=c(0,0,2,0), mar = c(5,4,0,2)+0.1)
					if (ncol(as.matrix(data)) > 2) title(curr_title, outer = TRUE) else title(curr_title)
					combiM = combMat(K,which(j==i[K]),K) %*% combiM 
					par(ask=TRUE)	
				}	

		}
		
	if(any(what == "entropy")) 
	  entPlot(z = output$MclustOutput$z, combiM = output$combiM, reg = reg, ...)

}

combiPlot <- function(data, z, combiM, ...)
{
  p <- ncol(as.matrix(data))
  if (p > 2) {
  		clPairs(data[,1:min(5,p)], classification = map(t(combiM %*% t(z))), ...) # To be tried
  	}
	else if (p == 2) {
		mclust2Dplot(data = data, parameters = NULL, classification = map(t(combiM %*% t(z))), what = "classification", ...)
	}
	else {
		mclust1Dplot(data = as.matrix(data), parameters = NULL, classification = map(t(combiM %*% t(z))), what = "classification", ...) # To be tried
	}
}

entPlot <- function(z, combiM, abc = c("standard", "normalized"), reg = 2, ...)
{
	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar))

	ent <- numeric()
	Kmax <- ncol(z)
	z0 <- z
	for (K in Kmax:1) 
		{
			z0 <- t(combiM[[K]] %*% t(z0))
			ent[K] <- -sum(xlog(z0))
		}

	if (any(abc == "normalized"))
		{
			mergedn <- numeric()
			z0 <- z
			for (K in (Kmax-1):1)
				{
					z0 <- t(combiM[[K+1]] %*% t(z0))
					mergedn[K] = sum(sapply(map(z0),function(x) any(which(as.logical(combiM[[K]][rowSums(combiM[[K]])==2,]))==x)))				}
		}
		
	if (Kmax == 2) reg = NULL
		
	if (any(abc == "standard"))
		{
			par(mfrow=c(1,2), oma=c(0,0,3,0), mar = c(4,4,0,1)+0.1)
			plot(1:Kmax, ent, xlab = "Number of clusters", ylab = "Entropy", ...)
			if (any(reg == 2)) 
				{	
					pcwsreg <- pcws2_reg(1:Kmax,ent)
					lines(1:pcwsreg$c, pcwsreg$a1*(1:pcwsreg$c) + pcwsreg$b1, lty = 2, col = "red")
					lines(pcwsreg$c:Kmax, pcwsreg$a2*(pcwsreg$c:Kmax) + pcwsreg$b2, lty = 2, col = "red")
				}
			if (any(reg == 3)) 
				{	
					pcwsreg <- pcws3_reg(1:Kmax,ent)
					lines(1:pcwsreg$c1, pcwsreg$a1*(1:pcwsreg$c1) + pcwsreg$b1, lty = 2, col = "blue")
					lines(pcwsreg$c1:pcwsreg$c2, pcwsreg$a2*(pcwsreg$c1:pcwsreg$c2) + pcwsreg$b2, lty = 2, col = "blue")
					lines(pcwsreg$c2:Kmax, pcwsreg$a3*(pcwsreg$c2:Kmax) + pcwsreg$b3, lty = 2, col = "blue")
				}
			plot(1:(Kmax-1), ent[2:Kmax]-ent[1:(Kmax-1)], xlab = "Number of clusters", ylab = "Difference in entropy", ...)
			title("Entropy plot", outer=TRUE)
			par(ask=TRUE)
		}
	if (any(abc == "normalized"))
		{
			par(mfrow=c(1,2), oma=c(0,0,3,0), mar = c(4,4,0,1)+0.1)
			plot(cumsum(c(0,mergedn)), ent, xlab = "Cumul. count of merged obs.", ylab = "Entropy", ...)
			if (any(reg == 2)) 
				{	
					X <- cumsum(c(0,mergedn))
					pcwsreg <- pcws2_reg(X,ent)
					lines(X[1:pcwsreg$c], pcwsreg$a1*(X[1:pcwsreg$c]) + pcwsreg$b1, lty = 2, col = "red")
					lines(X[pcwsreg$c:Kmax], pcwsreg$a2*(X[pcwsreg$c:Kmax]) + pcwsreg$b2, lty = 2, col = "red")
				}
			if (any(reg == 3)) 
				{	
					X <- cumsum(c(0,mergedn))
					pcwsreg <- pcws3_reg(X,ent)
					lines(X[1:pcwsreg$c1], pcwsreg$a1*(X[1:pcwsreg$c1]) + pcwsreg$b1, lty = 2, col = "blue")
					lines(X[pcwsreg$c1:pcwsreg$c2], pcwsreg$a2*(X[pcwsreg$c1:pcwsreg$c2]) + pcwsreg$b2, lty = 2, col = "blue")
					lines(X[pcwsreg$c2:Kmax], pcwsreg$a3*(X[pcwsreg$c2:Kmax]) + pcwsreg$b3, lty = 2, col = "blue")
				}
			plot(1:(Kmax-1), (ent[2:Kmax]-ent[1:(Kmax-1)])/mergedn, xlab = "Number of clusters", ylab = "Normalized difference in entropy", ...)
			title("Normalized entropy plot", outer=TRUE)
		}

}

# pcws2_reg computes the piecewise linear regression -- with two pieces -- to (x,y), for any possible change point and chooses the one leading to the smallest least-square error.

pcws2_reg <- function(x, y)
{
	
	C <- length(x)
	ssBest = Inf
	for (c in 2:(C-1))
		{
			x1 <- x[1:c]
			y1 <- y[1:c]
			x2 <- x[c:C]
			y2 <- y[c:C]
			
			a1 <- sum((x1-mean(x1))*(y1-mean(y1)))/sum((x1-mean(x1))^2)
			b1 <- -a1 * mean(x1) + mean(y1)

			a2 <- sum((x2-mean(x2))*(y2-mean(y2)))/sum((x2-mean(x2))^2)
			b2 <- -a2 * mean(x2) + mean(y2)
			
			ss <- sum((a1*x1+b1-y1)^2) + sum((a2*x2+b2-y2)^2)

			if (ss < ssBest) 
				{
					ssBest <- ss
					cBest <- c
					a1Best <- a1
					a2Best <- a2
					b1Best <- b1
					b2Best <- b2
				}
		}
	
	return(list(c=cBest, a1=a1Best, b1=b1Best, a2=a2Best, b2=b2Best, residuals = c(a1*x1+b1-y1,a2*x2+b2-y2)))
}

# pcws3_reg computes the piecewise linear regression -- with three pieces -- to (x,y), for any possible change points and chooses the ones leading to the smallest least-square error.

pcws3_reg <- function(x, y)
{
	
	C <- length(x)
	ssBest = Inf
	for (c1 in 2:(C-2))
		{
			for (c2 in (c1+1):(C-1))
				{
			
					x1 <- x[1:c1]
					y1 <- y[1:c1]
					x2 <- x[c1:c2]
					y2 <- y[c1:c2] 
					x3 <- x[c2:C]
					y3 <- y[c2:C] 
			
					a1 <- sum((x1-mean(x1))*(y1-mean(y1)))/sum((x1-mean(x1))^2)
					b1 <- -a1 * mean(x1) + mean(y1)

					a2 <- sum((x2-mean(x2))*(y2-mean(y2)))/sum((x2-mean(x2))^2)
					b2 <- -a2 * mean(x2) + mean(y2)
			
					a3 <- sum((x3-mean(x3))*(y3-mean(y3)))/sum((x3-mean(x3))^2)
					b3 <- -a3 * mean(x3) + mean(y3)

			ss <- sum((a1*x1+b1-y1)^2) + sum((a2*x2+b2-y2)^2) + sum((a3*x3+b3-y3)^2)

			if (ss < ssBest) 
				{
					ssBest <- ss
					c1Best <- c1
					c2Best <- c2
					a1Best <- a1
					b1Best <- b1
					a2Best <- a2
					b2Best <- b2
					a3Best <- a3
					b3Best <- b3
				}
		}
	}
	return(list(c1=c1Best, c2=c2Best, a1=a1Best, b1=b1Best, a2=a2Best, b2=b2Best, a3=a3Best, b3=b3Best, residuals = c(a1*x1+b1-y1,a2*x2+b2-y2,a3*x3+b3-y3)))
}

print.clustCombi <- function(x, ...)
{
  output <- x # Argh. Really want to use 'output'
	cat("\n EM/BIC Solution\n")
	cat(" --------------- \n\n")
	cat("Number of components: ", as.character(output$MclustOutput$G), "\n", sep = "") 

	cat("Model name: ", output$MclustOutput$parameters$var$modelName, "\n\n", sep="")
	for (K in 1:output$MclustOutput$G)
		{
			cat("Component num.", as.character(K),": ", "\n", sep="")
			cat("		proportion: ", sprintf(fmt = "%4.2f ", output$MclustOutput$parameters$pro[K]), "\n", sep="")
			if (output$Mclust$d == 1) cat("		mean: ", sprintf(fmt = "%4.2f  ", output$MclustOutput$parameters$mean[K]), "\n", sep="") else cat("		mean: ", sprintf(fmt = "%4.2f  ", output$MclustOutput$parameters$mean[,K]), "\n", sep="")
		}
	
	cat("\n Combining steps \n")
	cat(" --------------- \n\n")

	cl = paste(rep(" ", max(output$MclustOutput$G-4,0)), "Classes labels after this step", rep(" ", max(output$MclustOutput$G-4,0)), sep="")
	
	if (output$MclustOutput$G>4) for (K in 5:output$MclustOutput$G) cl = paste(" ", cl, " ", sep="")
	
	cat("  Step | Classes combined at this step | Classes labels after this step", "\n", sep="")
	cat("-------|-------------------------------|-------------------------------", "\n", sep="")
	curr = 1:output$MclustOutput$G

	cat("   0   |              ---              |", sprintf(fmt = "%2d  ", curr), "\n", sep="")
			
	for (K in 1:(output$MclustOutput$G-1))
		{
			Kp = output$MclustOutput$G - K + 1
			l1 = which(!output$combiM[[Kp-1]] %*% rep(1,Kp) == 1)
			l2 = (output$combiM[[Kp-1]] %*% curr)[l1] - curr[l1]

			nc1 = floor((7-nchar(as.character(K)))/2)
			nc2 = (7-nchar(as.character(K))) - nc1
			nc3 = floor((33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))))/2) 
			nc4 = 33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))) - nc3
			 
			curr <- output$combiM[[Kp-1]] %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(Kp-1-l1)))

			cat(rep(" ", nc1), as.character(K), rep(" ", nc2), "|", rep(" ", nc3), as.character(l1), "  &  ", as.character(l2), rep(" ", nc4), "|", sprintf(fmt = "%2d  ", curr), "\n", sep="")
		
		}

	cat("\n Classification for K classes: output$classification[[K]]\n")
	cat(" Combining matrix (K classes -> (K-1) classes): output$combiM[[K]]\n\n")
}
