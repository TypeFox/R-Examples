# computes the RegMCD estimator
rrest <- function(data, lambda=0.5, hp=0.75, thresh=0.0001, maxit=10, penalty="L2")
{

# this function was taken from croux/haesbroeck and computes the L2 solution
regL2 <- function(scattermatrix,lambda) 
{ 
	prov <- eigen(scattermatrix) 
	thetai <- (-prov$values+sqrt(prov$values^2+8*lambda))/(4*lambda) 
	wi <- prov$vectors %*% diag(thetai) %*% t(prov$vectors) 
	return(list(wi=wi)) 
	#return(list(wi=wi,w=solve(wi))) 
}

# computes lp norm of matrix m
lpnorm <- function(m, p="L2")
{
	if(p == "L2")
	{
		ret <- sum(m^2)
	}
	else
	{
		ret <- sum(abs(m))
	}
	ret
}

	## check requirements
	if(!is.data.frame(data) && !is.matrix(data))
		stop("parameter data has to be a matrix or data frame!")

	if(penalty == "L2")
		p_func <- regL2
	else
		p_func <- glasso

	data <- as.data.frame(data)
		
	n <- dim(data)[1]
	h <- floor(n*hp)
	p <- dim(data)[2]
	n_iter <- 0
	na_count <- c()
	
	if(hp < 1) ## if h < n we have to do the whole procedure
	{
		na_count <- -1
		v_loglik <- NA
		loglik_val <- NA
		while(is.na(v_loglik))
		{
			na_count <- na_count+1
			pcoutw <- pcout(data)$wfinal # use pcout for outlier detection, and focus on the inner 2 points
			rsam <- order(pcoutw,decreasing=TRUE)[1:trunc(h/2)]

			rsub <- data[rsam,]
			v_cov <- var(rsub)
			v_mean <- colMeans(rsub)
			v_gl <- p_func(v_cov, lambda)

			ld <- log(det(v_gl$wi))
			md <- mean(mahalanobis(rsub, center=v_mean, cov=v_gl$wi,inverted=TRUE))
			v_loglik <- ld - md - lambda*lpnorm(v_gl$wi, penalty)		
			loglik_val <- (h/2)*(ld-md)
		}

		while(n_iter < maxit)
		{
			n_iter <- n_iter + 1
			v_mahal <- mahalanobis(x=data, center=v_mean, cov=v_gl$wi, inverted=TRUE)
			v_sorted <- sort.int(v_mahal, index.return=TRUE)
			rsam <- v_sorted$ix[1:h]
			rsub <- data[rsam,]
			v_cov <- var(rsub)
			v_mean <- colMeans(rsub)
			v_gl <- p_func(v_cov, lambda)
			ld <- log(det(v_gl$wi))
			md <- mean(mahalanobis(rsub, center=v_mean, cov=v_gl$wi,inverted=TRUE))
			v_loglik_new <- ld - md - lambda*lpnorm(v_gl$wi, penalty)		
			loglik_val <- (h/2)*(ld-md)
		
			## this break statement assures that the algorithm does not terminate with error msg
			if(is.na( abs((v_loglik_new - v_loglik)/v_loglik) ))
			{
				#print("NA produced, while loop left!")
				break;
			}	
			if(as.logical(abs((v_loglik_new - v_loglik)/v_loglik) < thresh))
			{
				#print(paste("Break condition reached after", n_iter, "iterations.", sep=" "))
				break;
			}	
			v_loglik <- v_loglik_new
		}
	}
	else ## hp = 1, so glasso is just computed once
	{	
		rsam <- 1:nrow(data)
		rsub <- data[rsam,]
		v_cov <- var(rsub)
		v_mean <- colMeans(rsub)
		v_gl <- p_func(v_cov, lambda)
		ld <- log(det(v_gl$wi))
		md <- mean(mahalanobis(rsub, center=v_mean, cov=v_gl$wi,inverted=TRUE))
		v_loglik_new <- ld - md - lambda*sum((v_gl$wi)^2)
		loglik_val <- (h/2)*(ld-md)
		n_iter <- 1
		nacount <- 0
	}
	l <- list(mean=v_mean, covi_nocons=v_gl$wi, subset=rsam, objective=v_loglik_new, loglik=loglik_val, niter=n_iter)
	l
}
