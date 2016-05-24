## performs LDA with cov matrix computed by repeated rrest; computes bic value
rrlda <- function(x, grouping, prior=NULL, lambda=0.5, hp=0.75, nssamples=30, maxit=50, penalty="L2")
{

# computes repeated rrest and returns object with largest objective value
rep_rrest <- function(rep=10, data, lambda=0.5, hp=0.75, thresh=0.0001, maxit=10, penalty="L2")
{
	max_obj <- -Inf
	obj <- NULL
	
	
	for(i in 1:rep)
	{
		#print(paste("Lambda=", lambda, "; ", "Repetition: ", i, "\n", sep=""))
		mrr <- rrest(data=data, lambda=lambda, hp=hp, thresh=thresh, maxit=maxit, penalty=penalty)
		if(mrr$objective > max_obj)
		{
			max_obj <- mrr$objective
			obj <- mrr
		}
	}
	obj
}

# computes degrees of freedom for concentration matrix m, lambda and penalty p
degf <- function(m, lambda, p)
{
	if(p == "L2")
	{
		ret <- alt.df.L2glasso(solve(m), lambda)
	}
	else
	{
		ret <- sum(m[upper.tri(m,diag=T)]!=0)
	}
	ret
}

getMeans <- function(gr_means, est_mean, lev, grouping)
{
	m_final <- t(t(gr_means) + est_mean)
	
	ind_vector <- c()
	for(i in 1:length(lev))
	{
		ind_vector <- c(ind_vector, which(grouping == lev[i]))
	}
	m_final <- m_final[ind_vector, ]
	means <- unique(m_final)
	rownames(means) <- lev
	means
}

# centers the observations using L1 median
centerObservations <- function(x, k, grouping, lev)
{
	m <- t(x)
	for(i in 1:k)
	{
		if(!is.vector(x[grouping==lev[i],]))
		{
			m[,grouping==lev[i]] <- l1median_NLM(x[grouping==lev[i],])$par
		}
		else
		{
			m[,grouping==lev[i]] <- as.vector(x[grouping==lev[i],])
		}
	}
	m <- t(m)
	z <- x - m
	l <- list(z=z, m=m)
	l
}

# returns prior probs
getPriors <- function(pr, k, props)
{
	m_prior <- pr
	if(!is.null(m_prior)) 
	{
        if(any(m_prior < 0) || round(sum(m_prior), 5) != 1) 
		{
            stop("invalid prior")
		}
        if(length(m_prior) != k) 
		{
            stop("'prior' is of incorrect length")
		}
    }
	else
	{
		m_prior <- props
	}
	m_prior
}

# computes degrees of freedom for L2 penalty; computation intense!
df.L2glasso<-function(sigma,lambda)
{
	p=ncol(sigma)
	G=duplication.matrix(p)
	a1=kronecker(sigma,sigma)
	a2=a1+2*lambda*diag(rep(1,p*p))
	a1=t(G)%*%a1%*%G;
	a2=t(G)%*%a2%*%G;
	a=solve(a2)%*%a1
	df=p+sum(diag(a))
	return(df)
}

# alternative: computes degrees of freedom for L2 penalty
alt.df.L2glasso<-function(sigma,lambda)
{
	p=ncol(sigma)
	eigvalue=eigen(sigma)$val
	help.matrix=eigvalue%*%t(eigvalue)
	df=p+sum(lower.triangle(help.matrix/(help.matrix+2*lambda)))
	return(df)
}
	## check requirements
    if (is.null(dim(x))) 
        stop("'x' is not a matrix or data.frame")

	x <- as.matrix(x)
    if (any(!is.finite(x))) 
        stop("Infinite, NA or NaN values in 'x'")

	n <- nrow(x)
	p <- ncol(x)
	h <- floor(n*hp)
	
    if (n != length(grouping)) 
        stop("nrow(x) and length(grouping) are different")	

	if(!is.factor(grouping))
		grouping <- as.factor(grouping)

	lev = levels(grouping)
	k = nlevels(grouping)
	counts = table(grouping)
	proportions <- counts/n
	
	co <- centerObservations(x, k, grouping, lev)
	est <- rep_rrest(rep=nssamples, co$z, lambda=lambda, hp=hp, maxit=maxit, penalty=penalty)
	means <- getMeans(co$m, est$mean, lev, grouping)
	bic <- (-2)*est$loglik + log(h)*(k*p + degf(est$covi_nocons, lambda, penalty)) 
	
	cl <- match.call()	
    cl[[1L]] <- as.name("rrlda")
	 
	ret = structure(list(call=cl, prior=getPriors(prior, k, proportions), counts=counts, means=means, covi=est$covi_nocons, lev=lev, n=n, h=h, bic=bic, loglik=est$loglik, df=degf, subs=est$subset), class="rrlda")
	ret
}

# predicts group membership for new data x and object from rrlda
predict.rrlda <- function(object, x, ...)
{
	## check requirements
	if(class(object) != "rrlda")
		stop("object has to be of type rrlda")
		
	## get dimension	
	p = ncol(object$means)

	if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x))
		stop("'x' has to be a vector, matrix or data.frame")
	
	n <- 0
	if(is.vector(x))
	{
		if(length(x) != p)
			stop("'x' has to be of length ", p)

		y <- NULL
		for(i in 1:length(x))
			y <- cbind(y, x[i])
		x <- y
		n <- 1
	}
	else
	{
		if(ncol(x) != p)
			stop("'x' has to be of dimension ", p)
			
		x <- as.matrix(x)
		n <- nrow(x)
	}
		
	probs <- t(-2*object$means %*% object$covi %*% t(x) + diag(object$means %*% object$covi %*% t(object$means)) - 2*log(object$prior[[1]]))
	cl = object$lev[apply(probs, 1, which.min)]
	
	ret = list(class=cl, posterior=probs)
	ret
}
