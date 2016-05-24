frmhet.links <- function(link) 
{
	switch(link,
		logit = {
			linkfun <- function(mu) qlogis(mu)
			linkinv <- function(eta) plogis(eta)
			mu.eta <- function(eta) dlogis(eta)
			H1 <- function(dep) dep/(1-dep)
			G2 <- function(eta) exp(eta)
			g2 <- function(eta) exp(eta)
			H2 <- function(uu) log(uu)
			gd <- function(eta) exp(eta)*(1-exp(eta))/((1+exp(eta))^3)
			valideta <- function(eta) TRUE
		},

		probit = {
			linkfun <- function(mu) qnorm(mu)
			linkinv <- function(eta) {
				thresh <- -qnorm(.Machine$double.eps)
				eta <- pmin(pmax(eta, -thresh), thresh)
	 			pnorm(eta)
			}
			mu.eta <- function(eta) pmax(dnorm(eta),.Machine$double.eps)
			H1 <- NULL
			G2 <- NULL
			g2 <- NULL
			H2 <- NULL
			gd <- function(eta) -eta*dnorm(eta)
			valideta <- function(eta) TRUE
		},

		cauchit = {
		  	linkfun <- function(mu) qcauchy(mu)
			linkinv <- function(eta) {
				thresh <- -qcauchy(.Machine$double.eps)
				eta <- pmin(pmax(eta, -thresh), thresh)
				pcauchy(eta)
			}
			mu.eta <- function(eta) pmax(dcauchy(eta),.Machine$double.eps)
			H1 <- NULL
			G2 <- NULL
			g2 <- NULL
			H2 <- NULL
			gd <- function(eta) -2*eta/(pi*(eta^2+1)^2)
			valideta <- function(eta) TRUE
		},

		cloglog = {
			linkfun <- function(mu) log(-log(1 - mu))
			linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)),1-.Machine$double.eps),.Machine$double.eps)
			mu.eta <- function(eta) {
				eta <- pmin(eta, 700)
				pmax(exp(eta) * exp(-exp(eta)),.Machine$double.eps)
			}
			H1 <- function(dep) -log(1-dep)
			G2 <- function(eta) exp(eta)
			g2 <- function(eta) -exp(-eta)
			H2 <- function(uu) log(uu)
			gd <- function(eta) (exp(-exp(eta))*exp(eta))*(1-exp(eta))
			valideta <- function(eta) TRUE
		},

		loglog = {
			linkfun <- function(mu) -log(-log(mu))
			linkinv <- function(eta) exp(-exp(-eta))
			mu.eta <- function(eta) exp(-exp(-eta)-eta)
			H1 <- NULL
			G2 <- NULL
			g2 <- NULL
			H2 <- NULL
			gd <- function(eta) (exp(-exp(-eta))*exp(-eta))*(exp(-eta)-1)
			valideta <- function(eta) TRUE
		}
	)

	structure(list(linkfun=linkfun,linkinv=linkinv,mu.eta=mu.eta,H1=H1,G2=G2,g2=g2,H2=H2,gd=gd,valideta=valideta,name=link),class="link-glm")
}

frmhet.gi <- function(type,z,Hy,XB,link)
{
	if(type!="QMLxv")
	{
		G2 <- frmhet.links(link)$G2(XB)
		u <- Hy/G2-1
	}
	if(type=="QMLxv")
	{
		yhat <- frmhet.links(link)$linkinv(XB)
		g <- frmhet.links(link)$mu.eta(XB)
		u <- Hy-yhat
		u <- g*u/(yhat*(1-yhat))
	}

	gi <- t(z*u)

	return(list(gi=gi,u=u))
}

frmhet.Gn <- function(type,x,z,z.in,u,bet,p)
{
	N <- length(u)

	if(all(type!=c("GMMxv","LINxv","QMLxv"))) Gn <- -(1/N)*t(z)%*%diag(bet)%*%x
	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		k <- ncol(x)
		
		G11 <- -(1/N)*t(z)%*%diag(bet)%*%x

		G12 <- (1/N)*t(x)%*%diag(bet*p[k])%*%z.in
		G12[k,] <- G12[k,]-(1/N)*apply(u*z.in,2,sum)

		G21 <- matrix(0,nrow=ncol(z.in),ncol=k)
		G22 <- -(1/N)*t(z.in)%*%z.in

		Gn <- rbind(cbind(G11,G12),cbind(G21,G22))
	}

	return(Gn)
}

frmhet.est <- function(type,x,z,link,start,Hy,variance,var.type,var.cluster,gixv,vhat,...)
{
	if(any(type==c("GMMxv","QMLxv")))
	{
		z.in <- z
		z <- x
	}

	GMM.est <- T
	if(any(type==c("GMMx","GMMxv")) & !any(Hy==0))
	{
		results <- tryCatch(glm(Hy ~ x-1,family=Gamma(link=log),maxit=100),error=function(e) return(NULL))
		if(any(is.null(results))) converged <- F
		else
		{
			p <- results$coefficients
			XB <- results$linear.predictors
			converged <- results$converged*(1-results$boundary)
		}

		if(converged==T) GMM.est <- F
	}
	if(GMM.est==T)
	{
		GMMn <- function(p)
		{
			XB <- as.vector(x%*%p)
			gi <- frmhet.gi(type,z,Hy,XB,link)$gi

			gn <- as.matrix(apply(gi,1,mean))
			Qn <- t(gn)%*%S%*%gn

			return(Qn)
		}

		S <- diag(ncol(z))
		results <- nlminb(start=start,objective=GMMn,...)
		p <- results$par
		XB <- as.vector(x%*%p)

		if(type=="GMMz" & ncol(z)>ncol(x))
		{
			fi.inv <- frmhet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,T,gixv,vhat)$fi.inv
			if(!is.character(fi.inv))
			{
				S <- fi.inv
				results <- nlminb(start=start,objective=GMMn,...)
				p <- results$par
				XB <- as.vector(x%*%p)
			}

			Qn <- results$objective
		}

		converged <- ifelse(results$convergence==0,T,F)
	}

	ret.list <- list(p=p,XB=XB,converged=converged)
	if(type=="GMMz" & ncol(z)>ncol(x)) ret.list[["Qn"]] <- Qn

	if(variance==F | converged==F) return(ret.list)

	if(any(type==c("GMMxv","QMLxv"))) z <- z.in
	p.var <- frmhet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,F,gixv,vhat)$p.var
	ret.list[["p.var"]] <- p.var

	return(ret.list)
}

frmhet.var <- function(type,p,XB,x,z,link,Hy,var.type,id,step.one,gixv,vhat)
{
	N <- length(Hy)

	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		z.in <- z
		z <- x
	}
	else z.in <- z

	if(any(type==c("LINx","LINxv","LINz")))
	{
		u <- Hy-XB
		gi <- t(z*u)
		bet <- rep(1,N)
	}
	if(any(type==c("GMMx","GMMz","GMMxv","QMLxv")))
	{
		results <- frmhet.gi(type,z,Hy,XB,link)
		gi <- results$gi
		u <- results$u

		if(type!="QMLxv")
		{
			G2 <- frmhet.links(link)$G2(XB)
			g2 <- frmhet.links(link)$g2(XB)
			bet <- Hy*g2/(G2^2)
		}
		if(type=="QMLxv")
		{
			yhat <- frmhet.links(link)$linkinv(XB)
			g <- frmhet.links(link)$mu.eta(XB)
			bet <- (g^2)/(yhat*(1-yhat))
		}
	}

	if(any(type==c("GMMxv","LINxv","QMLxv"))) gi <- rbind(gi,gixv)

	if(var.type=="robust") fi <- (1/N)*gi%*%t(gi)

	if(var.type=="cluster")
	{
		fi <- 0

		for(j in unique(id))
		{
			Zi <- matrix(z[id==j,],ncol=ncol(z))
			ui <- u[id==j]
			zu <- t(Zi)%*%ui

			if(any(type==c("GMMxv","LINxv","QMLxv")))
			{
				Zi <- matrix(z.in[id==j,],ncol=ncol(z.in))
				vi <- vhat[id==j]
				zu2 <- t(Zi)%*%vi
				zu <- rbind(zu,zu2)
			}

			fi <- fi+zu%*%t(zu)
		}

		fi <- fi/N
	}

	fi.inv <- tryCatch(solve(fi),error=function(e) NaN)
	if(any(is.nan(fi.inv))) fi.inv <- "singular"

	if(step.one==T) return(list(fi.inv=fi.inv))

	Gn <- frmhet.Gn(type,x,z,z.in,u,bet,p)

	if(is.numeric(fi.inv))
	{
		sigma <- N*t(Gn)%*%fi.inv%*%Gn
		p.var <- tryCatch(solve(sigma),error=function(e) NaN)
		if(any(is.nan(p.var))) p.var <- "singular"
	}
	else p.var <- "singular"

	if(is.character(p.var))
	{
		if(ncol(z)>ncol(x)) p.var <- "singular"
		else
		{
			Gn.inv <- tryCatch(solve(Gn),error=function(e) NaN)
			if(any(is.nan(Gn.inv))) p.var <- "singular"
			else p.var <- (1/N)*Gn.inv%*%fi%*%t(Gn.inv)
		}
	}

	ret.list <- list(p.var=p.var)

	return(ret.list)
}

frmhet.table <- function(p,p.var,x.names,type,link,converged,N,var.type,adjust,k,J,dfJ)
{
	if(converged==T)
	{
		stars <- rep("",length(p))

		if(!is.character(p.var))
		{
			p.sd <- diag(p.var)^0.5
			z.ratio <- p/p.sd
			p.value <- 2*(1-pnorm(abs(z.ratio)))

			stars[p.value<=0.01] <- "***"
			stars[p.value>0.01 & p.value<=0.05] <- "**"
			stars[p.value>0.05 & p.value<=0.1] <- "*"

			p.sd <- formatC(p.sd,digits=6,format="f")
			z.ratio <- formatC(z.ratio,digits=3,format="f")
			p.value <- formatC(p.value,digits=3,format="f")
			stars <- format(stars,justify="left")
		}
		else
		{
			p.sd <- rep(".",length(p))
			z.ratio <- rep(".",length(p))
			p.value <- rep(".",length(p))
		}

		p <- formatC(p,digits=6,format="f")
		results <- data.frame(cbind(p,p.sd,z.ratio,p.value,stars),row.names=NULL)

		namcol <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
		dimnames(results) <- list(x.names,namcol)

		cat("\n")
		cat("*** Fractional",link,"regression model ***")
		cat("\n")
		cat("*** Estimator:",type)
		if(adjust!=0)
		{
			cat("\n")
			if(is.numeric(adjust)) cat("*** Adjustment:",adjust,"added to all observations")
			else cat("*** Adjustment: all boundary observations dropped")
		}
		cat("\n\n")
		if(all(type!=c("GMMxv","LINxv","QMLxv"))) print(results)
		if(any(type==c("GMMxv","LINxv","QMLxv")))
		{
 			print(results[1:k,])
			cat("\n")
			cat("Reduced form:")
			cat("\n") 
			print(results[-(1:k),])
		}
		cat("\n")
		cat("Note:",var.type,"standard errors")
		cat("\n\n")
		cat("Number of observations:",N,"\n")
		cat("\n")
		if(any(type==c("GMMz","LINz")) & dfJ>0)
		{
			p.value <- 1-pchisq(J,dfJ)
			cat("J test of overidentifying moment conditions:",J,"(p-value:",p.value,")")
		}
	}
	else cat("ALGORITHM DID NOT CONVERGE")
	cat("\n")
}

frmhet.pe.var <- function(x,npar,which.x,x.names,xvar.names,type,p,xbhat,g,link,p.var,smearing,what,N)
{
	if(smearing==F) gd <- frmhet.links(link)$gd(xbhat)
	if(smearing==T) 
	{
		gd <- rep(NA,N)
		for(j in 1:N) gd[j] <- mean(frmhet.links(link)$gd(xbhat[j]+what))
	}

	PE.sd <- matrix(NA,nrow=npar,ncol=npar)

	for(i in 1:npar)
	{
		for(j in 1:npar) PE.sd[i,j] <- mean(p[i]*gd*x[,j]+(i==j)*g)
	}

	PE.sd <- PE.sd%*%p.var%*%t(PE.sd)

	PE.sd <- diag(PE.sd)^0.5
	if(any(x.names=="INTERCEPT")) PE.sd <- PE.sd[-1]

	names(PE.sd) <- xvar.names
	PE.sd <- PE.sd[which.x]

	return(PE.sd)
}

frmhet.pe.table <- function(PE.p,PE.sd,PE.type,which.x,xvar.names,title,adjust,at)
{
	z.ratio <- PE.p/PE.sd
	p.value <- 2*(1-pnorm(abs(z.ratio)))

	stars <- rep("",length(PE.p))
	stars[p.value<=0.01] <- "***"
	stars[p.value>0.01 & p.value<=0.05] <- "**"
	stars[p.value>0.05 & p.value<=0.1] <- "*"

	PE.p <- formatC(PE.p,digits=4,format="f")
	PE.sd <- formatC(PE.sd,digits=4,format="f")
	z.ratio <- formatC(z.ratio,digits=3,format="f")
	p.value <- formatC(p.value,digits=3,format="f")
	stars <- format(stars,justify="left")

	results <- data.frame(cbind(PE.p,PE.sd,z.ratio,p.value,stars),row.names=NULL)

	namcol <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
	dimnames(results) <- list(which.x,namcol)

	cat("\n\n")
	if(PE.type=="APE") cat("*** Average partial effects",title[1])
	if(PE.type=="CPE") cat("*** Conditional partial effects",title[1])
	cat("\n\n")
	cat(title[2])
	cat("\n")
	cat(title[3])
	cat("\n")
	if(adjust!=0) cat(title[4])
	cat("\n\n")
	print(results)
	cat("\n")
	if(PE.type=="CPE")
	{
		cat("------------------")

		if(length(at)==1)
		{
			if(any(at==c("mean","median"))) cat("\nNote: covariates evaluated at",at,"(or mode, for dummies) values\n")
			else
			{
				names(at) <- xvar.names
				cat("\nNote: covariates evaluated at the following values:\n\n")
				print(at)
			}
		}
		else
		{
			names(at) <- xvar.names
			cat("\nNote: covariates evaluated at the following values:\n\n")
			print(at)
		}
	}
}

frmhet.tests.table <- function(test.which,S,Sp,ver,title1,title2)
{
	stars <- rep("",length(S))
	stars[Sp<=0.01] <- "***"
	stars[Sp>0.01 & Sp<=0.05] <- "**"
	stars[Sp>0.05 & Sp<=0.1] <- "*"

	S <- formatC(S,digits=3,format="f")
	Sp <- formatC(Sp,digits=3,format="f")
	stars <- format(stars,justify="left")

	results <- data.frame(cbind(ver,S,Sp,stars))
	namcol <- c("Version","Statistic","p-value","")

	results <- results[-1,]
	dimnames(results) <- list(1:nrow(results),namcol)

	cat("\n")
	cat("***",test.which,"test ***")
	cat("\n")
	cat(title1)
	cat("\n\n")
	cat("H0: ",title2)
	cat("\n")
	cat("\n")
	print(results,row.names=F)
	cat("\n")
}

frmhet <- function(y,x,z=x,var.endog,start,type="GMMx",link="logit",intercept=T,table=T,variance=T,var.type="robust",var.cluster,adjust=0,...)
{
	### 1. Error and warning messages

	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")
	if(any(y>1) | any(y<0)) stop("The dependent variable has values outside the unit interval")
	if((any(y==1) | any(y==0)) & adjust==0 & any(type==c("LINx","LINz","LINxv"))) stop("0/1 values for the response variable: LIN estimators require adjustment")
	if(all(link!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(link)," - link not recognised")
	if(all(type!=c("GMMx","GMMz","GMMxv","LINx","LINz","LINxv","QMLxv"))) stop(sQuote(type)," - type not recognised")
	if(any(type==c("GMMx","GMMz","GMMxv")) & all(link!=c("logit","cloglog"))) stop("type and link not compatible")
	if(any(type==c("GMMx","GMMz","GMMxv")) & any(y==1)) stop("estimator does not allow y = 1")
	if(!is.numeric(adjust) & adjust!="drop") stop("adjust not defined properly")
	if(any(type==c("GMMxv","LINxv","QMLxv")) & missing(var.endog)) stop(sQuote(type)," requires var.endog to be specified")
	if(all(type!=c("GMMxv","LINxv","QMLxv")) & !missing(var.endog)) stop("var.endog should not be specified for this estimator")
	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
	}
	if(all(var.type!=c("robust","cluster"))) stop(sQuote(var.type)," - var.type not recognised")
	if(var.type=="cluster" & missing(var.cluster)) stop("option cluster for covariance matrix but no var.cluster supplied")
	if(var.type=="robust" & !missing(var.cluster)) stop("option robust for covariance matrix but var.cluster supplied")

	if(!is.logical(intercept)) stop("non-logical value assigned to option intercept")
	if(!is.logical(table)) stop("non-logical value assigned to option table")
	if(!is.logical(variance)) stop("non-logical value assigned to option variance")

	### 2. Data and variables preparation

	y <- as.vector(y)

	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.data.frame(z)) z <- as.matrix(z)

	if(!is.matrix(x)) stop("x is not a matrix")
	if(!is.matrix(z)) stop("z is not a matrix")

	x.names <- dimnames(x)[[2]]
	z.names <- dimnames(z)[[2]]

	if(is.null(x.names)) stop("x has no column names")
	if(is.null(z.names)) stop("z has no column names")

	if(intercept==T)
	{
		x <- cbind(1,x)
		z <- cbind(1,z)
		x.names <- c("INTERCEPT",x.names)
		z.names <- c("INTERCEPT",z.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(z.names)!=length(unique(z.names))) stop("some instrument names in z are identical")
	if(identical(x.names,z.names) & any(type==c("GMMz","GMMxv","LINz","LINxv","QMLxv"))) stop("instruments and covariates are identical")
	if(!identical(x.names,z.names) & any(type==c("GMMx","LINx"))) stop("instruments should not be specified for this estimator")
	if(length(x.names)>length(z.names)) stop("number of instruments not enough")

	if(adjust=="drop")
	{
		x <- as.matrix(x[y>0 & y<1,])
		z <- as.matrix(z[y>0 & y<1,])
		if(var.type=="cluster") var.cluster <- var.cluster[y>0 & y<1]
		if(!missing(var.endog)) var.endog <- var.endog[y>0 & y<1]
		y <- y[y>0 & y<1]
	}

	if(any(length(y)!=c(nrow(x),nrow(z)))) stop("the number of observations for y, x and/or z is different")
	if(var.type=="cluster")
	{
		if(length(y)!=length(var.cluster)) stop("var.cluster does not have the appropriate dimension")
	}
	if(!missing(var.endog))
	{
		if(length(y)!=length(var.endog)) stop("var.endog does not have the appropriate dimension")
	}

	class <- "frmhet"

	### 3. Estimation

	if(any(type==c("GMMx","LINx"))) z <- x

	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		results <- lm(var.endog ~ 0+z)
		PIhat <- results$coefficients
		vhat <- results$residuals

		gixv <- t(z*vhat)
		x <- cbind(x,vhat)
	}
	else
	{
		gixv <- NA
		vhat <- NA
	}

	k <- ncol(x)
	kz <- ncol(z)
	N <- nrow(x)
	J <- NA
	dfJ <- kz-k

	if(any(type==c("GMMx","GMMz","GMMxv","QMLxv")))
	{
		if(any(type==c("GMMx","GMMz","GMMxv"))) Hy <- frmhet.links(link)$H1(y)
		if(type=="QMLxv") Hy <- y

		if(missing(start)) start <- rep(0,k)
		if(length(start)!=k) stop("start is not of the same dimension as the covariate vector (including vhat, in case of GMMxv or QMLxv)")
		
		results <- frmhet.est(type,x,z,link,start,Hy,variance,var.type,var.cluster,gixv,vhat,...)
		p <- results$p
		converged <- results$converged

		XB <- results$XB
		if(any(type=="GMMz") & dfJ>0)
		{
			Qn <- results$Qn
			J <- N*Qn
		}
	}

	if(any(type==c("LINx","LINz","LINxv")))
	{
		if(is.numeric(adjust))
		{
			y <- y+adjust
			if(any(y>=1) | any(y<=0)) stop("After adjustment, the dependent variable has values outside or at the boundaries of the unit interval")
		}

		Hy <- frmhet.links(link)$linkfun(y)

		if(any(type==c("LINx","LINxv")))
		{
			results <- lm(Hy ~ 0+x)
			p <- results$coefficients

			XB <- results$fitted.values
		}
		if(type==c("LINz"))
		{
			XZ <- t(x)%*%z
			p <- solve(XZ%*%t(XZ))%*%(XZ%*%(t(z)%*%Hy))
			if(kz>k)
			{
				u <- as.vector(Hy-x%*%p)

				gi <- t(z*u)
				fi <- (1/N)*gi%*%t(gi)
				fi.inv <- solve(fi)
 				p <- as.vector(solve(XZ%*%fi.inv%*%t(XZ))%*%(XZ%*%fi.inv%*%(t(z)%*%Hy)))
			}

			XB <- as.vector(x%*%p)

			if(kz>k) J <- N*t(Hy-XB)%*%z%*%fi.inv%*%t(z)%*%(Hy-XB)
		}

		converged <- T
		if(variance==T) results <- frmhet.var(type,p,XB,x,z,link,Hy,var.type,var.cluster,F,gixv,vhat)
	}

	if(variance==T & converged==T) p.var <- results$p.var
	else p.var <- "singular"

	x.names.in <- x.names

	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		p <- c(p,PIhat)
		x.names <- c(x.names,"vhat",paste("Z",z.names,sep="_"))
	}
	names(p) <- x.names

	if(table==T) frmhet.table(p,p.var,x.names,type,link,converged,N,var.type,adjust,k,J,dfJ)

	formula <- y ~ x - 1

	res <- list(class=class,formula=formula,type=type,link=link,adjust=adjust,p=p,Hy=Hy,xbhat=as.vector(XB),converged=converged,x.names=x.names.in)
	if(any(type==c("GMMz","LINz")) & kz>k) res[["J"]] <- J

	if(variance==T & converged==T)
	{
		if(is.character(p.var)) p.var <- matrix(NA,nrow=length(p),ncol=length(p))
		dimnames(p.var) <- list(x.names,x.names)
		res[["p.var"]] <- p.var
		res[["var.type"]] <- var.type
		if(var.type=="cluster") res[["var.cluster"]] <- var.cluster
	}

	### 4. Return results

	return(invisible(res))
}

frmhet.reset <- function(object,lastpower.vec=3,version="Wald",table=T,...)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an frmhet command")
	if(object$class!="frmhet") stop("object is not the output of an frmhet command")
	if(object$converged==0) stop("object is not the output of a successful (converged) frmhet command")
	if(any(lastpower.vec<2)) stop(sQuote(lastpower.vec)," - lastpower.vec contains elements lower than 2")
	if(all(object$type!=c("GMMx","LINx"))) stop("frmhet.reset is only implemented for GMMx and LINx estimators")
	if(all(version!="LM") & all(version!="Wald")) stop("test version not correctly specified")
	if(any(version=="LM") & object$type=="LINx") stop("LM version not implemented for LINx; choose Wald version")
	if(is.null(object$var.type)) stop("frmhet command was run with variance = F")

	### 2. Recovering definitions and estimates

	mf <- model.frame(object$formula)
	y <- model.response(mf)
	x <- model.matrix(object$formula)

	type <- object$type
	xbhat <- object$xbhat
	Hy <- object$Hy
	link <- object$link
	adjust <- object$adjust

	if(any(version=="Wald") | type=="GMMx")
	{
		var.type <- object$var.type
		if(var.type=="cluster") var.cluster <- object$var.cluster
	}

	title1 <- paste("Fractional",link,"regression model")
	title2 <- paste("Estimator:",type)
	if(adjust!=0)
	{
		if(is.numeric(adjust)) title3 <- paste("(adjustment:",adjust,"added to all observations)")
		else title3 <- "(adjustment: all boundary observations dropped)"
		title2 <- paste(title2,title3,sep=" ")
	}

	### 3. Test

	lastpower.vec <- round(lastpower.vec,0)

	N <- length(Hy)
	g <- frmhet.links(link)$mu.eta(xbhat)
	gx <- g*x

	xx.all <- as.matrix(xbhat^2)
	if(max(lastpower.vec)>2) for(i in 3:max(lastpower.vec)) xx.all <- cbind(xx.all,xbhat^i)

	ver <- NA
	S <- NA
	Sp <- NA

	for(m in lastpower.vec)
	{
		df <- m-1
		xx <- xx.all[,1:(m-1)]

		X <- cbind(x,xx)
		k <- ncol(X)

		if(any(version=="LM"))
		{
			name <- paste("LM(",m,")",sep="")
			ver <- c(ver,name)

			if(type=="GMMx")
			{
				gi <- frmhet.gi(type,X,Hy,xbhat,link)$gi
				gn <- as.matrix(apply(gi,1,mean))

				if(var.type=="robust") fi <- (1/N)*gi%*%t(gi)
				if(var.type=="cluster")
				{
					fi <- 0
					id <- var.cluster
					u <- Hy-xbhat

					for(j in unique(id))
					{
						Zi <- matrix(X[id==j,],ncol=k)
						ui <- u[id==j]
						zu <- t(Zi)%*%ui

						fi <- fi+zu%*%t(zu)
					}

					fi <- fi/N
				}

				fi.inv <- tryCatch(solve(fi),error=function(e) NaN)
				if(any(is.nan(fi.inv))) fi.inv <- "singular"

				if(!is.character(fi.inv))
				{
					Sj <- N*t(gn)%*%fi.inv%*%gn
					S <- c(S,Sj)
					Sp <- c(Sp,1-pchisq(Sj,df))
				}
				else
				{
					S <- c(S,NA)
					Sp <- c(Sp,NA)
				}
			}
		}

		if(any(version=="Wald"))
		{
			name <- paste("Wald(",m,")",sep="")
			ver <- c(ver,name)

			if(type=="GMMx")
			{
				results <- frmhet.est(type,X,X,link,c(object$p,rep(0,df)),Hy,T,var.type,var.cluster,NA,NA,...)
				converged <- results$converged
			}
			if(type=="LINx")
			{
				results <- lm(Hy ~ 0+X)
				converged <- T
			}

			if(converged==T)
			{
				if(type=="GMMx")
				{
					p1 <- results$p
					p.var1 <- results$p.var
				}
				if(type=="LINx")
				{ 
					p1 <- results$coefficients
					XB <- results$fitted.values
					p.var1 <- frmhet.var(type,p1,XB,X,X,link,Hy,var.type,var.cluster,F,NA,NA)$p.var
				}

				if(!is.character(p.var1))
				{
					p.n <- length(p1)
					p1 <- p1[(p.n-df+1):p.n]
					p.var1 <- p.var1[(p.n-df+1):p.n,(p.n-df+1):p.n]

					Sj <- t(p1)%*%solve(p.var1)%*%p1
					S <- c(S,Sj)
					Sp <- c(Sp,1-pchisq(Sj,df))
				}
				else
				{
					S <- c(S,NA)
					Sp <- c(Sp,NA)
				}
			}
			else
			{
				S <- c(S,NA)
				Sp <- c(Sp,NA)
			}
		}
	}

	if(any(!is.na(S)) & table==T) frmhet.tests.table("RESET",S,Sp,ver,title1,title2)
	if(all(is.na(S))) cat("RESET test could not be computed; either algorithm did not converge (Wald version) or covariance matrix is singular (Wald/LM versions)\n")

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- ver[-1]

	return(invisible(statistics))
}

frmhet.pe <- function(object,smearing=T,APE=T,CPE=F,at=NULL,which.x=NULL,table=T,variance=T)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an frmhet command")
	if(object$class!="frmhet") stop("object is not the output of an frmhet command")

	if(all(c(APE,CPE)==F)) stop("You must specify at least one option: APE and/or CPE")
	if(CPE==F & !is.null(at)) stop("option at is only required for cpe")

	if(object$converged==0) stop("object is not the output of a successful (converged) frmhet command")

	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
	}

	### 2. Recovering definitions and estimates and other definitions

	type <- object$type
	link <- object$link
	adjust <- object$adjust
	p <- object$p
	p.var <- object$p.var
	x <- model.matrix(object$formula)
	x.names <- object$x.names
	Hy <- object$Hy

	if(any(x.names=="INTERCEPT")) xvar.names <- x.names[-1]
	else xvar.names <- x.names

	k <- length(xvar.names)
	npar <- ncol(x)
	n <- nrow(x)

	if(any(type==c("QMLxv"))) pv <- p[npar]
	if(any(type==c("GMMxv","LINxv","QMLxv")))
	{
		npar <- npar-1
		p <- p[1:npar]
		p.var <- p.var[1:npar,1:npar]
	}

	if(all(type!=c("GMMxv","LINxv","QMLxv"))) xbhat <- object$xbhat
	if(any(type==c("GMMxv","LINxv","QMLxv"))) xbhat <- as.vector(x[,-(npar+1)]%*%as.matrix(p))

	if(is.null(which.x)) which.x <- xvar.names
	xw.names <- unique(c(xvar.names,which.x))
	if(!identical(xvar.names,xw.names)) stop("option which not appropriately defined")

	### 3. Average partial effects

	if(smearing==T) title1 <- "(conditional only on unobservables, based on the smearing estimator)"
	else title1 <- "(conditional on both observables and unobservables, with error term = 0)"

	title2 <- paste("Fractional",link,"regression model")
	title3 <- paste("Estimator:",type)
	title <- c(title1,title2,title3)
	if(adjust!=0)
	{
		if(is.numeric(adjust)) title4 <- paste("Adjustment:",adjust,"added to all observations")
		else title4 <- "Adjustment: all boundary observations dropped"
		title <- c(title,title4)
	}

	if(APE==T)
	{
		PE.type <- "APE"

		if(any(x.names=="INTERCEPT")) p.pe <- matrix(rep(p[-1],each=n),ncol=k)
 		else p.pe <- matrix(rep(p,each=n),ncol=k)
		dimnames(p.pe) <- list(NULL,xvar.names)

		if(smearing==F)
		{
			g <- frmhet.links(link)$mu.eta(xbhat)
			what <- NA
		}
		if(smearing==T)
		{
			if(type=="QMLxv") what <- x[,npar+1]*pv

			if(any(type==c("LINx","LINxv","LINz"))) what <- Hy-xbhat
			if(any(type==c("GMMx","GMMxv","GMMz")))
			{
				G2 <- frmhet.links(link)$G2(xbhat)
				uhat.star <- Hy/G2-1
				what <- frmhet.links(link)$H2(uhat.star+1)
			}

			g <- rep(NA,n)
			for(j in 1:n) g[j] <- mean(frmhet.links(link)$mu.eta(xbhat[j]+what))
		}

		PE.p <- as.matrix(p.pe[,which.x])*g
		PE.p <- apply(PE.p,2,mean)

		resAPE <- list(PE.p=PE.p)

		if(variance==T)
		{
			PE.sd <- frmhet.pe.var(x,npar,which.x,x.names,xvar.names,type,p,xbhat,g,link,p.var,smearing,what,n)
			resAPE[["PE.sd"]] <- PE.sd
		}

		if(table==T) frmhet.pe.table(PE.p,PE.sd,PE.type,which.x,xvar.names,title,adjust,at)
	}

	### 4. Conditional partial effects

	if(CPE==T)
	{
		PE.type <- "CPE"

		if(is.null(at)) at <- "mean"

		if(length(at)==1)
		{
			if((!any(at==c("mean","median")) & k!=1)) stop("at not appropriately specified")

			if(any(at==c("mean","median")))
			{
				if(at=="mean") xm <- apply(x,2,mean)
				if(at=="median") xm <- apply(x,2,median)

				xdum <- apply(x,2,function(a) all(a %in% c(0,1)))

				if(any(type==c("GMMxv","LINxv","QMLxv")))
				{
					xm <- xm[-npar]
					xdum <- xdum[-npar]
				}
				xm[xdum==T] <- round(xm,0)[xdum==T]
			}
			else
			{
				if(is.numeric(at))
				{
					if(any(x.names=="INTERCEPT")) xm <- c(1,at)
					else xm <- at
				}
				else stop("at not appropriately specified")
			}
		}
		else
		{
			if(length(at)!=k) stop("at not appropriately specified")
			else
			{
				if(any(x.names=="INTERCEPT")) xm <- c(1,at)
				else xm <- at
			}
		}

		xmbhat <- as.vector(xm%*%p)

		if(smearing==F) g <- frmhet.links(link)$mu.eta(xmbhat)

		if(smearing==T)
		{
			if(type=="QMLxv") what <- x[,npar]*pv
			if(any(type==c("LINx","LINxv","LINz"))) what <- Hy-xbhat
			if(any(type==c("GMMx","GMMxv","GMMz")))
			{
				G2 <- frmhet.links(link)$G2(xbhat)
				uhat.star <- Hy/G2-1
				what <- frmhet.links(link)$H2(uhat.star+1)
			}

			g <- mean(frmhet.links(link)$mu.eta(xmbhat+what))
		}

		PE.p <- p[which.x]*g

		resCPE <- list(PE.p=PE.p)

		if(variance==T)
		{
			PE.sd <- frmhet.pe.var(x,npar,which.x,x.names,xvar.names,type,p,xbhat,g,link,p.var,smearing,what,n)
			resCPE[["PE.sd"]] <- PE.sd
		}

		if(table==T) frmhet.pe.table(PE.p,PE.sd,PE.type,which.x,xvar.names,title,adjust,at)
	}

	### 5. Return results

	if(APE==T & CPE==T) return(invisible(list(ape=resAPE,cpe=resCPE)))
	if(APE==T & CPE==F) return(invisible(resAPE))
	if(APE==F & CPE==T) return(invisible(resCPE))
}
