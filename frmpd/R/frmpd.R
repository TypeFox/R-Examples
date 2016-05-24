frmpd.links <- function(link) 
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
		}
	)

	structure(list(linkfun=linkfun,linkinv=linkinv,mu.eta=mu.eta,H1=H1,G2=G2,g2=g2,H2=H2,gd=gd,valideta=valideta,name=link),class="link-glm")
}

frmpd.gi <- function(type,id,Ti,Hy,z,XB,link,at,at1)
{
	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpfe")))
	{
		XBexp <- exp(XB)

		if(type=="GMMc") u <- (XBexp[at1]/XBexp[at])*Hy[at]-Hy[at1]
		if(type=="GMMww")
		{
			u <- Hy/XBexp
			u <- u[at]-u[at1]
		}
		if(type=="GMMbgw")
		{
			Hy.m <- rep(as.vector(by(Hy[at],id,mean)),times=Ti)
			XBexp.m <- rep(as.vector(by(XBexp[at],id,mean)),times=Ti)
			u <- Hy[at]-(XBexp[at]/XBexp.m)*Hy.m
		}
			if(type=="GMMpfe")
		{
			HyXBexp.m <- rep(as.vector(by(Hy[at]/XBexp[at],id,mean)),times=Ti)
			u <- Hy[at]/(XBexp[at]*HyXBexp.m)-1
		}
	}
	if(any(type==c("GMMpre","GMMcre")))
	{
		G2 <- frmpd.links(link)$G2(XB)
		u <- Hy[at]/G2[at]-1
	}
	if(type=="QMLcre")
	{
		yhat <- frmpd.links(link)$linkinv(XB)
		g <- frmpd.links(link)$mu.eta(XB)
		u <- (Hy-yhat)*g/(yhat*(1-yhat))
	}

	gi <- t(z*u)

	return(list(gi=gi,u=u))
}

frmpd.Gn <- function(type,x.exogenous,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,p,z.in,u)
{
	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpfe")))
	{
		XBexp <- exp(XB)

		if(type=="GMMc") Gn <- t(z)%*%diag((XBexp[at1]/XBexp[at])*Hy[at])%*%(x[at1,]-x[at,])
		if(type=="GMMww")
		{
			uu <- -Hy/XBexp
			Gn <- t(z)%*%(diag(uu[at])%*%x[at,]-diag(uu[at1])%*%x[at1,])
		}
		if(type=="GMMbgw")
		{
			Hy.m <- rep(as.vector(by(Hy[at],id,mean)),times=Ti)
			XBexp.m <- rep(as.vector(by(XBexp[at],id,mean)),times=Ti)

			XBexpX.m <- matrix(NA,nrow=NT,ncol=k)
			for(j in 1:k) XBexpX.m[,j] <- rep(as.vector(by(XBexp[at]*x[at,j],id,mean)),times=Ti)

			Gn <- -t(z)%*%diag(XBexp[at]*Hy.m/XBexp.m)%*%x[at,]+t(z)%*%((XBexp[at]*Hy.m/((XBexp.m)^2))*XBexpX.m)
		}
		if(type=="GMMpfe")
		{
			HyXBexp.m <- rep(as.vector(by(Hy[at]/XBexp[at],id,mean)),times=Ti)

			HyXBexpX.m <- matrix(NA,nrow=NT,ncol=k)
			for(j in 1:k) HyXBexpX.m[,j] <- rep(as.vector(by((Hy[at]/XBexp[at])*x[at,j],id,mean)),times=Ti)
			Gn <- -t(z)%*%diag(Hy[at]/(XBexp[at]*HyXBexp.m))%*%x[at,]+t(z)%*%(Hy[at]/(XBexp[at]*(HyXBexp.m^2))*HyXBexpX.m)
		}
	}
	if(any(type==c("GMMpre","GMMcre")))
	{
		G2 <- frmpd.links(link)$G2(XB)
		g2 <- frmpd.links(link)$g2(XB)
		bet <- Hy*g2/(G2^2)
		Gn <- -t(z)%*%diag(bet[at])%*%x[at,]
	}
	if(type=="QMLcre")
	{
		yhat <- frmpd.links(link)$linkinv(XB)
		g <- frmpd.links(link)$mu.eta(XB)
		bet <- g^2/(yhat*(1-yhat))
		Gn <- -t(z)%*%diag(bet)%*%x

		if(x.exogenous==F)
		{
			G12 <- t(x)%*%diag(bet*p[k])%*%z.in
			G12[k,] <- G12[k,]-apply(u*z.in,2,sum)

			G21 <- matrix(0,nrow=ncol(z.in),ncol=k)
			G22 <- -t(z.in)%*%z.in

			Gn <- rbind(cbind(Gn,G12),cbind(G21,G22))
		}
	}

	Gn <- Gn/NT

	return(Gn)
}

frmpd.est <- function(type,x.exogenous,lags,id,Ti,Hy,x,z,link,var.type,start,at,at1,variance,NT,k,kz,gixv,vhat,bootstrap,...)
{
	if(type=="QMLcre" & x.exogenous==F)
	{
		z.in <- z
		z <- x
		kz <- ncol(z)
	}

	GMM.est <- T
	if(x.exogenous==T & ((any(type==c("GMMpre","GMMcre")) & !any(Hy==0) & lags==F) | type=="QMLcre"))
	{
		if(type=="QMLcre") results <- tryCatch(glm(Hy ~ x-1,family=quasibinomial(link=frmpd.links(link)),maxit=100),error=function(e) return(NULL))
		if(any(type==c("GMMpre","GMMcre"))) results <- tryCatch(glm(Hy ~ x-1,family=Gamma(link=log),maxit=100),error=function(e) return(NULL))

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
			gi <- frmpd.gi(type,id,Ti,Hy,z,XB,link,at,at1)$gi
			gn <- as.matrix(apply(gi,1,mean))
			Qn <- t(gn)%*%S%*%gn

			return(Qn)
		}

		S <- diag(kz)
		results <- nlminb(start=start,objective=GMMn,...)
		p <- results$par
		XB <- as.vector(x%*%p)

		if(type!="QMLcre" & kz>k)
		{
			fi.inv <- frmpd.var(type,x.exogenous,var.type,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,kz,T,gixv,vhat,p)$fi.inv
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

	ret.list <- list(p=p,converged=converged)
	if(type!="QMLcre" & kz>k) ret.list[["Qn"]] <- Qn

	if(variance==F | converged==F | bootstrap==T) return(ret.list)

	if(type=="QMLcre" & x.exogenous==F) z <- z.in
	p.var <- frmpd.var(type,x.exogenous,var.type,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,kz,F,gixv,vhat,p)$p.var
	ret.list[["p.var"]] <- p.var

	return(ret.list)
}

frmpd.var <- function(type,x.exogenous,var.type,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,kz,step.one,gixv,vhat,p)
{
	if(type=="QMLcre" & x.exogenous==F)
	{
		z.in <- z
		z <- x
		kz <- ncol(z)
	}
	else z.in <- z

	results <- frmpd.gi(type,id,Ti,Hy,z,XB,link,at,at1)
	u <- results$u

	if(var.type=="robust")
	{
		gi <- results$gi
		if(type=="QMLcre" & x.exogenous==F) gi <- rbind(gi,gixv)
		fi <- (1/NT)*gi%*%t(gi)
	}

	if(var.type=="cluster")
	{
		fi <- 0

		for(j in unique(id))
		{
			Zi <- matrix(z[id==j,],ncol=kz)
			ui <- u[id==j]
			zu <- t(Zi)%*%ui

			if(type=="QMLcre" & x.exogenous==F)
			{
				Zi <- matrix(z.in[id==j,],ncol=ncol(z.in))
				vi <- vhat[id==j]
				zu2 <- t(Zi)%*%vi
				zu <- rbind(zu,zu2)
			}

			fi <- fi+zu%*%t(zu)
		}

		fi <- fi/NT
	}

	fi.inv <- tryCatch(solve(fi),error=function(e) NaN)
	if(any(is.nan(fi.inv))) fi.inv <- "singular"

	if(step.one==T) return(list(fi.inv=fi.inv))

	Gn <- frmpd.Gn(type,x.exogenous,id,Ti,Hy,x,z,XB,link,at,at1,NT,k,p,z.in,u)

	if(is.numeric(fi.inv))
	{
		sigma <- NT*t(Gn)%*%fi.inv%*%Gn
		p.var <- tryCatch(solve(sigma),error=function(e) NaN)
		if(any(is.nan(p.var))) p.var <- "singular"
	}
	else p.var <- "singular"

	if(is.character(p.var))
	{
		if(k!=kz) p.var <- "singular"
		else
		{
			Gn.inv <- tryCatch(solve(Gn),error=function(e) NaN)
			if(any(is.nan(Gn.inv))) p.var <- "singular"
			else p.var <- (1/NT)*Gn.inv%*%fi%*%t(Gn.inv)
		}
	}

	ret.list <- list(p.var=p.var)

	return(ret.list)
}

frmpd.table <- function(p,p.var,x.names,x.exogenous,lags,type,link,converged,N.ini,N,NT.ini,NT,J,dfJ,k,var.type,bootstrap)
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
		cat("\n")
		cat("*** Exogeneity:",x.exogenous)
		cat("\n")
		cat("*** Use first lag of instruments:",lags)
		cat("\n\n")

		if(type!="QMLcre" | x.exogenous==T) print(results)
		else
		{
 			print(results[1:k,])
			cat("\n")
			cat("Reduced form:")
			cat("\n") 
			print(results[-(1:k),])
		}

		cat("\n")
		if(bootstrap==F) cat("Note:",var.type,"standard errors")
		if(bootstrap==T) cat("Note: bootstrap standard errors")
		cat("\n\n")
		cat("Number of observations (initial):",NT.ini,"\n")
		cat("Number of observations (for estimation):",NT,"\n")
		cat("Number of cross-sectional units (initial):",N.ini,"\n")
		cat("Number of cross-sectional units (for estimation):",N,"\n")
		cat("Average number of time periods per cross-sectional unit (initial):",NT.ini/N.ini,"\n")
		cat("Average number of time periods per cross-sectional unit (for estimation):",NT/N,"\n")
		cat("\n")
		if(type!="QMLcre" & dfJ>0)
		{
			p.value <- 1-pchisq(J,dfJ)
			cat("J test of overidentifying moment conditions:",J,"(p-value:",p.value,")")
		}
	}
	else cat("ALGORITHM DID NOT CONVERGE")
	cat("\n")
}

frmpd <- function(id,time,y,x,z,var.endog,x.exogenous=T,lags,start,type,GMMww.cor=T,link="logit",intercept=T,table=T,variance=T,var.type="cluster",tdummies=F,bootstrap=F,B=200,...)
{
	### 1. Error and warning messages

	if(missing(id)) stop("variable id is missing")
	if(missing(time)) stop("variable time is missing")
	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")

	if(missing(type)) stop("type must be specified")
	if(any(y>1) | any(y<0)) stop("The dependent variable has values outside the unit interval")
	if(!is.logical(x.exogenous)) stop("x.exogenous must be logical")
	if(!is.logical(GMMww.cor)) stop("GMMww.cor must be logical")
	if(!is.logical(intercept)) stop("intercept must be logical")
	if(!is.logical(table)) stop("table must be logical")
	if(!is.logical(variance)) stop("variance must be logical")
	if(!is.logical(tdummies)) stop("tdummies must be logical")
	if(all(link!=c("logit","probit","cloglog"))) stop(sQuote(link)," - link not recognised")
	if(all(type!=c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre","QMLcre"))) stop(sQuote(type)," - type not recognised")
	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre")) & all(link!=c("logit","cloglog"))) stop("type and link not compatible")
	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre")) & any(y==1)) stop("estimator does not allow y = 1")
	if(type=="QMLcre" & link!="probit") stop("type and link not compatible")
	if(!missing(z) & x.exogenous==T) stop("z must not be specified in case of exogeneity")
	if(missing(z) & x.exogenous==F) stop("z needs to be specified in case of endogeneity")
	if(type=="QMLcre" & x.exogenous==F & missing(var.endog)) stop("QMLcre under endogeneity requires var.endog to be specified")
	if((type!="QMLcre" | x.exogenous==T) & !missing(var.endog)) stop("var.endog should not be specified for this estimator")
	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
	}
	if(all(var.type!=c("robust","cluster"))) stop(sQuote(var.type)," - var.type not recognised")

	if(any(unlist(by(time,id,duplicated)))) stop("variable time has duplicated entries for the same id")

	### 2. Data and variables preparation - general

	id <- as.vector(id)
	time <- as.vector(time)
	y <- as.vector(y)
	if(missing(z)) z <- x

	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.data.frame(z)) z <- as.matrix(z)

	if(!is.matrix(x)) stop("x is not a matrix")
	if(!is.matrix(z)) stop("z is not a matrix")

	if(any(is.na(y))) stop("y has missing values")
	if(any(is.na(x))) stop("x has missing values")
	if(any(is.na(z))) stop("z has missing values")

	x.names <- dimnames(x)[[2]]
	z.names <- dimnames(z)[[2]]

	if(is.null(x.names)) stop("x has no column names")
	if(is.null(z.names)) stop("z has no column names")

	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpfe","GMMcre","QMLcre"))) intercept <- F
	if(intercept==T)
	{
		x <- cbind(1,x)
		z <- cbind(1,z)
		x.names <- c("INTERCEPT",x.names)
		z.names <- c("INTERCEPT",z.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(z.names)!=length(unique(z.names))) stop("some instrument names in z are identical")
	if(identical(x.names,z.names) & x.exogenous==F) stop("instruments and covariates are identical")
	if(!identical(x.names,z.names) & x.exogenous==T) stop("instruments should not be specified for this estimator")
	if(length(x.names)>length(z.names)) stop("number of instruments not enough")

	if(any(length(y)!=c(nrow(x),nrow(z)))) stop("the number of observations for y, x and/or z is different")
	if(!missing(var.endog))
	{
		if(length(y)!=length(var.endog)) stop("var.endog does not have the appropriate dimension")
	}

	### 3. Data and variables preparation - panel / estimator specifics

	if(missing(lags))
	{
		if(any(type==c("GMMc","GMMww"))) lags <- T
		if(any(type==c("GMMbgw","GMMpre","GMMpfe","GMMcre","QMLcre"))) lags <- F
	}
	if(!is.logical(lags)) stop("lags must be logical")
	if(any(type==c("GMMcre","QMLcre")) & lags==T) stop("GMMcre/QMLcre cannot be used with lagged instruments")

	N.ini <- length(unique(id))
	Ti.ini <- as.vector(by(id,id,length))
	NT.ini <- length(id)

	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpfe","GMMcre","QMLcre")))
	{
		for(j in 1:ncol(x))
		{
			xm <- rep(as.vector(by(x[,j],id,mean)),times=Ti.ini)
			x.dif <- x[,j]-xm

			if(all(x.dif==0)) stop("Time-invariant covariates not allowed")

			if(any(type==c("GMMcre","QMLcre")) & x.exogenous==T)
			{
				if(j==1) x.m <- cbind(xm)
				if(j>1) x.m <- cbind(x.m,xm)
			}
		}
		for(j in 1:ncol(z))
		{
			zm <- rep(as.vector(by(z[,j],id,mean)),times=Ti.ini)
			z.dif <- z[,j]-zm

			if(all(z.dif==0)) stop("Time-invariant instruments not allowed")

			if(any(type==c("GMMcre","QMLcre")) & x.exogenous==F)
			{
				if(j==1) z.m <- cbind(zm)
				if(j>1) z.m <- cbind(z.m,zm)
			}
		}
	}

	T.ord <- rep(NA,NT.ini)

	if(tdummies==T)
	{
		if(lags==T) out <- 2
		else out <- 1

		td <- matrix(0,nrow=NT.ini,ncol=length(unique(time))-out)
	}

	a <- 0
	for(j in min(time):max(time))
	{
		a <- a+1
		T.ord[time==j] <- a
		if(tdummies==T)
		{
			if(a>out)
			{
				td[time==j,a-out] <- 1
				if(a==out+1) td.names <- paste("time",j,sep=".")
				if(a>(out+1)) td.names <- c(td.names,paste("time",j,sep="."))
			}
		}
	}

	ord <- order(id,time)
	id <- id[ord]
	T.ord <- T.ord[ord]
	y <- y[ord]
	x <- cbind(x[ord,])
	if(tdummies==T) td <- cbind(td[ord,])
	z <- cbind(z[ord,])

	if(any(type==c("GMMc","GMMww")) | lags==T)
	{
		at <- rep(F,NT.ini)
		at1 <- rep(F,NT.ini)

		if(any(Ti.ini!=max(Ti.ini)))
		{
			keep <- rep(F,NT.ini)

			for(j in 1:NT.ini)
			{
				if(j==1 & (T.ord[2]-T.ord[1])==1 & id[1]==id[2]) at1[1] <- T
				if(j>1 & j<NT.ini)
				{
					if((T.ord[j]-T.ord[j-1])==1 & id[j]==id[j-1]) at[j] <- T
					if((T.ord[j+1]-T.ord[j])==1 & id[j]==id[j+1]) at1[j] <- T
				}
				if(j==NT.ini & (T.ord[NT.ini]-T.ord[NT.ini-1])==1 & id[NT.ini]==id[NT.ini-1]) at[NT.ini] <- T
			}

			keep <- at==T | at1==T
			id <- id[keep==T]
			y <- y[keep==T]
			x <- cbind(x[keep==T,])
			if(tdummies==T) td <- cbind(td[keep==T,])
			z <- cbind(z[keep==T,])
			at <- at[keep==T]
			at1 <- at1[keep==T]
		}
		else
		{
			at[T.ord!=1] <- T
			at1[T.ord!=max(Ti.ini)] <- T			
		}

		z.full <- z
		id.full <- id

		if(lags==T)
		{
			z <- cbind(z[at1,])
			id <- id[at1]
		}
		if(lags==F)
		{
			z <- cbind(z[at,])
			id <- id[at]
		}
	}
	else
	{
		at <- rep(T,NT.ini)
		at1 <- rep(T,NT.ini)

		z.full <- z
		id.full <- id
	}

	if(tdummies==T)
	{
		x <- cbind(x,td)
		z <- cbind(z,td[at,])
		z.full <- cbind(z.full,td)

		if(any(type==c("GMMcre","QMLcre")))
		{
			xa.names <- x.names
			za.names <- z.names
		}

		x.names <- c(x.names,td.names)
		z.names <- c(z.names,td.names)
	}
	else
	{
		xa.names <- x.names
		za.names <- z.names
	}

	if(any(type==c("GMMcre","QMLcre")))
	{
		if(any(Ti.ini!=max(Ti.ini)))
		{
			Ti.unique <- unique(Ti.ini)
			Ti.unique <- Ti.unique[Ti.unique!=1]
			Ti.unique <- sort(Ti.unique)

			T.dum <- matrix(NA,nrow=NT.ini,ncol=length(Ti.unique))
			a <- 0

			Ti.ini.exp <- rep(Ti.ini,Ti.ini)

			for(j in Ti.unique)
			{
				a <- a+1
				T.dum[,a] <- Ti.ini.exp==j
			}

			if(x.exogenous==T)
			{
				x.m.aug <- T.dum*x.m[,1]
				if(ncol(x.m)>1) for(j in 2:ncol(x.m)) x.m.aug <- cbind(x.m.aug,T.dum*x.m[,j])
				x.m <- cbind(T.dum,x.m.aug)
				x.m.names <- c(paste("INTERCEPT",Ti.unique,sep="_"),paste(rep(paste(xa.names,"mean",sep="_"),each=length(Ti.unique)),rep(Ti.unique,length(xa.names)),sep="_"))
				x.m <- cbind(x.m[Ti.ini.exp!=1,])
			}
			if(x.exogenous==F)
			{
				z.m.aug <- T.dum*z.m[,1]
				if(ncol(z.m)>1) for(j in 2:ncol(z.m)) z.m.aug <- cbind(z.m.aug,T.dum*z.m[,j])
				z.m <- cbind(T.dum,z.m.aug)
				z.m.names <- c(paste("INTERCEPT",Ti.unique,sep="_"),paste(rep(paste(za.names,"mean",sep="_"),each=length(Ti.unique)),rep(Ti.unique,length(za.names)),sep="_"))
				z.m <- cbind(z.m[Ti.ini.exp!=1,])
			}

			id <- id[Ti.ini.exp!=1]
			y <- y[Ti.ini.exp!=1]
			x <- cbind(x[Ti.ini.exp!=1,])
			z <- cbind(z[Ti.ini.exp!=1,])
			at <- at[Ti.ini.exp!=1]
			if(type=="QMLcre" & x.exogenous==F) var.endog <- var.endog[Ti.ini.exp!=1]
		}
		else
		{
			if(x.exogenous==T)
			{
				x.m <- cbind(1,x.m)
				x.m.names <- c("INTERCEPT_mean",paste(xa.names,"mean",sep="_"))
			}
			if(x.exogenous==F)	
			{
				z.m <- cbind(1,z.m)
				z.m.names <- c("INTERCEPT_mean",paste(za.names,"mean",sep="_"))
			}
		}

		if(x.exogenous==F)
		{
			x <- cbind(x,z.m)
			z <- cbind(z,z.m)
			x.names <- c(x.names,z.m.names)

			if(type=="QMLcre")
			{
				results <- lm(var.endog ~ 0+z)
				PIhat <- results$coefficients
				vhat <- results$residuals

				gixv <- t(z*vhat)

				x <- cbind(x,vhat)
				x1.names <- c(x.names,"vhat")
				x2.names <- paste("Z",c(z.names,z.m.names),sep="_")
				x.names <- c(x1.names,x2.names)
			}
		}
		else
		{
			x <- cbind(x,x.m)
			z <- x
			x.names <- c(x.names,x.m.names)

			gixv <- NA
			vhat <- NA
		}

		z.full <- z
	}
	else
	{
		gixv <- NA
		vhat <- NA
	}

	N <- length(unique(id))
	Ti <- as.vector(by(id,id,length))
	NT <- sum(at)

	k <- ncol(x)
	kz <- ncol(z)
	dfJ <- kz-k

	### 4. Estimation

	if(any(type==c("GMMc","GMMww","GMMbgw","GMMpre","GMMpfe","GMMcre"))) Hy <- frmpd.links(link)$H1(y)
	if(type=="QMLcre") Hy <- y

	if(type=="GMMww" & GMMww.cor==T)
	{
		x <- x-matrix(apply(x,2,mean),nrow=nrow(x),ncol=k,byrow=T)
		z <- z-matrix(apply(z.full,2,mean),nrow=nrow(z),ncol=kz,byrow=T)
	}

	if(missing(start)) start <- rep(0,k)
	if(length(start)!=k) stop("start is not of the same dimension as the covariate vector (includes all auxiliary parameters)")

	results <- frmpd.est(type,x.exogenous,lags,id,Ti,Hy,x,z,link,var.type,start,at,at1,variance,NT,k,kz,gixv,vhat,bootstrap,...)
	p <- results$p
	converged <- results$converged

	if(type!="QMLcre" & dfJ>0)
	{
		Qn <- results$Qn
		J <- NT*Qn
	}

	if(variance==T & converged==T)
	{
 		if(bootstrap==F) p.var <- results$p.var
		else
		{
			Ti.full <- as.vector(by(id.full,id.full,length))

			pboot <- matrix(NA,nrow=B,ncol=length(x.names))

			for(j in 1:B)
			{
				index.id <- sample(unique(id),N,replace=T)
				a <- 1
				aa <- 1

				for(jj in unique(id))
				{
					n.id <- sum(index.id==jj)
					if(n.id>=1)
					{
						bb.1 <- (a!=1)*sum(Ti.full[1:(a-1)])+1
						bb.2 <- sum(Ti.full[1:a])

						refa <- rep(bb.1:bb.2,n.id)
						refe <- rep(aa:(aa+n.id-1),each=Ti.full[a])

						if(aa==1)
						{
							ref <- refa
							id.B <- refe
						}
						else
						{
							ref <- c(ref,refa)
							id.B <- c(id.B,refe)
						}

						aa <- aa+n.id
					}

					a <- a+1
				}

				Hy.B <- Hy[ref]
				x.B <- cbind(x[ref,])
				z.B <- cbind(z.full[ref,])
				if(type=="QMLcre" & x.exogenous==F) var.endog.B <- var.endog[ref]
				at.B <- at[ref]
				at1.B <- at1[ref]

				if(any(type==c("GMMc","GMMww")) | lags==T)
				{
					if(lags==T)
					{
						z.B <- cbind(z.B[at1.B,])
						id.B <- id.B[at1.B]
					}
					if(lags==F)
					{
						z.B <- cbind(z.B[at.B,])
						id.B <- id.B[at.B]
					}
				}

				Ti.B <- as.vector(by(id.B,id.B,length))
				NT.B <- sum(at.B)

				if(type=="QMLcre" & x.exogenous==F)
				{
					results <- lm(var.endog.B ~ 0+z.B)
					PIres <- results$coefficients
				}

				results <- frmpd.est(type,x.exogenous,lags,id.B,Ti.B,Hy.B,x.B,z.B,link,var.type,start,at.B,at1.B,variance,NT.B,k,kz,NA,NA,bootstrap,...)
				if(results$converged==T)
				{
					cat("1")
					if(any(j==seq(50,100000,50))) cat("\n")

					pres <- results$p
					if(type=="QMLcre" & x.exogenous==F) pres <- c(pres,PIres)
					pboot[j,] <- pres
				}
				else cat("0")
			}

			p.var <- matrix(NA,nrow=length(x.names),ncol=length(x.names))
			diag(p.var) <- apply(pboot,2,var,na.rm=T)
		}
	}
	else p.var <- "singular"

	if(type=="QMLcre" & x.exogenous==F) p <- c(p,PIhat)

	if(table==T) frmpd.table(p,p.var,x.names,x.exogenous,lags,type,link,converged,N.ini,N,NT.ini,NT,J,dfJ,k,var.type,bootstrap)

	names(p) <- x.names
	res <- list(type=type,link=link,Hy=Hy,p=p,converged=converged)
	if(dfJ>0) res[["J"]] <- J

	if(variance==T & converged==T)
	{ 
		if(is.character(p.var)) p.var <- matrix(NA,nrow=length(p),ncol=length(p))
		dimnames(p.var) <- list(x.names,x.names)
		res[["p.var"]] <- p.var
		res[["var.type"]] <- var.type
	}

	### 5. Return results

	return(invisible(res))
}


