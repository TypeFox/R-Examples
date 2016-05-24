frm.links <- function(link) 
{
	switch(link,
		logit = {
			linkfun <- function(mu) qlogis(mu)
			linkinv <- function(eta) plogis(eta)
			mu.eta <- function(eta) dlogis(eta)
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
			gd <- function(eta) (exp(-exp(eta))*exp(eta))*(1-exp(eta))
			valideta <- function(eta) TRUE
		},

		loglog = {
			linkfun <- function(mu) -log(-log(mu))
			linkinv <- function(eta) exp(-exp(-eta))
			mu.eta <- function(eta) exp(-exp(-eta)-eta)
			gd <- function(eta) (exp(-exp(-eta))*exp(-eta))*(exp(-eta)-1)
			valideta <- function(eta) TRUE
		},

		stop(sQuote(link), " - link not recognised")
	)

	structure(list(linkfun=linkfun,linkinv=linkinv,mu.eta=mu.eta,gd=gd,valideta=valideta,name=link),class="link-glm")
}

frm.est <- function(y,x,link,method,variance,var.type,var.eim,var.cluster,dfc,...)
{
	if(method=="ML") results <- glm(y ~ x-1,family=binomial(link=frm.links(link)),...)
	if(method=="QML") results <- glm(y ~ x-1,family=quasibinomial(link=frm.links(link)),...)
	p <- results$coefficients
	xbhat <- results$linear.predictors
	yhat <- results$fitted.values
	converged <- results$converged*(1-results$boundary)
	if(method=="ML") LL <- logLik(results)

	ret.list <- list(p=p,yhat=yhat,xbhat=xbhat,converged=converged)
	if(method=="ML") ret.list[["LL"]] <- LL

	if(variance==F) return(ret.list)

	p.var <- frm.var(y,x,yhat,xbhat,link,var.type,var.eim,var.cluster,dfc)$var
	ret.list[["p.var"]] <- p.var

	return(ret.list)
}

frm.var <- function(y,x,yhat,xbhat,link,var.type,var.eim,var.cluster,dfc)
{
	n <- nrow(x)
	uhat <- y-yhat
	g <- frm.links(link)$mu.eta(xbhat)
	gd <- frm.links(link)$gd(xbhat)

	A <- 0
	B <- 0

	for(jj in 1:n)
	{
		xx <- x[jj,]%*%t(x[jj,])
		A1 <- xx*(g[jj]^2)/(yhat[jj]*(1-yhat[jj]))
		if(var.eim==T) A <- A+A1
		if(var.eim==F)
		{
			div <- (yhat[jj]*(1-yhat[jj]))
			A2 <- -xx*uhat[jj]*gd[jj]/div
			A3 <- xx*uhat[jj]*(g[jj]^2)/(div*yhat[jj])
			A4 <- -xx*uhat[jj]*(g[jj]^2)/(div*(1-yhat[jj]))

			A <- A+A1+A2+A3+A4
		}
		if(var.type=="robust") B <- B+xx*(uhat[jj]^2)*(g[jj]^2)/((yhat[jj]*(1-yhat[jj]))^2)
	}

	if(var.type=="cluster")
	{
		id <- var.cluster
		id.uni <- unique(id)

		for(j in id.uni)
		{
			Xi <- matrix(x[id==j,],ncol=ncol(x))
			yhati <- yhat[id==j]
			gi <- g[id==j]
			ui <- uhat[id==j]

			ugGi <- ui*gi/(yhati*(1-yhati))

			B <- B+t(Xi)%*%ugGi%*%t(ugGi)%*%Xi
		}
	}

	if(dfc==T)
	{
		if(any(var.type==c("standard","robust"))) df <- n/(n-1)
		if(var.type=="cluster") df <- length(id.uni)/(length(id.uni)-1)
	}
	else df <- 1

	A.inv <- solve(A)

	if(var.type=="standard") var <- df*A.inv
	if(var.type!="standard") var <- df*A.inv%*%B%*%A.inv


	return(list(var=var))
}

frm.table <- function(y,yhat,p,p.var,x.names,type,link,converged,var.type)
{
	if(converged==T)
	{
		R2 <- cor(y,yhat)^2

		if(type!="2P")
		{
			n <- length(y)

			p.sd <- diag(p.var)^0.5
			z.ratio <- p/p.sd
			p.value <- 2*(1-pnorm(abs(z.ratio)))

			stars <- rep("",length(p))
			stars[p.value<=0.01] <- "***"
			stars[p.value>0.01 & p.value<=0.05] <- "**"
			stars[p.value>0.05 & p.value<=0.1] <- "*"

			p <- formatC(p,digits=6,format="f")
			p.sd <- formatC(p.sd,digits=6,format="f")
			z.ratio <- formatC(z.ratio,digits=3,format="f")
			p.value <- formatC(p.value,digits=3,format="f")
			stars <- format(stars,justify="left")

			results <- data.frame(cbind(p,p.sd,z.ratio,p.value,stars),row.names=NULL)

			namcol <- c("Estimate","Std. Error","t value","Pr(>|t|)","")
			dimnames(results) <- list(x.names,namcol)

			cat("\n")
			if(type=="2Pbin") cat("*** Binary component of a two-part model -",link,"specification ***")
			if(type=="1P") cat("*** Fractional",link,"regression model ***")
			if(type=="2Pfrac") cat("*** Fractional component of a two-part model -",link,"specification ***")
			cat("\n\n")

			print(results)
			cat("\n")
			if(var.type!="standard")
			{
				cat("Note:",var.type,"standard errors")
				cat("\n\n")
			}
			cat("Number of observations:",n,"\n")
			cat("R-squared:",round(R2,3),"\n")
			cat("\n")
		}
		else
		{
			cat("\n")
	     		cat("*** Two-part model - binary",link[1],"+ fractional",link[2]," ***")
			cat("\n\n")
			cat("R-squared:",round(R2,3),"\n")
			cat("\n")
		}
	}
	else cat("ALGORITHM DID NOT CONVERGE OR STOPPED AT A BOUNDARY VALUE")
	cat("\n")
}

frm.lm <- function(y,yhat,gx,gz,type)
{
	u <- y-yhat
	w <- as.vector(sqrt(yhat*(1-yhat)))

	uw <- u/w
	gxw <- gx/w

	gzw <- gz/w

	if(type=="2Pbin")
	{
		gxzw <- cbind(gxw,gzw)

		res <- lm(uw ~ gxzw-1)

		LM <- t(res$fitted.values)%*%uw

	}
	else
	{
		res <- lm(gzw ~ gxw-1)
		ur <- as.matrix(res$residuals)

		n <- length(u)
		ones <- rep(1,n)

		uwr <- as.matrix(uw*ur)
		res <- lm(ones ~ uwr-1)
		uu <- res$residuals
		RSS <- sum(uu^2)

		LM <- n-RSS
	}

	return(list(LM=LM))
}

frm.tests.table <- function(test.which,S,Sp,ver,title1,title2=NA,n.ver=NA,test.ggoff=NA)
{
	stars <- rep("",length(S))
	stars[Sp<=0.01] <- "***"
	stars[Sp>0.01 & Sp<=0.05] <- "**"
	stars[Sp>0.05 & Sp<=0.1] <- "*"

	if(all(is.na(S))) stop("NO WALD/LR TEST HAS ACHIEVED CONVERGENCE. USE OTHER TEST VERSIONS INSTEAD")

	S <- formatC(S,digits=3,format="f")
	Sp <- formatC(Sp,digits=3,format="f")
	stars <- format(stars,justify="left")

	if(test.which!="GGOFF")
	{
		results <- data.frame(cbind(ver,S,Sp,stars))
		namcol <- c("Version","Statistic","p-value","")	}
	else
	{
		results <- data.frame(cbind(test.ggoff,ver,S,Sp,stars))
		namcol <- c("Test","Version","Statistic","p-value","")
	}

	results <- results[-1,]
	dimnames(results) <- list(1:nrow(results),namcol)

	cat("\n")
	cat("***",test.which,"test ***")
	cat("\n\n")
	cat("H0: ",title1)
	cat("\n")
	if(test.which!="P")
	{
		cat("\n")
		print(results,row.names=F)
	}
	else
	{
		cat("H1: ",title2)
		cat("\n\n")
		print(results[1:n.ver,],row.names=F)
		cat("\n")
		cat("H0: ",title2)
		cat("\n")
		cat("H1: ",title1)
		cat("\n\n")
		print(results[(n.ver+1):(2*n.ver),],row.names=F)
	}
	cat("\n")
}

frm.pe.var <- function(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb=NA,xbhatb=NA,gb=NA,linkb=NA,pb.var=NA,yhata=NA,yhatb=NA)
{
	gda <- frm.links(linka)$gd(xbhata)
	PE1.sd <- matrix(NA,nrow=npar,ncol=npar)

	if(type!="2P")
	{
		for(i in 1:npar)
		{
			for(j in 1:npar) PE1.sd[i,j] <- mean(pa[i]*gda*x[,j]+(i==j)*ga)
		}

		PE.sd <- PE1.sd%*%pa.var%*%t(PE1.sd)
	}
	if(type=="2P")
	{
		gdb <- frm.links(linkb)$gd(xbhatb)
		PE2.sd <- matrix(NA,nrow=npar,ncol=npar)

		for(i in 1:npar)
		{
			for(j in 1:npar)
			{
				PE1.sd[i,j] <- mean(pb[i]*gb*ga*x[,j]+(i==j)*ga*yhatb+pa[i]*gda*yhatb*x[,j])
				PE2.sd[i,j] <- mean(pa[i]*ga*gb*x[,j]+(i==j)*gb*yhata+pb[i]*gdb*yhata*x[,j])
			}
		}

		PE.sd <- PE1.sd%*%pa.var%*%t(PE1.sd)+PE2.sd%*%pb.var%*%t(PE2.sd)
	}

	PE.sd <- diag(PE.sd)^0.5
	if(any(x.names=="INTERCEPT")) PE.sd <- PE.sd[-1]

	names(PE.sd) <- xvar.names
	PE.sd <- PE.sd[which.x]

	return(PE.sd)
}

frm.pe.table <- function(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
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
	if(PE.type=="APE") cat("*** Average partial effects ***")
	if(PE.type=="CPE") cat("*** Conditional partial effects ***")
	cat("\n\n")
	cat(title)
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

frm <- function(y,x,x2=x,linkbin,linkfrac,type="1P",inflation=0,intercept=T,table=T,variance=T,var.type="default",var.eim=T,var.cluster,dfc=F,...)
{
	### 1. Error and warning messages

	if(missing(y)) stop("dependent variable is missing")
	if(missing(x)) stop("explanatory variables are missing")

	if(all(type!=c("1P","2Pbin","2Pfrac","2P"))) stop(sQuote(type)," - type not recognised")
	if(any(y>1) | any(y<0)) stop("The dependent variable has values outside the unit interval")
	if(all(inflation!=c(0,1))) stop(inflation," - value not recognised for inflation")
	if(length(inflation)>1) stop(inflation," - only a single value allowed for inflation")
	if(type!="1P" & !any(y==inflation)) stop("The dependent variable has no ",sQuote(inflation)," values")

	if(type=="2Pbin")
	{
		if(missing(linkbin)) stop("linkbin is missing")
		if(!missing(linkfrac)) warning(sQuote(type)," and",sQuote(linkfrac)," - type does not use linkfrac")

		if(all(linkbin!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkbin)," - linkbin not recognised")
	}
	if(any(type==c("1P","2Pfrac")))
	{
		if(missing(linkfrac)) stop("linkfrac is missing")
		if(!missing(linkbin)) warning(sQuote(type)," and",sQuote(linkbin)," - type does not use linkbin")

		if(all(linkfrac!=c("logit","probit","cauchit","cloglog","loglog"))) stop(sQuote(linkfrac)," - linkfrac not recognised")
	}

	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
	}

	if(all(var.type!=c("standard","robust","cluster","default"))) stop(sQuote(var.type)," - var.type not recognised")
	if(var.type=="cluster" & missing(var.cluster)) stop("option cluster for covariance matrix but no var.cluster supplied")

	if(!is.logical(intercept)) stop("non-logical value assigned to option intercept")
	if(!is.logical(table)) stop("non-logical value assigned to option table")
	if(!is.logical(variance)) stop("non-logical value assigned to option variance")
	if(!is.logical(var.eim)) stop("non-logical value assigned to option var.eim")
	if(!is.logical(dfc)) stop("non-logical value assigned to option dfc")

	### 2. Data and variables preparation

	if(any(type==c("1P","2Pfrac"))) x2 <- x

	if(is.data.frame(x)) x <- as.matrix(x)
	if(is.data.frame(x2)) x2 <- as.matrix(x2)

	if(!is.matrix(x)) stop("x is not a matrix")
	if(!is.matrix(x2)) stop("x2 is not a matrix")

	x.names <- dimnames(x)[[2]]
	x2.names <- dimnames(x2)[[2]]

	if(is.null(x.names)) stop("x has no column names")
	if(is.null(x2.names)) stop("x2 has no column names")

	if(intercept==T)
	{
		x <- cbind(1,x)
		x2 <- cbind(1,x2)

		x.names <- c("INTERCEPT",x.names)
		x2.names <- c("INTERCEPT",x2.names)
	}

	if(length(x.names)!=length(unique(x.names))) stop("some covariate names in x are identical")
	if(length(x2.names)!=length(unique(x2.names))) stop("some covariate names in x2 are identical")

	if(length(y)!=nrow(x)) stop("the number of observations for y and x are different")
	if(var.type=="cluster")
	{
		if(length(y)!=length(var.cluster)) stop("var.cluster does not have the appropriate dimension")
	}

	if(any(type==c("2Pbin","2P")))
	{
		if(inflation==0) yb <- y>0
		if(inflation==1) yb <- y==1
	}

	if(any(type==c("2Pfrac","2P")))
	{
		if(inflation==0)
		{
			yf <- y[y>0]
			x2f <- x2[y>0,]
			if(var.type=="cluster") var.cluf <- var.cluster[y>0]
		}
		if(inflation==1)
		{
			yf <- y[y<1]
			x2f <- x2[y<1,]
			if(var.type=="cluster") var.cluf <- var.cluster[y<1]
		}

		if(length(yf)!=nrow(x2f)) stop("the number of observations for y and x2 are different")
	}

	### 3. Estimation

	class <- "frm"

	if(any(type==c("2Pbin","2P")))
	{
		if(var.type=="default") var.ty <- "standard"
		else var.ty <- var.type

		method <- "ML"
		results <- frm.est(yb,x,linkbin,method,variance,var.ty,var.eim,var.cluster,dfc,...)
		p <- results$p
		if(variance==T) p.var <- results$p.var
		yhat1 <- results$yhat
		xbhat <- results$xbhat
		converged1 <- results$converged
		LL <- results$LL

		if(table==T) frm.table(yb,yhat1,p,p.var,x.names,"2Pbin",linkbin,converged1,var.ty)

		formula <- yb ~ x - 1
		names(p) <- x.names

		resBIN <- list(class=class,formula=formula,type=type,link=linkbin,method=method,p=p,yhat=yhat1,xbhat=xbhat,converged=converged1,LL=LL,x.names=x.names)
		if(variance==T)
		{ 
			dimnames(p.var) <- list(x.names,x.names)
			resBIN[["p.var"]] <- p.var
			resBIN[["var.type"]] <- var.ty
			resBIN[["var.eim"]] <- var.eim
			resBIN[["dfc"]] <- dfc
			if(var.type=="cluster") resBIN[["var.cluster"]] <- var.cluster
		}
	}

	if(any(type==c("1P","2Pfrac","2P")))
	{
		if(type=="1P")
		{
			yy <- y
			xx2 <- x2
			ty <- "1P"
			if(var.type=="cluster") var.clu <- var.cluster
		}
		else
		{
			yy <- yf
			xx2 <- x2f
			ty <- "2Pfrac"
			if(var.type=="cluster") var.clu <- var.cluf
		}

		if(var.type=="default") var.ty <- "robust"
		else var.ty <- var.type

		method <- "QML"
		results <- frm.est(yy,xx2,linkfrac,method,variance,var.ty,var.eim,var.clu,dfc,...)
		p <- results$p
		if(variance==T) p.var <- results$p.var
		yhat <- results$yhat
		xbhat <- results$xbhat
		converged2 <- results$converged

		if(table==T) frm.table(yy,yhat,p,p.var,x2.names,ty,linkfrac,converged2,var.ty)

		formula <- yy ~ xx2 - 1
		names(p) <- x2.names

		resFRAC <- list(class=class,formula=formula,type=type,link=linkfrac,method=method,p=p,yhat=yhat,xbhat=xbhat,converged=converged2,x.names=x2.names)
		if(variance==T)
		{ 
			dimnames(p.var) <- list(x2.names,x2.names)
			resFRAC[["p.var"]] <- p.var
			resFRAC[["var.type"]] <- var.ty
			resFRAC[["var.eim"]] <- var.eim
			resFRAC[["dfc"]] <- dfc
			if(var.type=="cluster") resFRAC[["var.cluster"]] <- var.clu
		}
	}

	if(type=="2P")
	{
		yhat2 <- frm.links(linkfrac)$linkinv(x2%*%p)
		yhat <- yhat1*yhat2

		converged <- converged1*converged2

		if(table==T) frm.table(y,yhat,NA,NA,NA,type,c(linkbin,linkfrac),converged)

		ybase <- y
		x2base <- x2
	}

	### 4. Return results

	if(type=="2Pbin") return(invisible(resBIN))
	if(any(type==c("1P","2Pfrac"))) return(invisible(resFRAC))
	if(type=="2P") return(invisible(list(resBIN=resBIN,resFRAC=resFRAC,class=class,type=type,ybase=ybase,x2base=x2base,yhat2P=yhat,converged=converged)))
}

frm.reset <- function(object,lastpower.vec=3,version="LM",table=T,...)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an frm command")
	if(object$class!="frm") stop("object is not the output of an frm command")
	if(object$type=="2P") stop("The RESET test is not applicable to a two-part model")
	if(object$converged==0) stop("object is not the output of a successful (converged) frm command")
	if(object$method!="ML" & any(version=="LR")) stop("LR tests require ML estimation")
	if(any(lastpower.vec<2)) stop(sQuote(lastpower.vec)," - lastpower.vec contains elements lower than 2")
	if(all(version!="LM") & all(version!="Wald") & all(version!="LR")) stop("test version not correctly specified")
	if(!is.logical(table)) stop("non-logical value assigned to option table")

	### 2. Recovering definitions and estimates

	mf <- model.frame(object$formula)
	y <- model.response(mf)
	x <- model.matrix(object$formula)

	yhat <- object$yhat
	xbhat <- object$xbhat
	method <- object$method
	link <- object$link
	type <- object$type

	if(any(version=="Wald"))
	{
		if(is.null(object$var.type)) stop("Wald test required but frm command was run with variance = F")
		var.type <- object$var.type
		var.eim <- object$var.eim
		dfc <- object$dfc
		if(var.type=="cluster") var.cluster <- object$var.cluster
	}
	if(any(version=="LR")) LL0 <- object$LL
	if(type=="1P") title <- paste("Fractional",link,"model")
	if(type=="2Pbin") title <- paste("Binary",link,"component of a two-part model")
	if(type=="2Pfrac") title <- paste("Fractional",link,"component of a two-part model")

	### 3. Tests

	lastpower.vec <- round(lastpower.vec,0)

	g <- frm.links(link)$mu.eta(xbhat)
	gx <- g*x

	z.all <- as.matrix(xbhat^2)
	if(max(lastpower.vec)>2) for(i in 3:max(lastpower.vec)) z.all <- cbind(z.all,xbhat^i)

	ver <- NA
	S <- NA
	Sp <- NA

	for(j in lastpower.vec)
	{
		df <- j-1
		z <- z.all[,1:(j-1)]
		gz <- g*z

		if(any(version=="LM"))
		{
			name <- paste("LM(",j,")",sep="")
			ver <- c(ver,name)

			results <- frm.lm(y,yhat,gx,gz,type)

			Sj <- results$LM
			S <- c(S,Sj)
			Sp <- c(Sp,1-pchisq(Sj,df))
		}
		if(any(version=="LR") | any(version=="Wald"))
		{
			if(any(version=="Wald")) results <- frm.est(y,cbind(x,z),link,method,T,var.type,var.eim,var.cluster,dfc,...)
			else results <- frm.est(y,cbind(x,z),link,method,F,...)

			if(any(version=="LR"))
			{
				name <- paste("LR(",j,")",sep="")
				ver <- c(ver,name)

				if(results$converged==T)
				{
					LL1 <- results$LL

					Sj <- 2*(LL1-LL0)
					S <- c(S,Sj)
					Sp <- c(Sp,1-pchisq(Sj,df))
				}
				else
				{
					S <- c(S,NA)
					Sp <- c(Sp,NA)
				}
			}
			if(any(version=="Wald"))
			{
				name <- paste("Wald(",j,")",sep="")
				ver <- c(ver,name)

				if(results$converged==T)
				{
					p1 <- results$p
					p.var1 <- results$p.var

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
		}
	}

	if(table==T) frm.tests.table("RESET",S,Sp,ver,title)

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- ver[-1]

	return(invisible(statistics))
}

frm.ggoff <- function(object,version="LM",table=T,...)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an frm command")
	if(object$class!="frm") stop("object is not the output of an frm command")
	if(object$type=="2P") stop("GGOFF tests are not applicable to a two-part model")
	if(object$converged==0) stop("object is not the output of a successful (converged) frm command")
	if(object$method!="ML" & any(version=="LR")) stop("LR tests require ML estimation")
	if(all(version!="LM") & all(version!="Wald") & all(version!="LR")) stop("test version not correctly specified")
	if(!is.logical(table)) stop("non-logical value assigned to option table")

	### 2. Recovering definitions and estimates

	mf <- model.frame(object$formula)
	y <- model.response(mf)
	x <- model.matrix(object$formula)

	yhat <- object$yhat
	xbhat <- object$xbhat
	method <- object$method
	link <- object$link
	type <- object$type

	if(any(version=="Wald"))
	{
		if(is.null(object$var.type)) stop("Wald test required but frm command was run with variance = F")
		var.type <- object$var.type
		var.eim <- object$var.eim
		dfc <- object$dfc
		if(var.type=="cluster") var.cluster <- object$var.cluster
	}
	if(any(version=="LR")) LL0 <- object$LL
	if(type=="1P") title <- paste("Fractional",link,"model")
	if(type=="2Pbin") title <- paste("Binary",link,"component of a two-part model")
	if(type=="2Pfrac") title <- paste("Fractional",link,"component of a two-part model")

	### 3. Tests

	ver <- NA
	S <- NA
	Sp <- NA
	test <- NA
	tests <- c("GOFF1","GOFF2","GGOFF")

	g <- frm.links(link)$mu.eta(xbhat)
	gx <- g*x
	if(link!="loglog") z1 <- yhat*log(yhat)/g
	if(link!="cloglog") z2 <- (1-yhat)*log(1-yhat)/g

	for(j in 1:3)
	{
		if((j==1 & link!="loglog") | (j==2 & link!="cloglog") | j==3)
		{
			if(j==1) z <- z1
			if(j==2) z <- z2
			if(j==3)
			{
				if(all(link!=c("loglog","cloglog"))) z <- cbind(z1,z2)
				if(link=="loglog") z <- z2
				if(link=="cloglog") z <- z1
			}
			z <- as.matrix(z)

			df <- ncol(z)

			if(any(version=="LM"))
			{
				test <- c(test,tests[j])
				ver <- c(ver,"LM")

				gz <- g*z
				results <- frm.lm(y,yhat,gx,gz,type)
				Sj <- results$LM
				S <- c(S,Sj)
				Sp <- c(Sp,1-pchisq(Sj,df))
			}
			if(any(version=="LR") | any(version=="Wald"))
			{
				if(any(version=="Wald")) results <- frm.est(y,cbind(x,z),link,method,T,var.type,var.eim,var.cluster,dfc,...)
				else results <- frm.est(y,cbind(x,z),link,method,F,...)

				if(any(version=="LR"))
				{
					test <- c(test,tests[j])
					ver <- c(ver,"LR")

					if(results$converged==T)
					{
						LL1 <- results$LL

						Sj <- 2*(LL1-LL0)
						S <- c(S,Sj)
						Sp <- c(Sp,1-pchisq(Sj,df))
					}
					else
					{
						S <- c(S,NA)
						Sp <- c(Sp,NA)
					}
				}
				if(any(version=="Wald"))
				{
					test <- c(test,tests[j])
					ver <- c(ver,"Wald")

					if(results$converged==T)
					{
						p1 <- results$p
						p.var1 <- results$p.var

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
			}
		}
	}

	if(table==T) frm.tests.table("GGOFF",S,Sp,ver,title,test.ggoff=test)

	### 4. Return results

	statistics <- S[-1]
	names(statistics) <- paste(test[-1],ver[-1],sep="-")

	return(invisible(statistics))
}

frm.ptest <- function(object1,object2,version="Wald",table=T)
{
	### 1. Error and warning messages

	if(missing(object1)) stop("object1 is missing")
	if(missing(object2)) stop("object2 is missing")

	if(is.null(object1$class)) stop("object1 not the output of an frm command")
	if(is.null(object2$class)) stop("object2 not the output of an frm command")

	if(object1$class!="frm") stop("object1 is not the output of an frm command")
	if(object2$class!="frm") stop("object2 is not the output of an frm command")

	if(any(object1$type==c("1P","2P")) & all(object2$type!=c("1P","2P"))) stop("object1 and object2 cannot be compared")
	if(object1$type=="2Pbin" & object2$type!="2Pbin") stop("object1 and object2 cannot be compared")
	if(object1$type=="2Pfrac" & object2$type!="2Pfrac") stop("object1 and object2 cannot be compared")

	if(object1$converged==0) stop("object1 is not the output of a successful (converged) frm command")
	if(object2$converged==0) stop("object2 is not the output of a successful (converged) frm command")

	if(all(version!="LM") & all(version!="Wald")) stop("test version not correctly specified")

	if(!is.logical(table)) stop("non-logical value assigned to option table")

	### 2. Variables and other information for the tests

	# 2.1. Model 1

	if(any(object1$type==c("1P","2Pbin","2Pfrac")))
	{
		link1 <- object1$link
		type1 <- object1$type

		mf <- model.frame(object1$formula)
		y1 <- as.vector(model.response(mf))
		x1 <- model.matrix(object1$formula)
		dimnames(x1)[[2]] <- object1$x.names

		yhat1 <- object1$yhat
		xbhat1 <- object1$xbhat

		g1 <- frm.links(link1)$mu.eta(xbhat1)

		gx1 <- g1*x1

		if(type1=="1P") title1 <- paste("Fractional",link1,"model")
		if(type1=="2Pbin") title1 <- paste("Binary",link1,"component of a two-part model")
		if(type1=="2Pfrac") title1 <- paste("Fractional",link1,"component of a two-part model")

		if(type1=="2Pbin") type.both <- "2Pbin"
		else type.both <- "others"
	}
	if(object1$type=="2P")
	{
		y1 <- object1$ybase

		link1a <- object1$resBIN$link
		type1a <- "2Pbin"

		x1a <- model.matrix(object1$resBIN$formula)
		yhat1a <- object1$resBIN$yhat
		xbhat1a <- object1$resBIN$xbhat
		g1a <- frm.links(link1a)$mu.eta(xbhat1a)
		dimnames(x1a)[[2]] <- object1$resBIN$x.names

		link1b <- object1$resFRAC$link
		type1b <- "2Pfrac"

		x1b <- object1$x2base
		xbhat1b <- x1b%*%object1$resFRAC$p
		yhat1b <- frm.links(link1b)$linkinv(xbhat1b)
		g1b <- frm.links(link1b)$mu.eta(xbhat1b)
		dimnames(x1b)[[2]] <- object1$resFRAC$x.names

		yhat1 <- object1$yhat2P
		g1ab <- as.vector(g1a*yhat1b)
		g1ba <- as.vector(g1b*yhat1a)

		gx1 <- cbind(g1ab*x1a,g1ba*x1b)

		title1 <- paste("Binary",link1a,"+ Fractional",link1b,"two-part model")

		type.both <- "others"
	}

	# 2.2. Model 2

	if(any(object2$type==c("1P","2Pbin","2Pfrac")))
	{
		link2 <- object2$link
		type2 <- object2$type

		mf <- model.frame(object2$formula)
		y2 <- as.vector(model.response(mf))
		x2 <- model.matrix(object2$formula)
		dimnames(x2)[[2]] <- object2$x.names

		yhat2 <- object2$yhat
		xbhat2 <- object2$xbhat
		g2 <- frm.links(link2)$mu.eta(xbhat2)

		gx2 <- g2*x2

		if(type2=="1P") title2 <- paste("Fractional",link2,"model")
		if(type2=="2Pbin") title2 <- paste("Binary",link2,"component of a two-part model")
		if(type2=="2Pfrac") title2 <- paste("Fractional",link2,"component of a two-part model")
	}
	if(object2$type=="2P")
	{
		y2 <- object2$ybase

		link2a <- object2$resBIN$link
		type2a <- "2Pbin"

		x2a <- model.matrix(object2$resBIN$formula)
		yhat2a <- object2$resBIN$yhat
		xbhat2a <- object2$resBIN$xbhat
		g2a <- frm.links(link2a)$mu.eta(xbhat2a)
		dimnames(x2a)[[2]] <- object2$resBIN$x.names

		link2b <- object2$resFRAC$link
		type2b <- "2Pfrac"

		x2b <- object2$x2base
		xbhat2b <- x2b%*%object2$resFRAC$p
		yhat2b <- frm.links(link2b)$linkinv(xbhat2b)
		g2b <- frm.links(link2b)$mu.eta(xbhat2b)
		dimnames(x2b)[[2]] <- object2$resFRAC$x.names

		yhat2 <- object2$yhat2P
		g2ab <- as.vector(g2a*yhat2b)
		g2ba <- as.vector(g2b*yhat2a)

		gx2 <- cbind(g2ab*x2a,g2ba*x2b)

		title2 <- paste("Binary",link2a,"+ Fractional",link2b,"two-part model")
	}

	### 3. Further error and warning messages

	if(!all.equal(y1,y2)) stop("The dependent variable is not the same in the two models")

	if(any(object1$type==c("2Pbin","2Pfrac")) | (object1$type=="1P" & object2$type=="1P"))
	{
		if(object1$link==object2$link)
		{
			x1.names <- dimnames(x1)[[2]]
			x2.names <- dimnames(x2)[[2]]
			x12.names <- c(x1.names,x2.names)
			x12.names <- unique(x12.names)
			x1.len <- length(x1.names)
			x2.len <- length(x2.names)
			x12.len <- length(x12.names)

			if(identical(x1.names,x2.names)) stop("object 1 and object 2 are based on the same link function and covariates")
			if(x1.len==x12.len & x2.len!=x12.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1.len!=x12.len & x2.len==x12.len) stop("object 1 is nested in object 2 - no need to use the P test")
		}
	}
	if(object1$type=="2P" & object2$type=="2P")
	{
		if(object1$resBIN$link==object2$resBIN$link & object1$resFRAC$link==object2$resFRAC$link)
		{
			x1a.names <- dimnames(x1a)[[2]]
			x2a.names <- dimnames(x2a)[[2]]
			x12a.names <- c(x1a.names,x2a.names)
			x12a.names <- unique(x12a.names)
			x1a.len <- length(x1a.names)
			x2a.len <- length(x2a.names)
			x12a.len <- length(x12a.names)

			x1b.names <- dimnames(x1b)[[2]]
			x2b.names <- dimnames(x2b)[[2]]
			x12b.names <- c(x1b.names,x2b.names)
			x12b.names <- unique(x12b.names)
			x1b.len <- length(x1b.names)
			x2b.len <- length(x2b.names)
			x12b.len <- length(x12b.names)

			if(identical(x1a.names,x2a.names) & identical(x1b.names,x2b.names)) stop("object 1 and object 2 are based on the same link function and covariates")
			if(x1a.len==x12a.len & x2a.len==x12a.len & x1b.len==x12b.len & x2b.len!=x12b.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1a.len==x12a.len & x2a.len==x12a.len & x1b.len!=x12b.len & x2b.len==x12b.len) stop("object 1 is nested in object 2 - no need to use the P test")
			if(x1a.len==x12a.len & x2a.len!=x12a.len & x1b.len==x12b.len & x2b.len==x12b.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1a.len!=x12a.len & x2a.len==x12a.len & x1b.len==x12b.len & x2b.len==x12b.len) stop("object 1 is nested in object 2 - no need to use the P test")
			if(x1a.len==x12a.len & x2a.len!=x12a.len & x1b.len==x12b.len & x2b.len!=x12b.len) stop("object 2 is nested in object 1 - no need to use the P test")
			if(x1a.len!=x12a.len & x2a.len==x12a.len & x1b.len!=x12b.len & x2b.len==x12b.len) stop("object 1 is nested in object 2 - no need to use the P test")
		}
	}

	### 4. Tests

	ver <- NA
	S <- NA
	Sp <- NA

	df <- 1
	for(j in 1:2)
	{
		if(j==1)
		{
			yj <- y1
			yhatj <- yhat1
			gxj <- gx1
			gzj <- yhat2-yhat1
		}
		if(j==2)
		{
			yj <- y2
			yhatj <- yhat2
			gxj <- gx2
			gzj <- yhat1-yhat2
		}

		if(any(version=="LM"))
		{
			ver <- c(ver,"LM")

			results <- frm.lm(yj,yhatj,gxj,gzj,type.both)
			Sj <- results$LM
			S <- c(S,Sj)
			Sp <- c(Sp,1-pchisq(Sj,df))
		}
		if(any(version=="Wald"))
		{
			ver <- c(ver,"Wald")

			yt <- yj-yhatj
			gxzj <- cbind(gxj,gzj)
			results <- lm(yt ~ gxzj-1)

			gzj.b <- results$coefficients[ncol(gxzj)]
			dfcc <- nrow(gxzj)/(nrow(gxzj)-ncol(gxzj))
			gzj.var <- dfcc*(solve(t(gxzj)%*%gxzj)%*%t(gxzj)%*%diag(results$residuals^2)%*%gxzj%*%solve(t(gxzj)%*%gxzj))[ncol(gxj)+1,ncol(gxj)+1]

			Sj <- as.vector(gzj.b/(gzj.var^0.5))
			S <- c(S,Sj)
			Spj <- 2*(1-pt(abs(Sj),nrow(gxzj)-ncol(gxzj)))
			Sp <- c(Sp,Spj)
		}
	}

	n.ver <- length(ver[-1])/2 
	if(table==T) frm.tests.table("P",S,Sp,ver,title1,title2,n.ver)

	### 5. Return results

	statistics=S[-1]
	ver <- ver[-1]
	names(statistics) <- paste(c(rep("H0-obj1",n.ver),rep("H0-obj2",n.ver)),ver,sep="-")

	return(invisible(statistics))
}

frm.pe <- function(object,APE=T,CPE=F,at=NULL,which.x=NULL,variance=T,table=T)
{
	### 1. Error and warning messages

	if(missing(object)) stop("object is missing")
	if(is.null(object$class)) stop("object is not the output of an frm command")
	if(object$class!="frm") stop("object is not the output of an frm command")

	if(!is.logical(APE)) stop("non-logical value assigned to option APE")
	if(!is.logical(CPE)) stop("non-logical value assigned to option CPE")
	if(!is.logical(variance)) stop("non-logical value assigned to option variance")
	if(!is.logical(table)) stop("non-logical value assigned to option table")

	if(all(c(APE,CPE)==F)) stop("You must specify at least one option: APE and/or CPE")
	if(CPE==F & !is.null(at)) stop("option at is only required for CPE")

	if(object$converged==0) stop("object is not the output of a successful (converged) frm command")

	if(object$type!="2P")
	{
		if(is.null(object$p.var)) stop("frm command was run with variance = F")
	}
	else
	{
		if(is.null(object$resBIN$p.var)) stop("frm command was run with variance = F")
	}

	if(table==T & variance==F)
	{
		variance <- T
		warning("option variance changed from F to T, as required by table=T")
	}

	### 2. Recovering definitions and estimates and other definitions

	type <- object$type

	if(type!="2P")
	{
		linka <- object$link
		pa <- object$p
		pa.var <- object$p.var
		x <- model.matrix(object$formula)
		x.names <- object$x.names

		if(type=="1P") title <- paste("Fractional",linka,"model")
		if(type=="2Pbin") title <- paste("Binary",linka,"component of a two-part model")
		if(type=="2Pfrac") title <- paste("Fractional",linka,"component of a two-part model")
	}
	if(type=="2P")
	{
		linka <- object$resBIN$link
		pa <- object$resBIN$p
		pa.var <- object$resBIN$p.var
		xa <- model.matrix(object$resBIN$formula)
		xa.names <- object$resBIN$x.names

		linkb <- object$resFRAC$link
		pb <- object$resFRAC$p
		pb.var <- object$resFRAC$p.var
		xb <- object$x2base
		xb.names <- object$resFRAC$x.names

		if(!identical(xa.names,xa.names)) stop("currently frm.pe requires both components of two-part models to use the same covariates")
		x <- xa
		x.names <- xa.names

		title <- paste("Binary",linka,"+ Fractional",linkb,"two-part model")
	}

	if(any(x.names=="INTERCEPT")) xvar.names <- x.names[-1]
	else xvar.names <- x.names

	k <- length(xvar.names)
	npar <- ncol(x)
	n <- nrow(x)

	if(is.null(which.x)) which.x <- xvar.names
	xw.names <- unique(c(xvar.names,which.x))
	if(!identical(xvar.names,xw.names)) stop("option which not appropriately defined")

	### 3. Average partial effects

	if(APE==T)
	{
		PE.type <- "APE"

		if(any(x.names=="INTERCEPT")) p.pe <- matrix(rep(pa[-1],each=n),ncol=k)
 		else p.pe <- matrix(rep(pa,each=n),ncol=k)
		dimnames(p.pe) <- list(NULL,xvar.names)

		if(type!="2P") xbhata <- object$xbhat
		if(type=="2P") xbhata <- object$resBIN$xbhat

		ga <- frm.links(linka)$mu.eta(xbhata)
		PE.p <- as.matrix(p.pe[,which.x])*ga

		if(type!="2P") PE.p <- apply(PE.p,2,mean)

		if(type=="2P")
		{
			yhata <- object$resBIN$yhat
			PEa.p <- as.matrix(p.pe[,which.x])*ga

			if(any(x.names=="INTERCEPT")) p.pe <- matrix(rep(pb[-1],each=n),ncol=k)
 			else p.pe <- matrix(rep(pb,each=n),ncol=k)
			dimnames(p.pe) <- list(NULL,xvar.names)

			xbhatb <- as.vector(x%*%pb)
			gb <- frm.links(linkb)$mu.eta(xbhatb)
			yhatb <- frm.links(linkb)$linkinv(xbhatb)
			PEb.p <- as.matrix(p.pe[,which.x])*gb

			PE.p <- apply(PEb.p*yhata+PE.p*yhatb,2,mean)
		}

		resAPE <- list(PE.p=PE.p)

		if(variance==T)
		{
			if(type!="2P") PE.sd <- frm.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var)
			if(type=="2P") PE.sd <- frm.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb,xbhatb,gb,linkb,pb.var,yhata,yhatb)
			resAPE[["PE.sd"]] <- PE.sd
		}

		if(table==T) frm.pe.table(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
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

		xbhata <- as.vector(xm%*%pa)
		ga <- frm.links(linka)$mu.eta(xbhata)
		PE.p <- pa[which.x]*ga

		if(type=="2P")
		{
			yhata <- frm.links(linka)$linkinv(xbhata)

			xbhatb <- as.vector(xm%*%pb)
			gb <- frm.links(linkb)$mu.eta(xbhatb)
			yhatb <- frm.links(linkb)$linkinv(xbhatb)
			PEb.p <- pb[which.x]*gb

			PE.p <- PEb.p*yhata+PE.p*yhatb
		}

		resCPE <- list(PE.p=PE.p)

		if(variance==T)
		{
			if(type!="2P") PE.sd <- frm.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var)
			if(type=="2P") PE.sd <- frm.pe.var(x,npar,which.x,x.names,xvar.names,type,pa,xbhata,ga,linka,pa.var,pb,xbhatb,gb,linkb,pb.var,yhata,yhatb)
			resCPE[["PE.sd"]] <- PE.sd
		}

		if(table==T) frm.pe.table(PE.p,PE.sd,PE.type,which.x,xvar.names,title,at)
	}

	### 5. Return results

	if(APE==T & CPE==T) return(invisible(list(ape=resAPE,cpe=resCPE)))
	if(APE==T & CPE==F) return(invisible(resAPE))
	if(APE==F & CPE==T) return(invisible(resCPE))
}
