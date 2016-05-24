least.rect <-
function(formula,data,conf.level=0.95,theo=1,adj=TRUE){
  if (missing(formula)) {stop("missing formula")}
  if (length(formula[[3]])==3) {
    if (formula[[3]][[1]]!=as.name("|")) {
	stop("incorrect formula")
    } else {
	formula[[3]][[1]] <- as.name("+")
    }
  }
  m <- match.call()
  m$formula <- formula
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$conf.level <- m$theo <- m$adj <- NULL
  mf <- na.omit(eval(m,parent.frame()))
  dname <- names(mf)
  lr <- function(x,y,conf.level=conf.level,theo=theo) {
    corr <- cor.test(x,y,method="pearson",conf.level=conf.level)
    r <- as.numeric(corr$estimate)
    k <- qt((1+conf.level)/2,length(x)-2)^2*(1-r^2)/(length(x)-2)
    b <- sd(y)/sd(x)*sign(cov(x,y))
    b.ci1 <- b*sqrt(1+2*k-sqrt((1+2*k)^2-1))
    b.ci2 <- b*sqrt(1+2*k+sqrt((1+2*k)^2-1))
    b.inf <- min(b.ci1,b.ci2)
    b.sup <- max(b.ci1,b.ci2)
    a <- mean(y)-b*mean(x)
    a.inf <- mean(y)-b.sup*mean(x)
    a.sup <- mean(y)-b.inf*mean(x)
    t.obs <- abs(b^2-theo^2)*sqrt(length(x)-2)/(2*b*theo*sqrt(1-r^2))
    p <- min(pt(t.obs,length(x)-2),pt(t.obs,length(x)-2,lower.tail=FALSE))*2
    conf.int <- matrix(c(a.inf,b.inf,a,b,a.sup,b.sup),nrow=2,dimnames=list(c("(Intercept)","x"),
	c("inf","coeff","sup")))
    conform <- data.frame("observed"=b,"theoretical"=theo,"Df"=length(x)-2,"t"=t.obs,"Pr(>|t|)"=p,
	row.names=" ",check.names=FALSE)
    p.corr <- as.numeric(corr$p.value)
    corr.tab <- data.frame("inf"=as.numeric(corr$conf.int[1]),"r"=r,"sup"=as.numeric(corr$conf.int[2]),
	"Df"=as.numeric(corr$parameter),"t"=as.numeric(corr$statistic),"Pr(>|t|)"=p.corr,row.names=" ",
	check.names=FALSE)
    coeffs <- c(a,b)
    names(coeffs) <- c("(Intercept)","x")
    fit <- a+b*x
    names(fit) <- 1:length(x)
    res <- -(b*x-y+a)/sqrt(1+b^2)
    names(res) <- 1:length(x)
    result <- list(coef=coeffs,res=res,fit=fit,conf.int=conf.int,comp=conform,corr=corr.tab)
    return(result)
  }
  if (ncol(mf)>2) {
    if (ncol(mf)==3) {
	if (length(formula[[3]][[3]])!=1) {stop("incorrect formula")}
	if (!is.factor(mf[,3])) {stop("incorrect factor")}
	mf <- droplevels(mf)
	if (nlevels(mf[,3])==1) {mf <- mf[,-3]}
    } else {stop("incorrect formula")}
  }
  if (ncol(mf)==2) {
    mod <- lr(mf[,2],mf[,1],conf.level=conf.level,theo=theo)
    coef <- mod$coef
    res <- mod$res
    fit <- mod$fit
    conf.int <- mod$conf.int
    comp <- mod$comp
    corr <- mod$corr
    multiple <- FALSE
  } else {
    lev <- levels(mf[,3])
    nlev <- length(lev)
    coef <- data.frame("(Intercept)"=integer(nlev),x=integer(nlev),row.names=lev,check.names=FALSE)
    colnames(coef)[2] <- dname[2]
    res <- NULL
    fit <- NULL
    conf.int <- array(0,c(2,3,2),dimnames=list(lev,c("inf","coeff","sup"),c("(Intercept)",dname[2])))
    comp <- data.frame(observed=integer(nlev),theoretical=rep(theo,nlev),Df=integer(nlev),t=integer(nlev),
	"Pr(>|t|)"=integer(nlev),row.names=lev,check.names=FALSE)
    corr <- data.frame(inf=integer(nlev),r=rep(theo,nlev),sup=integer(nlev),Df=integer(nlev),t=integer(nlev),
	"Pr(>|t|)"=integer(nlev),row.names=lev,check.names=FALSE)
    multiple <- TRUE
    for (i in 1:nlev) {
	tab <- subset(mf,as.numeric(mf[,3])==i)
	cl <- if (adj) {
	  1-((1-conf.level)/nlev)
	} else {
	  conf.level
	}
	mod <- lr(tab[,2],tab[,1],conf.level=cl,theo=theo)
	coef[i,] <- mod$coef
	res.temp <- mod$res
	names(res.temp) <- rep(lev[i],length(res.temp))
	res <- c(res,res.temp)
	fit.temp <- mod$fit
	names(fit.temp) <- rep(lev[i],length(fit.temp))
	fit <- c(fit,fit.temp)
	conf.int[i,,1] <- mod$conf.int[1,]
	conf.int[i,,2] <- mod$conf.int[2,]
	comp[i,] <- mod$comp[1,]
	corr[i,] <- mod$corr[1,]
    }
    m <- match.call()
    if (adj) {
	comp[,"Pr(>|t|)"] <- p.adjust(comp[,"Pr(>|t|)"],method="bonferroni")
	corr[,"Pr(>|t|)"] <- p.adjust(corr[,"Pr(>|t|)"],method="bonferroni")
    }
  }
  m[[1]] <- as.name("least.rect")
  result <- list(coefficients=coef,residuals=res,fitted.values=fit,call=m,model=mf,conf.level=conf.level,
    conf.int=conf.int,theo=theo,comp=comp,corr=corr,multiple=multiple,adj=adj)
  class(result) <- "least.rect"
  return(result)
}

print.least.rect <- function(x,...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\n")
}

summary.least.rect <- function(object,...) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nResiduals:\n")
  print(summary(object$residuals),digits=3)
  if (object$multiple && object$adj) {
    cat("\nCoefficients and Bonferroni-adjusted ",100*object$conf.level,"% confidence interval:\n",sep="")
  } else {
    cat("\nCoefficients and ",100*object$conf.level,"% confidence interval:\n",sep="")
  }
  if (!object$multiple) {
    print(object$conf.int,digits=5)
  } else {
    for (i in dimnames(object$conf.int)[[3]]) {
	cat(" ",i,"\n")
	print(object$conf.int[,,i])
    }
  }
  if (object$multiple && object$adj) {
    cat("\nEquality of the slope",ifelse(object$multiple,"s","")," to ",object$theo,
	" (Bonferroni-adjusted p-value):\n",sep="")
  } else {
    cat("\nEquality of the slope",ifelse(object$multiple,"s","")," to ",object$theo,":\n",sep="")
  }
  printCoefmat(object$comp,na.print="",P.values=TRUE,has.Pvalue=TRUE)
  if (object$multiple && object$adj) {
    cat("\nPearson's linear correlation coefficient",ifelse(object$multiple,"s",""),
	" (Bonferroni-adjusted ",100*object$conf.level,"% confidence interval and p-value):\n",sep="")
  } else {
    cat("\nPearson's linear correlation coefficient",ifelse(object$multiple,"s",""),
	" (",100*object$conf.level,"% confidence interval):\n",sep="")
  }
  printCoefmat(object$corr,na.print="",P.values=TRUE,has.Pvalue=TRUE)
  cat("\n")
}


