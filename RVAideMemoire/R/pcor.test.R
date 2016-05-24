# boot : boot
# pspearman : pspearman

pcor.test <-
function(x,y,z,semi=FALSE,conf.level=0.95,nrep=1000,method=c("pearson","spearman")) {
  method <- match.arg(method)
  if (is.list(z)) {
    if (is.null(names(z))) {
	call.z <-match.call()$z
	if (length(call.z)==(length(z)+1)) {
	  names(z) <- unlist(lapply(call.z,function(a) as.character(a)[length(as.character(a))]))[-1]
	} else {
	  names(z) <- paste0("V",1:length(z))
	}
    }
    z <- as.data.frame(do.call("cbind",z))
  }
  if (is.vector(z)) {
    z <- as.data.frame(z)
    colnames(z) <- as.character(match.call()$z)[length(as.character(match.call()$z))]
  }
  tab <- data.frame(x,y,z)
  tab <- tab[complete.cases(tab),]
  dname <- paste0(deparse(substitute(x))," and ",deparse(substitute(y)),", controlling for ",
    paste(colnames(tab)[-c(1,2)],collapse=", "))
  estimate <- pcor(x,y,z,semi=semi,method=method)
  result <- list(data.name=dname,alternative="two.sided")
  null <- 0
  k <- ncol(tab)-2
  n <- nrow(tab)
  if (method=="pearson") {
    ddl <- n-2-k
    names(ddl) <- "df"
    statistic <- estimate*sqrt(ddl/(1-estimate^2))
    names(statistic) <- "t"
    p <- pt(statistic,ddl)
    pval <- 2*min(p,1-p)
    result$method <- ifelse(!semi,"Pearson's product-moment partial correlation",
	"Pearson's product-moment semi-partial correlation")
    names(estimate) <- "cor"
    names(null) <- "correlation"
    if (n>(3+k)) {
	z.estimate <- atanh(estimate)
	sigma <- 1/sqrt(n-3-k)
	cint <- z.estimate+c(-1,1)*sigma*qnorm((1+conf.level)/2)
	cint <- tanh(cint)
	attr(cint,"conf.level") <- conf.level
	result$conf.int <- cint
    }
  } else {
    TIES <- (min(length(unique(x)),length(unique(y)))<n)
    s <- (n^3-n)*(1-estimate)/6
    if (!TIES) {
	s <- 2*round(s/2)
    }
    statistic <- s
    names(statistic) <- "S"
    p <- if (s>(n^3-n)/6) {
	pspearman::pspearman(s,ifelse(n<10,n,n-k),lower.tail=FALSE,
	  approximation=ifelse(n<10,"AS89","t"))
    } else {
	pspearman::pspearman(s,ifelse(n<10,n,n-k),lower.tail=TRUE,
	  approximation=ifelse(n<10,"AS89","t"))
    }
    pval <- min(2*p,1)
    result$method <- ifelse(!semi,"Spearman's rank partial correlation",
	"Spearman's rank semi-partial correlation")
    names(estimate) <- names(null) <- "rho"
    cor.fun <- function(data,ind) {
	pcor(data[ind,1],data[ind,2],data[ind,3:ncol(data)],method="spearman")
    }
    simul <- boot::boot(tab,cor.fun,R=nrep)
    cint <- .ci(simul$t,conf.level=conf.level)
    attr(cint,"conf.level") <- conf.level
    result$conf.int <- cint
  }
  result$statistic <- statistic
  if (method=="pearson") {result$parameter <- ddl}
  result$p.value <- pval
  result$estimate <- estimate
  result$null.value <- null
  class(result) <- "htest"
  return(result)
}

