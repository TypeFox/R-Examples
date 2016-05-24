multtest.gp <- function(tab,fac,test=c("param","perm","rank"),transform=c("none","sqrt","log"),
  add=0,p.method="fdr",ordered=TRUE,...) {
  test <- match.arg(test)
  transform <- match.arg(transform)
  tab <- as.data.frame(tab)
  nlev <- nlevels(fac)
  if (nlev<2) {stop("at least 2 groups are needed")}
  gp.mean <- as.matrix(t(aggregate(as.matrix(tab)~fac,FUN=mean,na.rm=TRUE))[-1,])
  colnames(gp.mean) <- paste0("Mean.",abbreviate(levels(fac),1))
  mode(gp.mean) <- "numeric"
  gp.se <- as.matrix(t(aggregate(as.matrix(tab)~fac,FUN=se))[-1,])
  colnames(gp.se) <- paste0("SE.",abbreviate(levels(fac),1))
  mode(gp.se) <- "numeric"
  if (test %in% c("param","perm")) {
    if (transform=="sqrt") {tab <- apply(tab,2,sqrt)}
    if (transform=="log") {tab <- apply(tab,2,function(x) log(x+add))}    
  }
  test <- if (test=="param") {
    if (nlev>2) {"anova"} else {"t"}
  } else if (test=="perm") {
    if (nlev>2) {"perm.anova"} else {"perm.t"}
  } else {
    if (nlev>2) {"kw"} else {"mww"}
  }
  test.f <- switch(test,anova=multtest.gp.anova,t=multtest.gp.t,perm.anova=multtest.gp.perm.anova,
    perm.t=multtest.gp.perm.t,kw=multtest.gp.kw,mww=multtest.gp.mww)
  tab.res <- test.f(tab,fac,...)
  tab.res$P.value <- p.adjust(tab.res$P.value,method=p.method)
  tab.res[,3] <- integer(ncol(tab))
  tab.res <- cbind(tab.res,gp.mean,gp.se)
  tab.res <- signif(tab.res,5)
  tab.res <- as.data.frame(tab.res)
  tab.res$V3 <- .psignif(tab.res$P.value)
  colnames(tab.res)[3] <- " "
  if (ordered) {tab.res <- tab.res[order(tab.res$P.value),]}
  res <- list(tab=tab.res,p.method=p.method,labs=levels(fac))
  class(res) <- c("multtest","multtest.gp","list")
  return(res)
}

multtest.gp.anova <- function(tab,fac,...) {
  nvar <- ncol(tab)
  lab <- colnames(tab)
  res <- data.frame(F=integer(nvar),P.value=integer(nvar),row.names=lab)
  for (i in 1:ncol(tab)) {
    x <- tab[,i]
    mod <- lm(x~fac)
    test <- anova(mod)
    res[i,] <- test[1,c("F value","Pr(>F)")]
  }
  return(res)
}

multtest.gp.t <- function(tab,fac,...) {
  nvar <- ncol(tab)
  lab <- colnames(tab)
  res <- data.frame(t=integer(nvar),P.value=integer(nvar),row.names=lab)
  for (i in 1:ncol(tab)) {
    x <- tab[,i]
    test <- t.test(x~fac,...)
    res[i,] <- c(test$statistic,test$p.value)
  }
  return(res)
}

multtest.gp.perm.anova <- function(tab,fac,...) {
  nvar <- ncol(tab)
  lab <- colnames(tab)
  res <- data.frame(F=integer(nvar),P.value=integer(nvar),row.names=lab)
  for (i in 1:ncol(tab)) {
    x <- tab[,i]
    test <- perm.anova(x~fac,...)
    res[i,] <- c(as.numeric(test[1,"F value"]),as.numeric(test[1,"Pr(>F)"]))
  }
  return(res)
}

multtest.gp.perm.t <- function(tab,fac,...) {
  nvar <- ncol(tab)
  lab <- colnames(tab)
  res <- data.frame(t=integer(nvar),P.value=integer(nvar),row.names=lab)
  for (i in 1:ncol(tab)) {
    x <- tab[,i]
    test <- perm.t.test(x~fac,...)
    res[i,] <- c(test$statistic,test$p.value)
  }
  return(res)
}

multtest.gp.kw <- function(tab,fac,...) {
  nvar <- ncol(tab)
  lab <- colnames(tab)
  res <- data.frame(Chi2=integer(nvar),P.value=integer(nvar),row.names=lab)
  for (i in 1:ncol(tab)) {
    x <- tab[,i]
    test <- kruskal.test(x~fac)
    res[i,] <- c(test$statistic,test$p.value)
  }
  return(res)
}

multtest.gp.mww <- function(tab,fac,...) {
  nvar <- ncol(tab)
  lab <- colnames(tab)
  res <- data.frame(W=integer(nvar),P.value=integer(nvar),row.names=lab)
  for (i in 1:ncol(tab)) {
    x <- tab[,i]
    test <- wilcox.test(x~fac,...)
    res[i,] <- c(test$statistic,test$p.value)
  }
  return(res)
}

print.multtest <- function(x,...) {
  print(x$tab)
  cat(paste("\nP value adjustment method:",x$p.method,"\n"))
}

plot.multtest.gp <- function(x,signif=FALSE,alpha=0.05,xlab="Group",ylab="Mean (+/- SE) value",
  titles=NULL,groups=NULL,...) {
  rows <- if (signif) {
    which(x$tab$P.value<=alpha)
  } else {
    1:nrow(x$tab)
  }
  tab2 <- x$tab[rows,]
  n <- length(rows)
  col.m <- 4:(3+length(x$labs))
  col.s <- (4+length(x$labs)):(3+2*length(x$labs))
  labs <- if (is.null(groups)) {x$labs} else {groups}
  par(mfrow=n2mfrow(n))
  for (i in 1:n) {
    m <- unlist(tab2[i,col.m])
    names(m) <- labs
    s <- unlist(tab2[i,col.s])
    ymin <- ifelse(any(m-s<0),1.3*min(m-s),0)
    ymax <- ifelse(any(m+s>0),1.3*max(m+s),0)
    g <- barplot(m,main=ifelse(is.null(titles),rownames(tab2)[i],titles[i]),xlab=xlab,ylab=ylab,
	ylim=c(ymin,ymax),...)
    arrows(g,m-s,g,m+s,code=3,angle=90,length=0.06)
    ypval <- ifelse(any(m+s>0),1.2*max(m+s),1.2*min(m-s))
    P <- ifelse(tab2$P.value[i]<0.0001,"P < 0.0001",paste("P =",tab2$P.value[i]))
    text(mean(g),ypval,P)
  }
}


