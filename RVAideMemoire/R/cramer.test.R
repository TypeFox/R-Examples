# boot : boot

cramer.test <-
function (x,y,nrep=1000,conf.level=0.95) {
  if (is.matrix(x)) {
    if (is.null(rownames(x))) {rownames(x) <- letters[1:nrow(x)]}
    if (is.null(colnames(x))) {colnames(x) <- LETTERS[1:ncol(x)]}
    var1 <- rep(colnames(x),colSums(x))
    var2 <- rep(rep(rownames(x),ncol(x)),as.vector(x))
  } else {
    var1 <- x
    var2 <- y
  }
  if (length(var1)!=length(var2)) {
    stop(paste("'",deparse(substitute(var1)),"' and '",deparse(substitute(var2)),"' lengths differ",sep=""))
  }
  data.name <- if (is.matrix(x)) {
    deparse(substitute(x))
  } else {
    paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="")
  }
  if (!is.factor(var1)) {var1 <- factor(var1)}
  if (!is.factor(var2)) {var2 <- factor(var2)}
  nul <- as.numeric(row.names(table(c(which(is.na(var1)), which(is.na(var2))))))
  var1.2 <- if (length(nul)>0) {
    var1[-nul]
  } else {
    var1
  }
  var2.2 <- if (length(nul)>0) {
    var2[-nul]
  } else {
    var2
  }
  if (any(tapply(var1.2,var1.2,function(x) length(x)/length(var1.2))<0.05) |
    any(tapply(var2.2,var2.2,function(x) length(x)/length(var2.2))<0.05)) {
    warning("at least 1 level contains less than 5% of total number of individuals")
  }
  tab.cont <- table(var1.2,var2.2)
  v <- sqrt(as.numeric(suppressWarnings(chisq.test(tab.cont)$statistic))/(sum(tab.cont)*(min(dim(tab.cont))-1)))
  names(v) <- "V"
  v.fun <- function(dat,ind) {
    cont <- table(dat[ind,1],dat[ind,2])
    sqrt(as.numeric(suppressWarnings(chisq.test(cont)$statistic))/(sum(cont)*(min(dim(cont))-1)))
  }
  simul <- boot::boot(data.frame(var1.2, var2.2),v.fun,R=nrep)
  cint <- .ci(simul$t,conf.level=conf.level)
  attr(cint, "conf.level") <- conf.level
  test <- suppressWarnings(chisq.test(tab.cont))
  nval <- 0
  names(nval) <- "association"
  result <- list(method="Cram\u00E9r's association coefficient",statistic=test$statistic,parameter=test$parameter,
    p.value=test$p.value,data.name=data.name,estimate=v,conf.level=conf.level,rep=nrep,conf.int=cint,
    alternative="two.sided",null.value=nval)
  class(result) <- "htest"
  return(result)
}
