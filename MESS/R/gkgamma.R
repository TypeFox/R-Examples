gkgamma <- function(x, conf.level = 0.95) {

  if (! (is.matrix(as.matrix(x)))) {
    stop("x must be a table or matrix")
  }
  x <- as.matrix(x)
  
  DNAME <- deparse(substitute(x))
  
  rows <- nrow(x)
  cols <- ncol(x)
  n <- sum(x)
  
  con <- x
  dis <- x
  
  rseq <- 1:rows
  cseq <- 1:cols
  
  for (i in 1:rows) {
    for (j in 1:cols) {            
      con[i,j] <- sum(x[rseq<i,cseq<j]) + sum(x[rseq>i,cseq>j])       
      dis[i,j] <- sum(x[rseq>i,cseq<j]) + sum(x[rseq>i,cseq<j])     
    }
  }  
  CC <- sum(x*con)/2
  DC <- sum(x*dis)/2
  
  ESTIMATE <- (CC-DC)/(CC+DC)
  
  se1 <- sum(x*(DC*con - CC*dis)^2)*16/(CC+DC)^4
  se0 <- (sum(x*(con-dis)^2) - (CC-DC)^2/n)*4/(CC+DC)^2
  
  STATISTIC <- ESTIMATE / se0
  
  PVAL <- 2*pnorm(-abs(STATISTIC))
  CINT <- ESTIMATE + c(-1,1)*qnorm(conf.level)*se1
  if (!is.null(CINT)) 
    attr(CINT, "conf.level") <- conf.level
  
  METHOD <- "Goodman-Kruskal's gamma for ordinal categorical data"
  names(STATISTIC) <- "Z"
  names(ESTIMATE) <- "Goodman-Kruskal's gamma"
#  names(CINT) <- "Confidence "
#  names(PARAMETER) <- "df"
#  if (any(E < 5) && is.finite(PARAMETER)) 
#    warning("Chi-squared approximation may be incorrect")
  structure(list(statistic = STATISTIC, 
#                 parameter = PARAMETER, 
                 p.value = PVAL,
                 estimate = ESTIMATE,
                 conf.int = CINT,
                 method = METHOD, 
                 data.name = DNAME, 
                 observed = x, 
                 se0 = se0,
                 se1 = se1),
                 class = "htest")
  
    }
