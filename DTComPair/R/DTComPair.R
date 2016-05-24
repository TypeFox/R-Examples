# --------------------------------------------------------
# Description: Functions for DTComPair-package
# Author: C. Stock
# Last modified: Feb 15, 2014
# --------------------------------------------------------


# --------------------------------------------------------
# tab.1test
# --------------------------------------------------------
tab.1test <-  function(d, y, data=NULL, testname, ...) {
  # check arguments
  if (missing(d)) stop("Disease status (d) is missing.")
  if (missing(y)) stop("Test result (y) is missing.") 
  if (missing(testname)) testname <- deparse(substitute(y))
  d <- eval(substitute(d), data, parent.frame())
  y <- eval(substitute(y), data, parent.frame())  
  if ((is.numeric(d)==FALSE) & (is.factor(d)==FALSE)) 
    stop("Disease status (d) must be an integer variable.")
  if ((is.numeric(y)==FALSE) & (is.factor(y)==FALSE))
    stop("Test result (y) must be a integer variable.")
  d <- as.integer(d); y <- as.integer(y); 
  if (identical(sort(unique(d)), as.integer(c(0,1)))==FALSE) 
    stop("Disease status (d) must be coded as 0 or 1.")
  if (identical(sort(unique(y)), as.integer(c(0,1)))==FALSE) 
    stop("Test result (y) must be coded as 0 or 1.")
  if ((length(d) != length(y))) stop("Vector lengths differ.") 
  # matrix
  tab.1t <- matrix(rep(0,9), nrow=3, 
                   dimnames=list(c("Test pos.","Test neg.","Total"),
                                 c("Diseased","Non-diseased","Total")))
  data.1t <- data.frame(cbind(d,y))
  tab.1t[1,1] <- nrow(subset(data.1t, d==1 & y==1))
  tab.1t[1,2] <- nrow(subset(data.1t, d==0 & y==1))
  tab.1t[1,3] <- tab.1t[1,1]+tab.1t[1,2]
  tab.1t[2,1] <- nrow(subset(data.1t, d==1 & y==0))
  tab.1t[2,2] <- nrow(subset(data.1t, d==0 & y==0))
  tab.1t[2,3] <- tab.1t[2,1]+tab.1t[2,2]
  tab.1t[3,1] <- tab.1t[1,1]+tab.1t[2,1]
  tab.1t[3,2] <- tab.1t[1,2]+tab.1t[2,2]
  tab.1t[3,3] <- tab.1t[1,3]+tab.1t[2,3]   
  # results
  results <- list(tab.1t, testname)
  names(results) <- c("tab.1test", "testname")    
  class(results) <- "tab.1test"
  return(results)
}

print.tab.1test <- function(x,...) {
  cat(paste("Binary diagnostic test '",x$testname,"'\n\n",sep=''))
  colnames(x$tab.1test)[2] <- "Non-diseased"
  print(x$tab.1test)
}


# --------------------------------------------------------
# read.tab.1test
# --------------------------------------------------------
read.tab.1test <-  function(a, b, c, d, testname, ...) {
  if (missing(testname)) testname <- "Noname"
  tab.1t <- matrix(rep(0,9), nrow=3, 
                   dimnames=list(c("Test pos.","Test neg.","Total"),
                                 c("Diseased","Non-diseased","Total")))
  tab.1t[1,1] <- a; tab.1t[1,2] <- b
  tab.1t[1,3] <- a+b
  tab.1t[2,1] <- c; tab.1t[2,2] <- d
  tab.1t[2,3] <- c+d
  tab.1t[3,1] <- a+c; tab.1t[3,2] <- b+d
  tab.1t[3,3] <- a+b+c+d
  # results
  results <- list(tab.1t, testname)
  names(results) <- c("tab.1test", "testname")    
  class(results) <- "tab.1test"
  return(results)
}


# --------------------------------------------------------
# acc.1test
# --------------------------------------------------------
acc.1test <-  function(tab, alpha, testname, ...) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.1test") stop("Table must be of class 'tab.1test'")
  if (missing(testname)) testname <- tab$testname
  tab <- tab[[1]]
  if (missing(alpha)) alpha <- 0.05
  # sensitivity and specificity
  sens.est <- tab[1,1]/tab[3,1]
  sens.se <- sqrt((tab[1,1]*tab[2,1])/(tab[3,1]^3))
  sens.lcl <- sens.est-qnorm(1-alpha/2)*sens.se
    if (sens.lcl<0) sens.lcl <- 0
  sens.ucl <- sens.est+qnorm(1-alpha/2)*sens.se
    if (sens.ucl>1) sens.ucl <- 1
  sensitivity <- c(sens.est,sens.se,sens.lcl,sens.ucl)
  names(sensitivity) <- c("est","se","lcl","ucl")
  spec.est <- tab[2,2]/tab[3,2]
  spec.se <- sqrt((tab[1,2]*tab[2,2])/(tab[3,2]^3))
  spec.lcl <- spec.est-qnorm(1-alpha/2)*spec.se
    if (spec.lcl<0) spec.lcl <- 0
  spec.ucl <- spec.est+qnorm(1-alpha/2)*spec.se 
    if (spec.ucl>1) spec.ucl <- 1
  specificity <- c(spec.est,spec.se,spec.lcl,spec.ucl)
  names(specificity) <- c("est","se","lcl","ucl")
  # predictive values
  ppv.est <- tab[1,1]/tab[1,3]
  ppv.se <- sqrt((tab[1,1]*tab[1,2])/(tab[1,3]^3))
  ppv.lcl <- ppv.est-qnorm(1-alpha/2)*ppv.se
    if (ppv.lcl<0) ppv.lcl <- 0
  ppv.ucl <- ppv.est+qnorm(1-alpha/2)*ppv.se
    if (ppv.ucl>1) ppv.ucl <- 1
  ppv <- c(ppv.est,ppv.se,ppv.lcl,ppv.ucl)
  names(ppv) <- c("est","se","lcl","ucl")
  npv.est <- tab[2,2]/tab[2,3]
  npv.se <- sqrt((tab[2,1]*tab[2,2])/(tab[2,3]^3))
  npv.lcl <- npv.est-qnorm(1-alpha/2)*npv.se
    if (npv.lcl<0) npv.lcl <- 0
  npv.ucl <- npv.est+qnorm(1-alpha/2)*npv.se
    if (npv.ucl>1) npv.ucl <- 1
  npv <- c(npv.est,npv.se,npv.lcl,npv.ucl)
  names(npv) <- c("est","se","lcl","ucl")
  # diagnostic likelihood ratios
  pdlr.est <- sens.est/(1-spec.est)
  pdlr.se.log <- sqrt(((1-sens.est)/(tab[1,1]))+
                      (spec.est/tab[1,2]))
  pdlr.lcl <- exp(log(pdlr.est)-qnorm(1-alpha/2)*pdlr.se.log)
  pdlr.ucl <- exp(log(pdlr.est)+qnorm(1-alpha/2)*pdlr.se.log)
  pdlr <- c(pdlr.est,pdlr.se.log,pdlr.lcl,pdlr.ucl)
  names(pdlr) <- c("est","se.ln","lcl","ucl")
  ndlr.est <- (1-sens.est)/spec.est
  ndlr.se.log <- sqrt((sens.est/tab[2,1])+((1-spec.est)/
                 tab[2,2]))
  ndlr.lcl <- exp(log(ndlr.est)-qnorm(1-alpha/2)*ndlr.se.log)
  ndlr.ucl <- exp(log(ndlr.est)+qnorm(1-alpha/2)*ndlr.se.log)
  ndlr <- c(ndlr.est,ndlr.se.log,ndlr.lcl,ndlr.ucl)
  names(ndlr) <- c("est","se.ln","lcl","ucl")
  # results
  results <- list(tab, sensitivity, specificity,
                  ppv, npv, pdlr, ndlr, alpha, testname)
  names(results) <- c("tab", "sensitivity", "specificity",
                      "ppv", "npv", "pdlr", "ndlr", "alpha", "testname")    
  class(results) <- "acc.1test"
  return(results)
}

print.acc.1test <- function(x,...) {
  cat(paste("Diagnostic accuracy of test '",x$testname,"'\n",sep=''))
  cat(paste("\n(Estimates, standard errors and ",
            100*(1-x$alpha),
            "%-confidence intervals)\n\n",sep=""))
  acc.mat1 <- matrix(data=c(x$sensitivity[1:4],
                            x$specificity[1:4],
                            x$ppv[1:4], x$npv[1:4]), 
                     nrow=4, ncol=4, byrow=TRUE,
                     dimnames = list(c("Sensitivity", "Specificity",
                                       "PPV", "NPV"),
                                     c("Est.", "SE", 
                                       "Lower CL", "Upper CL")))  
  print(acc.mat1)
  cat("\n")
  acc.mat2 <- matrix(data=c(x$pdlr[1:4],x$ndlr[1:4]), 
                     nrow=2, ncol=4, byrow=TRUE,
                     dimnames = list(c("PDLR ","NDLR "),
                                     c("Est.", "SE (log)", 
                                       "Lower CL", "Upper CL")))
  print(acc.mat2)   
}


# --------------------------------------------------------
# tab.paired
# --------------------------------------------------------
tab.paired <-  function(d, y1, y2, data=NULL, testnames, ...) {
  # check arguments
  if (missing(d)) stop("Disease status (d) is missing.")
  if (missing(y1)) stop("Test result (y1) is missing.") 
  if (missing(y2)) stop("Test result (y2) is missing.")
  if (missing(testnames)) testnames <- c(deparse(substitute(y1)),
                                         deparse(substitute(y2)))
  d <- eval(substitute(d), data, parent.frame())
  y1 <- eval(substitute(y1), data, parent.frame())
  y2 <- eval(substitute(y2), data, parent.frame())  
  if ((is.numeric(d)==FALSE) & (is.factor(d)==FALSE)) 
    stop("Disease status (d) must be a numeric or factor variable.")
  if ((is.numeric(y1)==FALSE) & (is.factor(y1)==FALSE))
    stop("Test result (y1) must be a numeric or factor variable.")
  if ((is.numeric(y2)==FALSE) & (is.factor(y2)==FALSE))
    stop("Test result (y2) must be a numeric or factor variable.")
  d <- as.integer(d); y1 <- as.integer(y1); y2 <- as.integer(y2);
  if (identical(sort(unique(d)), as.integer(c(0,1)))==FALSE) 
    stop("Disease status (d) must be coded as 0 or 1.")
  if (identical(sort(unique(y1)), as.integer(c(0,1)))==FALSE) 
    stop("Test result (y1) must be coded as 0 or 1.")
  if (identical(sort(unique(y2)), as.integer(c(0,1)))==FALSE) 
    stop("Test result (y2) must be coded as 0 or 1.")
  if ((length(d) != length(y1)) | (length(y1) != length(y2)) ) 
    stop("Vector lengths differ.")
  ## matrices
  data.paired <- data.frame(d,y1,y2)
  # diseased
  tab.d <- matrix(rep(0,9), nrow=3, dimnames=list(
                    c("Test2 pos.","Test2 neg.","Total"),
                    c("Test1 pos.","Test1 neg.","Total")))
  tab.d[1,1] <- nrow(subset(data.paired, d==1 & y1==1 & y2==1))
  tab.d[1,2] <- nrow(subset(data.paired, d==1 & y1==0 & y2==1))
  tab.d[1,3] <- tab.d[1,1]+tab.d[1,2]
  tab.d[2,1] <- nrow(subset(data.paired, d==1 & y1==1 & y2==0))
  tab.d[2,2] <- nrow(subset(data.paired, d==1 & y1==0 & y2==0))
  tab.d[2,3] <- tab.d[2,1]+tab.d[2,2]
  tab.d[3,1] <- tab.d[1,1]+tab.d[2,1]
  tab.d[3,2] <- tab.d[1,2]+tab.d[2,2]
  tab.d[3,3] <- tab.d[1,3]+tab.d[2,3]
  # non-diseased
  tab.nd <- matrix(rep(0,9), nrow=3, dimnames=list(
                    c("Test2 pos.","Test2 neg.","Total"),
                    c("Test1 pos.","Test1 neg.","Total")))
  tab.nd[1,1] <- nrow(subset(data.paired, d==0 & y1==1 & y2==1))
  tab.nd[1,2] <- nrow(subset(data.paired, d==0 & y1==0 & y2==1))
  tab.nd[1,3] <- tab.nd[1,1]+tab.nd[1,2]
  tab.nd[2,1] <- nrow(subset(data.paired, d==0 & y1==1 & y2==0))
  tab.nd[2,2] <- nrow(subset(data.paired, d==0 & y1==0 & y2==0))
  tab.nd[2,3] <- tab.nd[2,1]+tab.nd[2,2]
  tab.nd[3,1] <- tab.nd[1,1]+tab.nd[2,1]
  tab.nd[3,2] <- tab.nd[1,2]+tab.nd[2,2]
  tab.nd[3,3] <- tab.nd[1,3]+tab.nd[2,3]
  # results
  results <- list(tab.d, tab.nd, testnames)
  names(results) <- c("diseased","non.diseased","testnames")    
  class(results) <- "tab.paired"
  return(results)
}

print.tab.paired <- function(x,...) {
  cat("Two binary diagnostic tests (paired design)")
  cat("\n\n")
  cat("Test1: '",x$testnames[1],"'\n",
      "Test2: '",x$testnames[2],"'\n\n", sep="")
  cat("Diseased:\n") 
  print(x$diseased)
  cat("\n")
  cat("Non-diseased:\n")
  print(x$non.diseased)
  cat("\n")
}


# --------------------------------------------------------
# read.tab.paired
# --------------------------------------------------------
read.tab.paired <-  function(d.a, d.b, d.c, d.d,
                             nd.a, nd.b, nd.c, nd.d,
                             testnames, ...) {
  if (missing(testnames)) testnames <- c("Noname 1","Noname 2")
  # diseased
  tab.d <- matrix(rep(0,9), nrow=3, dimnames=list(
    c("Test2 pos.","Test2 neg.","Total"),
    c("Test1 pos.","Test1 neg.","Total")))
  tab.d[1,1] <- d.a; tab.d[1,2] <- d.b
  tab.d[1,3] <- d.a+d.b
  tab.d[2,1] <- d.c; tab.d[2,2] <- d.d
  tab.d[2,3] <- d.c+d.d
  tab.d[3,1] <- d.a+d.c; tab.d[3,2] <- d.b+d.d
  tab.d[3,3] <- d.a+d.b+d.c+d.d
  # non-diseased
  tab.nd <- matrix(rep(0,9), nrow=3, dimnames=list(
    c("Test2 pos.","Test2 neg.","Total"),
    c("Test1 pos.","Test1 neg.","Total")))
  tab.nd[1,1] <- nd.a; tab.nd[1,2] <- nd.b
  tab.nd[1,3] <- nd.a+nd.b
  tab.nd[2,1] <- nd.c; tab.nd[2,2] <- nd.d
  tab.nd[2,3] <- nd.c+nd.d
  tab.nd[3,1] <- nd.a+nd.c; tab.nd[3,2] <- nd.b+nd.d
  tab.nd[3,3] <- nd.a+nd.b+nd.c+nd.d   
  # results
  results <- list(tab.d, tab.nd, testnames)
  names(results) <- c("diseased","non.diseased","testnames")    
  class(results) <- "tab.paired"
  return(results)
}


# --------------------------------------------------------
# generate.paired
# --------------------------------------------------------
generate.paired <- function(tab, ...) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  testnames <- tab$testnames
  # generate dataframe
  df <- expand.grid(d=c(1,0), y1=c(1,0), y2=c(1,0))
  df <- df[with(df, order(-d,-y1,-y2)), ]
  n <- c(tab$diseased[1,1], tab$diseased[2,1],
         tab$diseased[1,2], tab$diseased[2,2],
         tab$non.diseased[1,1], tab$non.diseased[2,1],
         tab$non.diseased[1,2], tab$non.diseased[2,2])
  df <- as.data.frame(cbind(df,n))
  df <- df[rep(seq(dim(df)[1]), df$n),-4]
  rownames(df) <- NULL
  return(df)
}


# --------------------------------------------------------
# acc.paired
# --------------------------------------------------------
acc.paired <-  function(tab, alpha, ...) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") stop("Table must be of class 'tab.paired'")
  if (missing(alpha)) alpha <- 0.05
  # tables for each test
  test1 <- read.tab.1test(tab$diseased[3,1], tab$non.diseased[3,1],
                          tab$diseased[3,2], tab$non.diseased[3,2], 
                          testname=tab$testnames[1])
  test2 <- read.tab.1test(tab$diseased[1,3], tab$non.diseased[1,3],
                          tab$diseased[2,3], tab$non.diseased[2,3], 
                          testname=tab$testnames[2])
  # accuracy of each test
  acc.test1 <- acc.1test(test1)
  acc.test2 <- acc.1test(test2)
  # results
  results <- list(acc.test1, acc.test2)
  names(results) <- c("Test1","Test2")
  class(results) <- "acc.paired"
  return(results)
}

print.acc.paired <- function(x,...) {
  print(x[[1]]); 
  cat("\n----------------------------------------------------------\n")
  print(x[[2]])
}


# --------------------------------------------------------
# represent.long
# --------------------------------------------------------
represent.long <- function(d, y1, y2) {
  df <- data.frame(d, y1, y2)
  colnames(df) <- c("d","y1","y2")
  df$id <- 1:nrow(df)
  dt1 <- df[,c(4,1,2)]; colnames(dt1)[3] <- "y"
  dt1$x <- rep(1,nrow(dt1))
  dt2 <- df[,c(4,1,3)]; colnames(dt2)[3] <- "y"
  dt2$x <- rep(0,nrow(dt2))
  df <- as.data.frame(rbind(dt1, dt2))
  df <- df[with(df, order(id)), ]
  df <- df[,c(1,2,4,3)]
  row.names(df) <- 1:nrow(df)
  return(df)
}


# --------------------------------------------------------
# sesp.mcnemar
# --------------------------------------------------------
sesp.mcnemar <- function(tab) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  # accuracy
  acc <- acc.paired(tab)
  # sensitivity
  se.1 <- acc$Test1$sensitivity["est"]; se.2 <- acc$Test2$sensitivity["est"]
  names(se.1) <- NULL; names(se.2) <- NULL
  diff.sens <- (se.2-se.1); names(diff.sens) <- NULL
  b <- tab$diseased[1,2]; c <- tab$diseased[2,1]
  X2 <- ((b-c)^2)/(b+c)
  p.value <- 1-pchisq(X2, df=1)
  sensitivity <- list(se.1, se.2, diff.sens, X2, p.value)
  # specificity
  sp.1 <- acc$Test1$specificity["est"]; sp.2 <- acc$Test2$specificity["est"]
  names(sp.1) <- NULL; names(sp.2) <- NULL
  diff.spec <- (sp.2-sp.1); names(diff.spec) <- NULL
  b <- tab$non.diseased[1,2]; c <- tab$non.diseased[2,1]
  X2 <- ((b-c)^2)/(b+c)
  p.value <- 1-pchisq(X2, df=1)
  specificity <- list(sp.1, sp.2, diff.spec, X2, p.value)
  # results
  method <- "mcnemar"
  results <- list(sensitivity, specificity, method)
  names(results) <- c("sensitivity","specificity","method")
  names(results$sensitivity) <- c("test1","test2","diff","test.statistic","p.value")
  names(results$specificity) <- c("test1","test2","diff","test.statistic","p.value")
  return(results)
}


# --------------------------------------------------------
# sesp.exactbinom
# --------------------------------------------------------
sesp.exactbinom <- function(tab) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  # accuracy
  acc <- acc.paired(tab)
  # sensitivity
  se.1 <- acc$Test1$sensitivity["est"]; se.2 <- acc$Test2$sensitivity["est"]
  names(se.1) <- NULL; names(se.2) <- NULL
  diff.sens <- (se.2-se.1); names(diff.sens) <- NULL
  m <- tab$diseased[1,2] + tab$diseased[2,1]
  k <- min(tab$diseased[1,2], tab$diseased[2,1])
  csum <- 0; for (j in 0:k) csum <- csum+choose(m,j)
  p.value <- 2*csum*(0.5^m)
  sensitivity <- list(se.1,se.2,diff.sens,p.value)
  # specificity
  sp.1 <- acc$Test1$specificity["est"]; sp.2 <- acc$Test2$specificity["est"]
  names(sp.1) <- NULL; names(sp.2) <- NULL
  diff.spec <- (sp.2-sp.1); names(diff.spec) <- NULL
  m <- tab$non.diseased[1,2] + tab$non.diseased[2,1]
  k <- min(tab$non.diseased[1,2], tab$non.diseased[2,1])
  csum <- 0; for (j in 0:k) csum <- csum+choose(m,j)
  p.value <- 2*csum*(0.5^m)
  specificity <- list(sp.1,sp.2,diff.spec,p.value)
  # results
  method <- "exactbinom"
  results <- list(sensitivity,specificity,method) 
  names(results) <- c("sensitivity","specificity","method")
  names(results$sensitivity) <- c("test1","test2","diff","p.value")
  names(results$specificity) <- c("test1","test2","diff","p.value")
  return(results)
}


# --------------------------------------------------------
# pv.gs
# --------------------------------------------------------
pv.gs <- function(tab) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  # accurac
  acc <- acc.paired(tab)
  ## ppv
  ppv.1 <- acc$Test1$ppv["est"]; ppv.2 <- acc$Test2$ppv["est"]
  names(ppv.1) <- NULL; names(ppv.2) <- NULL
  diff.ppv <- abs(ppv.1-ppv.2); names(diff.ppv) <- NULL
  # proportion of positive tests of type 2
  z.bar <- sum(c(tab$diseased[1,3], tab$non.diseased[1,3])) /
           sum(c(tab$diseased[1,c(1,3)], tab$diseased[2,1],
                 tab$non.diseased[1,c(1,3)], tab$non.diseased[2,1]))
  # all positive  tests in diseased subjects / all positive tests
  d.bar <- sum(c(tab$diseased[1,c(1,1,2)], tab$diseased[2,1])) /
           sum(c(tab$diseased[1,c(1,1,2)], tab$diseased[2,1],
                 tab$non.diseased[1,c(1,1,2)], tab$non.diseased[2,1]))
  numerator <- (tab$diseased[1,1]*(1-2*z.bar) + 
                tab$diseased[1,2]*(1-z.bar) + 
                tab$diseased[2,1]*(0-z.bar))^2 
  denominator <- (1-d.bar)^2 * 
                 (tab$diseased[1,1]*(1-2*z.bar)^2 + 
                  tab$diseased[1,2]*(1-z.bar)^2 + 
                  tab$diseased[2,1]*(0-z.bar)^2 ) +
                 (0-d.bar)^2 * 
                  (tab$non.diseased[1,1]*(1-2*z.bar)^2 + 
                   tab$non.diseased[1,2]*(1-z.bar)^2 + 
                   tab$non.diseased[2,1]*(0-z.bar)^2 )
  t.ppv <- numerator/denominator
  p.value <- 1-pchisq(t.ppv, df=1) 
  ppv <- list(ppv.1, ppv.2, diff.ppv, t.ppv, p.value)
  ## npv
  npv.1 <- acc$Test1$npv["est"]; npv.2 <- acc$Test2$npv["est"]
  names(npv.1) <- NULL; names(npv.2) <- NULL
  diff.npv <- abs(npv.1-npv.2); names(diff.npv) <- NULL
  # proportion of negative tests of type 2
  z.bar <- sum(c(tab$diseased[2,3], tab$non.diseased[2,3])) /
           sum(c(tab$diseased[2,c(2,3)], tab$diseased[1,2],
                 tab$non.diseased[2,c(2,3)], tab$non.diseased[1,2]))
  # all negative tests in non-diseased subjects / all negative tests
  d.bar <- sum(c(tab$non.diseased[2,c(2,3)], tab$non.diseased[1,2])) /
           sum(c(tab$diseased[2,c(2,3)], tab$diseased[1,2],
                 tab$non.diseased[2,c(2,3)], tab$non.diseased[1,2]))
  numerator <- (tab$non.diseased[2,2]*(1-2*z.bar) + 
                tab$non.diseased[2,1]*(1-z.bar) + 
                tab$non.diseased[1,2]*(0-z.bar))^2 
  denominator <- (1-d.bar)^2 * 
                  (tab$non.diseased[2,2]*(1-2*z.bar)^2 + 
                   tab$non.diseased[2,1]*(1-z.bar)^2 + 
                   tab$non.diseased[1,2]*(0-z.bar)^2 ) +
                 (0-d.bar)^2 * 
                  (tab$diseased[2,2]*(1-2*z.bar)^2 + 
                   tab$diseased[2,1]*(1-z.bar)^2 + 
                   tab$diseased[1,2]*(0-z.bar)^2 )
  t.npv <- numerator/denominator
  p.value <- 1-pchisq(t.npv, df=1)  
  npv <- list(npv.1, npv.2, diff.npv, t.npv, p.value)
  # results
  method <- "generalized score statistic (gs)"
  results <- list(ppv,npv,method)
  names(results) <- c("ppv","npv","method")
  names(results$ppv) <- c("test1","test2","diff","test.statistic","p.value")
  names(results$npv) <- c("test1","test2","diff","test.statistic","p.value")
  return(results)
}


# --------------------------------------------------------
# pv.wgs
# --------------------------------------------------------
pv.wgs <- function(tab) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  acc <- acc.paired(tab)
  ## ppv
  ppv.1 <- acc$Test1$ppv["est"]; ppv.2 <- acc$Test2$ppv["est"]
  names(ppv.1) <- NULL; names(ppv.2) <- NULL
  diff.ppv <- abs(ppv.1-ppv.2); names(diff.ppv) <- NULL
  ppv.pooled <- (tab$diseased[1,1]*2 + tab$diseased[1,2] + tab$diseased[2,1]) /
    (tab$diseased[1,3] + tab$non.diseased[1,3] + tab$diseased[3,1] + tab$non.diseased[3,1])
  numerator <- diff.ppv**2
  c.p.ppv <- (tab$diseased[1,1]*(1-ppv.pooled)**2 + tab$non.diseased[1,1]*(ppv.pooled**2)) /  
    (tab$diseased[1,3] + tab$non.diseased[1,3] + tab$diseased[3,1] + tab$non.diseased[3,1]) 
  denominator <- (ppv.pooled*(1-ppv.pooled) - 2*c.p.ppv) * 
    (( 1/ (tab$diseased[1,3] + tab$non.diseased[1,3])) +
       ( 1/ (tab$diseased[3,1] + tab$non.diseased[3,1]))   )
  t.ppv <- numerator/denominator
  p.value <- 1-pchisq(t.ppv, df=1)  
  ppv <- list(ppv.1, ppv.2, diff.ppv, t.ppv, p.value)
  ## npv
  npv.1 <- acc$Test1$npv["est"]; npv.2 <- acc$Test2$npv["est"]
  names(npv.1) <- NULL; names(npv.2) <- NULL  
  diff.npv <- abs(npv.1-npv.2); names(diff.npv) <- NULL  
  npv.pooled <- (tab$non.diseased[2,2]*2 + tab$non.diseased[1,2] + tab$non.diseased[2,1]) /
    (tab$diseased[2,3] + tab$non.diseased[2,3] + tab$diseased[3,2] + tab$non.diseased[3,2])
  numerator <- diff.npv**2
  c.p.npv <- (tab$diseased[2,2]*(npv.pooled)**2 + tab$non.diseased[2,2]*(1-npv.pooled)**2) /  
    (tab$diseased[2,3] + tab$non.diseased[2,3] + tab$diseased[3,2] + tab$non.diseased[3,2]) 
  denominator <- (npv.pooled*(1-npv.pooled) - 2*c.p.npv) * 
    (( 1/ (tab$diseased[2,3] + tab$non.diseased[2,3])) +
       ( 1/ (tab$diseased[3,2] + tab$non.diseased[3,2])))
  t.npv <- numerator/denominator
  p.value <- 1-pchisq(t.npv, df=1)  
  npv <- list(npv.1, npv.2, diff.npv, t.npv, p.value)  
  # results
  method <- "weighted generalized score statistic (wgs)"
  results <- list(ppv,npv,method)
  names(results) <- c("ppv","npv","method")
  names(results$ppv) <- c("test1","test2","diff","test.statistic","p.value")
  names(results$npv) <- c("test1","test2","diff","test.statistic","p.value")
  return(results)
}


# --------------------------------------------------------
# pv.rpv
# --------------------------------------------------------
pv.rpv <- function(tab, alpha) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  if (missing(alpha)) alpha <- 0.05
  # pre-processing 
  N <-  tab$non.diseased[3,3] + tab$diseased[3,3]
  p1 <- tab$non.diseased[1,1] / N
  p2 <- tab$non.diseased[1,2] / N
  p3 <- tab$non.diseased[2,1] / N
  p4 <- tab$non.diseased[2,2] / N
  p5 <- tab$diseased[1,1] / N
  p6 <- tab$diseased[1,2] / N
  p7 <- tab$diseased[2,1] / N
  p8 <- tab$diseased[2,2] / N
  acc <- acc.paired(tab)
  # rppv
  ppv.1 <- acc$Test1$ppv["est"]; ppv.2 <- acc$Test2$ppv["est"]
  names(ppv.1) <- NULL; names(ppv.2) <- NULL
  rel.ppv <- ppv.1/ppv.2; names(rel.ppv) <- NULL
  sigma2.p <- (1/((p5+p7)*(p5+p6))) *
    (p6*(1-ppv.1) + p5*(ppv.1-ppv.2) +
       + 2*(p7+p3)*ppv.2*ppv.1 + p7*(1-3*ppv.2)  )
  se.log.rel.ppv <- sqrt(sigma2.p/N)
  lcl <- exp(log(rel.ppv) - qnorm(1-alpha/2)*se.log.rel.ppv)
  ucl <- exp(log(rel.ppv) + qnorm(1-alpha/2)*se.log.rel.ppv)
  t.ppv <- log(rel.ppv) / se.log.rel.ppv
  p.value <- 2*pnorm(-abs(t.ppv))
  ppv <- list(ppv.1, ppv.2, rel.ppv, se.log.rel.ppv, lcl, ucl, t.ppv, p.value)
  # rnpv
  npv.1 <- acc$Test1$npv["est"]; npv.2 <- acc$Test2$npv["est"]
  names(npv.1) <- NULL; names(npv.2) <- NULL  
  rel.npv <- npv.1/npv.2; names(rel.ppv) <- NULL
  sigma2.n <- (1/((p2+p4)*(p3+p4))) *
    ( npv.1*(-p3+p4-2*(p4+p8)*npv.2) + 
        (p2+p3) - npv.2*(p2-p4) )
  se.log.rel.npv <- sqrt(sigma2.n/N)
  lcl <- exp(log(rel.npv) - qnorm(1-alpha/2)*se.log.rel.npv)
  ucl <- exp(log(rel.npv) + qnorm(1-alpha/2)*se.log.rel.npv)
  t.npv <- log(rel.npv) / se.log.rel.npv
  p.value <- 2*pnorm(-abs(t.npv))
  npv <- list(npv.1, npv.2, rel.npv, se.log.rel.npv, lcl, ucl, t.npv, p.value)
  # results
  method <- "relative predictive values (rpv)"
  results <- list(ppv,npv,method,alpha)
  names(results) <- c("ppv","npv","method","alpha")
  names(results$ppv) <- c("test1","test2","rppv","se.log.rppv","lcl.rppv","ucl.rppv","test.statistic","p.value")
  names(results$npv) <- c("test1","test2","rnpv","se.log.rnpv","lcl.rnpv","ucl.rnpv","test.statistic","p.value")
  return(results)
}


# --------------------------------------------------------
# End 
# --------------------------------------------------------
