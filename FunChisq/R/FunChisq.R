# FunChisq.R -- statistical tests for nonparametric functional dependencies
#
# YZ, HZ, MS
# Modified: Feb 8, 2015; June 25, 2015
# Jan 22, 2016 MS:
#    1. introduced the type argument to indicate whether alternative hypothesis is functional
#       or non-constant function
#    2. introduced the log.p argument to obtain the log of p-value
# Jan 23, 2016 MS:
#    Added a return of an estimate of functional index between 0 and 1. It is asymetrical and
#      different from Cramer's V.
# Jan 30, 2016 MS:
#    Added a new argument "index.kind" to specify the function index kind: "unconditional"
#      or "conditional" given marginal of Y
# Feb 5, 2016 MS:
#    Renamed the "type" argument to "alternative"
# Feb 7, 2016 MS:
#   Handled the special case when the degrees of freedom for the normalized FunChisq are zero

# exact.functional.test <- ExactFunctionalTest

fun.chisq.test <- function (x, method="fchisq", alternative="non-constant", log.p=FALSE,
                            index.kind="unconditional")
{
  if(!is.matrix(x) && !is.data.frame(x)) stop("input x must be matrix or data frame\n")

  row.chisq.sum <- sum(apply(x, 1,
                             function(v){
                               row.sum <- sum(v)
                               expected <- row.sum / length(v)
                               if(row.sum>0) sum( (v - expected)^2 / expected)
                               else 0
                             }))

  col.sum <- apply(x, 2, sum)
  n <- sum(col.sum)

  if(alternative == "non-constant") { # non-constant functional chi-square
    col.expected <- n / ncol(x)
    col.chisq <- sum((col.sum - col.expected)^2 / col.expected)
    fun.chisq <- row.chisq.sum - col.chisq
    df <- nrow(x) * (ncol(x) - 1) - (ncol(x) - 1)

    if(index.kind == "conditional") { # functional index given the column marginal
      max.fun.chisq <- n * (ncol(x) - 1) - col.chisq
      estimate.label <- "conditional function index xi.f"
    } else if(index.kind == "unconditional") { # unconditional functional index
      max.fun.chisq <- n * (ncol(x) - 1)
      estimate.label <- "function index xi.f"
    } else {
      stop("unrecognized function index kind ", index.kind)
    }

  } else if(alternative == "all") { # functional chi-square
    fun.chisq <- row.chisq.sum
    df <- nrow(x) * (ncol(x) - 1)
    max.fun.chisq <- n * (ncol(x) - 1)
    estimate.label <- "function index xi.f"
  } else {
    stop("unknown alternative ", alternative)
  }

  if(max.fun.chisq > 0) {
    estimate <- sqrt(abs(fun.chisq) / max.fun.chisq)
  } else {
    estimate <- 0
  }

  names(estimate) <- estimate.label

  DNAME <- deparse(substitute(x))

  if(method=="default" || method=="fchisq") {
    method.text <- "Functional chi-square test"
    names(fun.chisq) <- "statistic"
    names(df) <- "parameter"
    p.value <- pchisq( fun.chisq, df = df, lower.tail=FALSE, log.p=log.p )
    return(structure(list( statistic=fun.chisq, parameter=df, p.value=p.value,
                           estimate = estimate, data.name=DNAME,
                           method = method.text),
                     class = "htest"))

  } else if(method=="normalized" || method=="nfchisq") {
    method.text <- "Normalized functional chi-square test"
    if(df > 0) {
      normalized <- as.numeric((fun.chisq-df)/sqrt(2*df))
    } else {
      normalized <- -Inf
    }
    p.value <- pnorm( normalized, lower.tail=FALSE, log.p=log.p )
    names(normalized) <- "statistic"
    names(df) <- "parameter"
    return(structure(list(statistic = normalized, parameter = df, p.value = p.value,
                          estimate = estimate, data.name = DNAME, method = method.text),
                     class = "htest"))

  } else if(method=="exact") {

    method.text <- "Exact functional test"
    if(sum(x%%1!=0)>=1) { # Check whether numbers in x are all integers
      stop("ERROR: Exact test requires integer contingency tables!", call. = TRUE)
    }
    ####

    ####
    #Hua added, Nov 13, 2014
    #Exact functional test
    if((sum(x) <= 200 || sum(x)/nrow(x)/ncol(x) <=5)
       && nrow(x)<=5 && ncol(x)<=5) {
      # p.value <- exact.functional.test(x)
      p.value <- ExactFunctionalTest(x)
      if(log.p) p.value <- log(p.value)
      names(fun.chisq) <- "statistic"
      return(structure(list(statistic = fun.chisq, p.value = p.value, estimate = estimate,
                            data.name = DNAME, method = method.text),
                       class = "htest"))
    } else {
      return(fun.chisq.test(x, method="fchisq", alternative=alternative, log.p=log.p,
                            index.kind=index.kind))
    }

    ####
  } else {
    stop("ERROR: unrecognized method argument", method)
  }
}

cp.fun.chisq.test <- function(x, method="fchisq", log.p=FALSE)
{
  if(mode(x)!="list" || length(x)<2 )
  {
    stop("only accept list of 2 or more matrices as input!")
  }

  finalStat <- 0
  finalDf <- 0

  for(i in 1:nrow(x[[1]]))
  {
    oneT <- c() # one table for each row
    for(j in 1:length(x))
    {
      oneT <- rbind(oneT, x[[j]][i,])
    }
    oneresult <- fun.chisq.test(oneT)
    finalStat <- finalStat + oneresult$statistic
    finalDf <- finalDf + oneresult$parameter
  }

  DNAME <- deparse(substitute(x))

  if(method=="default" || method=="fchisq") {

    names(finalStat) <- "statistic"
    names(finalDf) <- "parameter"
    p.value <- pchisq(finalStat, df = finalDf, lower.tail=FALSE, log.p=log.p)
    return(structure(list( statistic=finalStat, parameter=finalDf, p.value=p.value,
                           method = "Comparative functional chi-square test for heterogeneity",
                           data.name= DNAME),
                     class = "htest"))

  } else if(method=="normalized" || method=="nfchisq") {

      finalStat <- as.numeric((finalStat-finalDf)/sqrt(2*finalDf))
      names(finalStat) <- "statistic"
      names(finalDf) <- "parameter"
      return(structure(list(statistic = finalStat, parameter = finalDf,
                            p.value = pnorm( finalStat, lower.tail=FALSE, log.p=log.p),
                            method = "Nomalized comparative functional chi-square test for heterogeneity",
                            data.name= DNAME),
                       class = "htest"))

  } else {
    stop("method can only be \"default\", \"normalized\", or \"exact\".\n")
  }
}

# exact.functional.test <- function(x){
#  res <- .Call("ExactFunctionalTest", x, PACKAGE="FunChisq")
#  return (as.double(res))
#}
####
