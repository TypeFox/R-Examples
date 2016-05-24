library(gmp)

##
##' @title Test a unary (if unary=TRUE) or *binary* function
##' @param FUN a function, such as add.bigq() ...
##' @param x a list of "numbers"
##' @param out string determining output class; if "str", use characters, otherwise double
##' @return
##' @author Antoine Lucas (& Martin Maechler)
##' @examples test(as.bigq, 0)
test <- function(FUN, x, xlabs, out = "str", unary = FALSE)
{
  if(missing(xlabs))
    xlabs <- if(is.character(names(x))) names(x) else sapply(x, formatN)
  stopifnot(is.function(FUN), is.list(x),
	    (n <- length(x)) >= 1, length(xlabs) == n)
  if(out == "str") {
      sortie <- as.character
      res <- ""
      error <- "error"
  } else {
      sortie <- as.double
      res <- 0
      error <- NA
  }
  nr <- if(unary) 1 else n
## FIXME  res <- matrix(res, nr, n, dimnames = list(if(!unary) xlabs, xlabs))
## for now
  res <- matrix(res, nr, n, dimnames = list(NULL, xlabs))

  for(i in 1:nr)
    for(j in 1:n) {
      e <- if(unary) try(FUN(x[[j]]),        silent = TRUE) else
                     try(FUN(x[[i]],x[[j]]), silent = TRUE)
      if(inherits(e, "try-error"))
        e <- error
      else if(length(e) == 0)
        e <- numeric()

      ## ## now, for some functions also compute the corresponding numeric values
      ## d <- if(unary) try(FUN(as.numeric(x[[j]])),                    silent = TRUE) else
      ##                try(FUN(as.numeric(x[[i]]),as.numeric(x[[j]])), silent = TRUE)

      res[i,j] <- sortie(e)[1]
    }
  res ## for printing, the user may prefer as.data.frame(.)
}## end{test}


allfunctionid <- c("as.bigz","add.bigz","sub.bigz","mul.bigz",
		   "divq.bigz","div.bigz","mod.bigz","pow.bigz",
		   "inv.bigz", "gcd.bigz", "gcdex", "lcm.bigz",
		   "as.bigq",
		   "add.bigq","sub.bigq","div.bigq", "mul.bigq", "^.bigq",
		   "chooseZ",
		   "max.bigq","max.bigz","min.bigq","min.bigz")
unaryfunctionid <- c("log.bigz","log2.bigz","log10.bigz","c.bigz",
		     "isprime","nextprime", "factorialZ",
		     "sizeinbase","fibnum","fibnum2","lucnum","lucnum2",
		     "factorize","abs")
numericFunName <- function(gmpName) {
  if(gmpName != (r <- sub("[ZQ]$","", gmpName)) &&
     r!="as" && existsFunction(r)) # e.g. chooseZ
    return(r)
  if(gmpName != (r <- sub("\\.big[zq]$","", gmpName)) &&
     r!="as" && r!="sub" && existsFunction(r))
    return(r)
  ttt <- c("add" = "+",
           "sub" = "-",
           "mul" = "*",
           "pow" = "^",
           "div" = "/",
           "divq" = "%/%",
           "mod" = "%%")
  if(!is.na(t.r <- ttt[r]))
    t.r[[1L]]
  else ## return argument
    gmpName
}


options(width = 125)

sapply(allfunctionid,   numericFunName)
sapply(unaryfunctionid, numericFunName)


ex <- expression(23, "25", 2.3, -4, 4L, 0, as.bigz(34),
                 as.bigq(32,7), as.bigz(31,45), NULL,NA, -3L)## TODO:  as.bigz(3)^700
x <- lapply(ex, eval)
## Those "numbers" in x for which arithmetic should also work in double precision:
## not modulo-arithmetic, not larger than double.prec
useN <- sapply(x, function(.) is.null(.) || is.na(.) ||
               (is.finite(as.numeric(.)) && (!inherits(., "bigz") || is.null(modulus(.)))))
## names(x) <- sapply(ex, format)
## shorter & easier:
names(x) <- sapply(x, formatN)
str(x)
x. <- x[useN]
nx <- lapply(x., as.numeric)
gmp.NS <- asNamespace("gmp")# also get namespace *hidden* functions, i.e. methods:
for(fid in allfunctionid)
  {
    ##cat ("------------------------------------------\n", fid, "\n\n", sep="")
    cat ("------------------------------------------\n", fid," ", sep="")
    FUN <- get(fid, envir = gmp.NS, mode="function")
    rc   <- test(FUN, x )
    res  <- test(FUN, x. , out = "numeric")
    if((nfid <- numericFunName(fid)) != fid || existsFunction(nfid, where=baseenv())) {
      FUN <- get(nfid, envir = gmp.NS, mode="function")
      if(nfid != fid) cat("-> num.fn.:", nfid)
      cat("\n-> all.equal(target = res, current = F(<numeric x>)):\n")
      nres <- test(FUN, nx, out = "numeric")
      print(all.equal(res, nres))
    } else cat("\n\n")
    print(as.data.frame(rc)); cat("\n")
    ##    ^^^^^^^^^^^^^ (for now, to diminuish difference to last version )
  }

##==============================================================================

for(fid in unaryfunctionid)
  {
    cat ("------------------------------------------\n", fid, "\n\n", sep="")
    FUN <- get(fid, envir = gmp.NS, mode="function")
    print(as.data.frame(test(FUN, x, unary=TRUE)))
  }

##==============================================================================

###----------- matrix -----------------------------
x  <- matrix(1:6,3)
stopifnot(identical(as.bigz(x), matrix(as.bigz(as.vector(x)), 3)),
          dim(x) == 3:2,
          dim(x) == dim(ym   <- as.bigz(x, 6:1)),
          dim(x) == dim(ymr  <- as.bigz(x, 4:6)),
          dim(x) == dim(ymc  <- as.bigz(x, 4)),
          dim(x) == dim(ymq  <- as.bigq(x)),
          dim(x) == dim(y    <- as.bigq(x, 6:1))
          ,
          apply(ym,1,max) == 1:3,
          apply(ym,2,min) == c(1,1,0))

x %*% t(x)

ym %*% t(ym)
ym %*% t(ymr)
ymc %*% t(ymc)
ymq %*% t(ymq)
y %*% t(y)

dd <- dim(D  <- diag(1:4))
stopifnot(dd == dim(Dmq  <- as.bigq(D)),
          dd == dim(Dm   <- as.bigz(D,6:1)),
          dd == dim(Dmr  <- as.bigz(D,7)),
          dd == dim(Dmc  <- as.bigz(D,4)),
          TRUE)
solve(D)
solve(Dmq)
solve(Dmr)
try(solve(Dmc))# Error: argument has no inverse
try(solve(Dm)) # Error: System is singular

(D.D <- D %*% t(Dm))# now [>= Jan.2012] works too
stopifnot(identical(D.D, tcrossprod(D,Dm)))

##
## some specific tests

factorize("33162879029270137")

