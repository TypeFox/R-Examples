## cp.calc takes
##    the output from the stepwise() function,
##    the data matrix,
##    the name of the response variable
## and constructs a matrix including the Cp statistic.
## The Cp plot uses several columns from this matrix.

## usair.step <- stepwise(y=usair$lnSO2,
##                        x=usair[, c(2,9,10,5:7)],
##                        method="exhaustive",
##                        plot=FALSE, nbest=2)
## usair.cp <- cp.calc(usair.step, usair, "lnSO2")
## plot(cp ~ p, data=usair.cp, type="n")
## abline(a=0, b=1)
## text(x=usair.cp$p, y=usair.cp$cp, row.names(usair.cp))
##
## tmp <- (usair.cp$cp <= 10)
## plot(cp ~ p, data=usair.cp[tmp,], ylim=c(0,10), type="n")
## abline(a=0, b=1)
## text(x=usair.cp$p[tmp], y=usair.cp$cp[tmp], row.names(usair.cp)[tmp])

cp.calc <- function(sw, data, y.name) {
  tss <-
    ## if.R(s=var(data[[y.name]], SumSquares=TRUE),
    ## r=
    var(data[[y.name]]) * (length(data[[y.name]])-1)
    ##)
  r2 <- (tss-sw$rss)/tss
  n <- nrow(data)
  num.x <- max(sw$size)
  full.i <- match(TRUE, num.x == sw$size)
  names(full.i) <- dimnames(sw$which)[[1]][full.i]
  rss.full <- sw$rss[full.i]
  p <- sw$size+1
  r2.adj <- 1 - ((n-1)/(n-p))*(1-r2)
  msres.full <- rss.full/(n-(1+ncol(sw$which)))
  cp <- sw$rss/msres.full + 2*p - n
  aic <- msres.full * (cp + n)
  result <- data.frame(p=sw$size+1, cp=cp, aic=aic, rss=sw$rss,
                       r2=r2, r2.adj=r2.adj)
  result$xvars <-
    apply(sw$which, 1, function(x,n) paste(n[x], collapse=","),
          n=dimnames(sw$which)[[2]])
##   result$xvars.abbrev0 <-
##     apply(sw$which, 1,
##           function(x,n) paste(substring(n,1,1)[x],collapse=""),
##           n=dimnames(sw$which)[[2]])
  row.names(result) <-
    apply(sw$which, 1,
          function(x,n) {
            first.letter <- substring(n,1,1)
            paste(ifelse(x,first.letter,"."), collapse="")
          },
          n=dimnames(sw$which)[[2]])
  result$sw.names <- dimnames(sw$which)[[1]]
  unique.row.names <- (!duplicated(result$sw.names))  &
                      (substring(result$sw.names, 1, 1) != "0")
  result <- result[unique.row.names,]
  attr(result,"full.i") <- row.names(result)[full.i]
  attr(result,"y.name") <- y.name
  attr(result,"n") <- n
  attr(result,"tss") <- tss
  oldClass(result) <- c("cp.object","data.frame")
  result
}
## trace(cp.calc, exit=browser)

print.cp.object <- function(x, ...) {
  cat("response variable = ", attr(x,"y.name"), "\n",
      "total sum of squares = ", signif(attr(x,"tss"),7), "\n",
      "number of observations = ", attr(x,"n"), "\n",
      "full model is row ", attr(x,"full.i"), "\n", sep="")
  xx <- x
  for (a in c("full.i", "y.name", "n", "tss" ))
    attributes(x)[[a]] <- NULL
  NextMethod("print", digits=4)
  invisible(xx)
}


"[.cp.object" <- function(x, ..., drop = TRUE)
  "[.data.frame"(x, ..., drop=TRUE)
