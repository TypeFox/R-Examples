NormalAndTplot.htest <- function(mean0, type="hypothesis", xlim=NULL, mean1=NA, ...,
                                 xbar, sd, df, n, alpha.left, alpha.right,  ## ignored
                                 distribution.name, sub ## these input arguments will be ignored
                               ) {
  if (names(mean0$statistic) == "X-squared")
    stop("Not yet available for proportion tests", call.=FALSE)

  ## this function takes "htest" objects from base::t.test and TeachingDemos::z.test
  t.htest <- mean0    ## store "htest" object as t.htest
  mean0 <- t.htest$null.value  ## use mean0 for the null hypothesis mean
  tstat <- t.htest$statistic
  distribution.name <- names(tstat)
  if (!distribution.name %in% c("t","z")) stop("t or z tests only", call.=FALSE)
  estimate <- t.htest$estimate
  xbar <- switch(length(estimate),
                 estimate,
                 estimate[1] - estimate[2])
  parameter <- t.htest$parameter
  df <- parameter["df"]
  if (is.null(df) || is.na(df)) df <- Inf
  n <- parameter["n"]
  if (is.null(n) || is.na(n)) n <- 1
  stderr <- (xbar-mean0) / tstat
  ## if (n > 1) sd <- stderr*sqrt(n)
  switch(t.htest$alternative,
         two.sided={
           alpha.left <- (1-attr(t.htest$conf.int, "conf.level"))/2
           alpha.right <- alpha.left
           sided <- "both"
           if (is.null(xlim))
             xlim <- if (type=="hypothesis")
                       range(mean0, mean1, xbar, na.rm=TRUE) + c(-3,3)*stderr
                       ## 1.12*c(-1,1)*max(abs(c(xbar, diff(t.htest$conf.int))))
                     else
                       xbar + .8*c(-1,1)*diff(t.htest$conf.int)
         },
         less={
           if (type=="hypothesis") {
             alpha.left <- 1-attr(t.htest$conf.int, "conf.level")
             alpha.right <- 0
             sided <- "left"
             if (is.null(xlim))
               xlim <- range(mean0, mean1, xbar, na.rm=TRUE) + c(-3,3)*stderr
               ## xlim <- c(min(xbar, mean0- 3.5*stderr),
               ##           max(t.htest$conf.int[2], mean0 + 3.5*stderr, xbar))
           }
           else { ## type == "confidence"
             alpha.left <- 0
             alpha.right <- 1-attr(t.htest$conf.int, "conf.level")
             sided <- "right"
             if (is.null(xlim))
               xlim <- xbar + c(-1,1)*3.5*stderr
           }
         },
         greater={
           if (type=="hypothesis") {
             alpha.left <- 0
             alpha.right <- 1-attr(t.htest$conf.int, "conf.level")
             sided <- "right"
             if (is.null(xlim))
               xlim <- range(mean0, mean1, xbar, na.rm=TRUE) + c(-3,3)*stderr
               ## xlim <- c(min(t.htest$conf.int[1], mean0 - 3.5*stderr, xbar),
               ##           max(mean0 + 3.5*stderr, xbar))
           }
           else { ## type == "confidence"
             alpha.left <- 1-attr(t.htest$conf.int, "conf.level")
             alpha.right <- 0
             sided <- "left"
             if (is.null(xlim))
               xlim <- xbar + c(-1,1)*3.5*stderr
           }
         }
         )

  sub <- paste(t.htest$method, t.htest$data.name, sep=": ")
  number.vars <- switch(t.htest$method,
                        "One Sample z-test"=1,
                        "One Sample t-test"=1,
                        "Paired t-test"=2,
                        "Two Sample t-test"=2,
                        " Two Sample t-test"=2,
                         "Welch Two Sample t-test"=2)

  NormalAndTplot(mean0=mean0, xbar=xbar, mean1=mean1, sd=stderr*sqrt(n), df=df, n=n,
                 alpha.left=alpha.left, alpha.right=alpha.right,
                 xlim=xlim, distribution.name=distribution.name,
                 sub=sub, type=type, number.vars=number.vars, ...)
}
