`lmodel2` <-
    function(formula, data = NULL, range.y = NULL, range.x = NULL,
             nperm = 0)
###
### Bivariate model II regression.
### Regression methods: OLS, MA, SMA, RMA
###
### formula: A formula specifying the model, as in 'lm' and 'aov'.
###
### data   : A data frame containing the variables specified in the formula. 
### 
### Ranged major axis: range.y = "relative" : variable y has a true
###    zero (relative-scale variable), range.y = "interval" : variable
###    y possibly includes negative values (interval-scale variable)
###    range.x = "relative" : variable x has a true zero
###    (relative-scale variable) range.x = "interval" : variable x
###    possibly includes negative values (interval-scale variable)
###
###          Pierre Legendre, December 2007 - June 2008
{

    ## From the formula, find the variables and the number of observations 'n'
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    ## var.names = colnames(mf)
    y <- mf[,1]
    x <- mf[,2]
    n <- nrow(mf)

    ## Check the data
    RMA <- FALSE
    if( (length(range.y) != 0) && (length(range.x) != 0) ) {
        RMA <- TRUE
    } else {
        if(length(range.y) != 0)
            stop("Give a value (relative or interval) to parameter 'range.x' for ranging") 
        if(length(range.x) != 0)
            stop("Give a value (relative or interval) to parameter 'range.y' for ranging") 
   }
    if(!RMA)
        message("RMA was not requested: it will not be computed.",'\n')
    if(nperm < 0)
        stop("'nperm' cannot be negative")

    if(length(x) != n)
        stop("Vectors y and x are not of the same length")
    if(var(y) == 0)
        stop("Variance(y) = 0. The regression coefficients cannot be computed")
    if(var(x) == 0)
        stop("Variance(x) = 0. The regression coefficients cannot be computed")

    ## Common calculations
    ybar <- mean(y)
    xbar <- mean(x)
    yx <- cbind(y,x)
    cov.mat <- cov(yx)
    vary <- cov.mat[1,1]
    varx <- cov.mat[2,2]
    r <- cor(y,x)
    rsquare <- r^2

    info.slope <- 0
    info.CI <- 0
    ## 2-tailed t-value, alpha = 0.05, to compute 95% C.I.
    t <- qt(0.025, (n-2), lower.tail=FALSE)
    epsilon <- .Machine$double.eps          

    ## Ordinary least squares regression, OLS
    met <- "OLS"
    ## lm {stats}  Fitting linear model
    temp <- lm(y ~ x)          
    b.ols <- summary(temp)$coefficients[2,1]
    angle <- atan(b.ols)*180/pi
    P.param <- summary(temp)$coefficients[2,4]
    ## confint {stats}  Confidence intervals of model parameters
    temp2 <- confint(temp)     

    res1 <- summary(temp)$coefficients[1,1]
    res2 <- b.ols
    res3 <- angle
    res4 <- temp2[1,1]
    res5 <- temp2[1,2]
    res6 <- temp2[2,1]
    res7 <- temp2[2,2]
    res8 <- NA

    ## Major axis regression, MA
    met <- c(met,"MA")

    MA.res <- MA.reg(b.ols, r, rsquare, ybar, xbar, epsilon, 1)
    b.ma <- MA.res[2]
    if(rsquare <= epsilon) {
        info.slope <- 1
        warning("MA:  R-square = 0",'\n')
    }
    CL <- CLma(b.ma, n, r, rsquare, xbar, ybar, cov.mat, t, 1, epsilon)
    lambda1 <- CL[5]
    lambda2 <- CL[6]
    H <- CL[7]
    A <- CL[8]
    if(CL[9] == 1)
        info.CI <- 1

    res1 <- c(res1, MA.res[1])
    res2 <- c(res2, MA.res[2])
    res3 <- c(res3, MA.res[3])
    res4 <- c(res4, CL[1])
    res5 <- c(res5, CL[2])
    res6 <- c(res6, CL[3])
    res7 <- c(res7, CL[4])
    res8 <- c(res8, NA)

    ## Standard major axis regression, SMA
    met <- c(met,"SMA")

    SMA.res <- SMA.reg(b.ols, r, rsquare, ybar, xbar, epsilon)
    b.sma <- SMA.res[2]
    if(rsquare <= epsilon) {
       info.slope <- 1
       warning("SMA: R-square = 0",'\n')
    }
    CL <- CLsma(b.sma, n, r, rsquare, xbar, ybar, t, epsilon)
    res1 <- c(res1, SMA.res[1])
    res2 <- c(res2, SMA.res[2])
    res3 <- c(res3, SMA.res[3])
    res4 <- c(res4, CL[1])
    res5 <- c(res5, CL[2])
    res6 <- c(res6, CL[3])
    res7 <- c(res7, CL[4])
    res8 <- c(res8, NA)

    ## Ranged major axis regression, RMA
    yx.2 = yx
    if(RMA) {
       met = c(met,"RMA")
       
       ## Range y and x
       range.yx = apply(yx,2,range)
       if(range.y == "relative") {
           if(range.yx[1,1] < 0 )
               stop("Negative values in 'y'. Use 'interval' for ranging")
           range.yx[1,1] <- 0
       }
       if(range.x == "relative") {
           if(range.yx[1,2] < 0)
               stop("Negative values in 'x'. Use 'interval' for ranging")
           range.yx[1,2] <- 0
       }
       ## Apply ranging to y and x. The table with ranged values is 'yx.2'
       yx.1 <- sweep(yx, 2, range.yx[1,], "-")
       yx.2 <- sweep(yx.1, 2, (range.yx[2,] - range.yx[1,]), "/")
       ran.y <-  (range.yx[2,1] - range.yx[1,1])
       ran.x <- (range.yx[2,2] - range.yx[1,2])
       ratio <- ran.y / ran.x
       ## cat("ran.y =",ran.y,"ran.x =",ran.x,"ratio =",ratio,'\n')   # Debug

       ## Compute RMA regression

       ## Note: cor(yx.2[,1],yx.2[,2]) = cor(y,x); the r-squares are
       ## thus also equal
       temp.ranged <- lm(yx.2[,1] ~ yx.2[,2])
       b.ols.ranged <- summary(temp.ranged)$coefficients[2,1]
       RMA.res <- MA.reg(b.ols.ranged, r, rsquare, ybar, xbar, epsilon, ratio)
       b.rma <- RMA.res[2]
       if(rsquare <= epsilon) {
           info.slope <- 1
           warning("RMA: R-square = 0")
       }
       cov.rma <- cov(yx.2)
       CL <- CLma(b.rma, n, r, rsquare, xbar, ybar, cov.rma, t, ratio, epsilon)
       if(CL[9] == 1)
           info.CI <- 1

       res1 <- c(res1, RMA.res[1])
       res2 <- c(res2, RMA.res[2])
       res3 <- c(res3, RMA.res[3])
       res4 <- c(res4, CL[1])
       res5 <- c(res5, CL[2])
       res6 <- c(res6, CL[3])
       res7 <- c(res7, CL[4])
       res8 <- c(res8, NA)
   }

    ## Angle (degrees) between the two OLS regression lines (Numerical
    ## ecology 1998, eq. 10.5)
    sdy <- sqrt(vary)
    sdx <- sqrt(varx)
    theta <- 90 - sign(r)*( atan(r*sdx/sdy)*180/pi + atan(r*sdy/sdx)*180/pi )
    
    ## Informative messages
    if(nperm == 0)
        message("No permutation test will be performed")
    
    ## Permutation tests
    if((nperm > 0) & (rsquare > epsilon)) {
        ## requires vegan if permuted.index2 used in permutest
        ## require(vegan) || stop("requires package 'vegan'")
        res8 <- permutest.lmodel2(yx, yx.2, b.ols, b.ma, b.rma/ratio,
                          RMA, ratio, nperm, epsilon)
    }
    
    ## Output results
    if(H == 0)
        H <- NA
    reg.res <- data.frame(met,res1,res2,res3,res8)
    CI.res <- data.frame(met,res4,res5,res6,res7)
    colnames(reg.res) <- c("Method","Intercept","Slope",
                           "Angle (degrees)","P-perm (1-tailed)")
    colnames(CI.res) <- c("Method","2.5%-Intercept","97.5%-Intercept",
                          "2.5%-Slope","97.5%-Slope")

    out <- list(y=y, x=x, regression.results=reg.res,
                confidence.intervals=CI.res,
                eigenvalues=c(lambda1,lambda2), H=H, n=n, r=r,
                rsquare=rsquare, P.param=P.param, theta=theta,
                nperm=nperm, epsilon=epsilon, info.slope=info.slope,
                info.CI=info.CI, call = match.call())

    class(out) <- "lmodel2"
    out
}


`plot.lmodel2` <-
    function(x, method="MA", confidence = TRUE, centroid = FALSE,
             col = c("red", "gray60", "black"), main, ...)

    ## Specify as follows the regression result to be plotted:
    ##   method="OLS" for ordinary least-squares regression
    ##   method="MA"  for major axis regression
    ##   method="SMA" for standard major axis regression
    ##   method="RMA" for ranged major axis regression

{
    out <- x
    y <- out$y
    x <- out$x
    centr.y <- mean(y)
    centr.x <- mean(x)
    row <- which(out$regression.results == method)
    a <- out$regression.results[row,2]
    b <- out$regression.results[row,3]
    b1 <- out$confidence.intervals[row,4]
    a1 <- centr.y - b1*centr.x
    b2 <- out$confidence.intervals[row,5]
    a2 <- centr.y - b2*centr.x
    plot(x, y,  ...)
    if (centroid)
        points(centr.x, centr.y, pch=3, cex=2)
    if((row != 1) && (out$rsquare <= out$epsilon)) {
        warning("R-square = 0: model and C.I. not drawn for MA, SMA or RMA")
    } else {
        abline(a, b, col=col[1])
        if(confidence && !( is.na(a1) )) {
            abline(a1, b1, col = col[2])
            abline(a2, b2, col = col[2])
        }
    }
    if (missing(main))
        title(main = paste(method,"regression"))
    else
        title(main = main)
    ##     cat('\n','a =',a,'  b =',b,'\n')
    ##     cat('\n','a1 =',a1,'  b2 =',b2,'\n')
    ##     cat('\n','a2 =',a2,'  b1 =',b1,'\n')
    invisible(x)
}

`lines.lmodel2` <-
    function(x, method="MA", confidence = TRUE, col = c("red", "gray60"), ...)
## Similar to plot, but only adds line, optionally with CI
{
    out <- x
    y <- out$y
    x <- out$x
    centr.y <- mean(y)
    centr.x <- mean(x)
    row <- which(out$regression.results == method)
    a <- out$regression.results[row,2]
    b <- out$regression.results[row,3]
    b1 <- out$confidence.intervals[row,4]
    a1 <- centr.y - b1*centr.x
    b2 <- out$confidence.intervals[row,5]
    a2 <- centr.y - b2*centr.x
    if((row != 1) && (out$rsquare <= out$epsilon)) {
        warning("R-square = 0: model and C.I. not drawn for MA, SMA or RMA")
    } else {
        abline(a, b, col=col[1])
        if(confidence && !(is.na(a1))) {
            abline(a1, b1, col = col[2])
            abline(a2, b2, col = col[2])
        }
    }
    invisible(x)
}

`print.lmodel2` <-
    function(x, ...)
{
## Print the regression results
    cat("\nModel II regression\n\n")
    writeLines(strwrap(paste("Call:",
                             paste(deparse(x$call), collapse=" "),
                             "\n")))
    cat("\n")
    cat("n =",x$n,"  r =",x$r,"  r-square =",x$rsquare,'\n')
    cat("Parametric P-values:   2-tailed =",x$P.param,
        "   1-tailed =",x$P.param/2,'\n')
    cat("Angle between the two OLS regression lines = ",
        x$theta," degrees",sep="",'\n')
    if(x$info.slope == 1)
        cat("MA, SMA, RMA slopes = NA when the correlation is zero\n")
    if(x$info.CI == 1)
        cat("Confidence limit of slope = Inf when infinite (90 deg. angle)\n")   
    if(x$nperm > 0) {
        cat("\nPermutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign\n")
        cat("A permutation test of r is equivalent to a permutation test of the OLS slope\n")
        cat("P-perm for SMA = NA because the SMA slope cannot be tested\n")
    }
    if(is.na(x$confidence.intervals[2,2]) ||
       is.na(x$confidence.intervals[3,2])) {
        cat("\nConfidence interval = NA when the limits of the confidence interval\n")
        cat("cannot be computed. This happens when the correlation is 0\n")
        cat("or the C.I. includes all 360 deg. of the plane (H >= 1)\n")
    }
    if((x$nperm > 0) && (x$rsquare < x$epsilon))
        cat("\nR-square = 0: no permutation test was performed\n")
    cat("\nRegression results\n")
    print(x$regression.results)
    cat("\nConfidence intervals\n")
    print(x$confidence.intervals)
    cat("\nEigenvalues: ")
    cat(x$eigenvalues,"\n")
    cat("\nH statistic used for computing C.I. of MA: ")
    cat(x$H,"\n")
    cat("\n")
    invisible(x) 
}

`MA.reg` <-
    function(b.ols, r, rsquare, ybar, xbar, epsilon, ratio)
{
    if(rsquare > epsilon) {
        d <- (b.ols^2 - rsquare) /(rsquare * b.ols)
   b.ma <- 0.5*(d + sign(r)*sqrt(d^2+4))
   b.ma <- b.ma*ratio
   b0 <- ybar - b.ma*xbar
   angle <- atan(b.ma)*180/pi
    } else { 
        b0 <- NA
        b.ma <- NA
        angle <- NA
    }
    MA.res <- c(b0, b.ma, angle)
}

`SMA.reg` <-
    function(b.ols, r, rsquare, ybar, xbar, epsilon)
{
    if(rsquare > epsilon) { 
        b.sma <- sign(r) * b.ols/r 
        b0 <- ybar - b.sma*xbar
        angle <- atan(b.sma)*180/pi
    } else { 
        b0 <- NA
        b.sma <- NA
        angle <- NA
    }
    SMA.res <- c(b0, b.sma, angle)
}

`CLma` <-
    function(b.ma, n, r, rsquare, xbar, ybar, cov.mat, t, ratio, epsilon)
## MA confidence intervals for slope and intercept.
## Note: the confidence interval of the intercept is underestimated.
{
    H <- 0
    A <- NA
    info.CI <- 0
    ## Compute eigenvalues by eigen(covariance matrix), as in PCA
    mat.eig <- eigen(cov.mat)
    lambda1 <- mat.eig$values[1]
    lambda2 <- mat.eig$values[2]
#
    if(rsquare > epsilon) {
        b.ma <- b.ma/ratio
        if(lambda2 <= epsilon) {               # If eigenvalue2 = 0
            b1inf <-  b.ma
            b1sup <- b.ma
        } else {
            if( (lambda1-lambda2) < epsilon) {  # If equal eigenvalues
                tmp <- atan(b.ma)*180/pi
                b1inf <- tan((tmp-90)*pi/180)
                b1sup <- tan((tmp+90)*pi/180)
            } else {
                H <- t^2 / (((lambda1/lambda2)+(lambda2/lambda1)-2)*(n-2))
                if(H >= 1) {                     # If H >= 1
                    b1inf <- NA
                    b1sup <- NA
                } else {
                    A <- sqrt(H/(1-H))
                    if((b.ma*A) == -1) {          # angle = -90 deg., b1inf = -infinity
                        b1inf <- -Inf
                        info.CI <- 1
                    } else {
                        b1inf <- ratio * (b.ma-A) / (1+b.ma*A)
                    }
                    if((b.ma*A) == 1) {           # angle = +90 deg., b1sup = +infinity
                        b1sup <- Inf
                        info.CI <- 1
                    } else {
                        b1sup <- ratio * (b.ma+A) / (1-b.ma*A)
                    }
                }
            }
        }
        if((H == 0) || (H >= 1)) {             # H >= 1
            b0inf <- NA
            b0sup <- NA
        } else {
            if(xbar >= 0) {
                b0inf <- ybar - b1sup*xbar
                b0sup <- ybar - b1inf*xbar
            } else {
                b0inf <- ybar - b1inf*xbar
                b0sup <- ybar - b1sup*xbar
            }
        }
    } else {
        b1inf <- NA
        b1sup <- NA
        b0inf <- NA
        b0sup <- NA
    }
    CL <- c(b0inf, b0sup, b1inf, b1sup, lambda1, lambda2, H, A, info.CI)
}

`CLsma` <-
    function(b.sma, n, r, rsquare, xbar, ybar, t, epsilon)
## SMA confidence intervals for slope and intercept.
## C.I. of slope following Jolicoeur & Mosimann (1968), McArdle (1988).
## Note: the confidence interval of the intercept is underestimated.
{
    if(rsquare > epsilon) {
        B <- t^2 * (1-rsquare) / (n-2)
        if(b.sma > 0) {
            b1inf <- b.sma * (sqrt(B+1) - sqrt(B))
            b1sup <- b.sma * (sqrt(B+1) + sqrt(B))
        } else {
            b1inf <- b.sma * (sqrt(B+1) + sqrt(B))
            b1sup <- b.sma * (sqrt(B+1) - sqrt(B))
        }
        if(xbar >= 0) {
            b0inf <- ybar - b1sup*xbar
            b0sup <- ybar - b1inf*xbar
        } else {
            b0inf <- ybar - b1inf*xbar
            b0sup <- ybar - b1sup*xbar
        }      
    } else {
        b1inf <- NA
        b1sup <- NA
        b0inf <- NA
        b0sup <- NA
    }
    CL <- c(b0inf, b0sup, b1inf, b1sup)
}

## vegan defines generic permutest, but if we don't depend on vegan we
## need to use another name
`permutest.lmodel2` <-
    function(x, yx.2, b.ols, b.ma, b.rma, RMA, ratio, nperm, epsilon, ...)

### One-tailed permutation tests of OLS, MA, and RMA regression slopes.
### The SMA slope cannot be tested for significance.
### A permutation test of r is equivalent to a permutation test of the OLS slope.
### When abs(b.ma) > 1, the test is done on 1/b.ma ; likewise for b.rma.
{
    yx <- x
    nGE.ols <- 1
    ols.pos <- TRUE
    if(b.ols < 0)
        ols.pos <- FALSE
    y <- yx[,1]
    x <- yx[,2]

    nGE.ma <- 1
    ma.pos <- TRUE
    if(b.ma < 0)
        ma.pos <- FALSE
    isur.ma <- FALSE
    if(abs(b.ma) > 1) {
        isur.ma <- TRUE
        b.ma <- 1/b.ma
    }

    if(RMA) {
        nGE.rma <- 1
        rma.pos <- TRUE
        if(b.rma < 0)
            rma.pos <- FALSE
        isur.rma <- FALSE
        if(abs(b.rma) > 1) {
            isur.rma <- TRUE
            b.rma <- 1/b.rma
        }
        y.2 <- yx.2[,1]
        x.2 <- yx.2[,2]
    }

    ## Permutation test begins
    n <- length(y)
    for(i in 1:nperm) {
        ## Permutation, could use permuted.index2
        idx <- sample(n)
        y.per <- y[idx]
        ## OLS regression
        temp <- lm(y.per ~ x)                # lm {stats}  Fitting linear model
        b.ols.per <- summary(temp)$coefficients[2,1]
        if(ols.pos) {
            if(b.ols.per >= b.ols)
                nGE.ols <- nGE.ols+1
        } else {
            if(b.ols.per <= b.ols)
                nGE.ols <- nGE.ols+1
        }

        ## MA regression
        r.per <- cor(y.per,x)
        rsq.per <- r.per^2
        temp <- (rsq.per * b.ols.per)
        if(abs(temp) > epsilon) {
            d <- (b.ols.per^2 - rsq.per) / temp
            b.ma.per <- 0.5*(d + sign(r.per)*sqrt(d^2+4))
            if(isur.ma) b.ma.per <- 1/b.ma.per
            if(ma.pos) {
                if(b.ma.per >= b.ma)
                    nGE.ma <- nGE.ma+1
            } else {
                if(b.ma.per <= b.ma)
                    nGE.ma <- nGE.ma+1 }
        }

        ## RMA regression
        if(RMA) {
            ## Permutation: could use permuted.index2
            y.2.per <- y.2[idx]
            r.per <- cor(y.2.per,x.2)
            rsq.per <- r.per^2
            temp.ranged <- lm(y.2.per ~ x.2)  # lm {stats}  Fitting linear model
            b.ols.ranged <- summary(temp.ranged)$coefficients[2,1]
            temp <- (rsq.per * b.ols.ranged)
            if(abs(temp) > epsilon) {
                d <- (b.ols.ranged^2 - rsq.per) / temp
                b.rma.per <- 0.5*(d + sign(r.per)*sqrt(d^2+4))
                if(isur.rma)
                    b.rma.per <- 1/b.rma.per
                if(rma.pos) {
                    if(b.rma.per >= b.rma)
                        nGE.rma <- nGE.rma+1
                } else {
                    if(b.rma.per <= b.rma)
                        nGE.rma <- nGE.rma+1 }
            }
        }
    }   # End permutations
    P.ols <- nGE.ols/(nperm+1)
    P.ma <- nGE.ma /(nperm+1)
    if(RMA)
        P.rma <- nGE.rma/(nperm+1)

    if(RMA) {
        res8 <- c(P.ols, P.ma, NA, P.rma)
    } else {
        res8 <- c(P.ols, P.ma, NA)
    }
}
