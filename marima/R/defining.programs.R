##' @title define.model
##' 
##' @description Function to define multivariate arma model
##' (indicator form) for marima.
##'
##' @param kvar = dimension of time series
##' @param ar = autoregresssion definition. For example ar=c(1,2,12) will
##' generate autoregression at lags 1, 2 and 12.
##' @param ma = moving average definition. Works like ar. If ma=c(1,2) moving
##' average terms at lags 1 and 2 are defined.
##' @param rem.var = no. of variable(s) not to be considered in marima.
##' @param reg.var = no. of variable(s) that can only act as regression
##' variable(s) such as (typically) a socalled leading indicator.
##' @param indep = no. of variable(s) that are independent of the other
##' variables. indep=c(2,4) makes variables 2 and 4 independent of all
##' other variables. Variables 2 and 4 may influence other variables.
##' @param no.dep = sequence of pairs of variables. For example
##' no.dep=c(1,2,2,3) means that variable 2 is not allowed in model
##' for variable 1, and variable 3 is not allowed in model for variable 2.
##' @param ar.fill = sequence of triplets: c(dependent variable,independent
##' variable, lag). ar.fill=c(2,3,12): Insert ar-indicator for model for
##' dependent variable 2 and independent variable 3 at lag 12.
##' @param ar.rem = sequence of triplets c(dependent variable,independent
##' variable, lag). ar.rem=c(2,3,12): remove (if present) ar-indicator for model
##' for dependent variable 2 and independent variable 3 at lag 12.
##' @param ma.fill =sequence of triplets: c(dependent variable,independent
##' variable, lag). ma.fill=c(2,3,12): Insert ma-indicator for model for
##' dependent variable 2 and independent variable 3 at lag 12.
##' @param ma.rem = sequence of triplets c(dependent variable,independent
##' variable, lag). ma.rem=c(2,3,12): remove (if present) ma-indicator for model
##' for dependent variable 2 and independent variable 3 at lag 12.
##'
##' The various parameters may (in some cases) accomplish the same model
##' requirements. The routine define.model apply these input parameters
##' successively in the following order: 1) rem.var, 2) reg.var, 3) indep,
##' 4) no.dep, 5) ar.fill, 6) ar.rem, 7) ma.fil, 8) ma.rem
##'
##' The parameters ar.fill, ar.rem, ma.fill and ma.rem are applied last, and
##' in that order. They overwrite what previously has been defined.
##' @param print (!0/0) If !0 is used, the generated patterns of the
##' arma model and other informations are printed on the console.
##' If 0 is used, no printout of the arma patterns are given. 
##'
##' @return ar.pattern = a matrix polynomium with 1's and 0's defining
##' the autoregressive matrix polynomium to be fitted by marima (type=array
##' with dim=c(kvar,kvar,1+ar_order) (with leading unity matrix)).
##'
##' @return ma.pattern = a matrix polynomium with 1's and 0's defining
##' the moving average matrix polynomium to be fitted by marima (type=array
##' with dim=c(kvar,kvar,1+ma_order) (with leading unity matrix)).
##'
##' @examples
##' #
##' # Example 1: 3-variate arma model with ar-lags at 1 and 2, and an
##' # ma-term at lag 1. And var=3 is a regression variable (X-variable).
##' #
##'  Model1<-define.model(kvar=3,ar=c(1,2),ma=c(1),reg.var=3)
##'  short.form(Model1$ar.pattern)
##'  short.form(Model1$ma.pattern,leading=FALSE)
##' #
##' # The object Model1 contains the ar- and ma-pattern arrays as defined.
##' #
##' # Model1$ar.pattern and Model1$ma.pattern are used as input to
##' # marima in order to define the model to be estimated.
##' #
##' # Example 2: arma model with ar-lags at 1, 2 and 6, and var=3
##' #  regression variable (X-variable).
##' #
##'  Model2<-define.model(kvar=3,ar=c(1,2,6),ma=c(1),reg.var=3)
##' # Print the ar- and ma-polynomial patterns using
##'  short.form(Model2$ar.pattern,leading=FALSE)
##'  short.form(Model2$ma.pattern,leading=TRUE)
##' #
##' # Example 3: arma model with ar-lags at 1, 2 and 6, and reg.var=3
##' # (X-variable). ma-order=1. Finally (ar.fill=c(2,3,4) puts  a '1'
##' # for (dep-var=2,indep-var=3,ar-lag=4).
##' #
##' # If further modifications of the ar- or ma-patterns are needed, it
##' # can be accomplished before calling marima (Model3$ar.pattern and
##' # Model3$ma.pattern are arrays).
##' #
##'  Model3<-define.model(kvar=3,ar=c(1,2,6),ma=c(1),reg.var=3,ar.fill=c(2,3,4))
##'  short.form(Model3$ar.pattern)
##'  short.form(Model3$ma.pattern)
##' #
##'  Model4<-define.model(kvar=3,ar=c(1,2,6),ma=c(1),reg.var=3,
##'  ar.fill=c(2,3,4),indep=c(1))
##'  short.form(Model4$ar.pattern)
##'  short.form(Model4$ma.pattern,leading=FALSE)
##'
##' @export

define.model <- function(kvar = 1, ar = 0, ma = 0, rem.var = 0,
                         reg.var = 0, no.dep = NULL, print = 0,
                         ar.fill = NULL, ar.rem = NULL, ma.fill = NULL,
                         ma.rem = NULL, indep = NULL) {
    if (print != 0) {
        cat("Call = kvar,ar,ma,rem.var,reg.var,no.dep  \n",
            "kvar", kvar, "ar", ar, "ma", ma, "rem.var", rem.var,
              "reg.var", reg.var, "\n\n   no.dep=", no.dep, "\n")
    }

    ar <- sort(ar[ar >= 0])
    ma <- sort(ma[ma >= 0])
    L <- max(ar) + max(ma)
    M <- 1
    if (max(ar) >= 0) {
        M <- max(ar) + 1
    }

    if ((max(ar) + max(ma)) <= 0) {
        cat("No model specified in define.model: Returning. \n")
           empty.model <- list(ar.pattern=diag(kvar),ma.pattern=diag(kvar))
        return(empty.model)
    }

    ind.poly <- array(data = 0, dim = c(kvar, kvar, L))

    if (length(ar) > 0) {
        if (max(ar) > 0) {
            for (i in 1:length(ar)) {
                ind.poly[, , ar[i]] <- 1
            }
        }
    }

    if (max(ma) > 0) {
        if (length(ma) > 0) {
            for (i in 1:length(ma)) {
                ind.poly[, , (max(ar) + ma[i])] <- 1
            }
        }
    }

    if (print != 0) {
        cat("ind.poly=", ind.poly, "\n")
        cat("rem.var=", rem.var, "\n")
        cat("reg.var=", reg.var, "\n")
    }

    if (length(reg.var) > 0) {
        reg.var <- reg.var[reg.var > 0]
        reg.var <- reg.var[reg.var <= kvar]
        if (length(reg.var) > 0) {
            if (print != 0) {
                cat("ind.poly", ind.poly[reg.var, , ], "\n")
            }
            ind.poly[reg.var, , ] <- 0

            if (print != 0) {
                cat("ind.poly", ind.poly[reg.var, , ], "\n")
            }
            if (L > max(ar)) {
                M <- max(ar) + 1

                if (print != 0) {
                  cat("Model=(kvar,( ar,ma ),M,L,reg.var,rem.var",
                      "\n(dim(ind.poly)) \n", kvar, "(", c(ar), c(ma), ")",
                      M, L, reg.var, rem.var, "(", dim(ind.poly), ") \n")
                }
                ind.poly[, reg.var, c(M:L)] <- 0
            }
        }
    }

    if (length(rem.var) > 0) {
        rem.var <- rem.var[rem.var > 0 & rem.var <= kvar]
        if (length(rem.var) > 0) {
            ind.poly[rem.var, , ] <- 0
            ind.poly[, rem.var, ] <- 0
        }
    }

    for (i in 1:kvar) {
        if (i == 1) {
            indic <- c(ind.poly[i, , ])
        }
        if (i > 1) {
            indic <- c(indic, c(ind.poly[i, , ]))
        }
    }

    indic <- as.integer(indic)

    poly <- array(ind.poly, dim = c(kvar, kvar, (max(ar) + max(ma))))
    ar.poly <- array(diag(kvar), dim = (c(kvar, kvar, 1)))
    ma.poly <- ar.poly
    if (max(ar) > 0) {
        ar.poly <- check.one(poly[, , 1:max(ar)])
    }
    if (max(ma) > 0) {
        ma.poly <- check.one(poly[, , M:L])
    }

    if (!is.null(indep)) {
        indep <- indep[indep > 0 & indep <= kvar]
        L <- length(indep)
        if (L > 0) {
            for (j in 1:L) {
                J <- indep[j]
                for (i in 1:(dim(ar.poly)[3])) {
                  ar.poly[J, , i] <- 0
                  ar.poly[J, J, i] <- 1
                }
                for (i in 1:(dim(ma.poly)[3])) {
                  ma.poly[J, , i] <- 0
                  ma.poly[J, J, i] <- 1
                }
            }
        }
    }

    if (!is.null(no.dep)) {
        # cat('no.dep = ',no.dep,'\n')
        no.dep <- matrix(no.dep, nrow = 2)
        if (print != 0) {
            cat("no.dep=", no.dep, "\n")
        }
        L <- dim(no.dep)[2]
        # cat('L=',L,'\n')
        for (i in 1:L) {
            # cat('i=',i,'L=',L)
            da <- dim(ar.poly)[3]
            dm <- dim(ma.poly)[3]
            # cat('no.dep = ',no.dep,dim(no.dep),'da,ma=',da,dm,'L=',L)
            for (j in 1:L) {
                ar.poly[no.dep[1, j], no.dep[2, j], 1:da] <- 0
            }
            for (j in 1:L) {
                ma.poly[no.dep[1, j], no.dep[2, j], 1:dm] <- 0
            }
        }
    }

    if (!is.null(ar.fill)) {
        arfill <- matrix(ar.fill, nrow = 3)
        if (length(c(arfill)) != length(c(ar.fill))) {
            cat(ar.fill, "ar.fill error \n")
        }
        L <- dim(arfill)[2]
        for (j in 1:L) {
            ar.poly[arfill[1, j], arfill[2, j], (arfill[3, j] + 1)] <- 1
        }
    }

    if (!is.null(ar.rem)) {
        arrem <- matrix(ar.rem, nrow = 3)
        if (length(c(arrem)) != length(c(ar.rem))) {
            cat(ar.rem, "ar.rem error \n")
        }
        L <- dim(arrem)[2]
        for (j in 1:L) {
            ar.poly[arrem[1, j], arrem[2, j], (arrem[3, j] + 1)] <- 0
        }
    }

    if (!is.null(ma.fill)) {
        mafill <- matrix(ma.fill, nrow = 3)
        if (length(c(mafill)) != length(c(ma.fill))) {
            cat(ma.fill, "ma.fill error \n")
        }
        L <- dim(mafill)[2]
        for (j in 1:L) {
            ma.poly[mafill[1, j], mafill[2, j], (mafill[3, j] + 1)] <- 1
        }
    }

    if (!is.null(ma.rem)) {
        marem <- matrix(ma.rem, nrow = 3)
        if (length(c(marem)) != length(c(ma.rem))) {
            cat(ma.rem, "ma.rem error \n")
        }
        L <- dim(marem)[2]
        for (j in 1:L) {
            ma.poly[marem[1, j], marem[2, j], (marem[3, j] + 1)] <- 0
        }
    }

    results <- list(ar.pattern = ar.poly, ma.pattern = ma.poly)
    return(results)
}

##' @title short.form
##' 
##' @description Function to condensate (and/or) print out matrix
##' polynomium leaving out empty lag matrices and, if specified,
##' the leading (unity) matrix.
##'
##' @param poly matrix polynomium (0-1 array as construced by define.model,
##' for example, or array of reals as estimated by marima).
##' @param name character string used as header in output (default='lag').
##' @param leading TRUE/FALSE. If leading=FALSE the leading (unity matrix)
##' is to be left out/suppressed.
##' @param tail TRUE/FALSE. If TRUE and the ar/ma-model only consists
##' of coefficient matrice(s) where all coefficients (except the
##' leading unity matrix) are all zero a first order coefficient matrix
##' (being zero) is retained (in order to avoid a model containing only
##' the leading unity matrix).
##'
##' If tail=TRUE and the coefficients in the first coefficient matrix
##' (after the leading unity matrix) are all zero, the leading unity
##' matrix is always retained.
##' 
##' @param digits the number of digits retained by short.form (default=6).
##'
##' @examples
##' Model<-define.model(kvar=4,ar=c(1,2,4),ma=c(1),reg.var=4)
##' short.form(Model$ar.pattern)
##' short.form(Model$ma.pattern)
##' short.form(Model$ar.pattern,leading=FALSE)
##' short.form(Model$ar.pattern,leading=FALSE)
##' #
##' M<-define.model(kvar=4,ma=c(1))
##' short.form(M$ar.pattern)
##' short.form(M$ar.pattern,tail=TRUE)
##' short.form(M$ar.pattern,leading=FALSE,tail=TRUE)
##'
##' @export

short.form <- function(poly = NULL, name = "Lag=", leading = TRUE,
            tail = FALSE, digits = 6 )
{
    "[" <- function(x, ...) .Primitive("[")(x, ..., drop = FALSE)
    # cat('poly= \n') print(poly)
    poly <- check.one(poly)
    # ensure array has leading unity matrix
    kvar <- dim(poly)[1]
    LL   <- dim(poly)[3]
    
    out.names  <- paste("y" , c(1:kvar),        sep = "")
    inpx.names <- paste("x" , c(1:kvar),        sep = "")
    inpy.names <- paste("=y", c(1:kvar),        sep = "")
    inpy.names <- paste(inpx.names, inpy.names, sep = "")
    pol.names  <- paste(name, 0:LL     ,        sep = "")

    KK <- max(which(colSums(abs(poly), dims = 2) != 0))
    poly <- poly[, , 1:KK]
    # cat('KK=',KK,'\n')

    if (tail & KK == 1) {
        poly <- array(c(poly, 0 * poly), dim = c(dim(poly)[1:2], 2))
        dimnames(poly) <- list(out.names, inpy.names, pol.names[1:2])
        Poly <- round(poly, digits)
        return(Poly)
    }

    dimnames(poly) <- list(out.names, inpy.names, pol.names[1:KK])
    Lags <- sort(which(colSums(abs(poly), dims = 2) != 0))
    # cat('Lags=',Lags,'\n') print(poly)
    if (leading != TRUE & max(Lags) > 1) {
        Lags <- setdiff(Lags, 1)
    }

    Poly <- round(poly[, , Lags], digits)
    # cat('Poly=','\n') print(Poly) cat('Lags=',Lags,'\n')
    return(Poly)
}

##'@title define.dif
##' 
##' @description Function to generate and apply a differencing matrix polynomial
##' (autoregressive form) defined by a pattern.
##'
##' To be used before calling marima in order to difference the
##' timeseries before the marima analysis. The averages of the variables
##' in the time series are subtracted from the input series before
##' differencing.
##'
##' @param series     = kvar-variate timeseries (kvar by n matrix).
##' @param difference = 2 by L matrix defining L differencing operations.
##'
##' @return y.dif = the differenced timeseries (the complete part)
##' @return y.lost = the first observations lost because of differencing
##' @return dif.poly = differencing polynomial
##' array = c(kvar,kvar,...) holding the autoregressive representation
##' of the specified differencing
##' @return averages = the averages of the original series as they
##' were subtracted before differencing
##' @return dif.series = the differenced series (y.lost followed by y.dif)
##'
##' @examples
##'
##' # Generate Y=series with 4 variables for illustration:
##' set.seed(4711)
##' Y<-matrix(round(100*rnorm(40)+10),nrow=4)
##'
##' # Example 1: use of difference parameter: If
##' difference=c(2,1,2,1,3,12)
##' difference
##' # the variable 2 is differenced
##' # twice, and variable 3 is differenced once with lag=12.
##'
##' # Example 2:
##' poly <- define.dif(series=Y,difference=c(2,1,3,1,3,1))
##' poly
##' # Generates a (4-variate) polynomial differencing array (with a leading
##' # unity matrix corresponding to lag=0, and (in the example) differencing
##' # of variable 2 for lag 1 and variable 3 for lag 1 but twice. Afterwards
##' # the series Y is differenced accordingly. Results in poly$series and
##' # poly$dif.poly .
##'
##' # Example 3: Generation and application of multivariate differencing
##' # polynomial. Re-use the 4-variate time series and use the
##' # differencing polynomial (ar-form):
##' # var=1, dif=1, var=2, dif=6, and var=3 and 4, no differencing.
##' dif.y <-define.dif(Y,c(1,1, 2,6, 3,0, 4,0))
##' # Now dif.y contains the differenced series and the differencing
##' # polynomial. Print the generated polynomial in short form:
##' short.form(dif.y$dif.poly)
##' # Specifying no differencing (3,0 and 4,0) may be omitted:
##' dif.y <-define.dif(Y,c(1,1, 2,6))
##' dif.y
##'
##' # Example 4:
##' y<-matrix(round(rnorm(1200)*100+50),nrow=6)
##' library(marima)
##' difference<-c(3,2,4,0,5,0,6,7)
##' matrix(difference,nrow=2)
##' Y<-define.dif(y,difference=difference)
##' round(rowMeans(Y$dif.series),2)
##' round(Y$averages,2)
##'
##' @export

define.dif <- function(series = series, difference = NULL) {
    Names <- rownames(series)
    d <- dim(series)
    T <- 0
    if (d[1] > d[2]) {
        series <- t(series)
        T <- 1
        Names <- rownames(series)
    }
    d <- dim(series)
    kvar <- d[1]
    OMIT <- 0
    DM <- rep(0, kvar)
    if (!is.null(difference)) {
        D <- matrix(difference, nrow = 2)
        L <- dim(D)[2]
        Means <- rep(1, kvar)
        for (i in 1:L) {
            v <- D[1, i]
            if (v < 0 | v > kvar) {
                cat("difference = ", difference, "\n")
                cat("Wrong specification of differencing: var = ", v, "\n")
            }
            t <- D[2, i]
            # cat('Var,t =',v,t,D[1:2,i],' \n')
            if (v > 0 & v <= kvar & t > 0) {
                DM[v] <- 1
                OMIT = max(c(OMIT, t))
            }
            # cat('DM = ',DM,' \n')
        }
    }
    for (i in 1:kvar) {
        Means[i] <- mean(series[i, ])
        series[i, ] <- series[i, ] - Means[i] * DM[i]
    }

    dif.poly <- array(diag(kvar), dim = c(kvar, kvar, 1))
    if (!is.null(difference)) {
        # cat('kvar=',kvar,' \n')
        dif <- difference
        L <- length(c(dif))/2
        # cat('L=',L,'\n')
        dim(dif) <- c(2, L)
        # cat('dif =',dif,' \n')
        pol <- array(diag(kvar), dim = c(kvar, kvar, 1))
        for (i in 1:L) {
            k <- dif[1, i]
            d <- dif[2, i]
            # cat('k,d=',c(k,d),' \n')
            one <- one.poly(kvar = kvar, dif = c(k, d))
            LL <- dim(pol)[3] + dim(one)[3]
            # cat('LL =',LL,' \n')
            pol <- pol.mul(pol, one, LL)
        }
        for (i in 1:LL) {
            if (sum(colSums(matrix(pol[, , i]))) != 0)
                J <- i
        }
        dif.poly <- pol[, , (1:J)]
    }
    dif.series <- arma.filter(series, ar.poly = dif.poly)$residuals
    L <- dim(dif.series)

    difs <- D[1, ]
    # cat('difs =',difs,'\n')
    vars <- rep(1, kvar)
    # cat('reconstruct means:',difs,'\n')
    vars[difs] <- 0
    # cat('reconstruct vars :',vars,'\n')

    for (i in 1:length(vars)) {
        if (vars[i] == 1) {
            dif.series[i, ] <- dif.series[i, ] + Means[i]
        }
    }
    
    # cat('means=',Means) cat('vars =',vars)
    # if(omit!=TRUE){
    # cat('All data used. Differencing with respect to previous', 'values \n',
    # (negative time index) was done relative to the mean of the series.\n')}
    # if(omit==TRUE){dif.series<-dif.series[ ,(OMIT+1):L[2]]
    # cat('Using data from obs ',(OMIT+1),' to obs ', L[2], ' incl. \n')}
    
    rownames(dif.series) <- Names
    Or <- pol.order(dif.poly)
    # cat('Order=',Or,'\n')
    y.lost <- dif.series[, 1:Or]
    y.dif <- dif.series[, (Or + 1):dim(dif.series)[2]]

    # if (T==1){dif.series<-t(dif.series)}

    names(Means) <- Names

    return(list(y.dif = y.dif, y.lost = y.lost,
                dif.poly = dif.poly, averages = Means,
                dif.series = dif.series))
}

one.poly <- function(kvar = 1, dif = c(0, 0)) {
    d <- dif[2]
    k <- dif[1]
    pol <- array(0, dim = c(kvar, kvar, d + 1))
    pol[, , 1] <- diag(kvar)
    pol[k, k, (d + 1)] <- -1
    return(pol)
}

##' @title define.sum
##' 
##' @description Function to aggregate multivariate time series.
##' Reverse of function 'define.dif'.
##' 
##' @param series = series to be summed up.
##' @param difference = differencing pattern (see define.dif).
##' @param averages of the individual series that (usually) have
##' been subtracted when differencing the time series (if so,
##' the averages are supplied in the output from define.dif(...).
##'
##' @return sum.series = the summed series.
##'
##' @examples
##' set.seed(4711)
##' y<-round(matrix(100*rnorm(48),nrow=4))
##' difference=matrix(c(1,1, 1,1, 2,1, 3,6),nrow=2)
##' dy<-define.dif(y,difference)$dif.series
##' averages<-define.dif(y,difference)$averages
##' sum.y<-define.sum(dy,difference,averages)$series.sum
##' y
##' dy
##' averages
##' sum.y
##'
##' @export
##' 

define.sum <- function(series = NULL, difference = NULL, averages = 0) {
    d <- dim(series)
    if (d[1] > d[2]) {
        series <- t(series)
    }
    d <- dim(series)
    difference <- matrix(difference, nrow = 2)
    # cat('difference=',difference,'\n')
    s <- dim(difference)
    for (j in (1:s[2])) {
        J <- difference[1, j]
        D <- difference[2, j]
        # cat('Summation of var',J,'T=',D,'\n')
        series[J, ] <- ser.sum(series[J, ], s = D)
    }
    if (max(abs(averages)) > 0) {
        for (i in 1:d[1]) {
            series[i, ] <- series[i, ] + averages[i]
        }
    }
    return(list(series.sum = series))
}


##' @title season.lagging
##' 
##' @description Generate new time series with (seasonally) lagged variables
##' from lagging pattern.
##'
##' @param y = data series
##' @param lagging = lagging array array describing what to be added
##' to y: c(1,3,6) adds a new y3, using y1 lagged 6 time steps.
##' lagging<-matrix(c(1,3,6-1,2,4,12-1),nrow=3) adds two new variables
##' (y3 and y4) using y1 lagged 6-1 time steps and y2 lagged 12-1 time
##' steps.
##'
##' @return y.lagged = the part of the new series (including new
##' lagged variables) that can be entered into marima
##' @return y.future = the part of the new series (including new
##' lagged variables) that does not include future observation
##' @return y.lost = previous values of the time series that is
##' incomplete with respect to the new variables generated by lagging
##'
##' cbind(y.lost,y.lagged.y,y.future) is the complete series after
##' creation and addition of the lagged variables.
##'
##' @examples
##' set.seed(4711)
##' # generate bivariate time series
##' y<-round(matrix(10*rnorm(36),nrow=2))
##' y
##' # define new lagged variables (y3 and y4) with seasonalities 6 and 12
##' lagging <-c(1,3,(6-1), 2,4,(12-1)) #
##' season.lagging(y,lagging)
##'
##' @export

season.lagging <- function(y, lagging = NULL) {
    if (is.null(lagging)) {
        y.lagged <- y
        y.lost <- NULL
        y.future <- NULL
    }

    if (!is.null(lagging)) {
        lagging <- matrix(lagging, nrow = 3)
        cat(" lagging definitions = ", c(lagging), " \n")
        Ny <- dim(lagging)[2]
        k <- dim(y)[1]
        s <- dim(y)[2]
        if (s < k) {
            y <- t(y)
        }
        k <- dim(y)[1]
        s <- dim(y)[2]

        L <- max(lagging[3, ])
        if (L >= s) {
            stop("Returning from function saeson.lagging: \n","max lagging",
            L,
            "longer than or equal to length of time series (", s, ") \n")
        }

        cat(" max lagging = ", L, " (OK) \n")

        newy <- matrix(NA, ncol = (s + L), nrow = k)
        newy[1:k, 1:s] <- y
        newy <- rbind(newy, matrix(NA, ncol = (s + L), nrow = k))

        for (i in 1:Ny) {
            old <- lagging[1, i]
            new <- lagging[2, i]
            lag <- lagging[3, i]
            seq <- c(1:s)
            # cat('i=',i,'old=',old,'new=',new,'lag=',lag,'\n')
            newy[new, seq + lag] <- newy[old, seq]
            if (lag < L) {
                ext <- c(1:(L - lag))
                newy[new, s + lag + ext] <- newy[old, ext]
            }
        }

        rownames(newy) <- c(c(paste("y", c(1:k), sep = "")),
               c(paste("y", c(lagging[1, ]), "-",
                       c(lagging[3, ]), sep = "")))

        e <- dim(newy)[2]
        colnames(newy) <- paste("t=", c(1:e), sep = "")

        y.lost <- newy[, c(1:L)]
        y.future <- newy[, c((s + 1):(e))]
        y.lagged <- newy[, c((L + 1):(s))]
    }
    
    return(list(y.lagged = y.lagged,
                y.future =
                y.future,
                y.lost = y.lost))
}


ser.sum <- function(y, s = 1) {
    # y = time series (univariate) s = summation for time difference=s
    # (default s=1)
    z <- y
    for (i in 1:length(y)) {
        j <- i - s
        if (j > 0) {
            z[i] <- z[i] + z[j]
        }
    }
    return(z)
}

ser.dif <- function(y, s = 1) {
    # y = time series (univariate) s = differencing for time=s
    z <- y
    for (i in 1:length(y)) {
        j <- i - s
        if (j > 0) {
            z[i] <- y[i] - y[j]
        }
    }
    return(z)
}
