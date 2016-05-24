## Copyright (C) 1997-2002  Adrian Trapletti
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

##
## Mostly time series tests
##

runs.test <-
function (x, alternative = c("two.sided", "less", "greater"))
{
    if(!is.factor(x))
        stop("x is not a factor")
    if(any(is.na(x)))
        stop("NAs in x")
    if(length(levels(x)) != 2)
        stop("x does not contain dichotomous data")
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    n <- length(x)
    R <- 1 + sum(as.numeric(x[-1] != x[-n]))
    n1 <- sum(levels(x)[1] == x)
    n2 <- sum(levels(x)[2] == x)
    m <- 1 + 2*n1*n2 / (n1+n2)
    s <- sqrt(2*n1*n2 * (2*n1*n2 - n1 - n2) / ((n1+n2)^2 * (n1+n2-1)))
    STATISTIC <- (R - m) / s
    METHOD <- "Runs Test"
    if(alternative == "two.sided")
        PVAL <- 2 * pnorm(-abs(STATISTIC))
    else if(alternative == "less")
        PVAL <- pnorm(STATISTIC)
    else if(alternative == "greater")
        PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
    else stop("irregular alternative")
    names(STATISTIC) <- "Standard Normal"
    structure(list(statistic = STATISTIC,
                   alternative = alternative,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME),
              class = "htest")
}

bds.test <-
function(x, m = 3, eps = seq(0.5*sd(x),2*sd(x),length=4), trace = FALSE)
{
    if((NCOL(x) > 1) || is.data.frame(x))
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(m < 2)
        stop("m is less than 2")
    if(length(eps) == 0)
        stop("invalid eps")
    if(any(eps <= 0))
        stop("invalid eps")
    DNAME <- deparse(substitute(x))
    n <- length(x)
    k <- length(eps)
    cc <- double(m+1)
    cstan <- double(m+1)
    STATISTIC <- matrix(0,m-1,k)
    for(i in (1:k)) {
        res <- .C("bdstest_main",
                  as.integer(n),
                  as.integer(m),
                  as.vector(x, mode="double"),
                  as.vector(cc),
                  cstan = as.vector(cstan),
                  as.double(eps[i]),
                  as.integer(trace),
                  PACKAGE="tseries")
        STATISTIC[,i] <- res$cstan[2:m+1]
    }
    colnames(STATISTIC) <- eps
    rownames(STATISTIC) <- 2:m
    PVAL <- 2 * pnorm(-abs(STATISTIC))
    colnames(PVAL) <- eps
    rownames(PVAL) <- 2:m
    METHOD <- "BDS Test"
    PARAMETER <- list(m = 2:m, eps = eps)
    structure(list(statistic = STATISTIC,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME,
                   parameter = PARAMETER),
              class = "bdstest")
}

print.bdstest <-
function(x, digits = 4, ...)
{
    if(!inherits(x, "bdstest"))
        stop("method is only for bdstest objects")
    cat("\n\t", x$method, "\n\n")
    cat("data: ", x$data.name, "\n\n")
    if(!is.null(x$parameter)) {
        cat("Embedding dimension = ",
            format(round(x$parameter$m, digits)), sep = " ", "\n\n")
        cat("Epsilon for close points = ",
            format(round(x$parameter$eps, digits)), sep = " ", "\n\n")
    }
    if(!is.null(x$statistic)) {
        colnames(x$statistic) <-
            round(as.numeric(colnames(x$statistic)), digits)
        colnames(x$statistic) <- paste("[",colnames(x$statistic),"]")
        rownames(x$statistic) <-
            round(as.numeric(rownames(x$statistic)), digits)
        rownames(x$statistic) <- paste("[",rownames(x$statistic),"]")
        cat("Standard Normal = \n")
        print(round(x$statistic, digits))
        cat("\n")
    }
    if(!is.null(x$p.value)) {
        colnames(x$p.value) <-
            round(as.numeric(colnames(x$p.value)), digits)
        colnames(x$p.value) <- paste("[",colnames(x$p.value),"]")
        rownames(x$p.value) <-
            round(as.numeric(rownames(x$p.value)), digits)
        rownames(x$p.value) <- paste("[",rownames(x$p.value),"]")
        cat("p-value = \n")
        print(round(x$p.value, digits))
        cat("\n")
    }
    cat("\n")
    invisible(x)
}

adf.test <-
function(x, alternative = c("stationary", "explosive"),
         k = trunc((length(x)-1)^(1/3)))
{
    if((NCOL(x) > 1) || is.data.frame(x))
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(k < 0)
        stop("k negative")
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    k <- k+1
    x <- as.vector(x, mode="double")
    y <- diff(x)
    n <- length(y)
    z <- embed(y, k)
    yt <- z[,1]
    xt1 <- x[k:n]
    tt <- k:n
    if(k > 1) {
        yt1 <- z[,2:k]
        res <- lm(yt ~ xt1 + 1 + tt + yt1)
    }
    else
        res <- lm(yt ~ xt1 + 1 + tt)
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[2,1] / res.sum$coefficients[2,2]
    table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96),
                   c(3.95, 3.80, 3.73, 3.69, 3.68, 3.66),
                   c(3.60, 3.50, 3.45, 3.43, 3.42, 3.41),
                   c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12),
                   c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25),
                   c(0.80, 0.87, 0.90, 0.92, 0.93, 0.94),
                   c(0.50, 0.58, 0.62, 0.64, 0.65, 0.66),
                   c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33))
    table <- -table
    tablen <- dim(table)[2]
    tableT <- c(25, 50, 100, 250, 500, 100000)
    tablep <- c(0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99)
    tableipl <- numeric(tablen)
    for(i in (1:tablen))
        tableipl[i] <- approx(tableT, table[, i], n, rule=2)$y
    interpol <- approx(tableipl, tablep, STAT, rule=2)$y
    if(is.na(approx(tableipl, tablep, STAT, rule=1)$y))
        if(interpol == min(tablep))
            warning("p-value smaller than printed p-value")
        else
            warning("p-value greater than printed p-value")
    if(alternative == "stationary")
        PVAL <- interpol
    else if(alternative == "explosive")
        PVAL <- 1 - interpol
    else stop("irregular alternative")
    PARAMETER <- k-1
    METHOD <- "Augmented Dickey-Fuller Test"
    names(STAT) <- "Dickey-Fuller"
    names(PARAMETER) <- "Lag order"
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   alternative = alternative,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME),
            class = "htest")
}

white.test <- function(x, ...) UseMethod("white.test")

white.test.default <-
function(x, y, qstar = 2, q = 10, range = 4,
         type = c("Chisq","F"), scale = TRUE, ...)
{
    DNAME <- paste(deparse(substitute(x)),
                   "and",
                   deparse(substitute(y)))
    x <- as.matrix(x)
    y <- as.matrix(y)
    if(any(is.na(x))) stop("NAs in x")
    if(any(is.na(y))) stop("NAs in y")
    nin <- dim(x)[2]
    t <- dim(x)[1]
    if(dim(x)[1] != dim(y)[1])
        stop("number of rows of x and y must match")
    if(dim(x)[1] <= 0)
        stop("no observations in x and y")
    if(dim(y)[2] > 1)
        stop("handles only univariate outputs")
    if(!missing(type) && !is.na(pmatch(type, "chisq"))) {
        warning(paste("value 'chisq' for 'type' is deprecated,",
                      "use 'Chisq' instead"))
        type <- "Chisq"
    }
    else
        type <- match.arg(type)
    if(scale) {
        x <- scale(x)
        y <- scale(y)
    }
    xnam <- paste("x[,", 1:nin, "]", sep="")
    fmla <- as.formula(paste("y~",paste(xnam,collapse= "+")))
    rr <- lm(fmla)
    u <- residuals(rr)
    ssr0 <- sum(u^2)
    max <- range/2
    gamma <- matrix(runif((nin+1)*q,-max,max),nin+1,q)
    phantom <- (1+exp(-(cbind(rep.int(1,t),x)%*%gamma)))^(-1)
    phantomstar <- as.matrix(prcomp(phantom,scale=TRUE)$x[,2:(qstar+1)])
    xnam2 <- paste("phantomstar[,", 1:qstar, "]", sep="")
    xnam2 <- paste(xnam2,collapse="+")
    fmla <- as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
                                          xnam2,sep="+")))
    rr <- lm(fmla)
    v <- residuals(rr)
    ssr <- sum(v^2)
    if(type == "Chisq") {
        STAT <- t*log(ssr0/ssr)
        PVAL <- 1-pchisq(STAT,qstar)
        PARAMETER <- qstar
        names(STAT) <- "X-squared"
        names(PARAMETER) <- "df"
    }
    else if(type == "F") {
        STAT <- ((ssr0-ssr)/qstar)/(ssr/(t-qstar-nin))
        PVAL <- 1-pf(STAT,qstar,t-qstar-nin)
        PARAMETER <- c(qstar,t-qstar-nin)
        names(STAT) <- "F"
        names(PARAMETER) <- c("df1","df2")
    }
    else
        stop("invalid type")
    ARG <- c(qstar,q,range,scale)
    names(ARG) <- c("qstar","q","range","scale")
    METHOD <- "White Neural Network Test"
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME,
                   arguments = ARG),
              class = "htest")
}

white.test.ts <-
function(x, lag = 1, qstar = 2, q = 10, range = 4,
         type = c("Chisq","F"), scale = TRUE, ...)
{
    if(!is.ts(x))
        stop("method is only for time series")
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(lag < 1)
        stop("minimum lag is 1")
    if(!missing(type) && !is.na(pmatch(type, "chisq"))) {
        warning(paste("value 'chisq' for 'type' is deprecated,",
                      "use 'Chisq' instead"))
        type <- "Chisq"
    }
    else
        type <- match.arg(type)
    DNAME <- deparse(substitute(x))
    t <- length(x)
    if(scale) x <- scale(x)
    y <- embed(x, lag+1)
    xnam <- paste("y[,", 2:(lag+1), "]", sep="")
    fmla <- as.formula(paste("y[,1]~",paste(xnam,collapse= "+")))
    rr <- lm(fmla)
    u <- residuals(rr)
    ssr0 <- sum(u^2)
    max <- range/2
    gamma <- matrix(runif((lag+1)*q,-max,max),lag+1,q)
    phantom <- (1+exp(-(cbind(rep.int(1,t-lag),y[,2:(lag+1)])%*%gamma)))^(-1)
    phantomstar <- as.matrix(prcomp(phantom,scale=TRUE)$x[,2:(qstar+1)])
    xnam2 <- paste("phantomstar[,", 1:qstar, "]", sep="")
    xnam2 <- paste(xnam2, collapse="+")
    fmla <- as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
                                          xnam2,sep="+")))
    rr <- lm(fmla)
    v <- residuals(rr)
    ssr <- sum(v^2)
    if(type == "Chisq") {
        STAT <- t*log(ssr0/ssr)
        PVAL <- 1-pchisq(STAT,qstar)
        PARAMETER <- qstar
        names(STAT) <- "X-squared"
        names(PARAMETER) <- "df"
    } else if(type == "F") {
        STAT <- ((ssr0-ssr)/qstar)/(ssr/(t-lag-qstar))
        PVAL <- 1-pf(STAT,qstar,t-lag-qstar)
        PARAMETER <- c(qstar,t-lag-qstar)
        names(STAT) <- "F"
        names(PARAMETER) <- c("df1","df2")
    }
    else
        stop("invalid type")
    ARG <- c(lag,qstar,q,range,scale)
    names(ARG) <- c("lag","qstar","q","range","scale")
    METHOD <- "White Neural Network Test"
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME,
                   arguments = ARG),
              class = "htest")
}

terasvirta.test <- function(x, ...) UseMethod("terasvirta.test")

terasvirta.test.default <-
function(x, y, type = c("Chisq", "F"), scale = TRUE, ...)
{
    DNAME <- paste(deparse(substitute(x)),
                   "and",
                   deparse(substitute(y)))
    x <- as.matrix(x)
    y <- as.matrix(y)
    if(any(is.na(x))) stop("NAs in x")
    if(any(is.na(y))) stop("NAs in y")
    nin <- dim(x)[2]
    if(nin < 1)
        stop("invalid x")
    t <- dim(x)[1]
    if(dim(x)[1] != dim(y)[1])
        stop("number of rows of x and y must match")
    if(dim(x)[1] <= 0)
        stop("no observations in x and y")
    if(dim(y)[2] > 1)
        stop("handles only univariate outputs")
    if(!missing(type) && !is.na(pmatch(type, "chisq"))) {
        warning(paste("value 'chisq' for 'type' is deprecated,",
                      "use 'Chisq' instead"))
        type <- "Chisq"
    }
    else
        type <- match.arg(type)
    if(scale) {
        x <- scale(x)
        y <- scale(y)
    }
    xnam <- paste("x[,", 1:nin, "]", sep="")
    fmla <- as.formula(paste("y~",paste(xnam,collapse= "+")))
    rr <- lm(fmla)
    u <- residuals(rr)
    ssr0 <- sum(u^2)
    xnam2 <- NULL
    m <- 0
    for(i in (1:nin)) {
        for(j in (i:nin)) {
            xnam2 <- c(xnam2,paste("I(x[,",i,"]*x[,",j,"])",sep=""))
            m <- m+1
        }
    }
    xnam2 <- paste(xnam2,collapse="+")
    xnam3 <- NULL
    for(i in (1:nin)) {
        for(j in (i:nin)) {
            for(k in (j:nin)) {
                xnam3 <- c(xnam3, paste("I(x[,", i, "]*x[,", j, "]*x[,",
                                        k ,"])", sep=""))
                m <- m+1
            }
        }
    }
    xnam3 <- paste(xnam3,collapse="+")
    fmla <- as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
                                          xnam2,xnam3,sep="+")))
    rr <- lm(fmla)
    v <- residuals(rr)
    ssr <- sum(v^2)
    if(type == "Chisq") {
        STAT <- t*log(ssr0/ssr)
        PVAL <- 1-pchisq(STAT,m)
        PARAMETER <- m
        names(STAT) <- "X-squared"
        names(PARAMETER) <- "df"
    } else if(type == "F") {
        STAT <- ((ssr0-ssr)/m)/(ssr/(t-nin-m))
        PVAL <- 1-pf(STAT,m,t-nin-m)
        PARAMETER <- c(m,t-nin-m)
        names(STAT) <- "F"
        names(PARAMETER) <- c("df1","df2")
    }
    else
        stop("invalid type")
    METHOD <- "Teraesvirta Neural Network Test"
    ARG <- scale
    names(ARG) <- "scale"
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME,
                   arguments = ARG),
              class = "htest")
}

terasvirta.test.ts <-
function(x, lag = 1, type = c("Chisq", "F"), scale = TRUE, ...)
{
    if(!is.ts(x))
        stop("method is only for time series")
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(lag < 1)
        stop("minimum lag is 1")
    if(!missing(type) && !is.na(pmatch(type, "chisq"))) {
        warning(paste("value 'chisq' for 'type' is deprecated,",
                      "use 'Chisq' instead"))
        type <- "Chisq"
    }
    else
        type <- match.arg(type)
    DNAME <- deparse(substitute(x))
    t <- length(x)
    if(scale) x <- scale(x)
    y <- embed(x, lag+1)
    xnam <- paste("y[,", 2:(lag+1), "]", sep="")
    fmla <- as.formula(paste("y[,1]~",paste(xnam,collapse= "+")))
    rr <- lm(fmla)
    u <- residuals(rr)
    ssr0 <- sum(u^2)
    xnam2 <- NULL
    m <- 0
    for(i in (1:lag)) {
        for(j in (i:lag)) {
            xnam2 <- c(xnam2,paste("I(y[,",i+1,"]*y[,",j+1,"])",sep=""))
            m <- m+1
        }
    }
    xnam2 <- paste(xnam2,collapse="+")
    xnam3 <- NULL
    for(i in (1:lag)) {
        for(j in (i:lag)) {
            for(k in (j:lag)) {
                xnam3 <- c(xnam3, paste("I(y[,", i+1, "]*y[,", j+1,
                                        "]*y[,", k+1, "])", sep=""))
                m <- m+1
            }
        }
    }
    xnam3 <- paste(xnam3,collapse="+")
    fmla <- as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
                                          xnam2,xnam3,sep="+")))
    rr <- lm(fmla)
    v <- residuals(rr)
    ssr <- sum(v^2)
    if(type == "Chisq") {
        STAT <- t*log(ssr0/ssr)
        PVAL <- 1-pchisq(STAT,m)
        PARAMETER <- m
        names(STAT) <- "X-squared"
        names(PARAMETER) <- "df"
    }
    else if(type == "F") {
        STAT <- ((ssr0-ssr)/m)/(ssr/(t-lag-m))
        PVAL <- 1-pf(STAT,m,t-lag-m)
        PARAMETER <- c(m,t-lag-m)
        names(STAT) <- "F"
        names(PARAMETER) <- c("df1","df2")
    }
    else
        stop("invalid type")
    METHOD <- "Teraesvirta Neural Network Test"
    ARG <- c(lag,scale)
    names(ARG) <- c("lag","scale")
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME,
                   arguments = ARG),
              class = "htest")
}

jarque.bera.test <-
function(x)
{
    if((NCOL(x) > 1) || is.data.frame(x))
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    DNAME <- deparse(substitute(x))
    n <- length(x)
    m1 <- sum(x)/n
    m2 <- sum((x-m1)^2)/n
    m3 <- sum((x-m1)^3)/n
    m4 <- sum((x-m1)^4)/n
    b1 <- (m3/m2^(3/2))^2
    b2 <- (m4/m2^2)
    STATISTIC <- n*b1/6+n*(b2-3)^2/24
    names(STATISTIC) <- "X-squared"
    PARAMETER <- 2
    names(PARAMETER) <- "df"
    PVAL <- 1 - pchisq(STATISTIC,df = 2)
    METHOD <- "Jarque Bera Test"
    structure(list(statistic = STATISTIC,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME),
              class = "htest")
}

pp.test <-
function(x, alternative = c("stationary", "explosive"),
         type = c("Z(alpha)", "Z(t_alpha)"), lshort = TRUE)
{
    if((NCOL(x) > 1) || is.data.frame(x))
        stop("x is not a vector or univariate time series")
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- as.vector(x, mode="double")
    z <- embed(x, 2)
    yt <- z[,1]
    yt1 <- z[,2]
    n <- length(yt)
    tt <- (1:n)-n/2
    res <- lm(yt ~ 1 + tt + yt1)
    if(res$rank < 3)
        stop("Singularities in regression")
    res.sum <- summary(res)
    u <- residuals(res)
    ssqru <- sum(u^2)/n
    if(lshort)
        l <- trunc(4*(n/100)^0.25)
    else
        l <- trunc(12*(n/100)^0.25)
    ssqrtl <- .C("R_pp_sum",
                 as.vector(u, mode="double"),
                 as.integer(n),
                 as.integer(l),
                 ssqrtl=as.double(ssqru),
                 PACKAGE="tseries")$ssqrtl
    n2 <- n^2
    trm1 <- n2*(n2-1)*sum(yt1^2)/12
    trm2 <- n*sum(yt1*(1:n))^2
    trm3 <- n*(n+1)*sum(yt1*(1:n))*sum(yt1)
    trm4 <- (n*(n+1)*(2*n+1)*sum(yt1)^2)/6
    Dx <- trm1-trm2+trm3-trm4
    if(type == "Z(alpha)") {
        alpha <- res.sum$coefficients[3,1]
        STAT <- n*(alpha-1)-(n^6)/(24*Dx)*(ssqrtl-ssqru)
        table <- cbind(c(22.5, 25.7, 27.4, 28.4, 28.9, 29.5),
                       c(19.9, 22.4, 23.6, 24.4, 24.8, 25.1),
                       c(17.9, 19.8, 20.7, 21.3, 21.5, 21.8),
                       c(15.6, 16.8, 17.5, 18.0, 18.1, 18.3),
                       c(3.66, 3.71, 3.74, 3.75, 3.76, 3.77),
                       c(2.51, 2.60, 2.62, 2.64, 2.65, 2.66),
                       c(1.53, 1.66, 1.73, 1.78, 1.78, 1.79),
                       c(0.43, 0.65, 0.75, 0.82, 0.84, 0.87))
    }
    else if(type == "Z(t_alpha)") {
        tstat <-
            (res.sum$coefficients[3,1]-1)/res.sum$coefficients[3,2]
        STAT <- sqrt(ssqru)/sqrt(ssqrtl)*tstat-(n^3) /
            (4*sqrt(3)*sqrt(Dx)*sqrt(ssqrtl))*(ssqrtl-ssqru)
        table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96),
                       c(3.95, 3.80, 3.73, 3.69, 3.68, 3.66),
                       c(3.60, 3.50, 3.45, 3.43, 3.42, 3.41),
                       c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12),
                       c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25),
                       c(0.80, 0.87, 0.90, 0.92, 0.93, 0.94),
                       c(0.50, 0.58, 0.62, 0.64, 0.65, 0.66),
                       c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33))
    }
    else
        stop("irregular type")
    table <- -table
    tablen <- dim(table)[2]
    tableT <- c(25, 50, 100, 250, 500, 100000)
    tablep <- c(0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99)
    tableipl <- numeric(tablen)
    for(i in (1:tablen))
        tableipl[i] <- approx(tableT, table[, i], n, rule=2)$y
    interpol <- approx(tableipl, tablep, STAT, rule=2)$y
    if(is.na(approx(tableipl, tablep, STAT, rule=1)$y))
        if(interpol == min(tablep))
            warning("p-value smaller than printed p-value")
        else
            warning("p-value greater than printed p-value")
    if(alternative == "stationary")
        PVAL <- interpol
    else if(alternative == "explosive")
        PVAL <- 1 - interpol
    else stop("irregular alternative")
    PARAMETER <- l
    METHOD <- "Phillips-Perron Unit Root Test"
    names(STAT) <- paste("Dickey-Fuller", type)
    names(PARAMETER) <- "Truncation lag parameter"
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   alternative = alternative,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME),
              class = "htest")
}

po.test <-
function(x, demean = TRUE, lshort = TRUE)
{
    if(NCOL(x) <= 1)
        stop("x is not a matrix or multivariate time series")
    DNAME <- deparse(substitute(x))
    x <- as.matrix(x)
    dimx <- ncol(x)
    if(dimx > 6) stop("no critical values for this dimension")
    if(demean)
        res <- lm(x[,1]~x[,-1])
    else
        res <- lm(x[,1]~x[,-1]-1)
    z <- embed(residuals(res), 2)
    ut <- z[,1]
    ut1 <- z[,2]
    n <- length(ut)
    res <- lm(ut ~ ut1 - 1)
    if(res$rank < 1)
        stop("Singularities in regression")
    res.sum <- summary(res)
    k <- residuals(res)
    ssqrk <- sum(k^2)/n
    if(lshort)
        l <- trunc(n/100)
    else
        l <- trunc(n/30)
    ssqrtl <- .C("R_pp_sum",
                 as.vector(k, mode="double"),
                 as.integer(n),
                 as.integer(l),
                 ssqrtl=as.double(ssqrk),
                 PACKAGE="tseries")$ssqrtl
    alpha <- res.sum$coefficients[1,1]
    STAT <- n*(alpha-1)-0.5*n^2*(ssqrtl-ssqrk)/(sum(ut1^2))
    if(demean) {
        table <- cbind(c(28.32, 34.17, 41.13, 47.51, 52.17),
                       c(23.81, 29.74, 35.71, 41.64, 46.53),
                       c(20.49, 26.09, 32.06, 37.15, 41.94),
                       c(18.48, 23.87, 29.51, 34.71, 39.11),
                       c(17.04, 22.19, 27.58, 32.74, 37.01),
                       c(15.93, 21.04, 26.23, 31.15, 35.48),
                       c(14.91, 19.95, 25.05, 29.88, 34.20))
    }
    else {
        table <- cbind(c(22.83, 29.27, 36.16, 42.87, 48.52),
                       c(18.89, 25.21, 31.54, 37.48, 42.55),
                       c(15.64, 21.48, 27.85, 33.48, 38.09),
                       c(13.81, 19.61, 25.52, 30.93, 35.51),
                       c(12.54, 18.18, 23.92, 28.85, 33.80),
                       c(11.57, 17.01, 22.62, 27.40, 32.27),
                       c(10.74, 16.02, 21.53, 26.17, 30.90))
    }
    table <- -table
    tablep <- c(0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)
    PVAL <- approx(table[dimx-1,], tablep, STAT, rule=2)$y
    if(is.na(approx(table[dimx-1, ], tablep, STAT, rule=1)$y))
        if(PVAL == min(tablep))
            warning("p-value smaller than printed p-value")
        else
            warning("p-value greater than printed p-value")
    PARAMETER <- l
    METHOD <- "Phillips-Ouliaris Cointegration Test"
    if(demean)
        names(STAT) <- "Phillips-Ouliaris demeaned"
    else
        names(STAT) <- "Phillips-Ouliaris standard"
    names(PARAMETER) <- "Truncation lag parameter"
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME),
              class = "htest")
}

kpss.test <-
function(x, null = c("Level", "Trend"), lshort = TRUE)
{
    if((NCOL(x) > 1) || is.data.frame(x))        
        stop("x is not a vector or univariate time series")
    DNAME <- deparse(substitute(x))
    null <- match.arg(null)
    x <- as.vector(x, mode="double")
    n <- length(x)
    if(null == "Trend") {
        t <- 1:n
        e <- residuals(lm(x ~ t))
        table <- c(0.216, 0.176, 0.146, 0.119)
    }
    else if(null == "Level") {
        e <- residuals(lm(x ~ 1))
        table <- c(0.739, 0.574, 0.463, 0.347)
    }
    tablep <- c(0.01, 0.025, 0.05, 0.10)
    s <- cumsum(e)
    eta <- sum(s^2)/(n^2)
    s2 <- sum(e^2)/n
    if(lshort)
        l <- trunc(3*sqrt(n)/13)
    else
        l <- trunc(10*sqrt(n)/14)
    s2 <- .C("R_pp_sum",
             as.vector(e, mode="double"),
             as.integer(n),
             as.integer(l),
             s2=as.double(s2),
             PACKAGE="tseries")$s2
    STAT <- eta/s2
    PVAL <- approx(table, tablep, STAT, rule=2)$y
    if(is.na(approx(table, tablep, STAT, rule=1)$y))
        if(PVAL == min(tablep))
            warning("p-value smaller than printed p-value")
        else
            warning("p-value greater than printed p-value")
    PARAMETER <- l
    METHOD <- paste("KPSS Test for", null, "Stationarity")
    names(STAT) <- paste("KPSS", null)
    names(PARAMETER) <- "Truncation lag parameter"
    structure(list(statistic = STAT,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = METHOD,
                   data.name = DNAME),
              class = "htest")
}
