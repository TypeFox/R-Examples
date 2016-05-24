#### Les fonctions generiques

setGeneric(name="summary",def=function(object,...){standardGeneric("summary")})


setGeneric(name="effect.size",def=function(object,...){standardGeneric("effect.size")})

setGeneric(name="slidingchart",def=function(object,...){standardGeneric("slidingchart")})











### Les fonctions classiques


winvar.Z<-
function(tr=0.2){
if(tr>=0.5){return(0)}
else{
if(tr<=0){return(1)}
else{
f<-function(x){x^2*dnorm(x)}
integrate(f,lower=qnorm(tr),upper=qnorm(1-tr))$value+2*tr*qnorm(tr)^2
}
}
}

summaryhelper<-function(x,tr=0.2){
res<-numeric(11)
res[1]<-length(x)
res[2]<- mean(x)
res[3] <- median(x)
res[4] <- mean(x, tr=tr)
res[5] <- sd(x)
res[6] <- IQR(x)/1.349
res[7] <- mad(x)
res[8] <- mean(abs(x - median(x))) * sqrt(pi/2)
res[9] <- sqrt(winvar(x,tr=tr))/sqrt(winvar.Z(tr=tr))
res[10] <- min(x)
res[11] <- max(x)
names(res) <- c("n", "mean", "median", "trim", "sd", "IQR (*)", 
        "median ad (*)", "mean ad (*)", "sd(w)", "min", "max")
return(res)
}


paired.summary<-function(x,y,tr=0.2){
X<-rbind(summaryhelper(x,tr=tr),summaryhelper(y,tr=tr),summaryhelper(x-y,tr=tr),summaryhelper((x+y)/2,tr=tr))
rownames(X)<-c(paste(deparse(substitute(x)),"(x)",sep=""),paste(deparse(substitute(y)),"(y)",sep=""),"x-y","(x+y)/2")
X
}


winvar<-
function(x,tr=.2,na.rm=FALSE){
#
#  Compute the gamma Winsorized variance for the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
if(na.rm)x<-x[!is.na(x)]
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
y<-ifelse(y<=xbot,xbot,y)
y<-ifelse(y>=xtop,xtop,y)
winvar<-var(y)
winvar
}

wincor<-
function(x,y,tr=.2){
#   Compute the Winsorized correlation between x and y.
#
#   tr is the amount of Winsorization
#   This function also returns the Winsorized covariance
#

g<-floor(tr*length(x))
xvec<-winval(x,tr)
yvec<-winval(y,tr)
wcor<-cor(xvec,yvec)
wcov<-var(xvec,yvec)
if(sum(x==y)!=length(x)){
test<-wcor*sqrt((length(x)-2)/(1.-wcor^2))
sig<-2*(1-pt(abs(test),length(x)-2*g-2))
}
list(cor=wcor,cov=wcov,siglevel=sig)
}

winval<-
function(x,tr=.2){
#
#  Winsorize the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
#  This function is used by several other functions that come with this book.
#
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
winval<-ifelse(x<=xbot,xbot,x)
winval<-ifelse(winval>=xtop,xtop,winval)
winval
}


### Programmer une fonction t.test pour un objet de type paired

t.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("t.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("t.test", c(DATA, paired=TRUE,list(...)))
y$data.name <- DNAME
    y
}



#### Fonction globale yuen.test sur le modele de t.test

yuen1.test<-
function (x, tr = 0.2, conf.level = 0.95, mu = 0, alternative = c("two.sided", 
    "less", "greater")) 
{
    alternative <- match.arg(alternative)
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    alpha <- 1 - conf.level
    n <- length(x)
    g <- floor(tr * n)
    df <- n - 2 * g - 1
    sw <- sqrt(winvar(x, tr))
    se <- sw/((1 - 2 * tr) * sqrt(n))
    dif <- mean(x, tr)
    test <- (dif - mu)/se
    if (alternative == "less") {
        crit <- qt(1 - alpha, df)
        low <- -Inf
        up <- dif + crit * se
        pval <- pt(test, df)
    }
    if (alternative == "greater") {
        crit <- qt(1 - alpha, df)
        low <- dif - crit * se
        up <- Inf
        pval <- 1 - pt(test, df)
    }
    if (alternative == "two.sided") {
        crit <- qt(1 - alpha/2, df)
        low <- dif - crit * se
        up <- dif + crit * se
        pval <- 2 * (1 - pt(abs(test), df))
    }
    estimate <- dif
    tstat <- test
    cint <- c(low, up)
    method <- paste("One sample Yuen test, trim=", tr, sep = "")
    names(estimate) <- c("trimmed mean of x")
    dname <- deparse(substitute(x))
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- c("trimmed means")
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
        conf.int = cint, estimate = estimate, null.value = mu, 
        alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}


yuenp.test <-
function (x, y = NULL, tr = 0.2, alternative = c("two.sided", "less", "greater"), 
    mu = 0, conf.level = 0.95) 
{
    alternative <- match.arg(alternative)
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x, y)
    x <- x[xok]
    y <- y[yok]
    h1 <- length(x) - 2 * floor(tr * length(x))
    q1 <- (length(x) - 1) * winvar(x, tr)
    q2 <- (length(y) - 1) * winvar(y, tr)
    q3 <- (length(x) - 1) * wincor(x, y, tr)$cov
    df <- h1 - 1
    se <- sqrt((q1 + q2 - 2 * q3)/(h1 * (h1 - 1)))
    dif <- mean(x, tr) - mean(y, tr)
    test <- (dif - mu)/se
    alpha <- 1 - conf.level
    if (alternative == "less") {
        crit <- qt(1 - alpha, df)
        low <- -Inf
        up <- dif + crit * se
        pval <- pt(test, df)
    }
    if (alternative == "greater") {
        crit <- qt(1 - alpha, df)
        low <- dif - crit * se
        up <- Inf
        pval <- 1 - pt(test, df)
    }
    if (alternative == "two.sided") {
        crit <- qt(1 - alpha/2, df)
        low <- dif - crit * se
        up <- dif + crit * se
        pval <- 2 * (1 - pt(abs(test), df))
    }
    estimate <- dif
    cint <- c(low, up)
    tstat <- test
    method <- paste("Paired Yuen test, trim=", tr, sep = "")
    names(estimate) <- "trimmed mean of x - trimmed mean of y"
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- c("difference in trimmed means")
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
        conf.int = cint, estimate = estimate, null.value = mu, 
        alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}



yuen2.test<-
function (x, y = NULL, tr = 0.2, alternative = c("two.sided", "less", "greater"), 
    mu = 0, conf.level = 0.95) 
{
    alternative <- match.arg(alternative)
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for two-sample test")
    }
    xok <- complete.cases(x)
    yok <- complete.cases(y)
    x <- x[xok]
    y <- y[yok]

    h1 <- length(x) - 2 * floor(tr * length(x))
    h2 <- length(y) - 2 * floor(tr * length(y))

    d1 <- (length(x) - 1) * winvar(x, tr)/(h1*(h1-1))
    d2 <- (length(y) - 1) * winvar(y, tr)/(h2*(h2-1))

    df <- (d1+d2)^2/(d1^2/(h1-1)+d2^2/(h2-1))

    se <- sqrt(d1+d2)

    dif <- mean(x, tr) - mean(y, tr)

    test <- (dif - mu)/se
    alpha <- 1 - conf.level
    if (alternative == "less") {
        crit <- qt(1 - alpha, df)
        low <- -Inf
        up <- dif + crit * se
        pval <- pt(test, df)
    }
    if (alternative == "greater") {
        crit <- qt(1 - alpha, df)
        low <- dif - crit * se
        up <- Inf
        pval <- 1 - pt(test, df)
    }
    if (alternative == "two.sided") {
        crit <- qt(1 - alpha/2, df)
        low <- dif - crit * se
        up <- dif + crit * se
        pval <- 2 * (1 - pt(abs(test), df))
    }
    estimate <- c(mean(x, tr), mean(y, tr))
    cint <- c(low, up)
    dif <- dif
    tstat <- test
    method <- paste("Two-sample Yuen test, trim=", tr, sep = "")
    names(estimate) <- c("trimmed mean of x", "trimmed mean of y")
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- c("difference in trimmed means")
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
        conf.int = cint, estimate = estimate, null.value = mu, 
        alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}

yuen.t.test<-function(x,...) UseMethod("yuen.t.test")

yuen.t.test.default<-
function (x, y = NULL, tr=0.2, alternative = c("two.sided", "less", "greater"), 
    mu = 0, paired=FALSE,conf.level = 0.95,...){
if(paired){h<-yuenp.test(x=x, y = y, tr = tr,alternative = alternative, 
    mu = mu, conf.level = conf.level)
return(h)}

if(is.null(y)){h<-yuen1.test(x=x, tr = tr,alternative = alternative, 
    mu = mu, conf.level = conf.level)
return(h)}
h<-yuen2.test(x=x, y=y,tr = tr,alternative = alternative, 
    mu = mu, conf.level = conf.level)
return(h)
}



yuen.t.test.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1L]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    y <- do.call("yuen.t.test", c(DATA, list(...)))
    y$data.name <- DNAME
    if(length(y$estimate) == 2L)
        names(y$estimate) <- paste("trimmed mean in group", levels(g))
    y
}




yuen.t.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("yuen.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("yuen.t.test", c(DATA, paired=TRUE,list(...)))
y$data.name <- DNAME
    y
}



paired.effect.size<-
function (x, y,tr=0.2) 
{
    tt <- matrix(numeric(8), nrow = 2)
    tt[1, 1] <- (mean(x) - mean(y))/sqrt((sd(x)^2 + sd(y)^2)/2)
    tt[1, 2] <- (mean(x) - mean(y))/sd(x)
    tt[1, 3] <- (mean(x) - mean(y))/sd(y)
    d <- x - y
    tt[1, 4] <- mean(d)/sd(d)
    Vx <- winvar(x, tr = tr)
    Vy <- winvar(y, tr = tr)
    Vd <- winvar(d, tr = tr)
        std<-sqrt(winvar.Z(tr=tr))
    tt[2, 1] <- std * (mean(x, tr = tr) - mean(y, tr = tr))/sqrt((Vx + 
        Vy)/2)
    tt[2, 2] <- std * (mean(x, tr = tr) - mean(y, tr = tr))/sqrt(Vx)
    tt[2, 3] <- std * (mean(x, tr = tr) - mean(y, tr = tr))/sqrt(Vy)
    tt[2, 4] <- std * mean(d, tr = tr)/sqrt(Vd)
    colnames(tt) <- c("Average", "Single (x)", "Single (y)", 
        "Difference")
    rownames(tt) <- c("OLS", "Robust")
    tt
}


paired.plotCor<-function(df,condition1,condition2,groups=NULL,facet=TRUE,...){
plotP<-ggplot(data=df)+aes_string(x=condition1,y=condition2)
if(is.null(groups)){
plotP+geom_point()+coord_equal()+ geom_abline(intercept =0, slope =1)
}
else{
if(facet){
formula<-paste(groups,"~.",sep="")
plotP+geom_point()+coord_equal()+ geom_abline(intercept =0, slope =1)+facet_grid(formula)
}
else{
plotP+geom_point()+coord_equal()+ geom_abline(intercept =0, slope =1)+aes_string(colour=groups)+ scale_colour_discrete(name = groups)
}
}
}


###une fonction de test

# pitman.morgan.test<-function(x,...) UseMethod("pitman.morgan.test")

pitman.morgan.test.default<-function (x, y = NULL, alternative = c("two.sided", "less", "greater"), 
    ratio = 1, conf.level = 0.95,...) 
{
    alternative <- match.arg(alternative)
    if (!missing(ratio) && (length(ratio) != 1 || is.na(ratio))) 
        stop("' ratio ' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x, y)
    x <- x[xok]
    y <- y[yok]
    alpha = 1 - conf.level
    n <- length(x)
    df = n - 2
    r <- cor(x, y)
    Var1 <- var(x)
    Var2 <- var(y)
    w <- var(x)/var(y)
    stat.t <- ((w - ratio) * sqrt(n - 2))/sqrt(4 * (1 - r^2) * 
        w * ratio)
    if (alternative == "two.sided") {
        k <- qt(1 - alpha/2, df = n - 2)
        K <- 1 + (2 * (1 - r^2) * k^2)/(n - 2)
        low <- w * (K - sqrt(K^2 - 1))
        up <- w * (K + sqrt(K^2 - 1))
        pval <- 2 * (1 - pt(abs(stat.t), df = df))
    }
    if (alternative == "less") {
        k <- qt(1 - alpha, df = n - 2)
        K <- 1 + (2 * (1 - r^2) * k^2)/(n - 2)
        low <- 0
        up <- w * (K + sqrt(K^2 - 1))
        pval <- pt(stat.t, df = df)
    }
    if (alternative == "greater") {
        k <- qt(1 - alpha, df = n - 2)
        K <- 1 + (2 * (1 - r^2) * k^2)/(n - 2)
        low <- w * (K - sqrt(K^2 - 1))
        up <- Inf
        pval <- 1 - pt(stat.t, df = df)
    }
    estimate <- c(Var1, Var2)
    cint <- c(low, up)
    tstat <- stat.t
    method <- c("Paired Pitman-Morgan test")
    names(estimate) <- c("variance of x", "variance of y")
    names(tstat) <- "t"
    names(df) <- "df"
    names(ratio) <- c("ratio of variances")
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
        conf.int = cint, estimate = estimate, null.value = ratio, 
        alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}



var.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("var.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("pitman.morgan.test.default", c(DATA,list(...)))
y$data.name <- DNAME
    y
}

# pitman.morgan.test.paired <-
# function(x,...)
# {
#    if(!is.numeric(x[,1])){return("pitman.morgan.test is only suitable
# to numeric paired data")}
#    DATA <- x
#DNAME <- paste(names(DATA), collapse = " and ")
#    names(DATA) <- c("x", "y")
#    y <- do.call("pitman.morgan.test", c(DATA, list(...)))
# y$data.name <- DNAME
#    y
# }

winsor.cor.test<-function(x,...) UseMethod("winsor.cor.test")

winsor.cor.test.default <-
function(x,y,tr=0.2,alternative = c("two.sided", "less", "greater"),...)
{
alternative <- match.arg(alternative)

method <- paste("winsorized correlation, trim=",tr,sep="")   

NVAL <- 0
names(NVAL)<-"(winsorized) correlation"

DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

xx<-x
yy<-y
g<-floor(tr*length(xx))
xvec<-winval(xx,tr)
yvec<-winval(yy,tr)
wcor<-cor(xvec,yvec)


ESTIMATE <- c(cor = wcor)
PARAMETER <- c(df = length(xx)-2*g-2)
STATISTIC <- c(t=wcor*sqrt((length(xx)-2)/(1.-wcor^2)))

PVAL <-switch(alternative,
              "two.sided" = 2*(1-pt(abs(STATISTIC),length(xx)-2*g-2)),
              "greater" = 1-pt(STATISTIC,length(xx)-2*g-2),
              "less" = pt(STATISTIC,length(xx)-2*g-2))


RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = as.numeric(PVAL),
                 estimate = ESTIMATE,
                 null.value = NVAL,
                 alternative = alternative,
                 method = method,
                 data.name = DNAME)
    
    class(RVAL) <- "htest"
    RVAL

}


winsor.cor.test.paired <-
function(x,tr=0.2,alternative = c("two.sided", "less", "greater"),...)
{
alternative <- match.arg(alternative)

method <- paste("winsorized correlation, trim=",tr,sep="")   

NVAL <- 0
names(NVAL)<-"(winsorized) correlation"

DNAME <- paste(names(x), collapse = " and ")

xx<-x[,1]
yy<-x[,2]
g<-floor(tr*length(xx))
xvec<-winval(xx,tr)
yvec<-winval(yy,tr)
wcor<-cor(xvec,yvec)


ESTIMATE <- c(cor = wcor)
PARAMETER <- c(df = length(xx)-2*g-2)
STATISTIC <- c(t=wcor*sqrt((length(xx)-2)/(1.-wcor^2)))

PVAL <-switch(alternative,
              "two.sided" = 2*(1-pt(abs(STATISTIC),length(xx)-2*g-2)),
              "greater" = 1-pt(STATISTIC,length(xx)-2*g-2),
              "less" = pt(STATISTIC,length(xx)-2*g-2))


RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = as.numeric(PVAL),
                 estimate = ESTIMATE,
                 null.value = NVAL,
                 alternative = alternative,
                 method = method,
                 data.name = DNAME)
    
    class(RVAL) <- "htest"
    RVAL

}




wilcox.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("t.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("wilcox.test", c(DATA, paired=TRUE,list(...)))
y$data.name <- DNAME
    y
}





bonettseier.var.test<-function(x,...) UseMethod("bonettseier.var.test")


bonettseier.var.test.default<-
function (x, y = NULL, alternative = c("two.sided", "less", "greater"), 
    omega = 1, conf.level = 0.95,...) 
{
    alternative <- match.arg(alternative)
    if (!missing(omega) && (length(omega) != 1 || is.na(omega))) 
        stop("'omega' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x, y)
    x <- x[xok]
    y <- y[yok]
    alpha <- 1 - conf.level
    n <- length(x)
    tauX <- mean(abs(x - median(x)))
    tauY <- mean(abs(y - median(y)))
    ratio <- tauX/tauY
    muX <- mean(x)
    sdX <- sd(x)
    gammaX <- (sdX^2)/(tauX^2)
    deltaX <- (muX - median(x))/tauX
    var.ln.tauX <- (gammaX + deltaX^2 - 1)/n
    muY <- mean(y)
    sdY <- sd(y)
    gammaY <- (sdY^2)/(tauY^2)
    deltaY <- (muY - median(y))/tauY
    var.ln.tauY <- (gammaY + deltaY^2 - 1)/n
    dX <- abs(x - median(x))
    dY <- abs(y - median(y))
    corD <- cor(dX, dY)
    var.ln.ratio <- var.ln.tauX + var.ln.tauY - 2 * corD * sqrt(var.ln.tauX * 
        var.ln.tauY)
    stat.t <- (log(ratio) - log(omega))/sqrt(var.ln.ratio)
    if (alternative == "two.sided") {
        z <- qnorm(1 - alpha/2)
        low <- exp(log(ratio) - z * sqrt(var.ln.ratio))
        up <- exp(log(ratio) + z * sqrt(var.ln.ratio))
        pval <- 2 * (1 - pnorm(abs(stat.t)))
    }
    if (alternative == "less") {
        z <- qnorm(1 - alpha)
        low <- 0
        up <- exp(log(ratio) + z * sqrt(var.ln.ratio))
        pval <- pnorm(stat.t, lower.tail = TRUE)
    }
    if (alternative == "greater") {
        z <- qnorm(1 - alpha)
        low <- exp(log(ratio) - z * sqrt(var.ln.ratio))
        up <- Inf
        pval <- pnorm(stat.t, lower.tail = FALSE)
    }
    estimate <- c(tauX, tauY)
    tstat <- stat.t
    cint <- c(low, up)
    method <- c("Paired Bonett-Seier test")
    names(estimate) <- c("mean abs. dev. of x", " mean abs. dev. of y")
    names(tstat) <- "z"
    names(omega) <- c("ratio of means absolute deviations")
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = tstat, p.value = pval, conf.int = cint, 
        estimate = estimate, null.value = omega, alternative = alternative, 
        method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}



bonettseier.var.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("bonett.seier.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("bonettseier.var.test", c(DATA, list(...)))
y$data.name <- DNAME
    y
}



### Encore un autre

grambsch.var.test<-function(x,...) UseMethod("grambsch.var.test")


grambsch.var.test.default<-
function (x, y = NULL, alternative = c("two.sided", "less", "greater"),...) 
{
    alternative <- match.arg(alternative)
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
    xok <- yok <- complete.cases(x, y)
    x <- x[xok]
    y <- y[yok]
    S <- x + y
    D <- x - y
    Z <- (D - mean(D)) * (S - mean(S))
    statistique <- sqrt(length(Z)) * mean(Z)/sd(Z)
    if (alternative == "two.sided") {
        probabilite <- 2 * pnorm(abs(statistique), lower.tail = FALSE)
    }
    if (alternative == "less") {
        probabilite <- pnorm(statistique, lower.tail = TRUE)
    }
    if (alternative == "greater") {
        probabilite <- pnorm(statistique, lower.tail = FALSE)
    }
    tstat <- statistique
    method <- c("Paired Grambsch test")
    pval <- probabilite
    names(tstat) <- "z"
    omega <- 1
    names(omega) <- c("ratio of variances")
    rval <- list(statistic = tstat, p.value = pval, null.value = omega, 
        alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
}

grambsch.var.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("grambsch.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("grambsch.var.test", c(DATA, list(...)))
y$data.name <- DNAME
    y
}


### Et maintenant les fonctions graphiques

paired.plotBA<-
function (df, condition1, condition2, groups = NULL, facet = TRUE, 
    ...) 
{
    name.mean<-paste("(",condition1,"+",condition2,")/2",sep="")
    name.difference<-paste(condition2,"-",condition1,sep="")
    df2 <- df
    df2$M <- (df[, condition1] + df[, condition2])/2
    df2$D <- (df[, condition2] - df[, condition1])
    plotP <- ggplot(data = df2) + aes_string(x = "M", y = "D")
    if (is.null(groups)) {
        plotP + geom_point() + geom_abline(intercept = 0, slope = 0) + 
            geom_smooth(method = "lm", formula = "y~1") + xlab(name.mean) + 
            ylab(name.difference)
    }
    else {
        if (facet) {
            formula <- paste(groups, "~.", sep = "")
            plotP + geom_point() + geom_abline(intercept = 0, 
                slope = 0) + geom_smooth(method = "lm", formula = "y~1", 
                fullrange = TRUE) + xlab(name.mean) + ylab(name.difference) + 
                facet_grid(formula)
        }
        else {
            plotP + geom_point() + geom_abline(intercept = 0, 
                slope = 0) + geom_smooth(method = "lm", formula = "y~1", 
                fullrange = TRUE) + aes_string(colour = groups, 
                fill = groups) + xlab(name.mean) + ylab(name.difference)
        }
    }
}

paired.plotMcNeil<-
function(df, condition1, condition2, groups = NULL, subjects,
    facet = TRUE, ...) 
{

    if (is.null(groups)) {
        df2 <- unroll(df, condition1, condition2, subjects)
    }
    else {
        df2 <- unrollG(df, condition1, condition2, subjects, 
            groups)
    }
    df2[, "Subjects"] <- reorder(df2[, "Subjects"], df2[, "Measurements"])
df2[, "Conditions"]<-ordered(df2[, "Conditions"],levels=c(condition1,condition2))
 
    plotP <- ggplot(data = df2) + aes_string(x = "Measurements", 
        y = "Subjects", group = "Subjects")
    if (is.null(groups)) {
        plotP + geom_point() + aes_string(colour = "Conditions")+xlab("")+  theme(legend.key = element_rect())+ scale_colour_discrete(name = "")+ylab(subjects)
    }
    else {
        if (facet) {
            formula <- "Groups~."
            plotP + geom_point() + aes_string(colour = "Conditions") + 
                facet_grid(formula, scales = "free_y")+xlab("")+ theme(legend.key = element_rect())+ scale_colour_discrete(name = "")+ylab(subjects)
        }
        else {
            plotP + geom_point() + aes_string(colour = "Conditions") + 
                aes_string(shape = "Groups")+xlab("")+ theme(legend.key = element_rect())+ scale_colour_discrete(name = "")+ylab(subjects)+ scale_shape_discrete(name=groups)

        }
    }
}

paired.plotProfiles<-
function (df, condition1, condition2, groups = NULL,subjects, 
    facet = TRUE, ...) 
{
    if (is.null(groups)) {
        df2 <- unroll(df, condition1, condition2, subjects)
    }
    else {
        df2 <- unrollG(df, condition1, condition2, subjects, 
            groups)
    }
    plotP <- ggplot(data = df2) + aes_string(x = "Conditions", 
        y = "Measurements", group = "Subjects")
    if (is.null(groups)) {
        plotP + geom_boxplot(aes(group = NULL), width = 0.1) + 
            geom_line()+xlab("")+ylab("")+xlim(c(condition1,condition2))
    }
    else {
        if (facet) {
            formula <- "Groups~."
            plotP + geom_boxplot(aes(group = NULL), width = 0.1) + 
                geom_line() + facet_grid(formula) +xlab("")+ylab("")+xlim(c(condition1,condition2))
        }
        else {
            plotP + geom_boxplot(aes(group = NULL), width = 0.1) + 
                geom_line(aes_string(colour = "Groups"))+
		  aes_string(fill = "Groups")+xlab("")+ylab("")+
			scale_fill_discrete(name=groups)+ scale_colour_discrete(name=groups) +xlim(c(condition1,condition2))
        }
    }
}




paired.slidingchart<-
function (df, condition1, condition2, xlab = "", ylab = "", ...) 
{
    x <- df[, condition1]
    y <- df[, condition2]
    require(MASS)
    par(mar = c(2, 2, 5, 4))
    mX <- min(x)
    MX <- max(x)
    mY <- min(y)
    MY <- max(y)
    addX <- (max(x) - min(x))/25
    addY <- (max(y) - min(y))/25
    addXY <- min(addX, addY)
    minX <- mX - (MY - mY)/2
    maxX <- MX + (MY - mY)/2
    minY <- mY - (MX - mX)/2
    maxY <- MY
    eqscplot(x, y, xlim = c(minX - addX - addXY, maxX + addX + 
        addXY), ylim = c(minY - addY - addXY, maxY + addY + addXY), 
        axes = FALSE, xlab = "", ylab = "", ...)
    polygon(c(mX - addX, mX - addX, MX + addX, MX + addX), c(mY - 
        addY, MY + addY, MY + addY, mY - addY), lty = 3)
    rug(x, side = 3, ticksize = 0.02)
    rug(y, side = 4, ticksize = 0.02)
    lines(c(minX, (mX + MX)/2 + addXY) - addX, c((mY + MY)/2, 
        minY - addXY) - addY, col = "red", lwd = 0.25)
    lines(c((mX + MX)/2 - addXY, maxX) + addX, c(minY - addXY, 
        (mY + MY)/2) - addY, col = "blue", lwd = 0.25)
    axis(3, at = pretty(x))
    if (xlab == "") {
        xlab <- condition1
    }
    mtext(xlab, side = 3, line = 2)
    axis(4, at = pretty(y))
    if (ylab == "") {
        ylab <- condition2
    }
    mtext(ylab, side = 4, line = 2)
    xy <- (x + y)/2
    lines(c(max(mX - addX, mY - addY), min(MX + addX, MY + addY)), 
        c(max(mX - addX, mY - addY), min(MX + addX, MY + addY)), 
        lty = 2)
    xx <- mX - addX
    yy <- mY - addY
    c <- xx + yy
    points(c/2 - addXY + (x - y)/2, c/2 - addXY + (y - x)/2, 
        col = "red")
    x1 <- MX + addX
    y1 <- mY - addY
    b <- y1 - x1
    e <- b/2
    points(xy - e + addXY, xy + e - addXY, col = "blue")
}



#### divers
unroll<-function(df,c1,c2,subjects){
df2<-data.frame(c(df[,c1],df[,c2]),rep(c(c1,c2),rep(length(df[,c1]),2)),rep(df[,subjects],2))
names(df2)<-c("Measurements","Conditions","Subjects")
df2
}
unrollG<-function(df,c1,c2,subjects,groups){
df2<-data.frame(c(df[,c1],df[,c2]),rep(c(c1,c2),rep(length(df[,c1]),2)),rep(df[,subjects],2),rep(df[,groups],2))
names(df2)<-c("Measurements","Conditions","Subjects","Groups")
df2
}





### la tricherie
var1.test <-
function(x,ratio = 1,
         alternative = c("two.sided", "less", "greater"),
         conf.level = 0.95, ...)
{
    if (!((length(ratio) == 1L) && is.finite(ratio) && (ratio > 0)))
        stop("'ratio' must be a single positive number")

    alternative <- match.arg(alternative)

    if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
          (conf.level > 0) && (conf.level < 1)))
        stop("'conf.level' must be a single number between 0 and 1")

   
 DNAME <- deparse(substitute(x))

    if (inherits(x, "lm")) {
        DF.x <- x$df.residual
        V.x <- sum(x$residuals^2) / DF.x
    } else {
        x <- x[is.finite(x)]
        DF.x <- length(x) - 1L
        if (DF.x < 1L)
            stop("not enough 'x' observations")
        V.x <- var(x)
    }
    ESTIMATE <- V.x
    STATISTIC <- DF.x * ESTIMATE / ratio
    PARAMETER <- c("df" = DF.x)
    PVAL <- pchisq(STATISTIC, DF.x)
    if (alternative == "two.sided") {
        PVAL <- 2 * min(PVAL, 1 - PVAL)
        BETA <- (1 - conf.level) / 2
        CINT <- c(DF.x * ESTIMATE / qchisq(1 - BETA, DF.x),
                  DF.x * ESTIMATE / qchisq(BETA, DF.x))
    }
    else if (alternative == "greater") {
        PVAL <- 1 - PVAL
        CINT <- c(DF.x * ESTIMATE / qchisq(conf.level, DF.x), Inf)
    }
    else
        CINT <- c(0, DF.x * ESTIMATE / qchisq(1 - conf.level, DF.x))
    names(STATISTIC) <- "X-squared"
    names(ESTIMATE) <- names(ratio) <- "variance"
    attr(CINT, "conf.level") <- conf.level
    RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 conf.int = CINT,
                 estimate = ESTIMATE,
                 null.value = ratio,
                 alternative = alternative,
                 method = "One-sample variance test",
                 data.name = DNAME)
    attr(RVAL, "class") <- "htest"
    return(RVAL)
}




var2.test <-
function(x, y, ratio = 1,
         alternative = c("two.sided", "less", "greater"),
         conf.level = 0.95, ...)
{
    if (!((length(ratio) == 1L) && is.finite(ratio) && (ratio > 0)))
        stop("'ratio' must be a single positive number")

    alternative <- match.arg(alternative)

    if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
          (conf.level > 0) && (conf.level < 1)))
        stop("'conf.level' must be a single number between 0 and 1")

    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if (inherits(x, "lm") && inherits(y, "lm")) {
        DF.x <- x$df.residual
        DF.y <- y$df.residual
        V.x <- sum(x$residuals^2) / DF.x
        V.y <- sum(y$residuals^2) / DF.y
    } else {
        x <- x[is.finite(x)]
        DF.x <- length(x) - 1L
        if (DF.x < 1L)
            stop("not enough 'x' observations")
        y <- y[is.finite(y)]
        DF.y <- length(y) - 1L
        if (DF.y < 1L)
            stop("not enough 'y' observations")
        V.x <- var(x)
        V.y <- var(y)
    }
    ESTIMATE <- V.x / V.y
    STATISTIC <- ESTIMATE / ratio
    PARAMETER <- c("num df" = DF.x, "denom df" = DF.y)
    PVAL <- pf(STATISTIC, DF.x, DF.y)
    if (alternative == "two.sided") {
        PVAL <- 2 * min(PVAL, 1 - PVAL)
        BETA <- (1 - conf.level) / 2
        CINT <- c(ESTIMATE / qf(1 - BETA, DF.x, DF.y),
                  ESTIMATE / qf(BETA, DF.x, DF.y))
    }
    else if (alternative == "greater") {
        PVAL <- 1 - PVAL
        CINT <- c(ESTIMATE / qf(conf.level, DF.x, DF.y), Inf)
    }
    else
        CINT <- c(0, ESTIMATE / qf(1 - conf.level, DF.x, DF.y))
    names(STATISTIC) <- "F"
    names(ESTIMATE) <- names(ratio) <- "ratio of variances"
    attr(CINT, "conf.level") <- conf.level
    RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 conf.int = CINT,
                 estimate = ESTIMATE,
                 null.value = ratio,
                 alternative = alternative,
                 method = "F test to compare two variances",
                 data.name = DNAME)
    attr(RVAL, "class") <- "htest"
    return(RVAL)
}




var.test.default<-
function (x, y = NULL, ratio=1, alternative = c("two.sided", "less", "greater"),paired=FALSE,conf.level = 0.95,...){
if(paired){h<- pitman.morgan.test.default(x, y, alternative = alternative, 
    ratio = ratio, conf.level = conf.level) 
return(h)}

if(is.null(y)){h<-var1.test(x=x,alternative = alternative, 
    ratio=ratio, conf.level = conf.level)
return(h)}
h<-var2.test(x=x, y=y,alternative = alternative, 
    ratio=ratio, conf.level = conf.level)
return(h)
}
# scale comparison test
sandvikolsson.var.test<-function(x,...) UseMethod("sandvikolsson.var.test")


sandvikolsson.var.test.default<-
function (x, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, exact = NULL, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95, location=c("trim","median"),tr=0.1, ...)
{
    
if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
location<-match.arg(location)
if(location=="trim"){
X<-abs(x-mean(x,tr=tr))
Y<-abs(y-mean(y,tr=tr))
}
if(location=="median"){
X<-abs(x-median(x))
Y<-abs(y-median(y))
}
rval<-wilcox.test(X,Y,alternative = alternative,
            mu = mu, paired = TRUE, exact = exact, correct = correct,
            conf.int = conf.int , conf.level = conf.level ,...)
rval$data.name<-dname
rval$method<-"Sandvik-Olsson paired test for scale comparison"
    return(rval)
}


sandvikolsson.var.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("sandvikolsson.var.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("sandvikolsson.var.test", c(DATA, list(...)))
y$data.name <- DNAME
    y
}


### levene


levene.var.test<-function(x,...) UseMethod("levene.var.test")


levene.var.test.default<-
function (x, y = NULL,
       alternative = c("two.sided", "less", "greater"),
       mu = 0,conf.level = 0.95, location=c("trim","median"),tr=0.1,...)
{
    
if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
location<-match.arg(location)
if(location=="trim"){
X<-abs(x-mean(x,tr=tr))
Y<-abs(y-mean(y,tr=tr))
}
if(location=="median"){
X<-abs(x-median(x))
Y<-abs(y-median(y))
}

rval<-t.test(X,Y,alternative = alternative,
            mu = mu, paired = TRUE,conf.level = conf.level ,...)
rval$data.name<-dname
rval$method<-"Levene paired test for scale comparison"
    return(rval)
}


levene.var.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("levene.var.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("levene.var.test", c(DATA, list(...)))
y$data.name <- DNAME
    y
}



### imam


imam.var.test<-function(x,...) UseMethod("imam.var.test")


imam.var.test.default<-
function (x, y = NULL,
       alternative = c("two.sided", "less", "greater"),
       mu = 0,conf.level = 0.95, location=c("trim","median"),tr=0.1, ...)
{
    
if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
location<-match.arg(location)
if(location=="trim"){
X1<-abs(x-mean(x,tr=tr))
Y1<-abs(y-mean(y,tr=tr))
}
if(location=="median"){
X1<-abs(x-median(x))
Y1<-abs(y-median(y))
}
n<-length(x)
X1Y1<-rank(c(X1,Y1))
X2<-X1Y1[1:n]
Y2<-X1Y1[(n+1):(2*n)]
rval<-t.test(X2,Y2,alternative = alternative,
            mu = mu, paired = TRUE,conf.level = conf.level ,...)
rval$data.name<-dname
rval$method<-"Imam paired test for scale comparison"
    return(rval)
}


imam.var.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("imam.var.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("imam.var.test", c(DATA, list(...)))
y$data.name <- DNAME
    y
}






### mcculloch


mcculloch.var.test<-function(x,...) UseMethod("mcculloch.var.test")


mcculloch.var.test.default<-
function (x, y = NULL,
       alternative = c("two.sided", "less", "greater"),method= c("spearman","pearson", "kendall"),
       exact = NULL,conf.level = 0.95,continuity = FALSE, ...)
{
    
if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        stop("'y' is missing for paired test")
    }
D<-x-y
S<-x+y
method<-match.arg(method)
alternative<-match.arg(alternative)
rval<-cor.test(D,S,alternative = alternative,
            method = method,
exact = exact, conf.level = conf.level , continuity = continuity,...)
rval$data.name<-dname
rval$method<-"McCulloch paired test for scale comparison"
    return(rval)
}


mcculloch.var.test.paired <-
function(x,...)
{
    if(!is.numeric(x[,1])){return("mcculloch.var.test is only suitable to numeric paired data")}
    DATA <- x
DNAME <- paste(names(DATA), collapse = " and ")
    names(DATA) <- c("x", "y")
    y <- do.call("mcculloch.var.test", c(DATA, list(...)))
y$data.name <- DNAME
    y
}










