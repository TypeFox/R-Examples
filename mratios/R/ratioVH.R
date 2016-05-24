

sci.ratioVH <- function(formula, data, type = "Dunnett", base = 1, method = "Plug", 
    Num.Contrast = NULL, Den.Contrast = NULL, alternative = "two.sided", 
    conf.level = 0.95, names = TRUE) 
{
   
    method <- match.arg(method, choices = c("Plug", "Bonf", "Unadj"))
    alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))

    if (length(formula) != 3) {stop("Argument 'formula' mis-specified")}

    mf <- model.frame(formula, data)

    if (ncol(mf) != 2) 
     {stop("Specify one response variable and only one class variable in the formula")}

    if (is.numeric(mf[, 1]) == FALSE)
     {stop("Response variable must be numeric")}

    Response <- mf[, 1]
    Treatment <- droplevels(as.factor(mf[, 2]))
    varnames <- levels(Treatment)
    k <- length(varnames)
    splitdat <- split(Response, Treatment)
    ni <- as.numeric(lapply(splitdat, FUN = length))


   if(any(ni<2))
    {stop("the number of observations in each group should be at least 2")}


    if (is.null(Num.Contrast) == FALSE || is.null(Num.Contrast) == 
        FALSE) {
        if (is.null(Den.Contrast == TRUE) && is.null(Den.Contrast == 
            FALSE)) {
            stop("Num.Contrast is specified, but Den.Contrast is missing")
        }
        if (is.null(Den.Contrast == FALSE) && is.null(Den.Contrast == 
            TRUE)) {
            stop("Den.Contrast is specified, but Num.Contrast is missing")
        }
        if (is.null(Den.Contrast) == FALSE && is.null(Num.Contrast) == 
            FALSE) {
            if (nrow(Den.Contrast) != nrow(Num.Contrast)) {
                stop("number of rows in Num.Contrast and Den.Contrast is not the same")
            }
            if (ncol(Den.Contrast) != k || ncol(Num.Contrast) != 
                k) {
                stop("number of columns in Num.Contrast or Den.Contrast is not the same as number of groups")
            }
            NC0 <- apply(X = Num.Contrast, MARGIN = 1, function(x) {
                all(x == 0)
            })
            DC0 <- apply(X = Den.Contrast, MARGIN = 1, function(x) {
                all(x == 0)
            })
            if (any(c(NC0, DC0))) {
                cat("Warning: At least one row of the numerator or denominator contrast matrices is a vector with all components equal to zero", 
                  "\n")
            }
            Num.C <- Num.Contrast
            Den.C <- Den.Contrast
            type <- "User defined"
            if (is.null(rownames(Num.C)) && is.null(rownames(Den.C))) {
                compnames <- paste("C", 1:nrow(Num.C), sep = "")
            }
            else {
                if (any(rownames(Num.C) != rownames(Den.C))) {
                  compnames <- paste(rownames(Num.C), rownames(Den.C), 
                    sep = "/")
                }
                else {
                  compnames <- rownames(Num.C)
                }
            }
        }
    }
    else {
        type <- match.arg(type, choices = c("Dunnett", "Tukey", 
            "Sequen", "AVE", "GrandMean", "Changepoint", "Marcus", 
            "McDermott", "Williams", "UmbrellaWilliams"))
        if (names == TRUE) {
            names(ni) <- varnames
        }
        Cmat <- contrMatRatio(n = ni, type = type, base = base)
        Num.C <- Cmat$numC
        Den.C <- Cmat$denC
        compnames <- Cmat$rnames
    }
    out <- sci.ratioII(Response = Response, Treatment = Treatment, 
        Num.Contrast = Num.C, Den.Contrast = Den.C, alternative = alternative, 
        conf.level = conf.level, method = method)

    out$type <- type
    out$compnames <- compnames
    colnames(out$Num.Contrast) <- colnames(out$Den.Contrast) <- varnames
    rownames(out$conf.int) <- compnames
    if (type == "User defined")
    {rownames(out$estimate) <- compnames}



if(method=="Unadj")
{
methodname<-paste( round(conf.level*100,4), "% confidence intervals (heteroscedasticity)", sep="")
}
else{
methodname<-paste("Simultaneous ", round(conf.level*100,4), "% confidence intervals (heteroscedasticity)", sep="")
}

     out$methodname<-methodname


    class(out) <- "sci.ratio"
    return(out)
}



sci.ratioII <- 
function (Response, Treatment, Num.Contrast, Den.Contrast, alternative = "two.sided", 
    conf.level = 0.95, method = "Plug") 
{
    CMat <- Num.Contrast
    DMat <- Den.Contrast
    n.Treat <- tapply(Response, Treatment, length)
    Mean.Treat <- tapply(Response, Treatment, mean)
    Var.Treat <- tapply(Response, Treatment, var)

if(!is.numeric(conf.level) | length(conf.level)!=1 | conf.level<=0.5 | conf.level>=1)
 {stop("Argument 'conf.level' must be a single numeric value between 0.5 and 1")}

if(any( sqrt(Var.Treat) < 10 * .Machine$double.eps * abs(Mean.Treat))) 
 {warning("Data are essentially constant in a least one group")}

if(any( n.Treat < 2 )) 
 {warning("There are less than 2 observations in a least one group")}

#if(any( Mean.Treat<0 ))
# {warning("At least one sample mean is negative. Note, that with negative denominators in the ratio of interest, tests with one-sided alternatives may test the incorrect direction.")}

    MMH <- diag(Var.Treat/n.Treat)  # Diagonal matrix containing the variances divided by the sample sizes
    n.comp <- nrow(CMat)
    degree.f <- cpUAd <- cpBon <- cpMtI <- Cplug <- as.numeric(rep(NA,n.comp))
    gammaC.vec <- CMat %*% Mean.Treat/DMat %*% Mean.Treat

    for (i in 1:n.comp) {
      dfNum.i  <- (  sum(  ((CMat[i,] - gammaC.vec[i]*DMat[i,])^2)*(Var.Treat/n.Treat)  )  )^2
      dfDen.i  <-    sum(  ((CMat[i,] - gammaC.vec[i]*DMat[i,])^4)*((Var.Treat/n.Treat)^2)/(n.Treat - 1)  )
      degree.f[i] <- max(round(dfNum.i/dfDen.i, 4),2) # Minimal degrees of freedom = 2 to avoid adj. p.val. < raw p.val.
    }

    CorrMat.plug <- matrix(as.numeric(rep(NA, n.comp * n.comp)), nrow = n.comp)
    for (i in 1:n.comp) {
        for (j in 1:n.comp) {
            CorrMat.plug[i, j] <- (gammaC.vec[i] * DMat[i, ] - 
                CMat[i, ]) %*% MMH %*% (gammaC.vec[j] * DMat[j, 
                ] - CMat[j, ])/(sqrt((gammaC.vec[i] * DMat[i, 
                ] - CMat[i, ]) %*% MMH %*% (gammaC.vec[i] * DMat[i, 
                ] - CMat[i, ])) * sqrt((gammaC.vec[j] * DMat[j, 
                ] - CMat[j, ]) %*% MMH %*% (gammaC.vec[j] * DMat[j, 
                ] - CMat[j, ])))
        }
    }

    Quad.root <- function(Aj, Bj, Cj) {
        Discrimi <- Bj^2 - 4 * Aj * Cj
        if ((Aj > 0) & (Discrimi >= 0))
            Limit.s <- (-Bj + plus.minus * sqrt(Discrimi))/(2 * 
                Aj)
        else Limit.s <- as.numeric(NA)
        return(Limit.s)
    }


switch(method,

# Unadjusted CIs:

Unadj = {

for (i in 1:n.comp) {
    if (alternative == "two.sided") {
        side <- 2
        plus.minus <- c(-1, 1)
        cpUAd[i] <- qt(1 - (1 - conf.level)/(side), degree.f[i], lower.tail = TRUE)
    }
    if ((alternative == "less") | (alternative == "greater")) {
        side <- 1
        if (alternative == "less") 
            plus.minus <- 1
        else plus.minus <- -1
        cpUAd[i] <- qt(1 - (1 - conf.level)/(side), degree.f[i], lower.tail = TRUE)

    }
}
        UAdCL <- matrix(as.numeric(rep(NA, side * n.comp)), nrow = n.comp)
        for (j in 1:n.comp) {
            AjUAd <- (DMat[j, ] %*% Mean.Treat)^2 - (cpUAd[j]^2) * 
                DMat[j, ] %*% MMH %*% DMat[j, ]
            BjUAd <- -2 * ((CMat[j, ] %*% Mean.Treat) * (DMat[j, 
                ] %*% Mean.Treat) - (cpUAd[j]^2) * 
                CMat[j, ] %*% MMH %*% DMat[j, ])
            CjUAd <- (CMat[j, ] %*% Mean.Treat)^2 - (cpUAd[j]^2) * 
                CMat[j, ] %*% MMH %*% CMat[j, ]
            UAdCL[j, ] <- Quad.root(AjUAd, BjUAd, CjUAd)
        }
        sci.table <- data.frame(UAdCL)
	df <- degree.f; critp <- cpUAd
    },


# Bonferroni CIs

 Bonf = {

for (i in 1:n.comp)
 {
    if (alternative == "two.sided") {
        side <- 2
        plus.minus <- c(-1, 1)
        cpBon[i] <- qt(1 - (1 - conf.level)/(side * n.comp), degree.f[i], lower.tail = TRUE)

    }
    if ((alternative == "less") | (alternative == "greater")) {
        side <- 1
        if (alternative == "less") 
            plus.minus <- 1
        else plus.minus <- -1
        cpBon[i] <- qt(1 - (1 - conf.level)/(side * n.comp), degree.f[i], lower.tail = TRUE)
    }
 }
        BonCL <- matrix(as.numeric(rep(NA, side * n.comp)), nrow = n.comp)
        for (j in 1:n.comp) {
            AjBon <- (DMat[j,]%*%Mean.Treat)^2 - (cpBon[j]^2)*DMat[j,]%*%MMH%*%DMat[j,]
            BjBon <- -2*((CMat[j,]%*%Mean.Treat)*(DMat[j,]%*%Mean.Treat) - (cpBon[j]^2)*CMat[j,]%*%MMH%*%DMat[j,])
            CjBon <- (CMat[j,]%*%Mean.Treat)^2 - (cpBon[j]^2)*CMat[j,]%*%MMH%*%CMat[j,]
            BonCL[j,]  <- Quad.root(AjBon, BjBon,  CjBon)
        }
        sci.table <- data.frame(BonCL)
	df <- degree.f; critp <- cpBon
    },


 Plug = {

for (i in 1:n.comp) {
    if (alternative == "two.sided") {
        side <- 2
        plus.minus <- c(-1, 1)

        Cplug[i] <- qmvt(conf.level, df = as.integer(degree.f[i]), 
            corr = CorrMat.plug, delta = rep(0, n.comp), tail = "both", 
            abseps = 1e-05)$quantile
    }
    if ((alternative == "less") | (alternative == "greater")) {
        side <- 1
        if (alternative == "less") 
            plus.minus <- 1
        else plus.minus <- -1

        Cplug[i] <- qmvt(conf.level, df = as.integer(degree.f[i]), 
            corr = CorrMat.plug, delta = rep(0, n.comp), tail = "lower.tail", 
            abseps = 1e-05)$quantile
    }
}
        PlugCL <- matrix(as.numeric(rep(NA, side * n.comp)), nrow = n.comp)
        for (j in 1:n.comp) {
            AjPlug <- (DMat[j,]%*%Mean.Treat)^2 - (Cplug[j]^2)*DMat[j,]%*%MMH%*%DMat[j,]
            BjPlug <- -2*((CMat[j,]%*%Mean.Treat)*(DMat[j,]%*%Mean.Treat) - (Cplug[j]^2)*CMat[j,]%*%MMH%*%DMat[j,])
            CjPlug <- (CMat[j,]%*%Mean.Treat)^2 - (Cplug[j]^2)*CMat[j,]%*%MMH%*%CMat[j,]
            PlugCL[j,] <- Quad.root(AjPlug, BjPlug,  CjPlug)
        }
        sci.table <- data.frame(PlugCL)
	df <- as.integer(degree.f); critp <- Cplug
    })
    if (alternative == "two.sided") {
        names(sci.table) <- c("lower", "upper")
    }
    if (alternative == "less") {
        names(sci.table) <- c("upper")
    }
    if (alternative == "greater") {
        names(sci.table) <- c("lower")
    }

    if (any(is.na(sci.table)))
     {NSD <- TRUE}
    else
     {NSD <- FALSE}

    list(estimate = gammaC.vec, CorrMat.est = CorrMat.plug, Num.Contrast = CMat, 
        Den.Contrast = DMat, conf.int = sci.table, NSD = NSD, 
        method = method, alternative = alternative, conf.level = conf.level, df=df, quantile=critp)
}

#####################################################


simtest.ratioVH <- function(formula, data, type = "Dunnett", base = 1, alternative = "two.sided", 
    Margin.vec = NULL, FWER = 0.05, Num.Contrast = NULL, Den.Contrast = NULL, 
    names = TRUE) 
{
    alternative <- match.arg(alternative, choices = c("two.sided", 
        "less", "greater"))
    if (length(formula) != 3) {
        stop("formula mis-specified")
    }
    mf <- model.frame(formula, data)
    if (ncol(mf) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(mf[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    Response <- mf[, 1]
    Treatment <- as.factor(mf[, 2])
    varnames <- levels(Treatment)
    k <- length(varnames)
    splitdat <- split(Response, Treatment)
    ni <- as.numeric(lapply(splitdat, FUN = length))

   if(any(ni<2))
    {stop("the number of observations in each group should be at least 2")}

    if (is.null(Num.Contrast) == FALSE || is.null(Num.Contrast) == 
        FALSE) {
        if (is.null(Den.Contrast == TRUE) && is.null(Den.Contrast == 
            FALSE)) {
            stop("Num.Contrast is specified, but Den.Contrast is missing")
        }
        if (is.null(Den.Contrast == FALSE) && is.null(Den.Contrast == 
            TRUE)) {
            stop("Den.Contrast is specified, but Num.Contrast is missing")
        }
        if (is.null(Den.Contrast) == FALSE && is.null(Num.Contrast) == 
            FALSE) {
            if (nrow(Den.Contrast) != nrow(Num.Contrast)) {
                stop("number of rows in Num.Contrast and Den.Contrast should be the same")
            }
            if (ncol(Den.Contrast) != k || ncol(Num.Contrast) != 
                k) {
                stop("number of columns in Num.Contrast or Den.Contrast should be the same as number of groups")
            }
            NC0 <- apply(X = Num.Contrast, MARGIN = 1, function(x) {
                all(x == 0)
            })
            DC0 <- apply(X = Den.Contrast, MARGIN = 1, function(x) {
                all(x == 0)
            })
            if (any(c(NC0, DC0))) {
                cat("Warning: At least one row of the numerator or denominator contrast matrices is a vector with all components equal to zero", 
                  "\n")
            }
            Num.C <- Num.Contrast
            Den.C <- Den.Contrast
            type <- "User defined"
            if (is.null(rownames(Num.C)) && is.null(rownames(Den.C))) {
                compnames <- paste("C", 1:nrow(Num.C), sep = "")
            }
            else {
                if (any(rownames(Num.C) != rownames(Den.C))) {
                  compnames <- paste(rownames(Num.C), rownames(Den.C), 
                    sep = "/")
                }
                else {
                  compnames <- rownames(Num.C)
                }
            }
        }
    }
    else {
        type <- match.arg(type, choices = c("Dunnett", "Tukey", 
            "Sequen", "AVE", "GrandMean", "Changepoint", "Marcus", 
            "McDermott", "Williams", "UmbrellaWilliams"))
        if (names == TRUE) {
            names(ni) <- varnames
        }
        Cmat <- contrMatRatio(n = ni, type = type, base = base)
        Num.C <- Cmat$numC
        Den.C <- Cmat$denC
        compnames <- Cmat$rnames
    }


if(is.null(Margin.vec))
 {Margin.vec <- rep(1,nrow(Num.C))}
else
 {
  if(is.numeric(Margin.vec) && length(Margin.vec)<=nrow(Num.C))
   {Margin.vec <- cbind(Margin.vec,Num.C)[,1]}
  else{
    Margin.vec <- Margin.vec[1:nrow(Num.C)]
    warning( paste("Margin.vec has more elements than there are comparisons. Only the first ", nrow(Num.C)," elements are used!") )
    }
 }


#    if (is.null(Margin.vec)) {
#        Margin.vec <- rep(1, nrow(Num.C))
#    }
#    else {
#        if (is.numeric(Margin.vec) && length(Margin.vec) <= nrow(Num.C)) {
#            Margin.vec <- cbind(Margin.vec, Num.C)[, 1]
#        }
#        else {
#            stop("Margin.vec must be a single numeric value or numeric vector not longer than nrow of contrasts ")
#        }
#    }

    out <- simtest.ratioII(Response = Response, Treatment = Treatment, 
        alternative = alternative, Margin.vec = Margin.vec, FWER = FWER, 
        Num.Contrast = Num.C, Den.Contrast = Den.C)
    out$type <- type
    out$compnames <- compnames
    colnames(out$Num.Contrast) <- varnames
    colnames(out$Den.Contrast) <- varnames
    names(out$p.value.raw) <- compnames
    names(out$p.value.adj) <- compnames
    names(out$estimate) <- compnames
    names(out$teststat) <- compnames


    out$methodname<-"Tests for ratios of means assuming heterogeneous variances \n"


    class(out) <- "simtest.ratio"
    return(out)
}




simtest.ratioII <- 
function (Response, Treatment, alternative = "two.sided", Margin.vec = NULL, 
    FWER = 0.05, Num.Contrast, Den.Contrast) 
{
    
    CMat <- Num.Contrast
    DMat <- Den.Contrast
    n.Treat <- tapply(Response, Treatment, length)
    ybar.Treat <- tapply(Response, Treatment, mean)
    var.Treat <- tapply(Response, Treatment, var)


if(!is.numeric(FWER) | length(FWER)!=1 | FWER<=0 | FWER>=0.5)
 {stop("Argument 'FWER' must be a single numeric value between 0 and 0.5")}

if(any( sqrt(var.Treat) < 10 * .Machine$double.eps * abs(ybar.Treat))) 
 {warning("Data are essentially constant in a least one group")}

if(any( n.Treat < 2 )) 
 {warning("There are less than 2 observations in a least one group")}

if(any( ybar.Treat<0 ))
 {warning("At least one sample mean is negative. Note, that with negative denominators in the ratio of interest, tests with one-sided alternatives may test the incorrect direction.")}


    MMH <- diag(var.Treat/n.Treat)  # Diagonal matrix containing the variances divided by the sample sizes
    ncomp <- nrow(CMat)
    Ratio.Estimate <- Test.Stat <- P.raw <- P.adjusted <- d.freedom <- Critical.pt <- as.numeric(rep(NA, ncomp))
    for (i in 1:ncomp){
      dfNum.i  <- (  sum(  ((CMat[i,] - Margin.vec[i]*DMat[i,])^2)*(var.Treat/n.Treat)  )  )^2
      dfDen.i  <-    sum(  ((CMat[i,] - Margin.vec[i]*DMat[i,])^4)*((var.Treat/n.Treat)^2)/(n.Treat - 1)  )
      d.freedom[i] <- max(round(dfNum.i/dfDen.i,4),2) # min df = 2 to avoid adj. p.val. < raw p.val.
    }
    CorrMat.H0 <- matrix(as.numeric(rep(NA, ncomp * ncomp)), nrow = ncomp)
    for (i in 1:ncomp) {
        for (j in 1:ncomp) {
            CorrMat.H0[i, j] <- (Margin.vec[i] * DMat[i, ] - 
                CMat[i, ]) %*% MMH %*% (Margin.vec[j] * DMat[j, 
                ] - CMat[j, ])/(sqrt((Margin.vec[i] * DMat[i, 
                ] - CMat[i, ]) %*% MMH %*% (Margin.vec[i] * DMat[i, 
                ] - CMat[i, ])) * sqrt((Margin.vec[j] * DMat[j, 
                ] - CMat[j, ]) %*% MMH %*% (Margin.vec[j] * DMat[j, 
                ] - CMat[j, ])))
        }
    }
    for (i in 1:ncomp) {
        Ratio.Estimate[i] <- (CMat[i,]%*%ybar.Treat)/(DMat[i,]%*%ybar.Treat)
        Test.Stat[i] <- ((CMat[i,] - Margin.vec[i]*DMat[i,])%*%ybar.Treat)/
                        sqrt((CMat[i,] - Margin.vec[i]*DMat[i,])%*%MMH%*%(CMat[i,] - Margin.vec[i]*DMat[i,]))

        if (alternative=='two.sided'){ 
            P.adjusted[i] <- 1 - pmvt(lower=rep(-abs(Test.Stat[i]),ncomp), upper=rep(abs(Test.Stat[i]),ncomp), df=as.integer(d.freedom[i]),corr=CorrMat.H0, abseps = 0.00001)
            Critical.pt[i] <- qmvt(1-FWER, df=as.integer(d.freedom[i]),corr=CorrMat.H0, delta=0,tail='both' ,abseps = 0.00001)$quantile
            P.raw[i]  <-  2*pt(abs(Test.Stat[i]),d.freedom[i],lower.tail=FALSE)
        }
    
        if (alternative=='greater'){
            P.adjusted[i] <- 1 - pmvt(lower=rep(-Inf,ncomp), upper=rep(Test.Stat[i],ncomp), df=as.integer(d.freedom[i]),corr=CorrMat.H0, abseps = 0.00001)
            Critical.pt[i] <- qmvt(1-FWER,  df=as.integer(d.freedom[i]),corr=CorrMat.H0, delta=0,tail='lower.tail' ,abseps = 0.00001)$quantile
            P.raw[i]  <-  pt(Test.Stat[i],d.freedom[i],lower.tail=FALSE)
        }
        if (alternative=='less'){
            P.adjusted[i] <- 1 - pmvt(lower=rep(-Inf,ncomp), upper=rep(-Test.Stat[i],ncomp), df=as.integer(d.freedom[i]),corr=CorrMat.H0,abseps = 0.00001)
            Critical.pt[i] <-  qmvt(1-FWER,  df=as.integer(d.freedom[i]),corr=CorrMat.H0, delta=0,tail='lower.tail' ,abseps = 0.00001)$quantile
            P.raw[i]  <-  pt(Test.Stat[i],d.freedom[i],lower.tail=TRUE)
        }
    }
    return(list(estimate = Ratio.Estimate, teststat = Test.Stat, 
        Num.Contrast = Num.Contrast, Den.Contrast = Den.Contrast, 
        CorrMat = CorrMat.H0, critical.pt = Critical.pt, p.value.raw = P.raw, 
        p.value.adj = P.adjusted, Margin.vec = Margin.vec, alternative = alternative, 
        FWER = FWER, df=d.freedom, dfmvt=as.integer(d.freedom)))
}


