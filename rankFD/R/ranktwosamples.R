#' A function for analyzing two-sample problems
#' 
#' The rank.two.sample() function 
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'    contains the response variable and the right hand side contains the factor
#'    variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in 
#'    \code{formula}. The default option is \code{NULL}.
#' @param conf.level A number specifying the confidence level; the default is 0.95.
#' @param alternative Which alternative is considered? One of "two.sided", "less", "greater".
#' @param rounds Value specifying the number of digits the results are rounded to.
#' @param method specifying the method used for calculation of the confidence intervals.
#'    One of "logit", "probit", "normal", "t.app" and "permu".
#' @param plot.simci Logical, indicating whether or not confidence intervals
#'    should be plotted
#' @param info Logical. info = FALSE suppresses the output of additional information
#'    concerning e.g. the interpretation of the test results.
#' @param wilcoxon asymptotic or exact calculation of Wilcoxon test.
#' @param shift.int Logical, indicating whether or not shift effects should be considered.
#' @param nperm Number of permutations used, default is 10000.
#' 
#' 
#' @examples
#' data(Muco)
#' Muco2 <- subset(Muco, Disease != "OAD")
#' twosample <- rank.two.samples(HalfTime ~ Disease, data = Muco2, 
#'    alternative = "greater", method = "probit", wilcoxon="exact")
#' 
#' 
#' @export


rank.two.samples <- function (formula, data, conf.level = 0.95, 
alternative = c("two.sided", 
    "less", "greater"), rounds = 3, method = c("logit", "probit", 
    "normal", "t.app", "permu"), plot.simci = FALSE, info = TRUE, 
wilcoxon=c("asymptotic","exact"),shift.int=TRUE,
    nperm = 10000) 
{
    alpha <- 1 - conf.level
    if (alpha >= 1 || alpha <= 0) {
        stop("The confidence level must be between 0 and 1!")
        if (is.null(alternative)) {
            stop("Please declare the alternative! (two.sided, less, greater)")
        }
    }
    alternative <- match.arg(alternative)
    method <- match.arg(method)
wilcoxon <- match.arg(wilcoxon)
    if (length(formula) != 3) {
        stop("You can only analyse one-way layouts!")
    }
    dat <- model.frame(formula, droplevels(data))
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    response <- dat[, 1]
    factorx <- as.factor(dat[, 2])
    fl <- levels(factorx)
    a <- nlevels(factorx)
    if (a > 2) {
        stop("You want to perform a contrast test (the factor variable has more than two levels). Please use the function mctp!")
    }
    samples <- split(response, factorx)
    n <- sapply(samples, length)
    n1 <- n[1]
    n2 <- n[2]
    if (any(n == 1)) {
        warn <- paste("The factor level", fl[n == 1], "has got only one observation!")
        stop(warn)
    }
    N <- sum(n)
    cmpid <- paste("p(", fl[1], ",", fl[2], ")", sep = "")
    plotz <- 1
    rxy <- rank(c(samples[[1]], samples[[2]]))
    rx <- rank(c(samples[[1]]))
    ry <- rank(c(samples[[2]]))
    pl1 <- 1/n2 * (rxy[1:n1] - rx)
    pl2 <- 1/n1 * (rxy[(n1 + 1):N] - ry)
    pd <- mean(pl2)
    pd1 <- (pd == 1)
    pd0 <- (pd == 0)
    pd[pd1] <- 0.999
    pd[pd0] <- 0.001
    s1 <- var(pl1)/n1
    s2 <- var(pl2)/n2
    V <- N * (s1 + s2)
    singular.bf <- (V == 0)
    V[singular.bf] <- N/(2 * n1 * n2)
    switch(method, normal = {
        AsyMethod <- "Normal - Approximation"
        T <- sqrt(N) * (pd - 1/2)/sqrt(V)
        switch(alternative, two.sided = {
            text.Output <- paste("True relative effect p is less or equal than 1/2")
            p.Value <- min(2 - 2 * pnorm(T), 2 * pnorm(T))
            crit <- qnorm(1 - alpha/2)
            Lower <- pd - crit/sqrt(N) * sqrt(V)
            Upper <- pd + crit/sqrt(N) * sqrt(V)
        }, less = {
            text.Output <- paste("True relative effect p is less than 1/2")
            p.Value <- pnorm(T)
            crit <- qnorm(1 - alpha)
            Lower <- 0
            Upper <- pd + crit/sqrt(N) * sqrt(V)
        }, greater = {
            text.Output <- paste("True relative effect p is greater than 1/2")
            p.Value <- 1 - pnorm(T)
            crit <- qnorm(1 - alpha)
            Lower <- pd - crit/sqrt(N) * sqrt(V)
            Upper <- 1
        })
        data.info <- data.frame(Sample = fl, Size = n)
        Analysis <- data.frame(Effect = cmpid, Estimator = round(pd, 
            rounds), Lower = round(Lower, rounds), Upper = round(Upper, 
            rounds), T = round(T, rounds), p.Value = round(p.Value, 
            rounds))
        rownames(Analysis) <- 1
    }, t.app = {
        T <- sqrt(N) * (pd - 1/2)/sqrt(V)
        df.sw <- (s1 + s2)^2/(s1^2/(n1 - 1) + s2^2/(n2 - 1))
        df.sw[is.nan(df.sw)] <- 1000
        AsyMethod <- paste("Brunner - Munzel - T - Approx with", 
            round(df.sw, rounds), "DF")
        switch(alternative, two.sided = {
            text.Output <- paste("True relative effect p is less or equal than 1/2")
            p.Value <- min(2 - 2 * pt(T, df = df.sw), 2 * pt(T, 
                df = df.sw))
            crit <- qt(1 - alpha/2, df = df.sw)
            Lower <- pd - crit/sqrt(N) * sqrt(V)
            Upper <- pd + crit/sqrt(N) * sqrt(V)
        }, less = {
            text.Output <- paste("True relative effect p is less than 1/2")
            p.Value <- pt(T, df = df.sw)
            crit <- qt(1 - alpha, df = df.sw)
            Lower <- 0
            Upper <- pd + crit/sqrt(N) * sqrt(V)
        }, greater = {
            text.Output <- paste("True relative effect p is greater than 1/2")
            p.Value <- 1 - pt(T, df = df.sw)
            crit <- qt(1 - alpha, df = df.sw)
            Lower <- pd - crit/sqrt(N) * sqrt(V)
            Upper <- 1
        })
        data.info <- data.frame(Sample = fl, Size = n)
        Analysis <- data.frame(Effect = cmpid, Estimator = round(pd, 
            rounds), Lower = round(Lower, rounds), Upper = round(Upper, 
            rounds), T = round(T, rounds), p.Value = round(p.Value, 
            rounds))
        rownames(Analysis) <- 1
        result <- list(Info = data.info, Analysis = Analysis)
    }, logit = {
        AsyMethod <- "Logit - Transformation"
        logitf <- function(p) {
            log(p/(1 - p))
        }
        expit <- function(G) {
            exp(G)/(1 + exp(G))
        }
        logit.pd <- logitf(pd)
        logit.dev <- 1/(pd * (1 - pd))
        vd.logit <- logit.dev^2 * V
        T <- (logit.pd) * sqrt(N/vd.logit)
        switch(alternative, two.sided = {
            text.Output <- paste("True relative effect p is less or equal than 1/2")
            p.Value <- min(2 - 2 * pnorm(T), 2 * pnorm(T))
            crit <- qnorm(1 - alpha/2)
            Lower <- expit(logit.pd - crit/sqrt(N) * sqrt(vd.logit))
            Upper <- expit(logit.pd + crit/sqrt(N) * sqrt(vd.logit))
        }, less = {
            text.Output <- paste("True relative effect p is less than 1/2")
            p.Value <- pnorm(T)
            crit <- qnorm(1 - alpha)
            Lower <- 0
            Upper <- expit(logit.pd + crit/sqrt(N) * sqrt(vd.logit))
        }, greater = {
            text.Output <- paste("True relative effect p is greater than 1/2")
            p.Value <- 1 - pnorm(T)
            crit <- qnorm(1 - alpha)
            Lower <- expit(logit.pd - crit/sqrt(N) * sqrt(vd.logit))
            Upper <- 1
        })
        data.info <- data.frame(Sample = fl, Size = n)
        Analysis <- data.frame(Effect = cmpid, Estimator = round(pd, 
            rounds), Lower = round(Lower, rounds), Upper = round(Upper, 
            rounds), T = round(T, rounds), p.Value = round(p.Value, 
            rounds))
        rownames(Analysis) <- 1
        result <- list(Info = data.info, Analysis = Analysis)
    }, probit = {
        AsyMethod <- "Probit - Transformation"
        probit.pd <- qnorm(pd)
        probit.dev <- sqrt(2 * pi)/(exp(-0.5 * qnorm(pd) * qnorm(pd)))
        vd.probit <- probit.dev^2 * V
        T <- (probit.pd) * sqrt(N/vd.probit)
        switch(alternative, two.sided = {
            text.Output <- paste("True relative effect p is less or equal than 1/2")
            p.Value <- min(2 - 2 * pnorm(T), 2 * pnorm(T))
            crit <- qnorm(1 - alpha/2)
            Lower <- pnorm(probit.pd - crit/sqrt(N) * sqrt(vd.probit))
            Upper <- pnorm(probit.pd + crit/sqrt(N) * sqrt(vd.probit))
        }, less = {
            text.Output <- paste("True relative effect p is less than 1/2")
            p.Value <- pnorm(T)
            crit <- qnorm(1 - alpha)
            Lower <- 0
            Upper <- pnorm(probit.pd + crit/sqrt(N) * sqrt(vd.probit))
        }, greater = {
            text.Output <- paste("True relative effect p is greater than 1/2")
            p.Value <- 1 - pnorm(T)
            crit <- qnorm(1 - alpha)
            Lower <- pnorm(probit.pd - crit/sqrt(N) * sqrt(vd.probit))
            Upper <- 1
        })
        data.info <- data.frame(Sample = fl, Size = n)
        Analysis <- data.frame(Effect = cmpid, Estimator = round(pd, 
            rounds), Lower = round(Lower, rounds), Upper = round(Upper, 
            rounds), T = round(T, rounds), p.Value = round(p.Value, 
            rounds))
        rownames(Analysis) <- 1
        result <- list(Info = data.info, Analysis = Analysis)
    }, permu = {
   
Tperm=Tlogitperm=Tprobitperm=c()

ausgang = BMstat(samples[[1]],samples[[2]],n1,n2)
for(h in 1:nperm){
respperm=sample(response)
phelp=BMstat(respperm[1:n1],respperm[(n1+1):N],n1,n2)
Tperm[h] = phelp$T
Tlogitperm[h] = phelp$Logit
Tprobitperm[h] = phelp$Probit
}
p.PERM1 = mean(ausgang$T >= Tperm)
p.PERMLogit1 = mean(ausgang$Logit >= Tlogitperm)
p.PERMProbit1 = mean(ausgang$Probit >= Tprobitperm)
c1 = quantile(Tperm,(1-conf.level)/2)
c2 = quantile(Tperm,1-(1-conf.level)/2)

c1LOGIT = quantile(Tlogitperm,(1-conf.level)/2)
c2LOGIT = quantile(Tlogitperm,1-(1-conf.level)/2)
c1PROBIT = quantile(Tprobitperm,(1-conf.level)/2)
c2PROBIT = quantile(Tprobitperm,1-(1-conf.level)/2)
c1lower = quantile(Tperm,(1-conf.level))
c2upper = quantile(Tperm,1-(1-conf.level))
c1LOGITlower = quantile(Tlogitperm,(1-conf.level))
c2LOGITupper = quantile(Tlogitperm,1-(1-conf.level))
c1PROBITlower = quantile(Tprobitperm,(1-conf.level))
c2PROBITupper = quantile(Tprobitperm,1-(1-conf.level))
  

        switch(alternative, two.sided = {
            text.Output <- paste("True relative effect p is less or equal than 1/2")
            p.PERM <- min(2 - 2 * p.PERM1, 2 * p.PERM1)
            p.LOGIT <- min(2 - 2 * p.PERMLogit1, 2 * p.PERMLogit1)
            p.PROBIT <- min(2 - 2 * p.PERMProbit1, 2 *p.PERMProbit1)
        UntenRS <- pd - sqrt(ausgang$sdx/N) * c2
        ObenRS <- pd - sqrt(ausgang$sdx/N) * c1

        ULogitRS <- logit(pd) - ausgang$slogit/sqrt(N) * c2LOGIT
        OLogitRS <- logit(pd) - ausgang$slogit/sqrt(N) * c1LOGIT
        UntenLogitRS <- expit(ULogitRS)
        ObenLogitRS <- expit(OLogitRS)
        UProbitRS <- qnorm(pd) - ausgang$sprobit/sqrt(N) * c2PROBIT
        OProbitRS <- qnorm(pd) - ausgang$sprobit/sqrt(N) * c1PROBIT

        UntenProbitRS <- pnorm(UProbitRS)
        ObenProbitRS <- pnorm(OProbitRS)
        Statistic <- round(c(ausgang$T, ausgang$Logit, ausgang$Probit), rounds)
        Estimator <- round(rep(pd, 3), rounds)
        Lower <- round(c(UntenRS, UntenLogitRS, UntenProbitRS), 
            rounds)
        Upper <- round(c(ObenRS, ObenLogitRS, ObenProbitRS), 
            rounds)
        p.value <- c(p.PERM, p.LOGIT, p.PROBIT)
		Analysis <- data.frame(Estimator, Statistic, Lower, 
            Upper, p.value, row.names = c("id", "logit", "probit"))
        }, less = {
            text.Output <- paste("True relative effect p is less than 1/2")
            p.PERM = p.PERM1
            p.LOGIT <- p.PERMLogit1
            p.PROBIT <- p.PERMProbit1
        UntenRS <- 0
        ObenRS <- pd - sqrt(ausgang$sdx/N) * c1lower
        OLogitRS <- logit(pd) - ausgang$slogit/sqrt(N) * c1LOGITlower
        UntenLogitRS <-0
        ObenLogitRS <- expit(OLogitRS)
        OProbitRS <- qnorm(pd) - ausgang$sprobit/sqrt(N) * c1PROBITlower

        UntenProbitRS <- 0
        ObenProbitRS <- pnorm(OProbitRS)
        Statistic <- round(c(ausgang$T, ausgang$Logit, ausgang$Probit), rounds)
        Estimator <- round(rep(pd, 3), rounds)
        Lower <- round(c(UntenRS, UntenLogitRS, UntenProbitRS), 
            rounds)
        Upper <- round(c(ObenRS, ObenLogitRS, ObenProbitRS), 
            rounds)
        p.value <- c(p.PERM, p.LOGIT, p.PROBIT)
		Analysis <- data.frame(Estimator, Statistic, Lower, 
            Upper, p.value, row.names = c("id", "logit", "probit"))

        }, greater = {
            text.Output <- paste("True relative effect p is greater than 1/2")
            p.PERM = 1-p.PERM1
            p.LOGIT <-1- p.PERMLogit1
            p.PROBIT <- 1-p.PERMProbit1
        UntenRS <- pd - sqrt(ausgang$sdx/N) * c2upper
        ObenRS <- 1
	  ULogitRS <- logit(pd) - ausgang$slogit/sqrt(N) * c2LOGITupper
        UntenLogitRS <-expit(ULogitRS)
        ObenLogitRS <- 1
      

        UProbitRS <- qnorm(pd) - ausgang$sprobit/sqrt(N) * c2PROBITupper
        UntenProbitRS <- pnorm(UProbitRS)
	  ObenProbitRS <- 1
        Statistic <- round(c(ausgang$T, ausgang$Logit, ausgang$Probit), rounds)
        Estimator <- round(rep(pd, 3), rounds)
        Lower <- round(c(UntenRS, UntenLogitRS, UntenProbitRS), 
            rounds)
        Upper <- round(c(ObenRS, ObenLogitRS, ObenProbitRS), 
            rounds)
        p.value <- round(c(p.PERM, p.LOGIT, p.PROBIT), rounds)
		Analysis <- data.frame(Estimator, Statistic, Lower, 
            Upper, p.value, row.names = c("id", "logit", "probit"))


        })

        
        AsyMethod <- "Studentized Permutation Test (+ delta-method)"
        #cmpid <- c("id", "logit", "probit")
        data.info <- data.frame(Sample = fl, Size = n)
       result <- list(Info = data.info, Analysis = Analysis)
    })


#------------SHIFT EFFECTS--------------------------------#
HL.help=expand.grid(samples[[1]],samples[[2]])
HL=median(HL.help[,2]-HL.help[,1])
switch(wilcoxon,asymptotic={
Wilcox = wilcox_test(response~factorx,distribution="asymptotic",
alternative=alternative,
conf.int=TRUE,conf.level=(1 - alpha))
p.wilcox=pvalue(Wilcox)
Z.wilcox=statistic(Wilcox)
if(shift.int==TRUE){
shiftint=sort(-1*c(confint(Wilcox)$conf.int))
Lower.Shift=shiftint[1]
Upper.Shift=shiftint[2]
}
},
exact={
Wilcox = wilcox_test(response~factorx,distribution="exact",
alternative=alternative,conf.int=TRUE,conf.level=(1 - alpha))
p.wilcox=pvalue(Wilcox)
Z.wilcox=sum(ry)
if(shift.int==TRUE){
shiftint=sort(-1*c(confint(Wilcox)$conf.int))
Lower.Shift=shiftint[1]
Upper.Shift=shiftint[2]
}
})




if(shift.int==FALSE){
Lower.Shift = NA
Upper.Shift= NA
HL = NA
}

cmpidWilcoxon <- paste("delta","(",fl[2], "-", fl[1], ")", sep = "")


Wilcoxon.Test=data.frame(Effect = cmpid,Estimator=pd,
Statistic=Z.wilcox,p.Value=p.wilcox,Shift=cmpidWilcoxon, Hodges.Lehmann=HL,Lower=Lower.Shift,Upper=Upper.Shift)
 result <- list(Info = data.info, Analysis = Analysis, Wilcoxon=Wilcoxon.Test)

    if (plot.simci == TRUE) {
        text.Ci <- paste((1 - alpha) * 100, "%", "Confidence Interval for p")
        Lowerp <- "|"
        plot(rep(pd, plotz), 1:plotz, xlim = c(0, 1), pch = 15, 
            axes = FALSE, xlab = "", ylab = "")
        points(Lower, 1:plotz, pch = Lowerp, font = 2, cex = 2)
        points(Upper, 1:plotz, pch = Lowerp, font = 2, cex = 2)
        abline(v = 0.5, lty = 3, lwd = 2)
        for (ss in 1:plotz) {
            polygon(x = c(Lower[ss], Upper[ss]), y = c(ss, ss), 
                lwd = 2)
        }
        axis(1, at = seq(0, 1, 0.1))
        axis(2, at = 1:plotz, labels = cmpid, font = 2)
        box()
        title(main = c(text.Ci, paste("Method:", AsyMethod)))
    }
    if (info == TRUE) {
        cat("\n", "#------Nonparametric Test Procedures and Confidence Intervals for relative  effects-----#", 
            "\n", "\n", "-", "Alternative Hypothesis: ", text.Output, 
            "\n", "-", "Confidence level:", (1 - alpha) * 100, 
            "%", "\n", "-", "Method", "=", AsyMethod, "\n", "\n", 
            "#---------------------------Interpretation----------------------------------#", 
            "\n", "p(a,b)", ">", "1/2", ":", "b tends to be larger than a", 
            "\n", "#---------------------------------------------------------------------------#", 
            "\n", "\n")
    }
    #result$input <- input.list
    #result$text.Output <- text.Output
    #result$cmpid <- cmpid
    #result$AsyMethod <- AsyMethod
    #class(result) <- "ranktwosamples"
    return(result)
}
