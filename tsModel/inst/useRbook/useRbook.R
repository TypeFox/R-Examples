###############################################################################
## Fit seasonal models
## Copyright (C) 2004, Roger D. Peng <rpeng@jhsph.edu>
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
###############################################################################

## Fit a single city seasonal model using glm() and ns()

multiDFFit <- function(dfVec, city, ...) {
        if(!require(NMMAPSlite))
                stop("need 'NMMAPSlite' package")
        results <- vector("list", length = length(dfVec))
        dataframe <- readCity(city)

        results <- lapply(dfVec, function(d) {
                cat(d, "\n")
                try( fitCitySeason(dataframe, dfyr.Time = d[1], ...) )
        })
        cat("\n")
        invisible(results)
}

## Fit a single city seasonal model using glm() and ns()

## Note: Just running `fitCitySeason(city)', where `city' is the data
## frame for a particular city, should fit the usual NMMAPS model,

fitCitySeason <- function(data, pollutant = "l1pm10tmean", cause = "death",
                          season = c("none", "periodic", "factor2"),
                          tempModel = c("default", "rm7", "tempInt", "SeasonInt"),
                          dfyr.Time = 7, pdfyr.time = 0.15,
                          df.Temp = 6, df.Dew = 3,
                          df.Season = 1,  ## For "periodic" fits
                          obsThreshold = 50, extractors = NULL) {
        season <- match.arg(season)
        tempModel <- match.arg(tempModel)

        if(inherits(pollutant, "formula"))
                pollutant <- attr(terms(pollutant), "term.labels")

        ## Create some seasonal variables
        data <- setupSeason(data)
        
        ## Some exclusions
        is.na(data[, cause]) <- as.logical(data[, paste("mark", cause, sep = "")])

        ## Coerce to factor; very important!
        data$dow <- as.factor(data$dow)
        data$agecat <- as.factor(data$agecat)
        nAges <- length(levels(data[, "agecat"]))
        
        ## Specify/setup temperature variables
        temp.info <- setupTemp(data, df.Temp, tempModel)
        data <- temp.info$adj.data
        temp.f <- temp.info$temp.f
        
        ## Need to modify degrees of freedom based on missingness of data
        subform <- paste("~", paste(c("time", "time2", "time3", "agecat",
                                      "tmpd", "rmtmpd", "dptp", "rmdptp",
                                      temp.info$addedVars, cause,
                                      paste(pollutant, collapse = "+")),
                                    collapse = "+"))
        sub <- model.frame(as.formula(subform), data = data, na.action = na.pass)
        subset <- complete.cases(sub)
        df.Time <- numdf(subset, dfyr.Time)
        df.time <- df.Time * pdfyr.time
        data$subset <- subset
        nobs <- sum(subset) / nAges
        if(nobs < obsThreshold) 
                stop("Not enough observations: ", nobs, "\n")

        ## Set up smooth functions of time (with separate functions for
        ## the older two age categories
        smoothTime.info <- setupSmoothTime(data, df.Time, df.time)
        data <- smoothTime.info$adj.data
        time.f <- smoothTime.info$time.f

        ## Setup other formulas
        covariates.f <- paste(cause, "~ dow + agecat")
        weather.f <- paste(c(paste("ns(dptp,", df.Dew, ")"),
                             paste("ns(rmdptp,", df.Dew, ")")), 
                           collapse = " + ")

        if(season == "none") {
                poll.f <- paste(pollutant, collapse = " + ")
                form.str <- paste(c(covariates.f, time.f, temp.f,
                                    weather.f, "Season", poll.f), collapse = " + ")
        }
        else {
                pollSeason.f <-
                        switch(season,
                               periodic = {
                                       nfreq <- df.Season
                                       paste("harmonic(SeasonTime,", nfreq,
                                             ",365, intercept = TRUE):", pollutant,
                                             sep = "")
                               },
                               periodic2 = {
                                       nfreq <- df.Season
                                       paste("periodicBasis2(SeasonTime,", nfreq,
                                             ",365, intercept = TRUE):", pollutant,
                                             sep = "")
                               },
                               factor2 = paste("Season + Season", pollutant,
                               sep = ":"),
                               )
                form.str <- paste(covariates.f, time.f, temp.f,
                                  weather.f, pollSeason.f, sep = " + ")
        }
        form <- as.formula(form.str)

        ## Fit the model!
        fit <- glm(form, family = quasipoisson, subset = subset, data = data,
                   control = glm.control(epsilon = 1e-10, maxit = 1000),
                   na.action = na.exclude)

        ## Extract information from the fitted glm model object using the
        ## list of functions in `extractors'.  If no extractors are
        ## specified, just return the entire fitted model object.
        if(is.null(extractors))
                rval <- fit
        else 
                rval <- lapply(extractors, function(f) f(fit))
        invisible(rval)
}


periodicBasis2 <- function(x, nfreq, period, intercept = FALSE, ortho = TRUE) {
        pi <- base::pi  ## Just in case someone has redefined pi!
        stopifnot(nfreq > 0)
        x <- as.numeric(x)
        nax <- is.na(x)
        N <- seq(0, nfreq - 1)
        k <- 2^N * 2 * pi / period
        M <- outer(x, k)
        sinM <- apply(M, 2, sin)
        cosM <- apply(M, 2, cos)
        if(!intercept) { 
                basis <- cbind(sinM, cosM)

                if(ortho) 
                        basis <- qr.Q(qr(basis))
        }
        else {
                basis <- cbind(1, sinM, cosM)

                if(ortho) {
                        basis <- qr.Q(qr(basis))
                        basis[,1] <- rep.int(1, nrow(basis))
                }
        }
        basis
}

setupTemp <- function(dataframe, df.Temp, tempModel) {
        default.temp.f <- paste(c(paste("ns(tmpd,", df.Temp, ")"),
                                  paste("ns(rmtmpd,", df.Temp, ")")),
                                collapse = " + ")
        orig.namelist <- names(dataframe)
        
        if(tempModel == "default") 
                temp.f <- default.temp.f
        else if(tempModel == "rm7") {
                ## Create a running 7 day mean of tmpd
                tmpd <- dataframe[, "tmpd"]
                rm7tmpd <- filter(tmpd, filter = c(0, rep(1/7, 7)),
                                  sides = 1, circular = FALSE)
                dataframe[, "rm7tmpd"] <- as.numeric(rm7tmpd)

                ## Using 3 df here is arbitrary
                temp.f <- paste(default.temp.f, "ns(rm7tmpd, 3)", sep = " + ")
        }
        else if(tempModel == "tempInt") 
                temp.f <- paste("ns(tmpd,", df.Temp, "):ns(rmtmpd,", df.Temp, ")")
        else if(tempModel == "SeasonInt") 
                temp.f <- paste("(", default.temp.f, "):Season")
        else
                stop("Wrong temperature model specification")
        list(adj.data = dataframe, temp.f = temp.f,
             addedVars = setdiff(names(dataframe), orig.namelist))
}

## Return a modified data frame with two extra variables.  `SeasonTime'
## is the day number within each year.  `Season' is a factor
## indicating the season for each particular day.

setupSeason <- function(dataframe) {
        n <- as.numeric(table(dataframe$agecat)[1])
        year <- floor(dataframe[1:n, "date"] / 10000)
        nyeardays <- as.numeric(table(year))
        seas.time <- unlist(lapply(paste(1, nyeardays, sep = ":"),
                                   function(x) eval(parse(text = x))))
        SeasonTime <- rep(seas.time, 3)

        indlist <- genSeasonInd(dataframe)
        seas.ind <- with(indlist,
                         factor(winter * 1 + spring * 2 + summer * 3 + fall * 4,
                                labels = c("Winter", "Spring", "Summer", "Fall")))
        Season <- rep(seas.ind, 3)

        data.frame(dataframe, SeasonTime = SeasonTime, Season = Season)
}

## Re-create season indicators (taken from Aidan's `dates.R' file)

genSeasonInd <- function(dataframe) {
        n <- as.numeric(table(dataframe$agecat)[1])
        dates <- dataframe$date[1:n]
        year  <- floor(dates/10000)
        mhyr  <- dates %% 10000
        month <- floor(mhyr/100)
        day   <- mhyr %% 100
        rm(mhyr)
        
        winter <-  (( month == 12 & day > 21 ) | ( month ==  1)
                    | ( month ==  2 ) | ( month ==  3 & day <= 21))
        spring <- (( month ==  3 & day > 21 ) | ( month ==  4)
                   | ( month ==  5 ) | ( month ==  6 & day <= 21))
        summer <- (( month ==  6 & day > 21 ) | ( month ==  7)
                   | ( month ==  8 ) | ( month ==  9 & day <= 21))
        fall   <- (( month ==  9 & day > 21 ) | ( month == 10)
                   | ( month == 11 ) | ( month == 12 & day <= 21))
        list(winter = winter, spring = spring, summer = summer, fall = fall,
             all = rep(TRUE, length(winter)))
}

## This function modifies the original data frame so that the smooth
## functions of time for the older two age groups are setup correctly.
## The correct formula is also constructed.  This code was extracted
## from Aidan's original fitmodel2() function.

setupSmoothTime <- function(dataframe, df.Time, df.time) {
        df.Time <- round(df.Time)
        df.time <- round(df.time)
        subset <- dataframe$subset

        if(is.null(subset))
                stop(sQuote("dataframe"), " missing ", sQuote("subset"), " variable")

        Time <- dataframe[subset, "time"]

        if ( round(df.Time) == 1 | round(df.Time) == 0) {
                TIME <- Time
                Time <- matrix(NA, ncol = 1, nrow = nrow(dataframe))
                Time[subset, ] <- TIME
                dimnames(Time) <- list(NULL, paste("Time.", 1, sep = ""))
        }
        else {
                knots  <- quantile(Time, prob = (1:(round(df.Time)-1))/round(df.Time))
                TIME   <- ns(Time, knots = knots)
                Time   <- matrix(NA, ncol = ncol(TIME), nrow = nrow(dataframe))
                Time[subset, ] <- TIME
                dimnames(Time) <- list(NULL, paste("Time.", 1:ncol(TIME), sep = ""))
        }
        Age    <- dataframe[, "agecat"]
        Time2  <- dataframe[subset & (Age == levels(Age)[2]), "time"]

        if ( round(df.time) == 1 | round(df.time) == 0) {
                TIME2 <- Time2
                Time2  <- matrix(NA,ncol=1,nrow=nrow(dataframe))
                Time2[subset & Age == levels(Age)[2],]  <- TIME2
                Time2[subset & Age == levels(Age)[1],]  <- 0
                Time2[subset & Age == levels(Age)[3],]  <- 0
                dimnames(Time2) <- list(NULL,paste("Time2.",1,sep=""))
        }
        else {
                knots2 <- quantile(Time2, prob= (1:(round(df.time)-1))/round(df.time) )
                TIME2  <- ns(Time2, knots = knots2)
                Time2  <- matrix(NA,ncol = ncol(TIME2), nrow = nrow(dataframe))
                Time2[subset & Age == levels(Age)[2],]  <- TIME2
                Time2[subset & Age == levels(Age)[1],]  <- 0
                Time2[subset & Age == levels(Age)[3],]  <- 0
                dimnames(Time2) <- list(NULL, paste("Time2.", 1:ncol(TIME2), sep = ""))
        }
        Time3  <- dataframe[subset & (Age == levels(Age)[3]), "time"]
        
        if ( round(df.time) == 1 | round(df.time) == 0 ) {
                TIME3 <- Time3
                Time3  <- matrix(NA,ncol=1,nrow=nrow(dataframe))
                Time3[subset & Age== levels(Age)[3],]  <- TIME3
                Time3[subset & Age== levels(Age)[1],]  <- 0
                Time3[subset & Age== levels(Age)[2],]  <- 0
                dimnames(Time3) <- list(NULL,paste("Time3.",1,sep=""))
        }
        else {
                knots3 <- quantile(Time3, prob= (1:(round(df.time)-1))/round(df.time) )
                TIME3  <- ns(Time3,knots=knots3)
                Time3  <- matrix(NA,ncol=ncol(TIME3),nrow=nrow(dataframe))
                Time3[subset & Age== levels(Age)[3],]  <- TIME3
                Time3[subset & Age== levels(Age)[1],]  <- 0
                Time3[subset & Age== levels(Age)[2],]  <- 0
                dimnames(Time3) <- list(NULL,paste("Time3.",1:ncol(TIME3),sep=""))
        }
        time.f <- paste(paste("Time.", 1:ncol(Time), sep="",collapse="+"),
                        "+", paste("Time2.",1:ncol(Time2),sep="",collapse="+"),
                        "+", paste("Time3.",1:ncol(Time3),sep="",collapse="+"))
        
        list(adj.data = cbind(dataframe, Time, Time2, Time3),
             time.f = time.f)
}


############################################################################


## -----------------------------------------------------------------------
## Some support functions used in model fitting
## -----------------------------------------------------------------------

##--------------------------------------------
numdf <- function(usedata,num=5){
        n <- length(usedata)
        use <- usedata[1:(n/3)]        
        ll  <- round(length(use)/12)

        ## this is to eliminate the warning message the length of use is
        ## not a multiple of 12
        usenew <- use[1:(ll*12)]
        
        mat <- matrix(usenew, ncol=12,byrow=T)
        m   <- sum(ceiling(apply(mat,1,sum,na.rm=T)/12)) ##-365.25/12   
        df  <- round(12*m/365.25*num)   
        max(1,df)
}
##--------------------------------------------


################################################################################
################################################################################

## beta: vector of length N, or N x p matrix.
## cov: a vector of length N or a p x p x N array.

poolCoef <- function(b, cov = NULL, w = NULL, 
                     method = c("tlnise", "fixed"),
                     extractors = NULL, ...) {
        method <- match.arg(method)

        if(is.null(cov) && !is.list(b))
                stop("'b' must be list if 'cov' is 'NULL'")
        if(is.list(b)) {
                cov <- b$cov
                b <- b$beta
        }
        if(is.null(w))
                w <- rep(1, NROW(b))
        if(method == "tlnise") {
                stopifnot(exists("tlnise"))
                fit <- tlnise(b, cov, w, prnt = FALSE, ...)
                if(is.null(extractors))
                        pooled <- drop(fit$gamma[, c("est", "se")])
        }
        else if(method == "fixed") {
                stopifnot(NCOL(b) == 1)
                fit <- lm(b ~ w, weights = 1 / cov)
                if(is.null(extractors))
                        pooled <- summary(fit)$coefficients[1, 1:2]
        }   
        else
                stop("invalid 'method' specified")
        if(!is.null(extractors)) 
                pooled <- lapply(extractors, function(f) f(fit))
        pooled
}

coefSeasonal <- function(results, pollutant, method = "factor2",
                         Seasons = c("Winter", "Spring", "Summer", "Fall")) {
        if(method == "factor2") {
                betacov <- extractBetaCov(results, pollutant)
                coefmat <- poolCoef(betacov$beta, betacov$cov)
                rownames(coefmat) <- Seasons
                colnames(coefmat) <- c("Estimate", "Std. Error")
        }
        else
                stop("invalid 'method' specified")
        coefmat
}


## Return a n x 2 matrix where n is the number of seasons.  The first
## column is the pooled beta coefficient and the second column is the
## standard error.  `results' is a list of length equal to the number
## of cities.
extractBetaCov <- function(results, pollutant) {
        use <- !sapply(results, inherits, what = "try-error")
        results <- results[use]
        covlist <- lapply(results, function(x) {
                stopifnot(!is.null(x$cov))
                i <- grep(pollutant, rownames(x$cov))
                if(length(i) == 0)
                        stop(sQuote("grep"), " failed: ", sQuote(pollutant),
                             " did not match any row names")
                x$cov[i, i]
        })
        cov <- array(unlist(covlist), c(dim(covlist[[1]]), length(covlist)))

        beta <- t(sapply(results, function(x) {
                stopifnot(!is.null(x$coefficients))
                i <- grep(pollutant, rownames(x$coefficients))
                if(length(i) == 0)
                        stop(sQuote("grep"), " failed: ", sQuote(pollutant),
                             " did not match any coefficient names")
                x$coefficients[i, 1]
        }))
        if(NCOL(cov) == 1)
                cov <- as.vector(cov)
        list(beta = drop(beta), cov = cov, lcov = covlist)
}


LouisFormat <- function(x, type = c("stderr", "confint"), digits = 2) {
        type <- match.arg(type)

        if(!is.matrix(x))
                x <- as.matrix(x)

        switch(type,
               stderr = {
                       stopifnot(ncol(x) == 2)
                       minus <- ifelse(x[,1] < 0, "-", "")
                       
                       r <- paste("$",
                                  formatC(abs(x[,1]), digits = digits, format = "f"),
                                  "_", "{(",
                                  formatC(x[,2], digits = digits, format = "f"),
                                  ")}", "$", sep = "")
                       paste(minus, r, sep = "")
               },
               confint = {
                       stopifnot(ncol(x) == 3)  ## est, lo, hi
                       minus <- ifelse(x[,1] < 0, "-", "")
                       r <- paste("$", formatC(abs(x[,1]), digits = digits,
                                               format = "f"),
                                  "_", "{(", formatC(x[,2], digits = digits,
                                                     format = "f"),
                                  ",\\,", formatC(x[,3], digits = digits,
                                                  format = "f"),
                                  ")}", "$", sep = "")
                       paste(minus, r, sep = "")
               })
}



################################################################################
################################################################################

gibbspool <- function(b, v, scale = 1, maxiter = 2000, burn = 500,
                      pScale = 1e-5) {
        A <- diag(1, ncol(b)) * scale^2
        df0 <- 1
        S0 <- diag(pScale, ncol(b))
        mu <- rep(0, ncol(b))
        Sigma <- rowMeans(v, dims = 2)
        I <- diag(1, ncol(b))
        n <- nrow(b)
        if(maxiter < burn)
                burn <- 1
        
        results <- vector("list", length = maxiter)
        
        for(i in seq_len(maxiter)) {
                ## Sample beta
                beta <- lapply(seq_len(n), function(j) {
                        B <- Sigma %*% solve(v[,,j] + Sigma)
                        mvrnorm(1, mu + B %*% (b[j, ] - mu), (I - B) %*% Sigma)
                })
                beta <- do.call("rbind", beta)

                ## Sample mu
                B <- A %*% solve(Sigma / n + A)
                mu <- mvrnorm(1, B %*% colMeans(beta), (I - B) %*% A)

                ## Sample Sigma
                S1 <- 0

                for(j in seq_len(n)) 
                        S1 <- S1 + (beta[j, ] - mu) %o% (beta[j, ] - mu)
                Sigma <- riwish(n + df0, solve(S0 + S1))

                results[[i]] <- list(mu = mu, Sigma = Sigma)
        }
        if(burn > 0)
                results <- results[-seq_len(burn)]
        mu.mcmc <- t(sapply(results, "[[", "mu"))

        m <- cbind(est = colMeans(mu.mcmc),
                   se = apply(mu.mcmc, 2, sd))
        list(gamma = m, results = results)
}

## This function was taken from the MCMCpack package.

riwish <- function(v, S) {
        w <- rwish(v, S)
        solve(w)
}

rwish <- function (v, S){
        if (!is.matrix(S)) 
                S <- matrix(S)
        if (nrow(S) != ncol(S)) {
                stop("S not square in rwish()")
        }
        if (v < nrow(S)) {
                stop("v is less than the dimension of S in rwish()")
        }
        p <- nrow(S)
        CC <- chol(S)
        Z <- matrix(0, p, p)
        diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
        if (p > 1) {
                pseq <- 1:(p - 1)
                idx <- rep(p * pseq, pseq) + unlist(lapply(pseq, seq))
                Z[idx] <- rnorm(p * (p - 1)/2)
        }
        crossprod(Z %*% CC)
}


################################################################################
## MOM pooling

pooling <- function(estimate, var) {
    mu <- weighted.mean(estimate, 1 / var)
    nv <- max(var(estimate) - mean(var), 0)
    std <- sqrt(1 / sum(1 / (var + nv)))
    structure(list(mu = mu, std = std, 
                   df = length(estimate)-2, 
                   het = sqrt(nv)),
              class = "poolingResult")
}

print.poolingResult <- function(x, ...) {
    m <- with(x, matrix(c(mu, std, mu/std, 
                          pt(abs(mu/std),  df,
                             lower.tail = FALSE)),
                        byrow = TRUE, ncol = 4))
    colnames(m) <- c("Estimate", "Std. Error", 
                     "t value", "Pr(>|t|)")
    rownames(m) <- c("National avg.")
    printCoefmat(m)
    invisible(x)             
}
