logistic <- function(x) exp(x)/(1+exp(x))

#######################################
###### Create ictregj function  #######
####### FOR ONE STEP ESTIMATOR ########
#######################################

ictreg.joint <- function(formula, data = parent.frame(), treat = "treat", J, outcome = "outcome", outcome.reg = "logistic",
                             constrained = FALSE, 
                             maxIter = 1000) {

    ictreg.joint.call <- match.call()
    mf <- match.call(expand.dots = FALSE)

    ## make all other call elements null in mf <- NULL in next line
    mf$maxIter <- mf$J <- mf$treat <- mf$constrained <- mf$outcome <- mf$outcome.reg <- mf$bayes <- mf$priorscale <- NULL
    mf[[1]] <- as.name("model.frame")
    mf$na.action <- 'na.pass'
    mf <- eval.parent(mf)
    
    ## set up the data for all points 
    x.all <- model.matrix(attr(mf, "terms"), data = mf)
    y.all <- model.response(mf)

    ## list-wise missing deletion
    na.x <- apply(is.na(x.all), 1, sum)
    na.y <- is.na(y.all)
    na.o <- is.na(data[,paste(outcome)])

    t <- data[na.x == 0 & na.y == 0 & na.o == 0, paste(treat)]
    o <- data[na.x == 0 & na.y == 0 & na.o == 0, paste(outcome)] 
    if(class(t) == "factor") {
        
        levels(t) <- tolower(levels(t))
        
        if (length(which(levels(t) == "control")) == 1) {
            t <- relevel(t, ref = "control")
        } else {
            warning("Note: using the first level as the control condition, but it is not labeled 'control'.")
        }
        
        condition.labels <- levels(t)
        t <- as.numeric(t) - 1
        treatment.labels <- condition.labels[2:length(condition.labels)]
        control.label <- condition.labels[1]
        
    } else {
        condition.labels <- sort(unique(t))
        treatment.labels <- condition.labels[condition.labels != 0]
        control.label <- 0
    }
    
    ## list wise delete
    y.all <- y.all[na.x == 0 & na.y == 0 & na.o == 0]
    x.all <- x.all[na.x == 0 & na.y == 0 & na.o == 0, , drop = FALSE]
    
    ## so that the output data has the same dimension as x.all and y.all
    data.all <- data[na.x == 0 & na.y == 0 & na.o == 0, , drop = FALSE]

    data.treatment <- model.frame(formula = formula, data = subset(data.all, t == 1))
    data.control <- model.frame(formula = formula, data = subset(data.all, t == 0))

     # set up data objects for y and x for each group
    x.treatment <- model.matrix(attr(data.treatment, "terms"), data = data.treatment)
    y.treatment <- model.response(data.treatment)
    x.control <- model.matrix(attr(data.control, "terms"), data = data.control)
    y.control <- model.response(data.control)

    ## starting values
    coef.z.start <- 0
    if(constrained == FALSE){ ## unconstrained model
        start.fit.control <- lm(y.control ~ x.control - 1, data = data.control)
        start.fit.treat <- lm(y.treatment ~ x.treatment - 1, data = data.treatment)
        if(outcome.reg == "logistic"){
            start.fit.outcome <- glm(o ~ x.all + y.all - 1,
                                     family = binomial(logit), data = data.all)
            coef.outcome.draws <- rmvnorm(n = 1, coef(start.fit.outcome),
                                          vcov(start.fit.outcome))
            coef.outcome.start <- c(coef.outcome.draws, coef.z.start)
        }
        else if(outcome.reg == "linear"){
            start.fit.outcome <- lm(o ~ x.all + y.all - 1,
                                    data = data.all)
            coef.outcome.draws <- rmvnorm(n = 1, coef(start.fit.outcome),
                                          vcov(start.fit.outcome))
            coef.outcome.start <- c(coef.outcome.draws, coef.z.start)            
        }
        ## draw starting values
        coef.control.draws <- rmvnorm(n = 1, start.fit.control$coefficients,
                                      vcov(start.fit.control))
        coef.control.start <- c(coef.control.draws, coef.z.start)
        coef.treat.start <- rmvnorm(n = 1, start.fit.treat$coefficients,
                                    vcov(start.fit.treat))
        par <- c(coef.control.start, coef.treat.start, coef.outcome.start)
        if(outcome.reg == "linear"){
            sigma.start <- sigma <- summary(start.fit.outcome)$sigma
            par <- c(coef.control.start, coef.treat.start, coef.outcome.start, sigma.start)
        }
    } ## end unconstrained model
        else { # constrained model
        start.fit.control <- lm(y.control ~ x.control - 1, data = data.control)
        start.fit.treat <- lm(y.treatment ~ x.treatment - 1, data = data.treatment)
        if(outcome.reg == "logistic"){
            start.fit.outcome <- glm(o ~ x.all - 1 , family = binomial(logit),
                                     data = data.all)
            #print(start.fit.outcome)
            coef.outcome.draws <- rmvnorm(n = 1, start.fit.outcome$coefficients,
                                          vcov(start.fit.outcome))
            coef.outcome.start <- c(coef.outcome.draws, coef.z.start)
        }
        else if(outcome.reg == "linear"){
            start.fit.outcome <- lm(o ~ x.all - 1,
                                     data = data.all)
            coef.outcome.draws <- rmvnorm(n = 1, start.fit.outcome$coefficients,
                                          vcov(start.fit.outcome))
            coef.outcome.start <- c(coef.outcome.draws, coef.z.start)            
        }
        ## draw starting values
        coef.control.draws <- rmvnorm(n = 1, start.fit.control$coefficients,
                                      vcov(start.fit.control))
        coef.control.start <- c(coef.control.draws)
        coef.treat.start <- rmvnorm(n = 1, start.fit.treat$coefficients,
                                    vcov(start.fit.treat))
        par <- c(coef.control.start, coef.treat.start, coef.outcome.start)
        if(outcome.reg == "linear"){
            sigma.start <- sigma <- summary(start.fit.outcome)$sigma
            par <- c(coef.control.start, coef.treat.start, coef.outcome.start, sigma.start)
        }
    }

    obs.llik.binom.std <- function(par, J, y, treat, x, o, constrained = FALSE,
                                   outcome.reg = "logistic") {
        k <- ncol(x) ## number of covariates
        if (constrained == FALSE) { ## unconstrained model
            coef.h <- par[1:(k+1)]
            coef.g <- par[(k+2):(2*k+1)]
            x.h0 <- cbind(x,0)
            x.h1 <- cbind(x,1)
            hX0 <- logistic(x.h0 %*% coef.h) ## z = 0
            hX1 <- logistic(x.h1 %*% coef.h) ## z = 1
            gX <- logistic(x %*% coef.g)
            if (outcome.reg == "logistic") {
                coef.f <- par[(2*k + 2):(3*k + 3)]
                x.fy0 <- cbind(x, as.matrix(y), 0) #where z = 0
                x.fy1 <- cbind(x, as.matrix(y - treat), 1)
                fXy1 <- logistic(x.fy1 %*% coef.f)
                fXy0 <- logistic(x.fy0 %*% coef.f)
            } else if(outcome.reg == "linear") {
                coef.f <- par[(2*k + 2):(3*k + 3)]
                sigma <- par[3*k + 4]
                x.fy0 <- cbind(x, as.matrix(y), 0) #where z = 0
                x.fy1 <- cbind(x, as.matrix(y - treat), 1)
                fXy1 <- x.fy1 %*% coef.f
                fXy0 <- x.fy0 %*% coef.f                
            }
        } else{ ## constrained model
            coef.h <- par[1:(k)]
            coef.g <- par[(k+1):(2*k)]
            x.h0 <- cbind(x)
            x.h1 <- cbind(x)
            hX0 <- logistic(x.h0 %*% coef.h)
            hX1 <- logistic(x.h1 %*% coef.h)
            gX <- logistic(x %*% coef.g)
            if(outcome.reg == "logistic") {
                coef.f <- par[(2*k + 1):(3*k + 1)]
                x.fy0 <- cbind(x, 0) #where z = 0
                x.fy1 <- cbind(x, 1)
                fXy1 <- logistic(x.fy1 %*% coef.f)
                fXy0 <- logistic(x.fy0 %*% coef.f)
            } else if(outcome.reg == "linear") {
                coef.f <- par[(2*k + 1):(3*k + 1)]
                sigma <- par[3*k + 2]
                x.fy0 <- cbind(x, 0) #where z = 0
                x.fy1 <- cbind(x, 1)
                fXy1 <- x.fy1 %*% coef.f
                fXy0 <- x.fy0 %*% coef.f
            }
        }
        if(outcome.reg == "linear") {
            ind10 <- ((treat == 1) & (y == 0))
            ind1J1 <- ((treat == 1) & (y == (J+1)))
            ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
            nottreat <- (treat==0)
            if (sum(ind10) > 0) {
                p10 <- sum(log(1-gX[ind10])
                           + dbinom(x = y[ind10], size = J, prob = hX0[ind10], log = TRUE)
                           + dnorm(o[ind10], mean = fXy0[ind10], sd = sigma, log = TRUE))
            } else {
                p10 <- 0
            }
            if (sum(ind1J1) > 0) {
                p1J1 <- sum(log(gX[ind1J1])
                        + dbinom(y[ind1J1] - 1, size = J, prob = hX1[ind1J1], log = TRUE)
                            + dnorm(o[ind1J1], mean = fXy1[ind1J1], sd = sigma, log = TRUE))
            } else {
                p1J1 <- 0
            }
            if (sum(ind1y) > 0) {
                p1y <- sum(log(exp(log(gX[ind1y]) + dbinom(y[ind1y]-1, size = J, prob = hX1[ind1y],log = TRUE)
                                  + dnorm(o[ind1y], mean = fXy1[ind1y], sd = sigma, log = TRUE))
                               + exp(log(1-gX[ind1y]) + dbinom(y[ind1y], size = J, prob = hX0[ind1y], log = TRUE)
                                     + dnorm(o[ind1y], mean = fXy0[ind1y], sd = sigma, log = TRUE))))
            } else {
                p1y <- 0
            }
            if (sum(treat == 0) > 0) {
                p0y <- sum(log(exp(log(gX[nottreat]) + dbinom(y[nottreat], size = J, prob = hX1[nottreat], log = TRUE)
                                   + dnorm(o[nottreat], mean = fXy1[nottreat], sd = sigma, log = TRUE)) 
                               + exp(log(1-gX[nottreat]) + dbinom(y[nottreat], size = J, prob = hX0[nottreat], log = TRUE)
                                     + dnorm(o[nottreat], mean = fXy0[nottreat], sd = sigma, log = TRUE) ) )) ## control group
            } else {
                p0y <- 0
            }
        }
        else if(outcome.reg == "logistic") {
            ind10 <- ((treat == 1) & (y == 0))
            ind1J1 <- ((treat == 1) & (y == (J+1)))
            ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
            nottreat <- (treat==0)
            if (sum(ind10) > 0) {
                p10 <- sum(log(1-gX[ind10])
                           + dbinom(x = y[ind10], size = J, prob = hX0[ind10], log = TRUE)
                           + o[ind10] * log(fXy0[ind10]) + (1 - o[ind10]) * log(1-fXy0[ind10]))
            } else {
                p10 <- 0
            }
            if (sum(ind1J1) > 0) {
                p1J1 <- sum(log(gX[ind1J1])
                            + dbinom(y[ind1J1] - 1, size = J, prob = hX1[ind1J1], log = TRUE)
                            + o[ind1J1] * log(fXy1[ind1J1]) + (1 - o[ind1J1]) * log(1 - fXy1[ind1J1]))
            } else {
                p1J1 <- 0
            }
            if (sum(ind1y) > 0) {
                p1y <- sum(log(exp(log(gX[ind1y]) + dbinom(y[ind1y]-1, size = J, prob = hX1[ind1y],log = TRUE)
                                   + o[ind1y]*log(fXy1[ind1y]) + (1-o[ind1y])*log(1-fXy1[ind1y]))
                               + exp(log(1-gX[ind1y]) + dbinom(y[ind1y], size = J, prob = hX0[ind1y], log = TRUE)
                                     + o[ind1y]*log(fXy0[ind1y]) + (1-o[ind1y])*log(1-fXy0[ind1y]))))
            } else {
                p1y <- 0
            }
            if (sum(treat == 0) > 0) {
                p0y <- sum(log(exp(log(gX[nottreat]) + dbinom(y[nottreat], size = J, prob = hX1[nottreat], log = TRUE)
                                   + o[nottreat]*log(fXy1[nottreat]) + (1-o[nottreat])*log(1-fXy1[nottreat])) 
                               + exp(log(1-gX[!treat]) + dbinom(y[nottreat], size = J, prob = hX0[nottreat], log = TRUE)
                                     + o[nottreat]*log(fXy0[nottreat]) + (1-o[nottreat])*log(1-fXy0[nottreat]) ) )) ## control group
            } else {
                p0y <- 0
            }
        }
        return(p10+p1J1+p1y+p0y)
    }
    
    ## Estep (weights)
    Estep.binom.std <- function(par, J, y, treat, x, o, constrained = FALSE, outcome.reg = "logistic") {
        k <- ncol(x) ## number of covariates
        if (constrained == FALSE) { ## unconstrained model
            coef.h <- par[1:(k+1)]
            coef.g <- par[(k+2):(2*k+1)]
            x.h0 <- cbind(x,0)
            x.h1 <- cbind(x,1)
            hX0 <- logistic(x.h0 %*% coef.h) ## z = 0
            hX1 <- logistic(x.h1 %*% coef.h) ## z = 1
            gX <- logistic(x %*% coef.g)
            if (outcome.reg == "logistic") {
                coef.f <- par[(2*k + 2):(3*k + 3)]
                x.fy0 <- cbind(x, as.matrix(y), 0) #where z = 0
                x.fy1 <- cbind(x, as.matrix(y - treat), 1)
                fXy1 <- logistic(x.fy1 %*% coef.f)
                fXy0 <- logistic(x.fy0 %*% coef.f)
            } else if(outcome.reg == "linear") {
                coef.f <- par[(2*k + 2):(3*k + 3)]
                sigma <- par[3*k + 4]
                x.fy0 <- cbind(x, as.matrix(y), 0) #where z = 0
                x.fy1 <- cbind(x, as.matrix(y - treat), 1)
                fXy1 <- x.fy1 %*% coef.f
                fXy0 <- x.fy0 %*% coef.f
            }
        } else { ## constrained model
            coef.h <- par[1:(k)]
            coef.g <- par[(k+1):(2*k)]
            x.h0 <- cbind(x)
            x.h1 <- cbind(x)
            hX0 <- logistic(x.h0 %*% coef.h)
            hX1 <- logistic(x.h1 %*% coef.h)
            gX <- logistic(x %*% coef.g)
            if(outcome.reg == "logistic") {
                coef.f <- par[(2*k + 1):(3*k + 1)]
                x.fy0 <- cbind(x, 0) #where z = 0
                x.fy1 <- cbind(x, 1)
                fXy1 <- logistic(x.fy1 %*% coef.f)
                fXy0 <- logistic(x.fy0 %*% coef.f)
            } else if(outcome.reg == "linear") {
                coef.f <- par[(2*k + 1):(3*k + 1)]
                sigma <- par[3*k + 2]
                x.fy0 <- cbind(x, 0) #where z = 0
                x.fy1 <- cbind(x, 1)
                fXy1 <- x.fy1 %*% coef.f
                fXy0 <- x.fy0 %*% coef.f
            }
        }
        ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
        w <- rep(NA, length(y))
        if(outcome.reg == "logistic") {        
            w[ind] <- exp(log(o[ind]*fXy1[ind] + (1-o[ind])*(1-fXy1[ind]))
                          + log(gX[ind]) + dbinom((y-treat)[ind], size = J, prob = hX1[ind], log = TRUE)
                          - log(exp(log(o[ind]*fXy1[ind] + (1-o[ind])*(1-fXy1[ind]))+log(gX[ind])
                                    + dbinom((y-treat)[ind], size = J, prob = hX1[ind], log = TRUE))
                                + exp(log(o[ind]*fXy0[ind] + (1-o[ind])*(1-fXy0[ind]))
                                      + log(1-gX[ind])+dbinom(y[ind], size = J, prob = hX0[ind], log = TRUE))))
        }
        else if(outcome.reg == "linear") {
            w[ind] <- exp(dnorm(o[ind], mean = fXy1[ind], sd = sigma, log = TRUE)
                          + log(gX[ind]) + dbinom((y-treat)[ind], size = J, prob = hX1[ind], log = TRUE) 
                          - log(exp(dnorm(o[ind], mean = fXy1[ind], sd = sigma, log = TRUE)
                                    + log(gX[ind]) + dbinom((y-treat)[ind], size = J, prob = hX1[ind], log = TRUE))
                                + exp(dnorm(o[ind], mean = fXy0[ind], sd = sigma, log = TRUE)
                                      + log(1-gX[ind]) + dbinom(y[ind], size = J, prob = hX0[ind], log = TRUE))))

        }
        w[(treat == 1) & (y == 0)] <- 0
        w[(treat == 1) & (y == (J+1))] <- 1        
        return(w)
    }
    
    ## Mstep 1: g (sensitive item)
    
    wlogit.fit.std <- function(y, treat, x, w, par = NULL) {
        yrep <- rep(c(1,0), each = length(y))
        xrep <- rbind(x, x)
        wrep <- c(w, 1-w)
        return(glm(cbind(yrep, 1-yrep) ~ xrep - 1, weights = wrep, family = binomial(logit), start = par))       
    }
  
    ## Mstep 3: h (non-sensitive)
    wbinomial.fit <- function(y, treat, x, w, par = NULL, constrained = FALSE) {
        nPar <- ncol(x)
        yrep <- c(y - treat, y)
        xrep <- rbind(x, x)
        zrep <- rep(c(1, 0), each = length(y))
        wrep <- c(w, 1-w)
        yrep <- yrep[wrep > 0]
        xrep <- xrep[wrep > 0, ]
        zrep <- zrep[wrep > 0]
        wrep <- wrep[wrep > 0]
        if(constrained == FALSE){
            par <- par[1:(nPar+1)]
            fit <- glm(cbind(yrep, J - yrep) ~ xrep + zrep - 1, weights = wrep, family = binomial(logit), start = par)
        } else {
            par <- par[1:(nPar)]
            fit <- glm(cbind(yrep, J - yrep) ~ xrep - 1, weights = wrep, family = binomial(logit), start = par)
        }
        return(fit)
    }
        
    ## Mstep 2: f (outcome)
    wlogit.fit.outcome <- function(y, treat, x, w, o, par = NULL, constrained = FALSE) {            
        yrep <- c((y - treat), y)
        xrep <- rbind(x, x)
        zrep <- rep(c(1,0), each = length(y))
        wrep <- c(w, 1-w)
        orep <- c(o,o)
        nPar <- ncol(x)
        yrep <- yrep[wrep > 0]
        orep <- orep[wrep > 0]
        xrep <- xrep[wrep > 0, ]
        zrep <- zrep[wrep > 0]
        wrep <- wrep[wrep > 0]
        if(constrained == FALSE){
            if(outcome.reg == "logistic"){
                par <- par[(2*nPar + 2):(3*nPar + 3)]
                fit <- glm(cbind(orep, 1-orep) ~ xrep + yrep + zrep - 1, weights = wrep, family = binomial(logit),
                           start = par)
            } else if(outcome.reg == "linear"){
                par <- par[(2*nPar + 2):(3*nPar + 3)]
                fit <- lm(orep ~ xrep + yrep + zrep - 1, weights = wrep)
            }
        } else { #constrained
            if(outcome.reg == "logistic"){
                par <- par[(2*nPar + 1):(3*nPar + 1)]
                fit <- glm(cbind(orep, 1-orep) ~ xrep + zrep - 1, weights = wrep, family = binomial(logit), start = par)
            } else if(outcome.reg == "linear"){
                par <- par[(2*nPar + 1):(3*nPar + 1)]
                fit <- lm(orep ~ xrep + zrep - 1, weights = wrep)
            }
        }
        return(fit)
    }
    
########################################
########## NOW THE ALGORITHM ###########

    ## start off with an infinitely negative log likelihood
    pllik.const <- -Inf
    counter <- 0
    nPar <- ncol(x.all)
    
    ## calculate the log likelihood at the starting values
    llik.const <- obs.llik.binom.std(par, J = J, y = y.all, treat = t, o = o,
                                     x = x.all, constrained = constrained,
                                     outcome.reg = outcome.reg)
    #print(llik.const)
    while (((llik.const - pllik.const) > 10^(-4)) & (counter < maxIter)) {
        
        w <- Estep.binom.std(par, J, y = y.all, treat = t, x = x.all, o = o,
                             outcome.reg = outcome.reg, constrained = constrained)
        ## treatment 
        if(constrained == FALSE) {
            part <- par[(nPar+2):(2*nPar+1)]
        }
        else {
            part <- par[(nPar+1):(2*nPar)]
        }
        lfit <- wlogit.fit.std(y.all, treat = t, x.all, w, par = part)

        ## outcome
        if(outcome.reg == "logistic" | outcome.reg == "linear"){
            outfit <- wlogit.fit.outcome(y.all, treat = t, x.all, w, o = o, par = par,
                                         constrained = constrained)
        } 

        ## control 
        fit <- wbinomial.fit(y = y.all, treat = t, x = x.all, w = w, par = par,
                             constrained = constrained)
        outfit.sum <- summary(outfit)
        sigma <- outfit.sum$sigma*sqrt(outfit.sum$df[2]/(sum(outfit$weights)-outfit.sum$df[1]))
        par <- c(coef(fit), coef(lfit), coef(outfit), sigma) #control, treat, outcome, sigma
        pllik.const <- llik.const ## make old llik the new one
        
        ## calculate the new log likelihood 
        llik.const <- obs.llik.binom.std(par, J = J, y = y.all, treat = t, o = o,
                                         x = x.all, constrained = constrained,
                                         outcome.reg = outcome.reg)
        counter <- counter + 1 ## up the counter
        #print(counter)
        
        if (llik.const < pllik.const)
            warning("log-likelihood is not monotonically increasing.")  
        if(counter == (maxIter-1))
            warning("number of iterations exceeded maximum in ML")
        
        cat("llik:", llik.const, "\n")
        cat("llik diff:", llik.const - pllik.const, "\n")
        cat("par:", par, "\n")        
    }
    
    score <- function(par, J, y, treat, x, o, constrained = TRUE, outcome.reg = "logistic") {
        k <- ncol(x)
        if (constrained == FALSE) { ## unconstrained model
            coef.h <- par[1:(k+1)]
            coef.g <- par[(k+2):(2*k+1)]
            x.h0 <- cbind(x,0)
            x.h1 <- cbind(x,1)
            coef.f <- par[(2*k + 2):(3*k + 3)]
            x.fy0 <- cbind(x, as.matrix(y), 0) #where z = 0
            x.fy1 <- cbind(x, as.matrix(y - treat), 1)
            if (outcome.reg == "linear") {
                sigma <- par[3*k + 4]
            }
        } else { ## constrained model
            k <- ncol(x)
            coef.h <- par[1:(k)]
            coef.g <- par[(k+1):(2*k)]
            x.h0 <- x.h1 <- cbind(x)
            coef.f <- par[(2*k + 1):(3*k + 1)]
            x.fy0 <- cbind(x, 0) #where z = 0
            x.fy1 <- cbind(x, 1)
            if (outcome.reg == "linear") {
                sigma <- par[3*k + 2]
            }
        }
        hX0 <- logistic(x.h0 %*% coef.h) ## z = 0
        hX1 <- logistic(x.h1 %*% coef.h) ## z = 1
        gX <- logistic(x %*% coef.g)
        Xb1 <- x.fy1 %*% coef.f #this is x*beta
        Xb0 <- x.fy0 %*% coef.f
        if (outcome.reg == "logistic") {
            fXy1 <- exp(o * Xb1) / (1 + exp(Xb1)) 
            fXy0 <- exp(o * Xb0) / (1 + exp(Xb0)) 
            fprime1 <- as.vector((2*o - 1) * (exp(Xb1) / (1 + exp(Xb1))^2)) * x.fy1
            fprime0 <- as.vector((2*o - 1) * (exp(Xb0) / (1 + exp(Xb0))^2)) * x.fy0
        } else if (outcome.reg == "linear") {
            fXy1 <- dnorm(o, mean = Xb1, sd = sigma, log = TRUE)
            fXy0 <- dnorm(o, mean = Xb0, sd = sigma, log = TRUE)
            f.partial.beta1 <- as.vector(exp(fXy1) * (1/sigma^2) * (o - Xb1)) * x.fy1 #this is right, 1000 x 4, vector multi
            f.partial.beta0 <- as.vector(exp(fXy0) * (1/sigma^2) * (o - Xb0)) * x.fy0 
            f.partial.sigma1 <- exp(fXy1) * (1/(2*sigma^2)) * (-1 + (1/sigma^2) * (o - Xb1)^2)
            f.partial.sigma0 <- exp(fXy0) * (1/(2*sigma^2)) * (-1 + (1/sigma^2) * (o - Xb0)^2) 
            fprime1 <- cbind(f.partial.beta1, f.partial.sigma1)
            fprime0 <- cbind(f.partial.beta0, f.partial.sigma0)            
        }
        
        ind10 <- ((treat == 1) & (y == 0))
        ind1J1 <- ((treat == 1) & (y == (J+1)))
        ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
        nottreat <- !treat       
        gprime <- as.vector(gX) * as.vector((1 - gX)) * x
        
        if (outcome.reg == "logistic") {
            parttheta.1 <- as.vector(1 - fXy0) * x.fy0 
            parttheta.2 <- as.vector(1 - fXy1) * x.fy1 
            parttheta.3 <- (fprime1 * dbinom(y - treat, size = J, prob = hX1, log = FALSE) * as.vector(gX) + fprime0 * dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector(1 - gX)) / (as.vector(fXy1) * dbinom(y - treat, size = J, prob = hX1, log = FALSE) * as.vector(gX) + as.vector(fXy0) * dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector(1 - gX))
            parttheta.4 <- (fprime1 * dbinom(y, size = J, prob = hX1, log = FALSE) * as.vector(gX) + fprime0 * dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector(1 - gX)) / (as.vector(fXy1) * dbinom(y, size = J, prob = hX1, log = FALSE) * as.vector(gX) + as.vector(fXy0)* dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector(1 - gX))
        } else if (outcome.reg == "linear") {
            parttheta.1 <- cbind(as.vector((1/sigma^2) * (o - Xb0)) * x.fy0, (1/(2*sigma^2)) * (-1 + (1/sigma^2) * (o - Xb0)^2))
            parttheta.2 <- cbind(as.vector((1/sigma^2) * (o - Xb1)) * x.fy1, (1/(2*sigma^2)) * (-1 + (1/sigma^2) * (o - Xb1)^2))
            parttheta.3 <- (fprime1 * exp(dbinom(y - treat, size = J, prob = hX1, log = TRUE) + log(as.vector(gX))) + fprime0 * exp(dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector(1 - gX)))) / ((exp(as.vector(fXy1) + dbinom(y - treat, size = J, prob = hX1, log = TRUE) + log(as.vector(gX)))) + exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector(1 - gX))))

            parttheta.4 <- (fprime1 * exp(dbinom(y, size = J, prob = hX1, log = TRUE) + log(as.vector(gX))) + fprime0 * exp(dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector(1 - gX)))) / (exp(as.vector(fXy1) + dbinom(y, size = J, prob = hX1, log = TRUE) + log(as.vector(gX))) + exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector(1 - gX))))
        }
        parttheta.1r <- parttheta.1 * ind10
        parttheta.2r <- parttheta.2 * ind1J1
        parttheta.3r <- parttheta.3 * ind1y
        parttheta.4r <- parttheta.4 * nottreat ## Not defined where y = 4
        
        idx <- apply(is.na(parttheta.4r), 1, all)
        if (outcome.reg == "logistic" & constrained == FALSE) {
            cero <- rep(0, ncol(x) + 2)
        } else if (outcome.reg == "logistic" & constrained == TRUE) {
            cero <- rep(0, ncol(x) + 1)
        } else if (outcome.reg == "linear" & constrained == FALSE) {
            cero <- rep(0, ncol(x) + 3)
        } else if (outcome.reg == "linear" & constrained == TRUE){
            cero <- rep(0, ncol(x) + 2)
        }
        parttheta.4r[idx,] <- cero
        
        parttheta <- parttheta.1r + parttheta.2r + parttheta.3r + parttheta.4r
        
        partdel.1 <- -as.vector(gX) * x
        partdel.1r <- partdel.1 * ind10
        
        partdel.2 <- as.vector(1 - gX) * x
        partdel.2r <- partdel.2 * ind1J1

        if (outcome.reg == "logistic") {
            partdel.3 <- (as.vector(fXy1) * dbinom(y - treat, size = J, prob = hX1, log = FALSE) * gprime - as.vector(fXy0) * dbinom(y, size = J, prob = hX0, log = FALSE) * gprime) / (as.vector(fXy1) * dbinom(y - treat, size = J, prob = hX1, log = FALSE) * as.vector(gX) + as.vector(fXy0) * dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector(1 - gX))
            partdel.4 <- (as.vector(fXy1) * dbinom(y, size = J, prob = hX1) * gprime -  as.vector(fXy0) * dbinom(y, size = J, prob = hX0, log = FALSE) * gprime) / (as.vector(fXy1) * dbinom(y, size = J, prob = hX1) * as.vector(gX) + as.vector(fXy0) * dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector(1 - gX))
        } else if (outcome.reg == "linear") {
            partdel.3 <- (exp(as.vector(fXy1) + dbinom(y - treat, size = J, prob = hX1, log = TRUE)) * gprime - exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE)) * gprime) / (exp(as.vector(fXy1) + dbinom(y - treat, size = J, prob = hX1, log = TRUE) + log(as.vector(gX))) + exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector(1 - gX))))
            partdel.4 <- (exp(as.vector(fXy1) + dbinom(y, size = J, prob = hX1, log = TRUE)) * gprime - exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE)) * gprime) / (exp(as.vector(fXy1) + dbinom(y, size = J, prob = hX1, log = TRUE) + log(as.vector(gX))) + exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector(1 - gX))))
        }
        partdel.3r <- partdel.3 * ind1y
        partdel.4r <- partdel.4 * nottreat ## Not defined where y = 4
        idx <- apply(is.na(partdel.4r), 1, all)
        cero <- rep(0, ncol(x))  
        partdel.4r[idx,] <- cero
        
        partdel <- partdel.1r + partdel.2r + partdel.3r + partdel.4r
        
        hpyx0 <- as.vector(dbinom(y, size = J, prob = hX0, log = FALSE) * (y - J * hX0)) * x.h0 #h prime y = y, z = 0
        hpy1x1 <- as.vector(dbinom(y - treat, size = J, prob = hX1, log = FALSE) * ((y - treat) - J * hX1)) * x.h1 #hprime y = y-1, z = 1
        hpy1x0 <- as.vector(dbinom(y, size = J, prob = hX1, log = FALSE) * (y - J * hX1)) * x.h1 #hprime y = y, z = 1
                                        #last two are exactly the same thing, but for treat vs. control group
        partpsi.1 <- as.vector(y - J * hX0) * x.h0
        partpsi.1r <- partpsi.1 * ind10
    
        partpsi.2 <- as.vector((y - treat) - J * hX1) * x.h1
        partpsi.2r <- partpsi.2 * ind1J1

        if (outcome.reg == "logistic") {
            partpsi.3 <- (as.vector(fXy1) * hpy1x1 * as.vector(gX) + as.vector(fXy0) * hpyx0 * as.vector(1 - gX)) / (as.vector(fXy1) * dbinom(y - treat, size = J, prob = hX1, log = FALSE) * as.vector(gX) + as.vector(fXy0) * dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector((1 - gX)))
            partpsi.4 <- (as.vector(fXy1) * hpy1x0 * as.vector(gX) + as.vector(fXy0) * hpyx0 * as.vector(1 - gX)) / (as.vector(fXy1) * dbinom(y, size = J, prob = hX1, log = FALSE) * as.vector(gX) + as.vector(fXy0) * dbinom(y, size = J, prob = hX0, log = FALSE) * as.vector(1 - gX))
        } else if (outcome.reg == "linear") {
            partpsi.3 <- (exp(as.vector(fXy1) + log(as.vector(gX))) * hpy1x1 + hpyx0 * exp(log(as.vector(1 - gX)) + as.vector(fXy0))) / (exp(as.vector(fXy1) + dbinom(y - treat, size = J, prob = hX1, log = TRUE) + log(as.vector(gX))) + exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector((1 - gX)))))
            partpsi.4 <- (exp(as.vector(fXy1) + log(as.vector(gX))) * hpy1x0 + hpyx0 * exp(log(as.vector(1 - gX)) + as.vector(fXy0))) / (exp(as.vector(fXy1) + dbinom(y, size = J, prob = hX1, log = TRUE) + log(as.vector(gX))) + exp(as.vector(fXy0) + dbinom(y, size = J, prob = hX0, log = TRUE) + log(as.vector(1 - gX))))
        }
        partpsi.3r <- partpsi.3 * ind1y
        partpsi.4r <- partpsi.4 * nottreat ## Not defined where y = 4.
        
        idx <- apply(is.na(partpsi.4r), 1, all)
        if (constrained == FALSE) {
            cero <- rep(0, ncol(x) + 1)
        } else if (constrained == TRUE) {
            cero <- rep(0, ncol(x))
        }
        partpsi.4r[idx,] <- cero
        
        partpsi <- partpsi.1r + partpsi.2r + partpsi.3r + partpsi.4r

        ssize <- nrow(x.all)
        S <- cbind(partpsi, partdel, parttheta) #control, treat, outcome 
        #S <- cbind(partdel, partpsi, parttheta) #treat, control, outcome 
        info <- (t(S) %*% S) / ssize
        vcov <- solve(info) / ssize #Hessian / n
        ses <- sqrt(diag(vcov))
        returnses <- list(ses, vcov)
        return(returnses)
    }

    if (outcome.reg == "linear") {
        sigma <- par[3*nPar + 4]
    } else {
        sigma <- NA
    }
    ses <- score(par = par, J = J, y = y.all, treat = t, o = o,
                 x = x.all, outcome.reg = outcome.reg, constrained = constrained)
    
    if(constrained == FALSE) {
        par.control <- par[1:(nPar+1)]
        ses.control <- ses[[1]][1:(nPar+1)]
    } else {
        par.control <- par[1:(nPar)]
        ses.control <- ses[[1]][1:(nPar)]
    }
    
    if(constrained == FALSE) {
        par.treat <- par[(nPar+2):(2*nPar+1)]
        ses.treat <- ses[[1]][(nPar+2):(2*nPar+1)]
    } else {
        par.treat <- par[(nPar+1):(2*nPar)]
        ses.treat <- ses[[1]][(nPar+1):(2*nPar)]
    }
    
    if(constrained == FALSE) {
        par.outcome <- par[(2*nPar + 2):(3*nPar + 3)]
        ses.outcome <- ses[[1]][(2*nPar + 2):(3*nPar + 3)]
    } else {
        par.outcome <- par[(2*nPar + 1):(3*nPar + 1)]
        ses.outcome <- ses[[1]][(2*nPar + 1):(3*nPar + 1)]
    }
    

    return.object <- list(par, par.control, par.treat, par.outcome, ses.control, ses.treat, ses.outcome, sigma,
                          counter, ses, w, match.call(), treatment.labels, control.label,
                          "ml", "standard", colnames(x.all), FALSE, constrained, FALSE, llik.const, FALSE, J,
                          ses[[2]], x.all, y.all, t, outcome.reg)
    names(return.object) <- list("par", "par.control", "par.treat", "par.outcome",
                                 "se.control", "se.treat", "se.outcome", "sigma", "no.iterations", "ses", "w",
                                 "call", "treat.labels", "control.labels", "method",
                                 "design", "coef.names", "multi", "constrained", "overdispersed", "llik", "boundary", "J",
                                 "vcov", "x", "y", "treat", "outcome.reg")
    class(return.object) <- "ictreg.joint"
    return(return.object)
}


vcov.ictreg.joint <- function(object, ...){
    vcov <- object$vcov
    nPar <- length(object$coef.names)

    if(object$constrained == FALSE) { #Unconstrained
        ## PUT in order of treat, control, outcome
        vcov <- rbind(cbind(vcov[(nPar+2):(nPar*2 + 1),(nPar+2):(nPar*2 + 1)],
                             vcov[(nPar+2):(nPar*2 + 1), 1:(nPar + 1)] ,
                             vcov[(nPar+2):(nPar*2 + 1), (nPar*2+2):(nPar*3 + 3)]),
                       cbind(vcov[1:(nPar + 1), (nPar+2):(nPar*2 + 1)],
                             vcov[1:(nPar + 1), 1:(nPar + 1)] ,
                             vcov[1:(nPar + 1), (nPar*2+2):(nPar*3 + 3)]),
                       cbind(vcov[(nPar*2 + 2):(nPar*3 + 3), (nPar+2):(nPar*2 + 1)],
                             vcov[(nPar*2 + 2):(nPar*3 + 3),1: (nPar + 1)],
                             vcov[(nPar*2 + 2):(nPar*3 + 3),(nPar*2+2):(nPar*3 + 3)]))
        rownames(vcov)[1:(nPar)] <- colnames(vcov)[1:(nPar)] <- paste("sensitive.", object$coef.names, sep = "")
        rownames(vcov)[(nPar+1):(2*nPar)] <- colnames(vcov)[(nPar+1):(2*nPar)] <- paste("control.", object$coef.names, sep = "")
        rownames(vcov)[2*nPar + 1] <- colnames(vcov)[2*nPar + 1]  <- paste("control.sensitiveitem")
        rownames(vcov)[(2*nPar + 2):(3*nPar + 1)] <- colnames(vcov)[(2*nPar + 2):(3*nPar + 1)] <- paste("outcome.", object$coef.names, sep = "")
        rownames(vcov)[3*nPar + 2] <- colnames(vcov)[3*nPar + 2] <- paste("outcome.controlitems")
        rownames(vcov)[(3*nPar + 3)] <- colnames(vcov)[(3*nPar + 3)] <- paste("outcome.sensitiveitem")
     }
    else { #constrained
        vcov <- rbind(cbind(vcov[(nPar+1):(nPar*2),(nPar+1):(nPar*2)],
                             vcov[(nPar+1):(nPar*2), 1:nPar] ,
                             vcov[(nPar+1):(nPar*2), (nPar*2+1):(nPar*3 + 1)]),
                       cbind(vcov[1:nPar, (nPar+1):(nPar*2)],
                             vcov[1:nPar,1:nPar] ,
                             vcov[1:nPar, (nPar*2+1):(nPar*3 + 1)]),
                       cbind(vcov[(nPar*2+1):(nPar*3 + 1), (nPar+1):(nPar*2)],
                             vcov[(nPar*2+1):(nPar*3 + 1),1:nPar],
                             vcov[(nPar*2+1):(nPar*3 + 1),(nPar*2+1):(nPar*3 + 1)]))
        rownames(vcov)[1:(nPar)] <- colnames(vcov)[1:(nPar)] <- paste("sensitive.", object$coef.names, sep = "")
        rownames(vcov)[(nPar+1):(2*nPar)] <- colnames(vcov)[(nPar+1):(2*nPar)] <- paste("control.", object$coef.names, sep = "")
        rownames(vcov)[(2*nPar + 1):(3*nPar)] <- colnames(vcov)[(2*nPar + 1):(3*nPar)] <- paste("outcome.", object$coef.names, sep = "")
        rownames(vcov)[(3*nPar + 1)] <- colnames(vcov)[(3*nPar + 1)] <- paste("outcome.sensitiveitem")
    }
    return(vcov)
}

summary.ictreg.joint <- function(object, ...) {
    structure(object, class = c("summary.ictreg.joint", class(object)))
}

print.summary.ictreg.joint <- function(x, ...){
    ##print.summary.ictreg(x) ## Use GB function for treatment and control models.
                            ## Only works for constrained models. Need some changes for unconstrained.

    tb.sensitive <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
    colnames(tb.sensitive) <- c("Est.", "S.E.")
    rownames(tb.sensitive) <- c(x$coef.names)
    tb.sensitive[,1] <- x$par.treat
    tb.sensitive[,2] <- x$se.treat

    tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.control))
    colnames(tb.control) <- c("Est.", "S.E.")
    if(x$constrained == TRUE)
        rownames(tb.control) <- c(x$coef.names)
    else
        rownames(tb.control) <- c(x$coef.names, "Sensitive item")
    tb.control[,1] <- x$par.control
    tb.control[,2] <- x$se.control
    
    tb.outcome <- matrix(NA, ncol = 2, nrow = length(x$par.outcome))
    colnames(tb.outcome) <- c("Est.", "S.E.")
    if(x$constrained == TRUE)
        rownames(tb.outcome) <- c(x$coef.names, "Sensitive item")
    else
        rownames(tb.outcome) <- c(x$coef.names, "y", "Sensitive item")
    tb.outcome[,1] <- x$par.outcome
    tb.outcome[,2] <- x$se.outcome

    cat("\nJoint List Experiment Outcome Regression\n\n")
    
     cat("Sensitive item \n")
    print(as.matrix(round(tb.sensitive, 5)))
    
    cat("\n\nControl items \n")
    print(as.matrix(round(tb.control, 5)))
   
    cat("\n\nOutcome \n")
    print(as.matrix(round(tb.outcome, 5)))

    cat("\n")
    
    invisible(x)
}

print.ictreg.joint <- function(x, ...){
    cat("\nItem Count Technique Regression \n\nCall: ")
    
    dput(x$call)
    
    cat("\nCoefficient estimates\n")

    print(coef.ictreg.joint(x))
    
    treat.print <- c()
    for (i in 1:length(x$treat.labels)) {
        treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
        if (i != length(x$treat.labels))
            treat.print <- c(treat.print, " and ")
    }
    
    cat("Number of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
    cat(treat.print, sep ="")
    cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
    
    invisible(x)

}

coef.ictreg.joint <- function(object, ...){
    coef <- c(object$par.treat, object$par.control, object$par.outcome)
    names(coef) <- c(paste("sensitive.", object$coef.names, sep = ""),
                     paste("control.", object$coef.names, sep = ""),
                     paste("outcome.", object$coef.names, sep = ""),
                     "outcome.sensitiveitem") ## This will need changes for unconstrained

    return(coef)
}


# Users do prediction by setting either Z = 0 or Z = 1
# Users want the difference in prediction between Z = 0 and Z = 1
# In either case, users may want to compute the average quantities across all observations in the data frame

## predict.ictreg.joint produces predicted values, obtained by evaluating the regression function in the
## frame newdata (which defaluts to model.frame(object)). By using sensitive.value, users must set
## the value of z -- the latent response to the sensitive item -- to be either zero or
## one, depending on the prediction that the user requires.

## Two additional types of mean prediction are also available. The first, if a newdata.diff data frame
## is provided by the user, calculates the mean predicted values across two datasets, as well as the
## mean difference in predicted value. Standard errors and confidence intervals can also be added.
## Users may also set the logical sensitive.diff to TRUE and sensitive.value to "both", which will output 
## the mean predicted values across all observations for z = 0 as well as z = 1, in addition to the mean 
## difference in predicted value. Standard errors and confidence intervals can also be added. For
## difference predictions (sensitive.diff and avg.diff), the option avg must be set to TRUE.

## Users can also use the predict.sensitive = TRUE option to generate predictions of responses
## to the sensitive item. This uses

predict.ictreg.joint <- function(object, newdata, newdata.diff, se.fit = FALSE,
                                 interval = c("none","confidence"), level = .95,
                                 avg = FALSE, sensitive.value = c("0", "1", "both"),
                                 sensitive.diff = FALSE, return.draws = FALSE,
                                 predict.sensitive = FALSE, ...){
    if(missing(sensitive.value))
        stop("Must set a value of 0 or 1 for the sensitive item.")
    if(missing(interval)) interval <- "none"
    nx <- ncol(object$x)
    logistic <- function(object) exp(object)/(1+exp(object))
    
    ## Get vcov and coefficients from outcome model    
    var.matrix <- object$vcov
    if(object$constrained == FALSE) { #UNconstrained
        var.matrix <- var.matrix[(2*nx + 2):(3*nx + 3),(2*nx + 2):(3*nx + 3)]
    } else { #constrained
        var.matrix <- var.matrix[(nx*2 + 1):(nx*3 +1),(nx*2 + 1):(nx*3 +1)]
    }
    coef.matrix <- object$par.outcome
    
    ## Get (new) dataframe
    if (missing(newdata)) {  
        xvar <- object$x # if no new data specified, use original 
    } else {
        if(nrow(newdata)==0)
            stop("No data in the provided data frame.")
        xvar <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata) # if newdata given, paste it as xvar. 
    }
    ## Add newdata.diff if provided
    if (!missing(newdata.diff)) {
        data.list <- list(xvar,
                          model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                                       newdata.diff)) # has two elements, newdata and newdata.diff
    } else if (missing(newdata.diff)) {
        data.list <- list(xvar)
    }
    ## If user wants predictions for z = 1 and z = 0, or they need the average 
    ## difference in predictions between those two, then data.list should have
    ## two elements, first where z = 0 and second where z = 1
    if (sensitive.value == "both"){
        data.list <- list(xvar, xvar)
    }
    
    
    return.object <- c()
    
    ## Before adding z and/or y to data.list, generate predictions from
    ## sensitive item model.
    if (predict.sensitive == TRUE) {
        var.matrix.sens <- vcov(object)
        var.matrix.sens <- var.matrix.sens[(1):(nx),(1):(nx)]
        coef.matrix.sens <- object$par.treat
        draws.predict.sens <- draws.mean.sens <- list()
        mean.sens <- ci.lower.sens <- ci.upper.sens <- sd.sens <- c()
        for (i in 1: length(data.list)){
            draws.sens <- mvrnorm(n = 10000, coef.matrix.sens, var.matrix.sens)
            draws.predict.sens[[i]] <- logistic(as.matrix(data.list[[i]])%*% t(draws.sens))  
            draws.mean.sens[[i]] <- apply(draws.predict.sens[[i]], 2, mean)
            mean.sens[i] <- mean(draws.mean.sens[[i]])
            sd.sens[i] <- sd(draws.mean.sens[[i]])
            ci.lower.sens[i] <- quantile(draws.mean.sens[[i]], probs = .025)
            ci.upper.sens[i] <- quantile(draws.mean.sens[[i]], probs = .975)
        }
        if (missing(newdata.diff)){
            fit.matrix.sens <- as.data.frame(rbind(c(mean.sens[1], ci.lower.sens[1], ci.upper.sens[1])))
            names(fit.matrix.sens) <- c("fit", "lwr", "upr")
        }
        if (!missing(newdata.diff)) {
            fit.matrix.sens <- as.data.frame(rbind(c(mean.sens[1], ci.lower.sens[1], ci.upper.sens[1]),
                                                   c(mean.sens[2], ci.lower.sens[2], ci.upper.sens[2])))
            names(fit.matrix.sens) <- c("fit", "lwr", "upr")
            rownames(fit.matrix.sens) <- c("newdata", "newdata.diff")
        }
        return.object$fitsens <- fit.matrix.sens

        if(return.draws == TRUE & missing(newdata.diff)){
            return.object$draws.predict.sens <- draws.predict.sens[[1]]
            return.object$draws.mean.sens <- draws.mean.sens[[1]]

        }
        if (return.draws == TRUE & !missing(newdata.diff)) {
            return.object$draws.predict.sens <- list(draws.predict.sens[[1]], draws.predict.sens[[2]])
            return.object$draws.mean.sens <- list(draws.mean.sens[[1]], draws.mean.sens[[2]])
        }
    }
    
    ## If user wants predictions for z = 1 AND z = 0, or they need the average 
    ## difference in predictions between those two, or they provide newdata 
    ## and newdata.diff, need to add z and y to all those dataframes
    
    k <- 1   
    multi.return.object <- list()
    
    for (i in 1:length(data.list)){
        ## Add Y0 (unconstrained only) and z to x dataframe(s)
        if(object$constrained == FALSE){ #UNconstrained
            if(sensitive.value == "0"){
                data.list[[i]] <- cbind(data.list[[i]], object$y, 0)
            } else if (sensitive.value == "1") { # z = 1
                data.list[[i]] <- cbind(data.list[[i]], object$y - object$treat, 1) ##PROBLEM IS HERE
                                        # when given newdata with NAs:
                                        # object$y - object$treat drops some obs, xvar doesn't
                                        # because user gives full vector of observations.
            }
        } else { #constrained
            if(sensitive.value == "0") {
                data.list[[i]] <- cbind(data.list[[i]], 0)
            } else if (sensitive.value == "1") {
                data.list[[i]] <- cbind(data.list[[i]], 1)
            } else if (sensitive.value == "both") {
                data.list[[1]] <- cbind(data.list[[i]])
                data.list[[2]] <- cbind(data.list[[i]])
            }
        }
    } ## END LOOP
    
    if (sensitive.value == "both" & object$constrained == FALSE) {
        data.list[[1]] <- cbind(data.list[[1]], object$y, 0)
        data.list[[2]] <- cbind(data.list[[2]], object$y - object$treat, 1)
    } else if (sensitive.value == "both" & object$constrained == TRUE) {
        data.list[[1]] <- cbind(data.list[[1]], 0)
        data.list[[2]] <- cbind(data.list[[2]], 1)
    } ## Done with getting x matrix
    
    ## Get predictions with confidence intervals
    draws.predict <- draws.mean <- mean.obs <- lower.obs <- upper.obs <- list()
    mean <- ci.lower <- ci.upper <- sd <- c()
    if(object$outcome.reg == "logistic") {
        for(i in 1:length(data.list)){
            draws <- mvrnorm(n = 10000, coef.matrix, var.matrix)
            draws.predict[[i]] <- logistic(as.matrix(data.list[[i]]) %*% t(draws)) ##obs are rows, draws are columns
            draws.mean[[i]] <- apply(draws.predict[[i]], 2, mean) ## averaged over all observations -- 10,000 draws
            mean.obs[[i]] <- apply(draws.predict[[i]], 1, mean) ## averaged over all draws -- n observations
            lower.obs[[i]] <- apply(draws.predict[[i]], 1, quantile, probs = 0.025) ## averaged over all draws -- n observations
            upper.obs[[i]] <- apply(draws.predict[[i]], 1, quantile, probs = 0.975) ## averaged over all draws -- n observations
            mean[i] <- mean(draws.mean[[i]]) ## averaged over all draws
            sd[i] <- sd(draws.mean[[i]])
            ci.lower[i] <- quantile(draws.mean[[i]], probs = 0.025)
            ci.upper[i] <- quantile(draws.mean[[i]], probs = 0.975)
        }
    } else if (object$outcome.reg == "linear") {
        for(i in 1:length(data.list)){
            draws <- mvrnorm(n = 10000, coef.matrix, var.matrix)
            draws.predict[[i]] <- as.matrix(data.list[[i]]) %*% t(draws) ##obs are rows, draws are columns
            draws.mean[[i]] <- apply(draws.predict[[i]], 2, mean) ## averaged over all observations -- 10,000 draws
            mean.obs[[i]] <- apply(draws.predict[[i]], 1, mean) ## averaged over all draws -- n observations
            lower.obs[[i]] <- apply(draws.predict[[i]], 1, quantile, probs = 0.025) ## averaged over all draws -- n observations
            upper.obs[[i]] <- apply(draws.predict[[i]], 1, quantile, probs = 0.975) ## averaged over all draws -- n observations
            mean[i] <- mean(draws.mean[[i]]) ## averaged over all draws
            sd[i] <- sd(draws.mean[[i]])
            ci.lower[i] <- quantile(draws.mean[[i]], probs = 0.025)
            ci.upper[i] <- quantile(draws.mean[[i]], probs = 0.975)
        }
    }

    ## Up until here, all works with newdata and newdata.diff -- each return object has two elements,
    ## the first is for newdata and the second is for newdata.diff. Now it's about the reporting...

    ## When there is NO newdata:
    if (sensitive.diff == FALSE & avg == TRUE & missing(newdata.diff) & missing(newdata) & sensitive.value == "both") {#Return predictions for z = 1 & z = 0
        fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1]),
                                          c(mean[2], ci.lower[2], ci.upper[2])))
        names(fit.matrix) <- c("fit", "lwr", "upr")
        rownames(fit.matrix) <- c("Sensitive Item = 0", "Sensitive Item = 1")
    }
    if (sensitive.diff == FALSE & avg == TRUE & missing(newdata.diff) & missing(newdata) & sensitive.value == "1") {#Return predictions for z = 1
        fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1])))
        names(fit.matrix) <- c("fit", "lwr", "upr")
        rownames(fit.matrix) <- c("Sensitive Item = 1")
    }
    if (sensitive.diff == FALSE & avg == TRUE & missing(newdata.diff) & missing(newdata) & sensitive.value == "0") {#Return predictions for z = 0
        fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1])))
        names(fit.matrix) <- c("fit", "lwr", "upr")
        rownames(fit.matrix) <- c("Sensitive Item = 0")
    }
    
    ## When there IS newdata:
    if (sensitive.diff == FALSE & avg == TRUE & !missing(newdata) & missing(newdata.diff)) { ## Return predictions for newdata
        if (sensitive.value == "both") {
            fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1]),
                                              c(mean[2], ci.lower[2], ci.upper[2])))
            names(fit.matrix) <- c("fit", "lwr", "upr")
            rownames(fit.matrix) <- c("Sensitive Item = 0", "Sensitive Item = 1")
        } else if (sensitive.value == "1") {
            fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1])))
            names(fit.matrix) <- c("fit", "lwr", "upr")
            rownames(fit.matrix) <- c("Sensitive Item = 1")
        } else if (sensitive.value == "0") {
            fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1])))
            names(fit.matrix) <- c("fit", "lwr", "upr")
            rownames(fit.matrix) <- c("Sensitive Item = 0")   
        }
    }

    ## When there is newdata.diff:
    if (sensitive.diff == FALSE & avg == TRUE & !missing(newdata.diff)) {## Return predictions for newdata, newdata.diff, & diff
        ## Get difference
        diff <- draws.mean[[1]] - draws.mean[[2]] # newdata minus newdata.diff
        diff.mean <- mean(diff)
        diff.lower <- quantile(diff, probs = 0.025)
        diff.upper <- quantile(diff, probs = 0.975)
        diff.sd <- sd(diff)
        
        fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1]),
                                          c(mean[2], ci.lower[2], ci.upper[2]),
                                          c(diff.mean, diff.lower, diff.upper)))
        names(fit.matrix) <- c("fit", "lwr", "upr")
        rownames(fit.matrix) <- c("newdata", "newdata.diff", "newdata-newdata.diff")
    }
    
    if (sensitive.diff == TRUE & avg == TRUE) { ## Return predictions for z = 1, z = 0, and difference
        
        ## Get difference
        sens.diff <- draws.mean[[2]] - draws.mean[[1]] # z=1 minus z=0
        sens.diff.mean <- mean(sens.diff)
        sens.diff.lower <- quantile(sens.diff, probs = 0.025)
        sens.diff.upper <- quantile(sens.diff, probs = 0.975)
        sens.diff.sd <- sd(sens.diff)
        
        fit.matrix <- as.data.frame(rbind(c(mean[1], ci.lower[1], ci.upper[1]),
                                          c(mean[2], ci.lower[2], ci.upper[2]),
                                          c(sens.diff.mean, sens.diff.lower, sens.diff.upper)))
        names(fit.matrix) <- c("fit", "lwr", "upr")
        rownames(fit.matrix) <- c("Sensitive Item = 0", "Sensitive Item = 1", "Difference (Sensitive1 - Sensitive0)") 
    }
    
    if (sensitive.diff == FALSE & avg == FALSE) { ## Return all predictions for z = 1 and z = 0
        fit.matrix <- list(cbind(mean.obs[[1]], lower.obs[[1]], upper.obs[[1]]), ## z = 0
                           cbind(mean.obs[[2]], lower.obs[[2]], upper.obs[[2]])) ## z = 1
        names(fit.matrix) <- c("Z0", "Z1")
        colnames(fit.matrix[[1]]) <- colnames(fit.matrix[[2]]) <- c("fit", "lwr", "upr")
    }
    if (sensitive.diff == TRUE & avg == FALSE) { ## Return all predictions for z = 1 and z = 0 and the difference
        
        ## Get difference
        diff <- draws.predict[[2]] - draws.predict[[1]] ## obs are rows, draws are columns
        mean.diff <- apply(diff, 1, mean) ## averaged over all draws -- n observations
        lower.diff <- apply(diff, 1, quantile, probs = 0.025) ## averaged over all draws -- n observations
        upper.diff <- apply(diff, 1, quantile, probs = 0.975) ## averaged over all draws -- n observations
        sens.diff.sd <- apply(diff, 1, sd)
        
        fit.matrix <- list(cbind(mean.obs[[1]], lower.obs[[1]], upper.obs[[1]]), ## z = 0
                           cbind(mean.obs[[2]], lower.obs[[2]], upper.obs[[2]]), ## z = 1
                           cbind(mean.diff, lower.diff, upper.diff)) ## diff
    }
    
    return.object$fit <- fit.matrix
    if (se.fit == TRUE & sensitive.diff == FALSE) {
        return.object$se.fit <- c(sd[1], sd[2])
    } else if (se.fit == TRUE & sensitive.diff == TRUE) {
        return.object$se.fit <- c(sd[1], sd[2], sens.diff.sd)
    }
    if (return.draws == TRUE) {
        return.object$draws.predict <- draws.predict ## obs are rows, draws are columns
                                        # (can be list of 2 if sensitive.value = both) -- will be a list of 2 if
                                        # there is newdata.diff too...
        return.object$draws.mean <- draws.mean  ## averaged over all observations: vector of 10,000 draws
                                        # (can be list of 2 if sensitive.value = both)
        if (sensitive.diff == TRUE) {
            return.object$sens.diff <- sens.diff ## z=1 minus z=0 averaged over all obs: vector of 10,000 draws
        }
    }
    
    attr(return.object, "concat") <- TRUE
    
    class(return.object) <- "predict.ictreg.joint"
    return(return.object)
    
}


