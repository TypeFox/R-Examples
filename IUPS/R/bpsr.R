
##################################################
# Function to implement BPSR
##################################################
bpsr <- function (Y, t, X, K = 10000, S = 1000 ) {
     t1 <- proc.time()
     if ( K!= 10000 ) { K <- K }
     if ( S!= 1000 ) { S <- S }

     est <- rep(NA, 2)
     se <-  rep(NA, 2)

     N <- length(Y)
     D <- ncol(X)

     ### Estimate propensity score
     g <- glm( t ~ X - 1 , family=binomial(link="logit"))
     ps <- g$fitted

     ### Regression using estiamted propensitys scores
     att <- lm( Y ~ t + ps - 1 ) 
     est[1] <- summary(att)$coef[1,1]
     se[1] <-  summary(att)$coef[1,2]

     ### BPSR
     ET <- c( g$coef, summary(att)$coef[,1] )
     data <- list( "N", "D", "Y", "t", "X", "ET" )
     inits <- function() { list( "beta" = runif(1), "zeta" = runif(1),
                                 "theta" = g$coef, "sigma" = runif(1)) }
     parameters <- c("beta")
     bpsr <- jags (data, inits = NULL, parameters, model.file = modelpsr, n.iter = K )
     posterior <- bpsr$BUGSoutput$summar
     est[2] <- est[1]
     se[2] <- posterior[1,2]

     estimates <- t( rbind(est, se) )
     rownames(estimates) <- c( "Phat", "BPSR" )
     colnames(estimates) <- c( "est", "se" )

     time <- proc.time() - t1
     tab <- list( estimates = estimates, time = time, sims = K, posterior = posterior)
     return ( tab )
}
