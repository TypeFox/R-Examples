######################
# EM - NOW IT'S HSMM #
######################
require(mgcv)

hsmm <- function(x, od, od.par,
                 rd        = "nonp", 
                 rd.par    = list(np = matrix(0.1, nrow = 10, ncol = 2)),
                 pi.par    = c(0.5, 0.5),
                 tpm.par   = matrix(c(0, 1, 1, 0), 2),
                 M         = NA, 
                 Q.max     = 500, 
                 epsilon   = 1e-08,  
                 censoring = 1,
                 prt       = TRUE,
                 detailed  = FALSE,
                 r.lim     = c(0.01, 100), 
                 p.log.lim = c(0.001, 0.999),
                 nu.lim    = c(0.01, 100)){
  
  # Input checking:
  od.t <- c("bern", "norm", "pois", "t", "mvnorm")
  rd.t <- c("nonp", "geom", "nbinom", "log", "pois")
  
  # formatting, renaming, additional variables
  if (length(dim(x)) != 2) {
    inputData <- as.vector(x)
  }
  if (length(dim(x)) == 2) {
    inputData <- x
  }
  tau <- get.tau(inputData)
  censoring <- as.integer(censoring)
 
  # determine number of states
  J <- length(pi.par)

  # write all initial values into one list
  Para     <- list()
  Para$pi  <- pi.par 
  Para$tpm <- t(tpm.par) 
  Para$rd  <- rd.par 
  Para$od  <- od.par 

  # lower and upper bound for uniroot
  min.r     <- as.double(r.lim[1])     # bounds$min.r
  max.r     <- as.double(r.lim[2])     # bounds$max.r
  min.p.log <- as.double(p.log.lim[1]) # bounds$min.p.log
  max.p.log <- as.double(p.log.lim[2]) # bounds$max.p.log
  min.nu    <- as.double(nu.lim[1])    # bounds$min.nu
  max.nu    <- as.double(nu.lim[2])    # bounds$max.nu

  # lower bound for some variables >= 0
  lowerBound <- 1e-300

  # selection of the maximum runlength
  M <- as.integer(M)
  if (is.na(M)){
    if (rd == "nonp"){
      M <- as.integer(dim(rd.par$np)[1])
      } else {
      M <- as.integer(min(tau, 1000))
      }
    } # endif isna(M)
  
  # iteration variable of EM alg.
  Q    <- as.integer(1)   
  
  Is.solution.reached <- FALSE
  Some.var.undef      <- FALSE

  # initialize variables for parameters estimated by EM
  EM.Para <- list()
  llh     <- c()
  for (i in 1:(Q.max + 1)){
    EM.Para[[i]] <- Para
    llh[i]       <- 0
    } # endloop i in 1:Q.max

  # variables calculated by Forward-Backward alg.
  F    <- as.double(rep(0, times = J * tau))
  L    <- as.double(rep(0, times = J * tau))
  G    <- as.double(rep(0, times = J * tau))
  L1   <- as.double(rep(0, times = J * tau))
  N    <- as.double(rep(0, times = tau))
  Norm <- as.double(rep(0, times = J * tau))
  eta  <- as.double(rep(0, times = J * M))
  xi   <- as.double(rep(0, times = J * M))
    
  # main loop of EM
  max.iterations.reached <- FALSE
  error <- as.integer(0)
  while ((!max.iterations.reached) && (!(Is.solution.reached)) && (!(Some.var.undef)) && (error == 0)){
    # Store variables for calling FB
    FB.p.tpm      <- EM.Para[[Q]]$tpm
    dim(FB.p.tpm) <- c(J * J)
    FB.pi.ini     <- EM.Para[[Q]]$pi
    FB.d          <- get.d(rd, J, M, param = EM.Para[[Q]]$rd)
    FB.pdf        <- get.pdf(inputData, od, J, M, param = EM.Para[[Q]]$od)   

    # Call Forward-Backward Alg.
    FB.result <- FB(censoring, tau, J, M, FB.d, FB.p.tpm, FB.pi.ini, FB.pdf, F, L, G, L1, N, Norm, eta, xi, error)
    error     <- FB.result[[17]]
    
    if (error == 0){   
      # Save results of FB
      F    <- FB.result[[9]]
      L    <- FB.result[[10]]
      G    <- FB.result[[11]]
      L1   <- FB.result[[12]]
      N    <- FB.result[[13]]
      Norm <- FB.result[[14]]
      eta  <- FB.result[[15]]
      xi   <- FB.result[[16]]
  
      # Change results of FB to matrices
      dim(F)    <- c(tau, J)
      dim(L)    <- c(tau, J)
      dim(G)    <- c(tau, J)
      dim(L1)   <- c(tau, J)
      dim(Norm) <- c(tau, J)
      dim(eta)  <- c(M, J)
      dim(xi)   <- c(M, J)
  
      # Calculate llh
      llh[Q] <- sum(log(N[1:tau]))

      # Output
          if (prt){
            cat("iter =", Q - 1, "\n")
            cat("LL =", formatC(llh[Q], digits = log10(1/epsilon) + 2, format = "f"), "\n") 
            }
  
      # Check whether solution is reached
      k <- Q - 1
      if (k > 1){
        while ((k >= 1) && is.na(llh[k])){
          k <- k - 1
          } #endwhile k >= 1
        } # endif k > 1
      #Q >=3 for numerical reasons
      if ((Q >= 3) && (k >= 1) && (is.na(llh[Q]) == FALSE) && (is.na(llh[k]) == FALSE)){
        if ((llh[Q] - llh[k]) / abs(llh[k]) < epsilon){
          Is.solution.reached <- TRUE
          } # endif llh < epsilon
        } # endif Q >= 1 & k >= 1 ...
  
      # Maximum iterations reached?
      if (Q + 1 > Q.max + 1) {
        max.iterations.reached <- TRUE
      }
      
      #  Is solution reached?
      if ((!(Is.solution.reached) && (!max.iterations.reached))){
        # Reestimation pi
        for (j in 1:J){
          EM.Para[[Q + 1]]$pi[j] <- L[1, j]
          if (is.na(EM.Para[[Q + 1]]$pi[j])){
            EM.Para[[Q + 1]]$pi[j] <- lowerBound
            } else{
            if (EM.Para[[Q + 1]]$pi[j] < lowerBound){
              EM.Para[[Q + 1]]$pi[j] <- lowerBound
              } # endif EM.Para$pi
            } # endif is.na(EM.Para$pi)
        } # endloop j in 1:J
           
        # Reestimation tpm
        for (i in 1:J){
          r <- sum(L1[1:(tau - 1), i])
          for (j in 1:J){
            z <- sum(G[2:tau, j] * EM.Para[[Q]]$tpm[j, i] * F[1:(tau - 1), i])
            EM.Para[[Q + 1]]$tpm[j, i] = z / r;
            } # endloop j in 1:J
          } # endloop i in 1:j
  
        # Reestimation of parameters of rd #######################################################
        if (rd == "nonp"){
          # Reestimation d
          for (j in 1:J){
            if (censoring == 1){
              EM.Para[[Q + 1]]$rd$np[1:M, j] <- eta[1:M, j] / (sum(L1[1:(tau - 1), j]) + L[tau, j])
            } else{
              EM.Para[[Q + 1]]$rd$np[1:M, j] <- eta[1:M, j] / sum(L1[1:(tau - 1), j])
            }
            } # endloop j in 1:J
          } # endif rd == "non.parametric"
  
        if (rd == "geom"){
          # Reestimation p
          for (j in 1:J){
            w <- sum(eta[1:M, j])
            z <- sum(1:M * eta[1:M, j])
            EM.Para[[Q + 1]]$rd$p[j] <- w / z
            } # endloop j in 1:J
        } # endif rd == "geometric"
  
        if (rd == "nbinom"){
          # Reestimation r
          for (j in 1:J){
            z <- try(uniroot(f = r.nb.update, j = j, M = M, eta = eta, lower = min.r, upper = max.r)$root)
            if (inherits(z, "try-error")){
              Some.var.undef <- TRUE
              } else {
              EM.Para[[Q + 1]]$rd$r[j] <- z
            } # endif inherits.error
          } # endloop j in 1:J
  
          # Reestimation pi
          if (!(Some.var.undef))
            for (j in 1:J){
              w <- sum(eta[1:M, j])
              z <- sum(eta[1:M, j] * (1:M - 1 + EM.Para[[Q + 1]]$rd$r[j]))
              EM.Para[[Q + 1]]$rd$pi[j] <- EM.Para[[Q + 1]]$rd$r[j] * w / z
              } # endloop j in 1:J
          } # endif rd == "negative.binomial"
  
        if (rd == "log"){
          # Reestimation p
          for (j in 1:J){
            z <- try(uniroot(f = p.log.update, j = j, M = M, eta = eta, lower = min.p.log, upper = max.p.log)$root)
            if (inherits(z, "try-error")){
              Some.var.undef <- TRUE
              } else {
              EM.Para[[Q + 1]]$rd$p[j] <- z
              } # endif inherits.error
            } # endloop j in 1:J
         } # endif rd == "logarithmic"
  
  
  #      if (rd == "logarithmic.geometric"){
  #        # Reestimation p.log
  #        for (j in 1:J){
  #          result <- optim(fn = loggeom.update, par = c(qlogis(EM.Para[[Q]]$p.loggeom[j]), qlogis(EM.Para[[Q]]$theta[j])),
  #                          j = j, M = M, eta = eta)
  #          par <- result$par
  #          EM.Para[[Q + 1]]$rd$p.loggeom[j] <- plogis(par[1])
  #          EM.Para[[Q + 1]]$rd$theta[j]     <- plogis(par[2])
  #          }
  #        } # endif rd == "logarithmic.geometric"
  
        if (rd == "pois"){
          # Reestimation lambda
          for (j in 1:J){
            EM.Para[[Q + 1]]$rd$lambda[j] <- sum(eta[1:M, j] * (1:M - 1)) / sum(eta[1:M, j])
            } # endloop j in 1:J
          } # endif rd == "Poisson"
  
        # Reestimation of parameters of od #######################################################
        if (od == "bern"){
          # Number of possible outcomes of the Bernoulli distribution
          Y <- dim(od.par$b)[1]
          # Reestimation b
          for (j in 1:J){
            s <- sum(L[1:tau, j])
            z <- sum(L[inputData[1:tau] == 1, j])
            EM.Para[[Q + 1]]$od$b[j] <- z / s
            } # endloop j in 1:J
          } # endif od == "Bernoulli"
  
        if (od == "norm"){
          # Reestimation mu
          for (j in 1:J){
            r <- sum(L[1:tau, j] * inputData[1:tau])
            s <- sum(L[1:tau, j])
            EM.Para[[Q + 1]]$od$mean[j] = r / s;
            } # endloop j in 1:J
  
          # Reestimation var
          for (j in 1:J){
            r <- sum(L[1:tau, j] * (inputData[1:tau] - EM.Para[[Q + 1]]$od$mean[j]) * 
                 (inputData[1:tau] - EM.Para[[Q + 1]]$od$mean[j]))
            s <- sum(L[1:tau, j])
            EM.Para[[Q + 1]]$od$var[j] <- r / s;
            } #endloop j in 1:J
          } # endif od == "norm"
  
        if (od == "pois"){
          # Reestimation lambda.obs
          for (j in 1:J){
            r <- sum(L[1:tau, j] * inputData[1:tau])
            s <- sum(L[1:tau, j])
            EM.Para[[Q + 1]]$od$lambda[j] <- r / s;
            } # endloop j in 1:J
          } # endif od == "Poisson"
  
        if (od == "t"){                 
          for (j in 1:J){		
            U <- rep(0, times = tau)
            for (t in 1:tau){
              U[t] <- (EM.Para[[Q]]$od$df[j] + 1) / 
                      (EM.Para[[Q]]$od$df[j] + (inputData[t] - EM.Para[[Q]]$od$mean[j])^2 / EM.Para[[Q]]$od$var[j]);
              } # endloop t in 1:tau
            
            # Reestimation mu
            w <- sum(L[1:tau, j] * U[1:tau] * inputData[1:tau])
            z <- sum(L[1:tau, j] * U[1:tau])
            EM.Para[[Q + 1]]$od$mean[j] <- w / z;
  
            # Reestimation var
            w <- sum(L[1:tau, j] * U[1:tau] * (inputData[1:tau] - EM.Para[[Q + 1]]$od$mean[j])^2)
            z <- sum(L[1:tau, j])
            EM.Para[[Q + 1]]$od$var[j] <- w / z;
  
            # Reestimation nu
            w <- sum(L[1:tau, j] * (log(U[1:tau]) - U[1:tau])) 	
            z <- sum(L[1:tau, j])
            r <- try(exp(uniroot(f = nu.update, dblPara = w / z, lower = log(min.nu), upper = log(max.nu))$root))
            if (inherits(r, "try-error")){
              Some.var.undef <- TRUE
              } else {
              EM.Para[[Q + 1]]$od$df[j] <- r
              } # endif inherits.error
            } # endloop j in 1:J
          } # endif od == "Student.t"
    
          if (od == "mvnorm"){
            obs.dim <- dim(inputData)[1]
            # Reestimation mu
            for (j in 1:J){
              s <- sum(L[1:tau, j])
              for (k in 1:obs.dim) {
                r <- sum(L[1:tau, j] * inputData[k, 1:tau])
                EM.Para[[Q + 1]]$od$mean[k, j] = r / s;
                } # endloop k in 1:obs.dim
              } # endloop j in 1:J
    
            # Reestimation sigma
            for (j in 1:J){
              s <- sum(L[1:tau, j])   
              a <- array(0, dim = c(obs.dim, obs.dim, tau))
              for (u in 1:tau) {
                m <- (inputData[,u] - EM.Para[[Q + 1]]$od$mean[,j]) %*% 
                     t((inputData[,u] - EM.Para[[Q + 1]]$od$mean[,j]))
                   a[,,u] <- L[u, j] * array(m)
                }
              r <- rowSums(a, dims = 2)   
              y <- r / s
              for (k in 1:obs.dim) {
                for (l in 1:obs.dim) {
                  EM.Para[[Q + 1]]$od$sigma[k, l, j] <- y[k, l];
                  }
                }
              } #endloop j in 1:J
            } # endif od == "mvnorm"

          Q <- Q + 1
          
      } # Is solution reached or maximum iterations reached?
    } # Has an error occured?  

  } # while Q < Q.max and solution not reached
  
  # save results
  Para    <- EM.Para[[Q]]
  llh.out <- llh[Q]
  Q.used  <- Q - 1

  # reformat tpm
  Para$tpm <- t(Para$tpm)

  # control variable that contains additional information on the results
  ctrl                  <- list()
  ctrl$solution.reached <- Is.solution.reached
  ctrl$error            <- error
  
  # details on the iterations
  if (detailed){ 
    details <- list()
    for (i in 1:Q){
      details$logl[[i]] <- llh[i]
      details$pi[[i]]   <- EM.Para[[i]]$pi
      details$tpm[[i]]  <- t(EM.Para[[i]]$tpm)
      details$rd[[i]]   <- EM.Para[[i]]$rd
      details$od[[i]]   <- EM.Para[[i]]$od
      }
    ctrl$details <- details
    }
    
  # return results
  out <- list(call                = match.call(),
              iter                = Q.used, 
              logl                = llh.out,
              para                = Para,
              ctrl                = ctrl)
  class(out) <- "hsmm"
  return(out)
  }
