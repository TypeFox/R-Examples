

####################################
#### runMCMC_: internal
#### Call runMCMCBPcpp() in C++
####################################
runMCMC_ <- function(Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, famSig_, par1_, par2_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_){
  .Call( "runMCMCBPcpp", Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, famSig_, par1_, par2_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_, PACKAGE = "geoCount" )
}
runMCMCpartialPois_ <- function(Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, famT_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_){
  .Call( "runMCMCpartialPoiscpp", Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, famT_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_, PACKAGE = "geoCount" )
}
####################################
#### set up MCMCinput
####################################
MCMCinput <- function( run = 200, run.S = 1,
    rho.family = "rhoPowerExp", 
    Y.family = "Poisson", 
    priorSigma = "Halft", parSigma = c(1, 1),
    ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1) ){
    list( run=run, run.S=run.S, rho.family=rho.family, 
          Y.family=Y.family, priorSigma=priorSigma,
          parSigma=parSigma, ifkappa=ifkappa,
          scales=scales, phi.bound=phi.bound, initials=initials )  
    }
####################################
#### runMCMC(): base 
####################################
runMCMC <- function( Y, L=0, loc, X=NULL, 
    run = 200, run.S = 1,
    rho.family = "rhoPowerExp", 
    Y.family = "Poisson", 
    priorSigma = "Halft", parSigma = c(1, 1),
    ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1), 
    MCMCinput=NULL, partial = FALSE, famT=1 ){
      
if(!is.null(MCMCinput)){
  run <- MCMCinput$run; run.S <- MCMCinput$run.S
  rho.family <- MCMCinput$rho.family
  Y.family <- MCMCinput$Y.family
  priorSigma <- MCMCinput$priorSigma
  parSigma <- MCMCinput$parSigma
  ifkappa <- MCMCinput$ifkappa
  scales <- MCMCinput$scales; 
  phi.bound <- MCMCinput$phi.bound
  initials <- MCMCinput$initials
  }
  
  if(any(Y<=0)){
      message("\nY contains non-positive element and 0.1 is added to all elements.")
      Y <- Y+0.1
    } 
  Y <- matrix(Y,,1)
  if(any(L==0)){
    L <- matrix(rep(1,nrow(Y)),,1)
    message("\nL contains zero and L is set to 1 for all locations.")
    } else { L <- matrix(L,,1)}
  U <- loc2U(loc)
  D <- as.matrix(cbind( rep(1, nrow(Y)), X ))
  
  if(rho.family=="rhoPowerExp"){
      fam = 1
    } else if(rho.family=="rhoMatern"){
        fam = 2
      } else if(rho.family=="rhoSph"){
          fam = 3
      } else {
          warning( paste("\nrho.family=", rho.family, " doesn't exist! \nrho.family=rhoPowerExp is used!", sep="") )
          fam = 1
      }
        
  if(Y.family=="Poisson"){
      famY = 1
    } else if(Y.family=="Binomial"){
        famY = 2
      } else {
          warning( paste("\nY.family=", Y.family, " doesn't exist! \nY.family=Poisson is used!", sep="") )
          famY = 1
        }

par1 <- parSigma[1]
par2 <- parSigma[2]
if(priorSigma == "Halft"){
  famSig = 1
} else if(priorSigma == "InvGamma"){
  famSig = 2
} else if(priorSigma == "Reciprocal"){
  famSig = 3
} else {
  warning( paste("\npriorSigma=", priorSigma, " doesn't exist! \npriorSigma=Reciprocal is used!", sep="") )
  famSig = 3
}

  scale <- scales[1]; mscale <- scales[2]; sscale <- scales[3]; 
  ascale <- scales[4]; kscale <- scales[5]
  alow <- phi.bound[1]; aup <- phi.bound[2]
  mini <- matrix(initials[[1]],,1); sini <- initials[[2]]; 
  aini <- initials[[3]]; kini <- initials[[4]];
  if(ncol(D) != nrow(mini)) stop("The number of covariates is not equal to the number of coefficients!")
  
  message("### MCMC Starts!\n")
  t0 <- proc.time()
  if(partial){
    res <- runMCMCpartialPois_(Y, L, U, D, run, run.S, fam, famY, famT, ifkappa, scale, mscale, sscale, ascale, kscale, alow, aup, mini, sini, aini, kini)
  } else {
    res <- runMCMC_(Y, L, U, D, run, run.S, fam, famY, famSig, par1, par2, ifkappa, scale, mscale, sscale, ascale, kscale, alow, aup, mini, sini, aini, kini)
  }
  run.time <- proc.time() - t0
  message("### MCMC Done!\n")
  
  message("### MCMC Running Time: ")
  print(run.time)
  message("### MCMC Acceptance Rate: ")
  print(res$Acc)

  res[[1]] <- matrix(res[[1]], nrow(Y), , )
  if(nrow(mini)!=1) res[[2]] <- matrix(res[[2]], nrow(mini), , )
  return(res)
  }
####################################
#### runMCMC.multiChain(): multiple
# require pre-load library{multicore}
####################################
# runMCMC.multiChain <- function(Y, L=0, loc, X=NULL, 
#     run = 200, run.S = 1,
#     rho.family = "rhoPowerExp", Y.family = "Poisson", 
#     priorSigma = "Halft", parSigma = c(1, 1),
#     ifkappa = 0,
#     scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
#     phi.bound = c(0.005, 1), 
#     initials = list(c(1), 1.5, 0.2, 1), 
#     MCMCinput=NULL, partial = FALSE, famT=1,
#     n.chn = 2, n.cores = getOption("cores")) {
#       
# if(!is.null(MCMCinput)){
#   run <- MCMCinput$run; run.S <- MCMCinput$run.S
#   rho.family <- MCMCinput$rho.family
#   Y.family <- MCMCinput$Y.family
#   priorSigma <- MCMCinput$priorSigma
#   parSigma <- MCMCinput$parSigma
#   ifkappa <- MCMCinput$ifkappa
#   scales <- MCMCinput$scales; 
#   phi.bound <- MCMCinput$phi.bound
#   initials <- MCMCinput$initials
# }
# 
# ## set different starting points    
# s.ini <- 0.1*c( mean(initials[[1]]), initials[[2]], initials[[3]], initials[[4]])
# ini <- lapply( 1:n.chn, function(tt)
#           list( initials[[1]]+s.ini[1]*rnorm(length(initials[[1]])), 
#             abs(initials[[2]]+s.ini[2]*rnorm(1)), 
#             abs(initials[[3]]+s.ini[3]*rnorm(1)),
#             ifelse(ifkappa==0, initials[[4]],
#                   abs(initials[[4]]+s.ini[4]*rnorm(1)) )
#             )
#           )
# ## MCMC        
#   message("### multiChain Starts!\n")
#   t0 <- proc.time()
#   res.prl <- mclapply(1:n.chn, function(t) 
#               runMCMC(Y, L, loc, X, run, run.S, rho.family, Y.family, 
#                       priorSigma, parSigma,
#                       ifkappa, scales, phi.bound, ini[[t]], MCMCinput=NULL,
#                       partial, famT), 
#                       mc.cores = n.cores
#           )
#   run.time <- proc.time() - t0
#   message("### multiChain Done!\n")
#   
#   message("### multiChain Running Time: ")
#   print(run.time)
#   message("### MCMC Acceptance Rate: ")
#   print(sapply(res.prl, function(tt) tt$Acc))
#   
#   res.prl
# }

####################################
#### runMCMC.sf(): multiple-chain
# require pre-load library{snowfall}
####################################
runMCMC.sf <- function(Y, L=0, loc, X=NULL, 
    run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", 
    priorSigma = "Halft", parSigma = c(1, 1),
    ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1), 
    MCMCinput=NULL, partial = FALSE, famT=1,
    n.chn = 2, n.cores = getOption("cores"), cluster.type="SOCK") {
  
  
  if(!is.null(MCMCinput)){
    run <- MCMCinput$run; run.S <- MCMCinput$run.S
    rho.family <- MCMCinput$rho.family
    Y.family <- MCMCinput$Y.family
    priorSigma <- MCMCinput$priorSigma
    parSigma <- MCMCinput$parSigma
    ifkappa <- MCMCinput$ifkappa
    scales <- MCMCinput$scales; 
    phi.bound <- MCMCinput$phi.bound
    initials <- MCMCinput$initials
  }
  ## set different starting points    
  s.ini <- 0.1*c( mean(initials[[1]]), initials[[2]], initials[[3]], initials[[4]])
  ini <- lapply( 1:n.chn, function(tt)
    list( initials[[1]]+s.ini[1]*rnorm(length(initials[[1]])), 
          abs(initials[[2]]+s.ini[2]*rnorm(1)), 
          abs(initials[[3]]+s.ini[3]*rnorm(1)),
          ifelse(ifkappa==0, initials[[4]],
                 abs(initials[[4]]+s.ini[4]*rnorm(1)) )
    )
  )
  ## MCMC
  res.prl <- NULL
  t0 <- proc.time()
  if (requireNamespace("snowfall", quietly = TRUE)) {    
    snowfall::sfInit(parallel=TRUE, cpus= n.cores, type=cluster.type)
    snowfall::sfExportAll( except=NULL, debug=FALSE )
    snowfall::sfLibrary("geoCount", character.only= TRUE)
    snowfall::sfClusterSetupRNG()
    message("### multiChain Starts!\n")    
    res.prl <- snowfall::sfLapply(1:n.chn, function(t) 
      runMCMC(Y, L, loc, X, run, run.S, rho.family, Y.family, 
              priorSigma, parSigma,
              ifkappa, scales, phi.bound, ini[[t]], MCMCinput=NULL,
              partial, famT)
    )    
    message("### multiChain Done!\n")
    snowfall::sfStop()
  } else{
    stop("Please install and load {snowfall} first before using this function!")
  }
  
  run.time <- proc.time() - t0
  message("### multiChain Running Time: ")
  print(run.time)
  message("### MCMC Acceptance Rate: ")
  print(sapply(res.prl, function(tt) tt$Acc))
  
  res.prl
}

####################################
#### Prediction at new locations
####################################
predY <- function(res.m, loc, locp, X=NULL, Xp=NULL, Lp=0, k=1, 
                  rho.family="rhoPowerExp", Y.family="Poisson", 
                  parallel=NULL, n.cores = getOption("cores"), cluster.type="SOCK"){
  t0 <- proc.time()
  
  n <- nrow(loc); np <- nrow(locp); ns <- ncol(res.m$S)
  Sp.post <- Yp <- matrix(0,np,ns)
  
  if(any(Lp==0)){
    Lp <- matrix(rep(1,np),,1)
    message("\nLp contains zero and Lp is set to 1 for all locations.")
  } else { Lp <- matrix(Lp,,1)}
  
  Uxx <- loc2U(loc)
  Uyy <- loc2U(locp)
  Uyx <- locUloc(loc, locp)
  Dx <- cbind( rep(1, n), X )
  Dy <- cbind( rep(1, np), Xp )
  
  if( !is.null(res.m$k) ){
    k.post <- res.m$k
  } else  k.post <- rep(k, ns) 
  
  
  message("### Prediction Starts!\n")
  
  if(is.null(parallel)){
    res.prl <- lapply(1:ns, function(i){  
      x <-res.m$S[,i]; 
      s <- res.m$s[i]; a <- res.m$a[i]; k <- k.post[i]
      if(is.matrix(res.m$m)){
        m <- res.m$m[,i]
      } else m <- res.m$m[i]
      
      mu.x <- Dx%*%m; mu.y <- Dy%*%m
      
      if(rho.family=="rhoPowerExp"){
        Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
        Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
        Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
      } else if(rho.family=="rhoMatern"){
        Z.xx <- s^2* rhoMatern(Uxx, a, k)
        Z.yy <- s^2* rhoMatern(Uyy, a, k)
        Z.yx <- s^2* rhoMatern(Uyx, a, k)
      } else {
        stop(paste("rho.family=", rho.family, " is not appropriate!", sep=""))
      }
      E <- mu.y + Z.yx%*%solve(Z.xx)%*%(x-mu.x)
      V <- Z.yy - Z.yx%*%solve(Z.xx)%*%t(Z.yx)
      z <- rnorm(np)
      y <- E + chol(V)%*%z
      Sp.post[,i] <- y; 
      if(Y.family=="Poisson"){
        Yp[,i] <- rpois(np, Lp*exp(y))
      } else if(Y.family=="Binomial"){
        Yp[,i] <- rbinom(np, Lp, exp(y)/(1+exp(y)))
      }
      c(Sp.post[,i], Yp[,i])
    } 
    ) # end of lapply()
    # } else if(parallel=="multicore"){
    #   res.prl <- mclapply(1:ns, function(i){  
    # x <-res.m$S[,i]; 
    # s <- res.m$s[i]; a <- res.m$a[i]; k <- k.post[i]
    # if(is.matrix(res.m$m)){
    #   m <- res.m$m[,i]
    #   } else m <- res.m$m[i]
    # 
    # mu.x <- Dx%*%m; mu.y <- Dy%*%m
    # 
    #   if(rho.family=="rhoPowerExp"){
    #       Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
    #       Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
    #       Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
    #     } else if(rho.family=="rhoMatern"){
    #         Z.xx <- s^2* rhoMatern(Uxx, a, k)
    #         Z.yy <- s^2* rhoMatern(Uyy, a, k)
    #         Z.yx <- s^2* rhoMatern(Uyx, a, k)
    #       } else {
    #           cat("Notice: rho.family=", rho.family, " doesn't exist! rho.family=rhoPowerExp will be used.\n", sep="")
    #                 Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
    #                 Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
    #                 Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
    #         }
    # E <- mu.y + Z.yx%*%solve(Z.xx)%*%(x-mu.x)
    # V <- Z.yy - Z.yx%*%solve(Z.xx)%*%t(Z.yx)
    # z <- rnorm(np)
    # y <- E + chol(V)%*%z
    # Sp.post[,i] <- y; 
    # if(Y.family=="Poisson"){
    #   Yp[,i] <- rpois(np, Lp*exp(y))
    #   } else if(Y.family=="Binomial"){
    #           Yp[,i] <- rbinom(np, Lp, exp(y)/(1+exp(y)))
    #           }
    # c(Sp.post[,i], Yp[,i])
    # }, mc.cores = n.cores
    #   ) # end of mclapply()
  } else if (requireNamespace("snowfall", quietly = TRUE)) {
    
    snowfall::sfInit(parallel=TRUE, cpus= n.cores, type=cluster.type)
    snowfall::sfExportAll( except=NULL, debug=FALSE )
    snowfall::sfLibrary("geoCount", character.only= TRUE)
    snowfall::sfClusterSetupRNG()
    res.prl <- snowfall::sfLapply(1:ns, function(i){  
      x <-res.m$S[,i]; 
      s <- res.m$s[i]; a <- res.m$a[i]; k <- k.post[i]
      if(is.matrix(res.m$m)){
        m <- res.m$m[,i]
      } else m <- res.m$m[i]
      
      mu.x <- Dx%*%m; mu.y <- Dy%*%m
      
      if(rho.family=="rhoPowerExp"){
        Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
        Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
        Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
      } else if(rho.family=="rhoMatern"){
        Z.xx <- s^2* rhoMatern(Uxx, a, k)
        Z.yy <- s^2* rhoMatern(Uyy, a, k)
        Z.yx <- s^2* rhoMatern(Uyx, a, k)
      } else {
        cat("Notice: rho.family=", rho.family, " doesn't exist! rho.family=rhoPowerExp will be used.\n", sep="")
        Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
        Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
        Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
      }
      E <- mu.y + Z.yx%*%solve(Z.xx)%*%(x-mu.x)
      V <- Z.yy - Z.yx%*%solve(Z.xx)%*%t(Z.yx)
      z <- rnorm(np)
      y <- E + chol(V)%*%z
      Sp.post[,i] <- y; 
      if(Y.family=="Poisson"){
        Yp[,i] <- rpois(np, Lp*exp(y))
      } else if(Y.family=="Binomial"){
        Yp[,i] <- rbinom(np, Lp, exp(y)/(1+exp(y)))
      }
      c(Sp.post[,i], Yp[,i])
    }
    ) # end of sfLapply()
    snowfall::sfStop()    
    
  } else{
    stop("Please install and load {snowfall} first before using this function!")
  }

  
  run.time <- proc.time() - t0
  message("### Prediction Done!\n")
  message("### Prediction Running Time: ")
  print(run.time)
  res <- matrix(unlist(res.prl),,ns)
  list(latent.predict=res[1:np,], Y.predict=res[(np+1):(2*np),])
  
}
#################
#### END
#################






          
