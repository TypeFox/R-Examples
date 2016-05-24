# This code is useful for estimation of the COGARCH(p,q) model according the PseudoLogLikelihood
# maximization procedure developed in Iacus et al. 2015
#

PseudoLogLik.COGARCH <- function(yuima, start, method="BFGS", fixed = list(),
                     lower, upper, Est.Incr, call, grideq, ...){

  if(is(yuima,"yuima")){
    model <- yuima@model
    info <- model@info
    Data <- onezoo(yuima)
  }

  time <- index(Data)
  Obs <- as.numeric(as.matrix(Data)[,1])
  my.env <- new.env()
  param <- unlist(start)

  meas.par <- model@parameter@measure

  if(length(meas.par)==0 && Est.Incr=="IncrPar"){
    yuima.warn("The dimension of measure parameters is zero, yuima changes 'Est.Incr = IncrPar' into 'Est.Incr = Incr'")
    Est.Incr <- "Incr"
  }

  ar.names <- paste0(info@ar.par,c(1:info@q))
  ma.names <- paste0(info@ma.par,c(1:info@p))

  start.state <- paste0(paste0(info@Latent.var,0),c(1:info@q))

  e <- matrix(0, info@q,1)
  e[info@q,1] <- 1
  assign("e", e, envir = my.env)

  assign("start.state", start.state, envir = my.env)
  assign("q", info@q, envir = my.env)
  assign("p", info@p, envir = my.env)
  assign("B", matrix(0,info@q,info@q),  envir = my.env)
  # consider two cases:
      # 1) equally spaced grid
  if(grideq){
    assign("Deltat", tail(index(Data),1)/length(index(Data)), envir = my.env)
  #  assign("Deltat", diff(time)[1], envir = my.env)
      # 2) no-equally spaced grid
    assign("grideq", TRUE, envir = my.env)
  }else{
    assign("Deltat", diff(time), envir = my.env)
    assign("grideq", FALSE, envir = my.env)
  }
  assign("Obs", (diff(Obs))^2, envir = my.env)
  assign("nObs",length(Obs),envir = my.env)
  assign("ar.names", ar.names, envir = my.env)
  assign("ma.names", ma.names, envir = my.env)
  assign("loc.par",info@loc.par, envir = my.env)

  I <- diag(info@q)
  assign("I",I, envir = my.env)

  out<-NULL
  if(length(lower)==0 && length(upper)>0 && length(fixed)==0){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                   method = method, upper=upper, env = my.env,...)


  }

  if(length(lower)==0 && length(upper)==0 && length(fixed)>0){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                   method = method, fixed=fixed, env = my.env,...)

  }


  if(length(lower)>0 && length(upper)==0 && length(fixed)==0){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                   method = method, lower=lower, env = my.env,...)
    }

  if(length(lower)>0 && length(upper)>0 && length(fixed)==0){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                   method = method, upper = upper,
                   lower=lower, env = my.env,...)
  }


  if(length(lower)==0 && length(upper)>0 && length(fixed)>0){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                   method = method, upper = upper,
                   fixed = fixed, env = my.env,...)
  }

  if(length(lower)>0 && length(upper)==0 && length(fixed)>0){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                   method = method, lower = lower,
                   fixed = fixed, env = my.env,...)
  }


  if(length(lower)>0 && length(upper)>0 && length(fixed)>0){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                   method = method, lower = lower,
                   fixed = fixed, upper = upper,
                   env = my.env,...)
  }


  if(is.null(out)){
    out <- optim(par=param, fn=minusloglik.COGARCH1,
                 method = method, env = my.env,...)
  }

#                 control= list(maxit=100),
   # Write the object mle with result


  bvect<-out$par[ar.names]
  bq<-bvect[1]
  avect<-out$par[ma.names]
  a1<-avect[1]

  a0<-out$par[info@loc.par]

  if(length(meas.par)!=0){
    idx.dumm<-match(meas.par,names(out$par))
    out$par<-out$par[- idx.dumm]
  }

  idx.dumm1<-match(start.state,names(out$par))
  coef <- out$par[-idx.dumm1]

  vcov<-matrix(NA, length(coef), length(coef))
  names_coef<-names(coef)
  colnames(vcov) <- names_coef
  rownames(vcov) <- names_coef
  mycoef <- start
  # min <- out$value
  objFun <- "PseudoLogLik"
  min <- numeric()


  res<-new("cogarch.gmm", call = call, coef = coef, fullcoef = unlist(coef),
           vcov = vcov, min = min, details = list(),
           method = character(),
           model = model,
           objFun = objFun
  )


    return(res)
}

minusloglik.COGARCH1<-function(param,env){

#   assign("start.state", start.state, envir = my.env)
#   assign("q", info@q, envir = my.env)
#   assign("p", info@p, envir = my.env)
#   assign("time", time, envir = my.env)
#   assign("Obs", Obs, envir = my.env)
#   assign("nObs",length(Obs),envir = my.env)
#   assign("ar.names", ar.names, envir = my.env)
#   assign("ma.names", ma.names, envir = my.env)
#   assign("loc.par",info@loc.par, envir = my.env)
  a0<-param[env$loc.par]
  bq<-param[env$ar.names[env$q]]
  a1<-param[env$ma.names[1]]

  stateMean <- a0/(bq-a1)*as.matrix(c(1,numeric(length=(env$q-1))))

 param[env$start.state]<-stateMean
 state <- stateMean
#  state <- param[env$start.state]
  DeltaG2 <- env$Obs
  B <- env$B
  if(env$q>1){
    B[1:(env$q-1),] <- c(matrix(0,(env$p-1),1), diag(env$q-1))
  }
  B[env$q,] <- -param[env$ar.names[env$q:1]]
  a<-matrix(0,env$q,1)
  a[1:env$p,]<-param[env$ma.names]

  ta<-t(a)
  e <- env$e
  Btilde <-B+e%*%ta
  InvBtilde <- solve(Btilde)
  V1 <- a0+ta%*% state
  V <- V1[1]
  Deltat<- env$Deltat
  I <- env$I
#   VarDeltaG <- 0
   PseudologLik <- 0

#   DeltatB1 <- lapply(as.list(Deltat), function(dt,B){expm(B*dt)} , B)
# #   DeltatB <- lapply(as.list(Deltat), "*" , B)
# #   assign("DeltatB",as.list(DeltatB),.GlobalEnv)
#   outputB <- matrix(unlist(DeltatB1), ncol = env$q, byrow = TRUE)
  if(env$grideq){
  DeltatB1 <- expm(B*Deltat)

#   DeltatB2 <- lapply(as.list(Deltat), function(dt,B){expm(B*dt)} , Btilde)
#   #   DeltatB <- lapply(as.list(Deltat), "*" , B)
#   #   assign("DeltatB",as.list(DeltatB),.GlobalEnv)
#   outputB2 <- matrix(unlist(DeltatB2), ncol = env$q, byrow = TRUE)

  DeltatB2 <- expm(Btilde*Deltat)

#   DeltatB3 <- lapply(as.list(-Deltat), function(dt,B){expm(B*dt)} , Btilde)
#   outputB3 <- matrix(unlist(DeltatB3), ncol = env$q, byrow = TRUE)

  DeltatB3 <- expm(-Btilde*Deltat)

  dummyMatr <- ta%*%DeltatB2%*%InvBtilde%*%(I-DeltatB3)
  dummyeB1 <- e%*%ta%*%DeltatB1

  #aa <- .Call("myfun1", DeltatB1, state)
  PseudologLik <-.Call("pseudoLoglik_COGARCH1", a0, bq, a1, stateMean, Q=as.integer(env$q),
                         DeltaG2, Deltat, DeltatB1,  a, e,
                       V, nObs=as.integer(env$nObs-1),
                      dummyMatr, dummyeB1)
  #cat(sprintf("\n%.5f ", PseudologLik))
#
#
#   PseudologLikR<-0
#   V1 <- a0+ta%*% state
#   V <- V1[1]
#
#
#     for(i in c(1:(env$nObs-1))){
#
# #       cat(sprintf("\n dummy1R %.10f ",dummyMatr%*%(state-stateMean) ))
# #       d <- 0
# #       for(j in 1:2){
# #         d<- d+dummyMatr[1,j]*(state[j]-stateMean[j])
# #         cat(sprintf("\n %d dummy1R %.10f ",j,d ))
# #       }
#       VarDeltaG <- a0*Deltat*bq/(bq-a1)+ dummyMatr%*%(state-stateMean)
#
#         VarDeltaG <- VarDeltaG[1]
#         #state <- (I+DeltaG2[i]/V*e%*%ta)%*%DeltatB1%*%state+a0*DeltaG2[i]/V*e
#         state <- DeltatB1%*%state+DeltaG2[i]/V*dummyeB1%*%state+a0*DeltaG2[i]/V*e
#         #state <- DeltatB1%*%state+dummyeB1%*%state
# #         cat(sprintf("\n d1 %.10f d2 %.10f", DeltatB1%*%state, dummyeB1%*%state));
#         V <- a0+ta%*% state
#         V <- V[1]
#         if(is.nan(VarDeltaG))
#           VarDeltaG<- 10^(-6)
#         PseudologLikR<--0.5*(DeltaG2[i]/VarDeltaG+log(VarDeltaG)+log(2.*3.14159265))+PseudologLikR
# #         cat(sprintf("\n%.5f -  %.5f %.5f  -  %.5f",VarDeltaG, state[1], state[2],V ))
#         cat(sprintf("\n Part %.10f partial %.10f ", PseudologLikR, VarDeltaG))
#         if(is.nan(V))
#           V <- 10^(-6)
#
#      }
# #
#    cat(sprintf("\n%.5f -  %.5f",PseudologLikR, PseudologLik))
# # #
# #

  }else{
    DeltatB1 <- lapply(as.list(Deltat), function(dt,B){expm(B*dt)} , B)
    DeltatB2 <- lapply(as.list(Deltat), function(dt,B){expm(B*dt)} , Btilde)
    DeltatB3 <- lapply(as.list(-Deltat), function(dt,B){expm(B*dt)} , Btilde)

    for(i in c(1:(env$nObs-1))){
        VarDeltaG <- as.numeric(a0*Deltat[i]*bq/(bq-a1)+ta%*%DeltatB2[[i]]%*%InvBtilde%*%(I-DeltatB3[[i]])%*%(state-stateMean))
        state <- (I+DeltaG2[i]/V*e%*%ta)%*%DeltatB1[[i]]%*%state+a0*DeltaG2[i]/V*e
        V <- as.numeric(a0+ta%*% state)
        PseudologLik<--1/2*(DeltaG2[i]/VarDeltaG+log(VarDeltaG)+log(2*pi))+PseudologLik
      }
  }
#
#   PseudologLik <- 0
#
#   for(i in c(1:(env$nObs-1))){
#       VarDeltaG <- a0*Deltat*bq/(bq-a1)+ dummyMatr%*%(state-stateMean)
#       #state <- (I+DeltaG2[i]/V*e%*%ta)%*%DeltatB1%*%state+a0*DeltaG2[i]/V*e
#       state <- DeltatB1%*%state+DeltaG2[i]/V*dummyeB1%*%state+a0*DeltaG2[i]/V*e
#       V <- as.numeric(a0+ta%*% state)
#       PseudologLik<--1/2*(DeltaG2[i]/VarDeltaG+log(VarDeltaG)+log(2*pi))+PseudologLik
#     }

#   for(i in c(1:(env$nObs-1))){
#     VarDeltaG <- as.numeric(a0*Deltat[i]*bq/(bq-a1)+t(a)%*%DeltatB2[[i]]%*%InvBtilde%*%(I-DeltatB3[[i]])%*%(state-stateMean))



  #     state <- (I+DeltaG2[i]/V*e%*%t(a))%*%DeltatB1[[i]]%*%state+a0*DeltaG2[i]/V*e
#     V <- as.numeric(a0+t(a)%*% state)
#     PseudologLik<--1/2*(DeltaG2[i]/VarDeltaG+log(VarDeltaG)+log(2*pi))
#   }

#     dummyMatr <- ta%*%DeltatB2%*%InvBtilde%*%(I-DeltatB3)
#     dummyeB1 <- e%*%ta%*%DeltatB1
#     PseudologLik1 <- 0
#     for(i in c(1:(env$nObs-1))){
#       VarDeltaG <- as.numeric(a0*Deltat*bq/(bq-a1)+dummyMatr%*%(state-stateMean))
#       state <- DeltatB1%*%state+DeltaG2[i]/V*dummyeB1%*%state+a0*DeltaG2[i]/V*e
#       V <- as.numeric(a0+ta%*% state)
#       PseudologLik1 <- -1/2*(DeltaG2[i]/VarDeltaG+log(VarDeltaG)+log(2*pi))
#       if(is.finite(PseudologLik1)){
#         PseudologLik <- PseudologLik1 + PseudologLik
#       }
#     }

    minusPseudoLogLik <- -PseudologLik
   return(minusPseudoLogLik)
}

#   res<-.Call("pseudoLoglik_COGARCH", a0, bq, a1, stateMean, Q=as.integer(env$q),
#                        state, DeltaG2, Deltat, DeltatB, B, a, e,
#                        Btilde, InvBtilde, V, I, VarDeltaG,
#                        PseudologLik, nObs = as.integer(env$nObs-1), fn = quote(expm(x)) , rho= .GlobalEnv,
#                        PACKAGE = "yuima")
#
#   output <- matrix(unlist(res), ncol = env$q, byrow = TRUE)
#   res <- res

