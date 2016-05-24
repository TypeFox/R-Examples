"portest" <-
function (obj,lags=seq(5,30,5),order=0,test=c("PenaRodriguez","BoxPierce","LjungBox",
        "Hosking","LiMcLeod"),MonteCarlo=TRUE,Kernel=FALSE,nworkers=1,NREP=1000,
        InfiniteVarianceQ=FALSE,SquaredQ=FALSE,
        func=list(SimModel=NULL, FitModel=NULL),pkg=NULL,SetSeed=TRUE)  
{
    test <- match.arg(test)
     TestType <- "0"
    if (class(obj) == "ts" || class(obj) == "numeric" || class(obj) == 
        "matrix" || (class(obj)[1] == "mts" && class(obj)[2] == 
        "ts")) 
        TestType <- "1"
    if (class(obj) == "ar" || class(obj) == "arima0" || class(obj) == 
        "Arima" || class(obj) == "varest" || class(obj) == "FitAR" || 
        class(obj) == "FitFGN" || class(obj) == "garch" || 
        class(obj) == "fGARCH" || class(obj) == "list") 
        TestType <- "2"
    if (TestType == "0") 
        stop("obj must be class ar, arima0, Arima, varest, FitAR, 
             FitFGN, garch, fGARCH, ts, numeric, matrix, (mts ts), or list")
    if (TestType == "1") {
        res <- as.ts(obj)
        Order <- order
    }
    else {
        GetResid <- GetResiduals(obj)
        res <- GetResid$res
        Order <- GetResid$order
    }
    k <- NCOL(res)
    n <- NROW(res)
    res <- matrix(res,ncol=k,nrow=n)
    if (MonteCarlo == FALSE){ 
      if (test =="PenaRodriguez")
        return(gvtest(res, lags, Order, SquaredQ))
      else if (test =="BoxPierce")
        return(BoxPierce(res, lags, Order, SquaredQ))  
      else if (test =="LjungBox")
        return(LjungBox(res, lags, Order, SquaredQ)) 
      else if (test =="Hosking")
        return(Hosking(res, lags, Order, SquaredQ)) 
      else if (test =="LiMcLeod")
        return(LiMcLeod(res, lags, Order, SquaredQ)) 
   }           
    else {
          if (test =="PenaRodriguez"){
              if (Kernel==FALSE){
              ans <- gvtest(res, lags, Order, SquaredQ)
              obs.stat <- ans[, 2]
              }
              else
              {
                ans <- gvtest(res, lags, Order, SquaredQ,Kernel=TRUE)
                obs.stat <- ans[, 2]
              }
         }
         else if (test =="BoxPierce"){
             ans <- BoxPierce(res, lags, Order, SquaredQ)
             obs.stat <- ans[, 2]
        }
        else if (test =="LjungBox"){
            ans <- LjungBox(res, lags, Order, SquaredQ)
            obs.stat <- ans[, 2]
       }
       else if (test =="Hosking"){
           ans <- Hosking(res, lags, Order, SquaredQ)
           obs.stat <- ans[, 2]
      }
      else if (test =="LiMcLeod"){
           ans <- LiMcLeod(res, lags, Order, SquaredQ)
           obs.stat <- ans[, 2]
      }
      if (InfiniteVarianceQ){ 
          StableParameters <- matrix(fitstable(res), ncol = 4)
          ALPHA<-StableParameters[,1]
          BETA<-StableParameters[,2]
          GAMMA<-StableParameters[,3]
          DELTA<-StableParameters[,4]
      }
      else 
          StableParameters <- NA
      if (TestType == 1)  
           sigma <- matrix(stats::acf(res, lag.max = 1, plot = FALSE,type = "covariance")$acf[1, , ], k, k)
      else {        
                        if (all(class(obj) == "ar")) {
                            p <- obj$order
                            q <- 0
                            if (k == 1) 
                                Model <- 1
                            else 
                                 Model <- 2
                            sigma <- obj$var.pred
                            if (is.array(obj$ar)) {
                                arrayphi <- array(numeric(k * k * p), dim = c(k^2, p))
                                for (i in 1:p) arrayphi[, i] <- c(obj$ar[i, ,])
                                phi <- array(c(arrayphi), dim = c(k, k, p))
                            }
                            else phi <- obj$ar
                            theta <- NULL
                            if (!is.null(obj$x.intercept)){ 
                                constant <- obj$x.intercept
                                intercept <- TRUE
                            }
                            else {
                                 constant <- rep(0, k)
                                 intercept <- FALSE
                            }
                            trend <- rep(0, k)
                            demean <- obj$x.mean
                            if(all(demean==0)) 
                               Demean <- FALSE
                            else
                               Demean <- TRUE
                            }
                         else if (all(class(obj) == "varest")) {
                              sigma <- summary(obj)[[3]]
                              p <- obj$p
                              q <- 0
                              theta <- NULL
                              if (obj$type == "none") {
                                  Model <- "3A"
                                  Phi <- matrix(numeric(p * k^2), nrow = k, ncol = p * k)
                                  constant <- NA
                                  trend <- NA
                                  demean <- NA
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  phi <- array(Phi, dim = c(k, k, p))
                              }
                              else if (obj$type == "const") {
                                  Model <- "3B"
                                  Phi <- matrix(numeric(k + p * k^2), nrow = k,ncol = p * k + 1)
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  constant <- Phi[, p * k + 1]
                                  trend <- NA
                                  demean <- NA
                                  phi <- array(Phi[, -(p * k + 1)], dim = c(k, k, p))
                              }
                              else if (obj$type == "trend") {
                                  Model <- "3C"
                                  Phi <- matrix(numeric(k + p * k^2), nrow = k,ncol = p * k + 1)
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  constant <- NA
                                  trend <- Phi[, p * k + 1]
                                  demean <- NA
                                  phi <- array(Phi[, -(p * k + 1)], dim = c(k, k, p))
                               }
                               else if (obj$type == "both") {
                                  Model <- "3D"
                                  Phi <- matrix(numeric(2 * k + p * k^2), nrow = k, ncol = p * k + 2)
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  constant <- Phi[, p * k + 1]
                                  trend <- Phi[, p * k + 2]
                                  demean <- NA
                                  phi <- array(Phi[, -((p * k + 1):(p * k + 2))], dim = c(k, k, p))
                                }
                         }
                         else if (all(class(obj) == "arima0") || all(class(obj) == "Arima")) {
                                Model <- 4
                                pdq <- obj$arma
                                if(as.integer(pdq[3])!= 0 || as.integer(pdq[4]) != 0 || as.integer(pdq[7]) != 0) 
                                     stop("portes is applied only for nonseasonal arima models")
	                            p <- pdq[1]
	                            q <- pdq[2]
	                            d <- pdq[6]
                                phi <- theta <- NULL
                                if (p > 0 && q == 0)               
                                    phi <- as.vector(obj$coef[1:p])
                                else if (p > 0 && q >0){
                                     phi <- as.vector(obj$coef[1:p])
                                     theta <- as.vector(obj$coef[(p+1):(p+q)])
                                }
                                else if (p==0 && q>0)
                                     theta <- as.vector(obj$coef[1:q])
                                sigma <- obj$sigma2
                                if (length(obj$coef)==p+q+2) 
                                     constant <- as.vector(obj$coef[p+q+2])
                                else
                                     constant <- 0  
                                trend <- 0
                                if (d==0) 
                                     demean <- as.vector(obj$coef[p+q+1])
                                else
                                     demean <- 0   
                         }
                         else if (all(class(obj) == "FitAR")) {
                               Model <- 5
                               p <- length(obj$phiHat)
                               q <- 0
                               phi <- obj$phiHat
                               theta <- NULL
                               sigma <- obj$sigsqHat
                               constant <- 0
                               trend <- 0
                               demean <- obj$muHat
                         }
                         else if (all(class(obj) == "FitFGN")) {
                               Model <- 6
                               phi <- theta <- NULL
                               p <- 1
                               q <- 0
                               H <- obj$H
                               sigma <- sqrt(obj$sigsq)
                               constant <- 0
                               trend <- 0
                               demean <- obj$muHat
                         }
                         else if (all(class(obj) == "garch")) {
                               Model <- 7
	                           GARCHOrder <- as.vector(obj$order)
                               p2 <- GARCHOrder[2]
                               q2 <- GARCHOrder[1]
                               Beta <- 0
                               demean <- as.vector(obj$coef[1])
                               Alpha <- as.vector(obj$coef[2:(p2+1)])
                               if (q2 > 0)
                                  Beta <- as.vector(obj$coef[(p2+2):(p2+q2+1)])
                               p <- 0
                               q <- 0
                               phi <- NULL
                               theta <- NULL
                               sigma <- 0
                               constant <- 0
                               trend <- 0
                               cond.dist <-"norm"
                         }
                         else if (all(class(obj) == "fGARCH")) {
                               Model <- 8
		                       GARCHOrder <- as.vector(obj@fit$series$order)
	                           p1 <- GARCHOrder[1]
	                           q1 <- GARCHOrder[2]
	                           p2<- GARCHOrder[3]
	                           q2<- GARCHOrder[4]
		                       if (obj@fit$params$mu==0){
    		                      demean <- 0
    			                  if (p1>0 && q1 >0){
      			                     phi <- as.vector(obj@fit$par)[1:p1]
  				                     theta <- as.vector(obj@fit$par)[(p1+1):(p1+q1)]
   			                      }
    			                  else if (p1>0 && q1 == 0){
     				                  phi <- as.vector(obj@fit$par)[1:p1]
    				                  theta <- 0 
   			                       }
    			                   else if (p1 ==0 && q1 >0){
   				                      phi <- 0
    			                      theta <- as.vector(obj@fit$par)[1:q1]
  			                       }
    			                   else if (p1 ==0 && q1 == 0){
     				                  phi <- 0
     				                  theta <- 0 
 			                       }
		                         OMEGA <- as.vector(obj@fit$par)[p1+q1+1]
	                             Alpha <- as.vector(obj@fit$par)[(p1+q1+2):(p1+q1+p2+1)]
	                             Beta <- as.vector(obj@fit$par)[(p1+q1+p2+2):(p1+q1+p2+q2+1)]
                               }
                               else{
			                     demean <- as.vector(obj@fit$par)[1]
    			                 if (p1>0 && q1 >0){
    				                 phi <- as.vector(obj@fit$par)[2:(p1+1)]
    				                 theta <- as.vector(obj@fit$par)[(p1+2):(p1+q1+1)]
  		                         }
   			                     else if (p1>0 && q1 == 0){
     				                 phi <- as.vector(obj@fit$par)[2:(p1+1)]
     				                 theta <- 0 
   			                     }
    			                 else if (p1 ==0 && q1 >0){
    				                 phi <- 0
     				                 theta <- as.vector(obj@fit$par)[2:(q1+1)]
   		                         }
   			                     else if (p1 ==0 && q1 == 0){
    				                 phi <- 0
                                     theta <- 0 
                                 }
		                        OMEGA <- as.vector(obj@fit$par)[p1+q1+2]
		                        Alpha <- as.vector(obj@fit$par)[(p1+q1+3):(p1+q1+p2+2)]
		                        Beta <- as.vector(obj@fit$par)[(p1+q1+p2+3):(p1+q1+p2+q2+2)]
                              }
                             Delta <- obj@fit$params$delta
		                     SKEW <- obj@fit$params$skew
		                     SHAPE <- obj@fit$params$shape
	                         cond.dist <- obj@fit$params$cond.dist
		                     if (obj@fit$params$leverage==FALSE) 
  			                     Leverage <- NULL
		                     else
 			                     Leverage <- obj@fit$params$leverage
                             p <- p1
                             q <- q1
                             sigma <- 0
                             constant <- 0
                             trend <- 0  
                         }
                         else if (all(class(obj) == "list")) {
                               Model <- 9
                           stopifnot(!is.null(pkg)==TRUE)
                           pkg <<- as.name(pkg)
                           stopifnot(length(func)==2)
                 SimModel <- as.function(func[[1]])
                 FitModel <- as.function(func[[2]])
                         }
      }
      if (TestType == "1") {
          (OneMonteCarlo <- function() {
                if (InfiniteVarianceQ) 
                  rboot <- ts(rStable(n, ALPHA, BETA, GAMMA,DELTA))
                else
                  rboot <- res[sample(x=1:n,size=n,replace=TRUE, prob = NULL),]
              if (test =="PenaRodriguez")
              {
                if (Kernel==FALSE){  
                OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                  2]
                }
                else{
                  OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[, 
                                                                      2]
                }
              }
              else if (test =="BoxPierce")
                  OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                  2]
              else if (test =="LjungBox")
                  OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                  2]
              else if (test =="Hosking")
                  OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                  2]
              else if (test =="LiMcLeod")
                  OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                  2]
              return(OneSim.stat)
          })
      }
      else {
          (OneMonteCarlo <- function() {
                if (Model == 1) {
                  if (InfiniteVarianceQ) 
                     innov <- ts(rStable(n, ALPHA, BETA, GAMMA,DELTA))
                  else
                     innov <- sample(x=res,size=n,replace = TRUE, prob = NULL)
                  Sim.Data <- constant + stats::arima.sim(n = n, 
                    list(ar = phi, ma = theta), innov = innov, sd = sqrt(sigma), 
                    mean = demean)
                  FitSimModel <- stats::ar(Sim.Data, aic = FALSE, demean = Demean, order.max = p,method="yule-walker")
                  rboot <- FitSimModel$resid[-(1:p)]
                  if (test =="PenaRodriguez"){
                    if (Kernel==FALSE){
                    OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                    else
                    {
                      OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]
                    }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
                else if (Model == 2) {
                  if (InfiniteVarianceQ) 
                    Sim.Data <- varima.sim(phi = phi, theta = theta, sigma = sigma, 
                       n = n, constant = constant, trend = trend, demean = demean,
                       innov.dist="stable", StableParameters=StableParameters)
                  else{
                     innov <- res[sample(x=(1:(n-p)),size=n,replace=TRUE, prob = NULL),]
                     Sim.Data <- varima.sim(phi = phi, theta = theta, sigma = sigma,
                       n = n, constant = constant, trend = trend, demean = demean,
                       innov = innov)
                   }
                  FitSimModel <- stats::ar.ols(Sim.Data, aic = FALSE, demean =Demean, intercept=intercept,
                    order.max = p)
                  rboot <- ts(as.matrix(FitSimModel$resid)[-(1:p),])
                  if (test =="PenaRodriguez"){
                    if (Kernel==FALSE){
                    OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                    else {
                      OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]   
                    }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
                else if (Model == "3A" || Model == "3B" || Model == 
                  "3C" || Model == "3D") {
                  if (InfiniteVarianceQ) 
                    Sim.Data <- varima.sim(phi = phi, theta = theta, sigma = sigma, 
                      n = n, constant = constant, trend = trend, demean = demean, 
                      innov.dist="stable", StableParameters=StableParameters)
                  else{
                     innov <- res[sample(x=1:n,size=n,replace=TRUE, prob = NULL),]
                     Sim.Data <- varima.sim(phi = phi, theta = theta, sigma = sigma, 
                      n = n, constant = constant, trend = trend, demean = demean, 
                      innov = innov)
                  }
                  if (Model == "3A") 
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "none")
                  else if (Model == "3B") 
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "const")
                  else if (Model == "3C") 
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "trend")
                  else if (Model == "3D") 
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "both")
                  rboot <- resid(FitSimModel)
                  if (test =="PenaRodriguez")
                  {
                    if (Kernel==FALSE){ 
                    OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                    else 
                    {
                      OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]
                    }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
                else if (Model == 4) {
                     if (InfiniteVarianceQ) 
                       innov <- ts(rStable(n, ALPHA, BETA, GAMMA,DELTA))
                     else
                       innov <- sample(x=res,size=n,replace = TRUE, prob = NULL)
                      Sim.Data <- constant + stats::arima.sim(n = n, list(ar = phi, ma = theta), innov = innov,
                               sd = sqrt(sigma), mean = demean)
                  if (d>0)
                      Sim.Data <- stats::diffinv(Sim.Data,differences = d)[-(1:d)]
                  if (constant !=0)
                      FitSimModel <- forecast::Arima(Sim.Data, order = c(p,d,q),include.drift = TRUE)                       
                  else                             
                     FitSimModel <- stats::arima(Sim.Data, order = c(p,d,q))
                  rboot <- FitSimModel$resid
                  if (test =="PenaRodriguez")
                  {
                    if (Kernel==FALSE){
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                   else {
                     OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]
                   }
                  
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
                else if (Model == 5) {
                  if (InfiniteVarianceQ) 
                     innov <- ts(rStable(n, ALPHA, BETA, GAMMA,DELTA))
                  else
                     innov <- sample(x=res,size=n,replace = TRUE, prob = NULL)
                  Sim.Data <- stats::arima.sim(n = n, list(ar = phi, 
                    ma = theta), innov = innov, sd = sqrt(sigma), mean = demean)
                  pvec <- obj$pvec
                  if (obj$SubsetQ) {
                    if (obj$ARModel=="ARz")
                        rboot<-ts(FitAR::GetFitARz(Sim.Data, pvec=pvec, MeanValue=mean(Sim.Data))$res)
                    else
                        rboot<-ts(FitAR::GetFitARpLS(Sim.Data, pvec=pvec)$res) #let function do mean-correction
                        }
                  else #must be AR(p)
                    rboot<-ts(FitAR::GetFitARpLS(Sim.Data, pvec=pvec)$res)
                  if (test =="PenaRodriguez")
                  {
                    if (Kernel==FALSE){
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                   else {
                     OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]
                   }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
                else if (Model == 6) {
                  Sim.Data <- demean + FGN::SimulateFGN(n, H) * sigma
                  FitSimModel <- FGN::FitFGN(Sim.Data)
                  rboot <- ts(FitSimModel$res)
                  if (test =="PenaRodriguez")
                  {
                    if (Kernel==FALSE){
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                    else{
                      OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]   
                    }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
                else if (Model == 7) {
                 spec <- fGarch::garchSpec(model = list(mu=demean,alpha = Alpha, beta = Beta),cond.dist = cond.dist)
                  Sim.Data <- fGarch::garchSim(spec, n = n)
                  GARCH <- capture.output({
                        FitSimModel <- tseries::garch(Sim.Data,order = c(q2,p2))
                  })
                  rboot <- ts(FitSimModel$residuals[-(1:p2)])
                  if (test =="PenaRodriguez"){
                   if (Kernel==FALSE){
                    OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                   }
                   else {
                     OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]
                   }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
                else if (Model == 8) {
                  spec <- fGarch::garchSpec(model=list(mu=demean,omega=OMEGA,ar=phi,ma=theta,
                          alpha=Alpha,beta=Beta,delta=Delta,skew=SKEW,shape=SHAPE,leverage=Leverage),
                          cond.dist =cond.dist)
                  Sim.Data <- fGarch::garchSim(spec, n = n)
                  GARCH <- capture.output({
                  form <- as.formula(paste("~arma(", p, ",", q, ")+garch(", p2, ",", q2, ")"))
                   FitSimModel <- garchFit(formula = form, data = Sim.Data)
                  })
                  rboot <- ts(residuals(FitSimModel, standardize = TRUE))
                  if (test =="PenaRodriguez")
                  {
                    if (Kernel==FALSE){
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                   else {
                     OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2] 
                   }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]
                  return(OneSim.stat)
                }
               else if (Model == 9) {
                  Sim.Data <- SimModel(obj)
                  FitSimModel <- FitModel(Sim.Data)
                  rboot <- ts(FitSimModel$res) 
                  if (test =="PenaRodriguez")
                  {
                    if (Kernel==FALSE){
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[,2]
                    }
                    else {
                      OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ,Kernel=TRUE)[,2]  
                    }
                  }
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[,2]               
                  return(OneSim.stat)
               }
            })
       }
       if (nworkers == 1) {
        if (SetSeed) 
          set.seed(21597341)
         sim.stat <- replicate(NREP, OneMonteCarlo())
       }
       else {
            OneMonteCarlo <<- OneMonteCarlo
         if (all(class(obj) == "list")){
            SimModel <<- SimModel 
            FitModel <<- FitModel 
         }
            cl <- makeCluster(nworkers)
            if (SetSeed) 
              clusterSetRNGStream(cl, 21597341)
            if (all(class(obj) == "garch")||all(class(obj) == "fGARCH")){ 
                clusterEvalQ(cl, library("tseries"))
                clusterEvalQ(cl, library("fGarch"))
            }
            if (all(class(obj) == "list")){
              clusterExport(cl, list("SimModel","FitModel"))
              Package <- as.call(list(library, as.name(pkg)))
              clusterCall(cl,eval,Package,env = .GlobalEnv)
            }
            clusterExport(cl, list("GetResiduals","gvtest","LjungBox",
                "BoxPierce","Hosking","LiMcLeod","ImpulseVMA","InvertQ", 
                "OneMonteCarlo", "varima.sim","vma.sim", "ToeplitzBlock"))
           if (InfiniteVarianceQ){ 
                clusterEvalQ(cl, library("akima"))
                clusterExport(cl, list("fitstable", "rStable"))
           }
            sim.stat <- parSapply(cl, 1:NREP, function(j) (OneMonteCarlo()))
            stopCluster(cl)
       }
       pvalue <- numeric(length(lags))
        for (i in 1:length(lags)) {
            if (is.matrix(sim.stat)) 
                pvalue[i] <- (1 + sum(as.numeric(sim.stat[i, 
                  ] >= obs.stat[i])))/(NREP + 1)
            else pvalue[i] <- (1 + sum(as.numeric(sim.stat[i] >= 
                obs.stat[i])))/(NREP + 1)
       }
      summary <- matrix(c(lags,ans[,2],ans[,3],pvalue),ncol=4)
    dimnames(summary) <- list(rep("", length(lags)),c("Lags","Statistic","df","pvalue"))
  return(summary)
    }
}