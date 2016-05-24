
sparse.tscgm <- function(data=data, lam1=NULL, lam2=NULL, nlambda=NULL, model=c("ar1","ar2"),
            optimality=c("NULL","bic","bic_ext","bic_mod","aic","gic"),
            control=list())
{
   if(is.longitudinal(data)==TRUE)
  {
   n_obs = get.time.repeats(data)$repeats[[1]]
   tps = get.time.repeats(data)$time
   p_num = dim(data)[2]
   time = dim(data)[1]/n_obs
   xy.data <-array(NA, c(time,p_num,n_obs))
   data <- as.matrix(data)
   for(i in 1:n_obs){
     for(t in 1:time){
       cc <- 1+(t-1)*n_obs + (i-1)
       xy.data[t,,i] <- data[cc,]
     }
   }
 }
 else cat("Data format is not longitudinal.", "\n")

  model = match.arg(model)

  if (model=="ar1") {
   X <- xy.data[1:time-1,,,drop=FALSE]
   Y <- xy.data[2:time,,,drop=FALSE]
      }
  else if (model=="ar2") {
    t1=time-1
    t2=time-2
    Y <- round(xy.data[3:time,,,drop=FALSE],3)
    X1 <- round(xy.data[2:t1,,,drop=FALSE],3)
    X2<- round(xy.data[1:t2,,,drop=FALSE],3)
    X <- abind(X1, X2,along = 2 )

  }
  T <- dim(Y)[1]
  p <- dim(X)[2]
  n <- dim(Y)[3]
  q <- dim(Y)[2]
  xtyi <- array(NA, c(p,q,n))
  xtxi <- array(NA, c(p,p,n))
  ytyi <- array(NA, c(q,q,n))

  for(i in 1:n){
      XX <- X[,,i]
      YY <- Y[,,i]
      XX2 <- X[,,i]^2
      YY2 <- Y[,,i]^2
      xtyi[,,i]=crossprod(XX,YY)
      xtxi[,,i]=crossprod(XX)
      ytyi[,,i]=crossprod(YY)

    }
  xty=apply(xtyi, c(1,2), sum)
  xtx=apply(xtxi, c(1,2), sum)
  yty=apply(ytyi, c(1,2), sum)
  xtxt=apply(xtxi, c(1,2), sum)/(n*T)
  xtx2=(n*T)*colMeans(apply(XX2, c(1,2), sum))
  yty2=(n*T)*colMeans(apply(YY2, c(1,2), sum))

  if(is.null(lam1) | is.null(lam2) )  {
      SX <- xtx/(n*T)
      mSX <- glasso(SX,0.05,penalize.diagonal=FALSE)
      if(model == "ar1") {
          SX <- xtx/(n*T)
          mSX <- glasso(SX,0.05,penalize.diagonal=FALSE)
          SXi <- mSX$wi
          SS =(yty)/(n*T)
          SS = cov2cor(SS)
          SAs = xty/(n*T)
          SA = SAs %*% SXi
         }
      if(model == "ar2") {
          SX <- xtx/(n*T)
          mSX <- glasso(SX,0.05,penalize.diagonal=FALSE)
          SXi <- mSX$wi
          SS =(yty)/(n*T)
          SS = cov2cor(SS)
          SAs = xty/(n*T)
          SA =  SXi %*% SAs
         }
      lambda <-  lambda.seq(SS=SS,SA=SA, nlambda=nlambda)
      lam1 <- round(lambda$lambda1,3) 
      lam2 <- round(lambda$lambda2,3)
      lam2 <- round(lam2/max(lam2),3)    
    }

    optimality = match.arg(optimality)
    nobic = (length(lam1) + length(lam2) == 2)
    doms=(length(lam1)+length(lam2) > 2)
    if(!is.list(control))
	stop("control is not a list")
	setting <- list(maxit.out = 5, maxit.in = 50, tol.out = 1e-04, silent = TRUE)
	nmsSetting <- names(setting)
	setting[(nms <- names(control))] <- control
	if(length(noNms <- nms[!nms %in% nmsSetting]))
	warning("unknow names in control: ",paste(noNms,collapse =", "))
 	if(nobic != 2 & optimality == "nosel" )
	stop("Specify positive scalar values for the tuning parameters")
	if(doms > 2 & optimality != "nosel")
	stop("Specify vector of positive decreasing values for the tuning parameters")
	if(setting$maxit.out < 1 )
	stop("read the documentation for 'maxit.out' more carefully")
	if(setting$maxit.in < 1 )
	stop("read the documentation for 'maxit.in' more carefully")
	if(setting$tol.out <= 0)
	stop("read the documentation for 'tol.out' more carefully")
	if(!setting$silent %in% c("TRUE","FALSE"))
	stop("read the documentation for 'silent' more carefully")

    doNULL = nobic & (optimality == "NULL")
    dobic =  doms & (optimality == "bic")
    dobic1 =  doms & (optimality == "bic_ext")
    dobic2 =  doms & (optimality == "bic_mod")
    dobic3 =  doms & (optimality == "aic")
    dobic4 =  doms & (optimality == "gic")
    gamma = NULL
    theta = NULL


    if (doNULL) {
      lambda1 = NULL
      lambda2 = NULL

        tmp.out = compute.sparse.tscgm(X=X, Y=Y, lam1 = lam1, lam2 = lam2,
        optimality = "NULL", setting=setting)

      }
    else if (dobic) {
        tmp.out = sparse.tscgm.bic(X=X, Y=Y, lam1 = lam1, lam2 = lam2,
        optimality="bic", setting=setting)
    }
    else if (dobic1) {
        tmp.out = sparse.tscgm.bic(X=X, Y=Y, lam1 = lam1, lam2 = lam2,
        optimality="bic_ext", setting=setting)
    }
    else if (dobic2) {
        tmp.out = sparse.tscgm.bic(X=X, Y=Y, lam1 = lam1, lam2 = lam2,
         optimality="bic_mod", setting=setting)
    }
    else if (dobic3) {
        tmp.out = sparse.tscgm.bic(X=X, Y=Y, lam1 = lam1, lam2 = lam2,
        optimality="aic", setting=setting)
    }
    else if (dobic4) {
        tmp.out = sparse.tscgm.bic(X=X, Y=Y, lam1 = lam1, lam2 = lam2,
        optimality="gic", setting=setting)
    }
    gamma = tmp.out$gamma
    gamma = gamma*(1*(abs(gamma) > 0.01))

    theta = tmp.out$theta
    theta = theta*(1*(abs(theta) > 0.01))

    lam1.opt = tmp.out$lam1.opt
    lam2.opt = tmp.out$lam2.opt
    lam1.seq = tmp.out$lam1.seq
    lam2.seq = tmp.out$lam2.seq
    s.gamma  =  tmp.out$s.gamma
    s.theta = tmp.out$s.theta
    tun.ic = tmp.out$tun.ic
    min.ic = tmp.out$min.ic
    if(model=="ar1") {colnames(gamma) <- colnames(data)}
    else if(model=="ar2"){
     colnames(gamma) <- colnames(data)
     rownames(gamma) <- c(colnames(data), colnames(data))
     }
	   colnames(theta) <- rownames(theta) <- colnames(data)
	   ###
    return(list(gamma = gamma, theta = theta, lam1.opt=lam1.opt, lam2.opt=lam2.opt, lam1.seq=lam1.seq, lam2.seq=lam2.seq,
     min.ic=min.ic, tun.ic=tun.ic, s.gamma=s.gamma,  s.theta=s.theta))

}

