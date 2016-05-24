optim.relatedness <- function(obs,theta0=0,theta1=0.03,theta.tol=10^(-7),theta.step=NULL,max.bisect=15,
                              probs,var.list=NULL,init.alpha=10^c(-4,-6,-8,-10),init.keep=FALSE,
                              objFunction=c("T2","T1","C3","C2","C1"),collapse=FALSE,trace=FALSE,
                              solnp.ctrl=list(tol=10^(-9),rho=10,delta=min(init.alpha)*0.01,trace=FALSE)){
  if(!any(objFunction==c("T2","T1","C3","C2","C1"))) stop("Wrong objFunction is supplied. Use any of: \"T2\", \"T1\", \"C3\", \"C2\" or \"C1\")")
  if(sum(init.alpha)>1 | any(init.alpha<0) | any(init.alpha>1)) stop("init.alpha is a probability vector, i.e. 0<=alpha[i]<=1 and sum(alpha)<=1")
  objFunction <- objFunction[1]
  n <- (1+sqrt(1+8*sum(obs,na.rm=TRUE)))/2
  if(is.matrix(obs)) obs <- t(obs)[up.tri(obs)]
  if(collapse){ ## Collapse the M_{m/p} matrix into matrix of matching alleles
    warning("attempt to 'collapse' observations into mathcing alleles ignored. The options is not yet implemented.")
  }
  
  ## Declaration of the different objective functions
  C1 <- function(x,expected,observed,variance=NULL){
    EE <- x[1]*expected[[1]]
    for(r in 2:length(expected)) EE <- EE + x[r]*expected[[r]]
    sum(((observed-EE)^2/EE),na.rm=TRUE)
  }
  C2 <- function(x,expected,observed,variance=NULL){
    EE <- x[1]*expected[[1]]
    for(r in 2:length(expected)) EE <- EE + x[r]*expected[[r]]
    sqrt(sum((observed-EE)^2,na.rm=TRUE))
  }
  C3 <- function(x,expected,observed,variance=NULL){
    EE <- x[1]*expected[[1]]
    for(r in 2:length(expected)) EE <- EE + x[r]*expected[[r]]
    sum(abs(observed-EE)/EE,na.rm=TRUE)
  }
  T1 <- function(x,expected,observed,variance=NULL){
    EE <- x[1]*expected[[1]]
    for(r in 2:length(expected)) EE <- EE + x[r]*expected[[r]]
    sum((EE-observed)^2/diag(variance))
  }
  T2 <- function(x,expected,observed,variance=NULL){
    EE <- x[1]*expected[[1]]
    for(r in 2:length(expected)) EE <- EE + x[r]*expected[[r]]
    gginv <- eigen(variance,T,F,T)
    ivariance <- gginv$vec%*%diag(c(1/gginv$val[-length(gginv$val)],0))%*%t(gginv$vec)
    as.numeric(t(EE-observed)%*%ivariance%*%(EE-observed))
  }
  objFun <- get(objFunction)
  ## Now the objFun contains the function declaration for the object function needed in the optimisation steps

  tgrid <- seq(from=theta0, to=theta1, len=3) ## Initial theta grid
  val.df <- data.frame(theta=NA,value=NA)[0,] ## The output containing the theta values and 
  min.t <- rep(tgrid[2],2)

  ## GRID SEARCH
  if(!is.null(theta.step)){
    tgrid <- seq(from=theta0,to=theta1,by=theta.step)
    grid.search <- TRUE
  }
  else grid.search <- FALSE
  
  bb <- 0
  min.val <- Inf
  min.id <- 1
  if(sum(init.alpha)!=1) init.alpha <- c(1-sum(init.alpha),init.alpha)
  min.res <- init.alpha
  if(trace){
    if(grid.search) cat("Grid search... ")
    else cat("Bisectional search...\n")
  }
  while(bb < max.bisect){
    bb <- bb+1 ## number of bisection searches performed
    if(trace & !grid.search) cat(paste("Iteration: ",format(bb,width=2),"\n",sep=""))
    if(diff(range(tgrid))<theta.tol & !grid.search){
      if(trace){
        if(solnp.ctrl$trace) cat("\n")
        cat("Interval converged\n")
      }
      break
    }
    if(length(min.t)>10){
      if(abs(diff(min.t[length(min.t)-c(0,5)]))<(theta.tol)^2 & !grid.search){
        if(trace) cat("No change in theta for several iterations\n")
        break
      }
    }
    if(grepl("T",objFunction)){ ## Object function is either T2 or T1 - variances are needed..
      if(is.null(var.list)){
        if(trace) cat("Variances are being computed... Please wait")
        var.list <- dbVariance(probs=probs,theta=tgrid,n=1)
        names(var.list) <- paste(tgrid)
        if(trace) cat("  Done..!\n")
      }
      else{
        if(all(is.element(paste(tgrid),names(var.list))) & trace) cat("All needed variances are provided in the input...\n")
        else{
          if(trace) cat(paste("Missing variances (",sum(!is.element(paste(tgrid),names(var.list))),") are being computed... Please wait",sep=""))
          ttgrid <- tgrid[!is.element(paste(tgrid),names(var.list))]
          vvar.list <- dbVariance(probs,theta=ttgrid,n=1)
          if(length(ttgrid)==1) vvar.list <- list(vvar.list)
          names(vvar.list) <- ttgrid
          var.list <- c(var.list,vvar.list)
          var.list <- var.list[sort.list(names(var.list))]
          if(trace) cat("  Done..!\n")
        }
      }
      variances <- lapply(var.list,function(x,n) choose(n,2)*x$V1+6*choose(n,3)*x$V2+6*choose(n,4)*x$V3,n=n)
      variances <- variances[paste(tgrid)]
    }
    else{ ## Object function is either C1, C2 or C3
      var.list <- NULL
      variances <- replicate(length(tgrid), NULL, simplify=FALSE)
    }
    expects <- lapply(tgrid,function(t,n) list("UN"=dbExpect(probs=probs,theta=t,n=n,vector=TRUE,k=c(0,0,1)),
                                               "FC"=dbExpect(probs=probs,theta=t,n=n,vector=TRUE,k=c(0,1,3)/4),
                                               "AV"=dbExpect(probs=probs,theta=t,n=n,vector=TRUE,k=c(0,1,1)/2),
                                               "PC"=dbExpect(probs=probs,theta=t,n=n,vector=TRUE,k=c(0,1,0)),
                                               "FS"=dbExpect(probs=probs,theta=t,n=n,vector=TRUE,k=c(1,2,1)/4)),n=n)
    ## Part of the actual computations for each value of theta in the bisection or grid search 
    for(i in 1:length(tgrid)){
      t <- tgrid[i] ## current theta value under consideration
      solnpObjFun <- function(x) objFun(x,expected=expects[[i]],observed=obs,variance=variances[[i]])
      alpha <- min.res
      if(init.keep) alpha <- init.alpha
      est <- try(solnp(pars = alpha, fun=solnpObjFun, eqfun=sum, eqB=1, LB=rep(0,5), UB=rep(1,5), control=solnp.ctrl),silent=TRUE)
      if(length(est)==1){ ## 'est' contains the error message produced by the solnp function
        est <- list(pars=init.alpha,values=NA)
        val.df <- rbind(val.df,data.frame(theta=t,value=NA))
        warn.message <- paste("NAs were returned for theta =",format(t,digits=5),"by the solnp procedure. This could indicate numerical problems with the identified solution.",sep=" ")
        if(!grid.search) warn.message <- paste(warn.message,"\n  An attempt to overcome or approximate/solve the problem is to use a grid search instead using 'theta.step' to set step size.")
        warning(warn.message)
      }
      else{
        step.val <- est$values[length(est$values)]
        if(step.val<=min.val){
          min.val <- step.val
          min.res <- est$pars
          min.id <- i
          min.t <- c(min.t,t)
        }
        val.df <- rbind(val.df,data.frame(theta=t,value=step.val))
      }
    }
    ##
    if(grid.search){
      if(trace) cat("\n")
      break ## breaks out of bisectional search loop
    }
    ##
    if(min.id==1) tgrid <- seq(from=tgrid[1],to=tgrid[2],len=3)
    else if(min.id==length(tgrid)) tgrid <- seq(from=tgrid[min.id-1],to=tgrid[min.id],len=3)
    else{ ## The mid point t1 is the minima of [t0,t1,t2].. Update to search in [t0+0.5*l,t2-0.5*l], where l = t1-t0 = t2-t1
      tgrid.length <- c(tgrid[min.id]-min(tgrid),max(tgrid)-tgrid[min.id])
      tgrid <- sort(c(tgrid[1],tgrid[min.id],tgrid[length(tgrid)],tgrid[min.id]+0.25^bb*tgrid.length[1],tgrid[min.id]-0.25^bb*tgrid.length[2]))
    }
  }
  names(min.res) <- c("Unrelated","First-Cousins","Avuncular","Parent-child","Full-siblings")
  val.df <- val.df[!duplicated(val.df$theta),]
  res <- list(value=val.df[order(val.df$theta),],solution=c(theta=tgrid[min.id],min.res),var.list=var.list)
  attributes(res)$objFun <- objFunction
  attributes(res)$class <- "dbOptim"
  res
}
  
plot.dbOptim <- function(x,type="l",...){
  objFun <- attributes(x)$objFun
  ylabel <- switch(objFun, C1 = expression(C[1](theta)), C2 = expression(C[2](theta)),
                   C3 = expression(C[3](theta)), T1 = expression(T[1](theta)), T2 = expression(T[2](theta)))
  plot(value~theta,x$value,xlab=expression(theta),ylab=ylabel,type=type,...)
}

points.dbOptim <- function(x,type="p",...){
  points(value~theta,x$value,type=type,...)
}

lines.dbOptim <- function(x,type="l",...){
  points(value~theta,x$value,type=type,...)
}

print.dbOptim <- function(x,var.list=FALSE,...){
  if(var.list) print(x[c("value","solution","var.list")])
  else print(x[c("value","solution")])
}
