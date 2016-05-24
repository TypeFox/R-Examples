max_qEI <- function(model, npoints, lower, upper, crit="exact", minimization = TRUE, optimcontrol=NULL) {  
  if (is.null(optimcontrol$method)) optimcontrol$method <- "BFGS"
  optim.method <- optimcontrol$method
  d <- model@d 
  parinit <- optimcontrol$parinit
   if (crit == "CL") {
    res <- max_qEI.CL(model, npoints, optimcontrol$L, lower, upper, parinit=parinit,minimization = minimization, control=optimcontrol)
    res <- list(par = res$par, value = matrix(qEI(x = res$par, model = model)))
   } else if(crit=="exact") {
     EI.envir <- new.env()
     environment(qEI) <- environment(qEI.grad) <- EI.envir
     
     LOWER <- c(apply(matrix(lower,d,1),1,rep,npoints))
     UPPER <- c(apply(matrix(upper,d,1),1,rep,npoints))
      if (optim.method == "BFGS") {
        if (is.null(parinit))
          parinit <- array(NaN,dim=c(npoints,model@d,0))
        if (is.null(optimcontrol$nStarts))
          optimcontrol$nStarts <- 4
        if (is.null(optimcontrol$fastCompute))
          optimcontrol$fastCompute <- TRUE
        if (is.null(optimcontrol$sampleFun))
          optimcontrol$sampleFun <- sampleFromEI
        if (is.null(optimcontrol$gradNum))
          optimcontrol$gradNum <- FALSE
        if (is.null(optimcontrol$maxit))
          optimcontrol$maxit <- 100
        maxi <- 0
        startPoints <- optimcontrol$sampleFun(model=model,minimization=TRUE,n = (npoints*optimcontrol$nStarts),lower=lower,upper=upper)
        for (i in 1:(optimcontrol$nStarts)) {
          if (i > length(parinit[1,1,])){
            x <- startPoints[((i-1)*npoints+1):(i*(npoints)),]
          } else {
            x <- parinit[,,i]
          }
          if (!(optimcontrol$gradNum)) {
            o <- optim(c(x), qEI ,qEI.grad, control=list(trace = 0,REPORT=1,fnscale = -1,maxit=optimcontrol$maxit),method="L-BFGS-B",lower=LOWER,upper=UPPER,model = model,envir=EI.envir,fastCompute=optimcontrol$fastCompute,minimization = minimization)
          } else {
            o <- optim(par = c(x), fn = qEI ,gr=NULL, control=list(trace = 0,REPORT=1,fnscale = -1,maxit=optimcontrol$maxit),method="L-BFGS-B",lower=LOWER,upper=UPPER,model = model,envir=EI.envir,fastCompute=optimcontrol$fastCompute,minimization=minimization)
          }
          par <- matrix(o$par,ncol=d)
          value <- as.matrix(o$value)
          if (value>=maxi) {
            maxi <- value
            parmax <- par
          }
        }
        colnames(parmax) <- colnames(model@X)
        colnames(maxi) <- "EI"
        res <- list(par = parmax, value = maxi)
        #res <- max_qEI.BFGS(model,npoints,lower,upper, parinit = NULL,minimization = minimization, control = control)
      } else if (optim.method =="genoud") {
        if (is.null(optimcontrol$pop.size)) optimcontrol$pop.size <- 50*d
        if (is.null(optimcontrol$max.generations)) optimcontrol$max.generations <- 5
        if (is.null(optimcontrol$wait.generations)) optimcontrol$wait.generations <- 2
        if (is.null(optimcontrol$BFGSburnin)) optimcontrol$BFGSburnin <- 2
        if (is.null(optimcontrol$fastCompute)) optimcontrol$fastCompute <- TRUE
        if (is.null(optimcontrol$print.level)) optimcontrol$print.level <- 0
        domaine <- cbind(LOWER,UPPER)
        
        o <- genoud(qEI, nvars = (d*npoints), max = TRUE,
                    gr = qEI.grad,
                    pop.size = optimcontrol$pop.size, 
                    max.generations = optimcontrol$max.generations, 
                    wait.generations = optimcontrol$wait.generations, 
                    hard.generation.limit = TRUE, starting.values = c(optimcontrol$parinit), 
                    MemoryMatrix = TRUE, Domains = domaine, default.domains = 10, 
                    solution.tolerance = 1e-09, boundary.enforcement = 2, 
                    lexical = FALSE, gradient.check = FALSE, BFGS = TRUE, 
                    data.type.int = FALSE, hessian = FALSE, 
                    unif.seed = floor(runif(1, max = 10000)), 
                    int.seed = floor(runif(1, max = 10000)),
                    print.level = 0, share.type = 0, instance.number = 0, 
                    output.path = "stdout", output.append = FALSE, project.path = NULL, 
                    P1 = 50, P2 = 50, P3 = 50, P4 = 50, P5 = 50, P6 = 50, 
                    P7 = 50, P8 = 50, P9 = 0, P9mix = NULL, BFGSburnin = optimcontrol$BFGSburnin, 
                    BFGSfn = NULL, BFGShelp = NULL, cluster = FALSE, balance = FALSE, 
                    debug = FALSE, model = model,fastCompute = optimcontrol$fastCompute,minimization=minimization,envir=EI.envir)
        
        o$par <- matrix(o$par,ncol=d)
        colnames(o$par) <- colnames(model@X)
        o$value <- as.matrix(o$value)
        colnames(o$value) <- "EI"
        res <- list(par = o$par, value = o$value)
        #res <- max_qEI.gen(model,npoints,lower,upper, parinit = NULL,minimization = minimization, control = control)
      } else {
        stop(paste(paste("\'",optim.method,"\'",sep=""),"is unknown."))
      }
   } else {
     stop(paste(paste("\'",crit,"\'",sep=""),"is unknown."))
   }
   return(res)
}



max_qEI.CL <- function(model, npoints, L, lower, upper, parinit=NULL, minimization = TRUE, control=NULL) {
  KB <- FALSE
  n1 <- nrow(model@X)
  lengthParinit <- length(parinit[1,])
  if (is.null(L)) {
    liar <- minimization*min(model@y) + (!minimization)*max(model@y)
  } else if (L == "max") {
    liar <- max(model@y)
  } else if (L == "min") {
    liar <- min(model@y)
  } else if (L == "mean") {
    KB <- TRUE
  } else {
    if ((!is.numeric(L))||(length(L)!=1)) stop("control$L must be NULL, \"max\", \"min\", \"mean\" or a scalar specifying the plugin.")
    liar <- L
  }
  for (s in 1:npoints) {
    if (s > lengthParinit){
      startPoint <- NULL
    } else {
      startPoint <- parinit[,s]
    }
    oEGO <- max_EI(model=model, lower=lower, upper=upper, parinit=startPoint, minimization = minimization, control=c(control,list(print.level=0)))
  		model@X <- rbind(model@X, oEGO$par)
    if (KB) liar <- predict(object=model,newdata=oEGO$par,se.compute=FALSE,cov.compute=FALSE,light.return=TRUE,checkNames=FALSE)$mean
    model@y <- rbind(model@y, liar, deparse.level=0)   		
    model@F <- trendMatrix.update(model, Xnew=data.frame(oEGO$par))	
    if (model@noise.flag) {
	# heterogenous case : use 0 nugget for new points
	model@noise.var = c(model@covariance@nugget, 0)
    }	
    model <- computeAuxVariables(model)
  }
  return(list(par = model@X[(n1+1):(n1+npoints),, drop=FALSE], value = model@y[(n1+1):(n1+npoints),, drop=FALSE]))
}
