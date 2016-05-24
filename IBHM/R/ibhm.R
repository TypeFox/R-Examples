source('R/scals.R')
source('R/activs.R')
source('R/optim.R')
source('R/stopCriteria.R')
source('R/config.R')

loadModule("evaluator",what=T)

# IBHM method -------------------------------------------------------------

TrainIBHM <- function(x,y, config = ConfigureIBHM()  ){
  if(config$jit){    
    enableJIT(3)
  }  
  
  ctx <- CreateContext(x,y,config)
  
  ctx$m$par$w0 <- mean(y)  
  ctx$yr <- ctx$y-ctx$m$par$w0
  
  while(ctx$continue){
    ctx$iteration <- ctx$iteration + 1
  
    ctx <- ctx$sc$eval( DoIteration(ctx) )    
  }  
    
  
  ctx$log('Final parameter estimation ')
  ctx <- config$final.estimation(ctx, x = ctx$x.final, y = ctx$y.final, maxit = config$final.estimation.maxit)
  
  ctx$log('\nFinal RMSE: ',CalcRMSE(ctx))
  
  return(ctx$m)
}


DoIteration <- function(ctx){
  ctx$log('\nIteration ',ctx$iteration, ' RMSE: ', CalcRMSE(ctx))    
  
  ctx <- FindActiv( FindScal(ctx) ) 
  ctx$m$par$w <- append(ctx$m$par$w,1)  
  ctx <- OptimizeWeights(ctx)
  UpdateResidual(ctx)
}

CalcRMSE <- function(ctx){
  sqrt(mean((ctx$y-predict.IBHM(ctx$m,ctx$x))^2))
}



# Method runtime context --------------------------------------------------

CreateContext <- function(x,y, params){
  x <- as.matrix(x)
  y <- as.matrix(y)
  stopifnot(nrow(x)==nrow(y), nrow(x)>10, is.null(params$final.estimation.x) == is.null(params$final.estimation.y))  
  
  
  if(params$verbose){
    log <- function(...) cat(...,'\n')
  }else{
    log <- function(...) {}
  }
  
  x.m <- max(abs(x))
  params$scal$optim.params$par.min <- min(-1,-x.m)
  params$scal$optim.params$par.max <- max(1, x.m)
  
  params$activ$optim.params$par.min <- min(-1,-x.m)
  params$activ$optim.params$par.max <- max(1,x.m)
  
  append(params,
         list( x = x,
               y = y, 
               yr = y,               
               x.final = if(is.null(params$final.estimation.x)) x else params$final.estimation.x,
               y.final = if(is.null(params$final.estimation.y)) y else params$final.estimation.y, 
               m=CreateEmptyModel(x,y), 
               scal.y = list(),
               scal.evaluator = new(ScalEvaluator, y),
               continue=TRUE, 
               iteration=0,
               log = log
         )
  )
}


# IBHM model structure ----------------------------------------------------
CreateEmptyModel <- function(x=NULL,y=NULL){
  m <- list(params = list(
                w0=0,
                w=vector(),
                d=list(),
                a=vector(),
                b=vector()),
            scals=list(),            
            activs=list(),            
            train=list(x=x,y=y))
  class(m) <- 'IBHM'
  return(m)
}



# Standard S3 methods -----------------------------------------------------
predict.IBHM <- function(object,x = NULL,...){  
  m<-object
  if(is.null(x)){
    x <- m$train$x
  }else{
    x <- as.matrix(x)  
  }
  
  
  d <- dim(x)
  len <- length(m$scals)  
  y <- matrix(m$par$w0,d[[1]],1)
  
  
  if(len > 0){
    for(i in 1:len){
      act <- m$activs[[i]]
      a <- m$par$a[[i]]
      b <- m$par$b[[i]]
            
      scal <- m$scals[[i]]
      d <- m$par$d[[i]]
      
      w <- m$par$w[[i]]    
      
      y <- y + w * act( a*scal(x,d)+b)
    }
  }
  
  y
}


length.IBHM <- function(x){
  length(x$par$w)
}

ToString <- function(m){  
  f <- function(v){
    sprintf('%.2e',v)
  }
  
  str <- f(m$par$w0)
  len <- length(m$par$w)
  
  if(len>0){  
    for(i in 1:length(m)){
      w <- m$par$w[[i]]
      if(w>0) s<-'+' else s<-''
      
      act <- expr(m$activs[[i]])
      a <- f(m$par$a[[i]])
      b <- f(m$par$b[[i]])
      
      scal <- expr(m$scals[[i]])
      d <- paste(f(m$par$d[[i]]),collapse=' ')
      
      
      str <- paste(str,s,f(w),act,' (',a,'*',scal,'(x,[',d,'])','+',b,') ')
    }
  }
  
  str
  
}

summary.IBHM <- function(object,...){  
  y <- object$train$y
  x <- object$train$x
  
  se <- c((y - predict(object))^2)
  
  res <- list(model = ToString(object),
              model.size = length(object),
              TrainSize = nrow(y),
              TrainDim = ncol(x),
              mse = mean(se),
              sd.se = sd(se),
              rmse = sqrt(mean(se)),
              cor=cor(predict(object),y))
  class(res) <- 'summary.IBHM'
  res
}

print.summary.IBHM <- function(x,...){
  cat('Model equation: ', x$model,'\n',
      'Model size: ',x$model.size,'\n\n',
      'Train set dim: ', x$TrainDim, ' Train set size: ', x$TrainSize,'\n',      
      'MSE:  ', x$mse, ' Std. dev:', x$se.sd,'\n',
      'RMSE: ', x$rmse, '\n',
      'Pearson correlation coefficient: ', x$cor,'\n')
}


# Finding a scalarization function ----------------------------------------

FindScal <- function(ctx){
  ctx$log(' Finding next scal')
  
  nScals <- length(ctx$scal$candidates)  
  
  best.eval <- Inf
  best <- NULL
  
  for( i in  1:nScals){
    scal <- ctx$scal$candidates[[i]]    
    candidate <- OptimizeScal(scal,ctx)
    ctx$log('  ',expr(scal), ' eval : ', candidate$eval)
    
    if(candidate$eval < best.eval){
      best <- candidate
      best.eval <- candidate$eval
    }
  }
  
  ctx$log(' Best scal: ',expr(best$scal),' w.par:',best$w.par)
  
  k <- ctx$iteration
  ctx$m$scals[[k]]  <- best$scal
  ctx$m$par$d[[k]] <- best$d
  ctx$w <- ctx$wf(best$scal(ctx$x, best$d), best$w.par)   
  ctx$scal.evaluator$addScal(best$scal(ctx$x, best$d))
  
  ctx
}


OptimizeScal <- function(scal, ctx){
  
  goal <- function(par){              
    w.par <- par[[1]]
    zCandidate <- scal(ctx$x, par[-1])
    w <- ctx$wf(zCandidate, w.par)
              
    return(1 - ctx$scal.evaluator$evaluate(zCandidate, w))
  }
  n <- ScalParamDim(scal,ctx$x) + 1  
  res <- ctx$scal$optim(goal,n,ctx$scal$optim.params)    
  
  append(res, list(scal = scal))
}

EvaluateScal <- function(y,yh,w.par,ctx){
  w <- ctx$wf(yh, w.par)
  y.var <- weightedVar(y,w)
  yh.var <- weightedVar(yh,w) 
  
  
  if( is.finite(y.var) && y.var > 1e-5 && is.finite(yh) && yh > 1e-5){      
    indep <- 1.0
    for(ys in ctx$scal.y){          
      indep <- indep * (1-abs(weightedRho(yh,ys,w)))
    }
    
    
    cor.y <- weightedRho(y, yh,w)    
    
    if(is.nan(cor.y)){
      cor.y<-0
    }
    if(is.nan(indep)){
      indep <- 1.0  
    }
  
    return(abs(cor.y) * indep)    
  }else{
    return(0.0)
  }
}

# Finding an activation function #####################


FindActiv <- function(ctx){
  ctx$log(' Finding next activ')    
  
  yh <- RunLastScal(ctx)
  best <- NULL
  best.eval <- Inf
    
  for(i in 1:length(ctx$activ$candidates)){
    cand <- OptimizeActiv(ctx$activ$candidates[[i]], yh, ctx)        
    if(cand$eval < best.eval){
      best <- cand      
      best.eval <- cand$eval
    }
  }    
  
    
  k <- ctx$iteration  
  
  best.output  <- best$activ(best$a*yh+best$b)
  best.output.sd <- sd(best.output)
  
  ctx$m$activs[[k]] <- best$activ
  ctx$m$par$a[[k]]      <- best$a
  ctx$m$par$b[[k]]      <- best$b
  
  ctx$log(' Best activ: ', expr(best$activ), ' a:',best$a, ' b:',best$b,' alpha:',ctx$m$alpha[[k]],' beta:',ctx$m$beta[[k]])
  
  ctx
}


OptimizeActiv <- function(activ, yh, ctx){  
  goal <- function(par){          
    a <- par[[1]]
    b <- par[[2]]
    
    ya <- activ(a*yh+b)    
    
    eval <- 1-abs(EvaluateActiv(ctx$y,ya,ctx$w))     
    
    if(is.nan(eval)){  
      eval <- Inf        
    } 
    return(eval)
  }  
  
  res <- ctx$activ$optim(goal, ctx$activ$optim.params)  
      
  ctx$log('  ',expr(activ),' eval:',res$eval)
  append(res, list(activ=activ))
}

EvaluateActiv <- function(y,ya,w){
  val<-weightedR(y, ya,w)  
  if(is.nan(val)){
    val<-0
  }
  if(!is.finite(val)){
    val <- 0
  }
  return(val)
}

RunLastScal <-function(ctx){
  lastScalIdx <- length(ctx$m$scals)
  scal <- ctx$m$scals[[lastScalIdx]]  
  d <- ctx$m$par$d[[lastScalIdx]]
  
  scal(ctx$x,d)  
}


# Estimating regression weights -------------------------------------------

OptimizeWeights <- function(ctx){  
  m <- ctx$m
  w.idx <- length(m$par$w)
  
  goal <- function(w){            
    m$par$w[[w.idx]] <- w
    
    return(weightedMean((ctx$y - predict.IBHM(m,ctx$x))^2,ctx$w))            
  }
  
  res <- optim(1,goal, method="BFGS")  
  m$par$w[[w.idx]] <- res$par
  
  ctx$m <-m  
  ctx
}

OptimizeAllWeights <- function(ctx, x = ctx$x, y = ctx$y, maxit=100){  
  ctx$log(' Optimizing all weights')
  m <- ctx$m
  goal <- function(par){    
    m$par$w0 <- par[[1]]
    m$par$w <- par[2:length(par)]
    
    return(mean((y - predict.IBHM(m, x))^2))
  }
  
  res <- optim(c(m$par$w0,m$par$w),goal, method="BFGS", control=list(maxit=maxit))
  m$par$w0 <- res$par[[1]]
  m$par$w <- res$par[2:length(res$par)]
  
  ctx$m <-m
  
  ctx
}

OptimizeAllParams <- function(ctx,x = ctx$x, y = ctx$y, maxit=200){
  ctx$log(' Optimizing all parameters')
  m <- ctx$m
  goal <- function(par){    
    
    m$par <- relist(par, m$par)        
    return(mean((y - predict.IBHM(m,x))^2))
  }
            
  res <- optim(unlist(m$par),goal, method="BFGS", control=list(maxit=maxit, trace=F))  
  ctx$m$par <- relist(res$par, m$par)  
  ctx
}

UpdateResidual <- function(ctx){
  ctx$yr <- ctx$y - predict.IBHM(ctx$m, ctx$x)
  ctx$scal.evaluator$updateResidual(ctx$yr)
  ctx
}





