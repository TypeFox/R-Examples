cv.glm1path = function(object, block = NULL, best="min", plot=TRUE, prop.test=0.2, n.split = 10, seed=NULL, show.progress=FALSE, ...)
{

  tol=c(1.e-8,.Machine$double.eps)
  
  #extract some useful bits from the glm1path object
  n.rows     = length(object$y)
  lambdas    = object$lambdas
  n.lambda   = length(lambdas)

  if(is.null(seed)==FALSE)
    n.split = length(seed)
  
  #pre-set the size of some big fat uglies  
  ll.test    = array(NA, c(n.lambda,n.split), dimnames=list(lambdas,seed))
  ll.train   = array(NA, c(n.lambda,n.split), dimnames=list(lambdas,seed))
  phi.cv     = array(NA, c(n.lambda,n.split), dimnames=list(lambdas,seed))
  df.cv      = array(NA, c(n.lambda,n.split), dimnames=list(lambdas,seed))
  counter.cv = array(NA, c(n.lambda,n.split), dimnames=list(lambdas,seed))
  beta.cv    = array(NA, c(dim(object$X)[2],n.lambda,n.split), dimnames=list(dimnames(object$X)[[2]],lambdas,seed))

  # if a block argument has been specified as input, work out which obs belong in which block.
  if(is.null(block)==FALSE)
  {
    tb        = table(block)
    n.levels  = length(tb)
    blockIDs  = vector("list",n.levels)
    for(i.level in 1:n.levels)
      blockIDs[[i.level]] = which(block==names(tb)[i.level])
  }
  else
    n.levels = n.rows
  
  n.test = floor(n.levels*prop.test)
  for(i.split in 1:n.split)
  {
      if(is.null(seed)==FALSE) #set seed if it has been given
        set.seed(seed[i.split])
      is.test = sample(1:n.levels, n.test)
      if(is.null(block) == FALSE)
        is.test    = unlist(blockIDs[is.test]) # if block sampling, turn into a vector of row IDs
      y.test     = object$y[is.test]
      
      out       = glm1path(object$y[-is.test], object$X[-is.test,], lambdas=object$lambdas, family=object$family, penalty=object$penalty, ...)
      n.lambdai = length(out$logL) #in case it stops before if finishes the full range of lambdas
      ll.train[1:n.lambdai,i.split] = out$logL

#WOULD BE BETTER TO USE EVAL: like in anova.manyany
#      assign(as.character(object1$call[[3]]),yMat[,is.zeroton==FALSE]) 
#      assign(as.character(object2$call[[3]]),yMat[,is.zeroton==FALSE]) 
#        ft.1i=eval(object1$call)
#        ft.2i=eval(object2$call)
        
      # get test likelihood
# If I write a predict function, the following can be done generically and can apply across manyglm, manyany, etc...
      mu.test      = object$glm1.best$family$linkinv(object$X[is.test,] %*% out$all.coefficients)
      weights.test = object$glm1.best$weights[is.test]
      if( object$glm1.best$family$family=="binomial" || object$glm1.best$family$family=="poisson" || pmatch("Negative Binomial",object$glm1.best$family$family,nomatch=0)==1 )
        dev = NA #to save a little computation time, don't get dev when it is not needed by the aic function.
      else
        dev = object$glm1.best$family$dev(object$y[is.test], mu.test, weights.test) #only actually needed by gaussian, gamma, inverse.gaussian
      for(i.lambda in 1:n.lambdai)
        ll.test[i.lambda,i.split]  = -object$glm1.best$family$aic(y.test, 1, mu.test[,i.lambda], weights.test, dev)/2

#      ll.test[abs(ll.test[,i.split])<100*tol[1],i.split]=min(ll.test[,i.split],na.rm=T) #what is this line about???
      df.cv[1:n.lambdai,i.split] = apply(abs(out$all.coefficients)>tol[1],2,sum)
      beta.cv[,1:n.lambdai,i.split] = out$all.coefficients
      phi.cv[1:n.lambdai,i.split] <- out$phis
      counter.cv[1:n.lambdai,i.split] = out$counter

      if(show.progress)
      {
        flush.console()
        print(paste("train/test split number", i.split, "of", n.split, "completed"))
#        print( cbind( ll.test[,i.split], df.cv[,i.split] ) )
      }
  }
  ll.cv = apply(ll.test,1,mean)

## Find which model minimises predictive logL
  ll.min = max(ll.cv, na.rm=TRUE)
  id.min = which(ll.cv==ll.min)[1]

## Find se hence best model by 1 se rule (if it can be computed from replicate test/training splits)
  if(n.split>1)
  {
    ll.se = apply(ll.test,1,sd) * sqrt( n.test / n.levels ) #jackknife-style estimator
    llminusSE = ll.min - ll.se[id.min]
  }
  else
  {
    ll.se = NULL 
    llminusSE = ll.min
  }
  lam.1se = max( object$lambdas[ll.cv>llminusSE], na.rm=TRUE )
  id.1se = which(object$lambdas == lam.1se)[1]

  # choose which criterion to use for "best" model: best ll or 1se rule
  if(best=="1se")
    id.use = id.1se
  else
    id.use = id.min

  # add these funky new objects to the output object.
  object$ll.cv = ll.cv
  object$se = ll.se
  object$lambda = object$lambdas[id.use]
  beta.best   = object$all.coefficients[,id.use]
  object$coefficients = beta.best
  penalty.i = object$lambdas[id.use] * object$penalty

  best = glm1(object$y, object$X, penalty.i, family=object$family, b.init=beta.best, phi.init=object$phis[id.use])
  object$glm1.best = best  


  if(plot==TRUE)
  {
    lls = object$ll.cv[is.na(object$ll.cv)==FALSE]
    ses = object$se[is.na(object$ll.cv)==FALSE]
    dfs = object$df[is.na(object$ll.cv)==FALSE]
    n.lambdai = length(lls)
    y.lab = "log-likelihood (per test observation)"

    plot(dfs,lls,type="l",ylim=range(c(lls-ses,lls+ses)),xlab="Number of terms in model (df)",ylab=y.lab)
    polygon( c(dfs,dfs[n.lambdai:1]) , c(lls+ses,lls[n.lambdai:1]-ses[n.lambdai:1]), col="lightgrey", border="lightgrey" )
    lines(dfs,lls,type="l")
    points( dfs[id.min], lls[id.min],         col="red",  pch=19 )
    points( dfs[id.1se],  lls[id.1se], col="blue", pch=19 )
    legend("bottomleft",c("maximum predictive likelihood","best within 1 se of max"),col=c("red","blue"),pch=19)
  }
  return( object )
    
}
