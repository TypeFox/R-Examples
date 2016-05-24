getSort <- function(x){
  sortcode = apply(x,2,order)
#  sortcode = x
#  for (i in 1:ncol(x)){
#    sortcode[,i] = order(x[,i])
#  }
  sortcode
}

getSorted <- function(x,y,sortcodes=FALSE){
  if (sortcodes){
    sorted = apply(x,2,function(a) y[a])
  }else {
    sorted = apply(x,2,function(a) y[order(a)])
  }
  sorted
}

#do a parameterwise thing, in the vein of the additive isotone one
#now use multistep format
#UPDATE: add weights (does not work with snapping!)
#UPDATE2: add sign handling code
#UPDATE3: add huber

liso.backfit <- function(x,y,lambda=0,givebeta = FALSE,tol.target = 1e-04,weights = rep(1, length(y)), covweights = rep(1,ncol(x)),feed, trace=FALSE, monotone=TRUE, randomise=FALSE, huber=Inf){
  if (missing(feed)){
    feed  = multistep(coefc=rep(0, (nrow(x)-1) *ncol(x)), x=x, intercept=sum(weights * y)/sum(weights), residual=y-sum(weights *y)/sum(weights), sorter=getSort(x))
  }
#turn x into table of sortcodes
  if (length(lambda) > 1){
    result = NULL
    for (ilam in 1:(length(lambda))){
      if (trace) cat("Starting calc with lambda = " , lambda[ilam], "\n")
        if ((length(dim(covweights)) > 1) && (length(covweights) > ncol(x))){
          result[[ilam]] = liso.backfit(x, y, lambda[ilam], givebeta, tol.target, weights, covweights[ilam,], feed, trace, monotone, randomise)
        } else {
          result[[ilam]] = liso.backfit(x, y, lambda[ilam], givebeta, tol.target, weights, covweights, feed, trace, monotone, randomise)
        }
      feed = result[[ilam]]
    }
    return(result)
  }
  dims <- dim(x)
  nn <- dims[1]
  pp <- dims[2]
  if (length(monotone) == 1) monotone = rep(monotone, pp)
#duplicate the covweights if appropriate. Otherwise take the whole sequence
  if (length(covweights) == pp) covweights = c(covweights, covweights[which(monotone == 0)])
  vn <- dimnames(x)[[2]]
  meany <- sum(weights * y)/sum(weights) #feed$intercept
  y <- y - meany
  x.old <- x
  res <- y - scale(feed * x, scale=F)#feed$residual
  hres <- rep(0, length(y))
  x <- scale(x, center=F, scale = replace(monotone, which(monotone == 0), 1))
  x <- cbind(x, -x[, which(monotone==0)]) 
  sortedx <- apply(x,2,sort)
  fits <- as.array.multistep(feed, explode=T, signs = monotone) #matrix(rep(0,nn*pp),nrow = nn)
  oldrem = Inf
  currvar = 1
  roundoff = 0 # accumulated roundoff error
  pp = pp + sum(monotone == F)
  sortcode <- getSort(x)
  weightsgrid = apply(getSorted(sortcode, weights, T),2, cumsum)
    #we fudge things a bit then
#    sortedx = cbind(sortedx, -apply(sortedx,2,rev))
#    sortcode = cbind(sortcode, apply(sortcode,2,rev))

  dupx = apply(sortedx,2,duplicated)
  dupindex = apply(sortedx,2,function(a) any(duplicated(a)))
  indices = 1:pp
  steps = 0
  if (sum(weights != 1)){
    weightedliso = TRUE
  } else {
    weightedliso = FALSE
  }
  if (trace) print("Starting calc...")
  while (steps < 10000->MAXstep){
    res <- drop(res)
    if (randomise) {
      indices = sample(pp)
    }
    for (i in indices){
      work = res[sortcode[,i]] + fits[,i]
#if covweights is high enough, don't even bother
      if (covweights[i] <= 0){
        fit = rep(0, length(y))
#to avoid 0/0 problems...
        covweights[i] = -1
      } else {
        if (dupindex[i]){
          if (weightedliso){
            fit = pavareg(ave.sorted(work,sortedx[,i], weights = weights[sortcode[,i]]),lambda/covweights[i],long.out=F, weightsum = weightsgrid[,i], weights = weights[sortcode[,i]])
          } else {
            fit = pavareg(ave.sorted(work,sortedx[,i]),lambda/covweights[i],long.out=F  )
          }
        } else {
          if (weightedliso){
            fit = pavareg(work,lambda/covweights[i],long.out=F, weightsum = weightsgrid[,i], weights = weights[sortcode[,i]])
          } else {

            fit = pavareg(work,lambda/covweights[i],long.out=F )
          }
        }
      }
      fits[,i] = fit
      res[sortcode[,i]] = work - fit
    }
    if (huber < Inf){
      res = res + hres
      hres = sign(res) * ( pmax(abs(res) - huber, 0))
      res = res - hres
    }
    if (weightedliso){
      res <- res - (res %*% weights/sum(weights) -> rounderr)#scale(res,scale=FALSE)
      if (trace) roundoff = roundoff + abs(rounderr)  #attr(res, "scaled:center"))
    } else {
      res = res - (mean(res) -> rounderr)  #scale(res, scale=FALSE)
      if (trace) roundoff = roundoff + abs(attr(res, "scaled:center"))
    }
    meany = meany + rounderr


    steps = steps +1
    #print(fits)
     if (weightedliso) {
        score = sum( weights*(res)^2)/2 + sum((fits[nn,]-fits[1,])*lambda/covweights)
      } else {
        score = sum( (res)^2)/2 + sum((fits[nn,]-fits[1,])*lambda/covweights)
      }
      if (huber <Inf) score = score + huber*sum(weights * abs(hres))
      if (trace) cat("Step ", steps, "\t Score: ", score, "\n")
    
    if (sum(abs(oldrem - score)) < tol.target ) {
      beta = as.vector(apply(fits,2,diff))
        break
    }
    score -> oldrem
  }

    pp = length(monotone)
    sortcode <- getSort(x.old)#sortcode[,1:pp]
    sortedx <- apply(x.old, 2, sort)#sortedx[,1:pp]
    resbeta = beta[1:(pp*(nn-1))] 
    nmk = 1
#fold the betas back in
    for (k in 1:pp){
      if (monotone[k] == -1){
        resbeta[ ((k-1)*(nn-1) + 1):(k*(nn-1))] = -rev( resbeta[ ((k-1)*(nn-1) + 1):(k*(nn-1))] )
      } else if (monotone[k] == 0) {
        resbeta[ ((k-1)*(nn-1) + 1):(k*(nn-1))] = resbeta[ ((k-1)*(nn-1) + 1):(k*(nn-1))] - rev( beta[ pp * (nn-1)+ ((nmk-1)*(nn-1) + 1):(nmk*(nn-1))] )
        nmk = nmk +1
      }
    }

    beta = resbeta
  if (trace) cat("Roundoff error: ", roundoff, "\n")
  if (givebeta){
    beta
  }else{
    result = multistep(beta, intercept=meany, sortedx=sortedx, names = vn, steps = MAXstep - steps, residual = res, sorter=sortcode,  tolerance = tol.target, monotone = monotone, lambda= lambda, pinters = -colMeans(apply(rbind(0,matrix(beta, nn-1)),2, cumsum)), Rsq = 1 - sum(res^2) / sum(y^2) , huber.residual=hres)
    class(result) = c("lisofit", "multistep")
    result
  }
}

liso.adaptive = function(x,y,lambda=0,givebeta = FALSE, tol.target = 1e-04,weights = rep(1, length(y)), covweights = rep(1,ncol(x)),feed = multistep(coefc=rep(0, (nrow(x)-1) *ncol(x)), x=x, intercept=sum(weights * y)/sum(weights), residual=y-sum(weights *y)/sum(weights), sorter=getSort(x)), trace=FALSE, monotone=TRUE, randomise=FALSE, retries = 2, signfind = FALSE){
  if (length(lambda) == 1) lambda = rep(lambda,retries)
  for (itry in 1:retries) {
    cat("Run number ", itry, "\n")
    feed = liso.backfit( x,y,lambda[itry],givebeta=FALSE, tol.target,weights, covweights,feed, trace, monotone, randomise)
      if (signfind){
        covweights = as.vector(t(apply(matrix(as.vector(feed), ncol = ncol(x)),2,function(a) c(sum(pmax(a,0)), abs(sum(pmin(a,0)))))))
      } else {
        covweights = abs(feed)
      }
    print(covweights)
  }
  attr(feed, "covweights") = covweights
  feed
}

#UPDATE: Add weights
pavareg <- function(x,lambda,long.out=FALSE,plotting=FALSE,trace=FALSE, weightsum, weights){
  if ( missing(weightsum) || (weights == rep(1, length(x))) ){
    return(pavaregorig(x,lambda,long.out=FALSE,plotting=FALSE,trace=FALSE))
  }
  fit <- unreg <- pava(x,w=weights)
    deltas = diff(unreg)
    if (sum(deltas) < 1e-10){
      if (long.out) {
        return(list(unreg = unreg, lambda = lambda, y = rep(0,length(x))))
      } else {
        return( rep(0, length(x)))
      }
    }
    cross = min(which(unreg >= -1e-10))
    if ((cross == 1) || (( weights[1:(cross-1)] %*% unreg[1:(cross-1)]) >= -lambda)){
      if (long.out) {
        return(list(unreg = unreg, lambda = lambda, y = rep(0,length(x))))
      } else {
        return( rep(0, length(x)))
      }
    }
    changepts = which(deltas > 0)#change comes after changept
    blocks = length(changepts) +1
    nn = length(x)
    lareas = deltas[deltas > 0] * weightsum[changepts]
    uareas = rev(deltas[deltas > 0] * weightsum[nn] -lareas)
    lslice = min(which(cumsum(lareas) > lambda)) #so we know we're just before there
    uslice = min(which(cumsum(uareas) > lambda))
    fit[1:changepts[lslice]] = unreg[changepts[lslice]+1] - (sum(lareas[1:lslice]) - lambda) / weightsum[changepts[lslice]]
    fit[(changepts[blocks - uslice] + 1):length(x)] = unreg[changepts[blocks - uslice]] + (sum(uareas[1:uslice]) - lambda) / ( weightsum[nn] - weightsum[changepts[blocks - uslice]])
  
  if (long.out) {
    list(unreg = unreg, lambda = lambda, y = fit)
  } else {
    fit
  }
}


pavaregorig <- function(x,lambda,long.out=FALSE,plotting=FALSE,trace=FALSE){
  fit <- unreg <- pava(x)
    deltas = diff(unreg)
    if (sum(deltas) < 1e-10){
      if (long.out) {
        return(list(unreg = unreg, lambda = lambda, y = rep(0,length(x))))
      } else {
        return( rep(0, length(x)))
      }
    }
    cross = min(which(unreg >= -1e-10))
    if ((cross == 1) || (sum(unreg[1:(cross-1)]) >= -lambda)){
      if (long.out) {
        return(list(unreg = unreg, lambda = lambda, y = rep(0,length(x))))
      } else {
        return( rep(0, length(x)))
      }
    }
    changepts = which(deltas > 0)#change comes after changept
    blocks = length(changepts) +1
    lareas = deltas[deltas > 0] * changepts
    uareas = rev(deltas[deltas > 0] * length(x) -lareas)
    lslice = min(which(cumsum(lareas) > lambda)) #so we know we're just before there
    uslice = min(which(cumsum(uareas) > lambda))
    fit[1:changepts[lslice]] = unreg[changepts[lslice]+1] - (sum(lareas[1:lslice]) - lambda) / changepts[lslice]
    fit[(changepts[blocks - uslice] + 1):length(x)] = unreg[changepts[blocks - uslice]] + (sum(uareas[1:uslice]) - lambda) / ( length(x) - changepts[blocks - uslice])
  
  if (long.out) {
    list(unreg = unreg, lambda = lambda, y = fit)
  } else {
    fit
  }
}


#get the maximum value of lambda
liso.maxlamb <- function(x=NULL,y=NULL,monotone=TRUE, covweights=rep(1, ncol(x)), weights=rep(1, length(y))){
   sorted = getSorted(x,y)
  if (length(monotone) == 1) monotone = rep(monotone, ncol(sorted))
  if (length(covweights) == ncol(sorted)) covweights = c(covweights, covweights[which(monotone == 0)])
  for (i in 1:length(monotone)){
    if (monotone[i] == -1){
      sorted[,i] = rev(sorted[,i])
    } else if (monotone[i] == 0) {
      sorted = cbind(sorted, rev(sorted[,i]))
    }
  }

  if (sum((weights - rep(1, length(y)))^2) < 1e-16){
    sorted = sorted - mean(y)
    return(max(covweights * apply(sorted,2,function(a) max(cumsum(rev(a))))))
  } else {
    sortweights = getSorted(x, weights)
    sorted = sorted - sum(weights * y)/sum(weights) 
    res = 0
    for (i in 1:ncol(sorted))
        unreg <- pava(sorted[,i],w=sortweights[,i])
        cross = min(which(unreg >= -1e-10))
        if (cross == 1){
          next
        } else {
          res = max(res, -covweights[i] * sortweights[1:(cross-1),i] %*% unreg[1:(cross-1)])
        }
    }
    res
}

#incorporate huber
cv.liso <- function (x, y, K = 10, lambda = NULL, trace = FALSE, plot.it = FALSE, se = TRUE, weights = rep(1, length(y)), weightedcv = FALSE,huber=Inf,covweights= rep(1, ncol(x)), gridsize = 50, gridmin = 0.01, ...) {
    if (is.null(lambda)) lambda = liso.maxlamb(x,y, covweights = covweights)*seq.log(1,0.01,gridsize)
    all.folds <- cv.folds(length(y), K)
    residmat <- matrix(0, length(lambda), K)
    cat("\n Conducting", K, "fold cv with grid", lambda[1], "-", tail(lambda,1), "\n")
    for (i in seq(K)) {
        cat("\nFold: ", i, "\n")
        omit <- all.folds[[i]]
        fit <- liso.backfit(x[-omit, , drop = FALSE], y[-omit], lambda = lambda, trace = trace, weights = weights[-omit],huber=huber, covweights = covweights,...)
        fit <- sapply(fit, function(obj) obj*x[omit,])
        if (length(omit) == 1) 
            fit <- matrix(fit, nrow = 1)
    if (huber < Inf){
          tmp =(y[omit] - fit)^2
          ind = abs(y[omit] - fit) > huber
          tmp[ind] = (huber^2 + 2*huber*(abs(y[omit] - fit) - huber))[ind]

        if (weightedcv){
          tmp =  weights[omit]*tmp
          residmat[, i] <- colSums(tmp)
        } else {
          residmat[, i] <- colMeans(tmp)
        }
    } else {
     if (weightedcv){
        residmat[, i] <- apply(weights[omit]*(y[omit] - fit)^2, 2, sum)
        } else {
        residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
        }
    }

    }
    if (weightedcv){
    cv <- apply(residmat, 1, sum)/sum(weights)
    } else {
    cv <- apply(residmat, 1, mean)
    }
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object <- list(lambda = lambda, cv = cv, cv.error = cv.error, residmat = residmat, optimlam = lambda[which.min(cv)])
    if (plot.it) 
        plotCV(object, se = se)
    print("Optimum lambda = ")
    print(lambda[which.min(cv)])
    invisible(object)
}

plotCV = function(cv.object,se=TRUE, ...){
#  attach(cv.object)
      plot(cv.object$lambda, cv.object$cv, type = "b", ylim = range(cv.object$cv, cv.object$cv + cv.object$cv.error, cv.object$cv - cv.object$cv.error), ...)
    if(se)
      error.bars(cv.object$lambda, cv.object$cv + cv.object$cv.error, cv.object$cv - cv.object$cv.error, 
                 width = 1/length(cv.object$lambda))
#  detach(cv.object)
  
invisible()
}

liso.covweights = function(obj, signfind = FALSE){
      if (signfind){
        covweights = as.vector(t(apply(matrix(as.vector(obj), ncol = nrow(obj$params)),2,function(a) c(sum(pmax(a,0)), abs(sum(pmin(a,0)))))))
      } else {
        covweights = abs(obj)
      }
      covweights
}

#print liso fit (different from multistep)
print.lisofit <- function(x, ...){
  nonzeroes = which(x$params[,2]!=0)
  totalvar = abs(x)
  signs = sign(matrix(liso.covweights(x, T), nrow(x$params)))
  signs = signs[,1] - signs[,2]
  cat("\nLiso fit: n = ", x$n, ", p = ", nrow(x$params), ", Lambda = ", x$lambda ,"\n")
  cat("Fitted ", nrow(x$values), " steps in ", length(nonzeroes), "variables:\n")
  for (k in nonzeroes){
    if (is.null(dimnames(x$params)[[1]])){
    cat("-Variable ", k, ": ")
    } else {
    cat("-Variable ", k, " (", dimnames(x$params)[[1]][k],"): ")
    }
    cat(x$params[k,2], " steps, ")
    cat("totalvar = " , totalvar[k], ", ")
    if (signs[k] == 1){
      cat("increasing ")
    } else if (signs[k] == -1){
      cat("decreasing ")
    } else {
      cat("non-monotonic ")
    }
    if (signs[k] != x$monotone[k]) {
      cat("(estimated)")
    } else {
      cat("(specified)")
    }
    cat("\n")
  }
  cat("RSS = ", sum(x$res^2), ", R^2 = ", x$Rsq, "\n") 
}

# totally automated fit wrapper
liso <- function(x, y, adaptive = TRUE, lambda = NULL, monotone = NULL, control=list(cv = NULL, liso = NULL )){
  signfind=FALSE
  if (is.null(monotone)) {
    signfind=TRUE
    monotone=FALSE
  }
  if (is.null(lambda)){
print("Conducting crossvalidation")
    cvobj = do.call("cv.liso", c(list(x=x,y=y, monotone=monotone), control$cv))
    lambda1 = cvobj$optimlam
  } else {
    lambda1 = lambda[1]
  }
   print("Making initial fit")
    fitobj = do.call("liso.backfit", c(list(x=x,y=y, lambda=lambda1, monotone=monotone), control$liso)) 
  if (!adaptive){
    return(fitobj)
  } else {
    covweights = liso.covweights(fitobj, signfind)
    if (length(lambda) < 2){
    print("Conducting second crossvalidation")
    cvobj2 = do.call("cv.liso",  c(list(x=x,y=y, monotone=monotone, covweights=covweights), control$cv))
    lambda2 = cvobj2$optimlam
    } else {
      lambda2 = lambda[2]
    }

    print("Making second fit")
    fitobj2 = do.call("liso.backfit", c(list(x=x,y=y, lambda=lambda2, monotone=monotone, covweights=covweights), control$liso)) 
  }
  fitobj2
}
