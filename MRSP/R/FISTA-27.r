################################################################################
#####    FISTA: (F)ast (I)terative (S)hrinkage (T)hresholding (A)lgorithm  #####
#####    An extensible R implementation of the FISTA algorithm             #####
#####    of Beck & Teboulle (2009).                                        #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 24.01.2014, 15:43                                  #####
################################################################################


##### documentation for the main function "fista":
##### --------------------------------------------------------------------------
##### Arguments:
   ##
   ## coef.init:  initial coef object with which optimization begins, its size
   ##             and class define the exact method used. internally, coef is
   ##             always handled as a list.
   ## dat:        data list which contains all the data needed for the calls to
   ##             loglik, gradient or fistaProximal. for example, dat could be a
   ##             list with elements "y" and "x".
   ## offset:     offset vector. 
   ## eta.init:   initial linear predictor. its specific form depends on the
   ##             model used. for lms, glms etc., eta.init will be a vector.
   ##             for multinomial logit models with K categories, as implemented
   ##             for example in package MRSP, eta.init is an nobs x K-1 matrix.
   ## mu.init:    initial expected values. same size as eta.init.
   ## weights:    weight object, usually a vector or list.
   ## model:      a model object containing various functions required by fista.
      ## -------------------------------------------------------------------- ##
      ## necessary slots of "model":
      ##
      ## invlink:  a function with argument "eta" implementing the inverse link
      ##           function.
      ## link:     a function with argument "mu" implementing the link function.
      ## loglik:   a function with arguments "y", "mu" and "weights" that
      ##           returns the loglikelihood.
      ## gradient: a function with arguments "dat", "mu", "coef", "weights",
      ##           "penweights", "grpindex" and "Proximal.args" that gives the
      ##           gradient with  respect to all covariate elements in dat.
      ##           the output must be a list!
      ## sancheck: a function with arguments "mu", "eta" and "weights" that 
      ##           warns if fista converges towards a degenarated solution. 
      ## further slots can be specified for the calls to fistaProximal!
      ## -------------------------------------------------------------------- ##
   ## Proximal.args: a list with objects to pass on to fistaProximal.   
   ## Proximal.control: an object with control parameters for the calls to
   ##             fistaProximal.
   ## grpindex:   a vector or list specifying which covariates in dat form a
   ##             parameter group that will be penalized jointly by group
   ##             lasso type penalties. all predictors with the same number in
   ##             grpindex form one parameter group. grpindex is not used by
   ##             fista itself, but passed on to fistaProximal.
   ## penindex:   a vector or list specifying how each particular predictor
   ##             shall be treated. the entries must match the method definition
   ##             of fistaProximal for the class of the considered coef object.
   ##             its size and structure must also match that of grpindex in the
   ##             sense entries which have the same number in grpindex must also
   ##             have the same value in penindex. please note that no checks
   ##             of any kind are performed by fista itself for speed reasons.
   ## penweights.init: a list with weights for the penalty terms in "penalty".
   ##             its size and form depends on the needs of penalty and thus on
   ##             the class of coef.
   ## tuning:     list of tuning parameters used in the calls to penalty. 
   ##             ATTENTION: tuning must only contain true "lambda"-values. 
   ##             for example, for the elastic net penalty, it could be a list
   ##             with one element called lambda. the distributional factor 
   ##             alpha of the elastic net penalty should then be supplied as 
   ##             an entry in Proximal.args!
   ## max.iter:   maximum amount of iterations.
   ## rel.tol:    relative tolerance for determining convergence of fista.
   ## inverse.stepsize.init: initial inverse.stepsize. must be positive and
   ##             should be chosen small.
   ## adaptive.stepsize: if false, the stepsizes become smaller and can't
   ##             return to bigger values over iterations.
   ## ----------------------------------------------------------------------- ##
   ## generic functions which are called by fista and whose methods must be
   ## defined for each particular function/algorithm/problem that uses fista:
     ## fistaProximal: computes the proximal operator associated with the task
     ##               at hand. input: coef, tuning, penweights, penindex,
     ##               grpindex, Proximal.args, ...
     ##               output: the updated coef object.
     ## updateEta:   a function that updates the linear predictor eta. input:
     ##               dat, coef, offset, weights, ...
     ##               output: the updated eta object.
     ## penalty:      a function computing the value of the overall penalty term
     ##               of the considered model. input: coef, tuning, penweights,
     ##               penindex, grpindex, Proximal.args, ...
     ##               output: a numerical value
     ## updatePenweights: a function for updating the weights used in the
     ##               penalty terms. this is used for, e.g., adaptive lasso-type
     ##               penalties. input: penweights, coef, weights, penindex,
     ##               grpindex, Proximal.args, ...
     ##               output: an updated penweights object.
     ## --------------------------------------------------------------------- ##
##### --------------------------------------------------------------------------
   
fista <- function(coef.init, dat, offset, eta.init, mu.init, weights,
                  model, Proximal.args, Proximal.control = NULL,
                  grpindex, penindex, penweights.init, tuning,
                  inverse.stepsize.init = 1, control = fista.control(),
                  do.sancheck = TRUE, ...)
{
 ## since functions like lapply, Map and Reduce dont preserve the class of an
 ## object, we must store the appropriate fista generics locally.
 clcoef <- class(coef.init)
 fistaProximal <- getMethod("fistaProximal", signature(coef = clcoef))
 updateEta <- getMethod("updateEta", signature(coef = clcoef))
 penalty <- getMethod("penalty", signature(coef = clcoef))
 updatePenweights <- getMethod("updatePenweights", signature(coef = clcoef))
 lossapprox <- getMethod("lossapprox", signature(coefdiff = clcoef))
 proximterm <- getMethod("proximterm", signature(coefdiff = clcoef))
 useSPG <- getMethod("useSPG", signature(coef = clcoef))
 SPG <- getMethod("SPG", signature(coef = clcoef))
 SPGsmoothpen <- getMethod("SPGsmoothpen", signature(coef = clcoef))
 
 ## will SPG be used?
 doSPG <- useSPG(coef = coef.init, penweights = penweights.init, Proximal.args = Proximal.args)
 
 ## initializing
 coef.old1 <- coef.init
 coef <- coef.init
 ## internally, coef must be a list, therefore:
 if(!is.list(coef.init)){
  coef.old1 <- list(coef.old1)
  coef <- list(coef)
 }
 eta <- eta.init
 mu <- mu.init
 penweights <- penweights.init

 max.iter <- control@max.iter
 rel.tol <- control@rel.tol
 m <- control@m
 #if(m > 0){monotone <- F}else{monotone <- T}
 fista.alpha <- numeric(max.iter + 2)
 fista.alpha[1] <- 0
 fista.alpha[2] <- 1
 inverse.stepsize <- inverse.stepsize.init 
 
 update.mu <- model@invlink
 loglik <- model@loglik
 gradient <- model@gradient
 sancheck <- model@sancheck
 environment(sancheck) <- sys.frame(sys.nframe())
 #do.sancheck <- TRUE
 Proximal.args$do.sancheck <- do.sancheck

 ## initialize stuff for the computation of the BB-stepsizes
 if(doSPG){
  SPGlist <- SPG(coef=coef, tuning=tuning, penweights=penweights, grpindex=grpindex,
                 Proximal.args=Proximal.args)
  dualopt <- SPGlist$dualopt
 }
 
 if(control@BB.stepsize){
  BB.min <- control@BB.min
  BB.max <- control@BB.max
  BB.every <- control@BB.every

  eta.init <- updateEta(dat=dat, coef=coef, offset=offset, weights=weights)
  mu.init <- update.mu(eta.init)
  grad <- gradient(dat=dat, mu=mu.init, coef=coef, weights=weights,
                   penweights=penweights, grpindex=grpindex, Proximal.args=Proximal.args)
  if(doSPG){
   #SPGlist$grad <- lapply(SPGlist$grad, function(u){u/100/tuning[[7]]^2})                                                   ###### !
   totalgrad <- Map('-', grad, SPGlist$grad)                                       
  }else{
   totalgrad <- grad
  }
 }

 loglik.old <- loglik(y=dat$y, mu=mu, weights=weights) 
 pen.old <- penalty(coef=coef, tuning=tuning, penweights=penweights,
                    penindex=penindex, grpindex=grpindex, Proximal.args=Proximal.args)
 SPGpen.old <- ifelse(doSPG, SPGsmoothpen(coef=coef, dualopt=dualopt, penweights=penweights,
                      grpindex=grpindex, tuning=tuning, Proximal.args=Proximal.args), 0)
 fn.val.old <- -loglik.old + pen.old + SPGpen.old
 #fn.val.old <- 1e61

 ## a list that keeps the old function values:
 fn.list <- numeric(max.iter + 1)
 fn.list[1] <- fn.val.old

 d.fn <- 1
 fn.val <- 1e30
 fn.val.old2 <- 1e30
 worsen.count <- 0
 iter.count <- 0
 descent.iter.count <- 0
 
 best.iter <- 0
 best.coef <- coef
 best.pen <- pen.old
 best.l <- loglik.old
 best.l.approx <- loglik.old
 best.fn.val <- fn.val.old
 best.fn.val.approx <- fn.val.old
 best.eta <- eta
 best.mu <- mu
 best.inverse.stepsize <- inverse.stepsize
 best.penweights <- penweights
 best.d.fn <- d.fn
 best.proxim <- 0

 if(doSPG) dualopt.old <- dualopt

 ##### the main loop #####
 while(d.fn > control@rel.tol || fn.val > fn.val.old){# || iter.count < 5){
  iter.count <- iter.count + 1
  if(iter.count > control@max.iter){
   if(control@trace){warning(paste("Maximum number of iterations reached"))}
   break
  } 
  coef.old2 <- coef.old1    ## this is used for computations
  coef.old1 <- coef         ## this is used for (potential) convergence checks

  if(control@BB.stepsize) grad.old <- grad
  
  ## perform a sanity check to detect if the model is deteriorating towards a
  ## degenerate solution.
  if(Proximal.args$do.sancheck){
     sancheck(coef = coef, coef.old2 = coef.old2, mu = mu, eta = eta,
              weights = weights, Proximal.args = Proximal.args)
  }

  fista.beta <- (fista.alpha[iter.count] - 1)/fista.alpha[iter.count + 1]
  
  ## compute the search point and its gradient:
  ## only use the fista extrapolation formula if the line it extrapolates is a
  ## descent direction
  #if(fn.val.old < fn.val.old2){
   #descent.iter.count <- descent.iter.count + 1
   #fista.beta <- (fista.alpha[descent.iter.count] - 1)/fista.alpha[descent.iter.count + 1]
   #coefs <- Map('+', coef, lapply(Map('-', coef, coef.old2), function(u){fista.beta*u}))
   #etas <- updateEta(dat=dat, coef=coefs, offset=offset, weights=weights)
   #mus <- update.mu(etas)
   #grad <- gradient(dat=dat, mu=mus, weights=weights, Proximal.args=Proximal.args)
   #fista.alpha[descent.iter.count+2] <- (1 + sqrt(1 + 4*(fista.alpha[descent.iter.count+1])^2)) / 2
  #}else{
   coefs <- coef
   etas <- updateEta(dat=dat, coef=coefs, offset=offset, weights=weights)
   mus <- update.mu(etas)
   if(!control@BB.stepsize){
    grad <- gradient(dat=dat, mu=mus, coef=coefs, weights=weights,
                     penweights=penweights, grpindex=grpindex, Proximal.args=Proximal.args)
    if(doSPG){
     SPGlist <- SPG(coef=coefs, tuning=tuning, penweights=penweights, grpindex=grpindex,
                    Proximal.args=Proximal.args)
     dualopt <- SPGlist$dualopt
     #SPGlist$grad <- lapply(SPGlist$grad, function(u){u/100/tuning[[7]]^2})                                                   ###### !
     totalgrad <- Map('-', grad, SPGlist$grad)                                                                          
    }else{
     totalgrad <- grad
    }
    #inverse.stepsize <- 1
   }
  #}
  


  ## stuff for the stepsize search
  #SPGpen.old <- ifelse(doSPG, SPGsmoothpen(coef=coefs, dualopt=dualopt, penweights=penweights,
  #                 grpindex=grpindex, tuning=tuning, Proximal.args=Proximal.args), 0)
  #lS <- loglik(y=dat$y, mu=mus, weights=weights) 
  fn.val <- 1
  fn.val.approx <- 0
  scalefac <- 1

  #totaliter <- get("totaliter", envir = .GlobalEnv)
  ## stepsize search and main computations
  ## a note for the third line: the loss function is the negative loglik, so 
  ## subtracting the gradient of the loss means adding the gradient of loglik.
  while(fn.val > fn.val.approx & inverse.stepsize < 1e30){
   inverse.stepsize <- scalefac * inverse.stepsize
   coeft <- Map('+', coefs, lapply(totalgrad, function(u){u/inverse.stepsize}))
   tuning.scaled <- lapply(tuning, function(u){u/inverse.stepsize})
   coefprox <- fistaProximal(coef=coeft, tuning=tuning.scaled, penweights=penweights,
                            penindex=penindex, grpindex=grpindex,
                            Proximal.args=Proximal.args, Proximal.control=Proximal.control)
   etaprox <- updateEta(dat=dat, coef=coefprox, offset=offset, weights=weights)
   muprox <- update.mu(etaprox)
   lprox <- loglik(y=dat$y, mu=muprox, weights=weights)
   coefdiff <- Map('-', coefprox, coefs)
   #loss.approx <-  lossapprox(lS, totalgrad, coefdiff, Proximal.args)
   penprox <- penalty(coef=coefprox, tuning=tuning, penweights=penweights,
                      penindex=penindex, grpindex=grpindex, Proximal.args=Proximal.args)
   proxim <- proximterm(coefdiff, Proximal.args)
   if(proxim <= 1e-20) break
   if(doSPG){
    SPGlist <- SPG(coef=coefprox, tuning=tuning, penweights=penweights, grpindex=grpindex,                                         
                   Proximal.args=Proximal.args, doGrad=F)
    dualopt <- SPGlist$dualopt
   }
   SPGpen <- ifelse(doSPG, SPGsmoothpen(coef=coefprox, dualopt=dualopt, penweights=penweights,
                    grpindex=grpindex, tuning=tuning, Proximal.args=Proximal.args), 0) 
   fn.val <- -lprox + SPGpen + penprox                                               #####                                                        
   #fn.list[iter.count + 1] <- -loss.approx
   #fn.val.approx <- max(fn.list[seq(from=max(1, iter.count-m+1), to=iter.count)]) + inverse.stepsize/2e5 * proxim
   #fn.val.approx <- -loss.approx + (inverse.stepsize/2 * proxim) + SPGpen.old                                                      ####
   fn.val.approx <- max(fn.list[seq(from=max(1, iter.count-m+1), to=iter.count)]) - inverse.stepsize/2e5 * proxim
   scalefac <- 2
  } ## end while(fn.val > ...)

  if(inverse.stepsize > 1e31){warning(paste("FISTA couldnt find a feasible stepsize during iteration",
                             iter.count))
                             break}                         
  coef <- coefprox
  eta <- etaprox
  mu <- muprox
  #penprox <- penalty(coef=coefprox, tuning=tuning, penweights=penweights,                                                                           
  #                   penindex=penindex, grpindex=grpindex, Proximal.args=Proximal.args)                                                             

  #SPGpen <- ifelse(doSPG, SPGsmoothpen(coef=coef, dualopt=dualopt, penweights=penweights, grpindex=grpindex, tuning=tuning, Proximal.args=Proximal.args), 0)     
  #fn.val <- -lprox + penprox + SPGpen                                                                                                                 
  #fn.val.approx <- -loss.approx + (inverse.stepsize/2 * proxim) + penprox  + SPGpen.old                                                                                                         ####
  if(doSPG) SPGpen.old <- SPGpen                      #dualopt.old <- dualopt                                                                                  ####

   #print("coefprox")
   #print(coefprox)
   #print("gradient")
   #print(grad)
   #print("totalgrad")
   #print(lapply(totalgrad, function(u){u/inverse.stepsize}))
   #print("iter.count")
   #print(iter.count)
   #print("inverse.stepsize")
   #print(inverse.stepsize)
   #print("SPGpen")
   #print(SPGpen)
   #print("fn.val.approx")
   #print(fn.val.approx)
   #print("-loss.approx")
   #print(-loss.approx)
   #print("fn.val.old")
   #print(fn.val.old)
   #print("fn.val")
   #print(fn.val)
   #print("-lprox")
   #print(-lprox)
   #print("coef")
   #print(coef)
   #print("SPGpen")
   #print(SPGpen)

  fista.alpha[iter.count+2] <- (1 + sqrt(1 + 4*(fista.alpha[iter.count+1])^2)) / 2
  penweights <- updatePenweights(penweights=penweights, coef=coef,
                                  weights=weights, grpindex=grpindex,
                                  Proximal.args=Proximal.args)             

  d.fn <- abs(fn.val - fn.val.old)/(control@rel.tol/10 + abs(fn.val))

  if(fn.val < fn.val.old){
   best.iter <- iter.count
   best.coef <- coef
   best.eta <- eta
   best.mu <- mu
   best.fn.val <- fn.val
   #best.fn.val.approx <- fn.val.approx
   best.penweights <- penweights
   best.pen <- penprox
   best.proxim <- proxim
   best.l <- lprox
   #best.l.approx <- loss.approx
   best.d.fn <- d.fn
   best.inverse.stepsize <- inverse.stepsize
   worsen.count <- max(0, worsen.count - 0.5)
  }else{worsen.count <- worsen.count + 1}
  if(worsen.count >= control@worsen.max){
   if(control@trace){warning(paste("fista terminated after", control@worsen.max,
            "consecutive iterations in which the objective function worsened"))}
   #break
   m <- 0
   #control@BB.stepsize <- F
   slot(control, "BB.stepsize") <- F
  }

  ## computation of the BB-stepsize.
  ## a note for the second line: the BB formula uses the grad of the loss, but
  ## our 'grad' function computes the gradient of the loglikelihood, with
  ## the negative loglik being the loss function. therefore the sign switch.
  if(control@BB.stepsize){
   grad <- gradient(dat=dat, mu=mu, coef=coef, weights=weights,
                    penweights=penweights, grpindex=grpindex, Proximal.args=Proximal.args)
   grad.diff <- do.call("c", Map('-', grad.old, grad))
   if(doSPG){
    SPGlist <- SPG(coef=coef, tuning=tuning, penweights=penweights, grpindex=grpindex,
                  Proximal.args=Proximal.args)
    dualopt <- SPGlist$dualopt
    #SPGlist$grad <- lapply(SPGlist$grad, function(u){u/100/tuning[[7]]^2})                                                   ###### !   
    totalgrad <- Map('-', grad, SPGlist$grad)                                                                   
   }else{
    totalgrad <- grad
   }
   #grad.diff <- do.call("c", Map('-', grad.old, totalgrad))   
   coef.diff <- do.call("c", Map('-', coef, coef.old1))
   BB.num <- crossprod(coef.diff, grad.diff)
   if(BB.num <= 0){
    inverse.stepsize <- BB.min
   }else{
    BB.denum <- crossprod(coef.diff)
    inverse.stepsize <- max(BB.min, min(BB.max, BB.num/BB.denum))
   }
   #print(inverse.stepsize)
  }

  fn.val.old2 <- fn.val.old
  fn.val.old <- fn.val
  fn.list[iter.count + 1] <- fn.val#.approx # -loss.approx + ifelse(doSPG, SPGsmoothpen(coef=coef, dualopt=dualopt, penweights=penweights, grpindex=grpindex, tuning=tuning, Proximal.args=Proximal.args), 0)       ######
  #fn.list[iter.count + 1] <- fn.val
 } ## end while(d.fn > ...)
  
 ## coef was transformed to be a list. before fista returns its output, coef thus
 ## has to be transformed back into the form that coef.init had. if coef.init 
 ## was a vector or matrix, the artificial listversion of coef has length 1.
 ## Note: actually, this does more harm than good. defunct for now.
 #if(length(coef) == 1){best.coef <- best.coef[[1]]}
 
 class(best.coef) <- clcoef

 out <- list(coef = best.coef,
             inverse.stepsize = best.inverse.stepsize,
             penweights = best.penweights,
             eta = best.eta,
             mu = best.mu,
             loglikval = best.l,
             #loglikapprox = best.l.approx,
             penval = best.pen,
             proximval = best.proxim,
             offset = offset,
             fn.val = best.fn.val,
             #fn.val.approx = best.fn.val.approx,
             iter.count = iter.count,
             best.iter = best.iter,
             d.fn = best.d.fn,
             ridgestabil = Proximal.args$ridgestabil,
             ridgestabilrf = Proximal.args$ridgestabilrf)
 return(out)
}

fista <- cmpfun(fista)