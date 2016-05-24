
#########################################################
# Function to return gcv/ubre ...                      ##
#########################################################

gcv.ubre <- function(rho,G,gamma,env)
{  ## function to get GCV.UBRE value for optim()...
   if (length(rho)!= length(G$off)) stop (paste("length of rho and n.terms has to be the same"))
   sp <- exp(rho)
   b <- scam.fit(G=G, sp=sp, env=env)
   if (G$scale.known) #  value of Mallow's Cp/UBRE/AIC ....
      {  n <- nrow(G$X)
         gcv.ubre <- b$dev/n - G$sig2 +2*gamma*b$trA*G$sig2/n
      }  else   # value of GCV ...
           gcv.ubre <- b$gcv
   return(gcv.ubre)
}

#########################################################
## function to get the gradient of the gcv/ubre.....   ##
#########################################################

gcv.ubre.derivative <- function(rho,G, gamma,env, check.analytical=FALSE, del) 
{  ## function to return derivative of GCV or UBRE for optim...
   gcv.ubre_grad(rho, G, gamma,env,check.analytical, del)$gcv.ubre.rho
}


#############################################################################
## for nlm() function to get the gcv/ubre and gradient of the gcv/ubre.....##
#############################################################################

dgcv.ubre.nlm <- function(rho,G, gamma,env, check.analytical=FALSE, del)
{  ## GCV UBRE objective function for nlm
   gg <- gcv.ubre_grad(rho, G, gamma,env,check.analytical, del)
   attr(gg$gcv.ubre,"gradient") <- gg$gcv.ubre.rho
   gg$gcv.ubre
}



#######################################################
#### estimate.scam()....                             ##
#######################################################


estimate.scam <- function(G,optimizer,optim.method,rho, gamma,env,
                     check.analytical, del, devtol, steptol)
{  ## function to select smoothing parameter...
   if (!(optimizer %in% c("bfgs", "nlm", "optim","nlm.fd")) )
          stop("unknown outer optimization method")
   if (optimizer == "bfgs") ## minimize GCV/UBRE by BFGS...
      {  b <- bfgs_gcv.ubre(gcv.ubre_grad,rho=rho, G=G,gamma=gamma,env=env,
                         check.analytical=check.analytical, del=del)
         sp <- exp(b$rho)
         object <- b$object
         object$gcv.ubre <- b$gcv.ubre
         object$dgcv.ubre <- b$dgcv.ubre
         object$termcode <- b$termcode
         object$check.grad <- b$check.grad
         object$dgcv.ubre.check <- b$dgcv.ubre.check
         object$conv.bfgs <- b$conv.bfgs
         object$iterations <- b$iterations
      }
   else if (optimizer=="optim")  ## gr=gcv.ubre.derivative
           {  if (!(optim.method[1] %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")))
                 {   warning("unknown optim() method, `L-BFGS-B' were used")
                     optim.method[1] <- "L-BFGS-B"
                 }
              if (is.na(optim.method[2])) 
                 {   warning("the second parameter of optim.method argument is not supplied, 
                      finite-difference approximation of the gradient were used")
                     grr <- NULL
                 }  else if (!(optim.method[2] %in% c("fd","grad")))
                      {   warning("only `fd' and `grad' options are possible, finite-difference 
                              approximation of the gradient were used")
                          grr <- NULL
                      }  else if (optim.method[2] == "grad")
                                 grr <- gcv.ubre.derivative
                              else 
                                 grr <- NULL
              b <- optim(par=rho,fn=gcv.ubre, gr=grr, method=optim.method[1],G=G, gamma=gamma,env=env)
              sp <- exp(b$par)
              gcv.ubre <- b$value
              dgcv.ubre <- NULL
              iterations <- b$counts
              termcode <- b$convergence
              if (termcode == 0)
                     conv <- "Successful completion"
              else if (termcode == 1)  
                      conv <- "The iteration limit `maxit' had been reached"
              else if (termcode == 10)  
                      conv <- "Degeneracy of the Nelder-Mead simplex"
              else if (termcode == 51)  
                      conv <- "A warning from the `L-BFGS-B' method; see help for `optim' for further details"
              else if (termcode == 52)  
                      conv <- "An error from the `L-BFGS-B' method; see help for `optim' for further details"
           }
   else if (optimizer=="nlm.fd") ## nlm() with finite difference derivatives...
           {  b <- nlm(f=gcv.ubre, p=rho,iterlim=100, G=G, gamma=gamma,env=env)
           }
   else if (optimizer=="nlm")  ## nlm() with analytical derivatives...
           { b <- nlm(f=dgcv.ubre.nlm, p=rho,iterlim=100,G=G,gamma=1,env=env,
                     check.analytical=check.analytical, del=del)
           }
   if (optimizer== "nlm.fd" || optimizer== "nlm") 
      {   sp <- exp(b$estimate)
          gcv.ubre <- b$minimum
          dgcv.ubre <- b$gradient
          iterations <- b$iterations
          termcode <- b$code
          if (termcode == 1)
                 conv <- "Relative gradient is close to zero, current iterate is probably solution"
          else if (termcode == 2)  
                 conv <- "Successive iterates within tolerance, current iterate is probably solution"
          else if (termcode == 3)  
                 conv <- "Last global step failed to locate a point lower than `estimate'. Either 
                       `estimate' is an approximate local minimum of the function or 
                        `steptol' is too small"
          else if (termcode == 4)  
                 conv <- "Iteration limit exceeded"
          else if (termcode == 5)  
                 conv <- "Maximum step size `stepmax' exceeded five consecutive 
                         times. Either the function is unbounded below, becomes asymptotic 
                         to a finite value from above in some direction or stepmax is too small"   
      }
   ## fit the model using the optimal sp from "optim" or "nlm"...
   if (optimizer== "nlm.fd" || optimizer== "nlm" || optimizer== "optim")
      {  object <- scam.fit(G=G, sp=sp,env=env, devtol=devtol, steptol=steptol)
         object$gcv.ubre <- gcv.ubre
         object$dgcv.ubre <- dgcv.ubre 
         object$termcode <- termcode
         object$conv <- conv
         object$iterations <- iterations 
      }
   if (optimizer=="optim")
      {  object$optim.method <- rep(NA,2)
         object$optim.method[1] <- optim.method[1]
         if (!is.null(grr))
              object$optim.method[2] <- "grad"
      }
   object$sp <- sp
   object$q.f <- G$q.f
   object$p.ident <- G$p.ident
   object$S <- G$S
   object
}


