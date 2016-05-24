## function to check analytical gradient fo GCV/UBRE...


check.analytical <- function(object,data, del=1e-6){

    # require(mgcv)
     G<-gam(object$formula,object$family,data,fit=FALSE)
     n.terms <- length(G$smooth)  # number of smooths in the model
     n <- nrow(G$X)
     intercept <- G$intercept ## TRUE or FALSE
     n.pen <- length(G$S)
     Q <- penalty_pident(G)
     q.f <- rep(0,n.terms)
     for (i in 1:n.terms) {
             q.f[i] <- ncol(G$smooth[[i]]$S[[1]]) +1
     }
     G$S <- Q$S
     G$q.f <- q.f
     G$q0 <- G$off[1]-1  ## number of the parameters of the strictly parametric model
     G$p.ident <- Q$p.ident  # vector of 0's & 1's for the model parameters identification: 
     G$n.terms <- n.terms   ## number of the smooth terms in the SCAM
     G$weights <- object$prior.weights
     G$sig2 <- object$sig2
     G$scale.known <- object$scale.known

     ## Create new environments with `start' initially empty
     env <- new.env()
     assign("start",object$coefficients,envir=env)
     assign("dbeta.start",object$dbeta.rho,envir=env)
     assign("sp.last",object$sp,envir=env)


     sp1 <- rep(0,n.pen)
     sp <- object$sp
     dgcv.ubre.check <- rep(0,n.pen)
     for (j in 1:n.pen){
          sp1<- sp; sp1[j]<-sp[j]*exp(del)

           m1 <- scam.fit(G=G, sp=sp1,env=env,devtol=1e-8, steptol=1e-8)

          if (object$scale.known) gcv.ubre1 <- m1$deviance/n - object$sig2 +2*object$gamma*m1$trA*object$sig2/n
          else gcv.ubre1 <- m1$gcv
          dgcv.ubre.check[j] <- (gcv.ubre1-object$gcv)/del # finite difference derivative
       }

      check.grad <- 100*(object$dgcv.ubre - dgcv.ubre.check)/dgcv.ubre.check

   list(dgcv.ubre.fd=dgcv.ubre.check,check.grad=check.grad)
}






