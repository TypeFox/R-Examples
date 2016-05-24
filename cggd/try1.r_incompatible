  
  
  
cggd.lm <- function(x, y, beta0 = rep(0,2), kmax = 300, m = 20, eps1 = 1e-06, 
                      TRR=F, t0 = 1, 
                      alpha1 = 0, alpha2 = 0, w = 1, 
                      tau = 1, tautil = -1, eps = -1)
 {
    n <- length(y) # y is n by 1
    p <- ceiling(length(x)/n) # x is n by p
    if (length(beta0)!=p) beta0 <- rep(0,p) # O: make sure initial betas have correct dimensionality
    
    # kmax <- max(1,min(kmax, 10*n, 10*p, 1500))  O: try to limit kmax and m size (computational considerations space/time)
    # m <- max(2,min(m, 200, floor(5000/kmax)))
    
    if ((tau>1)||(tau<=0)) tau <- 1  # O: verify that tau is correct
      if (p==1) { # O: if one variable just directly do least squares and plot m straight line points
      olse <- sum(x*y)/sum(x^2)
      beta <- beta0 + ((0:(m-1))/m)*(olse - beta0)
      beta.tk <- c(beta0, olse)
      tk <-  -1
      a.set.tk <- rep(T,2)
      kmax <- 1
      gmax <- 0
      return(beta, beta.tk, tk, a.set.tk, kmax, m, p, gmax)
    }  # now p > 1
    eps <- eps * median( abs(t(x)%*%y)/apply(x^2,2,sum) ) # O: t-transpose, median of least squares solution for individual least squares
    
    # variable declaration, beta is dynamic 
    beta.tk <- matrix(beta0, p, kmax+1)
    tk <- rep(-1, kmax)   
    a.set.tk <- matrix(F, p, kmax+1) 
    
    # compute the active set of model 1 [O: initial model]
    a.set <- rep(F, p)
    g0 <- t(x)%*%(y - x%*%beta0)  # O: gradient at t0
    gmax <- max(abs(g0))
    a.set[ abs(g0) >= tau * gmax ]  <- T # active set for model 1 [O: a.set is the current active set]
    var.set <- a.set # cumulative union of a.set [O: ultimately will be closure of a.set] 
    a.set.tk[,1] <- a.set 
    k <- 1 
    success <- TRUE
    while (k <= kmax && success) { # compute the trajectory for kmax -1 models
      # step k 
      if (TRR==F) 
        # cggd.mgd.one.step <- function(x,y, a.set, beta0, t1=1, m=20, eps1=1e-08,
        #                 alpha1=0, alpha2=0, w=1, tau=1, tautil=-1, eps=-1)
        #   return(beta, a.set.init, a.set, g.init, g, g.abs, gmax, t1)
        #  O:
        #   beta = new beta is p x m+1 columns , final column is the t_k beta
        #   a.set.init = initial a.set for step k (from k-1 or inital)
        #   a.set = new a.set at t_k
        #   g.init = initial gradient ??
        #   g = final gradient  ????
        #   ??? need to check rest
        
        step.k<-cggd.mgd.one.step(x, y, a.set, beta.tk[,k], # [,k] is the kth column of beta.tk
                m=m, eps1=eps1, alpha1=alpha1, alpha2=alpha2, 
                w=w, tau=tau, tautil=tautil, eps=eps)
        success <- step.k$success        
      else {
        # cggd.trr.one.step <- function(x, y, a.set, beta0, 
        #                 t1 = 1, m = 20, eps1 = 1e-08,
        #                 t0 = 1, tau = 1, tautil = -1, eps = -1)
        if (k>1) t0 <- tk[k-1]
        step.k<-cggd.trr.one.step(x, y, a.set, beta.tk[,k], m=m, eps1=eps1, 
                                  t0=t0, tau=tau, tautil=tautil, eps=eps)
      }
      
      # update beta.tk, beta, var.set
      #beta.tk[,(k+1):(kmax+1)] <- matrix(step.k$beta[,m+1], p, kmax+1-k) # O: change all beta from this point on (has to do with vars going in and out)
      #browser()
      beta.tk[,(k+1):(kmax+1)] <- matrix(step.k$beta[,1], p, kmax+1-k) # O: change all beta from this point on (has to do with vars going in and out)
      
      b <- matrix(step.k$beta[a.set,], sum(a.set), m+1) # O: for plotting purposes only care about the a.set (of t_k-1) not entire p
      
      if (k==1) beta <- b[ ,1:m] # O: initial values of plotting varaible BETA
      
      if (k > 1) { # update beta and var.set  O: increasing size of beta with change in var.set size
        old.beta <- beta
        old.var.set <- var.set
        var.set[a.set] <- T
        beta <- matrix(beta0[var.set], sum(var.set), k*m) 
        beta[old.var.set[var.set], 1:(k*m-m)] <- old.beta
        beta[ ,(k*m-m+1):(k*m)] <- matrix(beta.tk[var.set,k],sum(var.set),m)
        beta[a.set[var.set], (k*m-m+1):(k*m)] <- b[,1:m]
        
        for (i in 1:length(var.set)) {
          if(var.set[i] != old.var.set[i]) cat("var ",i," different\n")
        }
      
      }
      
      # update a.set, a.set.tk, tk and k 
      
      cat("k =", k," model size = ",sum(a.set), " gmax = ", step.k$gmax, "\n")
      
      
      a.set <- step.k$a.set
      a.set.tk[,k+1] <- a.set
      tk[k] <- step.k$t1
      k <- k+1
      
      # if (sum(a.set) > n) kmax <- k-1   # O: if have n variables have perfect fit and exit (bioinformatics)
      
      
      if (step.k$t1 == -1) kmax <- k-1   # O: if previous step is end of algorithm then end loop
    } # end while (k <= kmax) 
    # kmax may have changed due to early stopping
    
    beta.tk <- beta.tk[,1:(kmax+1)]  # O: adjust to current value of kmax
    tk <- tk[1:kmax]
    a.set.tk <- a.set.tk[,1:(kmax+1)]
    gmax <- step.k$gmax       # O: max gradient of last iteration
    # return(beta, beta.tk, tk, a.set.tk, kmax, m, p, gmax)
    return(list(beta=beta, beta.tk=beta.tk, tk=tk, a.set.tk=a.set.tk, kmax=kmax, m=m, p=p, gmax=gmax))
}
  
  
  
 ##############################
 
 
 
 
 cggd.mgd.one.step <- function(x, y, a.set, beta_init,    t1 = 1, m = 20, eps1 = 1e-08,
                      alpha1 = 0, alpha2 = 0, w = 1, tau = 1, tautil = -1, eps = -1){
   
   
    
   # O: as computed before
   beta0 <- beta_init
   #beta <- beta0
   
   n <- length(y)
   p <- length(beta0)
   g0 <- t(x)%*%(y - x%*%beta0) 
   gmax <- max(abs(g0))
   
   
   # eps1 <- gmax * eps1  #O: another scaling
   # O: a.set is updated in the loop, add these vriables as default 
   a.set[ abs(g0) >= tau * gmax ]  <- T
   
   # O: cant be easier to exit than enter
   if (tautil >= tau) {
     tautil <- tau
     w <- 1
     alpha1 <-  -1
   }
   
   # O: m1 = number of points to examine for switching time
   
   m1 <- m * 10
   
   # q = number of active variables
   q <- sum(a.set)
   
   # the decomposition of P(t)GP(t)
   txa<-t(x[,a.set])
   txna<-t(x[,!a.set])
   
   svdG <- svd(txa %*% x[,a.set])  # svd of a q x q sized matrix
   tsvdGv <- t(svdG$v)
   Q1 <-  svdG$d  # diagonal, eigenvalues
   Q2 <-  Q1      # do twice for kernel calculations
   zeros <-  Q1 == 0 
   Q1[!zeros] <- Q1[!zeros]^alpha1
   if (alpha1==0) Q1[zeros]  <- 1
   Q2[!zeros] <- Q2[!zeros]^alpha2
   if (alpha2==0) Q2[zeros]  <- 1
   Q <- w*Q1 + (1-w)*Q2  # kernel calculation
   
   
   # note: c(M) returns a vector of a matrix
   
   
   # this is the kernel x g0
   Qg0 <-  (svdG$u)%*%(Q * c(t(svdG$v)%*%g0[a.set]))  
   
   # H = Q*G in eigen space
   H <- Q * (svdG$d)
   
   #FstarQg <- (svdG$u)%*%(c(tsvdGv%*%Qg0))
 get.switch.time <- function(terminate_value=1e9)
 {
 # This function finds the location of the next switch time
 # It works in two stages. First it uses bump following to find consecutive local minimas
 # of the switch time function, until it finds a negative minima.
 # Then it uses root finding to narrow in on the switch time.
 
 # looking for negative values of the switch.function
 
 boundary<-0
 find_minimum<-TRUE
 more<-TRUE
 
 
 while(more)
 {
 
   if(find_minimum) {
     res<-optim(boundary,cggd.mgd.switch.func.min,method="L-BFGS-B",lower=boundary,control=list(factor=1e2))
   }
   else { # maximize over the next hump
     res<-optim(boundary,cggd.mgd.switch.func.min,method="L-BFGS-B",lower=boundary,control=list(fnscale=-1,factor=1e2))       
   }
 
   if (!res$convergence) { return(list(switch.time=0,success=FALSE)) }
 	
 
   cat("res par =", res$par," res value = ",res$value," counts =", res$counts, " convergence =", res$convergence, " message = ",res$message,"\n")
   if(res$par >= terminate_value || res$value <0) { more<-FALSE }
 
   find_minimum <- !find_minimum
   last_boundary<-boundary
   boundary <- res$par
 }
 
 switch.time<-res$par
 #browser()
 # Find the exact switch time
 # can further optimize by figuering out it is in the m.in or m.out and only doing those
 # can also check which variables in particular are negative and only do it for them
 if( res$value <0  && res$par !=last_boundary) {
     r<-uniroot(cggd.mgd.switch.func.min,c(last_boundary,boundary),tol = .Machine$double.eps^0.5)
 
     switch.time <- r$root + r$estim.prec 
 }
 
 return(list(switch.time=switch.time,success=TRUE))
 }
 

   cggd.mgd.switch.func <- function(t,do_in=T, do_out=T)
     {
      # browser()
       time <- t
       beta <- matrix(beta0,p,1)

       if (q==1) {
         Fstar <-  - time*H
         Fstar <- cggd.Fstar(Fstar)
         beta[a.set,] <- beta[a.set,] + time*Fstar*Qg0
       } # end if (q==1)
       if (q>1) {
       # H matrix for every eigen value x time points calculate 
         Fstar <-  - H*time
         Fstar <- cggd.Fstar(Fstar)

         ### on the active set   beta(t) = beta0 + t F*(-t H) Q g0, H = QG (in eigen space)
         beta[a.set,] <- beta[a.set,] + (svdG$u)%*%(Fstar*c(tsvdGv%*%Qg0))*time # eq 6.65
       } # end if (q>1)

      #browser()
       ### g <- t(x)%*%(matrix(y,n,m1+1)-x%*%beta)
       C <- (matrix(y,n,1)-x%*%beta)
       AA<-abs(txa%*%C)  
       g.in.max <- apply(AA,2,max)

       # calculate the entrance time
       if(do_in)
       {
         BB <-abs(txna%*%C)
         g.out.max <- apply(BB,2,max)

         m.in <- tau*g.in.max - g.out.max 
       }
       else m.in <- 1e99

       # calculate the exit time    
       if (do_out) 
       {
         g.in.min <- apply(AA,2,min)
         m.out <- g.in.min - tautil*g.in.max
       }
       else m.out <- 1e99

       return(list(m.in=m.in,m.out=m.out,beta=beta))     
   }
   
   
   cggd.mgd.switch.func.min <- function(t)
    {
      a <- cggd.mgd.switch.func(t);
      m<-min(c(a$m.in,a$m.out))
      
      return(m)
    }
    
    cggd.mgd.switch.func.in <- function(t)
    {
      a <- cggd.mgd.switch.func(t,do_out=false);
      return(a$m.in)
    }
    cggd.mgd.switch.func.out <- function(t)
    {
       a <- cggd.mgd.switch.func(t,do_in=false);
       return(a$m.out)
    }
    
    ##############################
    
   
     
    

   
    
   switch.time <- get.switch.time ()  
   
   info.init <- get.switch.time (0)
   switch.info <- cggd.mgd.switch.func(switch.time)
   # browser()
   
   g.init <- t(x)%*%(y-x%*%info.init$beta) # the starting g (not sure why we need it)
   g <- t(x)%*%(y-x%*%switch.info$beta)  # the g at the switch time
   g.abs <- abs(g)
   gmax <- max(g.abs[a.set])
   a.set.init <- a.set
   a.set[g.abs >= tau*gmax] <- T
   a.set[ (g.abs <= tautil*gmax) ]  <- F  # & (abs(beta[,m+1]) <= eps) 
   return(list(success=switch.time$success,beta=switch.info$beta, a.set.init=a.set.init, a.set=a.set, g.init=g.init, g=g, g.abs=g.abs, gmax=gmax, t1=t1))
    
 }
 

 diabetes.ntgd.2a <- cggd.lm(x.noise,y,m=1,kmax=11,tau=3/4)
 
 