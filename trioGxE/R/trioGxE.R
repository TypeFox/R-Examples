trioGxE <- function(data,pgenos,cgeno,cenv,
                    penmod=c("codominant","dominant","additive","recessive"),
                    k=NULL,knots=NULL,sp=NULL,
                    lsp0=NULL,lsp.grid=NULL,
                    control=list(maxit=100,tol=10e-08,trace=FALSE),
                    testGxE=FALSE,return.data=TRUE,...){
  
  ## define other "second-" and "third-" level functions
  ## calculates W and W^-1: weight and inverse-weight matrices
  weight.inv.weight <- function(mu,n3,penmod){
    ## mu: current mu
    ## mu is a scalar value for mt1 and mt2, but it is a vector of length 2 for mt3
    ## all conditions are met; mu's are in the right range and available mating type
    
    ## initializing
    W3 <-inv.W3 <- initmat3
    diagvals <- inv.diagvals <- init.diagvals
    
    ## calculate weights
    if(penmod == "codominant"){
      ## calculate weights (for pseudo-replications) for mating types (1,0) and (1,2)
      diagvals[ind12] <- mu[ind12]*(1-mu[ind12])
      
      ## construct a weight matrix for mating type (1,1) only
      diagvals[ind31] <- (mu[ind31]+mu[ind32])*(1-(mu[ind31]+mu[ind32])) #**
      diagvals[ind32] <- mu[ind32]*(1-mu[ind32]) #**
      W3[diag.index] <- diagvals[!ind12] 
      
      off.diag <- mu[ind32]*(1-(mu[ind31]+mu[ind32]))
      W3[off.index31] <- W3[off.index32] <- off.diag
    }# if: penmod is codominant
    
    else{
      if(penmod  == "dominant"){
        diagvals[ind1] <- mu[ind1]*(1-mu[ind1])
        diagvals[ind31] <- (1-(mu[ind31]+mu[ind32]))*(mu[ind31]+mu[ind32])
        W3[diag.index] <- diagvals[!ind12]
        off.diag = 0
      }
      else if(penmod  == "recessive"){
        diagvals[ind2] <- mu[ind2]*(1-mu[ind2])
        diagvals[ind32] <- mu[ind32]*(1-mu[ind32]) #**
        W3[diag.index] <- diagvals[!ind12]
        off.diag = 0
      }
      else{#multiplicative
        diagvals[ind12] <- mu[ind12]*(1-mu[ind12]) 
        diagvals[ind31] <- (1-(mu[ind31]+mu[ind32]))*(mu[ind31]+mu[ind32])
        diagvals[ind32] <- mu[ind32]*(1-mu[ind32]) #**
        W3[diag.index] <- diagvals[!ind12] 
        
        off.diag <- mu[ind32]*(1-(mu[ind31]+mu[ind32]))
        W3[off.index31] <- W3[off.index32] <- off.diag
      }
    }# else: penmod is not codominant
    ## takes care of weight that are too small!!
    ## manual matrix inversion
    tol = 100*(.Machine$double.eps)
    ## checking needed! **?? 
    if(penmod %in% c("codominant","additive")){
      ## calculate W^(-1) for mating types (0,1) and (1,2)
      diagvals[abs(diagvals) < tol] <- tol
      diagvals[diagvals > (1/tol)] <- (1/tol)
      inv.diagvals[ind12] <- 1/diagvals[ind12]
      new.off.diag <- (-1)*off.diag
      det.mt3 <- (diagvals[ind32]*diagvals[ind31]-(new.off.diag^2))
      det.mt3[abs(det.mt3) < tol] <- tol ## determinant is too small
      det.mt3[det.mt3 > (1/tol)] <- (1/tol) ## determinant is too small
      inv.det.mt3 <- 1/det.mt3
      inv.W3[diag.index][which.31] <- inv.det.mt3*diagvals[ind32]
      inv.W3[diag.index][(which.31+1)] <- inv.det.mt3*diagvals[ind31]
      inv.W3[off.index31] <- inv.W3[off.index32] <- inv.det.mt3*new.off.diag
    }#if: penmod is codominant or multiplicative
    
    else{#penmod is neither codominant nor additive
      if(penmod == "dominant"){
        diagvals[ind1|ind31][abs(diagvals[ind1|ind31]) < tol] <- tol
        diagvals[ind1|ind31][diagvals[ind1|ind31] > (1/tol)] <- 1/tol
        inv.diagvals[ind1] <- 1/diagvals[ind1]
        inv.W3[diag.index][which.31] <- 1/diagvals[ind31]
      }
      else {#(penmod%in%c("rec","recessive"))
        diagvals[ind2|ind32][abs(diagvals[ind2|ind32]) < tol] <- tol
        diagvals[ind2|ind32][diagvals[ind2|ind32] > (1/tol)] <- 1/tol
        inv.diagvals[ind2] <- 1/diagvals[ind2]
        inv.W3[diag.index][(which.31+1)] <- 1/diagvals[ind32]
      }#recessive
    }  
    #assign("W3",W3,envir=.GlobalEnv) ## used in triogam.fit.outer() to calculate sqrt(W)
    #assign("inv.W3",inv.W3,envir=.GlobalEnv)
    wdiag <- diagvals
    inv.wdiag <- inv.diagvals
    #assign("wdiag",diagvals,envir=.GlobalEnv) ## used to calculate W^(1/2)
    #assign("inv.wdiag",inv.diagvals,envir=.GlobalEnv)
    res<- list(W3=W3,inv.W3=inv.W3,wdiag=diagvals,inv.wdiag=inv.diagvals)
    res
  }
  
  ## 
  W.sqrt.mat = function(mat){
    ## diag.index,off.index31,off.index32){
    ## mat = any block-diagonal matrix: each block is a 2x2 matrix
    ## initialized
    R = mat ## R will be a square-root of the matrix 'mat'
    diag.index1 <- 2*(1:(nrow(mat)/2))-1# indexes for the odd-number diagonal elements of the matrix
    
    ## extracting diagonal and off-diagonal elements
    diag1 <- mat[diag.index[diag.index1,]]
    diag2 <- mat[diag.index[(diag.index1+1),]]
    off.diag1 <- mat[off.index31]
    off.diag2 <- mat[off.index32]
    
    tau = diag1 + diag2 ## n3-by-1 vector
    det = diag1*diag2-off.diag1*off.diag2 ## ## n3-by-1 vector: all elements should be positive for a weight mat
    
    tol = 100*(.Machine$double.eps)
    if(all(det>=0)){
      det[det < tol] <- tol ## determinant is too small
      det[det > (1/tol)] <- (1/tol) ## determinant is too small
      
      s = sqrt(det) #take only positive value
      t2 = tau+2*s
      stopifnot(t2>0)
      t <- sqrt(t2)
      t[t < tol] <- tol ## t is too small
      t[t > (1/tol)] <- (1/tol) ## t is too small
      
      ## provided that t is not too big
      R[diag.index[diag.index1,]] <- (1/t)*(mat[diag.index[diag.index1,]]+s)
      R[diag.index[(diag.index1+1),]] <- (1/t)*(mat[diag.index[(diag.index1+1),]]+s)
      R[off.index31] <- (1/t)*off.diag1
      R[off.index32] <- (1/t)*off.diag2
      
      return(R)
    }
    
    else {
      neg.ind <- det<0  
      if(all(abs(det[neg.ind]) <= tol)) {
        det[neg.ind] <- tol
      }
      
      else 
        stop("Determinat of the weight matrix is not positive.")
    }
  }
  
  
  ## 
  triogam.fit <- function(y,X,H,S,lsp,n,n12,n3,mt,offsets,
                          mustart,start=NULL,etastart=NULL,weights,
                          null.coef,control,k,penmod,is.final=FALSE){#used in step-halving method
    triogam.fit.update <- function(X,S,L,Z,gamma,scale,n,n12,n3,penmod){
      ## something similar to magic() in mgcv
      ##L: smoothing parameter matrix with
      ## diag(c(rep(lsp[1],k1),rep(lsp[2],k2))),  
      ## routine to do the model fitting, rank-determination, score-calculation,
      ## and returning information needed to calculate derivatives for gcv/UBRE minimization
      ##
      ## H: for stabilizing S: [[not used at the moment]]
      ## Step 1: define y.tilde and X.tilde
      ## Define X.tilde, Y.tilde and S.tilde similarly to those defined in p.185
      ## These definitions already incorporate the constraints (for intercept)
      ## ** TO DO: haven't incorporated H matrix **
      #sqrt.W = mat.sqrt(W)
      diagvals <- y.tilde <- init.diagvals 
      X.tilde <- X
      
      diagvals[ind12] <- sqrt(wdiag[ind12]) ##updated W-values
      ## using the trick:
      ##matrix(c(1:4), nrow=4, ncol=5)
      ##      [,1] [,2] [,3] [,4] [,5]
      ##[1,]    1    1    1    1    1
      ##[2,]    2    2    2    2    2
      ##[3,]    3    3    3    3    3
      ##[4,]    4    4    4    4    4
      X.tilde[ind12,] <- matrix(diagvals[ind12],nrow=n12,ncol=ncol(X))*X[ind12,]
      y.tilde[ind12] <- diagvals[ind12]*Z[ind12]
      if(n3){
        sqrt.W3 <- initmat3
        if(penmod =="codominant"){
          ## Use explicit formula for 2x2 matrices (for f1, f2)
          ## REF: http://en.wikipedia.org/wiki/Square_root_of_a_2_by_2_matrix
          ## W = 
          ##     [,1] [,2]
          ##[1,] w_11 w_12
          ##[2,] w_21 w_22
          ## [diag.index][which.31]: w_11
          ## [diag.index][(which.31+1)]: w_22
          ## off.index31: w_12
          ## off.index32: w_21
          sqrt.W3 <- W.sqrt.mat(W3) ## not stable??
        }#if: penmod is codominant
        
        else{# penmod is not codominant
          sqrt.W3 <- sqrt(W3)
        }#else: penmod is not codominant
        X.tilde[!ind12,] <- sqrt.W3%*%X[!ind12,] #0.005 sec
        y.tilde[!ind12] <- sqrt.W3%*%Z[!ind12] #0.001 sec
      }
      
      S.tilde=L*S
      ##
      nS = nrow(S.tilde)
      #option1 <- function(){ ## Simon Wood's C function
      #myBmat = mat.sqrt(S.tilde) # 0.018sec (for 56 runs) CHECK IF LAPACK IS USED FOR THIS FUNCTION
      #er = .C("mroot", as.double(S.tilde),as.integer(nS),as.integer(nS))
      #myBmat = matrix(er[[1]],nS,nS)
      #}
      ## Calculate square root of the penalty matrix without using Simon Wood's C function
      #option2 <- function(){
      S.eig <- eigen(S.tilde)
      tol <- 100*.Machine$double.eps
      eig.vals <- S.eig$values
      eig.vals[abs(eig.vals/eig.vals[1])<tol] <- 0 ## replace 'too small' values with zero
      myBmat <- S.eig$vectors %*% diag(sqrt(eig.vals)) %*% solve(S.eig$vectors)
      #}
      # QR-decomposition of X.tilde
      qrx <- qr((X.tilde), LAPACK=TRUE)
      ## Copy out upper triangular factor R and unpivot the columns of it
      Q = qr.Q(qrx) #Q1
      R1 = qr.R(qrx)
      qrx.pivot = qrx$pivot
      myBmat = myBmat[qrx.pivot, qrx.pivot]
      S.tilde = S.tilde[qrx.pivot, qrx.pivot]
      
      ## approach in Wood (2008)
      RBmat <- rbind(R1,myBmat) # qr.R(qrx) returns k-by-k upper triangular matrix 
      ## R-Routine
      svd.RB <- La.svd(RBmat)
      U = svd.RB$u
      vt = svd.RB$vt
      d.values = (svd.RB$d)
      
      ## Determine the rank-deficiency of the least squares problem
      ## through the singular values svd.RB$d (numerical rank-deficiency)
      too.small.sv = ((d.values/d.values[1]) < sqrt(.Machine$double.eps))# changed     
      
      ## Any singular values that are "too small" compared to the largest singular value
      ## are to be removed along with the corresponding columns of U and rows of t(V)
      if(any(too.small.sv)){ 
        cat("Rank deficiency detected in the model matrix:",
            sum(too.small.sv),"singular values are too small and removed.\n")
        ##print(d.values/d.values[1])
        d.values = d.values[!too.small.sv]
        U = U[,!too.small.sv] 
        vt = vt[!too.small.sv,]
      }
      
      u1=U[1:nrow(R1),] #OK
      
      ## Define submatrix U1 of U such that R = U1%*%D%*%t(V)
      y1 = t(u1)%*%(t(Q)%*%y.tilde) #q-by-1 vector
      
      trA = sum(diag(crossprod(t(u1))))#trA = tr(u1%*%t(u1)) p.184
      delta = n-gamma*trA
      
      ## calculations of coef estimates and ubre score
      b = (t(vt)%*%diag(1/d.values))%*%y1
      
      ## un-shuffling beta's and penalty matrix and sqrt of it
      unpivot = order(qrx.pivot)
      
      ## WHY S.tilde???
      b = b[unpivot] ; myBmat = myBmat[unpivot,unpivot]; S.tilde=S.tilde[unpivot,unpivot]
      rm(unpivot)
      res = list(b=b,
                 qrx=qrx,
                 rS=myBmat,
                 delta=delta,
                 u1=u1,
                 vt=vt,
                 d.values=d.values,
                 X.tilde=X.tilde,
                 S.tilde=S.tilde)
      res
    }# function ends
    ## remove.linear.separation: indicator specifying whether or not to use the results for
    ## under the parameterization which results in linear-separation.
    remove.linear.separation = TRUE   
    
    ## tem.trace: trace UBRE
    ## L = projection matrix of (sp)'s
    scale = 1; gamma = 1;
    L = matrix(0, nrow=nrow(S), ncol=ncol(S))
    
    if(penmod =="codominant"){
      k1 <- k[1]
      k2 <- k[2]
      L[1:k1,1:k1] = exp(lsp[1])
      L[-c(1:k1),-c(1:k1)] = exp(lsp[2])
    }#penmod is codominant
    
    else{#penmod is not codominant
      k <- k[k!=0][1]
      L[1:k,1:k] = exp(lsp)
    }
    
    if(is.null(mustart)){
      mustart <- trio.mustart(y=y,mt=mt) #estimated from data
    }
    
    if(is.null(etastart))
      null.eta <- X %*% null.coef + as.numeric(offsets)
    else null.eta = etastart ##??
    
    ##---- P-IRLS iteration starts here. ----##
    maxit = control$maxit; tol=control$tol;#convergence control same as gam.fit
    old.pdev = (-2)*triogam.loglkhd(y=y,mt=mt,eta=null.eta,w=weights) + (t(null.coef)%*%((L*S)%*%null.coef))
    
    ##cat("null penalized deviance = ",old.pdev,"\n",sep="")
    mu <- mustart
    eta <- link(mu=mu, mt=mt) #initial additive predictor based on initial means
    Zpart <- rep(0,length(y)) ## W^(-1)%*%deriv.fcn
    
    eps <- 10 * .Machine$double.eps
    
    breakFlag <- FALSE
    for(iter in 1:maxit){
      ## calculate W and Z with the current coefficient estimates: get W^(k) and Z^(k)
      #print(system.time(( #0.112 second/run with SparseM
      weight.info = weight.inv.weight(mu=mu,n3=n3,penmod=penmod) #mt??
      W3 = weight.info$W3
      inv.W3 = weight.info$inv.W3
      wdiag=weight.info$wdiag
      inv.wdiag=weight.info$inv.wdiag    
      #)))
      
      deriv.vals <- deriv.fcn(y=y,mu=mu,mt=mt,w=weights)
      ind12 <- mt==1 |mt==2
      if(any(ind12))
        Zpart[ind12] <- inv.wdiag[ind12]*deriv.vals[ind12]
      else
        Zpart <- Zpart
      
      if(any(!ind12)) ## i.e., trios from Gp=3 are available
        Zpart[!ind12] <- inv.W3%*%deriv.vals[!ind12]
      else
        Zpart <- Zpart
      ##debugging
      Z <- (eta-offsets)+Zpart #pseudo-response ****
      
      ## Update:
      fit.res = triogam.fit.update(X=X,S=S,L=L,Z=Z,
                                   scale=scale,gamma=gamma,n=n,n12=n12,n3=n3,penmod=penmod)
      start = fit.res$b #beta^(t)
      rS = fit.res$rS #re-shuffled sqrt of sp*S...
      penalty = crossprod(rS%*%start) #penalty at t-th iteration
      
      eta <- drop(X%*%start)+offsets ## drop makes the matrix a vector
      mu <- linkinv(eta=eta, mt=mt)
      linear.separation <- (any(mu > 1 - eps) || any(mu < eps))
      if (linear.separation) {# *** needs to be worked on more carefully
        warning("fitted probabilities numerically 0 or 1 occurred.")
        if(remove.linear.separation){
          #cat("linear separation: ",linear.separation,"\n",sep="")
          #ubre <- Inf
          conv = FALSE
          break
        }
      }
      log.lkhd=triogam.loglkhd(y=y,mt=mt,eta=eta,w=weights)
      dev <- -2*log.lkhd # verify!!
      pdev = dev + penalty
      
      if(control$trace){
        cat("pdev = ", pdev, "; deviance =",dev, "; penalty =",penalty,
            "; Iteration -", iter, "\n")
        ##---- Check for convergence ----#
        ## taken from gam.fit3 of mgcv package:
        #********
        cat("|pdev - old.pdev| = ",abs(pdev - old.pdev),"\n")
        cat("linear separation: ",linear.separation,"\n",sep="")
      }
      
      #------- step halvig starts here-------#
      ## not implemented 
      #------- step-halving ends here -------#
      if(pdev == Inf){ #if pdev = Inf break
        breakFlag <- TRUE
        conv = FALSE
        warning("P-IRLS algorithm may have not converged: penalized deviance is Inf.\n" )
        break
      }
      
      if(abs(pdev - old.pdev)/(0.1 + abs(pdev))<tol){## ||iter>fixedSteps){#??
        if(max(abs(start-coefold))>tol*max(abs(start+coefold)/2)){
          ##assuming iter=1, it will not satisfy the condition
          old.pdev=pdev
          etaold=eta
          coef=coefold=start
        }
        
        else{
          conv = TRUE
          coef=start
          break
        }
      }#if(pdev - old.pdev)...)ends
      
      else{
        old.pdev=pdev
        etaold=eta
        coef=coefold=start
      }
    }#for(iter...) ends here P-IRLS convergence
    
    if(iter==maxit){ ## before it gets to MAXIT...needs to address linear.separation issue!
      conv = FALSE
      ubre = NA
      #warning("P-IRLS algorithm have not converged in ",iter," iterations.\n" )
      cat("P-IRLS algorithm have not converged in ",iter," iterations.\n" )
    }
    
    else #*** NOT WORKING WELL PROPERLY
    {
      if(linear.separation & remove.linear.separation){
        ubre=NA ## set it to be Inf so that it does not get used in the
      }
      else
        #        ubre = (1/(n-n3))*(dev -2*fit.res$delta*scale + (n-n3)*scale) ## since n=n+n3
        ubre = (1/n)*(dev -2*fit.res$delta*scale + (n)*scale) ## since n+n3 seems to be more appropriate?
    }
    ##cat("ubre:",ubre,"\n")
    ##cat("coef.est: \n\n",sep="");print(coef)
    
    ## added after the meeting on Oct 13, 2011
    if(is.final){
      full.fit <- list(ubre=ubre,
                       conv=conv,
                       log.lkgd=log.lkhd,
                       penalty=penalty,
                       deviance=dev,
                       pen.deviance = pdev,
                       coef=coef,
                       pirls.iter=iter,
                       qrx=fit.res$qrx,
                       u1=fit.res$u1,
                       vt=fit.res$vt,
                       d.values=fit.res$d.values,
                       X.tilde=fit.res$X.tilde,
                       S.tilde=fit.res$S.tilde)
    }
    
    else{
      full.fit <- list(ubre=ubre, conv=conv)
    }
    
    full.fit
  }
  
  ## data format:
  ## G, Gm, Gf: numeric, numbers of copies of the index-allele for the child, 
  ##            his/her mother, and father, respectively. 
  ## E: numeric, a continuous variable 
  ##    G      E       Gm Gf
  ## 1  1 -1.596706835  1  0
  ## 2  1 -1.100675492  1  0
  ## 3  1 -0.807925652  1  0
  ## penmod: assumed inheritance mode (default: codominant), 'dom' for dominant,
  ##         'add' for -log-additive and 'rec', for recessive model.
  ##         when 'penmod' is non-null, 'k' must be an integer.
  ## pgenos: names of columns corresponding to parental genotypes: can be a length-2 vector
  ## (p1, p2) where p1 and p2 are 0, 1 or 2 (number of the copies of the index allele) or
  ## a scalar of Gp, mating types (01,12 or 11)
  ## cgeno: name of the child genotype column
  ## k: c(k1,k2): numbers of knots for the GRR1 and GRR2 functions
  ## knots: position of the knots (list(xk1=xk1, xk2=xk2), pos1 and pos2 are vectors of positions)
  ## ...
  ## sp: a vector of two smoothing parameter values (sp1, sp2) controlling smoothness of 
  ##the GRR curves for G=1 vs. G=0 and G=2 vs. G=1; if provided, smoothing parameters are not estimated 
  ## weights = prior weights: do we need these?
  ## C: constraints other than centering sum!!
  ## 'risk.allele' removed (20100513)
  ## lsp.grid
  ## lsp0 = should be log(sp) estimates obtained from Linnea's algorithm. (WILL BE REMOVED LATER?)
  ## '...': further argument for passing on e.g., to triogam.fit: 
  ##        start=NULL,etastart=NULL,mustart=NULL,weights=1,
  ## FIRST INCORPORATE HER ALGORITHM INTO MINE!!
  ## testGxE: if TRUE, the function will return only 'coef' and frequentist 'Ve' and Bayesian 'Vf'
  object <- list()
  trace = control$trace
  
  ## Do basic error checking and create a data set with the binary response variable y
  ## introduced in Shin (2012), page 26
  triodat <- pre.trioGxE(data,pgenos,cgeno,cenv,testGxE)
  n=nrow(triodat) #n<-n+n(MT3); augmented data
  y=triodat[,"y"]; x=triodat[,"x"]; mt=triodat[,"mt"]
  ind31 <- mt==3.1; ind32 <- mt==3.2
  h.ind1 <- mt==1|ind31; h.ind2 <- mt==2|ind32
  
  penmod <- match.arg(penmod)
  ## ----------- Knot selection ---------------#
  if(penmod == "codominant") {
    if(is.null(knots)){
      if(is.null(k))## 'k' is also null
        k1 <- k2 <- 5
      
      else{## 'k' is not null
        stopifnot((length(k) == 2), all(k>2))
        k1 <- k[1]
        k2 <- k[2]
      }
      xu1 <- unique(x[h.ind1])
      xu2 <- unique(x[h.ind2])
      xk1 <- quantile(xu1,probs=seq(0,1,length = k1))
      xk2 <- quantile(xu2,probs=seq(0,1,length = k2))
    }# if: 'knots' is null
    
    else{#knots is not null
      stopifnot(is.list(knots), length(knots) == 2, length(knots[[1]])>=3, length(knots[[2]])>=3)
      
      if(!is.null(k)) {#k is not null, either
        warning("Both \'k\' and \'knots\' are provided: \'k\' will be ignored.")        
      }
      k1 <- length(knots[[1]])
      k2 <- length(knots[[2]])
    }
    
    knots <- list(xk1=xk1,xk2=xk2)
    k <- c(k1,k2)        
  }# codominant ends
  
  else{## not codominant: 
    ## 'k' is a scalar; 'knots' is a vector
    if(is.null(knots)){
      if(is.null(k))## 'k' is also null
        k <- 5
      else ## 'k' is not null
        stopifnot((length(k) == 1), k>=3)
      
      ## select the positions of knots under dom, add, rec penetrances
      if(penmod == "dominant"){
        k1 <- k
        k2 <- 0
        xk1 <- quantile(unique(x[h.ind1]),probs=c(0:(k-1))/(k-1))
        xk2 <- NULL
      }
      else if(penmod == "additive"){
        k1 <- k2 <- k
        xk1 <- xk2 <- quantile(unique(x),probs=c(0:(k-1))/(k-1))
      }
      else{#(penmod == "recessive") 
        k1 <- 0
        k2 <- k
        xk1 <- NULL
        xk2 <- quantile(unique(x[h.ind2]),probs=c(0:(k-1))/(k-1))
      }
    }# if: 'knots' is not provided
    
    else{#knots is not null: it must be a vector with at least three elements
      stopifnot(!is.vector(knots), length(knots)>=3)
      
      if(!is.null(k)) #k is not null, either
        warning("Both \'k\' and \'knots\' are provided: \'k\' will be ignored.")    
      
      if(penmod == "dominant"){
        xk1 <- knots
        xk2 <- NULL
        k1 <- length(knots)
        k2 <- 0
      }
      else if(penmod == "additive"){
        xk1<- xk2 <- knots
        k1 <- k2 <- length(knots)
      }
      else{#(penmod == "recessive") 
        xk1 <- NULL
        xk2 <- knots
        k1 <- 0
        k2 <- length(knots)
      }
    }
    ## update knots and k
    knots <- list(xk1=xk1,xk2=xk2)
    k <- c(k1,k2)    
  }# else: not codominant
  
  ## Incorporating the constraints for the intercepts: (reparameterize X and S)
  ## before...this had taken care of the dimension of the reduced matrix
  if(penmod== "additive"){
    bs.dim <- sum(k)/2
  }
  else
    bs.dim <- sum(k)
  
  X <- matrix(0,nrow=n,ncol=bs.dim) ## WHY DID I USE 0 instead of 1? 
  S <- matrix(0,nrow=bs.dim,ncol=bs.dim)
  
  ## Some of the following codes used for incorporating constraints are 
  ## similar to those in mgcv's 'smoothCon()' function.
  if(!(penmod == "codominant") ){
    if(penmod=="dominant"){
      ## Setting up design and penalty matrices based on the defined basis
      tem.X <- trio.Xmat(x,mt,k1,k2,xk1=xk1,xk2=xk2,cenv=cenv)#takes about 0.03 seconds
      ## pen.scale: TRUE if the user wants penalty matrices to be scaled, otherwise FALSE (see mgcv's gam())
      tem.S <- trio.Smat(xk1=xk1,xk2=xk2,pen.scale=TRUE,X=tem.X,mt=mt) 
      ## n.consh is the number of identifiabilty constraints on fh(e), following Wood (2006)
      n.cons1 = 1; n.cons2 = 0 #
      C1 <- matrix(colSums(tem.X),nrow=1)
      qrc1 <- qr(t(C1))
      ZSZ1 <- qr.qty(qrc1,tem.S)[(n.cons1 + 1):(k1),]
      S1 <- t(qr.qty(qrc1,t(ZSZ1))[(n.cons1 + 1):k1,])
      S[2:k1,2:k1] <- S1
      X1 <- t(qr.qy( qrc1,t(tem.X[,1:k1]) )[(n.cons1 + 1):k1,])
      X[,-1] <- X1
      X[h.ind1,1] <- 1
      
      qrc2 <- NULL
      qrc = list(qrc1=qrc1, qrc2=qrc2)
      rm(tem.X,tem.S,X1,C1,ZSZ1,qrc1,qrc2)
    }#if(penmod=="dominant")
    else if(penmod=="recessive"){
      ## Setting up design and penalty matrices based on the defined basis
      tem.X <- trio.Xmat(x,mt,k1,k2,xk1=xk1,xk2=xk2,cenv=cenv)#takes about 0.03 seconds      
      ## pen.scale: TRUE if the user wants penalty matrices to be scaled, otherwise FALSE (see mgcv's gam())
      tem.S <- trio.Smat(xk1=xk1,xk2=xk2,pen.scale=TRUE,X=tem.X,mt=mt) 
      ## n.consh is the number of identifiability constraints on fh(e) - following Wood (2006)
      n.cons1 <- 0; n.cons2 <- 1
      C2 <- matrix(colSums(tem.X),nrow=1)
      qrc2 <- qr(t(C2))
      ZSZ2 <- qr.qty(qrc2,tem.S)[(n.cons2 + 1):(k2),]
      S2 <- t(qr.qty(qrc2,t(ZSZ2))[(n.cons2 + 1):k2,])
      S[2:k2,2:k2] <- S2
      X2 <- t(qr.qy( qrc2,t(tem.X[,1:k2]) )[(n.cons2 + 1):k2,])
      X[,-1] <- X2
      X[h.ind2,1] <- 1
      
      qrc1 <- NULL
      qrc = list(qrc1=qrc1, qrc2=qrc2)
      rm(tem.X,tem.S,X2,C2,ZSZ2,qrc1,qrc2)
    }#else if(penmod=="recessive")
    else {#additive
      n.cons1 <- n.cons2 <- 1
      tem.X <- trio.Xmat.mul(x=x,xk=knots[[1]],cenv=cenv)
      tem.S <- trio.Smat.mul(xk=knots[[1]],X=tem.X)
      C <- matrix(colSums(tem.X),nrow=1)
      qrc <- qr(t(C))
      ZSZ <- qr.qty(qrc,tem.S)[(n.cons1+1):k1,]
      new.S <- t(qr.qty(qrc,t(ZSZ))[(n.cons1+1):k1,])
      S[(n.cons1+1):k1,(n.cons1+1):k1] <- new.S
      new.X <- t(qr.qy( qrc,t(tem.X[,1:k1]) )[(n.cons1+1):k1,])
      X[,-1] <- new.X
      X[,1] <- 1
      
      rm(tem.X,tem.S,new.X,new.S,C,ZSZ)
    }#else: additive
    n.cons <- c(n.cons1, n.cons2)
  }#if(!is.null(penmod))
  
  else{#if codominant
    ## SMOOTH 1
    ## Setting up design and penalty matrices based on the defined basis
    tem.X <- trio.Xmat(x,mt,k1,k2,xk1=xk1,xk2=xk2,cenv=cenv)#takes about 0.03 seconds
    ## pen.scale: TRUE if the user wants penalty matrices to be scaled, otherwise FALSE (see mgcv's gam())
    tem.S <- trio.Smat(xk1=xk1,xk2=xk2,pen.scale=TRUE,X=tem.X,mt=mt) 
    n.cons <- c(1,1)
    n.cons1 <- n.cons[1]
    if(trace){
      cat("using the centering constraints for GRR1 and GRR2. \n")
    }
    C1 <- matrix(colSums(tem.X[h.ind1,1:k1]),nrow=1)
    qrc1 <- qr(t(C1))
    ## t(Q[,-m])%*%S = (t(Q)%*%X)[-m,]
    ZSZ1 <- qr.qty(qrc1,tem.S[1:k1,1:k1])[(n.cons1 + 1):(k1),]
    ## t( ( t(Q)%*% t(t(Q[,-m]%*%S)) )[-m,] ) = t(Q[,-m])%*%S %*% Q[,-m]
    S1 <- t(qr.qty(qrc1,t(ZSZ1))[(n.cons1 + 1):k1,])
    S[2:k1,2:k1] <- S1
    ## X%*%Q[,-m] = t((Q%*%X)[-m,]) 
    X1 <- t(qr.qy( qrc1,t(tem.X[h.ind1,1:k1]) )[(n.cons1 + 1):k1,])
    X[h.ind1,-c(n.cons1,(k1+1):(k1+k2))] <- X1
    X[h.ind1,1] <- 1
    
    ## SMOOTH 2
    n.cons2 <- n.cons[2]
    C2 <- matrix(colSums(tem.X[h.ind2,(k1+1):(k1+k2)]),nrow=1)
    qrc2 <- qr(t(C2));
    ZSZ2 <- qr.qty(qrc2,tem.S[(k1+1):(k1+k2),(k1+1):(k1+k2)])[(n.cons2 + 1):(k2),]
    S2 <- t(qr.qty(qrc2,t(ZSZ2))[(n.cons2 + 1):k2,])
    S[(k1+2):(k1+k2),(k1+2):(k1+k2)] <- S2
    X2 <- t(qr.qy(qrc2,t(tem.X[h.ind2,(k1+1):(k1+k2)]))[(n.cons2 + 1):k2, ])
    X[h.ind2, -c((k1+n.cons2),1:k1)] <- X2
    X[h.ind2, (k1+1)] <- 1
    
    qrc = list(qrc1=qrc1,qrc2=qrc2)
    rm(tem.X,X1,X2,C1,C2,qrc1,qrc2,ZSZ1,ZSZ2)#,tem.S)
  }#if(is.null(penmod))
  
  ## indicator
  const.ind <- rep(TRUE,length(mt)) #all
  ## offsetss for MT3
  offsets <- rep(0,n)
  offsets[mt==3.1] <- log(2)
  offsets[mt==3.2] <- -log(2)
  ind1 <- mt==1; ind2 <- mt==2; ind12 <- ind1|ind2; ind31 <- mt==3.1; ind32 <- mt==3.2
  #assign("ind1",mt==1,envir=.GlobalEnv)
  #assign("ind2",mt==2,envir=.GlobalEnv)
  #assign("ind12",(mt==1|mt==2),envir=.GlobalEnv)
  #assign("ind31",mt==3.1,envir=.GlobalEnv)
  #assign("ind32",mt==3.2,envir=.GlobalEnv)
  n12 = sum(ind12); n3=sum(ind31)
  
  ## to save some time for calculating W,W^(-1) and sqrt(W)
  ## initialize
  init.diagvals = rep(0,sum(const.ind))
  #assign("init.diagvals",rep(0,sum(const.ind)),envir=.GlobalEnv)
  
  n3.double <- n3*2## twice the number of trios with Gp=11
  initmat3 <- (matrix(0,nrow=(n3.double),ncol=(n3.double)))
  #assign("initmat3",(matrix(0,nrow=(n3.double),ncol=(n3.double))),envir=.GlobalEnv)
  diag.index <- matrix(NA,nrow=(n3.double),ncol=2)
  off.index31 <- off.index32 <- matrix(NA,nrow=n3,ncol=2) #half the size of diag.index
  
  ## fill in the diagonal elements
  diag.index[,1] <- diag.index[,2] <- 1:(n3.double)
  #assign("diag.index",diag.index,envir=.GlobalEnv)
  
  ## fill in the off-diagonal elements
  which.31 <- ((1:n3)*2-1)
  #assign("which.31",((1:n3)*2-1),envir=.GlobalEnv) ## odd number
  off.index31[,1] <- off.index32[,2] <- which.31
  off.index31[,2] <- off.index32[,1] <- (which.31+1)
  
  #assign("off.index31",off.index31,envir=.GlobalEnv)
  #assign("off.index32",off.index32,envir=.GlobalEnv)
  #}#if(is.null(penmod))
  
  ## ---------- Fitting begins here ---------- ##
  start = NULL; mustart = NULL; etastart = NULL; weights = 1
  ## FOR NOW THE USER CANNOT PASS THESE VALUES
  ## start: starting value for the parameter vectors
  ## etastart: starting value for the vector of addtive predictors 
  ## mustart: starting value for the vector of means
  ## weights: an optional vector of 'prior weights' to be used in the fitting process
  
  if(is.null(sp)){
    ##** SMOOTHING PARAMETER ESTIMATION **##
    
    ## grid-points are selected using truncated normal distributions based on Duke's estimates;
    ## we assume that the optimum smoothing parameter values are close to the ones estimated by 
    ## Duke's conditional likelihoods (2007)
    if(is.null(lsp.grid)){
      ##based on lsp0
      if(is.null(lsp0)){
        lsp0 = c(0,0) ## initial
        if(any(k==0))## trioplot.res needs a vector for k
          tem.k = rep(k[k!=0][1],2)
        else
          tem.k = k
        
        trioplot.res = trioplot(data,pgenos=pgenos,cgeno=cgeno,
                                cenv=cenv,knots=NULL,k=tem.k,sp=NULL)        
        
        if( !is.null( trioplot.res$gamfit1 ) & !is.null( trioplot.res$gamfit2 ) ) {
          lsp0=log(c(trioplot.res$gamfit1$sp,trioplot.res$gamfit2$sp))
        }
        
        else if( !is.null( trioplot.res$gamfit1 ) & is.null( trioplot.res$gamfit2 ) ) {
          lsp0[1] <- log(trioplot.res$gamfit1$sp)
        }
        
        else if( is.null( trioplot.res$gamfit1 ) & !is.null( trioplot.res$gamfit2 ) ) {
          lsp0[2] <- log(trioplot.res$gamfit1$sp)
        }
        
        else lsp0=lsp0
      }#if(is.null(lsp0))
      
      else{ #if lsp0!=NULL
        trioplot.res=NULL
      }
      
      ## setting up the grids based on lsp0
      n.lsp.grid = 4 ## quartiles
      lsp.prob = 1/n.lsp.grid*(1:(n.lsp.grid-1)) ## whatever stepsize+1 if grid.step=even
      ## exp(min.lsp) and exp(max.lsp) are effectively 0 and infinity, respectively.
      min.lsp = -20; max.lsp = 20 
      
      ## lsp.grid will be a list with two vectors if codominant; otherwise, a vector
      if(penmod == "codominant"){
        ## GRR1
        if(lsp0[1] <= min.lsp)
          lsp1 = min.lsp
        else if(lsp0[1] >= max.lsp)
          lsp1 = max.lsp
        else
          lsp1 = lsp0[1]
        ## sd
        sd1 = max( c( ( max.lsp-lsp1 ), ( lsp1-min.lsp ) ) )
        
        ## grid-points (six or five)
        if((lsp1 == min.lsp) | (lsp1 == max.lsp))
          lsp.grid1 = sort(c(qtnorm(lsp.prob, mean=lsp1, sd=sd1, lower=min.lsp, upper=max.lsp),
                             min.lsp,max.lsp)) ## five grid points
        else 
          lsp.grid1 = sort(c(qtnorm(lsp.prob, mean=lsp1, sd=sd1, lower=min.lsp, upper=max.lsp),
                             min.lsp,max.lsp,lsp1)) ## six grid points
        
        ## GRR2
        if(lsp0[2] <= min.lsp)
          lsp2 = min.lsp
        else if(lsp0[2] >= max.lsp)
          lsp2 = max.lsp
        else
          lsp2 = lsp0[2]
        
        ## sd2
        sd2 = max( c( ( max.lsp-lsp2 ), ( lsp2-min.lsp ) ) )
        
        if((lsp2 == min.lsp) | (lsp2 == max.lsp))
          lsp.grid2 = sort( c(qtnorm(lsp.prob, mean=lsp2, sd=sd2, lower=min.lsp, upper=max.lsp),
                              min.lsp,max.lsp) )
        
        else
          lsp.grid2 = sort( c(qtnorm(lsp.prob, mean=lsp2, sd=sd2, lower=min.lsp, upper=max.lsp),
                              min.lsp,max.lsp,lsp2) )
        
        lsp <- c(lsp1,lsp2)
        sd <- c(sd1,sd2)## may not be necessary
        lsp.grid <- list(lsp.grid1=lsp.grid1, lsp.grid2=lsp.grid2)
      }#if:codominant
      
      else{# not codominant
        if(penmod=="dominant"){
          if(lsp0[1] <= min.lsp)
            lsp = min.lsp
          else if(lsp0[1] >= max.lsp)
            lsp = max.lsp
          else
            lsp = lsp0[1]
        }
        else if(penmod=="additive") {
          ## take the estimate estimated with a larger sample size
          if(sum(ind1) > sum(ind2)){ #n1 > n2, so take the first element
            if(lsp0[1] <= min.lsp)
              lsp = min.lsp
            else if(lsp0[1] >= max.lsp)
              lsp = max.lsp
            else
              lsp = lsp0[1]
          }
          
          else {#n1 <= n2, so take the second element
            if(lsp0[2] <= min.lsp)
              lsp = min.lsp
            else if(lsp0[2] >= max.lsp)
              lsp = max.lsp
            else
              lsp = lsp0[2]
          }
        }#else if: additive
        else {#(penmod=="recessive"){
          ## sp2
          ## mean
          if(lsp0[2] <= min.lsp)
            lsp = min.lsp
          else if(lsp0[2] >= max.lsp)
            lsp = max.lsp
          else
            lsp = lsp0[2]
        }#else: recessive
        sd = max( c( ( max.lsp-lsp ), ( lsp-min.lsp ) ) )
        
        ## grid-points (six or five)
        if((lsp == min.lsp) | (lsp == max.lsp)){
          lsp.grid = sort(c(qtnorm(lsp.prob, mean=lsp, sd=sd, lower=min.lsp, upper=max.lsp),
                            min.lsp,max.lsp)) ## five grid points
        }
        else{
          lsp.grid = sort(c(qtnorm(lsp.prob, mean=lsp, sd=sd, lower=min.lsp, upper=max.lsp),
                            min.lsp,max.lsp,lsp)) ## six grid points
        }
      }#else: not codominant
    }#is.null(lsp.grid)
    
    else{#lsp.grid is provided by the user
      if(penmod == "codominant"){
        stopifnot(is.list(lsp.grid), length(lsp.grid)==2)
      }
      else{#penmod is not codominant
        ## take vectors only: !is.list(lsp.grid)
        stopifnot(!is.list(lsp.grid),is.vector(lsp.grid))
      }      
    }#else (i.e., !if(is.null(lsp.grid))) ends
    
    ## grid search bigins here: at the end, we will have 'sp'
    null.coef=rep(0,bs.dim)
    if(penmod == "codominant" ){      
      ubre.val1 <- conv1 <- rep(NA,(length(lsp.grid1)))
      ubre.val2 <- conv2 <- rep(NA,(length(lsp.grid2)))
      
      ## sp1: fixed sp2 and search over lsp.grid1
      lsp2 <-  lsp0[2]
      for(i in 1:length(lsp.grid1)){
        lsp1 = lsp.grid1[i]
        if(trace)
          cat("\niteration - ",i,": sp1 is ",exp(lsp1),",for a fixed sp2: ",exp(lsp2),"\n",sep="")
        fit.res = triogam.fit(y=y,X=X,S=S,lsp=c(lsp1,lsp2),
                              n=sum(const.ind),n12=n12,n3=n3,mt=mt,offsets=offsets,
                              start=start,mustart=mustart,etastart=etastart,
                              weights=weights,null.coef=null.coef,control=control,k=k,
                              penmod=penmod)
        
        ubre.val1[i] = fit.res$ubre
        conv1[i]=fit.res$conv
        ## resulting sp1-estimate
      }#for-loop ends here
      lsp1 = lsp.grid1[which.min(ubre.val1)]
      if(!all(conv1)){
        warning("P-IRLS Algorithm has not converged for all smoothing-parameter grid-points for \'sp1\'; 
                will try another other grid searcg with a differet value of \'sp2\' .")
        for(i in 1:length(lsp.grid1)){
          lsp1 = lsp.grid1[i]; lsp2 = lsp.grid2[length(lsp.grid2)]
          if(trace)
            cat("\niteration - ",i,": sp1 is ",exp(lsp1),",for a fixed sp2: ",exp(lsp2),"\n",sep="")
          fit.res = triogam.fit(y=y,X=X,S=S,lsp=c(lsp1,lsp2),
                                n=sum(const.ind),n12=n12,n3=n3,mt=mt,offsets=offsets,
                                start=start,mustart=mustart,etastart=etastart,
                                weights=weights,null.coef=null.coef,control=control,k=k,
                                penmod=penmod)
          
          ubre.val1[i] = fit.res$ubre
          ## resulting sp1-estimate
        }#for-loop ends here
        
        lsp1 = lsp.grid1[which.min(ubre.val1)]
      }#if(!all(conv1))-ends
      ubre = min(ubre.val1)
      
      ## estimating sp2: fixing sp1
      for(i in 1:length(lsp.grid2)){
        lsp2 = lsp.grid2[i]
        if(trace)
          cat("\niteration - ",i,": sp2 is ",exp(lsp1),",for a fixed sp1: ",exp(lsp2),"\n",sep="")
        
        fit.res = triogam.fit(y=y,X=X,S=S,lsp=c(lsp1,lsp2),
                              n=n,n12=n12,n3=n3,mt=mt,offsets=offsets,
                              start=start,mustart=mustart,etastart=etastart,
                              weights=weights,null.coef=null.coef,control=control,
                              k=k,penmod=penmod)
        ubre.val2[i] = fit.res$ubre
        
        conv2[i]=fit.res$conv
      }#for(i in 1:length(lsp.grid2))-loop ends
      
      lsp2 = lsp.grid2[which.min(ubre.val2)]
      if(!all(conv2)) 
        stop("P-IRLS Algorithm has not converged: No output will be produced.\n")
      ubre = min(ubre.val2)       
      lsp <- c(lsp1,lsp2)
      ubre.val=list(ubre.val1=ubre.val1,ubre.val2=ubre.val2) ##just for calcultion check
      
      #debugging
{
        conv.list = list(conv1=conv1, conv2=conv2)
      }
  }#penmod codominant
    else{#penmod is not codominant
      ubre.val <- conv <- rep(NA,length(lsp.grid))
      for(i in 1:length(lsp.grid)){
        if(trace)
          cat("\niteration - ",i,": sp is ",exp(lsp),"\n",sep="")
        
        fit.res = triogam.fit(y=y,X=X,S=S,lsp=lsp.grid[i],
                              n=n,n12=n12,n3=n3,mt=mt,offsets=offsets,
                              start=start,mustart=mustart,etastart=etastart,
                              weights=weights,null.coef=null.coef,control=control,
                              k=k,penmod=penmod)
        ubre.val[i] = fit.res$ubre
        ## debugging
{
  conv[i]=fit.res$conv
}
      }#for-loop ends
      lsp = lsp.grid[which.min(ubre.val)]
      ubre = min(ubre.val)
      
      ##debugging
{
        conv.list = list(conv=conv)
      } 
      ubre.val <- list(ubre.val=ubre.val)## make it a list-object
      
    }#else: not codominant
    
    ## sp is obtained by exponentiating 
    sp <- exp(lsp)
    
    ## When doing the permutation test, when the UBRE values are all Inf,
    ## we discard the data:
    if ( testGxE ){
      ## testing
      if(penmod=="codominant"){
        test.ind <- c(length(ubre.val[[1]]),length(ubre.val[[2]]))>0
        for(i in which(test.ind)){
          if ( all( ubre.val[[i]] == Inf) ){
            cat("\n\n All ubre values for either f1 or f2 are Inf. \n\n")
            ret = NULL
            return(ret)     
            break
            ## return nothing!
          }        
        }
      }
      else{
        test.ind <- (length(ubre.val[[1]]) > 0)
        for(i in which(test.ind)){
          if ( all( ubre.val[[1]] == Inf) ){
            cat("\n\n All ubre values for f are Inf. \n\n")
            ret = NULL
            return(ret)     
            break
            ## return nothing!
          }       
        }
      }
    }
    
    ## do fitting again...for getting the beta coef est. given the sp est minimizing ubre:
    ## initialization:
    null.coef=rep(0,bs.dim)
    ## final fit is done with the sp minimizing sp (i.e., the one obtained above)
    final.fit = triogam.fit(y=y,X=X,S=S,lsp=log(sp),
                            n=sum(const.ind),n12=n12,n3=n3,mt=mt,offsets=offsets,
                            start=start,mustart=mustart,etastart=etastart,
                            weights=weights,null.coef=null.coef,
                            control=control,k=k,penmod=penmod,is.final=TRUE)
    sp.user=FALSE
  }#is.null(sp) ends
  
  else{ #i.e., !is.null(sp)
    ## initialization:
    trioplot.res=NULL
    ubre.val = NULL #just for returning
    lsp=log(sp) ## the provided sp
    if(trace)
    {
      cat("\n","sp is ",exp(lsp[1]),",",exp(lsp[2]),"\n",sep="")
    }
    
    null.coef=rep(0,bs.dim)
    final.fit = triogam.fit(y=y,X=X,S=S,lsp=lsp,
                            n=sum(const.ind),n12=n12,n3=n3,mt=mt,offsets=offsets,
                            start=start,mustart=mustart,etastart=etastart,weights=weights,
                            null.coef=null.coef,control=control,k=k,penmod=penmod,is.final=TRUE)
    sp.user=TRUE
  }#else: i.e., !is.null(sp) ends
  
  ## SVD results to be used for calculating
  ## Bayesian covariance matrix for the parameters (p.185 & pp.189--196)
  ## Calculate the other quantities of interest: trace(A),Ve,Vp,the
  ## frequentist- and Bayesian- variance and covariance matrix from
  ## the final fit
  
  R = qr.R(final.fit$qrx)
  Q = qr.Q(final.fit$qrx)
  u1 = final.fit$u1 #(K1+K2)-by-(K1+K2) matrix
  sqrt.A = t(Q%*%u1)
  A= t(sqrt.A)%*%sqrt.A
  
  vt = final.fit$vt
  d.values = final.fit$d.values
  piv.order = order(final.fit$qrx$pivot)
  
  ## Vp: calculate Bayesian posterior covariance matrix for the parameters (p185,p.189)
  ## Ve: calculate frequentist covariance matrix for the parameter estimates (p.189)
  ## sqrt(W)X = QR --> X'WX = (sqrt(W)X)'(sqrt(W)X) = R'Q'QR = R'R
  RtR <- crossprod(R)
  
  Vp <- crossprod((1/d.values)*vt)
  Vp.inverse.debug <- crossprod((d.values)*vt)
  
  ## Calculate the matrix F mapping the un-penalized estimates to the penalized
  ## ones (see Wood 2006,page 171)
  Fmat = Vp%*%RtR
  Ve <- Fmat%*%Vp 
  
  ## inversing permutation applied due to pivoted QR-factorization on sqrt(W)*X
  ## Covarinace matrix
  Vp <- Vp[piv.order,piv.order] #re-shuffle (compare with (X'WX+sp*S)^(-1)X'WX(X'WX+sp*S)^(-1))
  X.tilde = final.fit$X.tilde; S.tilde = final.fit$S.tilde
  
  ## inverse of var-cov matrix for smoothing coefficients
  if(penmod == "codominant")
    param.terms = c(1,(k1+1)) ## parametric terms
  else
    param.terms=1
  inverse.Vp = t( X.tilde ) %*% X.tilde + S.tilde
  
  Ve <- Ve[piv.order,piv.order] #re-shuffle 
  Fmat = Fmat[piv.order,piv.order] #re-shuffle
  
  coef <- final.fit$coef
  
  ## gathering the results
  if(penmod == "codominant")
    names(coef)=c("(Intercept)",paste("s","(",cenv,1,").",c(1:(k[1]-n.cons1)),sep=""),
                  "(Intercept)",paste("s","(",cenv,2,").",c(1:(k[2]-n.cons2)),sep=""))
  else{
    names(coef) <- c("(Intercept)",paste("s","(",cenv,").",c(1:(k[k!=0][1]-1)),sep=""))
  }
  
  
    object$coefficients <- coef
    object$control <- control
    object$data <- data ## original data
    ##edf: diagonal elements of F matrix mapping the unpenalized estimates to 
    ##penalized ones (Wood 2006,page 171)
    object$edf = diag(Fmat)
    object$Gp = attr(triodat,"Gp")
    object$lsp0 <- lsp0
    object$lsp.grid <- lsp.grid
    object$penmod <- penmod
    object$pirls.iter <- final.fit$prils.iter 
    object$qrc <- qrc
    object$smooth = list(model.mat = X, pen.mat = S, bs.dim = k, knots = knots)
    object$sp <- sp
    object$sp.user <- sp.user
    object$terms <- list(cgeno=cgeno,pgenos=pgenos,cenv=cenv)
    object$triodata <- triodat
    object$ubre=final.fit$ubre
    object$ubre.val <- ubre.val ## may not be necessary
    object$Vp <- Vp
  
  
  class(object) <- "trioGxE"
  object
}
