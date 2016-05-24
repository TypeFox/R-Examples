## 15.09.2005 alexandros 


setGeneric("csi", function(x, y, kernel="rbfdot",kpar=list(sigma=0.1), rank, centering = TRUE, kappa =0.99 ,delta = 40 ,tol = 1e-4) standardGeneric("csi"))
setMethod("csi",signature(x="matrix"),
function(x, y, kernel="rbfdot",kpar=list(sigma=0.1), rank, centering = TRUE, kappa =0.99 ,delta = 40 ,tol = 1e-5)
  {
    ## G,P,Q,R,error1,error2,error,predicted.gain,true.gain 
     ## INPUT
    ## x  : data
    ## y  : target vector n x d
    ## m  : maximal rank
    ## kappa : trade-off between approximation of K and prediction of y (suggested: .99)
    ## centering : 1 if centering, 0 otherwise (suggested: 1)
    ## delta : number of columns of cholesky performed in advance (suggested: 40)
    ## tol : minimum gain at iteration (suggested: 1e-4)
    
    ## OUTPUT
    ## G : Cholesky decomposition -> K(P,P) is approximated by G*G'
    ## P : permutation matrix
    ## Q,R : QR decomposition of G (or center(G) if centering)
    ## error1 : tr(K-G*G')/tr(K) at each step of the decomposition
    ## error2 : ||y-Q*Q'*y||.F^2 / ||y||.F^2 at each step of the decomposition
    ## predicted.gain : predicted gain before adding each column
    ## true.gain : actual gain after adding each column

    
    n <- dim(x)[1]
    d <- dim(y)[2]
    if(n != dim(y)[1]) stop("Labels y and data x dont match") 
    
    if(!is(kernel,"kernel"))
      {
        if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
        kernel <- do.call(kernel, kpar)
      }
    if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")
    
    m <- rank
    ## make sure rank is smaller than n 
    m <- min(n-2,m)

    G <- matrix(0,n,min(m+delta,n))    ## Cholesky factor
    diagK <- rep(drop(kernel(x[1,],x[1,])),n)
    P <- 1:n ## pivots
    Q <- matrix(0,n,min(m+delta,n))   ## Q part of the QR decomposition
    R <- matrix(0,min(m+delta,n),min(m+delta,n))   ## R part of the QR decomposition
    traceK <- sum(diagK)
    lambda <- (1-kappa)/traceK 
    if (centering) y <- y - (1/n) * t(matrix(colSums(y),d,n)) 
    sumy2 <- sum(y^2)
    mu <- kappa/sumy2
    error1 <- traceK
    error2 <- sumy2
    predictedgain <- truegain <- rep(0,min(m+delta,n))

    k <- 0            # current index of the Cholesky decomposition
    kadv <- 0         # current index of the look ahead steps
    Dadv <- diagK

    D <- diagK

    ## makes sure that delta is smaller than n - 2
    delta <- min(delta,n - 2)
    ## approximation cost cached quantities
    A1 <- matrix(0,n,1)
    A2 <- matrix(0,n,1)
    A3 <- matrix(0,n,1)
    GTG <- matrix(0,m+delta,m+delta)
    QTy <- matrix(0,m+delta,d)
    QTyyTQ <- matrix(0,m+delta,m+delta)




    ## first performs delta steps of Cholesky and QR decomposition
    if(delta > 0)
      for (i in 1:delta) {
        kadv <- kadv + 1
        ## select best index
        diagmax <- Dadv[kadv]
        jast <- 1
        for (j in 1:(n-kadv+1)) {
          if (Dadv[j+kadv-1] > diagmax/0.99){
            diagmax <- Dadv[j+kadv-1]
            jast <- j
          }
        } 
        if (diagmax < 1e-12){
          kadv <- kadv - 1 ## all pivots are too close to zero, stops 
          break  ## this can only happen if the matrix has rank less than delta
        }
        else{
          jast <- jast + kadv-1
          
          ## permute indices 
          P[c(kadv,jast)]          <- P[c(jast,kadv)] 
          Dadv[c(kadv, jast)]      <- Dadv[c(jast, kadv)] 
          D[c(kadv, jast)]         <- D[c(jast, kadv)] 
          A1[c(kadv, jast)]         <- A1[c(jast, kadv)] 
          G[c(kadv, jast),1:kadv-1] <- G[c(jast,kadv),1:kadv-1]
          Q[c(kadv, jast),1:kadv-1] <- Q[c(jast, kadv),1:kadv-1]
          
          ## compute new Cholesky column
          G[kadv,kadv] <- Dadv[kadv]
          G[kadv,kadv] <- sqrt(G[kadv,kadv])
          newKcol <- kernelMatrix(kernel, x[P[(kadv+1):n],,drop = FALSE],x[P[kadv],,drop=FALSE])
          G[(kadv+1):n,kadv]<- (1/G[kadv,kadv])*(newKcol - G[(kadv+1):n,1:kadv-1,drop=FALSE] %*% t(G[kadv,1:kadv-1,drop=FALSE]))
          
          ## update diagonal
          Dadv[(kadv+1):n] <-  Dadv[(kadv+1):n] - G[(kadv+1):n,kadv]^2
          Dadv[kadv] <- 0

          ## performs QR
          if (centering)
            Gcol <- G[,kadv,drop=FALSE] - (1/n) * matrix(sum(G[,kadv]),n,1)
          else
            Gcol <- G[,kadv, drop=FALSE]
          
          R[1:kadv-1,kadv] <- crossprod(Q[,1:kadv-1, drop=FALSE], Gcol)
          Q[,kadv] <- Gcol - Q[,1:kadv-1,drop=FALSE] %*% R[1:kadv-1,kadv,drop=FALSE]
          R[kadv,kadv] <- sqrt(sum(Q[,kadv]^2))
          Q[,kadv] <- Q[,kadv]/drop(R[kadv,kadv])
          
          ## update cached quantities
          if (centering)
            GTG[1:kadv,kadv] <- crossprod(G[,1:kadv], G[,kadv])
          else
            GTG[1:kadv,kadv] <- crossprod(R[1:kadv,1:kadv], R[1:kadv,kadv])
          

          GTG[kadv,1:kadv] <- t(GTG[1:kadv,kadv])
          QTy[kadv,] <- crossprod(Q[,kadv], y[P,,drop = FALSE])
          QTyyTQ[kadv,1:kadv] <- QTy[kadv,,drop=FALSE] %*% t(QTy[1:kadv,,drop=FALSE])
          QTyyTQ[1:kadv,kadv] <- t(QTyyTQ[kadv,1:kadv])
          
          ## update costs
          A1[kadv:n] <- A1[kadv:n] + GTG[kadv,kadv] * G[kadv:n,kadv]^2 
          A1[kadv:n] <- A1[kadv:n] + 2 * G[kadv:n,kadv] *(G[kadv:n,1:kadv-1] %*% GTG[1:kadv-1,kadv,drop=FALSE])
        } 
      }

    ## compute remaining costs for all indices
    A2 <- rowSums(( G[,1:kadv,drop=FALSE] %*% crossprod(R[1:kadv,1:kadv], QTy[1:kadv,,drop=FALSE]))^2)
    A3 <- rowSums((G[,1:kadv,drop=FALSE] %*% t(R[1:kadv,1:kadv]))^2)

    ## start main loop
    while (k < m){
      k <- k +1

      ## compute the gains in approximation for all remaining indices
      dJK <- matrix(0,(n-k+1),1)

      for (i in 1:(n-k+1)) {
        kast <- k+i-1
        
        if (D[kast] < 1e-12)
          dJK[i] <- -1e100   ## this column is already generated by already
                             ## selected columns -> cannot be selected
        else {
          dJK[i] <- A1[kast]
          
          if (kast > kadv)
            ## add eta
            dJK[i] <- dJK[i] + D[kast]^2 - (D[kast] - Dadv[kast])^2
          dJK[i] <- dJK[i] / D[kast]
        } 
      }
      dJy <- matrix(0,n-k+1,1)
      
      if (kadv > k){
        for (i in 1:(n-k+1)) {
          kast <- k+i-1
          if (A3[kast] < 1e-12)
            dJy[i] <- 0
          else
            dJy[i] <- A2[kast] / A3[kast]
        } 
      } 

   
      ## select the best column
      dJ <- lambda * dJK + mu * dJy   
      diagmax <- -1
      jast <- 0
      
      for (j in 1:(n-k+1)) {
        if (D[j+k-1] > 1e-12)
          if (dJ[j] > diagmax/0.9){
            diagmax <- dJ[j]
            jast <- j
          }
      } 
  
  
      if (jast==0) {
        ## no more good indices, exit
        k <- k-1 
        break 
      }  
      
      jast <- jast + k - 1
      predictedgain[k] <- diagmax
      
      ## performs one cholesky + QR step:
      ## if new pivot not already selected, use pivot
      ## otherwise, select new look ahead index that maximize Dadv
      
      if (jast > kadv){
        newpivot <- jast
        jast <- kadv + 1
      }
      else{
        a <- 1e-12
        b <- 0
        for (j in 1:(n-kadv)) {
          if (Dadv[j+kadv] > a/0.99){
            a <- Dadv[j+kadv]
            b <- j+kadv
          } 
        }

        if (b==0)
          newpivot <- 0
        else
          newpivot <- b
      }
      
      
      if (newpivot > 0){
        ## performs steps
        kadv <- kadv + 1
        ## permute
        P[c(kadv, newpivot)]  <- P[c(newpivot, kadv)] 
        Dadv[c(kadv, newpivot)]  <- Dadv[c(newpivot, kadv)] 
        D[c(kadv, newpivot)]  <- D[c(newpivot, kadv)] 
        A1[c(kadv, newpivot)]  <- A1[c(newpivot, kadv)] 
        A2[c(kadv, newpivot)]  <- A2[c(newpivot, kadv)] 
        A3[c(kadv, newpivot)]  <- A3[c(newpivot, kadv)] 
        G[c(kadv, newpivot),1:kadv-1] <- G[c(newpivot, kadv),1:kadv-1]
        Q[c(kadv, newpivot),1:kadv-1] <- Q[ c(newpivot, kadv),1:kadv-1]
    
        ## compute new Cholesky column

        G[kadv,kadv] <- Dadv[kadv]
        G[kadv,kadv] <- sqrt(G[kadv,kadv])
        newKcol <- kernelMatrix(kernel,x[P[(kadv+1):n],,drop=FALSE],x[P[kadv],,drop=FALSE])
        G[(kadv+1):n,kadv] <- 1/G[kadv,kadv]*( newKcol - G[(kadv+1):n,1:kadv-1,drop=FALSE]%*%t(G[kadv,1:kadv-1,drop=FALSE]))
        

        ## update diagonal
        Dadv[(kadv+1):n] <-  Dadv[(kadv+1):n] - G[(kadv+1):n,kadv]^2
        Dadv[kadv] <- 0
    
        ## performs QR
        if (centering)
          Gcol <- G[,kadv,drop=FALSE] - 1/n * matrix(sum(G[,kadv]),n,1 )
        else
          Gcol <- G[,kadv,drop=FALSE]
    
        R[1:kadv-1,kadv] <- crossprod(Q[,1:kadv-1], Gcol)
        Q[,kadv] <- Gcol - Q[,1:kadv-1, drop=FALSE] %*% R[1:kadv-1,kadv, drop=FALSE]
        R[kadv,kadv] <- sum(abs(Q[,kadv])^2)^(1/2)
        Q[,kadv] <- Q[,kadv] / drop(R[kadv,kadv])
        
        ## update the cached quantities
        if (centering)
          GTG[k:kadv,kadv] <- crossprod(G[,k:kadv], G[,kadv])
        else
          GTG[k:kadv,kadv] <- crossprod(R[1:kadv,k:kadv], R[1:kadv,kadv])
        
        GTG[kadv,k:kadv] <- t(GTG[k:kadv,kadv])
        QTy[kadv,] <- crossprod(Q[,kadv], y[P,,drop =FALSE])
        QTyyTQ[kadv,k:kadv] <- QTy[kadv,,drop = FALSE] %*% t(QTy[k:kadv,,drop = FALSE])
        QTyyTQ[k:kadv,kadv] <- t(QTyyTQ[kadv,k:kadv])
        
        ## update costs
        A1[kadv:n] <- A1[kadv:n] + GTG[kadv,kadv] * G[kadv:n,kadv]^2 
        A1[kadv:n] <- A1[kadv:n] + 2 * G[kadv:n,kadv] * (G[kadv:n,k:kadv-1,drop = FALSE] %*% GTG[k:kadv-1,kadv,drop=FALSE])
        A3[kadv:n] <- A3[kadv:n] + G[kadv:n,kadv]^2 * sum(R[k:kadv,kadv]^2)
        temp <- crossprod(R[k:kadv,kadv,drop = FALSE], R[k:kadv,k:kadv-1,drop = FALSE])
        A3[kadv:n] <- A3[kadv:n] + 2 *  G[kadv:n,kadv] * (G[kadv:n,k:kadv-1] %*% t(temp))
        temp <- crossprod(R[k:kadv,kadv,drop = FALSE], QTyyTQ[k:kadv,k:kadv,drop = FALSE])
        temp1 <- temp %*% R[k:kadv,kadv,drop = FALSE] 
        A2[kadv:n] <- A2[kadv:n] + G[kadv:n,kadv,drop = FALSE]^2 %*% temp1
        temp2 <- temp %*% R[k:kadv,k:kadv-1]
        A2[kadv:n] <- A2[kadv:n] + 2 * G[kadv:n,kadv] * (G[kadv:n,k:kadv-1,drop=FALSE] %*% t(temp2))
      } 
      ## permute pivots in the Cholesky and QR decomposition between p,q
      p <- k
      q <- jast
      if (p < q){
    
        ## store some quantities
        Gbef <- G[,p:q]
        Gbeftotal <- G[,k:kadv]
        GTGbef <- GTG[p:q,p:q]
        QTyyTQbef <- QTyyTQ[p:q,k:kadv]
        Rbef <- R[p:q,p:q]
        Rbeftotal <- R[k:kadv,k:kadv]
        tempG <- diag(1,q-p+1,q-p+1)
        tempQ <- diag(1,q-p+1,q-p+1)
    
        for (s in seq(q-1,p,-1)) {

          ## permute indices
          P[c(s, s+1)]  <- P[c(s+1, s)] 
          Dadv[c(s, s+1)]  <- Dadv[c(s+1, s)] 
          D[c(s, s+1)]  <- D[c(s+1, s)] 
          A1[c(s, s+1)]  <- A1[c(s+1, s)] 
          A2[c(s, s+1)]  <- A2[c(s+1, s)] 
          A3[c(s, s+1)]  <- A3[c(s+1, s)] 
          G[c(s, s+1),1:kadv] <- G[c(s+1,s), 1:kadv]
          Gbef[c(s, s+1),] <- Gbef[c(s+1, s),]
          Gbeftotal[c(s, s+1),] <- Gbeftotal[c(s+1, s),]
          Q[c(s, s+1),1:kadv] <- Q[c(s+1, s) ,1:kadv]
      
          ## update decompositions
          res <- .qr2(t(G[s:(s+1),s:(s+1)]))
          Q1 <- res$Q
          R1 <- res$R
          G[,s:(s+1)] <- G[,s:(s+1)] %*% Q1
          G[s,(s+1)] <- 0
          R[1:kadv,s:(s+1)] <- R[1:kadv,s:(s+1)] %*% Q1
          res <- .qr2(R[s:(s+1),s:(s+1)])
          Q2 <- res$Q
          R2 <- res$R
          R[s:(s+1),1:kadv] <- crossprod(Q2, R[s:(s+1),1:kadv])
          Q[,s:(s+1)] <- Q[,s:(s+1)] %*% Q2
          R[s+1,s] <- 0
          
          ## update relevant quantities
          if( k <= (s-1) && s+2 <= kadv)
            nonchanged <- c(k:(s-1), (s+2):kadv)
          if( k <= (s-1) && s+2 > kadv)
            nonchanged <- k:(s-1)
          if( k > (s-1) && s+2 <= kadv)
            nonchanged <- (s+2):kadv
      
          GTG[nonchanged,s:(s+1)] <- GTG[nonchanged,s:(s+1)] %*% Q1
          GTG[s:(s+1),nonchanged] <- t(GTG[nonchanged,s:(s+1)])
          GTG[s:(s+1),s:(s+1)] <- crossprod(Q1, GTG[s:(s+1),s:(s+1)] %*% Q1)
          QTy[s:(s+1),] <- crossprod(Q2, QTy[s:(s+1),])
          QTyyTQ[nonchanged,s:(s+1)] <- QTyyTQ[nonchanged,s:(s+1)] %*% Q2
          QTyyTQ[s:(s+1),nonchanged] <- t(QTyyTQ[nonchanged,s:(s+1)])
          QTyyTQ[s:(s+1),s:(s+1)] <- crossprod(Q2, QTyyTQ[s:(s+1),s:(s+1)] %*% Q2)
          tempG[,(s-p+1):(s-p+2)] <- tempG[,(s-p+1):(s-p+2)] %*% Q1  
          tempQ[,(s-p+1):(s-p+2)] <- tempQ[,(s-p+1):(s-p+2)] %*% Q2
        } 
        
        ## update costs
        tempG <- tempG[,1]
        tempGG <- GTGbef %*% tempG 
        A1[k:n] <- A1[k:n] - 2 * G[k:n,k] * (Gbef[k:n,] %*% tempGG)                # between p and q -> different

        if(k > (p-1) )
          kmin <- 0
        else kmin <- k:(p-1)

        if((q+1) > kadv)
          qmin <- 0
        else qmin <- (q+1):kadv
    
        A1[k:n] <- A1[k:n] - 2 * G[k:n,k] * (G[k:n,kmin,drop=FALSE] %*% GTG[kmin,k,drop=FALSE])         # below p
        A1[k:n] <- A1[k:n] - 2 * G[k:n,k] * (G[k:n,qmin,drop=FALSE] %*% GTG[qmin,k,drop=FALSE])   # above q
        tempQ <- tempQ[,1]
        temp <- G[k:n,qmin,drop=FALSE] %*% t(R[k,qmin,drop=FALSE])
        temp <- temp + G[k:n,kmin,drop=FALSE] %*% t(R[k,kmin,drop=FALSE])
        
        temp <- temp + Gbef[k:n,] %*% crossprod(Rbef, tempQ)
        A3[k:n] <- A3[k:n] - temp^2 
        A2[k:n] <- A2[k:n] + temp^2 * QTyyTQ[k,k]
        temp2 <- crossprod(tempQ,QTyyTQbef) %*% Rbeftotal
        A2[k:n] <- A2[k:n] - 2 * temp * (Gbeftotal[k:n,,drop=FALSE] %*% t(temp2))
      }
      else
        { 
          ## update costs
          A1[k:n] <- A1[k:n] - 2 * G[k:n,k] * (G[k:n,k:kadv,drop=FALSE] %*% GTG[k:kadv,k,drop=FALSE]) 
          A3[k:n]<- A3[k:n] - (G[k:n,k:kadv,drop=FALSE] %*% t(R[k,k:kadv,drop=FALSE]))^2
          temp <- G[k:n,k:kadv,drop=FALSE] %*% t(R[k,k:kadv,drop=FALSE])
          A2[k:n] <- A2[k:n] + (temp^2) * QTyyTQ[k,k]
          temp2 <- QTyyTQ[k,k:kadv,drop=FALSE] %*% R[k:kadv,k:kadv,drop=FALSE]
          A2[k:n] <- A2[k:n] - 2 * temp * (G[k:n,k:kadv,drop=FALSE] %*% t(temp2))
        } 
      
      ## update diagonal and other quantities (A1,B1)
      D[(k+1):n] <-  D[(k+1):n] - G[(k+1):n,k]^2
      D[k] <- 0 
      A1[k:n] <- A1[k:n] + GTG[k,k] * (G[k:n,k]^2)

      ## compute errors and true gains
      temp2 <-  crossprod(Q[,k], y[P,])
      temp2 <- sum(temp2^2)
      temp1 <- sum(G[,k]^2)
      truegain[k] <- temp1 * lambda + temp2 * mu
      error1[k+1] <- error1[k] - temp1 
      error2[k+1] <- error2[k] - temp2
      
      if (truegain[k] < tol)
        break  
    }
    
    
    ## reduce dimensions of decomposition
    G <- G[,1:k,drop=FALSE]
    Q <- Q[,1:k,drop=FALSE]
    R <- R[1:k,1:k,drop=FALSE]


    ## compute and normalize errors
    error <- lambda * error1 + mu * error2
    error1 <- error1 / traceK
    error2 <- error2 / sumy2

    repivot <- sort(P, index.return = TRUE)$ix
    
    return(new("csi",.Data=G[repivot, ,drop=FALSE],Q= Q[repivot,,drop = FALSE], R = R, pivots=repivot, diagresidues = error1, maxresiduals = error2, truegain = truegain, predgain = predictedgain))
    
  })

## I guess we can replace this with qr()
.qr2 <- function(M)
  { 
    ## QR decomposition for 2x2 matrices
    Q <- matrix(0,2,2)
    R <- matrix(0,2,2)
    x <- sqrt(M[1,1]^2 + M[2,1]^2)
    R[1,1] <- x
    Q[,1] <- M[,1]/x
    R[1,2] <- crossprod(Q[,1],  M[,2])
    Q[,2] <- M[,2] - R[1,2] * Q[,1]
    R[2,2] <- sum(abs(Q[,2])^2)^(1/2)
    Q[,2] <- Q[,2] / R[2,2]
    
    return(list(Q=Q,R=R))
  }
