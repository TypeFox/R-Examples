
#source R codes:

#histw()

#insert()

#obj.fcn()

#uwham()
#uwham.phi()

#uwham.boot()

############################################

histw <- function(x, w, xaxis, xmin, xmax, ymax, bar=TRUE, add=FALSE, col="black", dens=TRUE) {
     #x: data vector
     #w: inverse weight
     #xaxis: vector of cut points
     #xmin,xmax: the range of x coordinate
     #ymax: the maximum of y coordinate

     #bar: bar plot (if TRUE) or line plot
     #add: if TRUE, the plot is added to an existing plot
     #col: color of lines
     #dens: if TRUE, the histogram has a total area of one

     nbin <- length(xaxis)
     xbin <- cut(x, breaks=xaxis, include.lowest=T, labels=1:(nbin-1))

     y <- tapply(w, xbin, sum)
     y[is.na(y)] <- 0
     y <- y/sum(w)
     if (dens) y <- y/ (xaxis[-1]-xaxis[-nbin])

    if (!add) {
     plot.new()
     plot.window(xlim=c(xmin,xmax), ylim=c(0,ymax))
     axis(1, pos=0)
     axis(2, pos=xmin)
    }

     if (bar==1) {
        rect(xaxis[-nbin], 0, xaxis[-1], y)
     }
     else {
        xval <- as.vector(rbind(xaxis[-nbin],xaxis[-1])) 
        yval <- as.vector(rbind(y,y))
        lines(c(min(xmin,xaxis[1]), xval, max(xmax,xaxis[length(xaxis)])), 
              c(0,yval,0), lty="11", lwd=2, col=col)
     }
     invisible()
}


insert <- function(x, d, x0=0) {
  #This inserts a value x0 at d-th position of x
  
  if (d==1)
   c(x0, x)
  else
   c(x[1:(d-1)], x0, x[-(1:(d-1))])
}

obj.fcn <- function(ze, logQ, size, base) {

   #ze:   a vector of log-normalizing constants (or free energies)
   #logQ: log of the unnormalized density ratios (over the baseline)
   #size: the individual sample sizes for the distributions
   #base: the baseline index

   N <- dim(logQ)[2]
   rho <- size/N

   ze <- insert(ze, base)
 
   #
   Qnorm <- exp(logQ-ze) *rho
   Qsum <- apply(Qnorm, 2, sum)

   val <- sum(log(Qsum)) /N +sum(ze*rho)
   
   W <- t(Qnorm[-base, ,drop=FALSE])/Qsum   #remove one column
   grad <- -apply(W, 2, sum) /N +rho[-base]

   O <- t(W)%*%W /N
   hess <- -O + diag( apply(W, 2, sum), nrow=length(grad) ) /N

   list(value=val, gradient=grad, hessian=hess)
}

uwham <- function(label=NULL, logQ, size=NULL, base=NULL, init=NULL, fisher=TRUE) {

   #label: a vector of labels, indicating which observation is obtained from which distribution
   #logQ: N x M matrix of log unnormalized densities (or negative potential energies), where 
   #      N is the total sample size, i.e., sum(size),
   #      M is the number of distributions for which free energies are to be computed
   #size: a vector of length M, giving the individual sample sizes for the distributions
   #base: the baseline index, between 1 to M, for the distribution whose free energy is set to 0
   #init: a vector of length M, giving the initial values of free energies
   #fisher: logical; if NULL, no variance estimation; 
   #        if TRUE, variance estimation is based on Fisher information

   N <- dim(logQ)[1]
   M <- dim(logQ)[2]

   # compute size if needed
   if (is.null(size)) {
      if (is.null(label)) {
         stop("either label or size must be provided")

      } else {
         size <- c( tapply(1:N, factor(label, levels=1:M), length) )
         size[is.na(size)] <- 0
      }
   }

   #logical, indicating the distributions with observations
   sampled <- size>0

   # check size and logQ
   if (N!=sum(size))
      stop("inconsistent sum(size) and dim(logQ)[1]")

   if (M!=length(size))
      stop("inconsistent length(size) and dim(logQ)[2]")

   #check size and label
   if (!is.null(label)) {
      if ( any(tapply(1:N, label, length)!=size[sampled]) )
         stop("inconsistent label and size[sampled]")
   } else {
      if (!is.null(fisher))
       if (fisher==FALSE) {
          label <- rep((1:M)[sampled], times=size[sampled])
          print("assume that observations are ordered by thermodynamic state if fisher=FALSE and label=NULL")
       }
   }

   #m is the number of distributions from which observations are simulated
   m <- sum(sampled)   #m<=M
   rho <- size[sampled]/N

   #set base or init if not provided
   if (is.null(base)) 
      base <- (1:M)[sampled][1]
   else if (!sampled[base])
      stop("observations from the baseline are required")

   if (is.null(init))
      init <- rep(0,M)

   #the baseline index, between 1 to m, corresponding to the sampled distributions
   base0 <- (1:m)[ as.logical(insert(rep(0, M-1), base, 1)[sampled]) ]

   #log of unnormalized density ratios over the baseline
   logQ <- t(logQ - logQ[,base])

   #use trust
   out <- trust(obj.fcn, init[sampled][-base0], rinit=1, rmax=100, iterlim=1000,
                         logQ=logQ[sampled,], size=size[sampled], base=base0)

   ze0 <- insert(out$argument, base0)

   ze <- rep(0,M)
   ze[sampled] <- ze0

   Qnorm <- exp(logQ-ze)
   Qsum <- apply(Qnorm[sampled,] *rho, 2, sum)

   W <- t(Qnorm) /Qsum
   z <- apply(W, 2, mean)
  
   #all elements should be equal to 1
   check <- z[sampled]

   if (m<M) {
      ze[!sampled] <- log(z[!sampled])
      W[,!sampled] <- t(t(W[,!sampled])/z[!sampled])
   }

   #variance estimation 
  if (is.null(fisher)){

   list(ze=ze,
        W=W, check=check,
        out=out, 
        size=size, base=base)

  } else {
   O <- t(W)%*%W /N

   D <- matrix(0, M,M)
   D[,sampled] <- t( t(O[,sampled])*rho )

   H <- D - diag(1, nrow=M)
   H <- H[-base,-base]

   if (fisher) {
    G <- O - D[,sampled]%*%O[sampled,]
    G <- G[-base, -base]

    iHG <- -O + rep(1,M)%*%O[base, ,drop=F]
    iHG <- iHG[-base,-base]

    Ve <- iHG%*%solve(t(H)) /N

   } else {
    C <- matrix(0, m,M)
    for (j in 1:M)
       C[,j] <- tapply(W[,j], as.factor(label), mean) 
    
    G <- O - t(C)%*%diag(rho)%*%C
    G <- G[-base, -base]
    Ve <- solve(H, G)%*%solve(t(H)) /N
 
    #Equivalent
    #R <- matrix(0, N,M)
    #for (j in 1:M)
    #   R[,j] <- tapply(W[,j], as.factor(label), mean)[label] 
    # 
    #in.fcn <- solve(H, t(W[,-base] - R[,-base]))
    #Ve <- in.fcn%*%t(in.fcn) /N^2
   }

   #variance vector 
   ve <- insert(diag(Ve), base)

   #variance-covariance matrix
   Ve <- apply( apply(Ve, 2, insert, base), 1, insert, base )

   list(ze=ze, ve=ve, Ve=Ve,
        W=W, check=check,
        out=out, 
        label=label, size=size, base=base)
  }
}

uwham.phi <- function(phi, state, out.uwham, fisher=TRUE) {

   N <- dim(out.uwham$W)[1]
   M <- dim(out.uwham$W)[2]

   label <- out.uwham$label
   size <- out.uwham$size
   base <- out.uwham$base

   sampled <- size>0

   m <- sum(sampled)
   rho <- size[sampled]/N

   #
   L <- length(state)
   sampled2 <- c(sampled, rep(FALSE, L))

   W.phi <- out.uwham$W[,state, drop=FALSE] *phi
   phi.bar <- apply(W.phi, 2, mean)

   W <- cbind(out.uwham$W, W.phi)
   rm(W.phi)

   #variance estimation
  if (is.null(fisher)){

   list(phi=phi.bar)

  } else {
   O <- t(W)%*%W /N

   D <- matrix(0, M+L,M+L)
   D[,sampled2] <- t( t(O[,sampled2])*rho )

   H <- D - diag(1, nrow=M+L)
   H <- H[-base,-base]

   if (fisher) {
    G <- O - D[,sampled2]%*%O[sampled2,]
    G <- G[-base, -base]

    iHG <- -O + rep(1,M+L)%*%O[base, ,drop=F]
    iHG <- iHG[-base,-base]

    Ve <- iHG%*%solve(t(H)) /N

   } else {
    C <- matrix(0, m,M+L)
    for (j in 1:(M+L))
       C[,j] <- tapply(W[,j], as.factor(label), mean) 
    
    G <- O - t(C)%*%diag(rho)%*%C
    G <- G[-base, -base]
    Ve <- solve(H, G)%*%solve(t(H)) /N
 
    #Equivalently
    #R <- matrix(0, N,M+L)
    #for (j in 1:(M+L))
    #   R[,j] <- tapply(W[,j], as.factor(label), mean)[label] 
    # 
    #in.fcn <- solve(H, t(W[,-base] - R[,-base]))
    #Ve <- in.fcn%*%t(in.fcn) /N^2
   }

   Ve <- apply( apply(Ve, 2, insert, base), 1, insert, base )
   Ve <- Ve[c(state,M+1:L), c(state,M+1:L)] 

   mat <- cbind(-diag(phi.bar, nrow=L), diag(1, nrow=L))
   phi.V <- mat%*%Ve%*%t(mat)
   phi.v <- diag(phi.V)

   list(phi=phi.bar, phi.V=phi.V, phi.v=phi.v)
  }
}

uwham.boot <- function(proc.type, block.size, boot.size, seed=0, label=NULL, logQ, size=NULL, base=NULL, init=NULL, phi=NULL, state=NULL) {

   #proc.type: type of simulation, 
   #           "indep" for simulation of independent chains, 
   #           "parallel" for parallel tempering, or 
   #           "serial" for serial tempering 
   #block.size: recycled to be a vector of block sizes if proc.type="indep" or 
   #            the first element is treated as a single block size if proc.type="parallel" or "serial"
   #boot.size: the number of bootstrap replications
   #seed: seed for random number generation
   #label: a vector of labels, indicating which observation is obtained from which distribution

   N <- dim(logQ)[1]
   M <- dim(logQ)[2]

   # compute size if needed
   if (is.null(size)) {
      if (is.null(label)) {
         stop("either label or size must be provided")

      } else {
         size <- c( tapply(1:N, factor(label, levels=1:M), length) )
         size[is.na(size)] <- 0
      }
   }

   #logical, indicating the distributions with observations
   sampled <- size>0

   # check size and logQ
   if (N!=sum(size))
      stop("inconsistent sum(size) and dim(logQ)[1]")

   if (M!=length(size))
      stop("inconsistent length(size) and dim(logQ)[2]")

   #check size and label
   if (!is.null(label)) {
      if ( any(tapply(1:N, label, length)!=size[sampled]) )
         stop("inconsistent label and size[sampled]")
   } else {
      if (proc.type=="serial")
         stop("label is required if proc.type='serial'")
   }

   #m is the number of distributions from which observations are simulated
   m <- sum(sampled)

   #set base or init if not provided
   if (is.null(base)) {
      base <- (1:M)[sampled][1]
   } else if (!sampled[base]) {
      stop("observations from the baseline are required")
   }

   if (is.null(init))
      init <- rep(0,M)

   # set block.size
   if (proc.type=="indep") {
      block.size <- rep(block.size, length=m)

   } else if (proc.type=="parallel") {
      n <- N/m

      if (any(size[sampled]!=n))
         stop("equal sample sizes are required if proc.type='parallel'")

      block.size <- block.size[1]
      ind <- matrix(1:n, nrow=block.size)
      off <- rep((1:m)*n-n, each=n)

   } else if (proc.type=="serial") {
      block.size <- block.size[1]
      ind <- matrix(1:N, nrow=block.size)

   } else {
      stop("proc.type must be 'indep', 'parallel', or 'serial'")
   }

   #
   set.seed(seed)

   ans <- matrix(0, boot.size, M)

   if (!is.null(phi))
      ans2 <- matrix(0, boot.size, length(state))  

   for (i in 1:boot.size) {

     if (proc.type=="indep") {
       start <- 0
       sam <- NULL

       for (j in 1:m) {
          n <- size[sampled][j]
          ind <- matrix(1:n, nrow=block.size[j])

          sam2 <- sample(1:(n/block.size[j]), n/block.size[j], replace=T)
          sam2 <- c( ind[, sam2] ) 

          sam <- c(sam, start+sam2)
          start <- n+start
       } 
     } else if (proc.type=="parallel") {
       sam <- sample(1:(n/block.size), n/block.size, replace=T)
       sam <- c( ind[, sam] ) 
       sam <- off + rep(sam, times=m) 

     } else if (proc.type=="serial") {
       sam <- sample(1:(N/block.size), N/block.size, replace=T)
       sam <- c( ind[, sam] ) 

       size2 <- tapply(1:N, label[sam], length)
       if (length(size2)<m)
          stop("No observation is resampled from one or more thermodynamic states")
       else
          size[sampled] <- size2
     }

     out <- uwham(logQ=logQ[sam,], size=size, base=base, init=init, fisher=NULL)
     ans[i,] <- out$ze

     if (!is.null(phi))
        ans2[i,] <- uwham.phi(phi=phi[sam], state=state, out.uwham=out, fisher=NULL)$phi
   }

   boot.ze <- apply(ans, 2, mean)
   boot.ve <- apply(ans, 2, var)

   if (!is.null(phi)) {
      boot.phi <- apply(ans2, 2, mean)
      boot.phi.v <- apply(ans2, 2, var)

      list(ze=boot.ze, ve=boot.ve, 
           phi=boot.phi, phi.v=boot.phi.v)
   } else
      list(ze=boot.ze, ve=boot.ve)
}

#save(list=ls(), file="UWHAM.R")
