##  calculates the kernel maximum mean discrepancy for samples from two distributions
## author: alexandros karatzoglou

setGeneric("kmmd",function(x,...) standardGeneric("kmmd"))
setMethod("kmmd", signature(x = "matrix"),
          function(x, y, kernel="rbfdot",kpar="automatic", alpha = 0.05, asymptotic = FALSE,  replace = TRUE, ntimes = 150, frac = 1,  ...)
          {
            x <- as.matrix(x)
            y <- as.matrix(y)

             res <- new("kmmd")


            if(is.character(kernel)){
              kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot","matrix"))
              
              if(kernel == "matrix")
                if(dim(x)[1]==dim(x)[2])
                  return(kmmd(x= as.kernelMatrix(x), y = y, Kxy = as.kernelMatrix(x)%*%y, alpha = 0.05, asymptotic = FALSE,  replace = TRUE, ntimes = 100, frac = 1,  ...))
                 else
                  stop(" kernel matrix not square!")
              
              if(is.character(kpar))
                if((kernel == "tanhdot" || kernel == "vanilladot" || kernel == "polydot"|| kernel == "besseldot" || kernel== "anovadot"|| kernel=="splinedot") &&  kpar=="automatic" )
                  {
                    cat (" Setting default kernel parameters ","\n")
                    kpar <- list()
                  }
            }

            if (!is.function(kernel))
              if (!is.list(kpar)&&is.character(kpar)&&(kernel == "laplacedot"|| kernel=="rbfdot")){
                kp <- match.arg(kpar,"automatic")
                if(kp=="automatic")
                  kpar <- list(sigma=sigest(rbind(x,y),scaled=FALSE)[2])
                cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n")
                
              }
            if(!is(kernel,"kernel"))
              {
                if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
                kernel <- do.call(kernel, kpar)
              }
            
            if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

            m <- dim(x)[1]
            n <- dim(y)[1]

            N <- max(m,n)
            M <- min(m,n)

            Kxx <- kernelMatrix(kernel,x)
            Kyy <- kernelMatrix(kernel,y)
            Kxy <- kernelMatrix(kernel,x,y)

            resmmd <- .submmd(Kxx, Kyy, Kxy, alpha) 

            H0(res) <- (resmmd$mmd1 > resmmd$D1) 
            Radbound(res) <- resmmd$D1
            Asymbound(res) <- 0
            mmdstats(res)[1] <- resmmd$mmd1
            mmdstats(res)[2] <- resmmd$mmd3
            
            if(asymptotic){
              boundA <- .submmd3bound(Kxx, Kyy, Kxy, alpha, frac, ntimes, replace)

              AsympH0(res) <- (resmmd$mmd3 > boundA) 
              Asymbound(res) <- boundA
            }

            kernelf(res) <- kernel
            return(res)
          })



setMethod("kmmd",signature(x="list"),
               function(x, y, kernel="stringdot",kpar=list(type="spectrum",length=4), alpha = 0.05, asymptotic = FALSE,  replace = TRUE, ntimes = 150, frac = 1,  ...)
  {
    
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  Kxx <- kernelMatrix(kernel,x)
  Kyy <- kernelMatrix(kernel,y)
  Kxy <- kernelMatrix(kernel,x,y)
  
  ret <- kmmd(x=Kxx,y = Kyy,Kxy=Kxy, alpha=alpha, asymptotic= asymptotic, replace = replace, ntimes = ntimes, frac= frac)

  kernelf(ret) <- kernel
    
  return(ret)

})



setMethod("kmmd",signature(x="kernelMatrix"), function (x, y, Kxy, alpha = 0.05, asymptotic = FALSE,  replace = TRUE, ntimes = 100, frac = 1,  ...)
          {
            res <- new("kmmd")
            resmmd <- .submmd(x, y, Kxy, alpha) 
            H0(res) <- (resmmd$mmd1 > resmmd$D1) 
            Radbound(res) <- resmmd$D1
            Asymbound(res) <- 0
            mmdstats(res)[1] <- resmmd$mmd1
            mmdstats(res)[2] <- resmmd$mmd3
            
            if(asymptotic){
              boundA <- .submmd3bound(x, y, Kxy, alpha, frac, ntimes, replace)

              AsympH0(res) <- (resmmd$mmd1 > boundA) 
              Asymbound(res) <- boundA
            }
            kernelf(res) <- " Kernel matrix used as input."
            return(res)
           
          })

        
.submmd <- function(Kxx,Kyy, Kxy, alpha)
{

  m <- dim(Kxx)[1]
  n <- dim(Kyy)[1]

  N <- max(m,n)
  M <- min(m,n)

  sumKxx <- sum(Kxx)

  if(m!=n)
     sumKxxM <- sum(Kxx[1:M,1:M])
  else
     sumKxxM <- sumKxx

  dgxx <- diag(Kxx)

  sumKxxnd <- sumKxx - sum(dgxx)
  R <- max(dgxx)
  RM <- max(dgxx[1:M])
  hu <- colSums(Kxx[1:M,1:M]) - dgxx[1:M]

  sumKyy <- sum(Kyy)
  if(m!=n)
    sumKyyM <- sum(Kyy[1:M,1:M])
  else
    sumKyyM <- sumKyy

  dgyy <- diag(Kyy)

  sumKyynd <- sum(Kyy) - sum(dgyy)
  R <- max(R,dgyy)
  RM <- max(RM,dgyy[1:M]) # RM instead of R in original
  hu <- hu + colSums(Kyy[1:M,1:M]) - dgyy[1:M]
  
  sumKxy <- sum(Kxy)
  if (m!=n)
    sumKxyM <- sum(Kxy[1:M,1:M])
  else
    sumKxyM <- sumKxy
    
  dg <- diag(Kxy) # up to M only
  hu <- hu - colSums(Kxy[1:M,1:M]) - colSums(t(Kxy[1:M,1:M])) + 2*dg # one sided sum

  mmd1 <- sqrt(max(0,sumKxx/(m*m) + sumKyy/(n*n) - 2/m/n* sumKxy))
  mmd3 <- sum(hu)/M/(M-1)
  D1 <- 2*sqrt(RM/M)+sqrt(log(1/alpha)*4*RM/M)
  
  return(list(mmd1=mmd1,mmd3=mmd3,D1=D1))
}


.submmd3bound <- function(Kxx,Kyy, Kxy, alpha, frac, ntimes, replace)
  {
    ## implements the bootstrapping approach to the MMD3 bound by shuffling
    ## the kernel matrix
    ##  frac   : fraction of data used for bootstrap
    ##  ntimes : how many times MMD is to be evaluated

    m <- dim(Kxx)[1]
    n <- dim(Kyy)[1]
  
    M <- min(m,n)
    N <- max(m,n)
  
  poslabels <- 1:m
  neglabels <- (m+1):(m+n)
  
  ## bootstrap
  bootmmd3 <- rep(0,ntimes)
  
   for (i in 1:ntimes)
     {
       nsamples <- ceiling(frac*min(m,n))
       xinds <- sample(1:m,nsamples,replace=replace)
       yinds <- sample(1:n,nsamples,replace=replace)
       newlab <- c(poslabels[xinds],neglabels[yinds])
       samplenew <- sample(newlab, length(newlab), replace=FALSE)
       xinds <- samplenew[1:nsamples]
       yinds <- samplenew[(nsamples+1):length(samplenew)]
       
       newm <- length(xinds)
       newn <- length(yinds)
       newM <- min(newm,newn)
    
       ##get new kernel matrices (without concat to big matrix to save memory)
       xind1 <- xinds[xinds<=m]
       xind2 <- xinds[xinds>m]- m
       yind1 <- yinds[yinds<=m]
       yind2 <- yinds[yinds>m]-m
       
       ##Kxx (this should be implemented with kernelMult for memory efficiency)
       nKxx <- rbind(cbind(Kxx[xind1,xind1],Kxy[xind1,xind2]), cbind(t(Kxy[xind1,xind2]),Kyy[xind2,xind2]))
       dgxx <- diag(nKxx)
       hu <- colSums(nKxx[1:newM,1:newM]) - dgxx[1:newM]   # one sided sum
       rm(nKxx)

       #Kyy
       nKyy <- rbind(cbind(Kxx[yind1,yind1],Kxy[yind1,yind2]), cbind(t(Kxy[yind1,yind2]), Kyy[yind2,yind2]))
       dgyy <- diag(nKyy)
       hu <- hu + colSums(nKyy[1:newM,1:newM]) - dgyy[1:newM]
       rm(nKyy)

       ## Kxy
       nKxy <- rbind(cbind(Kxx[yind1,xind1],Kxy[yind1,xind2]), cbind(t(Kxy[xind1,yind2]),Kyy[yind2,xind2]))
       dg <- diag(nKxy)
       hu <- hu - colSums(nKxy[1:newM,1:newM]) - colSums(t(nKxy[1:newM,1:newM])) + 2*dg
       rm(nKxy)
    
       ## now calculate mmd3
       bootmmd3[i] <- sum(hu)/newM/(newM-1)
     }
  
    
    bootmmd3 <- sort(bootmmd3, decreasing=TRUE);
    aind <- floor(alpha*ntimes) ## better less than too much (-> floor);
    
  ## take threshold in between aind and the next smaller value:
    bound <- sum(bootmmd3[c(aind,aind+1)])/2;
    return(bound)

  }


setMethod("show","kmmd",
function(object){
 
  cat("Kernel Maximum Mean Discrepancy object of class \"kmmd\"","\n","\n")

    show(kernelf(object))

  if(is.logical(object@H0)){
    cat("\n")
    cat("\n","H0 Hypothesis rejected : ", paste(H0(object)))
    cat("\n","Rademacher bound : ", paste(Radbound(object)))
  }

  cat("\n") 

  if(Asymbound(object)!=0){
    cat("\n","H0 Hypothesis rejected (based on Asymptotic bound): ", paste(AsympH0(object)))
    cat("\n","Asymptotic bound : ", paste(Asymbound(object)))
  }

  cat("\n","1st and 3rd order MMD Statistics : ", paste( mmdstats(object)))
  cat("\n")
})
