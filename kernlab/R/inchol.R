
setGeneric("inchol", function(x, kernel="rbfdot",kpar=list(sigma=0.1), tol= 0.001, maxiter = dim(x)[1], blocksize = 50, verbose = 0) standardGeneric("inchol"))
setMethod("inchol",signature(x="matrix"),
function(x, kernel="rbfdot",kpar=list(sigma=0.1), tol= 0.001, maxiter = dim(x)[1], blocksize = 50, verbose = 0)
{

  ## 
  ## Description:
  ## 
  ## Find the incomplete Cholesky decomposition of the kernel matrix
  ## 
  ## Parameters:
  ## 
  ## data      : 
  ## kernel    : kernlab object
  ## tol       : algo stops when remaining pivots < tol 
  ## max.iter  : maximum number of colums in Tk
  ##
  ## Return:
  ## 
  ## Tk     : K \approx Tk * Tk'
  ## pivots : Indices on which we pivoted
  ## diag.residues : Residuals left on the diagonal
  ## maxresiduals  : Residuals we picked for pivoting
  ## 
  ## Authors : S.V.N. Vishwanathan / Alex Smola
  ## R Version : Alexandros Karatzoglou
  

  ## For aggressive memory allocation
  BLOCKSIZE <- blocksize

if(!is.matrix(x))
	stop("x must be a matrix")

if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
        kernel <- do.call(kernel, kpar)
    }

     if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")
## Begin initialization

  Tk <- T <- pivots <- maxresiduals <- padded.veck <- matrix(0,0,0)
  counter <- 0

## End initialization

  m <- dim(x)[1] 

## Compute the diagonal of kernel matrix

  diag.residues <- matrix(0, m, 1)

  for (i  in  1:m) 
    diag.residues[i] <- kernel(x[i,],x[i,])
  
## Choose first pivot

  residue <- max(diag.residues)
  index <- which.max(diag.residues == residue)

  dota <- rowSums(x^2)

  while( residue > tol && counter < maxiter )
    {
      ## Aggressively allocate memory
      if(counter %% BLOCKSIZE  == 0)
        {
          Tktmp <- matrix(0, m, dim(Tk)[2] + BLOCKSIZE)
          Tktmp[1:m > 0, 1:(dim(Tk)[2] + BLOCKSIZE) <= dim(Tk)[2]] <- Tk
          Tk <- Tktmp

          Ttmp <- matrix(0, dim(T)[1]+BLOCKSIZE, BLOCKSIZE+counter)
          ind <- 1:(dim(T)[1]+BLOCKSIZE) <= dim(T)[1]
          ind2 <- 1:(BLOCKSIZE + counter) <= counter
          Ttmp[ind , ind2] <- T
          Ttmp[ind == FALSE, ind2 == FALSE] <- diag(1, BLOCKSIZE)
          T <- Ttmp
          
          padded.veck.tmp <- matrix(0,dim(padded.veck)[1]+BLOCKSIZE)
          padded.veck.tmp[1:(dim(padded.veck)[1]+BLOCKSIZE) <= dim(padded.veck)[1]] <- padded.veck
          padded.veck <-  padded.veck.tmp
    
          pivots.tmp <- matrix(0, dim(pivots)[1]+BLOCKSIZE)
          pivots.tmp[1:(dim(pivots)[1] + BLOCKSIZE)<= dim(pivots)[1]] <- pivots
          pivots <- pivots.tmp
     
          maxresiduals.tmp <- matrix(0,dim(maxresiduals)[1]+BLOCKSIZE)
          maxresiduals.tmp[1:(dim(maxresiduals)[1]+BLOCKSIZE) <= dim(maxresiduals)[1]] <- maxresiduals
          maxresiduals <- maxresiduals.tmp

          if(counter == 0)
            t <- rep(0,BLOCKSIZE)
          else
            t <- rep(0,length(t)+BLOCKSIZE)
        } 
    
      veck <- kernelFast(kernel, x, x[index, ,drop=FALSE],dota)
  
      if (counter == 0)
        {
          ## No need to compute t here
          tau <- sqrt(veck[index])
          
          ## Update T
          T[1, 1] <- tau
          ## Compute the update for Tk
          update <- veck/tau
        }
      else
        {
          padded.veck[1:counter] <- veck[pivots[1:counter]]

          ## First compute t
          ##          t <- t(crossprod(padded.veck,backsolve(T,diag(1,nrow=dim(T)[1]))))
          ## cat("T: ",dim(T), " p:",length(padded.veck),",\n")

          t[1:counter] <- backsolve(T, k=counter, padded.veck, transpose = TRUE)

          ## Now compute tau
          tau <- as.vector(sqrt(veck[index] - crossprod(t)))

          ## Update T

          T[1:counter, counter+1] <- t[1:counter]
          T[counter + 1, counter + 1] <- tau

          ## Compute the update for Tk
          update <- (1/tau) * (veck - Tk %*% t)
        }  

      ## Update Tk
      Tk[,counter + 1] <- update
      
      ## Update diagonal residuals
      diag.residues <- diag.residues - update^2			
      
      ## Update pivots
      pivots[counter + 1]  <- index

      
      ## Monitor residuals
      maxresiduals[counter + 1] <- residue 
      
      ## Choose next candidate
      residue <- max( diag.residues )
      index <- which.max(diag.residues)

      ## Update counter
      counter <- counter + 1
      
      ## Report progress to the user
      if(counter%%blocksize == 0 && (verbose == TRUE))
        cat("counter = ",counter," ", "residue = ", residue, "\n")
    } 


  ## Throw away extra columns which we might have added
  Tk <- Tk[, 1:counter] 
  
  pivots <- pivots[1:counter]

  maxresiduals <- maxresiduals[1:counter]

  return(new("inchol",.Data=Tk,pivots=pivots,diagresidues = diag.residues, maxresiduals = maxresiduals))
  
})
