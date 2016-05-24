## $Id: baseline.fillPeaks.R 170 2011-01-03 20:38:25Z bhm $
baseline.fillPeaks <- function(spectra, lambda, hwi, it, int){
  ## Iterative baseline correction algorithm based on mean suppression
  ## By Kristian Hovde Liland

  # INPUT:
  # spectra - rows of spectra
  # lambda  - 2nd derivative penalty of primary smoothing
  # hwi     - half width of local window
  # it      - number of iterations in suppression loop
  # int     - number of bucket intervals to make of data (or vektor of bucket boundaries)
  #
  # OUTPUT:
  # baseline  - estimated baseline
  # corrected - baseline corrected spectra
  
  # Initialization
  np <- dim(spectra)
  if (missing(int)) int <- np[1]-1
  baseline  <- matrix(0,np[1],np[2])
  
  # Sparse empty matrix (m x m)
  speye <- as.matrix.csr(0,np[2],np[2])
  
  # Diagonal sparse matrix (m x m)
  diag(speye) <- 1
  D <- diff(speye,differences=2)

  # ------==== S1: Smoothing ====------
  Yorig   <- spectra
  if(lambda > 0){
    U       <- chol(speye+10^lambda*t(D)%*%D)
    spectra <- t(backsolve(U, t(spectra)))
  }
  
  # Exponential decrease in interval width
  if(it != 1){
    d1 <- log10(hwi)
    d2 <- 0
    w <- ceiling((10)^c(d1+(0:(it-2))*(d2-d1)/(floor(it)-1), d2))
  } else {
    w <- hwi
  }

  # Compute bucket locations
  if(length(int)==1){
    lims   <- seq(from = 1, to = np[2], length = int + 1)
  } else {
    lims  <- int
    int <- length(int)-1
  }
  lefts  <- ceiling(lims[-(int+1)])
  rights <- floor(lims[-1])
  minip  <- round((lefts + rights)/2)
  
  # Iterate through spectra
  for(s in 1:np[1]){

    # ------==== S2: Subsampling ====------
    xx <- numeric(int)
    for (i in 1:int) xx[i] <- min(spectra[s,lefts[i]:rights[i]])
    
    # ------==== S3: Suppression ====------
    for(k in 1:it){
      # Current window width
      w0 <- w[k]
      
      # Point-wise iteration to the right
      for(i in 2:(int-1)){
        # Interval cut-off close to edges
        v <- min(c(i-1,w0,int-i))
        
        # Baseline suppression
        a <- mean(xx[(i-v):(i+v)])
        xx[i] <- min(a,xx[i])
      }
      
      # Point-wise iteration to the left
      for(i in 2:(int-1)){
        j <- int-i+1
        # Interval cut-off close to edges
        v <- min(c(i-1,w0,int-i))
        
        # Baseline suppression
        a <- mean(xx[(j-v):(j+v)])
        xx[j] <- min(a,xx[j])
      }
    }
    
    # Prepare minimum vector
    minip[1] <- 1
    minip[int] <- np[2]
    
    # ------==== S4: Stretch ====------
    xxx <- approx(minip, xx, 1:np[2])$y
    baseline[s,] <- xxx
  }
  list(baseline = baseline, corrected = Yorig - baseline)
}
