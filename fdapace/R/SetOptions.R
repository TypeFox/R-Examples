#' Set the PCA option list
#'
#' @param y A list of \emph{n} vectors containing the observed values for each individual.
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' See '?FPCAfor more details. Usually users are not supposed to use this function directly.
#'


SetOptions = function(y, t, optns){

  methodMuCovEst = optns[['methodMuCovEst']]
  userBwMu =optns[['userBwMu']];                
  methodBwMu =optns[['methodBwMu']]; 
  userBwCov =optns[['userBwCov']];            
  methodBwCov =optns[['methodBwCov']];
  kFoldMuCov = optns[['kFoldMuCov']]
  methodSelectK =optns[['methodSelectK']];  
  FVEthreshold =optns[['FVEthreshold']];
  fitEigenValues <- optns[['fitEigenValues']];
  maxK =optns[['maxK']];                
  dataType =optns[['dataType']];          
  error =optns[['error']];
  nRegGrid =optns[['nRegGrid']];              
  methodXi =optns[['methodXi']];
  shrink =optns[['shrink']]
  kernel =optns[['kernel']];            
  numBins =optns[['numBins']];
  yname =optns[['yname']];
  rho =optns[['rho']];                
  diagnosticsPlot =optns[['diagnosticsPlot']];
  verbose =optns[['verbose']];   
  userMu =optns[['userMu']];                  
  #methodMu =optns[['methodMu']];
  outPercent =optns[['outPercent']];  
  userCov =optns[['userCov']];
  userSigma2 = optns[['userSigma2']]
  rotationCut =optns[['rotationCut']];    
  useBinnedData =optns[['useBinnedData']];
  useBinnedCov = optns[['useBinnedCov']]
  lean = optns[['lean']]

  if(is.null(userBwMu)){ # bandwidth choice for mean function is using CV or GCV
    userBwMu = 0;   
  }
  if(is.null(methodBwMu)){ # bandwidth choice for mean function is GCV if userBwMu = 0
    methodBwMu = 'GMeanAndGCV';  
  }
  if(is.null(userBwCov)){ # bandwidth choice for covariance function is CV or GCV
    userBwCov = 0; 
  }
  if(is.null(methodBwCov)){  # bandwidth choice for covariance function is GCV if userBwCov = c(0,0)
    methodBwCov = 'GMeanAndGCV';
  }
  #if(is.null(ngrid1)){ # number of support points for the covariance surface 
  #  ngrid1 = 30;
  #}
  if (is.null(kFoldMuCov)) { # CV fold for covariance smoothing
    kFoldMuCov <- 10L
  } else {
    kFoldMuCov <- as.integer(kFoldMuCov)
  }
  if(is.null(methodSelectK)){ # the method of choosing the number of principal components K
    #  TODO : Possibly have user-defined selection methods for the # of FPCs and we keep
    # an internal FVE-based method for us
    methodSelectK = "FVE";
  }
  if(is.null(FVEthreshold)){  # Default Value for the Fraction-of-Variance-Explained
     FVEthreshold = 0.9999;
  }
  if(is.null(maxK)){ # maximum number of principal components to consider
    maxK = min(20, length(y)-1);   
  }
  if(is.null(dataType)){ #do we have dataType or sparse functional data
    dataType = IsRegular(t);    
  }
  if (is.null(fitEigenValues)) {
    fitEigenValues <- FALSE
  }
  if (fitEigenValues && dataType == 'Dense') {
    stop('Fit method only apply to sparse data')
  }
  if(is.null(error)){ # error assumption with measurement error
      error = TRUE;    
  }
  if(is.null(nRegGrid)){ # number of support points in each direction of covariance surface 
    if(dataType == 'Dense' || dataType == 'DenseWithMV'){
      tt = unlist(t)
      nRegGrid = length(unique(signif(tt[!is.na(tt)],6)));
    } else { # for Sparse and p>>n
      nRegGrid = 51;
    }    
  }
  methodNames = c("IN", "CE");
  if(!is.null(methodXi) && !(methodXi %in% methodNames)){
    cat(paste('methodXi', methodXi, 'is unrecognizable! Reset to automatic selection now!\n')); 
    methodXi = NULL; 
  }   
  if(is.null(methodXi)){ # method to estimate the PC scores
    if(dataType == 'Dense'){
      methodXi = "IN";
    } else if(dataType == 'Sparse'){
      methodXi = "CE";
    } else if(dataType == 'DenseWithMV'){
      methodXi = "IN";
    } else { # for dataType = p>>n
      methodXi = "IN";
    }
  }
   if(is.null(shrink)){ 
     # apply shrinkage to estimates of random coefficients (dataType data
     # only)
     shrink = FALSE;
   }
   if(shrink == TRUE && (error != TRUE || methodXi != "IN")){ 
     # Check for valid shrinkage choice
     cat('shrinkage method only has effects when methodXi = "IN" and error = TRUE! Reset to shrink = FALSE now!\n');
     shrink = FALSE      
   }
  if(is.null(kernel)){ # smoothing kernel choice
    if(dataType == "Dense"){
      kernel = "epan";   # kernel: Epanechnikov
    }else{
      kernel = "gauss";  # kernel: Gaussian
    }
  }
  kernNames = c("rect", "gauss", "epan", "gausvar", "quar");
  if(!(kernel %in% kernNames)){ # Check suitability of kernel
    cat(paste('kernel', kernel, 'is unrecognizable! Reset to automatic selection now!\n')); 
    kernel = NULL; 
  }  
  if(is.null(kernel)){ # smoothing kernel choice
    if(dataType %in% c( "Dense", "DenseWithMV")){
      kernel = "epan";   # kernel: Epanechnikov
    }else{
      kernel = "gauss";  # kernel: Gaussian
    }
  }
  if(is.null(yname)){ # name of the variable analysed
    yname = as.character(substitute(y))      
  }
  if(maxK > (nRegGrid-2)){ # check if a reasonable number of eigenfunctions is requested
    cat(paste("maxK can only be less than or equal to", nRegGrid-2,"! Reset to be", nRegGrid-2, "now!\n"));
    maxK = nRegGrid -2;
  }
  if(is.numeric(methodSelectK)){
    FVEthreshold <- 1 # disable FVE selection.
    if(methodSelectK > (nRegGrid-2)){ # check if a reasonable number of eigenfunctions is requested
      cat(paste("maxK can only be less than or equal to", nRegGrid-2,"! Reset to be", nRegGrid-2, "now!\n"));
      maxK = nRegGrid -2;
    }else if(methodSelectK <= 0){ # check if a positive number of eigenfunctions is requested
      cat("methodSelectK must be a positive integer! Reset to BIC now!\n");
      methodSelectK = "BIC"
      FVEthreshold = 0.95;
    }
  }
  if(is.null(diagnosticsPlot)){ # make corrplot
    diagnosticsPlot = FALSE;
  }
  if(is.null(rho)){ # truncation threshold for the iterative residual that is used
    if (!is.null(userSigma2)) { # no regularization if sigma2 is specified
      rho <- 'no'
    } else {
      rho <- 'cv'
    }
  }
  if(is.null(verbose)){ # display diagnostic messages
    verbose = FALSE;
  }  
  if(is.null(userMu)){ # user-defined mean functions valued at distinct input time points
    userMu <- NULL
  }
  if(is.null(userCov)){
    userCov <- NULL
  }
  #if(is.null(methodMu)){ # method to estimate mu
  #  methodMu <- 'PACE'
  #}
  if(is.null(outPercent)){ 
    outPercent <- c(0,1)
  }  
  if(is.null(rotationCut)){ 
    rotationCut <- c(0.25,.75)
  } 
  # if(error == FALSE && (methodSelectK == "AIC" || methodSelectK == "BIC")){ # Check suitability of information criterion
  #  cat('When assume no measurement error, cannot use "AIC" or "BIC". Reset to "BIC" now!\n')
  #  methodSelectK = "BIC" 
  #}
  if(!is.null(numBins)){ 
    if(numBins < 10){   # Check suitability of number of bins
      cat("Number of bins must be at least +10!!\n");
      numBins = NULL;
    }
  }
  if(is.null(useBinnedData)){ 
    useBinnedData = 'AUTO';
  }
  if (is.null(useBinnedCov)) {
    useBinnedCov <- TRUE
    if (  ( 128 > length(y) ) && ( 3 > mean ( unlist( lapply( y, length) ) ) )){
      useBinnedCov <- FALSE
    } 
  }
  if(is.null(lean)){ 
    lean = FALSE;
  }
  if(is.null(methodMuCovEst)){
    if (dataType == 'Sparse'){
      methodMuCovEst = 'smooth'; #In the odd case that somehow we use this...
    } else {
      methodMuCovEst = 'cross-sectional';
    }
  }
  # if (!all.equal(outPercent, c(0, 1)) && methodMuCovEst == 'cross-sectional') {
    # stop('outPercent not supported for cross-sectional covariance estimate')
  # }
    
  retOptns <- list(userBwMu = userBwMu, methodBwMu = methodBwMu, userBwCov = userBwCov, methodBwCov = methodBwCov,
          kFoldMuCov = kFoldMuCov, methodSelectK = methodSelectK, FVEthreshold = FVEthreshold,
          fitEigenValues = fitEigenValues, maxK = maxK, dataType = dataType, error = error, shrink = shrink,
          nRegGrid = nRegGrid, rotationCut = rotationCut, methodXi = methodXi, kernel = kernel, 
          lean = lean, diagnosticsPlot = diagnosticsPlot, numBins = numBins, useBinnedCov = useBinnedCov, 
          yname = yname,  rho = rho, verbose = verbose, userMu = userMu, userCov = userCov, methodMuCovEst = methodMuCovEst,
          userSigma2 = userSigma2, outPercent = outPercent, useBinnedData = useBinnedData)

  invalidNames <- !names(optns) %in% names(retOptns)
  if (any(invalidNames)) {
    stop(sprintf('Invalid option names: %s',
                 paste0(names(optns)[invalidNames], collapse=', ')))
  }
  return( retOptns )
}
