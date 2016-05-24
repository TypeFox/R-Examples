#' Check option format
#'
#' Check if the options structure is valid and set the ones that are NULL
#' 
#' @param t is a n-by-1 list of vectors 
#' @param optns is an initialized option list
#' @param n is a total number of sample curves

CheckOptions = function(t,optns,n){
  
    
  if( !(  any(optns[['useBinnedData']] == c('FORCE','AUTO','OFF')) )){ 
    # Force, turn off or automatically decide about the use of bin data
    stop("FPCA is aborted because the argument: useBinnedData is invalid!\n"); 
  }
  if(  !( (length(optns[['userBwMu']])==1) &&  is.numeric(optns[['userBwMu']]) && (0<=optns[['userBwMu']]) ) ){ 
    # bandwidth Bhoice for mean function is using CV or GCV
    stop("FPCA is aborted because the argument: userBwMu is invalid!\n"); 
  }
  if( !(  any(optns[['methodBwMu']] == c('CV','GCV','GMeanAndGCV')) )){ 
    # bandwidth choice for mean function is GCV if userBwMu = 0
    stop("FPCA is aborted because the argument: methodBwMu is invalid!\n"); 
  }
  if(!(length(optns[['userBwCov']])==1) &&  is.numeric(optns[['userBwCov']]) && (all(optns[['userBwCov']]>=0))){ 
    # bandwidth choice for covariance function is CV or GCV
    stop("FPCA is aborted because the argument: userBwCov is invalid!\n"); 
  }
  if( !(  any(optns[['methodBwCov']] == c('CV','GCV','GMeanAndGCV') ) )){ 
    # bandwidth choice for covariance function is GCV if userBwCov = c(0,0)
    stop("FPCA is aborted because the argument: methodBwCov is invalid!\n");  
  }

  if (is.nan(optns[['kFoldMuCov']]) || optns[['kFoldMuCov']] < 2) {
    stop('Invalid `kFoldMuCov` option')
  }

  if( !(any(optns[['methodSelectK']] == c('FVE','AIC','BIC')))){
    if(is.numeric(optns[['methodSelectK']])){
      if(as.integer(optns[['methodSelectK']]) != optns[['methodSelectK']] || optns[['methodSelectK']] <= 0 || optns[['methodSelectK']] > n){
        stop('FPCA is aborted because the argument: methodSelectK is invalid!\n')
      }
    } else {
      stop('FPCA is aborted because the argument: methodSelectK is invalid!\n')
    }
  }
  if(  ( (length(optns[['FVEthreshold']])==1) &&  is.numeric(optns[['FVEthreshold']]) ) ){
    if (!( (0<=optns[['FVEthreshold']]) && (optns[['FVEthreshold']]<=1) ) ){  
      # the Fraction-of-Variance-Explained
      stop("FPCA is aborted because the argument: FVEthreshold is invalid!\n"); 
    } 
  }
  if( !( (length(optns[['maxK']])==1) &&  is.numeric(optns[['maxK']]) && (1<=optns[['maxK']]) && (optns[['maxK']]<=n) )){  
    # maximum number of principal components to consider
    stop("FPCA is aborted because the argument: maxK is invalid!\n");   
  } 
  #if( !is.null(optns$numComponents) ) {
  #  if( !( (length(optns$numComponents)==1) &&  is.numeric(optns$numComponents) && (1<=optns$numComponents) && (optns$numComponents<=n) )){  
  #    # maximum number of principal components to return
  #    stop("FPCA is aborted because the argument: numComponents is invalid!\n");   
  #  }
  #} 
  if( !( is.null(optns[['dataType']]) || any(optns[['dataType']]==c("Sparse","DenseWithMV","Dense","p>>n")) )){ 
    #do we have regualr or sparse functional data
    stop("FPCA is aborted because the argument: dataType is invalid!\n");     
  }   
  if( ( is.null(optns[['dataType']])  )){ 
    optns[['dataType']] = IsRegular(t)
  }   
  if(!is.logical(optns[['error']])){ 
    # error assumption with measurement error 
    stop("FPCA is aborted because the error option is invalid!\n");   
  }
  if( !( (length(optns[['nRegGrid']])==1) &&  is.numeric(optns[['nRegGrid']]) && (1<=optns[['nRegGrid']]) && (optns[['nRegGrid']]>optns[['maxK']]) ) ){
    # number of support points in each direction of covariance surface  
    stop("FPCA is aborted because the argument: nRegGrid is invalid!\n");    
  }
  if( !(any(optns[['methodXi']] == c('CE','IN')))){ 
    #method to estimate the PC scores
    stop("FPCA is aborted because the argument: methodXi is invalid!\n");   
  }
  if (optns[['methodXi']] == 'IN' && optns[['dataType']] != 'Dense') {
    stop("integration method can only be applied on dense data now!")
  }
  if(!(any(optns[['kernel']] == c('epan','gauss','rect','quar','gausvar')))){ 
    #method to estimate the PC scores
    stop("FPCA is aborted because the argument: kernel is invalid!\n");   
  }
  if( !( ( is.numeric(optns[['numBins']]) && (optns[['numBins']]>1)) || is.null(optns[['numBins']]) )  ){  
    # Check suitability of number of bins
    stop("FPCA is aborted because the argument: numBins is invalid!\n");   
  }
  if( ( ( optns[['useBinnedData']] == 'FORCE') &&  is.null(optns[['numBins']]) ) ){  
    # Check that we have a number of the bins if we force binning
    stop("FPCA is aborted because the argument: numBins is NULL but you FORCE binning!\n");   
  }
  if(!is.character(optns[['yname']])){ 
    # name of the variable analysed     
    stop("FPCA is aborted because the argument: yname is invalid!\n");  
  }
  if(!is.logical(optns[['diagnosticsPlot']])){ 
    # make diagnosticsPlot 
    stop("FPCA is aborted because the argument: diagnosticsPlot is invalid!\n");    
  }
  if(!(any(optns[['rho']] == c('cv-random', 'cv', 'none', 'no')))){ 
    # truncation threshold for the iterative residual that is used 
    stop("FPCA is aborted because the argument: rho is invalid!\n");     
  }
  if(!is.logical(optns[['verbose']])){ 
    # display diagnostic messages
    stop("FPCA is aborted because the argument: verbose is invalid!\n");    
  }
  if (! (  is.null(optns[['userMu']]) || 
        (is.list(optns[['userMu']]) && is.vector(optns[['userMu']][['t']]) &&  is.vector(optns[['userMu']][['mu']]) &&
         ( length(optns[['userMu']][['t']]) ==  length(optns[['userMu']][['mu']]) ) ))){      
    # display diagnostic messages
    stop("FPCA is aborted because the argument: userMu is invalid!\n");     
  }
  if (! ( is.null(optns[['userCov']]) ||
        ( is.list(optns[['userCov']]) && is.vector(optns[['userCov']][['t']]) &&  is.matrix(optns[['userCov']][['cov']]) &&
          (length(optns[['userCov']][['t']]) ==  ncol(optns[['userCov']][['cov']]) ) && ( isSymmetric(optns[['userCov']][['cov']]) ) ) ) ){
    # display diagnostic messages
    stop("FPCA is aborted because the argument: userCov is invalid! (eg. Check if 'cov' is symmetric and 't' is of appropriate size.)\n");
    return(TRUE);
  }
  if (!is.null(optns[['userSigma2']])) {
    if (!(is.numeric(optns[['userSigma2']]) && 
          length(optns[['userSigma2']]) == 1 && 
          optns[['userSigma2']] >= 0)) {
      stop('userSigma2 invalid.')
    }
    
    if (optns[['userSigma2']] == 0 && optns[['error']]) {
      stop('userSigma2 specified to be 0 but error = TRUE. If no measurement error is assumed then use error = FALSE.')
    }
  }
  #if(!(any(optns[['methodMu']] == c('PACE','RARE','CrossSectional')))){ 
  #  # user-defined mean functions
  #  stop("FPCA is aborted because the argument: methodMu is invalid!\n");     
  #}
  if( !( (length(optns[['outPercent']])==2) &&  is.numeric(optns[['outPercent']]) && all(0<=optns[['outPercent']]) && all(optns[['outPercent']]<=1) )){ 
    # display diagnostic messages
    stop("FPCA is aborted because the argument: outPercent is invalid!\n");    
  }
  if( !( (length(optns[['rotationCut']])==2) &&  is.numeric(optns[['rotationCut']]) && all(0<=optns[['rotationCut']]) && all(optns[['rotationCut']]<=1) )){ 
    # display diagnostic messages
    stop("FPCA is aborted because the argument: rotationCut is invalid!\n");    
  }
  if(is.logical(optns[['userCov']])){ 
    # display diagnostic messages
    stop("FPCA is aborted because the argument: userCov is invalid!\n");     
  }
  if(!(any(optns[['methodMuCovEst']] == c('smooth', 'cross-sectional')))){
    stop("FPCA is aborted because the argument: methodMuCovEst is invalid!\n");    
  }

}

