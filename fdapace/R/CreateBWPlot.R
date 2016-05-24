#' Functional Principal Component Analysis Bandwidth Diagnostics plot
#' 
#' This function by default creates the mean and first principal modes of variation plots for
#' 50%, 75%, 100%, 125% and 150% of the defined bandwidth choices in the fpcaObj provided as input.
#' If provided with a derivative options object (?FPCAder) it will return the 
#' differentiated mean and first two principal modes of variations for 50%, 75%, 100%, 125% and 150% of the defined bandwidth choice.
#'
#' @param fpcaObj An FPCA class object returned by FPCA().
#' @param derOptns A list of options to control the derivation parameters; see ?FPCAder. If NULL standard diagnostics are returned 
#' @param bwMultipliers A vector of multipliers that the original 'bwMu' and 'bwCov' will be multiplied by. (default: c(0.50, 0.75, 1.00, 1.25, 1.50))
#' - default: NULL
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res1 <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=FALSE))
#' CreateBWPlot(res1)
#' @export

CreateBWPlot <-function(fpcaObj, derOptns = NULL, bwMultipliers = NULL){ 
  
  if(class(fpcaObj) != 'FPCA'){
    stop("Input class is incorrect; CreateDiagnosticsPlot() is only usable from FPCA objects.")
  }
  oldPar <- par()
  if( is.null(bwMultipliers)){
    bwMultipliers = c(0.50, 0.75, 1.00, 1.25, 1.50)
    # This is knowingly wasteful as 1.00 is already computed; not having it would perplex the code 
    # unnecessarily though.
  }
 
  M <- length(bwMultipliers)
  
  if(is.null(derOptns) || !is.list(derOptns)){
  
    if(fpcaObj$optns$lean){
      stop("FPCA bandwidth diagnostics are not available for lean FPCA objects.")
    }

    newFPCA <- function(mlt){
      optnsNew = fpcaObj$optns; 
      optnsNew[c('userBwMu', 'userBwCov')] = mlt * unlist(fpcaObj[c('bwMu', 'bwCov')])
      return( FPCA(y= fpcaObj$input$y, t= fpcaObj$input$t, optnsNew) )
    }

    yy = lapply( bwMultipliers, function(x)  tryCatch( newFPCA(x), error = function(err) {                                                
                                                     #warning('Probable invalid bandwidth. Try enlarging the window size.')
                                                     stop( paste( collapse =' ', c('Multiplier :', x, 'failed to return valid FPCA object. Change multiplier values.'))) 
                                                      return(NA)}))
   # if( any(is.na(yy))){
   #   warning( paste( collapse =' ', c('Multipliers :', bwMultipliers[is.na(yy)], 'fail to return valid FPCA objects.')))
   #   bwMultipliers = bwMultipliers[!is.na(yy)]  
   #   M = length(bwMultipliers)
   #   yy[[is.na(yy)]] = NULL
   # }

    par(mfrow=c(1,3))

    Z = rbind( sapply(1:M, function(x) yy[[x]]$mu));
    matplot(x = fpcaObj$workGrid, lty= 1, type='l',  Z, ylab= expression(paste(collapse = '', mu, "(s)")), xlab = 's')
    grid(); legend('topright', lty = 1, col=1:M, legend = apply( rbind( rep('bwMu: ',M), round( digits = 2, bwMultipliers * fpcaObj$bwMu)), 2, paste, collapse = ''))

    Z = rbind(sapply(1:M, function(x) yy[[x]]$phi[,1]));
    matplot(x = fpcaObj$workGrid, lty= 1, type='l',  Z,   ylab= expression(paste(collapse = '', phi[1], "(s)")),  xlab = 's')
    grid(); legend('topright', lty = 1, col=1:M, legend = apply( rbind( rep('bwCov: ',M), round( digits = 2, bwMultipliers * fpcaObj$bwCov)), 2, paste, collapse = ''))
 
    Z = rbind(sapply(1:M, function(x) yy[[x]]$phi[,2]));
    matplot(x = fpcaObj$workGrid, lty= 1, type='l',  Z,   ylab= expression(paste(collapse = '', phi[2], "(s)")),  xlab = 's')
    grid(); legend('topright', lty = 1, col=1:M, legend = apply( rbind( rep('bwCov: ',M), round( digits= 2, bwMultipliers * fpcaObj$bwCov)), 2, paste, collapse = ''))
    

    } else {
    
    derOptns <- SetDerOptions(fpcaObj,derOptns = derOptns) 
    p <- derOptns[['p']]
    method <- derOptns[['method']]
    bw <- derOptns[['bw']]
    kernelType <- derOptns[['kernelType']]
    k <- derOptns[['k']]
    if(p==0){
      stop("Derivative diagnostics are inapplicable when p = 0")
    }
    
    yy = lapply( bwMultipliers *  bw, function(x) FPCAder(fpcaObj, list(p=p, method = method, kernelType = kernelType, k = k, bw = x)))
    
    par(mfrow=c(1,3))
    
    Z = rbind(sapply(1:M, function(x) yy[[x]]$muDer));
    matplot(x = fpcaObj$workGrid, lty= 1, type='l',  Z, ylab= expression(paste(collapse = '', 'd', mu, "/ds")), 
            main= substitute(paste("Derivatives of order ", p, " of ", mu)), xlab = 's')
    grid(); legend('topright', lty = 1, col=1:M, legend = apply( rbind( rep('bw: ',M), round( digits=2, bwMultipliers * bw)), 2, paste, collapse = ''))
    
    Z = rbind(sapply(1:M, function(x) yy[[x]]$phiDer[,1]));
    matplot(x = fpcaObj$workGrid, lty= 1, type='l',  Z,   ylab= expression(paste(collapse = '', 'd', phi[1], "/ds")), 
            main= substitute(paste("Derivatives of order ", p, " of ", phi[1])), xlab = 's')
    grid(); legend('topright', lty = 1, col=1:M, legend = apply( rbind( rep('bw: ',M), round(digits= 2, bwMultipliers * bw)), 2, paste, collapse = ''))
    
    Z = rbind(sapply(1:M, function(x) yy[[x]]$phiDer[,2]));
    matplot(x = fpcaObj$workGrid, lty= 1, type='l',  Z, ylab= expression(paste(collapse = '', 'd', phi[2], "/ds")), 
            main= substitute(paste("Derivatives of order ", p, " of ", phi[2])), xlab = 's')
    grid(); legend('topleft', lty = 1, col=1:M, legend = apply( rbind( rep('bw: ',M),  round( bwMultipliers * bw, digits=2) ), 2, paste, collapse = ''))
    
  }
  suppressWarnings(par(oldPar))
}

