#' Take derivative of an FPCA object
#' 
#' @param fpcaObj A object of class FPCA returned by the function FPCA().   
#' @param derOptns A list of options to control the derivation parameters specified by \code{list(name=value)}. See `Details'. (default = NULL)
#'
#' @details Available derivation control options are 
#' \describe{
#' \item{p}{The order of the derivatives returned (default: 0, max: 2)}
#' \item{bw}{Bandwidth for smoothing the derivatives (default: p * 0.1 * S)}
#' \item{kernelType}{Smoothing kernel choice; same available types are FPCA(). default('epan')}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' derRes <- FPCAder(res)
#' @export


FPCAder <-  function (fpcaObj, derOptns = list(p=1)) {

  derOptns <- SetDerOptions(fpcaObj,derOptns = derOptns)
  p <- derOptns[['p']]
  method <- derOptns[['method']]
  bw <- derOptns[['bw']]
  kernelType <- derOptns[['kernelType']]
  k <- derOptns[['k']]

  obsGrid <- fpcaObj$obsGrid
  workGrid <- fpcaObj$workGrid
  
  if (!class(fpcaObj) %in% 'FPCA'){
    stop("FPCAder() requires an FPCA class object as basic input")
  }

  if( ! (p %in% c(1, 2))){
    stop("'FPCAder()' is requested to use a derivative order other than p = {1, 2}!")
  } 
  
  if (p == 2) {
    warning('Second derivative is experimental only.')
  } 

  #if( k > SelectK( fpcaObj, FVEthreshold=0.95, criterion='FVE')$k ){
  #    warning("Potentially you use too many components to estimate derivatives. \n  Consider using SelectK() to find a more informed estimate for 'k'.");
  #}

   #fittedCurves <- fitted(fpcaObj, k = k, derOpts = NULL)
   #A = MakeFPCAInputs(yVec=fittedCurves, tVec= rep(fpcaObj$workGrid))
   #fcfpcaObj <- FPCA(A$Ly, t= A$Lt, list(lean=TRUE))
   #xin = unlist(fpcaObj$input$t );    
   #yin = unlist(fpcaObj$input$y)[order(xin)];
   #xin = sort(xin);    
   #win = rep(1, length(xin));
   
   # muDer  = Lwls1D(bw, kernelType, npoly = p+1, nder = p, xin = xin, yin= yin, xout = obsGrid, win = win) 
   
   muDer <- Lwls1D(bw, kernelType, rep(1, length(workGrid)), workGrid, fpcaObj$mu, workGrid, p+0, nder= p)
   phiDer <- apply(fpcaObj$phi, 2, function(phi) Lwls1D(bw, kernelType, rep(1, length(workGrid)), workGrid, phi, workGrid, p+0, nder= p))

 # muDer2<- fpcaObj$mu
 # phiDer2 <- fpcaObj$phi
 # for (i in seq_len(p)) {
 #    # derivative
 #    muDer2 <- getDerivative(y = muDer2, t = obsGrid, ord=p)
 #    phiDer2 <- apply(phiDer2, 2, getDerivative, t= workGrid, ord=p)
 #    # smooth
 #    muDer2 <- getSmoothCurve(t=obsGrid, 
 #                            ft= muDer2,
 #                            GCV = TRUE,
 #                            kernelType = kernelType, mult=2)
 #    phiDer2 <- apply(phiDer2, 2, function(x)
 #                      getSmoothCurve(t=workGrid, ft=x, GCV =TRUE, kernelType = kernelType, mult=1))
 # }

  fpcaObj <- append(fpcaObj, list(muDer = muDer, phiDer = phiDer, derOptns = derOptns))
  class(fpcaObj) <- c(class(fpcaObj), 'FPCAder')
  return(fpcaObj)
}
