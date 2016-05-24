#' Functional Principal Component Analysis
#' 
#' FPCA for dense or sparse functional data. 
#' 
#' @param y A list of \emph{n} vectors containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
#' @param t A list of \emph{n} vectors containing the observation time points for each individual corresponding to y.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#'
#' @details Available control options are 
#' \describe{
#' \item{userBwCov}{The bandwidth value for the smoothed covariance function; positive numeric - default: determine automatically based on 'methodBwCov'}
#' \item{methodBwCov}{The bandwidth choice method for the smoothed covariance function; 'GMeanAndGCV','CV','GCV' - default: 'GMeanAndGCV'')}
#' \item{userBwMu}{The bandwidth value for the smoothed mean function (using 'CV' or 'GCV'); positive numeric - default: determine automatically based on 'methodBwMu'}
#' \item{methodBwMu}{The bandwidth choice method for the mean function; 'GMeanAndGCV','CV','GCV' - default: 'GMeanAndGCV''} 
#' \item{dataType}{The type of design we have (usually distinguishing between sparse or dense functional data); 'Sparse', 'Dense', 'DenseWithMV', 'p>>n' - default:  determine automatically based on 'IsRegular'}
#' \item{diagnosticsPlot}{Make diagnostics plot (design plot, mean, scree plot and first k (<=3) eigenfunctions); logical - default: FALSE}
#' \item{error}{Assume measurement error in the dataset; logical - default: TRUE}
#' \item{fitEigenValues}{Whether also to obtain a regression fit of the eigenvalues - default: FALSE}
#' \item{FVEthreshold}{Fraction-of-Variance-Explained threshold used during the SVD of the fitted covar. function; numeric (0,1] - default: 0.9999}
#' \item{kernel}{Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "gauss"; dense data are assumed noise-less so no smoothing is performed. }
#' \item{kFoldMuCov}{The number of folds to be used for mean and covariance smoothing. Default: 10}
#' \item{lean}{If TRUE the 'inputData' field in the output list is empty. Default: FALSE}
#' \item{maxK}{The maximum number of principal components to consider; positive integer smaller than 128 - default: min(20, N-1), N:# of curves}
#' \item{methodXi}{The method to estimate the PC scores; 'CE' (Condit. Expectation), 'IN' (Numerical Integration) - default: 'CE' for sparse data, 'IN' for dense data.}
#' \item{methodMuCovEst}{The method to estimate the mean and covariance in the case of dense functional data; 'cross-sectional', 'smooth' - default: 'cross-sectional'}
#' \item{nRegGrid}{The number of support points in each direction of covariance surface; numeric - default: 51}
#' \item{numBins}{The number of bins to bin the data into; positive integer > 10, default: NULL}
#' \item{methodSelectK}{The method of choosing the number of principal components K; 'FVE','AIC','BIC', or a positive integer as specified number of components: default 'FVE')}
#' \item{shrink}{Whether to use shrinkage method to estimate the scores in the dense case (see Yao et al 2003) - default FALSE}
#' \item{outPercent}{A 2-element vector in [0,1] indicating the outPercent data in the boundary - default (0,1)}
#' \item{rho}{The truncation threshold for the iterative residual. 'cv': choose rho by leave-one-observation out cross-validation; 'no': no regularization - default "cv".}
#' \item{rotationCut}{The 2-element vector in [0,1] indicating the percent of data truncated during sigma^2 estimation; default  (0.25, 0.75))}
#' \item{useBinnedData}{Should the data be binned? 'FORCE' (Enforce the # of bins), 'AUTO' (Select the # of  bins automatically), 'OFF' (Do not bin) - default: 'AUTO'}
#' \item{useBinnedCov}{Whether to use the binned raw covariance for smoothing; logical - default:TRUE}
#' \item{userCov}{The user-defined smoothed covariance function; list of two elements: numerical vector 't' and matrix 'cov', 't' must cover the support defined by 'y' - default: NULL}
#' \item{userMu}{The user-defined smoothed mean function; list of two numerical vector 't' and 'mu' of equal size, 't' must cover the support defined 'y' - default: NULL}
#' \item{userSigma2}{The user-defined measurement error variance. A positive scaler. If specified then no regularization is used (rho is set to 'no', unless specified otherwise). Default to `NULL`}
#' \item{verbose}{Display diagnostic messages; logical - default: FALSE}
#' }
#' @return A list containing the following fields:
#' \item{sigma2}{Variance for measure error.}
#' \item{lambda}{A vector of length \emph{K} containing eigenvalues.}
#' \item{phi}{An nWorkGrid by \emph{K} matrix containing eigenfunctions, supported on workGrid.}
#' \item{xiEst}{A \emph{n} by \emph{K} matrix containing the FPC estimates.} 
#' \item{xiVar}{A list of length \emph{n}, each entry containing the variance estimates for the FPC estimates.}
#' \item{obsGrid}{The (sorted) grid points where all observation points are pooled.}
#' \item{mu}{A vector of length nObsGrid containing the mean function estimate.}
#' \item{workGrid}{A vector of length nWorkGrid. The internal regular grid on which the eigen analysis is carried on.}
#' \item{smoothedCov}{A nWorkGrid by nWorkGrid matrix of the smoothed covariance surface.}
#' \item{fittedCov}{A nWorkGrid by nWorkGrid matrix of the fitted covariance surface, which is guaranteed to be non-negative definite.}
#' \item{optns}{A list of actually used options.}
#' \item{bwMu}{The selected (or user specified) bandwidth for smoothing the mean function.}
#' \item{bwCov}{The selected (or user specified) bandwidth for smoothing the covariance function.}
#' \item{rho}{A regularizing scalar for the measurement error variance estimate.}
#' \item{cumFVE}{A vector with the percentages of the total variance explained by each FPC. Increase to almost 1.}
#' \item{FVE}{A percentage indicating the total variance explained by chosen FPCs with corresponding 'FVEthreshold'.}
#' \item{criterionValue}{A scalar specifying the criterion value obtained by the selected number of components with specific methodSelectK: FVE,AIC,BIC values or NULL for fixedK.}
#' \item{inputData}{A list containting the original 'y' and 't' lists used as inputs to FPCA. NULL if 'lean' was specified to be TRUE.}
#' 
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' CreateCovPlot(res, 'Fitted')
#' @references
#' \cite{Yao, F., Mueller, H.G., Clifford, A.J., Dueker, S.R., Follett, J., Lin, Y., Buchholz, B., Vogel, J.S. (2003). "Shrinkage estimation for functional principal component scores, with application to the population kinetics of plasma folate." Biometrics 59, 676-685. (Shrinkage estimates for dense data)}
#' 
#' \cite{Yao, Fang, Hans-Georg Mueller, and Jane-Ling Wang. "Functional data analysis for sparse longitudinal data." Journal of the American Statistical Association 100, no. 470 (2005): 577-590. (Sparse data FPCA)}
#'
#' \cite{Liu, Bitao, and Hans-Georg Mueller. "Estimating derivatives for samples of sparsely observed functions, with application to online auction dynamics." Journal of the American Statistical Association 104, no. 486 (2009): 704-717. (Sparse data FPCA)}
#'
#' \cite{Castro, P. E., W. H. Lawton, and E. A. Sylvestre. "Principal modes of variation for processes with continuous sample curves." Technometrics 28, no. 4 (1986): 329-337. (Dense data FPCA)}
#' @export

FPCA = function(y, t, optns = list()){
  
  # Check the data validity for further analysis
  CheckData(y,t)
  
  # Force the data to be list of numeric members and handle NA's
  #y <- lapply(y, as.numeric) 
  #t <- lapply(t, as.numeric)
  #t <- lapply(t, signif, 14)
  #inputData <- list(y=y, t=t);

  inputData <- HandleNumericsAndNAN(y,t);
  y <- inputData$y;
  t <- inputData$t;

  # Set the options structure members that are still NULL
  optns = SetOptions(y, t, optns);
  
  # Check the options validity for the PCA function. 
  numOfCurves = length(y);
  CheckOptions(t, optns,numOfCurves)

  # Bin the data
  if ( optns$useBinnedData != 'OFF'){ 
      BinnedDataset <- GetBinnedDataset(y,t,optns)
      y = BinnedDataset$newy;
      t = BinnedDataset$newt; 
  }

  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new
  # data
  # regGrid: the grid of time points for which the smoothed covariance
  # surface assumes values
  # cutRegGrid: truncated grid specified by optns$outPercent for the cov
  # functions
  obsGrid = sort(unique( c(unlist(t))));
  regGrid = seq(min(obsGrid), max(obsGrid),length.out = optns$nRegGrid);
  outPercent <- optns$outPercent
  buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
  rangeGrid <- range(regGrid)
  minGrid <- rangeGrid[1]
  maxGrid <- rangeGrid[2]
  cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * outPercent[1] -
                        buff & 
                        regGrid < minGrid + diff(rangeGrid) * outPercent[2] +
                        buff]

## Mean function
  # If the user provided a mean function use it
  userMu <- optns$userMu
  if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
    smcObj <- GetUserMeanCurve(optns, obsGrid, regGrid, buff)
    smcObj$muDense = ConvertSupport(obsGrid, regGrid, mu = smcObj$mu)
  } else if (optns$methodMuCovEst == 'smooth') { # smooth mean
    smcObj = GetSmoothedMeanCurve(y, t, obsGrid, regGrid, optns)
  } else if (optns$methodMuCovEst == 'cross-sectional') { # cross-sectional mean
    ymat = List2Mat(y,t)
    smcObj = GetMeanDense(ymat, obsGrid, optns)
  }
# mu: the smoothed mean curve evaluated at times 'obsGrid'
  mu <- smcObj$mu

## Covariance function and sigma2
  if (!is.null(optns$userCov) && optns$methodMuCovEst != 'smooth') { 
      scsObj <- GetUserCov(optns, obsGrid, cutRegGrid, buff, ymat)
  } else if (optns$methodMuCovEst == 'smooth') {
# smooth cov and/or sigma2
    scsObj = GetSmoothedCovarSurface(y, t, mu, obsGrid, regGrid, optns,
                                     optns$useBinnedCov) 
  } else if (optns$methodMuCovEst == 'cross-sectional') {
    scsObj = GetCovDense(ymat, mu, optns)
    scsObj$smoothCov = ConvertSupport(obsGrid, cutRegGrid, Cov =
                                      scsObj$smoothCov)
    scsObj$outGrid <- cutRegGrid
  }
  sigma2 <- scsObj[['sigma2']]
  # workGrid: possibly truncated version of the regGrid
  workGrid <- scsObj$outGrid

  
  # convert mu to truncated workGrid
  muWork <- ConvertSupport(obsGrid, toGrid = workGrid, mu=smcObj$mu)
  
  # Get the results for the eigen-analysis
  eigObj = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, workGrid, optns, muWork = muWork)

  # Truncated obsGrid, and observations. Empty observation due to truncation has length 0.
  truncObsGrid <- obsGrid
  if (!all(abs(optns$outPercent - c(0, 1)) < .Machine$double.eps * 2)) {
    truncObsGrid <- truncObsGrid[truncObsGrid >= min(workGrid) - buff &
                                 truncObsGrid <= max(workGrid) + buff]
    tmp <- TruncateObs(y, t, truncObsGrid)
    y <- tmp$y
    t <- tmp$t
  }

  # convert phi and fittedCov to obsGrid.
  muObs <- ConvertSupport(obsGrid, truncObsGrid, mu=mu)
  phiObs <- ConvertSupport(workGrid, truncObsGrid, phi=eigObj$phi)
  if (optns$methodXi == 'CE') {
    CovObs <- ConvertSupport(workGrid, truncObsGrid, Cov=eigObj$fittedCov)
  }

  # Get scores  
  if (optns$methodXi == 'CE') {
    if (optns$rho != 'no') { 
      if( length(y) > 2048 ){
        randIndx <- sample( length(y), 2048)
        rho <- GetRho(y[randIndx], t[randIndx], optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
      } else {
        rho <- GetRho(y, t, optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
      }
      sigma2 <- rho
    }
    scoresObj <- GetCEScores(y, t, optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
  } else if (optns$methodXi == 'IN') {
    ymat = List2Mat(y,t)
    scoresObj <- GetINScores(ymat, t, optns, muObs, eigObj$lambda, phiObs, sigma2)
  }

  if (optns$fitEigenValues) {
    fitLambda <- FitEigenValues(scsObj$rcov, workGrid, eigObj$phi, optns$maxK)
  } else {
    fitLambda <- NULL
  }

  # Make the return object by MakeResultFPCA
  ret <- MakeResultFPCA(optns, smcObj, muObs, scsObj, eigObj, 
                        inputData = inputData, 
                        scoresObj, truncObsGrid, workGrid, 
                        rho = if (optns$rho =='cv') rho else NULL, 
                        fitLambda=fitLambda)
  
  # select number of components based on specified criterion
  if(ret$optns$lean == TRUE){
    selectedK <- SelectK(fpcaObj = ret, criterion = optns$methodSelectK, FVEthreshold = optns$FVEthreshold,
                         y = y, t = t)
  } else {
    selectedK <- SelectK(fpcaObj = ret, criterion = optns$methodSelectK, FVEthreshold = optns$FVEthreshold)
  }
  
  ret <- append(ret, list(selectK = selectedK$k, criterionValue = selectedK$criterion))
  class(ret) <- 'FPCA'
  ret <- SubsetFPCA(fpcaObj = ret, k = ret$selectK)
  
  # Make a quick diagnostics plot     
  if(optns$diagnosticsPlot){
    CreateDiagnosticsPlot(ret);
  }

  return(ret); 
}

