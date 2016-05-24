#
# Various post-estimation functions
#

#'@title 	Total, Direct, and Indirect Effects for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the standard generic function \code{effects} computes total, direct, 
#'and indirect effects for a matrixpls results according to the method described in Fox (1980).
#'
#'Adapted from the \code{\link[sem]{effects}} function of the \code{sem} package
#'
#'@template postestimationFunctions
#'
#'@return A list with \code{Total}, \code{Direct}, and \code{Indirect} elements.
#'
#'@references
#'
#'Fox, J. (1980) Effect analysis in structural equation models: Extensions and simplified methods of computation. \emph{Sociological Methods and Research}
#'9, 3--28.
#'
#'@importFrom stats effects
#'
#'@method effects matrixpls
#'
#'@S3method effects matrixpls

effects.matrixpls <- function(object,  ...) {
  
  A <- attr(object,"inner")
  endog <- rowSums(attr(object,"model")$inner)!=0 
  
  I <- diag(endog)
  AA <- - A
  diag(AA) <- 1
  Total <- solve(AA) - I
  Indirect <-  Total - A
  result <- list(Total=Total[endog, , drop = FALSE], Direct=A[endog, , drop = FALSE], Indirect=Indirect[endog, , drop = FALSE])
  class(result) <- "matrixplseffects"
  result
}

#'@S3method print matrixplseffects

print.matrixplseffects <- function(x, ...){
  cat("\n Total Effects (column on row)\n")
  Total <- x$Total
  Direct <- x$Direct
  Indirect <- x$Indirect
  select <- !(apply(Total, 2, function(x) all( x == 0)) & 
                apply(Direct, 2, function(x) all( x == 0)) & 
                apply(Indirect, 2, function(x) all( x == 0)))
  print(Total[, select, drop = FALSE], ...)
  cat("\n Direct Effects\n")
  print(Direct[, select, drop = FALSE], ...)
  cat("\n Indirect Effects\n")
  print(Indirect[, select, drop = FALSE], ...)
  invisible(x)
}

#'@title Residual diagnostics for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for generic function \code{residuals} computes the residual
#'covariance matrix and various fit indices. 
#'
#'@details The residuals can be
#'either observed empirical residuals as presented by Lohmöller (1989, ch 2.4) or 
#'model implied residuals as used by Henseler et al. (2014). Two versions of the SRMR index are
#'provided, the traditional SRMR that includes all residual covariances, and the version
#'proposed by Henseler et al. (2014) where the within-block residual covariances are ignored.
#'
#'@template postestimationFunctions
#'
#'@param observed If \code{TRUE} (default) the observed residuals from the outer model regressions
#'(indicators regressed on composites) are returned. If \code{FALSE}, the residuals are calculated
#'by combining  \code{inner}, \code{reflective}, and \code{formative} as a simultaneous equations
#'system and subtracting the covariances implied by this system from the observed covariances.
#'The error terms are constrained to be uncorrelated and covariances between exogenous observed
#'values are fixed at their sample values.
#'
#'@return A list with three elements: \code{inner}, \code{outer}, and \code{indices} elements
#' containing the residual covariance matrix of regressions of composites on other composites,
#' the residual covariance matrix of indicators on composites, and various indices
#' calculated based on the residuals.
#'
#'@references
#'
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial
#'least squares.} Heidelberg: Physica-Verlag.
#'
#'Henseler, J., Dijkstra, T. K., Sarstedt, M., Ringle, C. M., Diamantopoulos, A., Straub, D. W., …
#'Calantone, R. J. (2014). Common Beliefs and Reality About PLS Comments on Rönkkö and Evermann 
#'(2013). Organizational Research Methods, 17(2), 182–209. doi:10.1177/1094428114526928
#'
#'@export
#'
#'@method residuals matrixpls
#'
#'@S3method residuals matrixpls

residuals.matrixpls <- function(object, ..., observed = TRUE) {
  
  RMS <- function(num) sqrt(sum(num^2)/length(num))
  
  
  
  # Equation numbers in parenthesis refer to equation number in Lohmoller 1989
  
  if(observed){
    nativeModel <- attr(object,"model")
    
    W <- attr(object,"W")
    S <- attr(object,"S")
    
    # Number of reflective indicators
    reflectiveIndicators<- rowSums(nativeModel$reflective)>0
    k <- sum(reflectiveIndicators)
    
    # Number of endog LVs
    endog <- rowSums(nativeModel$inner)>0
    h <- sum(endog)	
    
    
    
    # Factor loading matrix
    P <- nativeModel$reflective
    
    P[P==1] <- object[grepl("=~", names(object), fixed=TRUE)]
    
    B <- attr(object,"inner")
    
    # Lohmoller is not clear whether R should be based on the estimated betas or calculated scores
    # The scores are used here because this results in less complex code
    
    R <- attr(object,"C")
    R_star <- (B %*% R %*% t(B))[endog,endog] # e. 2.99
    
    
    # Model implied indicator covariances
    H <- P %*% R %*% t(P) # eq 2.96
    
    I <- diag(S)
    H2 <- (I * H)  %*% MASS::ginv(I * S) # eq 2.9	
    
    F <- P %*% B %*% R %*% t(B) %*% t(P) # eq 2.104
    
    F2 <- (I * F) %*% MASS::ginv(I * S) # eq 2.105
    
    r2 <- r2(object)
    
    
    # C in Lohmoller 1989 is different from the matrixpls C
    # residual covariance matrix of indicators
    
    C <- S-H
    
    Q <- (W %*% S %*% t(W))[endog,endog] - R_star
    
    
    indices <- list(Communality = psych::tr(H2)/k, # eq 2.109
                    Redundancy = psych::tr(F2)/k,  # eq 2.110
                    SMC = sum(r2)/h,         # eq 2.111
                    "RMS outer residual covariance" = RMS(C[lower.tri(C)]), # eq 2.118
                    "RMS inner residual covariance" = RMS(C[lower.tri(Q)]) # eq 2.118
                    )
  }
  else{
    C <- S-stats::fitted(object)
    indices <- list()
  }
  
  C_std <- diag(1/diag(S)) %*% C
  
  indices <- c(indices,list(
                  
                  # SRMR as calculated in SEM. (Hu and Bentler, 1999, p. 3)
                  
                  SRMR = sqrt(sum(C_std[lower.tri(C_std)]^2)/length(C[lower.tri(C_std, diag=TRUE)])),
                  
                  # SRMR calculated ignoring within block residuals from Henserler et al 2014.
                  
                  SRMR_Henseler = RMS(C_std[nativeModel$reflective %*% t(nativeModel$reflective)==0]))
  )
  
  if(observed){
    result<- list(inner = Q, outer = C, indices = indices)
  }
  else{
    result<- list(outer = C, indices = indices)
  }
  
  class(result) <- "matrixplsresiduals"
  result
}

#'@S3method print matrixplsresiduals

print.matrixplsresiduals <- function(x, ...){
  if("outer" %in% names(x)){
    cat("\n Inner model (composite) residual covariance matrix\n")
    print(x$inner, ...)
    cat("\n Outer model (indicator) residual covariance matrix\n")
    print(x$outer, ...)
  }
  else{
    cat("\n Model implied residual covariance matrix\n")
    print(x$outer, ...)
  }
  cat("\n Residual-based fit indices\n")
  print(data.frame(Value = unlist(x$indices)), ...)
}

#'@title Model implied covariance matrix baesd on matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for generic function \code{fitted} computes the model implied
#'covariance matrix by combining  \code{inner}, \code{reflective}, and \code{formative} as a 
#'simultaneous equations system. The error terms are constrained to be uncorrelated and 
#'covariances between exogenous variables are fixed at their sample values. Defining a
#'composite as dependent variable in both inner and formative creates an impossible model
#'and results in an errors
#'
#'@template postestimationFunctions
#'
#'@return a matrix containing the model implied covariances.
#'
#'@export
#'
#'@method fitted matrixpls
#'
#'@S3method fitted matrixpls

fitted.matrixpls <- function(object, ...) {
  
  # Equation numbers in parenthesis refer to equation number in Lohmoller 1989
  
  nativeModel <- attr(object,"model")
  
  # Check that the matrices form a valid model
  
  if(any(rowSums(nativeModel$formative) > 0 & rowSums(nativeModel$formative) > 0))
    stop("Cannot calculate model implied covariance matrix. A composite is a dependent variable in both formative and inner matrices resulting in an impossible model")
  
  S <- attr(object,"S")
  IC <- attr(object,"IC")
  C <- attr(object,"C")
  B <- attr(object,"inner")
  F <- attr(object,"formative")
  L <- attr(object,"reflective")
  
  
  # Matrices containing all regressions and covariances
  # indicators first, then composites
  
  fullB <- rbind(cbind(matrix(0,nrow(S),nrow(S)),L),
                 cbind(F,B))
  
  fullC <- rbind(cbind(S,t(IC)),
                 cbind(IC,C))
  
  exog <- rowSums(fullB) == 0
  
  # Add indicator errors and composite errors
  
  e <- c(diag(S) - rowSums(L * t(IC)),
         1-r2(object))
  
  fullC <- rbind(cbind(fullC,matrix(0, nrow(fullC), sum(!exog))),
                 cbind(matrix(0, sum(!exog), nrow(fullC)), diag(e)[!exog,!exog]))
  
  fullB <- cbind(fullB,diag(length(e))[,! exog])
  fullB <- rbind(fullB, matrix(0,sum(!exog), ncol(fullB)))
  
  # Update exog because fullB is now larger
  exog <- rowSums(fullB) == 0
  
  #
  # Derive the implied covariance matrix using the Bentler-Weeks model
  #
  # Bentler, P. M., & Weeks, D. G. (1980). Linear structural equations with latent variables. Psychometrika, 45(3), 289–308. doi:10.1007/BetaF02293905
  #
  
  # beta (regression paths between dependent variables)
  
  beta <- fullB[!exog,!exog]
  
  
  # gamma (regression paths between dependent and independent  variables)
  
  gamma <- fullB[!exog,exog]
  
  # exogenous variable covariance matrix
  Phi <- fullC[exog,exog]
  
  #
  # Calculate sigma
  #
  
  # Matrix orders
  
  # number of observed dependent
  p = sum(!exog[1:nrow(S)])
  # number of observed independent
  q = nrow(S)-p
  
  # number of independent variables
  n = nrow(Phi)
  # number of dependent variables
  m = nrow(beta)
  
  # Number of observed variables
  r = p + q
  # Number of variables
  s = m + n
  
  # The matrices
  
  Beta <- matrix(0,s,s)
  Gamma <- matrix(0,s,n)
  G <- matrix(0,r,s)
  
  # Identity matrix
  
  I<-diag(nrow(Beta))
  
  # Populate the parameter matrices
  
  Gamma[1:m,] <- gamma
  Gamma[(m+1):s,] <- diag(n)
  
  Beta[1:m,1:m]<-beta
  
  # Populate the selection matrix (columns of all variables  selected to rows of observed variables)
  observed <- NULL
  if(p>0){
    observed <- 1:p
  }
  if(q>0){
    observed <- c(observed,(1:q)+n)
  }
  
  G <- diag(s)[observed,]
  
  # G has dependent variables followed by independent variables. Reorder to match the variables
  # in S
  
  G <- G[order(exog[1:nrow(S)]),]
  
  # Calculate the model implied covariance matrix and return
  
  Sigma = G %*% solve(I-Beta) %*% Gamma %*% Phi %*% t(Gamma) %*% t(solve(I-Beta))%*% t(G)
  rownames(Sigma) <- colnames(Sigma) <- rownames(S)
  
  Sigma
}

#'@S3method print matrixplsresiduals

print.matrixplsresiduals <- function(x, ...){
  cat("\n Inner model (composite) residual covariance matrix\n")
  print(x$inner, ...)
  cat("\n Outer model (indicator) residual covariance matrix\n")
  print(x$outer, ...)
  cat("\n Residual-based fit indices\n")
  print(data.frame(Value = unlist(x$indices)), ...)
}

#'@title R2	for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{r2} computes the squared multiple correlation (R2)
#'for composites predicted by other composites in the model.
#'
#'@template postestimationFunctions
#'
#'@return A named numberic vector containing the R2 values.
#'
#'@export
#'

r2 <- function(object, ...){
  UseMethod("r2")
}

#'@S3method r2 matrixpls

r2.matrixpls <- function(object, ...){
  
  r2 <- rowSums(attr(object,"inner") * attr(object,"C"))
  names(r2) <- colnames(attr(object,"model")$inner)
  class(r2) <- "matrixplsr2"
  r2
}

#'@S3method print matrixplsr2

print.matrixplsr2 <- function(x, ...){
  cat("\n Inner model squared multiple correlations (R2)\n")
  print.table(x, ...)
}

#'@title Goodness of Fit indices for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{gof} computes the Goodness of Fit index for matrixpls results.
#'
#'@template postestimationFunctions
#'
#'@return The Goodness of Fit index.
#'
#'@references
#'
#'Henseler, J., & Sarstedt, M. (2013). Goodness-of-fit indices for partial least squares path modeling. \emph{Computational Statistics}, 28(2), 565–580. doi:10.1007/s00180-012-0317-1
#'
#'
#'@export
#'

gof <- function(object, ...){
  UseMethod("gof")
}

#'@S3method gof matrixpls

gof.matrixpls <- function(object, ...) {
  
  nativeModel <- attr(object,"model")
  IC <- attr(object,"IC")
  S <- attr(object,"S")
  
  exogenousLVs <- rowSums(nativeModel$inner) == 0
  
  
  IC_std <- IC %*% (diag(1/sqrt(diag(S))))
  
  result <- sqrt(mean(IC_std[t(nativeModel$reflective)==1]^2) * 
                   mean(r2(object)[!exogenousLVs]))
  
  # TODO: Relative gof
  # http://www.stat.wmich.edu/wang/561/codes/R/CanCor.R
  
  class(result) <- "matrixplsgof"
  result
}

#'@S3method print matrixplsgof

print.matrixplsgof <- function(x, digits=getOption("digits"), ...){
  cat("\n Absolute goodness of fit:", round(x, digits = digits))
  cat("\n")
}

#'@title Factor loadings matrix from matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{loadings} computes standardized factor loading matrix
#'from the model results.
#'
#'@template postestimationFunctions
#'
#'@return A matrix of estimated factor loadings.
#'
#'@importFrom stats loadings
#'
#'@export
#'
loadings.matrixpls <- function(object, ...){
  UseMethod("loadings")
}

#
#'@S3method loadings matrixpls

loadings <- function(object, ...) {
  nativeModel <- attr(object,"model")
  IC <- attr(object,"IC")
  
  #Standardize
  S <- attr(object,"S")
  IC_std <- IC %*% (diag(1/sqrt(diag(S))))
  
  res <- nativeModel$reflective
  res[res==1] <- t(IC_std)[res==1]
  res
}

#'@title Composite Reliability indices for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{cr} computes Composite Reliability 
#'indices for the model using the formula presented by Fornell and Larcker (1981).
#'
#'@template postestimationFunctions
#'
#'@return A named numeric vector containing the Composite Reliability indices.
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'Aguirre-Urreta, M. I., Marakas, G. M., & Ellis, M. E. (2013). Measurement of composite reliability in research using partial least squares: Some issues and an alternative approach. \emph{The DATA BASE for Advances in Information Systems}, 44(4), 11–43. \href{http://doi.org/10.1145/2544415.2544417}{DOI:10.1145/2544415.2544417}
#'
#'@export
#'

cr <- function(object, ...){
  UseMethod("cr")
}

#'@S3method cr matrixpls

cr.matrixpls <- function(object, ...) {
  
  loadings <- loadings(object)
  reflectiveModel <- attr(object, "model")$reflective
  
  result <- unlist(lapply(1:ncol(loadings), function(col){	
    useLoadings <- abs(loadings[,col][reflectiveModel[,col]==1])
    
    crvalue <- sum(useLoadings)^2/
      (sum(useLoadings)^2 + sum(1-useLoadings^2))
    crvalue
  }))
  
  class(result) <- "matrixplscr"
  names(result) <- colnames(loadings)
  result
}

#'@S3method print matrixplscr

print.matrixplscr <- function(x, ...){
  cat("\n Composite Reliability indices\n")
  print.table(x, ...)
}

#'@title Predict method for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{predict} predict.
#'Predicts the reflective indicators of endogenous latent variables using
#'estimated model and data for the indicators of exogenous latent variables
#'
#'@template postestimationFunctions
#'
#'@param newdata A data frame or a matrix containing data used for prediction.
#'
#'@return a matrix of predicted values for reflective indicators of endogenous latent variables.
#'
#'@references
#'
#'Wold, H. (1974). Causal flows with latent variables: Partings of the ways in the light of NIPALS modelling. \emph{European Economic Review}, 5(1), 67–86. doi:10.1016/0014-2921(74)90008-7
#'
#'
#'@export
#'
#'@method predict matrixpls
#'
#'@S3method predict matrixpls


predict.matrixpls <- function(object, newdata, ...){
  
  nativeModel <- attr(object,"model")
  exog <- rowSums(nativeModel$inner)==0
  W <- attr(object,"W")
  inner <- attr(object,"inner")
  
  reflective <- nativeModel$reflective
  reflective[which(reflective==1)] <- object[grep("=~",names(object))]
  
  # Reorder the variables in newdata
  data <- NULL
  
  for(name in colnames(W)){
    if(! name %in% colnames(newdata)){
      data <- cbind(data,NA)
    }
    else{
      data <- cbind(data,newdata[,name])
    }
  }
  
  colnames(data) <- colnames(W)
  
  LVScores <- as.matrix(data) %*% t(W)
  
  # Predict endog LVs usign reduced from equations
  Beta <- inner[! exog, ! exog]
  Gamma <- inner[! exog, exog]
  
  LVScoresEndo <- (LVScores[,exog] %*% t(Gamma)) %*% solve(diag(nrow(Beta))- t(Beta))
  
  LVScores[,! exog] <- LVScoresEndo
  
  LVScores %*% t(reflective)
}

#'@title Average Variance Extracted indices for matrixpls results
#'
#'@description
#'
#'The \code{matrixpls} method for the generic function \code{ave} computes Average Variance Extracted 
#'indices for the model using the formula presented by Fornell and Larcker (1981).
#'
#'@template postestimationFunctions
#'
#'@return A list containing the Average Variance Extracted indices in the first position and the differences
#'between aves and largest squared correlations with other composites in the second position.
#'
#'
#'@references
#'
#'Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of marketing research}, 18(1), 39–50.
#'
#'@export
#'

ave <- function(object, ...){
  UseMethod("ave")
}


#'@S3method ave matrixpls

ave.matrixpls <- function(object, ...) {
  
  loadings <- loadings(object, standardized = TRUE)
  reflectiveModel <- attr(object, "model")$reflective
  
  aves <- unlist(lapply(1:ncol(loadings), function(col){
    useLoadings <- loadings[,col][reflectiveModel[,col]==1]
    
    sum(useLoadings^2)/
      (sum(useLoadings^2) + sum(1-useLoadings^2))
    
  }))
  
  names(aves) <- colnames(loadings)
  
  C <- attr(object,"C")
  C2 <- C^2
  diag(C2) <- 0
  maxC <- apply(C2,1,max)
  
  aves_correlation <- aves - maxC
  
  names(aves_correlation) <- colnames(loadings)
  
  result = list(ave = aves,
                ave_correlation = aves_correlation)
  
  class(result) <- "matrixplsave"
  result
}

#'@S3method print matrixplsave

print.matrixplsave <- function(x, ...){
  cat("\n Average Variance Extracted indices\n")
  print.table(x$ave, ...)
  cat("\n AVE - largest squared correlation\n")
  print.table(x$ave_correlation, ...)
}

#'@title Heterotrait-monotrait ratio 
#'
#'@description 
#'
#'The \code{matrixpls} method for the generic function \code{htmt} computes Heterotrait-monotrait ratio 
#'for the model using the formula presented by Henseler et al (2014).
#'
#'@template postestimationFunctions
#'
#'@return Heterotrait-monotrait ratio as a scalar
#'
#'
#'@references
#'
#'Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for assessing discriminant validity in variance-based structural equation modeling. \emph{Journal of the Academy of Marketing Science}, 43(1), 115–135.
#'
#'
#'@export
#'

htmt <- function(object, ...){
  
  S <- attr(object,"S")
  reflective <- attr(object,"reflective")!=0
  
  # Mean within- and between-block correlations
  
  meanBlockCor <- (t(reflective) %*% (S*lower.tri(S)) %*% reflective) / 
    (t(reflective) %*% lower.tri(S) %*% reflective)
  
  # Choose the blocks that have at least two reflective indicators
  i <- which(colSums(reflective) > 1)
  meanBlockCor <- meanBlockCor[i,i]
  
  # Calculate the indices
  htmt <- meanBlockCor*lower.tri(meanBlockCor) /
    sqrt(diag(meanBlockCor) %o% diag(meanBlockCor))
  
  htmt
}

#'@title Summary of model fit of PLS model
#'
#'@description
#'
#'\code{fitSummary} computes a list of statistics
#'that are commonly used to access the overall fit of the PLS model.
#'
#'@template postestimationFunctions
#
#'@return A list containing a selection of fit indices.
#'
#'@export
#'

fitSummary <- function(object){
  
  # Returns the minumum or NA if all are NA
  
  m <- function(x){
    x <- stats::na.exclude(x)
    if(length(x) == 0 ) NA
    else min(x)
  }
  
  ret <- list("Min CR" = m(cr(object)),
              "Min AVE" = m(ave(object)$ave),
              "Min AVE - sq. cor" = m(ave(object)$ave_correlation),
              "Goodness of Fit" = gof(object),
              SRMR = stats::residuals(object)$indices$srmr,
              "SRMR (Henseler)" = stats::residuals(object)$indices$srmr_Henseler)
  
  class(ret) <- "matrixplfitsummary"
  ret
}

#'@S3method print matrixplfitsummary

print.matrixplfitsummary <- function(x, ...){
  cat("\n Summary indices of model fit\n")
  print(data.frame(Value = unlist(x)), ...)
}
