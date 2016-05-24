#' @name AI_algorithm
#' @title Use the AI-algortihm (average ingormation; Gilmour et al 1995) to obtain REML-estimates in a mixed model containing a random genetic effect and iid errors.   
#' @description This function is used for mixed-model based estimation of narrow-sense heritability. 
#'              The input consists of a data.frame containing the individuals as rows, and with columns indicating their genotype, phenotype and possible covariates.   
#'              The data-frame can contain replicates of the same genotype. In addition, a kinship matrix is required, defining the genetic relatedness between the genotypes. 
#'              The vector of genetic effects is multivariate normal with...
#' @param w weigts 1/r_i , i=1...n
#' @param Y  gsdfgsdfg
#' @param K  sdfgsdfgsdfgs
#' @param eps sdfgsdfgsg
#' @param max.iter sdfgsdfgsdfg  
#' @return a list with the following components:  a vector varcomp, containing respectively the estimated genetic- and residual variance; 
#'         a 2 x 2 matrix inv.ai, which is the inverse of the average information matrix; n.iter, the number of iterations until convergence; eps, the numerical precision
#'  
#' @author Willem Kruijer \email{willem.kruijer@@wur.nl}
#' 
#' @details
#' @examples
#' @export
#
#
AI_algorithm_weights <- function(Y,K,Dm,eps = 0.000000001,max.iter=100,
                                 escape=0.3,fix.h2=FALSE,h2=0.5) {
#
# We use the AI-alogortihm as described in Yang,  Hong Lee,Goddard and Visscher 2011 (GCTA paper)
#
# Model: a vector y (the second column of Y) has a normal distribution with covariance \sigma_A^2 K + \sigma_E^2 Dm,
#        where K is the kinship matrix.   
#        We will refer to Dm as the weight matrix. In case of a single observation per genotype, it is the identity
#        In case of a balanced CRD with r replicates, it is the identity divided by r
#        In the more general case, it is the covariance matrix of the BLUEs of the genotypes, by divided sigma^2, the estimated residual variance
#        The preceding quantities are obtained from fitting a linear model without an intercept, a fixed factor genotype (number of levels equal to the number of genotypes)
#        and possibly extra (fixed) terms (i.e. the intercept is included in the factor genotype)

  if (escape < 0 | escape >=0.5) {escape <- 0}

  if (!(class(K) %in% c('matrix','data.frame'))) {stop("K should be of class matrix or data.frame")}
  if (class(K)=='data.frame') {class(K) <- 'matrix'}
  if (is.null(colnames(K)) | is.null(rownames(K))) {stop("K should have column and row names")}

  if (!(class(Dm) %in% c('matrix','data.frame'))) {stop("Dm should be of class matrix or data.frame")}
  if (class(Dm)=='data.frame') {class(Dm) <- 'matrix'}
  if (is.null(colnames(Dm)) | is.null(rownames(Dm))) {stop("Dm should have column and row names")}

  if (class(Y)!='data.frame') {stop("Y should be of class data.frame")}
  if (ncol(Y) < 1) {stop("Y should be a data.frame with at least two columns: \n The first column of Y should contain the genotype labels, and should be of class factor or character. \n The second column should contain the phenotypic values.")}

  if (!(class(Y[,1]) %in% c('factor','character'))) {stop("The first column of Y should contain the genotype labels, and should be of class factor or character.")}
  if (sum( unlist(lapply(Y,class)) %in% c('integer','numeric')) < (ncol(Y) - 1)) {stop("The first column of Y should contain the genotype labels; \n All other columns should be of type numeric or integer.")}

  if (max(table(Y[,1])) > 1) {stop("The first column of Y should contain the genotype labels, and each genotype should occur at most once.")}

  # If there are no covariates, add an intercept
  if (ncol(Y)==2) {Y <- cbind(Y,rep(1,nrow(Y)))}

  var.cols <- apply(as.data.frame(Y[,-(1:2)]),2,var)

  # If there is no intercept (i.e. the original data-frame did have covariaates, but no intercept), 
  #   add an intercept to the data.frame
  if (sum(var.cols==0)==0) {Y <- cbind(Y,rep(1,nrow(Y)))}

  if (sum(var.cols==0)>1) {
    first.intercept <- which(var.cols==0)[1] + 2
    Y <- Y[,-((which(var.cols==0)[1] + 2)[-1])]
  }

  Y <- Y[!is.na(Y[,2]),]
  Y <- Y[!is.na(Y[,1]),]
  names(Y)[1:2] <- c('genotype','value')
  if (ncol(Y) > 2) {
    Y <- Y[apply(as.data.frame(Y[,-(1:2)]),1,function(x){sum(is.na(x))})==0,]
  }

  if (!all(Y$genotype %in% colnames(K))) {
    no.kinship <- unique(Y$genotype[!(Y$genotype %in% colnames(K))])
    Y <- Y[Y$genotype %in% colnames(K),]
    if (nrow(Y)==0) {stop("Kinship information is not available for any of the genotypes in Y")}
    warning(paste("K does not contain kinship information for genotypes",no.kinship,"; these genotypes are removed from Y."))
  }

  if (!all(Y$genotype %in% colnames(Dm))) {
    no.Dm <- unique(Y$genotype[!(Y$genotype %in% colnames(Dm))])
    Y <- Y[Y$genotype %in% colnames(Dm),]
    if (nrow(Y)==0) {stop("Weights are not available for any of the genotypes in Y")}
    warning(paste("The Dm (weight) matrix does not contain information for genotypes",no.Dm,"; these genotypes are removed from Y."))
  }

  ########

  y <- as.matrix(Y[,2])
  X <- as.matrix(Y[,-(1:2)])

  K          <-  K[as.character(Y$genotype),as.character(Y$genotype)]
  Dm         <- Dm[as.character(Y$genotype),as.character(Y$genotype)]

  Dm_inv     <- ginv(Dm)

  N          <- nrow(Y)

  n.iter     <- 0

  converge <- TRUE
      
  ########
  
  if (!fix.h2) {
    
    # varcomp has two components: the genetic variance sigma_g^2 (varcomp[1]) and the residual variance sigma_e^2 (varcomp[2])
    varcomp.old <- c(0,0)
    # initialize
    varcomp <- rep(var(as.numeric(y)) / 2,2) # starting values genetic- and residual variance

    # WE DON'T use a first EM-step as described in Yang,  Hong Lee,Goddard and Visscher 2011 (GCTA paper)
    # ... and directly start with the AI-algorthm

    while (sum((varcomp.old - varcomp)^2) > eps & n.iter < max.iter) {
  
      Vinv   <- ginv(varcomp[1] * K + varcomp[2] * Dm)
      P      <- Vinv - (Vinv %*% X %*% ginv(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv)
      
      deltaL <- -0.5 * c(sum(diag(P %*% K)) - t(y) %*% P %*% K %*% P %*% y,sum(diag(P %*% Dm)) - t(y) %*% P %*% Dm %*% P %*% y)
      
      AI     <- matrix(0,2,2)
      AI[1,2]<- AI[2,1] <- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% Dm %*% P %*% y)
      AI[1,1]<- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% K %*% P %*% y)
      AI[2,2]<- 0.5 * as.numeric(t(y) %*% P %*% Dm %*% P %*% Dm %*% P %*% y)
  
      varcomp.old <- varcomp
      varcomp     <- varcomp + as.numeric(ginv(AI) %*% deltaL)
  
      # reset negative values
      #if (sum(varcomp < 0) ==2) {rep(var(as.numeric(y)) / 2,2)}
      if (sum(varcomp < 0) ==2) {
        guess <- runif(n=1)      
        varcomp <- c(guess * var(as.numeric(y)), (1-guess) * var(as.numeric(y)) ) 
      }
  
      if (escape > 0) { 
        guess <- runif(1,0,escape)
        if (varcomp[1] < 0) {varcomp <- c(guess * var(as.numeric(y)), (1-guess) * var(as.numeric(y)) ) }
        if (varcomp[2] < 0) {varcomp <- c((1-guess) * var(as.numeric(y)), guess * var(as.numeric(y)) ) }
      }
      if (escape==0) {
        if (varcomp[1] < 0) {varcomp <- c(0.001 * var(as.numeric(y)), 0.999 * var(as.numeric(y)) ) }
        if (varcomp[2] < 0) {varcomp <- c(0.999 * var(as.numeric(y)), 0.001 * var(as.numeric(y)) ) }
      }
      
      n.iter <- n.iter + 1
  
      cat('Iteration: ',n.iter,'\t',varcomp,'\n')
    }

    if (sum((varcomp.old - varcomp)^2) > eps & n.iter == max.iter) {
      warning('No convergence in AI-algorithm.')
      converge <- FALSE
    }
    
    # FINAL ESTIMATE AI-MATRIX
    Vinv   <- ginv(varcomp[1] * K + varcomp[2] * Dm)
    P      <- Vinv - (Vinv %*% X %*% ginv(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv)    
    AI     <- matrix(0,2,2)
    AI[1,2]<- AI[2,1] <- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% Dm %*% P %*% y)
    AI[1,1]<- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% K %*% P %*% y)
    AI[2,2]<- 0.5 * as.numeric(t(y) %*% P %*% Dm %*% P %*% Dm %*% P %*% y)
    
  } else {

    # profile likelihood, see e.g. Crainiceanu Ruppert (2004), bottom of p.3  (applied here without the lambda parametrization)
    VunscaledInv <- ginv(h2 * K + (1-h2) * Dm)
    betaHat      <- ginv(t(X) %*% VunscaledInv %*% X) %*% t(X) %*% VunscaledInv %*% y
    sigma2       <- (1 / ncol(VunscaledInv)) * t(y - X %*% betaHat) %*% VunscaledInv %*% (y - X %*% betaHat)
    varcomp <- c(h2 * sigma2, (1-h2) * sigma2)

    Vinv   <- ginv(varcomp[1] * K + varcomp[2] * Dm)
    P      <- Vinv - (Vinv %*% X %*% ginv(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv)
      
    deltaL <- -0.5 * c(sum(diag(P %*% K)) - t(y) %*% P %*% K %*% P %*% y,sum(diag(P %*% Dm)) - t(y) %*% P %*% Dm %*% P %*% y)
      
    AI     <- matrix(0,2,2)
    AI[1,2]<- AI[2,1] <- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% Dm %*% P %*% y)
    AI[1,1]<- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% K %*% P %*% y)
    AI[2,2]<- 0.5 * as.numeric(t(y) %*% P %*% Dm %*% P %*% Dm %*% P %*% y)    

  }

  ######################################
  
  # now redefine AI, to get the AI-matrix as defined by Gilmour et al 1995 (we denote their P with Pnew, to distinguish it from
  # the P used above, where we followed the notation of Yang,  Hong Lee,Goddard and Visscher 2011 (GCTA paper)

  W  <- cbind(X,diag(N))

  tX_Z <- t(X) %*% Dm_inv #apply(Z,2,sum)

  Cm <- rbind(cbind(t(X) %*% Dm_inv %*% X,tX_Z),cbind(t(tX_Z),Dm_inv + (varcomp[2]/varcomp[1]) * ginv(K)))

  Pnew <- Dm_inv - Dm_inv %*% W %*% ginv(Cm) %*% t(W) %*% Dm_inv
  
  ###########
  
  AInew      <- AI
  AInew[1,2] <- AInew[2,1] <- 0.5 * as.numeric(t(y) %*% Pnew %*% K %*% Pnew %*% y) / (varcomp[2])^2
  AInew[1,1] <- 0.5 * as.numeric(t(y) %*% Pnew %*% K %*% Pnew %*% K %*% Pnew %*% y) / (varcomp[2])
  AInew[2,2] <- 0.5 * as.numeric(t(y) %*% Pnew %*% y) / (varcomp[2])^3

  ############

  V <- varcomp[1] * K + varcomp[2] * Dm

  # log(det(V)) ==> 2*sum(log(diag(chol(V))))

  loglik <-  -0.5 * (2*sum(log(diag(chol(V))))  + log(det(t(X) %*% Vinv %*% X)) +  t(y) %*% P %*% y)   

  return(list(varcomp=varcomp,inv.ai=ginv(AI),inv.ai.new=ginv(AInew),n.iter=n.iter,eps=eps,loglik=loglik,converge=converge))
}


