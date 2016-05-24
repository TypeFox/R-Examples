#' @name AI_algorithm
#' @title Use the AI-algortihm (average ingormation; Gilmour et al 1995) to 
#'        obtain REML-estimates in a mixed model containing a random genetic effect and i.i.d. errors.   
#' @description This function is used for mixed-model based estimation of narrow-sense heritability using individual plant or plot data. 
#'              The input consists of a data.frame containing the individuals as rows, and with columns indicating their genotype, phenotype and possible covariates.   
#'              A kinship matrix is required, defining the genetic relatedness between the genotypes. 
#' @param Y  gsdfgsdfg
#' @param K  sdfgsdfgsdfgs
#' @param eps sdfgsdfgsg
#' @param max.iter sdfgsdfgsdfg  
#' @return a list with the following components:  a vector varcomp, containing respectively the estimated genetic- and residual variance; 
#'         a 2 x 2 matrix inv.ai, which is the inverse of the average information matrix; n.iter, the number of iterations until convergence; eps, the numerical precision
#'  
#' @author Willem Kruijer \email{willem.kruijer@@wur.nl}
#' 
#' @details ...
#' @examples ...
#' @export              return(list(varcomp=varcomp,inv.ai=ginv(AI),inv.ai.new=ginv(AInew),n.iter=n.iter,eps=eps,loglik=loglik,converge=converge))
AI_algorithm <- function(Y,K,eps = 0.000000001,max.iter=25,escape=0.1,fix.h2=FALSE,h2=0.5) {   # X=matrix(rep(1,nrow(Y)))

  # Y=data.frame(genotype=her.frame$genotype,VALUE)
  #require(MASS)
  # eps : numerical error margin, used as convergence criterion
  # max.iter : maximal number of iterations

  # Y must be
  # - a vector (class numeric or integer), with

  #if (class(Y)=='integer') {class(Y)=='numeric'}

  #if (!(class(Y) %in% c('data.frame','numeric'))) {stop("Y should be of class data.frame, numeric or integer.")}

  #if (class(Y)=='numeric') {
  #  if (is.null(names(Y))) {stop("If Y is a vector of class numeric or integer, it should have names, i.e. the genotype names used in the kinship matrix.")}
  #  Y <- data.frame(genotype=names(Y),value=Y)
  #}

  #if (class(Y)=='data.frame') {
  #}

  #if (class(Y)=='matrix') {class(Y) <- 'data.frame'}
  #Y=data.frame(genotype=her.frame$genotype,value=her.frame$dat) ; eps = 0.000000001; max.iter=25

  # Y=data.frame(genotype=her.frame$genotype,VALUE) ; eps = 0.000000001 ; max.iter=25

  # escape : either zero, or 

  if (escape < 0 | escape >=0.5) {escape <- 0}

  if (class(Y)!='data.frame') {stop("Y should be of class data.frame")}
  if (!(class(K) %in% c('matrix','data.frame'))) {stop("K should be of class matrix or data.frame")}
  if (class(K)=='data.frame') {class(K) <- 'matrix'}
  if (is.null(colnames(K)) | is.null(rownames(K))) {stop("K should have column and row names")}

  if (class(Y)!='data.frame') {stop("Y should be of class data.frame")}
  if (ncol(Y) < 2) {stop("Y should be a data.frame with at least two columns: \n The first column of Y should contain the genotype labels, and should be of class factor or character. \n The second column should contain the phenotypic values")}

  if (!(class(Y[,1]) %in% c('factor','character'))) {stop("The first column of Y should contain the genotype labels, and should be of class factor or character.")}
  if (sum( unlist(lapply(Y,class)) %in% c('integer','numeric')) < (ncol(Y) - 1)) {stop("The first column of Y should contain the genotype labels; \n All other columns should be of type numeric or integer.")}

  # If there are no covariates, add an intercept
  if (ncol(Y)==2) {Y <- cbind(Y,rep(1,nrow(Y)))}

  var.cols <- apply(as.data.frame(Y[,-(1:2)]),2,var)

  # If there is no intercept (i.e. the original data-frame did have covariaates, but no intercept), 
  #   add an intercept to the data.frame
  if (sum(var.cols==0)==0) {Y <- cbind(Y,rep(1,nrow(Y)))}

  if (sum(var.cols==0)>1) {
    first.intercept <- which(var.cols==0)[1] + 2
    #Y <- Y[,-(first.intercept)]
    Y <- Y[,-((which(var.cols==0)[1] + 2)[-1])]
  }

  Y <- Y[!is.na(Y[,2]),]
  Y <- Y[!is.na(Y[,1]),]
  names(Y)[1:2] <- c('genotype','value')
  if (ncol(Y) > 2) {
    Y <- Y[apply(as.data.frame(Y[,-(1:2)]),1,function(x){sum(is.na(x))})==0,]
  }
  # other names ... ?

  if (!all(Y$genotype %in% colnames(K))) {
    no.kinship <- unique(Y$genotype[!(Y$genotype %in% colnames(K))])
    Y <- Y[Y$genotype %in% colnames(K),]
    if (nrow(Y)==0) {stop("Kinship information is not available for none of the individuals in Y")}
    warning(paste("K does not contain kinship information for genotypes",no.kinship,"; these genotypes are removed from Y."))
  }

  ########

  y <- as.matrix(Y[,2])
  X <- as.matrix(Y[,-(1:2)])

  K.original <- K
  K          <- K[as.character(Y$genotype),as.character(Y$genotype)]

  N   <- nrow(Y)
  
  n.iter <- 0
  
  converge <- TRUE
    
  ########
  

  if (!fix.h2) {
    
    # varcomp has two components: the genetic variance sigma_g^2 (varcomp[1]) and the residual variance sigma_e^2 (varcomp[2])
    varcomp.old <- c(0,0)
    # initialize
    varcomp <- rep(var(as.numeric(y)) / 2,2) # starting values genetic- and residual variance
      
    # Use the AI-alogortihm as described in Yang,  Hong Lee,Goddard and Visscher 2011 (GCTA paper)
  
    # initial step, following Yang et al (2011; GCTA paper) : update genetic variance following the EM-algorithm
    Vinv           <- ginv(varcomp[1] * K + varcomp[2] * diag(N))
    # before 17-1:
    #P              <- Vinv - (Vinv %*% (X %*% t(X)) %*% Vinv) / as.numeric(t(X) %*% Vinv %*% X)
    P              <- Vinv - (Vinv %*% X %*% ginv(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv)
  
    varcomp.old[1] <- varcomp[1]
    varcomp[1]     <- (1/N) * ((varcomp[1])^2 * as.numeric(t(y) %*% P %*% K %*% P %*% y) + sum(diag(varcomp[1] * diag(N) - (varcomp[1])^2 * P %*% K)) )
    varcomp[2]     <- max(0,var(as.numeric(y)) - varcomp[1])

    # reset negative values
    if (sum(varcomp < 0) > 0) {
      guess <- runif(n=1)      
      varcomp <- c(guess * var(as.numeric(y)), (1-guess) * var(as.numeric(y)) ) 
    }

    #if (varcomp[1] < 0) {varcomp <- c(0.1 * var(as.numeric(y)), 0.9 * var(as.numeric(y)) ) }
    #if (varcomp[2] < 0) {varcomp <- c(0.9 * var(as.numeric(y)), 0.1 * var(as.numeric(y)) ) }  
      
    # continue with the AI-algorthm
    
    while (sum((varcomp.old - varcomp)^2) > eps & n.iter < max.iter) {
  
      Vinv           <- ginv(varcomp[1] * K + varcomp[2] * diag(N))
      P              <- Vinv - (Vinv %*% X %*% ginv(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv)
      deltaL         <- -0.5 * c(sum(diag(P %*% K)) - t(y) %*% P %*% K %*% P %*% y, sum(diag(P)) - t(y) %*% P %*% P %*% y )
      
      AI     <- matrix(0,2,2)
      AI[1,2]<- AI[2,1] <- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% P %*% y)
      AI[1,1]<- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% K %*% P %*% y)
      AI[2,2]<- 0.5 * as.numeric(t(y) %*% P %*% P %*% P %*% y)
  
      varcomp.old <- varcomp
      varcomp     <- varcomp + as.numeric(ginv(AI) %*% deltaL)
  
      # reset negative values
      if (sum(varcomp < 0) ==2) {rep(var(as.numeric(y)) / 2,2)}
  
  
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
    Vinv           <- ginv(varcomp[1] * K + varcomp[2] * diag(N))
    P              <- Vinv - (Vinv %*% X %*% ginv(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv)    
    AI     <- matrix(0,2,2)
    AI[1,2]<- AI[2,1] <- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% P %*% y)
    AI[1,1]<- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% K %*% P %*% y)
    AI[2,2]<- 0.5 * as.numeric(t(y) %*% P %*% P %*% P %*% y)

  } else {            
    varcomp <- c(h2 * var(as.numeric(y)), (1-h2) * var(as.numeric(y)))
    
    Vinv   <- ginv(varcomp[1] * K + varcomp[2] * diag(N))
    P              <- Vinv - (Vinv %*% X %*% ginv(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv)
      
    deltaL <-  -0.5 * c(sum(diag(P %*% K)) - t(y) %*% P %*% K %*% P %*% y, sum(diag(P)) - t(y) %*% P %*% P %*% y )
      
    AI     <- matrix(0,2,2)
    AI[1,2]<- AI[2,1] <- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% P %*% y)
    AI[1,1]<- 0.5 * as.numeric(t(y) %*% P %*% K %*% P %*% K %*% P %*% y)
    AI[2,2]<- 0.5 * as.numeric(t(y) %*% P %*% P %*% P %*% y)
    
  }

  ######################################
  
  # now redefine AI, to get the AI-matrix as defined by Gilmour et al 1995 (we denote their P with Pnew, to distinguish it from
  # the P used above, where we followed the notation of Yang,  Hong Lee,Goddard and Visscher 2011 (GCTA paper)

  Z <- matrix(0,nrow(y),ncol(K.original))
  for (i in 1:nrow(y)) {Z[i,which(rownames(K.original)==as.character(Y$genotype[i]))] <- 1}

  W  <- cbind(X,Z)
  tX_Z <- t(X) %*% Z 
  Cm <- rbind(cbind(t(X) %*% X,tX_Z),cbind(t(tX_Z),t(Z) %*% Z + (varcomp[2]/varcomp[1]) * GINV(K.original)))
  Pnew <- diag(N) - W %*% GINV(Cm) %*% t(W)

  AInew      <- AI
  AInew[1,2] <- AInew[2,1] <- 0.5 * as.numeric(t(y) %*% Pnew %*% K %*% Pnew %*% y) / (varcomp[2])^2
  AInew[1,1]<- 0.5 * as.numeric(t(y) %*% Pnew %*% K %*% Pnew %*% K %*% Pnew %*% y) / (varcomp[2])
  AInew[2,2]<- 0.5 * as.numeric(t(y) %*% Pnew %*% y) / (varcomp[2])^3

  ## 
  
  V <- varcomp[1] * K + varcomp[2] * diag(N)

  loglik <-  -0.5 * (2*sum(log(diag(chol(V))))  + log(det(t(X) %*% Vinv %*% X)) +  t(y) %*% P %*% y)   


  #return(list(varcomp=varcomp,inv.ai=ginv(AI),inv.ai.new=ginv(AInew),n.iter=n.iter,eps=eps,loglik=loglik,converge=converge))
# 17 sept. 2014: ginv ==> ginv
return(list(varcomp=varcomp,inv.ai=GINV(AI),inv.ai.new=GINV(AInew),n.iter=n.iter,eps=eps,loglik=loglik,converge=converge))
}


