#' @name narrow
#' @title Title.
#' @description Description.
#' @param data.vector :  vector of phenotypic values
#' @param geno.vector :  vector (character or factor) of genotypes
#' @param K           :  kinship matrix, of class matrix. Column-names should contain all occurring genotypes (?)
#' @param covariates  :  additional covariates (which should not include an intercept, which is included automatically)
#' @param average.number.of.replicates :
#' @param alpha : confidence level, for the (1-alpha) confidence intervals
#' @return blabla
#' @author Willem Kruijer \email{willem.kruijer@@wur.nl}
#' @details blabla.
#' @examples blabla.
#' @export
#' 
#' 

marker_h2_means <- function(data.vector,geno.vector,K,Dm=NULL,alpha=0.05,
                           eps=1e-06,max.iter=100,fix.h2=FALSE,h2=0.5,grid.size=99){

# OUTPUT :
# heritability, genetic - and residual variance

    if (is.null(rownames(K)) | is.null(colnames(K))) {
      stop('K should have row- and column-names corresponding to the levels of geno.vector')
    }
  
  
    if (is.null(Dm)) {
      Dm <- diag(nrow(K))
      rownames(Dm) <- colnames(Dm) <- colnames(K)
      cat('Warning: it is assumed that each phenotypic values corresponds to an observation on a single individual.','\n')
      cat('If each phenotypic value is the mean of observations on r genetically identical individuals in a completely randomized design,','\n')
      cat('Dm should be diag(1/r,nrow=nrow(K)). More generally, in other designs, Dm should be the covariance matrix of the genotypic means. ','\n')
    } 
    
    if (is.null(rownames(Dm)) | is.null(colnames(Dm))) {
      stop('Dm should have row- and column-names corresponding to the levels of geno.vector')
    }
    
    her.frame            <- data.frame(dat=data.vector,genotype=geno.vector)
    her.frame$genotype   <- as.character(her.frame$genotype)

    her.frame     <- her.frame[!is.na(her.frame$dat),]

    # re-scale the kinship matrix, for those genotypes for which phenotypic data are available
    K <- K[unique(her.frame$genotype),unique(her.frame$genotype)]
    K <- K / KinshipTransform(K)

    #if (is.null(Dm)) {Dm <- K} else {
    Dm <- Dm[unique(her.frame$genotype),unique(her.frame$genotype)]

    VALUE <- as.data.frame(her.frame$dat)
    names(VALUE)[1] <- 'value'
    
    ######################
            
    reml.obj <- AI_algorithm_weights(Y=data.frame(genotype=her.frame$genotype,VALUE),K=K,Dm=Dm,eps=eps,max.iter=max.iter,fix.h2=fix.h2,h2=h2)

    if (reml.obj$converge == FALSE) {
      h2.grid     <- (1:grid.size)/(grid.size+1)
      Y <- rep(NA,grid.size)
      for (i in 1:grid.size) {
        Y[i]           <- AI_algorithm_weights(Y=data.frame(genotype=her.frame$genotype,VALUE),K=K,Dm=Dm,eps=eps,max.iter=max.iter,fix.h2=TRUE,h2=h2.grid[i])$loglik
      }
      reml.obj <- AI_algorithm_weights(Y=data.frame(genotype=her.frame$genotype,VALUE),K=K,Dm=Dm,eps=eps,max.iter=max.iter,fix.h2=TRUE,h2=min(h2.grid[which.max(Y)]))
    
    } else{
      Y <- numeric(0)    
    }

    reml.obj$gammas        <- c(reml.obj$varcomp[1],1,reml.obj$varcomp[2])    
    names(reml.obj$gammas) <- c("giv(genotype, var = TRUE).giv","R!variance","R!units.var")
    reml.obj$gammas.type   <- c(2,1,2)
    names(reml.obj$gammas.type) <- names(reml.obj$gammas)
    reml.obj$ai            <- reml.obj$inv.ai[lower.tri(reml.obj$inv.ai,diag=TRUE)] #reml.obj$inv.ai
    reml.obj$ai            <- c(reml.obj$ai[1],0,0,reml.obj$ai[2],0,reml.obj$ai[3])

    var.comp <- reml.obj$varcomp

    inv.ai <- reml.obj$inv.ai

    h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[2])
    st.error1 <- as.numeric(pin(reml.obj,h2 ~ V1 / (V1 + V3))[2])
    conf.int1 <- c(h2.estimate - qnorm(1-alpha/2)*st.error1,h2.estimate + qnorm(1-alpha/2)*st.error1)

    st.error2 <- as.numeric(pin(reml.obj,h2 ~ log(V3 / V1))[2])
    conf.int2 <- c(1 / (1 + exp(log(var.comp[2] /var.comp[1]) + qnorm(1-alpha/2)*st.error2)),1 / (1 + exp(log(var.comp[2]/var.comp[1]) - qnorm(1-alpha/2)*st.error2)))

return(list(va=var.comp[1],ve=var.comp[2],h2=h2.estimate,conf.int1=conf.int1,conf.int2=conf.int2,
            inv.ai=reml.obj$inv.ai,loglik=reml.obj$loglik,loglik.vector=Y))

}

# inv.ai.varcomp=reml.obj$inv.ai.new