##
##  PURPOSE:   Re-labeling of the MCMC output.
##             * method for objects of class GLMM_MCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   26/02/2010
##
##  FUNCTION:  NMixRelabel.GLMM_MCMC (26/02/2010) 
##
## ======================================================================

## *************************************************************
## NMixRelabel.GLMM_MCMC
## *************************************************************
NMixRelabel.GLMM_MCMC <- function(object, type = c("mean", "weight", "stephens"), par,
                                  prob = c(0.025, 0.5, 0.975), keep.comp.prob = FALSE, info, silent = FALSE, ...)
{
  thispackage <- "mixAK"

  silent <- as.logical(silent[1])
  if (is.na(silent)) silent <- FALSE
  
  if (!object$dimb) stop("No random effects in object, nothing to re-label.")
  LTp <- object$dimb * (object$dimb + 1)/2
  I   <- object$Cpar$I
  l_alpha <- sum(object$p) + sum(object$fixed.intercept)
  
  
  ##### Determine re-labeling algorithm to use and additional parameters
  ##### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RAlg <- NMixRelabelAlgorithm(type=type, par=par, dim=object$dimb)
  object$relabel_b <- RAlg$relabel                                     ## resulting re-labeling


  ##### Parameters of MCMC
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  if (is.null(object$K_b) | is.null(object$w_b) | is.null(object$mu_b) | is.null(object$Li_b) | is.null(object$Q_b) | is.null(object$Sigma_b)){
    stop("object does not contain sampled values of mixture related parameter(s)")
  }  
  if (object$R["Rc"] & is.null(object$sigma_eps)){
    stop("object does not contain sampled values of sigma(eps)")
  }
  if (l_alpha & is.null(object$alpha)){
    stop("object does not contain sampled values of fixed effects")
  }  
  
  keepMCMC <- length(object$w_b) / object$K_b[1]
  if (missing(info)) info <- keepMCMC
  if (info <= 0 | info > keepMCMC) info <- keepMCMC

  
  ##### Some input checks
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  if (object$prior.b$priorK != "fixed") stop("only implemented for models with a fixed number of mixture components")


  ##### Needed length or Pr_b_b
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  lPr_b <- object$Cpar$I * object$K_b[1] * keepMCMC
  
  ##### Perform re-labeling
  ##### ++++++++++++++++++++++++++++++++++++++++++++++
  l_nchange <- ifelse(RAlg$Ctype <= 2, 1, RAlg$relabel$par$maxiter)
  
  MCMC <- .C("GLMM_NMixRelabel",
             type           = as.integer(RAlg$Ctype),
             iparam         = as.integer(RAlg$iparam),
             nonSilent      = as.integer(!silent),
             Y_c            = as.double(object$Cpar$Y_c),
             Y_d            = as.integer(object$Cpar$Y_d),
             R_cd           = as.integer(object$Cpar$R_cd),
             dist           = as.integer(object$Cpar$dist),
             I              = as.integer(object$Cpar$I),
             n              = as.integer(object$Cpar$n),
             X              = as.double(object$Cpar$X),
             Z              = as.double(object$Cpar$Z),
             p_fI_q_rI      = as.integer(object$Cpar$p_fI_q_rI),
             shiftScale_b   = as.double(c(object$scale.b$shift, object$scale.b$scale)),
             keepMCMC       = as.integer(keepMCMC),
             info           = as.integer(info),
             tune_scale_b   = as.double(object$Cpar$tune_scale_b),
             chsigma_eps    = if (object$R["Rc"]) as.double(t(object$sigma_eps)) else as.double(0),
             distribution_b = as.integer(object$Cpar$priorInt_b["distribution"]),
             K_b            = as.integer(object$K_b[1]),
             chw_b          = as.double(t(object$w_b)),
             chmu_b         = as.double(t(object$mu_b)),
             chQ_b          = as.double(t(object$Q_b)),
             chSigma_b      = as.double(t(object$Sigma_b)),
             chLi_b         = as.double(t(object$Li_b)),
             chdf_b         = if (object$Cpar$priorInt_b["distribution"] == 1) as.double(t(object$df_b)) else as.double(0),
             chalpha        = if (l_alpha) as.double(t(object$alpha)) else as.double(0),             
             chorder_b      = integer(object$K_b[1] * keepMCMC),
             chrank_b       = integer(object$K_b[1] * keepMCMC),
             b              = as.double(t(object$state.first.b$b)),
             r_b            = integer(object$Cpar$I),
             naccept_b      = integer(object$Cpar$I),
             pm_w_b         = double(object$K_b[1]),
             pm_mu_b        = double(object$dimb * object$K_b[1]),
             pm_Q_b         = double(LTp * object$K_b[1]),
             pm_Sigma_b     = double(LTp * object$K_b[1]),
             pm_Li_b        = double(LTp * object$K_b[1]),
             sum_Ir_b       = integer(object$Cpar$I * object$K_b[1]),
             hatPr_b_b      = double(object$Cpar$I * object$K_b[1]),
             Pr_b_b         = double(lPr_b),
             hatPr_obs      = double(object$Cpar$I * object$K_b[1]),
             Pr_obs         = double(lPr_b),             
             iter_relabel   = as.integer(0),
             nchange        = integer(l_nchange),
             err            = as.integer(0),
             PACKAGE = thispackage)
  if (MCMC$err) stop("Something went wrong.")             


  ##### New chains for order and rank (corresponding to newly labeled sample)
  ##### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  object$order_b <- matrix(as.numeric(MCMC$chorder_b + 1), ncol=object$K_b[1], byrow=TRUE)
  colnames(object$order_b) <- paste("order", 1:object$K_b[1], sep="")  
  MCMC$chorder_b <- NULL
  
  object$rank_b <- matrix(as.numeric(MCMC$chrank_b+ 1), ncol=object$K_b[1], byrow=TRUE)
  colnames(object$rank_b) <- paste("rank", 1:object$K_b[1], sep="")                          
  MCMC$chrank_b <- NULL
  

  ##### Clustering based on posterior P(alloc = k | y) or on P(alloc = k | theta, b, y) 
  ##### +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (object$K_b[1] == 1){
    object$poster.comp.prob_u <- object$poster.comp.prob_b <- object$poster.comp.prob <- matrix(1, nrow = object$Cpar$I, ncol = 1)
  }else{

    ### Using mean(I(r=k))
    MCMC$sum_Ir_b <- matrix(MCMC$sum_Ir_b, ncol = object$K_b[1], nrow = object$Cpar$I, byrow = TRUE)
    Denom <- apply(MCMC$sum_Ir_b, 1, sum)       ### this should be a vector of length I with all elements equal to the number of saved MCMC iterations 
    object$poster.comp.prob_u <- MCMC$sum_Ir_b / matrix(rep(Denom, object$K_b[1]), ncol = object$K_b[1], nrow = object$Cpar$I)

    ### Using mean(P(r=k | theta, b, y))
    object$poster.comp.prob_b <- matrix(MCMC$hatPr_b_b, ncol = object$K_b[1], nrow = object$Cpar$I, byrow = TRUE)

    ### Using mean(P(r=k | theta, y))
    object$poster.comp.prob <- matrix(MCMC$hatPr_obs, ncol = object$K_b[1], nrow = object$Cpar$I, byrow = TRUE)    
  }  

  ##### Individual sampled values of P(alloc = k | theta, b, y)
  ##### and related quantiles
  ##### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  object$comp.prob_b <- matrix(MCMC$Pr_b_b, ncol = object$K_b[1] * object$Cpar$I, nrow=keepMCMC, byrow=TRUE)
  colnames(object$comp.prob_b) <- paste("P(", rep(1:object$Cpar$I, each=object$K_b[1]), ",", rep(1:object$K_b[1], object$Cpar$I), ")", sep="")
  
  if (length(prob)){
    qq <- apply(object$comp.prob_b, 2, quantile, prob=prob)
    if (length(prob) == 1){
      object$quant.comp.prob_b <- list(matrix(qq, ncol=object$K_b[1], byrow=TRUE))      
    }else{
      object$quant.comp.prob_b <- list()
      for (i in 1:length(prob)){
        object$quant.comp.prob_b[[i]] <- matrix(qq[i,], ncol=object$K_b[1], byrow=TRUE)
      }  
    }  
    names(object$quant.comp.prob_b) <- paste(prob*100, "%", sep="")
  }
  if (!keep.comp.prob) object$comp.prob_b <- NULL
    
  ##### Individual sampled values of P(alloc = k | theta, y), i.e., random effects directly integrated out
  ##### and related quantiles
  ##### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  object$comp.prob <- matrix(MCMC$Pr_obs, ncol = object$K_b[1] * object$Cpar$I, nrow=keepMCMC, byrow=TRUE)
  colnames(object$comp.prob) <- paste("P(", rep(1:object$Cpar$I, each=object$K_b[1]), ",", rep(1:object$K_b[1], object$Cpar$I), ")", sep="")
  
  if (length(prob)){
    qq <- apply(object$comp.prob, 2, quantile, prob=prob)
    if (length(prob) == 1){
      object$quant.comp.prob <- list(matrix(qq, ncol=object$K_b[1], byrow=TRUE))      
    }else{
      object$quant.comp.prob <- list()
      for (i in 1:length(prob)){
        object$quant.comp.prob[[i]] <- matrix(qq[i,], ncol=object$K_b[1], byrow=TRUE)
      }  
    }  
    names(object$quant.comp.prob) <- paste(prob*100, "%", sep="")
  }
  if (!keep.comp.prob) object$comp.prob <- NULL
  
  ##### Posterior means for mixture components
  ##### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  object$poster.mean.w_b <- as.numeric(MCMC$pm_w_b)
  names(object$poster.mean.w_b) <- paste("w", 1:object$K_b[1], sep="")

  object$poster.mean.mu_b <- matrix(MCMC$pm_mu_b, nrow=object$K_b[1], ncol=object$dimb, byrow=TRUE)
  rownames(object$poster.mean.mu_b) <- paste("j", 1:object$K_b[1], sep="")
  colnames(object$poster.mean.mu_b) <- paste("m", 1:object$dimb, sep="")

  object$poster.mean.Q_b <- object$poster.mean.Sigma_b <- object$poster.mean.Li_b <- list()
  for (j in 1:object$K_b[1]){
    tmpQ <- matrix(0, nrow=object$dimb, ncol=object$dimb)
    tmpQ[lower.tri(tmpQ, diag=TRUE)] <- MCMC$pm_Q_b[((j-1)*LTp+1):(j*LTp)]
    tmpQ[upper.tri(tmpQ, diag=FALSE)] <- t(tmpQ)[upper.tri(t(tmpQ), diag=FALSE)]
    object$poster.mean.Q_b[[j]] <- tmpQ
    
    tmpSigma <- matrix(0, nrow=object$dimb, ncol=object$dimb)
    tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- MCMC$pm_Sigma_b[((j-1)*LTp+1):(j*LTp)]
    tmpSigma[upper.tri(tmpSigma, diag=FALSE)] <- t(tmpSigma)[upper.tri(t(tmpSigma), diag=FALSE)]
    object$poster.mean.Sigma_b[[j]] <- tmpSigma
    
    tmpLi <- matrix(0, nrow=object$dimb, ncol=object$dimb)
    tmpLi[lower.tri(tmpLi, diag=TRUE)] <- MCMC$pm_Li_b[((j-1)*LTp+1):(j*LTp)]
    object$poster.mean.Li_b[[j]] <- tmpLi      
  }
  names(object$poster.mean.Q_b) <- names(object$poster.mean.Sigma_b) <- names(object$poster.mean.Li_b) <- paste("j", 1:object$K_b[1], sep="")
  
  return(object)  
}
