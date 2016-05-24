#' informed Genome-Wide Association Study
#'
#' Perform an informed Genome-Wide Association Study (iGWAS). This is based on a current study and a prior study. The goal is to discover significant SNPs in the current study using hypothesis testing. The prior study is used to improve power. P-values and sample sizes are used for both the current and the prior studies.
#'
#' This method computes the p-value weights based on the prior p-values, and uses them in multiple testing (p-value weighting) for the current p-values.  The p-value weighting method (e.g. Unweighted, Bayes) and the multiple testing adjustment (e.g. Bonferroni, Benjamini-Hochberg) can be specified independently.
#'
#'  For more details, see the paper "Optimal Multiple Testing Under a
#'  Gaussian Prior on the Effect Sizes", by Dobriban, Fortney, Kim and Owen,
#'   \url{http://arxiv.org/abs/1504.02935}
#'
#'@param P_current P-values in the current study, a numeric vector of length J, with entries between 0 and 1
#'@param N_current sample size in the current study, a positive integer (or vector of length J)
#'@param P_prior P-values in the prior study, a numeric vector of length J, with entries between 0 and 1
#'@param N_prior sample size in the current study, a positive integer (or vector of length J)
#' @param q (optional) uncorrected level at which tests should be performed. Default \code{q = 0.05}
#' @param weighting_method (optional) weighting method used. Available methods: \code{c("unweighted",
#' "bayes", "spjotvoll", "exponential")}. The default is \code{"bayes"}.
#' @param p_adjust_method (optional) adjustment method for multiple testing used. The available methods are
#' \code{"genome-wide"} and those from the \code{p.adjust} function in the \code{stats} package.:
#' \code{c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "fdr", "none")}.
#' \code{"genome-wide"} tests all hypotheses at the genome-wide level \code{5*10^-8}.
#' The default is \code{"genome-wide"}.
#' @param sides (optional) The prior p-values must be one or two-sided: sides = 1 or 2. Default \code{sides = 2}
#' @param phi (optional) dispersion factor used to multiply all standard errors. Default \code{phi = 1}.
#' Used only for Bayes weights.
#' @param  beta (optional) weights are proportional to \code{exp(mu*beta)}, default \code{beta=2}.
#' Used only for Exponential weights.
#' @param  UB_exp (optional) upper bound on the weights (default \code{UB = Inf}).
#'  Used only for Exponential weights.
#' @param  figure (optional) \code{figure = "TRUE"} creates a manhattan plot of the weighted and unweighted p-values.
#' Possible values \code{c("TRUE","FALSE")}, default \code{"FALSE"}
#' @param  GWAS_data_frame (optional) when \code{figure = "TRUE"}, this is the parameter used to create a
#' Manattan plot. it must be a data frame with columns containing \code{c("CHR","BP")}.
#' These parameters are passed to the qqman package for plotting manhattan plots. Default is NA.
#'
#'
#'@return A list containing:
#'
#'\code{sig_ind}: A vector of 0-1s indicating the significant tests (1-s)
#'
#'\code{num_sig}: The number of significant tests. Equals \code{sum(sig_ind)}
#'
#'\code{w}: The computed p-value weights
#'
#'\code{P_w}: The weighted p-values. These equal \code{P_current/w}
#'
#' @family p-value weighting
#'@export
#'
iGWAS <-
  function(P_current, N_current, P_prior, N_prior, q = 0.05, weighting_method = "bayes",
           p_adjust_method = "genome-wide", sides = 2, phi = 1, beta = 2, UB_exp = Inf,
           figure = "FALSE", GWAS_data_frame=NA) {
    # Error checking: stop if the variables are not in range
    {
      if (any(P_current > 1) |
          any(P_current < 0) | any(P_prior > 1) | any(P_prior < 0)) {
        stop("P-values must be between 0 and 1")
      }
      if (any(N_current < 1) | any(N_prior < 1)) {
        stop("Sample sizes must be at least 1")
      }
      if ((q <= 0) |
          (q >= 1)) {
        stop("Level q at which tests will be performed must be in (0,1)")
      }
      if (phi < 0) {
        stop("Dispersion parameter phi for Bayes weights must be non-negative")
      }
      if (UB_exp <= 1) {
        stop("Upper bound UB for exponential weights must be greater than 1")
      }
      if (!((sides == 1) | (sides == 2))) {
        stop("The prior p-values must be one or two-sided: sides = 1 or 2")
      }
    }
    #this is needed to avoid the NOTE: iGWAS: no visible binding for global variable ???qqman???
    qqman = NULL
    # Define auxiliary variables to compute weights
    {
      J <- length(P_current)
      if (length(N_current) == 1) {
        N_current <- N_current * rep(1,J)
      }
      if (length(N_prior) == 1) {
        N_prior <- N_prior * rep(1,J)
      }
      T_prior <- qnorm(P_prior / sides)
      mu <- T_prior * sqrt(N_current / N_prior)
      sigma <- phi * sqrt(N_current / N_prior)
    }

    # Compute weights: switch according to 'weighting_method'
    {
      w_methods = c("unweighted","bayes","spjotvoll","exponential")
      #for genome-wide significance level, need to adjust q
      if (p_adjust_method == "genome-wide") {
        sig_level = 5 * 10 ^(-8)
        q = sig_level*J
      }
      switch(
        weighting_method,
        unweighted = {
          w <- rep(1,J)
        },
        bayes = {
          #note: the weight functions use q/J instead of q
          if (phi>0) {
          res <- bayes_weights(mu, sigma, q / J)
          w <- res$w
          } else  { #phi = 0
            epsilon = 1e-3;
            w <- spjotvoll_weights(-abs(mu)-epsilon, q / J)
          }
        },
        spjotvoll = {
          w <- spjotvoll_weights(mu, q / J)
        },
        exponential = {
          w <- exp_weights(mu, beta, UB_exp)
        },
        {
          cat(c("Available methods:", methods))
          stop("Method must be one of the available methods")
        }
      )
    }

    #Perform weighted multiple testing
    {
      P_weighted <- P_current / w
      #if the adjustment method is genome-wide, compare all p-values against Genome-wide significance
      #threshold 5*10^{-8}
      if (p_adjust_method == "genome-wide") {
        sig_level = 5 * 10 ^(-8)
        sig_ind <- (P_weighted < sig_level)
      } else {
      #if the adjustment method is something else, run p.adjust
        P_w_adjusted <- p.adjust(P_weighted, p_adjust_method)
        sig_ind <- (P_w_adjusted < q)
      }
      cat(
        c(
          "Number of significant tests using", weighting_method, "weights and", p_adjust_method, "correction:" , sum(sig_ind), "\n"
        )
      )
    }

    #plot
    if (figure=="TRUE"){
      #library(qqman)
      requireNamespace(qqman, quietly = TRUE)
      GWAS_data_frame$CHR <- as.numeric(GWAS_data_frame$CHR)
      GWAS_data_frame$BP <- as.numeric(GWAS_data_frame$BP)
      GWAS_data_frame$P <- P_current
      GWAS_data_frame_wt <- GWAS_data_frame
      GWAS_data_frame_wt$P <- P_weighted

      suppressWarnings(qqman::manhattan(GWAS_data_frame_wt, ylim=c(0,29), col = c("chartreuse", "chartreuse"), suggestiveline=F))
      par(new=T)
      suppressWarnings(qqman::manhattan(GWAS_data_frame, ylim=c(0,29), col = c("black", "black"), suggestiveline=F))
    }
    results <- list(
      "w" = w, "P_w" = P_weighted, "sig_ind" = sig_ind, "num_sig" = sum(sig_ind)
    )
    return(results)
  }
