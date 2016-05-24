#' Replicate meta-analysis results and summaries from MetaWin 2.0.
#'
#' Replicate meta-analysis results and summaries from Rosenberg's et al. (2000)
#' software 'MetaWin' 2.0.  Currently only replicates moderator analyses and not
#' meta-regressions.  
#'
#' @param model A two-sided linear formula object describing the model, with the
#'    response (effect sizes) on the left of a ~ operator and the moderator 
#'    variables, separated by +, :, * operators, on the right.  NOTE:  
#'    MetaWin was limited to analyses with a single moderator variable.  
#'    This function currently supports only categorical moderators.
#' @param weights A vector of effect size variances that will be used as 
#'    weights for the meta-analysis.
#' @param effects_model The default is \code{"random"}, which specifies a 
#'    random-effects meta-analysis.  Other options include \code{"fixed"} 
#'    which presents fixed-effect analyses.
#' @param data An optional data frame containing the variables named in model
#'    and weights.   
#' @param bootstraps The number of bootstraps used to estimate confidence 
#'    intervals.  As with 'MetaWin' 2.0, the default is 999.
#'
#' @return NULL
#'
#' @references Rosenberg, M.S., Adams, D.C., and Gurevitch, J. 2000. 
#'    MetaWin: Statistical Software for Meta-Analysis. Sinauer Associates 
#'    Sunderland, Massachusetts. 
#'
#' @importFrom metafor rma
#' @importFrom stats anova lm as.formula pchisq qt
#' @export replicate_MetaWin2.0

replicate_MetaWin2.0 <- function(model, 
                                 weights, 
                                 effects_model = "random", 
                                 data, 
                                 bootstraps = 999) {
  
  cat("\n=== START of Rosenberg et al. (2002) MetaWin 2.0 output ===\n\n")
  
  effects_model <- switch(effects_model,
                          random = "DL",
                          fixed = "FE"
                    )
 
  if(length(all.vars(model)) == 1) {
    rma_pooled <- rma(model, vi = weights, data = data, method = effects_model)
    effect.pooled <- anova(do.call("lm", 
                  list(model, data, NULL, (1.0/(weights + rma_pooled$tau2)))))
    cat(paste0("\nEstimate of pooled variance: ", 
               round(rma_pooled$tau2, 6), "\n"))
    SUMMARY.RESULTS(rma_pooled, effect.pooled)
    
    
    cat("\n=== END of Rosenberg et al. (2002) MetaWin 2.0 output ===\n\n")
    return(NULL)
  }
  
  # get tau and pooled effects from model
  fixedModel <- as.formula(paste(all.vars(model)[1], 
                                 "~", 
                                 all.vars(model)[2], 
                                 "-1"))
  rma_results <- rma(fixedModel, vi = weights, data = data, method = effects_model)
  rma_pooled <- rma(as.formula(paste(all.vars(model)[1], "~ 1")), 
                    tau2 = rma_results$tau2, vi = weights, data = data, 
                    method = effects_model)

  # get model sums of squares from lm based on metafor's tau square
  effects.results <- anova(do.call("lm", list(model, 
                                              data, 
                                              NULL, 
                                              (1.0/(weights + rma_results$tau2)))))
  effect.pooled <- anova(do.call("lm", list(as.formula(paste(all.vars(model)[1], "~ 1")), 
                                            data, 
                                            NULL, 
                                            (1.0/(weights + rma_results$tau2)))))
  
  cat(paste0("\nEstimate of pooled variance: ", round(rma_results$tau2, 6), "\n"))
  CATEGORICAL.RESULTS(rma_results, effects.results, data)
  SUMMARY.RESULTS(rma_pooled, effect.pooled)
  
  cat("\n\n=== END of Rosenberg et al. (2002) MetaWin 2.0 output ===\n")
  return(NULL)
}


# a few replicate_MetaWin2.0 helpers
SUMMARY.RESULTS <- function (rma.results, effects.results) {
  cat("\nSUMMARY RESULTS\n\n")
  cat(sprintf(c("%-20s", "%-5s", "%-10s\n"), c("Heterogeneity", "df", "Prob(Chi-Square)")))
  cat("------------------------------------------------\n")
  cat(sprintf(c("%-20s", "%-5s", "%-10s\n\n"), 
              c(
                paste0("Qtotal  ", round(effects.results[1,2], 6)), 
                effects.results[1,1], 
                round(1.0 - pchisq(effects.results[1,2], df=effects.results[1,1]), 5))
  )
  )
  cat(sprintf(c("%-23s", "%-23s", "%-23s", "%-23s\n"), c("Mean Effect Size", "95% CI", "Bootstrap CI", "Bias CI")))
  cat("--------------------------------------------------------------------------------------------\n")
  cat(sprintf(c("%-23s", "%-23s", "%-23s", "%-23s\n\n"), 
              c(
                paste0("E++  ", round(rma.results$b, 6)), 
                paste0(round(t_lCI(rma.results$b, rma.results$se, effects.results[1,1]), 6), " to ", round(t_uCI(rma.results$b, rma.results$se, effects.results[1,1]), 6)),
                paste0(round(rma.results$ci.lb, 6), " to ", round(rma.results$ci.ub, 6)),
                paste0(round(rma.results$ci.lb, 6), " to ", round(rma.results$ci.ub, 6))
              )
  ))
  cat(paste0("Sqrt Pooled Variance = ", round(sqrt(rma.results$tau2), 6), "\n"))
  cat(sprintf(c("%-25s", "%-20s\n"), 
              c(
                paste0(" Mean Study Variance = ", round(mean(rma.results$vi), 6)),
                paste0(" Ratio = ", round((sqrt(rma.results$tau2))/mean(rma.results$vi), 6))
              )
  )
  )
  cat("\n------------------------------------------------\n")
}


CATEGORICAL.RESULTS <- function (rma.results, effects.results, theData) {
  cat("\nCATEGORICAL RESULTS\n\n")
  cat("--Heterogeneity--\n")
  cat(sprintf(c("%-10s", "%-10s", "%-10s\n"), c("Class", "#Studies", "PooledVar")))
  cat("-------------------------------\n")
  cat("MJL NOTE: I don't know how 'PooledVar' is calculated...\n\n\n\n")
  
  
  cat(sprintf(c("%-10s", "%-10s","%-15s","%-20s", "%-20s\n"), c("Model", "df", "Q", "Prob(Chi-Square)", "Prob(Rand)")))
  cat("---------------------------------------------------------------------\n")
  effectsRange <- nrow(effects.results) 
  cat(sprintf(c("%-10s", "%-10s","%-15s","%-20s", "%-20s\n"), 
              c("Between",
                effects.results[1,1],
                round(effects.results[1,2], 6),
                round(1.0 - pchisq(effects.results[1,2], df=effects.results[1,1]), 5),
                round(1.0 - pchisq(effects.results[1,2], df=effects.results[1,1]), 5)
              )
  )
  )
  cat(sprintf(c("%-10s", "%-10s","%-15s","%-15s", "%-10s\n"), 
              c("Within",
                effects.results[2,1],
                round(effects.results[2,2], 6),
                round(1.0 - pchisq(effects.results[2,2], df=effects.results[2,1]), 5),
                ""
              )
  )
  )
  cat("---------------------------------------------------------------------\n")
  cat(sprintf(c("%-10s", "%-10s","%-15s","%-15s", "%-10s\n"), 
              c("Total",
                sum(effects.results[1:effectsRange,1]),
                round(sum(effects.results[1:effectsRange,2]), 6),
                round(1.0 - pchisq(sum(effects.results[1:effectsRange,2]), df=sum(effects.results[1:effectsRange,1])), 5),
                ""
              )
  )
  ) 
  
  cat("\n\n--Mean Effect Sizes--\n")
  cat(sprintf(c("%-10s", "%-10s","%-10s","%-5s", "%-25s", "%-25s", "%-25s\n"), c("Class", "#Studies", "E+", "df", "95% CI", "Bootstrap CI",  "Bias CI")))
  cat("----------------------------------------------------------------------------------------------------------------\n")
  
  for(i in 1:length(levels(theData$Tree))) {
    cat(sprintf(c("%-10s", "%-10s","%-10s","%-5s", "%-25s", "%-25s", "%-25s\n"), 
                c(
                  levels(theData$Tree)[i], 
                  table(theData$Tree)[i],
                  round(rma.results$b[i], 6),
                  table(theData$Tree)[i] - 1,
                  paste0(round(t_lCI(rma.results$b[i], rma.results$se[i], table(theData$Tree)[i] - 1), 6), " to ", round(t_uCI(rma.results$b[i], rma.results$se[i], table(theData$Tree)[i] - 1), 6)),
                  paste0(round(rma.results$ci.lb, 6)[i], " to ", round(rma.results$ci.ub, 6)[i]),
                  paste0(round(rma.results$ci.lb, 6)[i], " to ", round(rma.results$ci.ub, 6)[i])
                )
    ))
  }
  
}

t_lCI <- function(E, var, df) {
  return(E - qt(1.0 - 0.024991, df) * var)
}

t_uCI <- function(E, var, df) {
  return(E + qt(1.0 - 0.024991, df) * var)
}
