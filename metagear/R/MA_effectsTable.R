#' Generate an ANOVA-like effects table for a meta-analysis.
#'
#' Generates an ANOVA-like effects table that summarizes the within and between-study 
#' homogeneity tests (Q-tests), as well as moderator level Q-tests as originally
#' described by Hedges and Olkin (1985; p. 156).  
#'
#' @param model A two-sided linear formula object describing the model, with the
#'    response (effect sizes) on the left of a ~ operator and the moderator 
#'    variables, separated by +, :, * operators, on the right.  
#' @param dataFrame An optional data frame containing the variables named in the
#'    model. 
#' @param weights A column label from data.frame of variances to be used as weights.
#' @param tau2_model Specifies the between-study variance estimator.  Default 
#'    is \code{"DL"}.
#'
#' @return NULL
#' 
#' @references Hedges, L.V., and I. Olkin. 1985. Statistical methods for 
#'    meta-analysis. Academic Press, New York, USA.
#'
#' @importFrom metafor rma
#' @importFrom stats anova lm pchisq
#' @export MA_effectsTable

MA_effectsTable <- function(model,
                            dataFrame,
                            weights,
                            tau2_model = "DL") {

  # get tau from metafor
  rma.results <- rma(model, 
                     vi = weights, 
                     data = dataFrame, 
                     method = tau2_model)

  # get model sums of squares from lm based on metafor's tau
  effects.results <- anova(lm(model, 
                              weights = 1.0/(weights + rma.results$tau2),
                              data = dataFrame))

  # get summary of the Overall Model
  printModelSummary(effects.results)
  printEffectTestsSummary(effects.results)
  
  return (NULL)
}

printModelSummary <- function(ANOVA) {
  effectsRange <- nrow(ANOVA) - 1 #ANOVA used to be effects.results
  model.summary <- data.frame(
    SOURCE = c("model", "residual error", "total"), 
    Q = c(
      sum(ANOVA[1:effectsRange, 2]),
      ANOVA[effectsRange + 1, 2],
      sum(ANOVA[1:effectsRange, 2]) + ANOVA[effectsRange + 1, 2]
    ), 
    DF = c(
      sum(ANOVA[1:effectsRange, 1]),
      ANOVA[effectsRange + 1, 1],
      sum(ANOVA[1:effectsRange, 1]) + ANOVA[effectsRange + 1, 1]
    ), 
    P = c(
      1.0 - pchisq(sum(ANOVA[1:effectsRange, 2]), 
                   df = sum(ANOVA[1:effectsRange, 1])),
      1.0 - pchisq(ANOVA[effectsRange + 1,2], 
                   df = ANOVA[effectsRange + 1,1]),
      1.0 - pchisq(sum(ANOVA[1:effectsRange, 2]) + ANOVA[effectsRange + 1, 2], 
                   df = sum(ANOVA[1:effectsRange, 1]) + ANOVA[effectsRange + 1, 1])
    )
  )
  print(model.summary, row.names = FALSE, digits = 5)
}

# get summary of the Effect Tests
printEffectTestsSummary <- function(ANOVA) {
  effectsRange <- nrow(ANOVA) - 1
  mainEffects <- rownames(ANOVA[1:effectsRange, ])
  model.summary <- data.frame(
    SOURCE = mainEffects, 
    Q = ANOVA[1:effectsRange, 2], 
    DF = ANOVA[1:effectsRange, 1], 
    P = 1.0 - pchisq(ANOVA[1:effectsRange, 2], ANOVA[1:effectsRange, 1])
  )
  print(model.summary, row.names = FALSE, digits = 5)
}

