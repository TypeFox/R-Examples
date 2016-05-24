#' Pick-A-Point Technique
#'
#' Probe moderation effect using the Pick-A-Point technique
#'
#' @aliases spotlight
#' @param model Regression model (lm, glm, list).
#' @param dv Dependent variable (character).
#' @param iv Independent variable (character).
#' @param mod Moderator variable(s) (character or character vector).
#' @param points List of points to test for each moderator variable (list).
#' @param method Method to use. Possible values are: \code{"meansd", "percentiles"}, \code{method="meansd"} by default.
#' @param alpha Alpha level to use (numeric).
#' @param yas Show y (or conditional effect) as: \code{"none", "ratio","probability","percentage"}, \code{yas="none"} by default.
#'
#' @return A list with the elements
#'
#' @examples
#' \dontrun{
#' myModel <- lm('dv ~ iv + mod', data=someData)
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD')
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD', alpha=.01)
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD', method='percentiles')
#' pickapoint(myModel, dv='DV', iv='IV', mod='MOD', points=c(1,2,3))
#' }
#' @rdname pickapoint
#'
#' @references
#' Spiller, S. A., Fitzsimons, G. J., Lynch, J. G., Jr, & McClelland, G. H. (2013). Spotlights, floodlights, and the magic number zero: Simple effects tests in moderated regression. Journal of Marketing Research, 50(2), 277-288.
#'
#' Aiken, L. S., & West, S. G. (1991). Multiple regression: Testing and interpreting interactions. Thousand Oaks, CA: Sage Publications.
#'
#' @export

pickapoint <- function(model,dv,iv,mod,points,method,alpha, yas) UseMethod("pickapoint")
pickapoint <- spotlight <- function(model,dv,iv,mod,points,method="meansd",alpha=.05, yas='none'){
  qrange <- c(.10,.25,.50,.75,.90)
  papret <- list()

  #run checks on params
  if(!is.character(dv) || !is.character(iv) || !is.character(mod)){
    stop('params dv/iv/mod incorrectly specified')
  }
  if(!missing(points)) {
    if(!is(points,'list')){
      stop('param points has to be a list')
    }
  } else{
    points <- list()
  }
  if(!(method == 'meansd' || method == 'percentiles')) {
    stop('param method has to be either "meansd" or "percentiles"')
  }
  if(!is.numeric(alpha)){
    stop('param alpha has to be numeric')
  }
  if(is(model,'list')){
    if(is.null(model$coefficients)){
      stop("list model does not contain 'coefficients'")
    } else{
      beta.hat <- model$coefficients
    }
    if(is.null(model$vcov)){
      stop("list model does not contain 'vcov'")
    } else{
      cov <- model$vcov
    }
  } else if(is(model,'glm')){
    if(model$family$link == 'log' | model$family$link == 'logit'){
      beta.hat <- coef(model)
      cov <- vcov(model)
      papret$link <- model$family$link
    } else{
      stop('this function currently supports only glm models with link == log or logit')
    }
  } else if(is(model,'lm')){
    beta.hat <- coef(model)
    cov <- vcov(model)
  } else {
    stop('this function currently supports only the lm, glm(link == log or logit) & custom list objects')
  }

  if(missing(yas)){
    papret$yas <- 'none'
  } else {
    if(yas == 'ratio'){
      papret$yas <- yas
      papret$sline <- 1
    } else if(yas == 'prob' | yas == 'probability'){
      papret$yas <- 'prob'
      papret$sline <- 0.5
    } else if(yas == 'percent' | yas == 'percentage'){
      papret$yas <- 'percent'
      papret$sline <- 0
    } else {
      stop("only 'ratio' or 'prob' ('probability') or 'percent' ('percentage') is supported in the current version")
    }
  }
  yasfunc <- function(pap,value){
    if(pap$yas == 'ratio') return(exp(value))
    else if(pap$yas == 'prob') return(exp(value) / (1 + exp(value)))
    else if(pap$yas == 'percent') return(100*(exp(value) - 1))
    else return(value)
  }

  data <- model$model
  papret$method <- method
  papret$iv <- iv
  papret$dv <- dv
  papret$mod <- mod

  df <- nrow(data) - (length(beta.hat) + 1)
  tcrit <- qt(p=alpha/2,df=df,lower.tail = FALSE)
  modpprefvals <- vector('list')
  modppmatrix <- vector('list')

  interactionterms <- vector()
  for (i in 1:length(mod)){
    if(length(points[[mod[i]]]) > 0) {
      modpprefvals[[mod[i]]] <- points[[mod[i]]]
    } else if(length(unique(data[[mod[i]]])) == 2){
      #binary treatment for binary variables
      modpprefvals[[mod[i]]] <- unique(data[[mod[i]]])
    } else if (method == 'percentiles'){
      #calculate values for different percentiles
      modpprefvals[[mod[i]]] <- quantile(data[[mod[i]]], qrange)
    } else{
      #calculate values for mean +/-1sd
      modpprefvals[[mod[i]]] <- c(mean(data[[mod[i]]]) - sd(data[[mod[i]]]), mean(data[[mod[i]]]), mean(data[[mod[i]]]) + sd(data[[mod[i]]]) )
    }
    modppmatrix[[mod[i]]] <- c(1:length(modpprefvals[[mod[i]]]))
    interactionterms <- c(interactionterms, paste(papret$iv,papret$mod[i],sep=':'))
  }
  modpp <- expand.grid(modppmatrix)
  modvals <- vector()
  for(i in 1:ncol(modpp)){
    modvals <- cbind(modvals,modpprefvals[[i]][modpp[,i]])
  }

  papret$outcome <- vector()
  for(i in 1:nrow(modvals)){
    cy <- beta.hat[papret$iv]
    cyse <- cov[papret$iv, papret$iv]
    for(j in 1:ncol(modvals)) {
      cy <- cy + (beta.hat[interactionterms[j]] * modvals[i,j])
      cyse <- cyse + (modvals[i,j]^2*cov[interactionterms[j], interactionterms[j]] + 2*modvals[i,j]*cov[papret$iv, interactionterms[j]])
    }
    #calculate se and ci
    cyse <- sqrt(cyse)
    culci <- cy + tcrit*cyse
    cllci <- cy - tcrit*cyse
    ct <- cy/cyse

    #recalculate if yas func
    if(papret$yas != 'none') {
      cy <- yasfunc(papret, cy)
      cyse <- sqrt(cy^2 * cyse)
      culci <- yasfunc(papret, culci)
      cllci <- yasfunc(papret, cllci)
    }

    #cyse <- yasfunc(papret$yas, cyse)
    cp <- 2*pt(-abs(ct),df=df)
    papret$outcome <- rbind(papret$outcome,c(modvals[i,],cy,cyse,ct,cp,cllci,culci))
  }

  colnames(papret$outcome) <- c(mod,'Effect','SE','t','p','llci','ulci')
  rownames(papret$outcome) <- rep('',nrow(papret$outcome))
  papret$call <- match.call()
  class(papret) <- "pickapoint"
  papret
}

