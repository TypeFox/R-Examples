#' @name writeJagsFormula
#' @export
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats update
#' @importFrom utils tail
#'  
#' @title Write the JAGS Formula for a Hyde Node
#' @description Based on the information provided about the node,
#'   an appropriate JAGS model is written in text.  This is combined with 
#'   the other node models to generate the complete network.
#'   
#' @param fit a model object
#' @param nodes a vector of node names, usually passed from \code{network$nodes}
#' @param ... Additional arguments to be passed to other methods
#' 
#' @details Methods for different model objects can be written so that 
#'   this function can be extended as desired.
#'   
#'   The resulting formulas are based on the coefficient matrix of the fitted
#'   model, and the returned result is the JAGS code representing the 
#'   regression equation of the model.
#'   
#'   In the \code{writeJagsFormula.glm} method, appropriate transformations
#'   exist for the following combinations:
#'   \enumerate{
#'     \item family = binomial; link = logit
#'     \item family = poisson; link = log
#'     \item family = gaussian; link = identity (calls \code{writeJagsFormula.lm})
#'   }
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @seealso \code{\link{writeJagsModel}}, \code{\link{writeNetworkModel}}
#' 
#' @examples
#' data(PE, package="HydeNet")
#' fit <- lm(d.dimer ~ pregnant + pe, data=PE)
#' writeJagsFormula(fit, nodes=c("d.dimer", "pregnant", "pe"))
#' 
#' fit.glm <- glm(death ~ pe + treat, data=PE, family="binomial")
#' writeJagsFormula(fit.glm, nodes=c("death", "pe", "treat"))

writeJagsFormula <- function(fit, nodes, ...) UseMethod("writeJagsFormula")

#' @rdname writeJagsFormula
#' @export
#' 

writeJagsFormula.cpt <- function(fit, nodes, ...)
{
  form <- paste0(utils::tail(names(dimnames(fit)), 1),
                 " ~ ",
                 paste0(names(dimnames(fit))[-length(names(dimnames(fit)))],
                        collapse = " + "))
  rToJags(stats::as.formula(form)) 
}

#' @rdname writeJagsFormula
#' @export
#' @importFrom stringr str_trim
#' 

writeJagsFormula.glm <- function(fit, nodes, ...){
  if (fit$family$family == "gaussian" & fit$family$link == "identity") 
    return(writeJagsFormula.lm(fit))
  
  mdl <- pixiedust::dust(fit, descriptors = c("term", "term_plain", "level")) %>%
    as.data.frame(sprinkled = FALSE)
  
  regex <- factorRegex(fit)
  
  mdl <- makeJagsReady(mdl, regex, nodes) %>%
    mutate(term_plain = gsub(":", "*", term_plain))
  
  #* rhs = right hand side
  rhs <- paste(round(mdl$estimate, getOption("Hyde_maxDigits")), 
               ifelse(is.na(mdl$term_plain), "", "*"),
               ifelse(is.na(mdl$term_plain), "", mdl$term_plain), 
               collapse=" + ")
  
  #* Binomial Proportion
  if (fit$family$family == "binomial" & fit$family$link == "logit"){
    rhs <- paste0("ilogit(", rhs, ")")
  }
  
  #* Poisson Regression
  if (fit$family$family == "poisson" & fit$family$link == "log"){
    rhs <- paste("exp(", rhs, ")")  
  }
  
  out_fm <- paste0(as.character(fit$call$formula)[2], " ~ ", rhs)
  rToJags(stats::as.formula(out_fm)) 
}

#' @rdname writeJagsFormula
#' @export
#' @importFrom stringr str_trim
#' 

writeJagsFormula.lm <- function(fit, nodes, ...){
  mdl <- pixiedust::dust(fit, descriptors = c("term", "term_plain", "level")) %>%
    as.data.frame(sprinkled = FALSE)
  
  regex <- factorRegex(fit)
  
  mdl <- makeJagsReady(mdl, regex, nodes) %>%
    mutate(term_plain = gsub(":", "*", term_plain))
  
  #* rhs = right hand side
  rhs <- paste(round(mdl$estimate, getOption("Hyde_maxDigits")), 
               ifelse(is.na(mdl$term_plain), "", "*"),
               ifelse(is.na(mdl$term_plain), "", mdl$term_plain), 
               collapse=" + ")
  
  out_fm <- paste0(as.character(fit$call$formula)[2], " ~ ", rhs)
  rToJags(stats::as.formula(out_fm)) 
}

#' @rdname writeJagsFormula
#' @export
#' @import nnet
#'

writeJagsFormula.multinom <- function(fit, nodes, ...){
#   mdl <- broom::tidy(fit, exponentiate=FALSE)[, c("y.level", "term", "estimate")] 
#   
#   regex <- factorRegex(fit)
#   
#   mdl <- makeJagsReady(mdl, regex)
#   mdl <- dplyr::arrange(mdl, y.level, term_name)
#   
#   right_side <- function(l, m=mdl)
#   {
#     m <- m[m$y.level == l, ]
#     paste(round(m$estimate, getOption("Hyde_maxDigits")), 
#           ifelse(m$jagsVar == "(Intercept)", "", "*"),
#           ifelse(m$jagsVar == "(Intercept)", "", m$jagsVar), 
#           collapse=" + ")
#   }
#   
#   sapply(unique(as.character(mdl$y.level)), right_side)
         
  
  
  if (is.null(fit$model)) fit <- stats::update(fit, model=TRUE)
  fm <- as.character(fit$call$formula)
  out_fm <- paste0("pi.", fm[2])
  fm <- stringr::str_trim(unlist(strsplit(fm[-(1:2)], "[+]")))
  
  fm
  
  fm <- unlist(sapply(fm, function(x){
    if (! x %in% names(attributes(fit$terms)$dataClasses)){
      return(NULL)
    }
    if (attributes(fit$terms)$dataClasses[x] == "factor"){
      return(paste0("(", x, " == ", 2:nlevels(fit$model[, x]), ")"))
    }
    else return(x)
  }))
  
  fm <- lapply(1:nrow(stats::coef(fit)), 
               function(r){
                 if (is.null(fm)) stats::coef(fit)[r, 1]
                 else paste0(round(stats::coef(fit)[r, 1], getOption("Hyde_maxDigits")),
                             " + ", 
                             paste(round(stats::coef(fit)[r, -1], getOption("Hyde_maxDigits")), 
                                   fm, sep="*", collapse=" + "))})
  
  fm <- sapply(fm, function(x) paste0("exp(", x, ") / (1 + ", 
                                      paste(sapply(fm, function(x) paste0("exp(", x, ")")), collapse=" + "),
                                      ")"))
  fm <- c(paste0("1 - (", paste(fm, collapse=" + "), ")"), fm)
  fm <- paste0(out_fm, "[", 1:length(fm), "] <- ", fm, collapse="; ")
  
  return(fm) 
}

#' @rdname writeJagsFormula
#' @export

writeJagsFormula.xtabs <- function(fit, ...){
  fm <- attributes(fit)$call$formula
  out_fm <- paste(fm[2], fm[1])
  
  pi <- fit/sum(fit)
  names(pi) <- 1:length(pi)
  pi <- paste0("pi.", 
               fm[2], 
               "[", 
               names(pi), 
               "] <- ", 
               round(pi, getOption("Hyde_maxDigits")), 
               collapse="; ")
  
  return(pi)
}

utils::globalVariables(c("term_plain"))
