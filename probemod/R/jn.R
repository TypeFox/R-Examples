#' Johnson-Neyman Technique
#'
#' Probe moderation effect using the Johnson-Neyman technique
#'
#' @aliases floodlight
#' @param model Regression model (lm, glm, list).
#' @param dv Dependent variable (character).
#' @param iv Independent variable (character).
#' @param mod Moderator variable(s) (character or character vector).
#' @param mrange Range of values that jn should examine for moderator variable. Uses the current range of moderator values by default (numeric vector).
#' @param alpha Alpha level to use (numeric).
#' @param yas Show y (or conditional effect) as: \code{"none", "ratio","probability","percentage"}, \code{yas="none"} by default.
#'
#' @return A list with the elements
#'
#' @examples
#' \dontrun{
#' myModel <- lm('DV ~ IV + MOD', data=someData)
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD')
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD', alpha=.01)
#' plot(jnresults)
#' }
#' @rdname jn
#'
#' @references
#' Spiller, S. A., Fitzsimons, G. J., Lynch, J. G., Jr, & McClelland, G. H. (2013). Spotlights, floodlights, and the magic number zero: Simple effects tests in moderated regression. Journal of Marketing Research, 50(2), 277-288.
#'
#' Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and multilevel regression: Inferential and graphical techniques. Multivariate Behavioral Research, 40(3), 373-400.
#'
#' @export

jn <- function(model,dv,iv,mod,mrange,alpha, yas) UseMethod("jn")
jn <- floodlight <- function(model,dv,iv,mod,mrange,alpha =.05, yas='none'){
  jnret <- list()

  #run checks on params
  if(!is.character(dv) || !is.character(iv) || !is.character(mod)){
    stop('params dv/iv/mod incorrectly specified')
  }
  if(!is.numeric(alpha)){
    stop('param alpha has to be numeric')
  }
  if(is(model,'list')){
    if(is.null(model$coefficients)){
      stop("list model does not contain 'coefficients'")
    } else {
      beta.hat <- model$coefficients
    }
    if(is.null(model$vcov)){
      stop("list model does not contain 'vcov'")
    } else {
      cov <- model$vcov
    }
  } else if(is(model,'glm')){
    if(model$family$link == 'log' | model$family$link == 'logit'){
      beta.hat <- coef(model)
      cov <- vcov(model)
      jnret$link <- model$family$link
    } else {
      stop('this function currently supports only glm models with link == log or logit')
    }
  } else if(is(model,'lm')){
    beta.hat <- coef(model)
    cov <- vcov(model)
  } else{
    stop('this method currently supports only the lm, glm(yas == log or logit) & custom list objects')
  }

  jnret$plot <- list()
  jnret$plot$sline <- 0
  if(missing(yas)){
    jnret$yas <- 'none'
  } else {
    if(yas == 'ratio'){
      jnret$yas <- yas
      jnret$plot$sline <- 1
    } else if(yas == 'prob' | yas == 'probability'){
        jnret$yas <- 'prob'
        jnret$plot$sline <- 0.5
    } else if(yas == 'percentage' | yas == 'percent'){
      jnret$yas <- 'percent'
      jnret$plot$sline <- 0
    } else {
      stop("only 'ratio' or 'prob' ('probability') or 'percent' ('percentage') is supported in the current version")
    }
  }
  yasfunc <- function(jn,value){
    if(jn$yas == 'ratio') return(exp(value))
    else if(jn$yas == 'prob') return(exp(value) / (1 + exp(value)))
    else if(jn$yas == 'percent') return(100*(exp(value) - 1))
    else return(value)
  }

  data <- model$model
  jnret$iv <- iv
  jnret$dv <- dv
  jnret$mod <- mod
  jnret$n <- nrow(model$model)
  df <- nrow(data) - (length(beta.hat) + 1)
  tcrit <- qt(p=alpha/2,df=df,lower.tail = FALSE)

  if(length(jnret$mod) > 1){
    jnret$error <- 'The current version supports only one moderator'
  } else {
    interactionterm <- paste(jnret$iv,jnret$mod,sep=':')

    jna <- tcrit^2*cov[interactionterm, interactionterm] - beta.hat[interactionterm]^2
    jnb <- 2 * (tcrit^2*cov[jnret$iv, interactionterm] - beta.hat[iv]*beta.hat[interactionterm])
    jnc <- tcrit^2*cov[jnret$iv, jnret$iv] - beta.hat[iv]^2

    jnret$plot$summary <- list()
    jnret$plot$summary$upper <- list()
    jnret$plot$summary$lower <- list()
    rs <- c( (-jnb - sqrt(jnb^2 - 4*jna*jnc)) / (2 * jna), (-jnb + sqrt(jnb^2 - 4*jna*jnc)) / (2 * jna))
    jnret$plot$summary$upper$x <- max(rs)
    jnret$plot$summary$upper$y <- beta.hat[jnret$iv] + beta.hat[interactionterm]*jnret$plot$summary$upper$x
    jnret$plot$summary$upper$se <- sqrt(cov[jnret$iv, jnret$iv] + jnret$plot$summary$upper$x^2*cov[interactionterm, interactionterm] + 2*jnret$plot$summary$upper$x*cov[jnret$iv, interactionterm])
    jnret$plot$summary$upper$t <- jnret$plot$summary$upper$y/jnret$plot$summary$upper$se
    jnret$plot$summary$upper$p <- round(2*pt(-abs(jnret$plot$summary$upper$t),df=df),digits=4)
    jnret$plot$summary$upper$llci <- jnret$plot$summary$upper$y - tcrit*jnret$plot$summary$upper$se
    jnret$plot$summary$upper$ulci <- jnret$plot$summary$upper$y + tcrit*jnret$plot$summary$upper$se

    jnret$plot$summary$lower$x <- min(rs)
    jnret$plot$summary$lower$y <- beta.hat[jnret$iv] + beta.hat[interactionterm]*jnret$plot$summary$lower$x
    jnret$plot$summary$lower$se <- sqrt(cov[jnret$iv, jnret$iv] + jnret$plot$summary$lower$x^2*cov[interactionterm, interactionterm] + 2*jnret$plot$summary$lower$x*cov[jnret$iv, interactionterm])
    jnret$plot$summary$lower$t <- jnret$plot$summary$lower$y/jnret$plot$summary$lower$se
    jnret$plot$summary$lower$p <- round(2*pt(-abs(jnret$plot$summary$lower$t),df=df),digits=4)
    jnret$plot$summary$lower$llci <- jnret$plot$summary$lower$y - tcrit*jnret$plot$summary$lower$se
    jnret$plot$summary$lower$ulci <- jnret$plot$summary$lower$y + tcrit*jnret$plot$summary$lower$se

    if(jnret$yas != 'none') {
      jnret$plot$summary$upper$llci <- yasfunc(jnret,jnret$plot$summary$upper$llci)
      jnret$plot$summary$upper$ulci <- yasfunc(jnret,jnret$plot$summary$upper$ulci)
      jnret$plot$summary$upper$y <- yasfunc(jnret, jnret$plot$summary$upper$y)
      jnret$plot$summary$upper$se <- sqrt(jnret$plot$summary$upper$y^2 * jnret$plot$summary$upper$se)

      jnret$plot$summary$lower$llci <- yasfunc(jnret,jnret$plot$summary$lower$llci)
      jnret$plot$summary$lower$ulci <- yasfunc(jnret,jnret$plot$summary$lower$ulci)
      jnret$plot$summary$lower$y <- yasfunc(jnret, jnret$plot$summary$lower$y)
      jnret$plot$summary$lower$se <- sqrt(jnret$plot$summary$lower$y^2 * jnret$plot$summary$lower$se)
    }

    if(jnret$plot$summary$upper$x > max(data[[jnret$mod]]) & jnret$plot$summary$lower$x < min(data[[jnret$mod]])) {
      warning('Conditional effects beyond the range of moderator values')
    }

    jnret$printsummary <- list()
    jnret$printsummary$x <- seq(min(data[[jnret$mod]]):max(data[[jnret$mod]]))
    jnret$printsummary$y <- beta.hat[jnret$iv] + beta.hat[interactionterm]*jnret$printsummary$x
    jnret$printsummary$se <- sqrt(cov[jnret$iv, jnret$iv] + jnret$printsummary$x^2*cov[interactionterm, interactionterm] + 2*jnret$printsummary$x*cov[jnret$iv, interactionterm])
    jnret$printsummary$llci <- jnret$printsummary$y - tcrit*jnret$printsummary$se
    jnret$printsummary$ulci <- jnret$printsummary$y + tcrit*jnret$printsummary$se
    jnret$printsummary$t <- jnret$printsummary$y/jnret$printsummary$se
    jnret$printsummary$p <- round(2*pt(-abs(jnret$printsummary$t),df=df),digits=4)
    if(jnret$yas != 'none') {
      jnret$printsummary$ulci <- yasfunc(jnret,jnret$printsummary$ulci)
      jnret$printsummary$llci <- yasfunc(jnret,jnret$printsummary$llci)
      jnret$printsummary$y <- yasfunc(jnret, jnret$printsummary$y)
      jnret$printsummary$se <- sqrt(jnret$printsummary$y^2 * jnret$printsummary$se)
    }
    jnret$plot$data <- list()
    #conditional effect line
    if(missing(mrange)){
      jnret$plot$data$x <- seq(min(data[[jnret$mod]]), max(data[[jnret$mod]]), length.out = 1000)
    } else {
      jnret$plot$data$x <- seq(min(mrange), max(mrange), length.out = 1000)
    }
    jnret$plot$data$y <- beta.hat[jnret$iv] + beta.hat[interactionterm]*jnret$plot$data$x
    jnret$plot$data$se <- sqrt(cov[jnret$iv, jnret$iv] + jnret$plot$data$x^2*cov[interactionterm, interactionterm] + 2*jnret$plot$data$x*cov[jnret$iv, interactionterm])
    jnret$plot$data$llci <- jnret$plot$data$y - tcrit*jnret$plot$data$se
    jnret$plot$data$ulci <- jnret$plot$data$y + tcrit*jnret$plot$data$se
    if(jnret$yas != 'none') {
      jnret$plot$data$llci <- yasfunc(jnret,jnret$plot$data$llci)
      jnret$plot$data$ulci <- yasfunc(jnret,jnret$plot$data$ulci)
      jnret$plot$data$y <- yasfunc(jnret, jnret$plot$data$y)
      jnret$plot$data$se <- sqrt(jnret$plot$data$y^2 * jnret$plot$data$se)
    }

    jnret$plot$signintervals <- vector()
    trackstartpoly <- 0
    trackstartvalpoly <- vector()
    for(i in 1:length(jnret$plot$data$x)) {
      if(jnret$plot$data$ulci[i] >= jnret$plot$sline & jnret$plot$data$llci[i] >= jnret$plot$sline) {
        if(trackstartpoly == 0) {
          trackstartpoly <- 1
          trackstartvalpoly <- i
        }
      } else if(jnret$plot$data$ulci[i] <= jnret$plot$sline & jnret$plot$data$llci[i] <= jnret$plot$sline) {
        if(trackstartpoly == 0) {
          trackstartpoly <- 1
          trackstartvalpoly <- i
        }
      } else if(trackstartpoly == 1){
        trackstartpoly <- 0
        jnret$plot$signintervals <- rbind(jnret$plot$signintervals,c(trackstartvalpoly,i-1))
        trackstartvalpoly <- vector()
      }
    }
    #wrap up poly calulations
    if(trackstartpoly == 1){
      trackstartpoly <- 0
      jnret$plot$signintervals <- rbind(jnret$plot$signintervals,c(trackstartvalpoly,i))
      trackstartvalspoly <- vector()
    }
  }
  jnret$call <- match.call()
  class(jnret) <- "jn"
  jnret
}
