#' Mixed Effects Logistic Regression for RR Data
#' 
#' Uses the package \code{\link{lme4}} to fit a generalized linear mixed model (GLMM) with an adjusted link funciton.
#' 
#' @param formula two-sided formula including random and fixed effects (see below or \code{\link{glmer}} for details)
#' @param data an optional data frame with variables named in formula
#' @param model type of RR design. Only 1-group RR designs are supported at the moment (i.e., \code{"Warner"}, \code{"FR"}, \code{"UQTknown"}, \code{"Crosswise"}, \code{"Kuk"}, \code{"Mangat"}, \code{"custom"}). See \code{\link{RRuni}} or \code{vignette(RRreg)} for details.
#' @param p randomization probability 
#' @param ... further arguments passed to \code{\link{glmer}}
#' @return an object of class \code{glmerMod}
#' @details 
#' Some examples for formula:
#' \itemize{
#'  \item{random intercept:   }{ \code{response ~ covariate + (1 | group)}}
#'  \item{random slope:   }{ \code{response ~ covariate + (0 + covariate | group)}}
#'  \item{both random slope and intercept:   }{ \code{response ~ covariate +(covariate | group)}}
#'  \item{level-2 predictor (must have constant values within groups!):   }{ \code{response ~ lev2 + (1|group)}}
#' }
#' @references  van den Hout, A., van der Heijden, P. G., & Gilchrist, R. (2007). The Logistic Regression Model with Response Variables Subject to Randomized Response. Computational Statistics & Data Analysis, 51, 6060–6069.
#' @examples 
#' # generate data with a level-1 predictor 
#' d <- data.frame(group=factor(rep(LETTERS[1:20],each=50)), 
#'                 cov=rnorm(20*50))
#' # generate dependent data based on logistic model (random intercept):
#' d$true <- simulate(~  cov + (1|group), newdata=d,
#'                      family=binomial(link="logit"),
#'                      newparams=list(beta=c("(Intercept)"=-.5, cov=1),
#'                                     theta=c("group.(Intercept)"=.8)))[[1]]
#' # scramble responses using RR:
#' model <- "FR"
#' p <- c(true0=.1, true1=.2)
#' d$resp <- RRgen(model="FR", p=p, trueState=d$true)$response
#' # fit model:
#' mod <- RRmixed(resp ~  cov +(1|group), data=d, model="FR", p=p)
#' summary(mod)
#' @import lme4
#' @export
RRmixed <- function(formula, data, model, p, ...){
  
  ######### check if model is allowed
  model <- match.arg(model, modelnames())
  if(is2group(model) | isContinuous(model) )
     stop("Only one-group, dichotomous RR models allowed at the moment (see ?RRmixed)")
  
  ######### get link function
  p <- getPW(model, p)[2,]
  
  ######### fit model using lme4
  mod <- glmer(formula, data, family=binomial(link=RRloglink(p)), ...)

  return(mod)
}




################ LINK FUNCTION #############################

# p = c(0,1)  = probability to respond 1 given a true state of 0 or 1 respectively
# gegeben: p = c(p(1|0), p(1|1))
# gesucht: p(resp==1) = c+d*pi
# p(1) = p(1|0)*(1-pi)     + p(1|1)*pi
# p(1) = p(1|0)- p(1|0)*pi + p(1|1)*pi
# p(1) = p(1|0) + (p(1|1)-p(1|0))*pi
RRloglink <- function(p=c(0,1))
{
  c <- p[1]
  d <- p[2] - p[1]
  linkfun <- function(mu) log((mu-c)/(c+d-mu))
  linkinv <- function(eta) c+d/(1+exp(-eta)) # (c+(c+d)* exp(eta))/(exp(eta)+1)
  mu.eta <- function(eta) (exp(eta)* d)/(exp(eta)+1)^2
  valideta <- function(eta) TRUE
  link <- paste0("logRR(", paste0(p,collapse=","), ")")
  simulate <- 
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}

# Basic checks for link function:
# vv <- logRR(p=c(.1,.9))
# ## check invertibility
# xx <- seq(-30,30, length=100)
# plot(xx, sapply(xx, function(x) vv$linkfun(vv$linkinv(x))), type="l") 
# abline(0,1, col="red")
# library("numDeriv")
# all.equal(grad(vv$linkinv,1),vv$mu.eta(1))  ## check derivative
# plot(xx, sapply(xx, vv$mu.eta), type="l")
# lines(xx, sapply(xx, function(x) grad(vv$linkinv,x)), col="red")