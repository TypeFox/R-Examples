#' Create a new Stochastic Mortality Model
#' 
#' Initialises a StMoMo object which represents a Generalised 
#' Age-Period-Cohort Stochastic Mortality Model.
#' 
#' @details
#' 
#' R implementation of the family of Generalised Age-Period-Cohort 
#' stochastic mortality models.  This family of models encompasses
#' many models proposed in the literature including the well-known
#' Lee-Carter model, CBD model and APC model. 
#' 
#' \code{StMoMo} defines an abstract representation of a Generalised 
#' Age-Period-Cohort (GAPC) Stochastic model that fits within the 
#' general class of generalised non-linear models defined as follows 
#' \deqn{D_{xt} \sim Poisson(E_{xt}\mu_{xt}), D_{xt} \sim 
#' Binomial(E_{xt},q_{xt})} 
#' \deqn{\eta_{xt} = \log \mu_{xt}, \eta_{xt} = \mathrm{logit}\,
#' q_{xt}} \deqn{\eta_{xt} = \alpha_x + \sum_{i=1}^N \beta_x^{(i)}\kappa_t^{(i)}
#' + \beta_x^{(0)}\gamma_{t-x}} \deqn{v: \{\alpha_{x}, \beta_x^{(1)},..., 
#' \beta_x^{(N)}, \kappa_t^{(1)},..., \kappa_t^{(N)}, \beta_x^{(0)}, \gamma_{t-x}\} \mapsto 
#' \{\alpha_{x}, \beta_x^{(1)},..., \beta_x^{(N)}, \kappa_t^{(1)},..., \kappa_t^{(N)}, 
#' \beta_x^{(0)}, \gamma_{t-x}\},} where  
#' \itemize{
#'  \item \eqn{\alpha_x} is a static age function; 
#'  \item \eqn{\beta_x^{(i)}\kappa_t^{(i)}, i = 1,..N}, are age/period terms; 
#'  \item \eqn{\beta_x^{(0)}\gamma_{t-x}} is the age/cohort term; and 
#'  \item \eqn{v} is a  function defining the identifiability constraints of 
#'  the model. 
#' }
#' Most Stochastic mortality models proposed in the literature can be cast to
#' this representation (See Hunt and Blake (2015)).
#' 
#' Parametric age functions should be scalar functions of the form 
#' \code{f <- function(x, ages)} taking a scalar age \code{x} and a vector 
#' of model fitting \code{ages} (see examples below).
#' 
#' Do to limitation of functions \code{\link[gnm]{gnm}} within package 
#' \pkg{gnm}, which is used for fitting \code{"StMoMo"} objects to data 
#' (see \code{\link{fit.StMoMo}}), models combining parametric and 
#' non-parametric age-modulating functions are not supported at the moment.
#' 
#' @seealso \code{\link{fit.StMoMo}}, \code{\link{lc}}, \code{\link{cbd}}, 
#' \code{\link{apc}}, \code{\link{rh}}, \code{\link{m6}}, \code{\link{m7}}, 
#' \code{\link{m8}}
#' 
#' @param link defines the link function and random component associated with 
#'   the mortality model. \code{"log"} would assume that deaths follow a 
#'   Poisson distribution and use a log link while \code{"logit"} would assume 
#'   that deaths follow a Binomial distribution and a logit link.
#'   
#' @param staticAgeFun logical value indicating if a static age function 
#'   \eqn{\alpha_x} is to be included.
#'   
#' @param periodAgeFun  a list of length \eqn{N} with the definitions of the 
#'   period age modulating parameters \eqn{\beta_x^{(i)}}. Each entry can take 
#'   values: \code{"NP"} for non-parametric age terms, \code{"1"} for 
#'   \eqn{\beta_x^{(i)}=1} or a predefined parametric function of 
#'   age (see details). Set this to \code{NULL} if there are no period terms 
#'   in the model.
#'   
#' @param cohortAgeFun defines the cohort age modulating parameter 
#'   \eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for non-parametric 
#'   age terms, \code{"1"} for \eqn{\beta_x^{(0)}=1}, a predefined parametric 
#'   function of age (see details) or \code{NULL} if there is no cohort effect. 
#'   
#' @param  constFun function defining the identifiability constraints of the 
#'   model. It must be a function of the form 
#'   \code{constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)} taking a set
#'   of fitted model parameters and returning a list 
#'   \code{list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
#'   of the model parameters with the identifiability constraints applied. If 
#'   omitted no identifiability constraints are applied to the model.
#'   
#' @return A list with class \code{"StMoMo"} with components:
#'   
#'   \item{link}{a character string defining the link function of the model.}
#'   
#'   \item{staticAgeFun}{a logical value indicating if the model has a static 
#'   age function.}
#'   
#'   \item{periodAgeFun}{a list defining the period age modulating parameters.}
#'   
#'   \item{cohortAgeFun}{an object defining the cohort age modulating 
#'   parameters.}
#'   
#'   \item{constFun}{a function defining the identifiability constraints.}
#'   
#'   \item{N}{an integer specifying The number of age-period terms in the 
#'   model.}
#'   
#'   \item{textFormula}{a character string of the model formula.}
#'   
#'   \item{gnmFormula}{a formula that can be used for fitting the model with 
#'   package \pkg{gnm}.}
#'   
#' @references 
#' 
#' Plat, R. (2009). On stochastic mortality modeling. Insurance: Mathematics 
#'   and Economics, 45(3), 393-404. 
#' 
#' Hunt, A., & Blake, D. (2015). On the Structure and Classification of 
#' Mortality Models Mortality Models. Pension Institute Working Paper.
#' \url{http://www.pensions-institute.org/workingpapers/wp1506.pdf}.
#' 
#'   
#' @examples
#' #Lee-Carter model
#' constLC <- function(ax, bx, kt, b0x, gc, wxt, ages) {
#'  c1 <- mean(kt[1, ], na.rm = TRUE)
#'  c2 <- sum(bx[, 1], na.rm = TRUE)
#'  list(ax = ax + c1 * bx, bx = bx / c2, kt = c2 * (kt - c1))
#' }
#' LC <- StMoMo(link = "log", staticAgeFun = TRUE, periodAgeFun = "NP",
#'              constFun = constLC)
#' 
#' plot(fit(LC, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'          ages = EWMaleData$ages, years = EWMaleData$years,
#'          ages.fit = 55:89))
#'
#' 
#' #CBD model   
#' f2 <- function(x, ages) x - mean(ages)
#' CBD <- StMoMo(link = "logit", staticAgeFun = FALSE,
#'               periodAgeFun = c("1", f2))
#' plot(fit(CBD, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'          ages = EWMaleData$ages, years = EWMaleData$years, 
#'          ages.fit = 55:89))
#' 
#' #Reduced Plat model (Plat, 2009)
#' f2 <- function(x, ages) mean(ages) - x
#' constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages) {
#'  nYears <- dim(wxt)[2]
#'  x <- ages
#'  t <- 1:nYears
#'  c <- (1 - tail(ages, 1)):(nYears - ages[1])
#'  xbar <- mean(x)
#'  #nsum g(c)=0, nsum cg(c)=0, nsum c^2g(c)=0
#'  phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
#'  phi <- coef(phiReg)
#'  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
#'  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
#'  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
#'  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
#'  #nsum kt[i, ] = 0
#'  ci <- rowMeans(kt, na.rm = TRUE)
#'  ax <- ax + ci[1] + ci[2] * (xbar - x)
#'  kt[1, ] <- kt[1, ] - ci[1]
#'  kt[2, ] <- kt[2, ] - ci[2]
#'  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
#' }
#' PLAT <- StMoMo(link = "log", staticAgeFun = TRUE,
#'                periodAgeFun = c("1", f2), cohortAgeFun = "1",
#'                constFun = constPlat)
#' 
#' plot(fit(PLAT, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
#'          ages = EWMaleData$ages, years = EWMaleData$years, 
#'          ages.fit = 55:89))
#' 
#' #Models not supported
#' \dontrun{
#' MnotSup1 <- StMoMo(periodAgeFun = c(f2, "NP"))
#' MnotSup1 <- StMoMo(periodAgeFun = f2, cohortAgeFun = "NP")
#' }
#' @export
StMoMo  <- function(link = c("log", "logit"), staticAgeFun = TRUE, 
                    periodAgeFun = 'NP', cohortAgeFun = NULL, 
                    constFun = function(ax, bx, kt, b0x, gc, wxt, ages) 
                      list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)) {
  
  
  #---------------------------------------------------------------------
  # Construct the model formula
  #---------------------------------------------------------------------
  
  #link
  link <- match.arg(link)
  if (link == "logit") {
    gnmFormula <- "D/E ~ -1 + offset(o)" 
    textFormula <- "logit q[x,t]"
  } else if (link == "log") {
    gnmFormula <- "D ~ -1 + offset(log(E)) + offset(o)"
    textFormula <- "log m[x,t]"
  }
  
  #static life table  
  termSep <- " = "
  if (staticAgeFun == TRUE) {
    gnmFormula <- paste(gnmFormula, " + factor(x)", sep = "")    
    textFormula <- paste(textFormula, termSep, "a[x]", sep = "")    
    termSep <- " + "
  }
  
  #age-period terms
  periodAgeFun <- c(periodAgeFun)
  N <- length(periodAgeFun)
  if (N>0) {            
    for (i in 1:N) {
      if (is.function(periodAgeFun[[i]])) {
        gnmFormula <- paste(gnmFormula, " + B", i, ":factor(t)", sep = "")
        textFormula <- paste(textFormula, termSep, "f", i, "[x] k", i, "[t]", 
                             sep = "")
      } else if (periodAgeFun[[i]] == "1") {
        gnmFormula <- paste(gnmFormula, " + factor(t)", sep = "")      
        textFormula <- paste(textFormula, termSep, "k", i, "[t]", sep = "")
      } else if ( periodAgeFun[[i]] == "NP" ) {
        ind <- i:N
        inst <- 1 + sum(periodAgeFun[-ind] == periodAgeFun[[i]]) 
        gnmFormula <- paste(gnmFormula, 
                            " + Mult(factor(x), factor(t), inst = ", inst,
                            ")", sep = "")   
        textFormula <- paste(textFormula, termSep, "b", i, "[x] k", i, "[t]", 
                             sep = "")
      } else {
        stop("Not appropriate period age modulating function definition")
      }      
      termSep <- " + "
    }    
    totalP <- sum(sapply(periodAgeFun, is.function))
    totalNP <- sum(periodAgeFun == "NP")
  } else {
    totalP <- 0
    totalNP <- 0    
  } 
  
  #age-cohort term
  if (!is.null(cohortAgeFun)) {
    if (is.function(cohortAgeFun)) {
      gnmFormula <- paste(gnmFormula, " + B0:factor(c)", sep = "")
      textFormula <- paste(textFormula, termSep, "f0[x] g[t-x]", sep = "")
    } else if (cohortAgeFun == "1") {
      gnmFormula <- paste(gnmFormula, " + factor(c)", sep = "")      
      textFormula <- paste(textFormula, termSep, "g[t-x]", sep = "")
    } else if (cohortAgeFun == "NP") {
      gnmFormula <- paste(gnmFormula, " + Mult(factor(x), factor(c))", 
                          sep = "")   
      textFormula <- paste(textFormula, termSep, "b0[x] g[t-x]", sep = "")
    } else {
      stop("Not appropriate cohort age modulating function definition")
    }
    totalP <- totalP + is.function(cohortAgeFun)
    totalNP <-totalNP + sum(list(cohortAgeFun) == "NP")
  }
  
  #Check constraint on number of simulataneous parametric and non-parametric 
  #age functions
  if (totalP >= 1 && totalNP >= 1) {
    stop("Models combining parametric and non-parametric age-modulating functions are not supported at the moment.")
  }
  
  out <- list(link = link, staticAgeFun = staticAgeFun, 
              periodAgeFun = periodAgeFun, cohortAgeFun = cohortAgeFun, 
              constFun = constFun, N = N, textFormula = textFormula,
              gnmFormula = gnmFormula)
  class(out) <- "StMoMo"
  out
}


#' Print an object of class \code{"StMoMo"}
#' 
#' \code{print} method for class \code{"StMoMo"}.
#' 
#' \code{print.StMoMo} prints a description of the Stochastic Mortality
#' Model 
#' 
#' @usage 
#' \method{print}{StMoMo}(x, ...)
#' 
#' @param x an object of class \code{"StMoMo"}.
#' @param ... arguments to be passed to or from other methods.
#' @export 
#' @method print StMoMo
print.StMoMo <- function(x, ...) {
  if (x$link == "logit") {
    cat("Binomial model with predictor: ")
  } else {
    cat("Poisson model with predictor: ")
  }  
  cat("")
  cat(x$textFormula)
}

