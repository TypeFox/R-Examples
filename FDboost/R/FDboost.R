#################################################################################
#' Model-based Gradient Boosting for Functional Response
#' 
#' Gradient boosting for optimizing arbitrary loss functions, where component-wise models 
#' are utilized as base-learners in the case of functional responses. 
#' Scalar responses are treated as the special case where each functional response has 
#' only one observation. 
#' This function is a wrapper for \code{mboost}'s \code{\link[mboost]{mboost}} and its 
#' siblings to fit models of the general form 
#' \deqn{\xi(Y_i(t) | X_i = x_i) = \sum_{j} h_j(x_i, t), i = 1, ..., N,} 
#' with a functional (but not necessarily continuous) response \eqn{Y(t)}, 
#' transformation function \eqn{\xi}, e.g., the expectation, the median or some quantile, 
#' and partial effects \eqn{h_j(x_i, t)} depending on covariates \eqn{x_i}  
#' and the current index of the response \eqn{t}. The index of the response can 
#' be for example time.  
#' Possible effects are, e.g., a smooth intercept \eqn{\beta_0(t)}, 
#' a linear functional effect \eqn{\int x_i(s)\beta(s,t)ds}, 
#' potentially with integration limits depending on \eqn{t}, 
#' smooth and linear effects of scalar covariates \eqn{f(z_i)} or \eqn{z_i \beta(t)}. 
#' 
#' @param formula a symbolic description of the model to be fit. 
#' Per default no intercept is added, only a smooth offset, see argument \code{offset}. 
#' To add a smooth intercept, use 1, e.g. \code{y ~ 1} for a pure intercept model. 
#' @param timeformula one-sided formula for the specification of the effect over the index of the response. 
#' For functional response \eqn{Y_i(t)} typically use \code{~ bbs(t)} to obtain smooth 
#' effects over \eqn{t}. 
#' In the limiting case of \eqn{Y_i} being a scalar response, 
#' use \code{~ bols(1)}, which sets up a base-learner for the scalar 1. 
#' Or use \code{timeformula = NULL}, then the scalar response is treated as scalar. 
#' @param id defaults to NULL which means that all response trajectories are observed
#' on a common grid allowing to represent the response as a matrix. 
#' If the response is given in long format for observation-specific grids, \code{id} 
#' contains the information which observations belong to the same trajectory and must 
#' be supplied as a formula, \code{~ nameid}, where the variable \code{nameid} should 
#' contain integers 1, 2, 3, ..., N. 
#' @param numInt integration scheme for the integration of the loss function.
#' One of \code{c("equal", "Riemann")} meaning equal weights of 1 or 
#' trapezoidal Riemann weights.
#' Alternatively a vector of length \code{nrow(response)} containing  
#' positive weights can be specified.
#' @param data a data frame or list containing the variables in the model.
#' @param weights (1) a numeric vector of weights for observational units, 
#' i.e. \code{length(weights)} has to be \code{nrow(response)},
#' (2) alternatively weights can be specified for single observations then
#' \code{length(weights)} has to be \code{nrow(response)}*\code{ncol(response)}
#' per default weights is constantly 1. 
#' @param offset_control parameters for the estimation of the offset, 
#' defaults to \code{o_control(k_min = 20, silent = TRUE)}, see \code{\link{o_control}}.  
#' @param offset a numeric vector to be used as offset over the index of the response (optional).
#' If no offset is specified, per default \code{offset = NULL} which means that a 
#' smooth time-specific offset is computed and used before the model fit to center the data. 
#' If you do not want to use a time-specific offset, set \code{offset = "scalar"} to get an overall scalar offset, 
#' like in \code{mboost}.
#' @param check0 logical, for response in matrix form, i.e. response that is observed on a common grid, 
#' check the fitted effects for the sum-to-zero constraint 
#' \eqn{h_j(x_i)(t) = 0} for all \eqn{t} and give a warning if it is not fulfilled. Defaults to \code{FALSE}. 
#' @param ... additional arguments passed to \code{\link[mboost]{mboost}}, 
#' including, \code{family} and \code{control}.
#' 
#' @details In matrix representation of functional response and covariates each row 
#' represents one functional observation, e.g. \code{Y[i,t_g]} corresponds to \eqn{Y_i(t_g)}, 
#' giving a <number of curves> by <number of evaluations> matrix. 
#' For the model fit, the matrix of the functional
#' response evaluations \eqn{Y_i(t_g)} are stacked internally into one long vector. 
#' 
#' If it is possible to represent the model as a generalized linear array model 
#' (Currie et al., 2006), the array structure is used for an efficient implementation, 
#' see \code{\link[mboost]{mboost}}. This is only possible if the design 
#' matrix can be written as the Kronecker product of two marginal design 
#' matrices yielding a functional linear array model (FLAM), 
#' see Brockhaus et al. (2015) for details. 
#' The Kronecker product of two marginal bases is implemented in R-package mboost 
#' in the function \code{\%O\%}, see \code{\link[mboost]{\%O\%}}. 
#' 
#' When \code{\%O\%} is called with a specification of \code{df} in both base-learners, 
#' e.g. \code{bbs(x1, df = df1) \%O\% bbs(t, df = df2)}, the global \code{df} for the 
#' Kroneckered base-learner is computed as \code{df = df1 * df2}. 
#' And thus the penalty has only one smoothness parameter lambda resulting in an isotropic penalty. 
#' A Kronecker product with anisotropic penalty is \code{\%A\%}, allowing for different 
#' amount of smoothness in the two directions, see \code{\link{\%A\%}}. 
#' If the formula contains base-learners connected by \code{\%O\%} or \code{\%A\%}, 
#' those effects are not expanded with \code{timeformula}, allowing for model specifications 
#' with different effects in time-direction.   
#' 
#' If the response is observed on curve-specific grids it must be supplied  
#' as a vector in long format and the argument \code{id} has  
#' to be specified (as formula!) to define which observations belong to which curve.  
#' In this case the base-learners are built as row tensor-products of marginal base-learners, 
#' see Scheipl et al. (2015) and Brockhaus et al. (2016), for details on how to set up the effects. 
#' The row tensor product of two marginal bases is implemented in R-package mboost 
#' in the function \code{\%X\%}, see \code{\link[mboost]{\%X\%}}. 
#' 
#' A scalar response can be seen as special case of a functional response with only
#' one time-point, and thus it can be represented as FLAM with basis 1 in 
#' time-direction, use \code{timeformula = ~bols(1)}. In this case, a penalty in the 
#' time-direction is used, see Brockhaus et al. (2015) for details.  
#' Alternatively, the scalar response is fitted as scalar response, like in the function
#' \code{\link[mboost]{mboost}} in package mboost. 
#' The advantage of using \code{FDboost} in that case 
#' is that methods for the functional base-learners are available, e.g. \code{plot}. 
#' 
#' The desired regression type is specified by the \code{family}-argument, 
#' see the help-page of \code{\link[mboost]{mboost}}. For example a mean regression model is obtained by  
#' \code{family = Gaussian()} which is the default or median regression 
#' by \code{family = QuantReg()}; 
#' see \code{\link[mboost]{Family}} for a list of implemented families. 
#' 
#' With \code{FDboost} the following covariate effects can be estimated by specifying 
#' the following effects in the \code{formula}
#' (similar to function \code{\link[refund]{pffr}} in R-package \code{\link[refund]{refund}}). 
#' The \code{timeformula} is used to expand the effects in \code{t}-direction. 
#' \itemize{
#' \item Linear functional effect of scalar (numeric or factor) covariate \eqn{z} that varies 
#'   smoothly over \eqn{t}, i.e. \eqn{z_i \beta(t)}, specified as
#'   \code{~ bolsc(z)}, see \code{\link{bolsc}}, 
#'   or for a group effect with mean zero use \code{~ brandomc(z)}.  
#' \item Nonlinear effects  of a scalar covariate that vary smoothly over \eqn{t}, 
#'   i.e. \eqn{f(z_i, t)}, specified as \code{bbsc(z)}, 
#'   see \code{\link{bbsc}}. 
#' \item (Nonlinear) effects of scalar covariates that are constant 
#'   over \eqn{t}, e.g. \eqn{f(z_i)}, specified as \code{~ c(bbs(z))}, 
#'   or \eqn{\beta z_i}, specified as \code{~ c(bols(z))}.
#' \item Interaction terms between two scalar covariates, e.g. \eqn{z_i1 zi2 \beta(t)}, 
#'   are specified as \code{bolsc(z1) \%Xc\% bolsc(z2)} and  
#'   an interaction \eqn{z_i1 f(zi2, t)} as \code{bolsc(z1) \%Xc\% bbsc(z2)}, as 
#'   \code{\%Xc\%} applies the sum-to-zero constraint to the desgin matrix of the tensor product 
#'   built by \code{\%Xc\%}, see  \code{\link{\%Xc\%}}, EXPERIMENTAL!
#' \item Function-on-function regression terms of functional covariates \code{x}, 
#'   e.g. \eqn{\int x_i(s)\beta(s,t)ds}, specified as \code{~ bsignal(x, s = s)}, 
#'   using P-splines, see \code{\link{bsignal}}. 
#'   Terms given by \code{\link{bfpc}} provide FPC-based effects of functional 
#'   covariates, see \code{\link{bfpc}}. 
#' \item Function-on-function regression terms of functional covariates \code{x} 
#'   with integration limits \eqn{[l(t), u(t)]} depending on \eqn{t},  
#'   e.g. \eqn{\int_[l(t), u(t)] x_i(s)\beta(s,t)ds}, specified as 
#'   \code{~ bhist(x, s = s, time = t, limits)}. The \code{limits} argument defaults to
#'   \code{"s<=t"} which yields a hisotircal effect with limits \eqn{[min(t),t]}, 
#'    see \code{\link{bhist}}.
#' \item Concurrent effects of functional covariates \code{x}
#'   measured on the same grid as the response, i.e., \eqn{x_i(s)\beta(t)}, 
#'   are specified as \code{~ bconcurrent(x, s = s, time = t)}, 
#'   see \code{\link{bconcurrent}}. 
#' \item Interaction effects can be estimated as tensor product smooth, e.g., 
#'   \eqn{ z \int x_i(s)\beta(s,t)ds} as \code{~ bsignal(x, s = s) \%X\% bolsc(z)}
#' \item For interaction effects with historical functional effects, e.g., 
#'   \eqn{ z_i \int_[l(t),u(t)] x_i(s)\beta(s,t)ds} the base-learner 
#'   \code{bhistx} should be used instead of \code{bhist}, 
#'   e.g. \code{~ bhistx(x, limits) \%X\% bolsc(z)}, see \code{\link{bhistx}}.
#' \item Generally, the \code{c()}-notation can be used to get effects that are 
#'   constant over the index of the functional response. 
#' \item If the \code{formula} in \code{FDboost} contains base-learners connected by 
#' \code{\%O\%} or \code{\%A\%}, those effects are not expanded with \code{timeformula}, 
#' allowing for model specifications with different effects in time-direction.  
#' } 
#'  
#' In order to obtain a fair selection of base-learners, the same  degrees of freedom (df) 
#' should be specified for all baselearners. It is recommended to use 
#' a rather small number of df for all base-learners. 
#' It is not possible to specify df larger than the rank of the design matrix.
#' For base-learners with rank-deficient penalty it is not possible to specify df smaller than the 
#' rank of the null space of the penalty (e.g., in \code{bbs} unpenalized part of P-splines). 
#' The df of the base-learners in an FDboost-object can be checked using \code{extract(object, "df")}, 
#' see \code{\link[mboost]{extract}}.  
#' 
#' The most important tuning parameter of component-wise gradient boosting 
#' is the number of boosting iterations. It is recommended to use the number of 
#' boosting iterations as only tuning parameter, 
#' fixing the step-length at a small value (e.g. nu = 0.1). 
#' Note that the default number of boosting iterations is 100 which is arbitrary and in most 
#' cases not adequate (the optimal number of boosting iterations can considerably exceed 100). 
#' The optimal stopping iteration can be determined by resampling methods like
#' cross-validation or bootstrapping, see the function \code{\link{cvrisk.FDboost}} which searches 
#' the optimal stopping iteration on a grid, which in many cases has to be extended.  
#' 
#' @return An object of class \code{FDboost} that inherits from \code{mboost}.
#' Special \code{\link{predict.FDboost}}, \code{\link{coef.FDboost}} and 
#' \code{\link{plot.FDboost}} methods are available. 
#' The methods of \code{\link[mboost]{mboost}} are available as well, 
#' e.g., \code{\link[mboost]{extract}}. 
#' 
#' The \code{FDboost}-object is a named list containing: 
#' \item{...}{all elements of an \code{\link[mboost]{mboost}-object}}
#' \item{yname}{the name of the response}
#' \item{ydim}{dimension of the response matrix, if the response is represented as such}
#' \item{yind}{the observation (time-)points of the response, i.e. the evaluation points, 
#'       with its name as attribute}
#' \item{data}{the data that was used for the model fit}
#' \item{id}{the id variable of the response}
#' \item{predictOffset}{the function to predict the smooth offset}
#' \item{offsetVec}{the offset for one trajectory if the response is represented as matrix and 
#' otherwise the offset for all trajectories}
#' \item{offsetFDboost}{offset as specified in call to FDboost} 
#' \item{offsetMboost}{offset as given to mboost}
#' \item{call}{the call to \code{FDboost}}
#' \item{callEval}{the evaluated function call to \code{FDboost} without data}
#' \item{timeformula}{the time-formula}
#' \item{formulaFDboost}{the formula with which \code{FDboost} was called}
#' \item{formulaMboost}{the formula with which \code{mboost} was called within \code{FDboost}}
#' 
#' @author Sarah Brockhaus, Torsten Hothorn
#' 
#' @seealso Note that \link{FDboost} calls \code{\link[mboost]{mboost}} directly.  
#' See, e.g., \code{\link[FDboost]{bsignal}} and \code{\link[FDboost]{bbsc}} 
#' for possible base-learners. 
#' 
#' @keywords models, nonlinear 
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300. 
#' 
#' Brockhaus, S., Melcher, M., Leisch, F. and Greven, S. (2016): 
#' Boosting flexible functional regression models with a high number of functional historical effects,  
#' under revision.  
#' 
#' Currie, I.D., Durban, M. and Eilers P.H.C. (2006):  
#' Generalized linear array models with applications to multidimensional smoothing. 
#' Journal of the Royal Statistical Society, Series B-Statistical Methodology, 68(2), 259-280.
#' 
#' Scheipl, F., Staicu, A.-M. and Greven, S. (2015):  
#' Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 24(2), 477-501.
#' 
#' @examples 
#' ######## Example for function-on-scalar-regression 
#' data("viscosity", package = "FDboost") 
#' ## set time-interval that should be modeled
#' interval <- "101"
#' 
#' ## model time until "interval" and take log() of viscosity
#' end <- which(viscosity$timeAll == as.numeric(interval))
#' viscosity$vis <- log(viscosity$visAll[,1:end])
#' viscosity$time <- viscosity$timeAll[1:end]
#' # with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))
#' 
#' ## fit median regression model with 100 boosting iterations,
#' ## step-length 0.4 and smooth time-specific offset
#' ## the factors are coded such that the effects are zero for each timepoint t
#' ## no integration weights are used!
#' mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df = 2) + bolsc(T_A, df = 2),
#'                timeformula = ~ bbs(time, df = 4),
#'                numInt = "equal", family = QuantReg(),
#'                offset = NULL, offset_control = o_control(k_min = 9),
#'                data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#' 
#' \dontrun{
#'   ## find optimal mstop over 5-fold bootstrap, small number of folds for example
#'   ## do the resampling on the level of curves
#'   set.seed(123)
#'   val1 <- validateFDboost(mod1, folds = cv(rep(1, length(unique(mod1$id))), B = 5), 
#'                         grid = 1:500)
#'   ## plot(val1)
#'   mstop(val1)
#'   mod1[mstop(val1)]
#' 
#'   ## find the optimal mstop over 5-fold bootstrap
#'   ## using the function cvrisk; be careful to do the resampling on the level of curves
#'   folds1 <- cvLong(id = mod1$id, weights = model.weights(mod1), B = 5)
#'   cvm1 <- cvrisk(mod1, folds = folds1, grid = 1:500)
#'   ## plot(cvm1)
#'   mstop(cvm1)
#' }
#' 
#' ## look at the model
#' summary(mod1)
#' coef(mod1)
#' ## plot(mod1)
#' ## plotPredicted(mod1, lwdPred = 2)
#' 
#' ######## Example for scalar-on-function-regression 
#' data("fuelSubset", package = "FDboost")
#' 
#' ## center the functional covariates per observed wavelength
#' fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
#' fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)
#' 
#' ## to make mboost:::df2lambda() happy (all design matrix entries < 10)
#' ## reduce range of argvals to [0,1] to get smaller integration weights
#' fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) / 
#'                                           (max(uvvis.lambda) - min(uvvis.lambda) ))
#' fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) / 
#'                                           (max(nir.lambda) - min(nir.lambda) )) 
#' 
#' ## possibility 1: as FLAM model with the scalar response 
#' ## adds a penalty over the index of the response
#' ## thus, mod2f and mod2 have different panlties 
#' mod2f <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
#'                + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE), 
#'                timeformula = ~bols(1), data = fuelSubset, control = boost_control(mstop = 200))
#' 
#' ## possibility 2: with scalar response 
#' mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
#'                + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE), 
#'                timeformula = NULL, data = fuelSubset, control = boost_control(mstop = 200)) 
#'     
#' ## bootstrap to find optimal mstop takes some time
#' ## set.seed(123)      
#' ## folds2 <- cv(weights = model.weights(mod2), B = 10)     
#' ## cvm2 <- cvrisk(mod2, folds = folds2, grid = 1:1000)
#' ## mstop(cvm2)
#' mod2[327]
#'                
#' summary(mod2) 
#' ## plot(mod2)
#' 
#' ## Example for function-on-function-regression 
#' if(require(fda)){
#' 
#'   data("CanadianWeather", package = "fda")
#'   CanadianWeather$l10precip <- t(log(CanadianWeather$monthlyPrecip))
#'   CanadianWeather$temp <- t(CanadianWeather$monthlyTemp)
#'   CanadianWeather$region <- factor(CanadianWeather$region)
#'   CanadianWeather$month.s <- CanadianWeather$month.t <- 1:12
#'   
#'   ## center the temperature curves per time-point
#'   CanadianWeather$temp <- scale(CanadianWeather$temp, scale = FALSE)
#'   
#'   ## fit model with cyclic splines over the year
#'   mod3 <- FDboost(l10precip ~ bols(region, df = 2.5, contrasts.arg = "contr.dummy") 
#'                    + bsignal(temp, month.s, knots = 11, cyclic = TRUE, 
#'                            df = 2.5, boundary.knots = c(0.5,12.5), check.ident = FALSE), 
#'                    timeformula = ~ bbs(month.t, knots = 11, cyclic = TRUE, 
#'                                   df = 3, boundary.knots = c(0.5, 12.5)), 
#'                    offset = "scalar", offset_control = o_control(k_min = 5), 
#'                   data = CanadianWeather)  
#'  ## find the optimal mstop over 5-fold bootstrap 
#'  ## using the function cvrisk; be careful to do the resampling on the level of curves
#'  ## set.seed(123)
#'  ## folds3 <- cvLong(id = mod3$id, weights = model.weights(mod3), B = 5)
#'  ## cvm3 <- cvrisk(mod3, folds = folds3, grid = 1:500)
#'  ## mstop(cvm3)
#'  mod3[64]
#'    
#'  summary(mod3)
#'  ## plot(mod3, pers = TRUE)
#' }
#' 
#' ######## Example for functional response observed on irregular grid
#' ######## Delete part of observations in viscosity data-set
#' data("viscosity", package = "FDboost")
#' ## set time-interval that should be modeled
#' interval <- "101"
#' 
#' ## model time until "interval" and take log() of viscosity
#' end <- which(viscosity$timeAll == as.numeric(interval))
#' viscosity$vis <- log(viscosity$visAll[,1:end])
#' viscosity$time <- viscosity$timeAll[1:end]
#' # with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))
#' 
#' ## only keep a quarter of the observation points
#' set.seed(123)
#' selectObs <- sort(sample(x = 1:(64*46), size = 64*46/4, replace = FALSE))
#' dataIrregular <- with(viscosity, list(vis = c(vis)[selectObs], 
#'                                       T_A = T_A, T_C = T_C,  
#'                                       time = rep(time, each = 64)[selectObs], 
#'                                       id = rep(1:64, 46)[selectObs]))
#' 
#' ## fit median regression model with 100 boosting iterations,
#' ## step-length 0.4 and smooth time-specific offset
#' ## the factors are in effect coding -1, 1 for the levels
#' ## no integration weights are used!
#' mod4 <- FDboost(vis ~ 1 + bols(T_C, contrasts.arg = "contr.sum", intercept = FALSE)
#'                 + bols(T_A, contrasts.arg = "contr.sum", intercept=FALSE),
#'                 timeformula = ~ bbs(time, lambda = 100), id = ~id, 
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = dataIrregular, control = boost_control(mstop = 100, nu = 0.4))
#' summary(mod4)
#' ## plot(mod4)
#' ## plotPredicted(mod4, lwdPred = 2)
#' 
#' ## Find optimal mstop, small grid/low B for a fast example
#' folds4 <- cvLong(id = mod4$id, weights = model.weights(mod4), B = 3)
#' cvm4 <- cvrisk(mod4, folds = folds4, grid = 1:50)
#' mstop(cvm4)
#' ## val4 <- validateFDboost(mod4, folds = cv(rep(1, length(unique(mod4$id))), B = 3), 
#' ##                         grid = 1:50)
#' 
#' ## Be careful if you want to predict newdata with irregular response,  
#' ## as the argument index is not considered in the prediction of newdata. 
#' ## Thus, all covariates have to be repeated according to the number of observations 
#' ## in each response trajectroy. 
#' ## Predict four response curves with full time-observations 
#' ## for the four combinations of T_A and T_C. 
#' newd <- list(T_A = factor(c(1,1,2,2), levels = 1:2, 
#'                         labels = c("low", "high"))[rep(1:4, length(viscosity$time))], 
#'              T_C = factor(c(1,2,1,2), levels = 1:2, 
#'                         labels = c("low", "high"))[rep(1:4, length(viscosity$time))], 
#'              time = rep(viscosity$time, 4))
#'              
#' ## pred <- predict(mod4, newdata = newd)
#' ## funplot(x = rep(viscosity$time, 4), y = pred, id = rep(1:4, length(viscosity$time)))
#'                   
#'                 
#' @export
#' @import methods Matrix mboost  
#' @importFrom grDevices heat.colors rgb
#' @importFrom graphics abline barplot contour legend lines matplot par persp plot points
#' @importFrom utils getS3method 
#' @importFrom stats approx as.formula coef complete.cases fitted formula lm median model.matrix model.weights na.omit predict quantile sd terms.formula variable.names 
#' @importFrom gamboostLSS GaussianLSS GaussianMu GaussianSigma make.grid cvrisk.mboostLSS
#' @importFrom splines bs splineDesign
#' @importFrom mgcv gam s
#' @importFrom zoo na.locf
#' @importFrom MASS Null
#' @importFrom parallel mclapply
#' @importFrom refund fpca.sc
FDboost <- function(formula,          ### response ~ xvars
                    timeformula,      ### time
                    id = NULL,          ### id variable if response is in long format
                    numInt = "equal",   ### option for approximation of integral over loss
                    data,             ### list of response, time, xvars
                    weights = NULL,   ### optional
                    offset = NULL,    ### optional
                    offset_control = o_control(), ### optional specification of offset model
                    check0 = FALSE,    ### check sum-to-zero-constraint of the fitted effects?
                    #offset_control = list(k_min = 20, silent = TRUE),
                    ...)              ### goes directly to mboost
{
  dots <- list(...)

  ## insert the id variable into the formula, to treat it like the other variables
  if(!is.null(id)){
    stopifnot(class(id) == "formula")
    tf <- terms.formula(formula, specials = c("c"))
    trmstrings <- attr(tf, "term.labels")
    if(length(trmstrings) > 0){
      ## insert id at end of each base-learner
      trmstrings2 <- paste(substr(trmstrings, 1 , nchar(trmstrings)-1), ", index=", id[2],")", sep = "")
      ## insert into the other base-learners of the tensor-product as well
      for(i in 1:length(trmstrings)){
        if(grepl( "%X", trmstrings2[i])){
          temp <- unlist(strsplit(trmstrings2[i], "%X"))
          temp1 <- temp[-length(temp)]
          ## http://stackoverflow.com/questions/2261079
          ## delete all trailing whitespace
          trim.trailing <- function (x) sub("\\s+$", "", x) 
          temp1 <- trim.trailing(temp1)
          temp1 <- paste(substr(temp1, 1 , nchar(temp1)-1), ", index=", id[2],")", sep = "")
          trmstrings2[i] <- paste0(paste0(temp1, collapse = " %X"), " %X", temp[length(temp)]) 
        } 
        ## do not add index to base-learners bhistx()
        if( grepl("bhistx", trmstrings[i]) ) trmstrings2[i] <- trmstrings[i] 
      }
      ### <FIXME> do not add an index if an index is already part of the formula
      #trmstrings2[grepl("index", trmstrings)] <- trmstrings[grepl("index", trmstrings)]
      trmstrings <- trmstrings2
    } 
    xpart <- paste(as.vector(trmstrings), collapse = " + ")
    if(xpart != ""){
      if(any(substr(tf[[3]], 1, 1) == "1")) xpart <- paste0("1 + ", xpart)
    }else{
      xpart <- 1
    }
    formula <- as.formula(paste(tf[[2]], " ~ ", xpart))
    #print(formula)
    nameid <- paste(id[2])
    id <- data[[nameid]]
  }else{
    nameid <- NULL
  }
  
  ### save formula of FDboost
  formulaFDboost <- formula
  
  stopifnot(class(formula) == "formula")
  if(!is.null(timeformula)) stopifnot(class(timeformula) == "formula")
  
  ### extract response; a numeric matrix or a vector
  yname <- all.vars(formula)[1]
  response <- data[[yname]]
  data[[yname]] <- NULL
  
  ### for scalar response ~bols(1) or NULL
  scalarResponse <- FALSE
  scalarNoFLAM <- FALSE
  if(is.null(timeformula) || timeformula == ~bols(1)){
    
    scalarResponse <- TRUE
    if(is.null(timeformula)) scalarNoFLAM <- TRUE
    
    if(grepl("df", formula[3]) | !grepl("lambda", formula[3]) ){
      timeformula <- ~bols(ONEtime, intercept = FALSE, df = 1)
    }else{
      timeformula <- ~bols(ONEtime, intercept = FALSE)
    }
    
    ## <FIXME> would be nice, as then no penalization in dummy-direction happens
    ## timeformula <- ~bols(ONEtime, lambda = 0)
    
    data$ONEtime <- 1
    response <- matrix(response, ncol = 1)
  }
  
  ### <TODO> find a reasonable way to warn the user in the case of non-identifiable models?
  #   else{   
  #     ### <FIXME> option to switch off those checks?    
  #     ## in the case of functional response: check formula to use identifiable base-learners
  #     tf <- terms.formula(formula, specials = c("c"))
  #     trmstrings <- attr(tf, "term.labels")    
  #     if(length(trmstrings) > 0){
  #       if(any(grepl("bbs(", trmstrings, fixed = TRUE))){
  #         warning("Use bbsc() instead of bbs() with functional response, to get an identifiable model.")
  #       }      
  #       if(any(grepl("bols(", trmstrings, fixed = TRUE))){
  #         temp <- trmstrings[grepl("bols(", trmstrings, fixed = TRUE)]
  #         #temp[[1]]
  #         warning("Use bolsc() instead of bols() with functional response, to get an identifiable model.")
  #       }      
  #       if(any(grepl("brandom(", trmstrings, fixed = TRUE))){
  #         warning("Use brandomc() instead of brandom() with functional response, to get an identifiable model.")
  #       }      
  #       if(any(grepl("bspatial(", trmstrings, fixed = TRUE))){
  #         warning("Use bbsc() instead of bspatial() with functional response, to get an identifiable model.")
  #       }      
  #       if(any(grepl("bmono(", trmstrings, fixed = TRUE))){
  #         warning("With functional response, base-learner bmono() yields generally not an identifiable effect.")
  #       }      
  #       if(any(grepl("bmrf(", trmstrings, fixed = TRUE))){
  #         warning("With functional response, base-learner bmrf() yields generally not an identifiable effect.")
  #       }
  #     }    
  #   }
  
  ### extract time;
  # <FixMe> Only keep first variable, 
  # otherwise it is impossible to have a vector for knots
  yind <- all.vars(timeformula)[[1]]
  stopifnot(length(yind) == 1)
  nameyind <- yind
  assign(yind, data[[yind]])
  time <- data[[yind]]
  stopifnot(is.numeric(time))
  data[[yind]] <- NULL
  attr(time, "nameyind") <- nameyind
  
  ### extract covariates
  # data <- as.data.frame(data)
  allCovs <- unique(c(nameid, all.vars(formula)))
  if(length(allCovs) > 1){
    data <- data[allCovs[!allCovs %in% c(yname, nameyind)] ]
    if( any(is.na(names(data))) ) data <- data[ !is.na(names(data)) ]
  }else data <- list(NULL)  # <SB> intercept-model without covariates
        
  ### get covariates that are modeled constant over time
  # code of function pffr() of package refund
  tf <- terms.formula(formula, specials = c("c"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]], simplify = FALSE) 
  #ugly, but getTerms(formula)[-1] does not work for terms like I(x1:x2) 
  frmlenv <- environment(formula)
  where.c <- attr(tf, "specials")$c - 1    # indices of scalar offset terms
  
  #transform: c(foo) --> foo
  if(length(where.c)){ 
    trmstrings[where.c] <- sapply(trmstrings[where.c], function(x){
      sub("\\)$", "", sub("^c\\(", "", x)) #c(BLA) --> BLA
    })
    blconstant <- "bols(ONEtime, intercept = FALSE)"
  }
  assign("ONEtime", rep(1.0, length(time)))
    
  if(is.null(id)){
    ### check dimensions
    ### response has trajectories as rows
    stopifnot(is.matrix(response))
    # <SB> dataframe is list and can contain time-points of functional covariates of arbitrary length
    # if (nrow(data) > 0) stopifnot(nrow(response) == nrow(data))
    nr <- nrow(response)
    stopifnot(ncol(response) == length(time))
    nc <- ncol(response)
    dresponse <- as.vector(response) # column-wise stacking of response 
    nobs <- nr # number of observed trajectories
    ## check wether time variable is used in other base-learners
    ## only check in regular response case, as for irregular response, the problem cannot occur
    #if(nameyind %in% allCovs){
    #  warning("Do not use the same variable t as time-variable in y(t) and in the base-learners, e.g. as x(t).")
    #}
  }else{
    stopifnot(is.null(dim(response))) ## stopifnot(is.vector(response))
    # check length of response and its time and index
    stopifnot(length(response) == length(time) & length(response) == length(id))
    
    if(any(is.na(response))) warning("For non-grid observations the response should not contain missing values.")
    if( !all(sort(unique(id)) == 1:length(unique(id))) ) stop("id has to be integers 1, 2, 3,..., N.")
    
    nr <- length(response) # total number of observations
    nc <- length(unique(id)) # number of trajectories
    dresponse <- as.vector(response) # column-wise stacking of response  
    nobs <- length(unique(id)) # number of observed trajectories
  }
  
  ### save original dimensions of response
  ydim <- dim(response)
  
  ### variable to fit smooth intercept
  assign("ONEx", rep(1.0, nobs))
  
  #   ### FIXME: plausibility check: 
  #   # for constrained effects in the formula the model should include an intercept
  #   grep("bbsc", trmstrings) + grep("brandomc", trmstrings)
  
  
  ### compose mboost formula

  ## get formula over time
  tfm <- paste(deparse(timeformula), collapse = "")
  tfm <- strsplit(tfm, "~")[[1]]
  tfm <- strsplit(tfm[2], "\\+")[[1]]
  
  cfm <- paste(deparse(formula), collapse = "") 
  cfm <- strsplit(cfm, "~")[[1]]
  #xfm <- strsplit(cfm[2], "\\+")[[1]]
  xfm <- trmstrings
  ### replace "1" with intercept base learner
  if ( any( gsub(" ", "", strsplit(cfm[2], "\\+")[[1]]) ==  "1")){
    ## if(length(time) == 1) warning("Smooth intercept may not be sensitive for scalar response.")
    if( grepl("lambda", deparse(timeformula)) || 
         ( grepl("bols", deparse(timeformula)) &  !grepl("df", deparse(timeformula))) ){
      xfm <- c("bols(ONEx, intercept = FALSE, lambda = 0)", xfm)
    } else{
      xfm <- c("bols(ONEx, intercept = FALSE, df = 1)", xfm)
    }
    
    ## for FLAM model with %O% use anisotropic Kronecker product for not penalizing in direction of ONEx
    ## use %A0%, as smooth intercept has smooting parameter 0 in 1-direction 
    if(is.null(id)){
      xfm[[1]] <- paste(xfm[[1]], "%A0%", tfm)
    }
    
    where.c <- where.c + 1
  }
  
  ## check that the timevariable in timeformula and in the bhistx-base-learners have the same name
  if(any(grepl("bhistx", trmstrings))){
    for(j in 1:length(trmstrings)){
      if(any(grepl("bhistx", trmstrings[j]))){
        temp_name <- all.vars(formula(paste("~", trmstrings[[j]])[1]))[1]
        if(getTimeLab(data[[temp_name]]) != nameyind){
          stop("The timeLab of the hmatrix-object in bhistx(), '", getTimeLab(data[[temp_name]]),
               "', must be euqal to the name of the time-variable in timeformula, '", nameyind, "'.")
        }
      }
    }
  }
    
  yfm <- strsplit(cfm[1], "\\+")[[1]]
  
  ## set up formula for effects constant in time
  if(length(where.c) > 0){
    ### set c_df to the df/lambda in timeformula
    ##<FIXME> Does this make sense for bols() base-learner?
    if( grepl("lambda", tfm) || 
          ( grepl("bols", tfm) &  !grepl("df", tfm)) ){
      c_lambda <- eval(parse(text = paste(tfm, "$dpp(rep(1.0,", length(time), "))$df()", sep = "")))["lambda"]
      cfm <- paste("bols(ONEtime, intercept = FALSE, lambda = ", c_lambda ,")")
    } else{
      c_df <- eval(parse(text=paste(tfm, "$dpp(rep(1.0,", length(time), "))$df()", sep = "")))["df"]
      cfm <- paste("bols(ONEtime, intercept = FALSE, df = ", c_df ,")")
    }
  }

  # expand formula as Kronecker or tensor product 
  if(is.null(id)){
    tmp <- outer(xfm, tfm, function(x, y) paste(x, y, sep = "%O%"))
  }else{
    # expand the bl according to id
    if(grepl("ONEx", xfm[[1]])) xfm[[1]] <- paste(substr(xfm[[1]], 1 , nchar(xfm[[1]])-1), ", index=", nameid, ")", sep = "")
    xfm <- paste(substr(xfm, 1 , nchar(xfm)-1), ")", sep = "") # , index=id is done in the beginning
    tmp <- outer(xfm, tfm, function(x, y) paste(x, y, sep = "%X%"))
    ## <FIXME> use the following specification for irregular response -> adapt coef() and plot()-function 
    ## if(grepl("ONEx", tmp[[1]])) tmp[[1]] <- tfm ## smooth intercept is just time-formula
  }

  # do not expand an effect bconcurrent() or bhist() with timeformula
  if( length(c(grep("bconcurrent", tmp), grep("bhis", tmp)) ) > 0 ) 
    tmp[c(grep("bconcurrent", tmp), grep("bhist", tmp))] <- xfm[c(grep("bconcurrent", tmp), grep("bhist", tmp))]  
  
  ## do not expand effects in formula including %A% with timeformula
  if( length(grep("%A%", xfm)) > 0 ) 
    tmp[grep("%A%", xfm)] <- xfm[grep("%A%", xfm)]
  
  ## do not expand effects in formula including %A0% with timeformula
  if( length(grep("%A0%", xfm)) > 0 ) 
    tmp[grep("%A0%", xfm)] <- xfm[grep("%A0%", xfm)]
  
  ## do not expand effects in formula including %O% with timeformula
  if( length(grep("%O%", xfm)) > 0 ) 
    tmp[grep("%O%", xfm)] <- xfm[grep("%O%", xfm)]
  
  ## expand with a constant effect in t-direction 
  if(length(where.c) > 0){
    tmp[where.c] <- outer(xfm[where.c], cfm, function(x, y) paste(x, y, sep = "%O%"))
  } 
  
  ## for scalar response without FLAM-model do not use the Kronecker product
  if(scalarNoFLAM){
    tmp <- xfm
  } 
  
  
  ## put together the model formula 
  xpart <- paste(as.vector(tmp), collapse = " + ")
  fm <- as.formula(paste("dresponse ~ ", xpart))
  
  ### expand weights for observations
  if (is.null(weights)) weights <- rep(1, nr)
  w <- weights
  if(is.null(id)){
    if (length(w) == nr) w <- rep(w, nc) # expand weights if they are only on the columns
    if(length(w) != nc*nr) stop("Dimensions of weights do not match the dimensions of the response.") # check dimensions of w  
  }

  ## save the integration weights as data_weights
  ## per default the data_weights are all 1 
  data_weights <- 1
  
  ### multiply integration weights numInt to weights and w
  if(is.numeric(numInt)){
    if(length(numInt) != length(time)) stop("Length of integration weights and time vector are not equal.")
    weights <- weights * numInt
    data_weights <- numInt
    if(!is.null(ydim)){ ## only blow up for array model
      w <- rep(weights, each = nr)
      data_weights <- rep(data_weights, each = nr)
    }
  }else{
    if(!numInt %in% c("equal", "Riemann")) warning("argument numInt is ignored as it is neither numeric nor one of (\"equal\", \"Riemann\")")
    if(numInt == "Riemann"){ 
      data_weights <- as.vector(integrationWeights(X1 = response, time, id = id))
      w <- w * data_weights
    }
  }
  
  ### set weights of missing values to 0
  if(sum(is.na(dresponse)) > 0){
    #warning(paste("The response contains", sum(is.na(dresponse)) ,"missing values. The corresponding weights are set to 0."))
    w[which(is.na(dresponse))] <- 0
  }
  
  ### offset == "scalar", or offset = numeric of length 1, or scalar response
  ### -> use one scalar/user-specified offset like in mboost
  
  ## remember the offset-specification of FDboost
  offsetFDboost <- offset
  
  if( (!is.null(ydim) && ydim[2] == 1) || # scalar response
      !is.null(offset) && length(offset) == 1 ){  # offset == "scalar" / offset = numeric of length 1
    
    if( !is.null(offset) && !is.numeric(offset) && offset != "scalar" ){
      stop("User-specified offset must be numeric or 'scalar' to get a scalar offset as in mboost.")
    } 
    
    # use one constant offset in mboost(), as default in mboost
    if(!is.null(offset) && offset == "scalar"){
      offsetVec <- NULL
      predictOffset <- NULL
      offset <- NULL 
    }else{ # use user-specified offset of length 1 or offset = NULL for scalar response 
      offsetVec <- offset
      tempOffset <- if(length(offset) == 1) offset else NULL
      predictOffset <- function(time) tempOffset
    }

  ### specify time-specific offset 
  }else{
    
    ## offset for regular and irregular data: handling of missings is different!
    if(is.null(id)){
      
      ## per default add smooth time-specific offset 
      if(is.null(offset) && dim(response)[2] > 1 && 
         any(colMeans(response, na.rm = TRUE) > .Machine$double.eps *10^10)){
        message("Use a smooth offset.") 
        ### <FixMe> is the use of family@offset correct?
        #meanY <- colMeans(response, na.rm = TRUE)
        if(! "family" %in% names(dots) ){ # get the used family
          myfamily <-  Gaussian()
        } else myfamily <- dots$family
        offsetFun <- myfamily@offset
        meanY <- c()
        # do a linear interpolation of the response to prevent bias because of missing values
        # only use responses with less than 90% missings for the calculation of the offset
        # only use response curves whose weights are not completely 0 (important for resampling methods)
        meanNA <- apply(response, 1, function(x) mean(is.na(x)))
        responseInter <- t(apply(response[meanNA < 0.9 & rowSums(matrix(w, ncol = nc)) != 0 , ], 1, function(x) 
          approx(time, x, rule = offset_control$rule, xout = time)$y))
        # check whether first or last columns of the response contain solely NA
        # then use the values of the next column
        if(any(apply(responseInter, 2, function(x) all(is.na(x)) ) )){
          warning("Column of interpolated response contains nothing but NA.")
          allNA <- apply(responseInter, 2, function(x) all(is.na(x)) ) 
          allNAlower <- allNAupper <- allNA
          allNAupper[1:round(ncol(responseInter)/2)] <- FALSE # missing columns with low index 
          allNAlower[round(ncol(responseInter)/2):ncol(responseInter)] <- FALSE # missing columns with high index 
          responseInter[, allNAlower] <- responseInter[, max(which(allNAlower))+1]
          responseInter[, allNAupper] <- responseInter[, min(which(allNAlower))-1]
        }
        
        for(i in 1:nc){
          try(meanY[i] <- offsetFun(responseInter[,i], 1*!is.na(responseInter[,i]) ), silent = TRUE)
        }
        # meanY <- sapply(1:nc, function(i) offsetFun(responseInter[,i], 1*!is.na(responseInter[,i])))
        
        if( is.null(meanY) ||  any(is.na(meanY)) ){
          warning("Mean offset cannot be computed by family@offset(). Use a weighted mean instead.")
          meanY <- c()
          for(i in 1:nc){
            meanY[i] <- Gaussian()@offset(responseInter[,i], 1*!is.na(responseInter[,i]) )
          }
        }
        rm(responseInter, meanNA)
        ### <FixMe> is the computation of k ok? 
        if(!offset_control$cyclic){
          modOffset <- try( gam(meanY ~ s(time, bs = "ad", 
                                          k = min(offset_control$k_min, round(length(time)/2))  ),
                                knots = offset_control$knots), 
                            silent = offset_control$silent )
        }else{ # use cyclic splines
          modOffset <- try( gam(meanY ~ s(time, bs = "cc", 
                                          k = min(offset_control$k_min, round(length(time)/2))  ),
                                knots = offset_control$knots), 
                            silent = offset_control$silent )
        }
        
        if(any(class(modOffset) == "try-error")){
          warning(paste("Could not fit the smooth offset by adaptive splines (default), use a simple spline expansion with 5 df instead.",
                        if(offset_control$cyclic) "This offset is not cyclic!"))
          if(round(length(time)/2) < 8) warning("Most likely because of too few time-points.")
          modOffset <- lm(meanY ~ bs(time, df = 5))
        } 
        offsetVec <- modOffset$fitted.values 
        predictOffset <- function(time){
          ret <- as.numeric(predict(modOffset, newdata = data.frame(time = time)))
          names(ret) <- NULL
          ret      
        } 
        offset <- as.vector(matrix(offsetVec, ncol = ncol(response), nrow = nrow(response), byrow = TRUE))    
      }else{ 
        ### scalar response or mean-centered response -> one constant offset value is used
        if(dim(response)[2] == 1 |  all(colMeans(response, na.rm = TRUE) < .Machine$double.eps *10^10)){ 
          offsetVec <- offset
          predictOffset <- offset 
        }else{ 
          # expand the offset to the long vector like dresponse
          if(length(offset) != nc) stop("Dimensions of offset and response do not match.")
          offsetVec <- offset
          offset <- as.vector(matrix(offset, ncol = ncol(response), nrow = nrow(response), byrow = TRUE))
          ### <FixMe> use a more sophisticated model to estimate the time-specific offset? 
          modOffset <- lm(offsetVec ~ bs(time, df = length(offsetVec)-2))
          predictOffset <- function(time){
            ret <- as.numeric(predict(modOffset, newdata = data.frame(time = time)))
            names(ret) <- NULL
            ret      
          }      
        }
      }  ## end else{} for no smooth time-specific offset for regular response
      
      # for irregular data the model for the smooth offset is computed on the available data
    }else{
      
      ## compute a time-specific smooth offset for irregular data 
      if(is.null(offset)){
        # only use response curves whose weights are not completely 0 (important for resampling methods)
        # do this by setting responses to NA whose weight is zero
        responseW <- response
        responseW[w == 0] <- NA
        message("Use a smooth offset for irregular data.") 
        ### <FixMe> is the computation of k ok? 
        if(!offset_control$cyclic){
          modOffset <- try( gam(responseW ~ s(time, bs = "ad", 
                                              k = min(offset_control$k_min, round(length(time)/10))  ),
                                knots = offset_control$knots), 
                            silent = offset_control$silent )
        }else{ # use cyclic splines
          modOffset <- try( gam(responseW ~ s(time, bs = "cc", 
                                              k = min(offset_control$k_min, round(length(time)/10))  ),
                                knots = offset_control$knots), 
                            silent = offset_control$silent )
        }
        
        if(any(class(modOffset) == "try-error")){
          warning(paste("Could not fit the smooth offset by adaptive splines (default), use a simple spline expansion with 5 df instead.",
                        if(offset_control$cyclic) "This offset is not cyclic!"))
          if(round(length(time)/2) < 8) warning("Most likely because of too few time-points.")
          modOffset <- lm(response ~ bs(time, df = 5))
        } 
        
        offsetVec <- as.numeric(predict(modOffset, newdata = data.frame(time = time)))
        predictOffset <- function(time){
          ret <- as.numeric(predict(modOffset, newdata = data.frame(time = time)))
          names(ret) <- NULL
          ret      
        } 
        offset <- offsetVec
        
        # no time-specific offset -> constant offset is estimated within mboost()
      }else{
        stop("User specified offset must be of length 1 for irregularly observed response.") 
      }
      

      
    } # end else from if(is.null(id))
    
  }
    
  offsetMboost <- offset
  
  if (length(data) > 0) {
    ### mboost isn't happy with nrow(data) == 0
    ret <- mboost(fm, data = data, weights = w, offset = offset, ...) 
  } else {
    ret <- mboost(fm, weights = w, offset = offset, ...)
  }
  
  # check sum-to-zero constraints for the fitted effects
  # for models with more than one effect and a regular response
  # not for scalar response
  if(check0 && length(ret$baselearner) > 1 && is.null(id) && dim(response)[2] != 1){
    
    # do not check the smooth intercept
    if(any( gsub(" ", "", strsplit(cfm[2], "\\+")[[1]]) ==  "1")){
      effectsToCheck <- 2:length(ret$baselearner)
    }else{
      effectsToCheck <- 1:length(ret$baselearner)
    }
    # predict each effect separately
    pred <- predict(ret, which = effectsToCheck)
    
    # check weather each effect is zero per time-point
    meanPerTime <- apply(pred, 2, function(x){
      tapply(x, rep(1:nc, each = nobs), mean) # compute mean per time-point
    })
    if(!all(meanPerTime < .Machine$double.eps *10^10)){
      temp <- colSums(meanPerTime > .Machine$double.eps *10^10) != 0
      message("The effects ", names(temp)[temp], " do not sum to zero per time-point.")
    }
  }
  
  ### assign new class (e.g. for specialized predictions)
  class(ret) <- c("FDboost", class(ret))
  if(!is.null(id)) class(ret) <- c("FDboostLong", class(ret))
  if(scalarResponse) class(ret) <- c("FDboostScalar", class(ret))
  ## generate an id-variable for a regular response 
  if(is.null(id)) id <- rep(1:ydim[1], times = ydim[2])
  
  ### reset weights for cvrisk etc., expanding works OK in bl_lin_matrix!
  # ret$"(weights)" <- weights
  # <SB> do not reset weights as than the integration weights are lost
  
  ret$yname <- yname
  ret$ydim <- ydim  
  ret$yind <- time
  ret$data <- data
  ret$id <- id
  attr(ret$id, "nameid") <- nameid
  
  # if the offset is just an integer the prediction gives back this integer
  ret$predictOffset <- predictOffset
  if(is.null(offset)) ret$predictOffset <- function(time) ret$offset
  
  # <FIXME> do no longer use offsetVec
  # offsetVec is an integer if no smooth offset was calculated
  ret$offsetVec <- offsetVec
  if(is.null(offset)) ret$offsetVec <- ret$offset
  ret$offsetFDboost <- offsetFDboost # offset as specified in call to FDboost 
  ret$offsetMboost <- offsetMboost # offset as given to mboost 
      
  # save the call
  ret$call <- match.call()
  
  # save the evaluated call
  ret$callEval <- ret$call
  ret$callEval[-1] <- lapply(ret$call[-1], function(x){  
    eval(x, parent.frame(3)) # use the environment of the call to FDboost()
  })
  ret$callEval$data <- NULL # do not save data and weights in callEval to save memory
  ret$callEval$weights <- NULL
  
  # save formulas as character strings to save memory
  ret$timeformula <- paste(deparse(timeformula), collapse = "")
  if(scalarNoFLAM) ret$timeformula <- ""
  ret$formulaFDboost <- paste(deparse(formulaFDboost), collapse = "")
  ret$formulaMboost <- paste(deparse(fm), collapse = "")
  
#   ret$timeformula <- timeformula
#   ret$formulaFDboost <- formulaFDboost
#   ret$formulaMboost <- fm
  
  ret
}

