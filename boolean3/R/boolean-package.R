## ----------------------------------------------------------------------------
## This file is part of boolean3
##
## Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
##
## boolean3 represents a substantial re-write of the original boolean package
## developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
## was developed under the direction of Bear Braumoeller and with support from
## The Ohio State University's College of Social and Behavioral Sciences.
##
## boolean3 and is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program.  If not, see <http://www.gnu.org/licenses/>.
##
## ----------------------------------------------------------------------------

##' @import stats rgenoud numDeriv optimx parallel mvtnorm rlecuyer lattice
NULL

##' 
##' Boolean binary response models are a family of partial-observability binary
##' response models designed to permit researchers to model causal complexity,
##' or multiple causal ``paths'' to a given outcome.
##'
##' Boolean permits estimation of Boolean binary response models (see
##' Braumoeller 2003 for derivation), which are a family of
##' partial-observability n-variate models designed to permit researchers to
##' model causal complexity, or multiple causal ``paths'' to a given
##' outcome. The various ``paths'' are modeled as latent dependent variables
##' that are multiplied together in a manner determined by the logic of their
##' (Boolean) interaction. If, for example, we wanted to model a situation in
##' which diet OR smoking causes heart failure, we would use one set of
##' independent variables (caloric intake, fat intake, etc.) to predict the
##' latent probability of diet-related coronary failure (\eqn{y1^*}), use
##' another set of variables (cigarettes smoked per day, exposure to second-hand
##' smoke, etc.)  to predict the latent probability of smoking-related coronary
##' failure (\eqn{y2^*}), and model the observed outcome (\eqn{y}, or coronary
##' failure) as a function of the Boolean interaction of the two: \eqn{\Pr(y=1) =
##' 1-([1-y1^*] \times [1-y2^*])}. Independent variables that have an impact on both
##' latent dependent variables can be included in both paths. Any combination of
##' any number of ANDs and ORs can be modeled, though the procedure becomes
##' exponentially more data-intensive as the number of interactions increases.
##' 
##' \tabular{ll}{
##' Package: \tab boolean3 \cr
##' Type: \tab Package\cr
##' Version: \tab 3.1.6\cr
##' Date: \tab 2014-11-15\cr
##' License: \tab GPL-3 \cr
##' LazyLoad: \tab yes\cr
##' }
##' 
##' @name boolean3-package
##' @aliases boolean3
##' @docType package
##' @title Modeling Causal Complexity
##' @note \code{boolean3} was developed by Jason W. Morgan under the direction
##' of Bear Braumoeller with support from The Ohio State University's College of
##' Social and Behavioral Sciences. The package represents a significant
##' re-write of the original \code{boolean} implementation developed by Bear
##' Braumoeller, Ben Goodrich, and Jacob Kline. Please see the release notes and
##' accompanying documentation for details regarding the many changes made in
##' this version.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @references Braumoeller, Bear F. (2003) ``Causal Complexity and the Study of
##' Politics.'' \emph{Political Analysis} 11(3): 209--233.
##' @keywords package
##' @seealso See \code{\link{boolprep}} for model setup and
##' \code{\link{boolean}} for estimation.
##' @examples
##' ## Generate some fake data
##' require(mvtnorm)
##' set.seed(12345)
##' N  <- 2000
##' Df <- cbind(1, rmvnorm(N, mean=rep(0, 5)))
##'
##' ## Set coefficients
##' beta.a <- c(-2.00, 0.33, 0.66, 1.00)
##' beta.b <- c(0.00, 1.50, -0.25)
##'
##' ## Generate path probabilities following a normal model.
##' y.a <- as.vector(pnorm(tcrossprod(beta.a, Df[, 1:4])))
##' y.b <- as.vector(pnorm(tcrossprod(beta.b, Df[, c(1, 5, 6)])))
##'
##' ## AND and OR-model
##' or <- function(x, y) { x + y - x * y }
##' and <- function(x, y) { x * y }
##' y.star.OR  <- or(y.a, y.b)
##' y.star.AND <- and(y.a, y.b)
##'
##' ## Observed responses
##' y.OR <- rbinom(N, 1, y.star.OR)
##' y.AND <- rbinom(N, 1, y.star.AND)
##'
##' ## Set up data.frame for estimation
##' Df <- cbind(1, Df)
##' Df <- as.data.frame(Df)
##' Df[,1] <- y.OR
##' Df[,2] <- y.AND
##' names(Df) <- c("y.OR", "y.AND", "x1", "x2", "x3", "x4", "x5")
##' 
##' ## Before estimating, boolean models need to be specified using the
##' ## boolprep function.
##' 
##' ## OR model, specified to use a probit link for each causal path. This
##' ## matches the data generating process above.
##' mod.OR <- boolprep(y.OR ~ (a | b), a ~ x1 + x2 + x3, b ~ x4 + x5,
##'                    data = Df, family=list(binomial("probit")))
##' 
##' ## IF you really want to, it's also possible to specify a different 
##' ## link function for each causal path. These should be in the same 
##' ## order as specified in the model formula.
##' mod.OR2 <- boolprep(y.OR ~ (a | b), a ~ x1 + x2 + x3, b ~ x4 + x5,
##'                     data = Df, family=list(binomial("probit"),
##'                                    binomial("logit")))
##' 
##' ## Fit the prepared model using the nlminb optimizer (the default).
##' (fit.OR <- boolean(mod.OR, method="nlminb", control=list(trace=1)))
##' 
##' ## Multiple optimizers can be specified in a single call to boolean. 
##' ## Here we fit with the nlm and nlminb optimizers.
##' (fit1.OR <- boolean(mod.OR, method=c("nlm", "nlminb")))
##' 
##' ## Re-fit, with BFGS and a higher maximum number of iterations. All 
##' ## of the options that go along with nlm(), optim(), and genoud() should 
##' ## be transparently passable through boolean.
##' (fit2.OR <- boolean(mod.OR, method="BFGS", control = list(maxit = 500)))
##' 
##' ## Induce a convergence warning message.
##' (fit3.OR <- boolean(mod.OR, method="BFGS", control = list(maxit = 5)))
##'
##' \dontrun{
##' ## Now estimate model with genoud, a genetic optimizer. This has the 
##' ## capability of using multiple processors via parallel.
##' (fit6.OR <- boolean(mod.OR, method="genoud",
##'                     cluster=c("localhost", "localhost"),
##'                     print.level=2))
##' 
##' ## The default SANN optimizer is also available.
##' (fit7.OR <- boolean(mod.OR, method="SANN"))
##' }
##' 
##' ## The fit is stored as "model.fit", within the boolean object.
##' str(fit.OR$model.fit)
##' 
##' ## Create a summary object, saving and printing it. Then take a look at 
##' ## the objects stored in the summary object.
##' (smry <- summary(fit.OR))
##' str(smry)
##' 
##' ## Extract log-likelihood and coefficient vector.
##' logLik(fit.OR)
##' coef(fit.OR)
##' 
##' \dontrun{
##' ## Display the contours of the likelihood given a change the value of 
##' ## the coefficients. Despite the function name, these are not true 
##' ## profile likelihoods as they hold all other coefficients fixed at
##' ## their MLE.
##' (prof <- boolprof(fit.OR))
##' 
##' ## Extract the plots for x1_a and x4_b.
##' plot(prof, y = c("x1_a", "x4_b"))
##' plot(prof, y = c(1, 3), scales = list(y = list(relation = "free")))
##' 
##' ## You can also use variable or index matching with boolprof to select 
##' ## particular covariates of interest.
##' boolprof(fit.OR, vars = c(1, 3))
##' boolprof(fit.OR, vars = c("x1_a", "x4_b"))
##' 
##' ## Plots of predicted probabilities are available through boolprob.
##' ## With boolprob, either vars or newdata *must* be specified.
##' boolprob(fit.OR, vars = c("x1_a", "x4_b"))
##' boolprob(fit.OR, vars = c(2, 3, 4, 6))
##' 
##' ## Specifying conf.int = TRUE produces simulated confidence intervals. 
##' ## The number of samples to pull from the distribution of the estimated
##' ## coefficients is controlled by n; n=100 is default. This can take a
##' ## while.
##' (prob <- boolprob(fit.OR, vars = c(2, 3, 4, 6), n = 1000,
##'                   conf.int = TRUE))
##' 
##' ## Choose a different method estimate upon which to base the estimates.
##' (prob <- boolprob(fit1.OR, method="nlm", vars=c(2, 3, 4, 6), n=1000,
##'                   conf.int=TRUE))
##' 
##' ## As with the other components of the model, you can extract the 
##' ## predicted probabilities.
##' str(prob)
##' prob$est
##'
##' ## Bootstrapping is also possible, and is the recommended method of 
##' ## making inferences. The boolbool function uses a simple sampling scheme: 
##' ## it resamples repeatedly from the provided data.frame. More the complex
##' ## data structures with, for example, clustering, need to be dealt with
##' ## manually.
##' (bs <- boolboot(fit, n=10))
##'
##' ## boolboot supports bootstrapping across multiple processors.
##' (bs <- boolboot(fit, n=100, cluster=c("localhost", "localhost")))
##' 
##' }
NULL

