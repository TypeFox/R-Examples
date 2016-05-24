#' Phylogenetic Generalised Linear Mixed Model for Community Data
#'
#' This function performs Generalized Linear Mixed Models for binary
#' and continuous phylogenetic data, estimating regression
#' coefficients with approximate standard errors. It is modeled after
#' \code{\link[lme4:lmer]{lmer}} but is more general by allowing
#' correlation structure within random effects; these correlations can
#' be phylogenetic among species, or any other correlation structure,
#' such as geographical correlations among sites. It is, however, much
#' more specific than \code{\link[lme4:lmer]{lmer}} in that it can
#' only analyze a subset of1 the types of model designed handled by
#' \code{\link[lme4:lmer]{lmer}}. It is also much slower than
#' \code{\link[lme4:lmer]{lmer}} and requires users to specify
#' correlation structures as covariance
#' matrices. \code{communityPGLMM} can analyze models in Ives and
#' Helmus (2011). It can also analyze bipartite phylogenetic data,
#' such as that analyzed in Rafferty and Ives (2011), by giving sites
#' phylogenetic correlations.
#' @param formula a two-sided linear formula object describing the
#' fixed-effects of the model; for example, \code{Y ~ X}.
#' @param data a \code{\link{data.frame}} containing the variables
#' named in formula. The data frame should have long format with
#' factors specifying species and sites. \code{communityPGLMM} will
#' reorder rows of the data frame so that species are nested within
#' sites. Please note that calling
#' \code{\link{as.data.frame.comparative.comm}} will return your
#' \code{comparative.comm} object into this format for you.
#' @param family either \code{gaussian} for a Linear Mixed Model, or
#' \code{binomial} for binary dependent data.
#' @param sp a \code{\link{factor}} variable that identifies species
#' @param site a \code{\link{factor}} variable that identifies sites
#' @param random.effects a \code{\link{list}} that contains, for
#' non-nested random effects, lists of triplets of the form
#' \code{list(X, group = group, covar = V)}. This is modeled after the
#' \code{\link[lme4:lmer]{lmer}} formula syntax \code{(X | group)}
#' where \code{X} is a variable and \code{group} is a grouping
#' factor. Note that \code{group} should be either your \code{sp} or
#' \code{site} variable specified in \code{sp} and \code{site}. The
#' additional term \code{V} is a covariance matrix of rank equal to
#' the number of levels of group that specifies the covariances among
#' groups in the random effect \code{X}. For nested variable random
#' effects, \code{random.effects} contains lists of quadruplets of the
#' form \code{list(X, group1 = group1, covar = V, group2 = group2)}
#' where \code{group1} is nested within \code{group2}.
#' @param REML whether REML or ML is used for model fitting. For the
#' generalized linear mixed model for binary data, these don't have
#' standard interpretations, and there is no log likelihood function
#' that can be used in likelihood ratio tests.
#' @param s2.init an array of initial estimates of s2 for each random
#' effect that scales the variance. If s2.init is not provided for
#' \code{family="gaussian"}, these are estimated using in a clunky way
#' using \code{\link{lm}} assuming no phylogenetic signal.  A better
#' approach is to run \code{link[lme4:lmer]{lmer}} and use the output
#' random effects for \code{s2.init}. If \code{s2.init} is not
#' provided for \code{family="binomial"}, these are set to 0.25.
#' @param B.init initial estimates of \eqn{B}{B}, a matrix containing
#' regression coefficients in the model for the fixed effects. This
#' matrix must have \code{dim(B.init)=c(p+1,1)}, where \code{p} is the
#' number of predictor (independent) variables; the first element of
#' \code{B} corresponds to the intercept, and the remaining elements
#' correspond in order to the predictor (independent) variables in the
#' formula.  If \code{B.init} is not provided, these are estimated
#' using in a clunky way using \code{\link{lm}} or \code{\link{glm}}
#' assuming no phylogenetic signal.  A better approach is to run
#' \code{\link[lme4:lmer]{lmer}} and use the output fixed effects for
#' \code{B.init}.
#' @param reltol a control parameter dictating the relative tolerance
#' for convergence in the optimization; see \code{\link{optim}}.
#' @param maxit a control parameter dictating the maximum number of
#' iterations in the optimization; see \code{\link{optim}}.
#' @param tol.pql a control parameter dictating the tolerance for
#' convergence in the PQL estimates of the mean components of the
#' binomial GLMM.
#' @param maxit.pql a control parameter dictating the maximum number
#' of iterations in the PQL estimates of the mean components of the
#' binomial GLMM.
#' @param verbose if \code{TRUE}, the model deviance and running
#' estimates of \code{s2} and \code{B} are plotted each iteration
#' during optimization.
#' @param ... additional arguments to summary and plotting functions
#' (currently ignored)
#' @details The vignette 'pez-pglmm-overview' gives a gentle
#' introduction to using PGLMMS. For linear mixed models (\code{family
#' = 'gaussian'}), the function estimates parameters for the model of
#' the form, for example,
#' 
#' \deqn{Y = \beta_0 + \beta_1x + b_0 + b_1x}{y = beta_0 + beta_1x + b_0 + b_1x}
#' \deqn{b_0 ~ Gaussian(0, \sigma_0^2I_{sp})}{b_0 ~ Gaussian(0, sigma_0^2I_(sp))}
#' \deqn{b_1 ~ Gaussian(0, \sigma_0^2V_{sp})}{b_0 ~ Gaussian(0, sigma_0^2V_(sp))}
#' \deqn{\eta ~ Gaussian(0,\sigma^2)}{e ~ Gaussian(0,sigma^2)}
#' 
#' where \eqn{\beta_0}{beta_0} and \eqn{\beta_1}{beta_1} are fixed
#' effects, and \eqn{V_{sp}}{V_(sp)} is a variance-covariance matrix
#' derived from a phylogeny (typically under the assumption of
#' Brownian motion evolution). Here, the variation in the mean
#' (intercept) for each species is given by the random effect
#' \eqn{b_0}{b_0} that is assumed to be independent among
#' species. Variation in species' responses to predictor variable
#' \eqn{x}{x} is given by a random effect \eqn{b_0}{b_0} that is
#' assumed to depend on the phylogenetic relatedness among species
#' given by \eqn{V_{sp}}{V_(sp_}; if species are closely related,
#' their specific responses to \eqn{x}{x} will be similar. This
#' particular model would be specified as
#'
#' \code{re.1 <- list(1, sp = dat$sp, covar = diag(nspp))}
#' \code{re.2 <- list(dat$X, sp = dat$sp, covar = Vsp)}
#' \code{z <- communityPGLMM(Y ~ X, data = data, family = "gaussian", random.effects = list(re.1, re.2))}
#' 
#' The covariance matrix covar is standardized to have its determinant
#' equal to 1. This in effect standardizes the interpretation of the
#' scalar \eqn{\sigma^2}{sigma^2}. Although mathematically this is
#' not required, it is a very good idea to standardize the predictor
#' (independent) variables to have mean 0 and variance 1. This will
#' make the function more robust and improve the interpretation of the
#' regression coefficients. For categorical (factor) predictor
#' variables, you will need to construct 0-1 dummy variables, and
#' these should not be standardized (for obvious reasons).
#'
#' For binary generalized linear mixed models (\code{family =
#' 'binomial'}), the function estimates parameters for the model of
#' the form, for example,
#'
#' \deqn{y = \beta_0 + \beta_1x + b_0 + b_1x}{y = beta_0 + beta_1x + b_0 + b_1x}
#' \deqn{Y = logit^{-1}(y)}{Y = logit^(-1)(y)}
#' \deqn{b_0 ~ Gaussian(0, \sigma_0^2I_{sp})}{b_0 ~ Gaussian(0, sigma_0^2I_(sp))}
#' \deqn{b_1 ~ Gaussian(0, \sigma_0^2V_{sp})}{b_0 ~ Gaussian(0, sigma_0^2V_(sp))}
#'
#' where \eqn{\beta_0}{beta_0} and \eqn{\beta_1}{beta_1} are fixed
#' effects, and \eqn{V_{sp}}{V_(sp)} is a variance-covariance matrix
#' derived from a phylogeny (typically under the assumption of
#' Brownian motion evolution).
#' 
#' \code{z <- communityPGLMM(Y ~ X, data = data, family =
#' 'binomial', random.effects = list(re.1, re.2))}
#' 
#' As with the linear mixed model, it is a very good idea to
#' standardize the predictor (independent) variables to have mean 0
#' and variance 1. This will make the function more robust and improve
#' the interpretation of the regression coefficients. For categorical
#' (factor) predictor variables, you will need to construct 0-1 dummy
#' variables, and these should not be standardized (for obvious
#' reasons).
#' @return an object of class \code{communityPGLMM}
#' \item{formula}{the formula for fixed effects}
#' \item{data}{the dataset}
#' \item{family}{either \code{gaussian} or \code{binomial} depending on the model fit}
#' \item{random.effects}{the list of random effects}
#' \item{B}{estimates of the regression coefficients}
#' \item{B.se}{approximate standard errors of the fixed effects regression coefficients}
#' \item{B.cov}{approximate covariance matrix for the fixed effects regression coefficients}
#' \item{B.zscore}{approximate Z scores for the fixed effects regression coefficients}
#' \item{B.pvalue}{approximate tests for the fixed effects regression coefficients being different from zero}
#' \item{ss}{random effects' standard deviations for the covariance matrix \eqn{\sigma^2V}{sigma^2 V} for each random effect in order. For the linear mixed model, the residual variance is listed last}
#' \item{s2r}{random effects variances for non-nested random effects}
#' \item{s2n}{random effects variances for nested random effects}
#' \item{s2resid}{for linear mixed models, the residual vairance}
#' \item{logLIK}{for linear mixed models, the log-likelihood for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalised linear mixed models}
#' \item{AIC}{for linear mixed models, the AIC for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalised linear mixed models}
#' \item{BIC}{for linear mixed models, the BIC for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalised linear mixed models}
#' \item{REML}{whether or not REML is used (\code{TRUE} or \code{FALSE})}
#' \item{s2.init}{the user-provided initial estimates of \code{s2}}
#' \item{B.init}{the user-provided initial estimates of \code{B}}
#' \item{Y}{the response (dependent) variable returned in matrix form}
#' \item{X}{the predictor (independent) variables returned in matrix form (including 1s in the first column)}
#' \item{H}{the residuals. For the generalized linear mixed model, these are the predicted residuals in the \eqn{logit^{-1}}{logit -1} space}
#' \item{iV}{the inverse of the covariance matrix for the entire system (of dimension (nsp*nsite) by (nsp*nsite))}
#' \item{mu}{predicted mean values for the generalized linear mixed model. Set to NULL for linear mixed models}
#' \item{sp, sp}{matrices used to construct the nested design matrix}
#' \item{Zt}{the design matrix for random effects}
#' \item{St}{diagonal matrix that maps the random effects variances onto the design matrix}
#' \item{convcode}{the convergence code provided by \code{\link{optim}}}
#' \item{niter}{number of iterations performed by \code{\link{optim}}}
#' @note These function \emph{do not} use a
#' \code{\link{comparative.comm}} object, but you can use
#' \code{\link{as.data.frame.comparative.comm}} to
#' create a \code{data.frame} for use with these functions. The power
#' of this method comes from deciding your own parameters parameters
#' to be determined (the data for regression, the random effects,
#' etc.), and it is our hope that this interface gives you more
#' flexibility in model selection/fitting.
#' @author Anthony R. Ives, cosmetic changes by Will Pearse
#' @references Ives, A. R. and M. R. Helmus. 2011. Generalized linear
#' mixed models for phylogenetic analyses of community
#' structure. Ecological Monographs 81:511-525.
#' @references Rafferty, N. E., and A. R. Ives. 2013. Phylogenetic
#' trait-based analyses of ecological networks. Ecology 94:2321-2333.
#' @rdname pglmm
#' @name pglmm
#' @aliases communityPGLMM
#' @export
#' @examples
#' ## Structure of examples:
#' # First, a (brief) description of model types, and how they are specified
#' # - these are *not* to be run 'as-is'; they show how models should be organised
#' # Second, a run-through of how to simulate, and then analyse, data
#' # - these *are* to be run 'as-is'; they show how to format and work with data
#'
#' \dontrun{
#' #########################################################
#' #First section; brief summary of models and their use####
#' #########################################################
#' ## Model structures from Ives & Helmus (2011)
#' # dat = data set for regression (note: *not* an comparative.comm object)
#' # nspp = number of species
#' # nsite = number of sites
#' # Vphy = phylogenetic covariance matrix for species
#' # Vrepul = inverse of Vphy representing phylogenetic repulsion
#'
#' # Model 1 (Eq. 1)
#' re.site <- list(1, site = dat$site, covar = diag(nsite))
#' re.sp.site <- list(1, sp = dat$sp, covar = Vphy, site = dat$site) # note: nested
#' z <- communityPGLMM(freq ~ sp, data = dat, family = "binomial", sp
#' = dat$sp, site = dat$site, random.effects = list(re.site,
#' re.sp.site), REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' 
#' # Model 2 (Eq. 2)
#' re.site <- list(1, site = dat$site, covar = diag(nsite))
#' re.slope <- list(X, sp = dat$sp, covar = diag(nspp))
#' re.slopephy <- list(X, sp = dat$sp, covar = Vphy)
#' z <- communityPGLMM(freq ~ sp + X, data = dat, family = "binomial",
#' sp = dat$sp, site = dat$site, random.effects = list(re.site,
#' re.slope, re.slopephy), REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' # Model 3 (Eq. 3)
#' re.site <- list(1, site = dat$site, covar = diag(nsite))
#' re.sp.site <- list(1, sp = dat$sp, covar = Vrepul, site = dat$site) # note: nested
#' z <- communityPGLMM(freq ~ sp*X, data = dat, family = "binomial",
#' sp = dat$sp, site = dat$site, random.effects = list(re.site,
#' re.sp.site), REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' ## Model structure from Rafferty & Ives (2013) (Eq. 3)
#' # dat = data set
#' # npp = number of pollinators (sp)
#' # nsite = number of plants (site)
#' # VphyPol = phylogenetic covariance matrix for pollinators
#' # VphyPlt = phylogenetic covariance matrix for plants
#' 
#' re.a <- list(1, sp = dat$sp, covar = diag(nspp))
#' re.b <- list(1, sp = dat$sp, covar = VphyPol)
#' re.c <- list(1, sp = dat$sp, covar = VphyPol, dat$site)
#' re.d <- list(1, site = dat$site, covar = diag(nsite))
#' re.f <- list(1, site = dat$site, covar = VphyPlt)
#' re.g <- list(1, site = dat$site, covar = VphyPlt, dat$sp)
#' #term h isn't possible in this implementation, but can be done with
#' available matlab code
#' 
#' z <- communityPGLMM(freq ~ sp*X, data = dat, family = "binomial",
#' sp = dat$sp, site = dat$site, random.effects = list(re.a, re.b,
#' re.c, re.d, re.f, re.g), REML = TRUE, verbose = TRUE, s2.init=.1)
#' }
#' 
#' #########################################################
#' #Second section; detailed simulation and analysis #######
#' #NOTE: this section is explained and annotated in #######
#' #      detail in the vignette 'pez-pglmm-overview'#######
#' #      run 'vignette('pez-pglmm-overview') to read#######
#' #########################################################
#' # Generate simulated data for nspp species and nsite sites
#' nspp <- 15
#' nsite <- 10
#' 
#' # residual variance (set to zero for binary data)
#' sd.resid <- 0
#' 
#' # fixed effects
#' beta0 <- 0
#' beta1 <- 0
#' 
#' # magnitude of random effects
#' sd.B0 <- 1
#' sd.B1 <- 1
#' 
#' # whether or not to include phylogenetic signal in B0 and B1
#' signal.B0 <- TRUE
#' signal.B1 <- TRUE
#' 
#' # simulate a phylogenetic tree
#' phy <- rtree(n = nspp)
#' phy <- compute.brlen(phy, method = "Grafen", power = 0.5)
#' 
#' # standardize the phylogenetic covariance matrix to have determinant 1
#' Vphy <- vcv(phy)
#' Vphy <- Vphy/(det(Vphy)^(1/nspp))
#' 
#' # Generate environmental site variable
#' X <- matrix(1:nsite, nrow = 1, ncol = nsite)
#' X <- (X - mean(X))/sd(X)
#' 
#' # Perform a Cholesky decomposition of Vphy. This is used to
#' # generate phylogenetic signal: a vector of independent normal random
#' # variables, when multiplied by the transpose of the Cholesky
#' # deposition of Vphy will have covariance matrix equal to Vphy.
#' 
#' iD <- t(chol(Vphy))
#' 
#' # Set up species-specific regression coefficients as random effects
#' if (signal.B0 == TRUE) {
#' 		b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0)
#' } else {
#' 		b0 <- beta0 + rnorm(nspp, sd = sd.B0)
#' }
#' if (signal.B1 == TRUE) {
#' 		b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1)
#' } else {
#' 		b1 <- beta1 + rnorm(nspp, sd = sd.B1)
#' }
#' 
#' # Simulate species abundances among sites to give matrix Y that
#' # contains species in rows and sites in columns
#' y <- rep(b0, each=nsite)
#' y <- y + rep(b1, each=nsite) * rep(X, nspp)
#' y <- y + rnorm(nspp*nsite) #add some random 'error'
#' Y <- rbinom(length(y), size=1, prob=exp(y)/(1+exp(y)))

#' y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
#' ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
#' e <- rnorm(nspp * nsite, sd = sd.resid)
#' y <- y + matrix(e, nrow = nspp, ncol = nsite)
#' y <- matrix(y, nrow = nspp * nsite, ncol = 1)
#' 
#' Y <- rbinom(n = length(y), size = 1, prob = exp(y)/(1 + exp(y)))
#' Y <- matrix(Y, nrow = nspp, ncol = nsite)
#' 
#' # name the simulated species 1:nspp and sites 1:nsites
#' rownames(Y) <- 1:nspp
#' colnames(Y) <- 1:nsite
#'
#' par(mfrow = c(3, 1), las = 1, mar = c(2, 4, 2, 2) - 0.1)
#' matplot(t(X), type = "l", ylab = "X", main = "X among sites")
#' hist(b0, xlab = "b0", main = "b0 among species")
#' hist(b1, xlab = "b1", main = "b1 among species")
#'
#' #Plot out; you get essentially this from plot(your.pglmm.model)
#' image(t(Y), ylab = "species", xlab = "sites", main = "abundance",
#' col=c("black","white"))
#' 
#' # Transform data matrices into "long" form, and generate a data frame
#' YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
#'
#' XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
#' nspp * nsite, ncol = 1)
#' 
#' site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
#' 1)), nrow = nspp * nsite, ncol = 1)
#' sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp),
#' nrow = nspp * nsite, ncol = 1)
#' 
#' dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))
#' 
#' # Format input and perform communityPGLMM()
#' # set up random effects
#' 
#' # random intercept with species independent
#' re.1 <- list(1, sp = dat$sp, covar = diag(nspp))
#' 
#' # random intercept with species showing phylogenetic covariances
#' re.2 <- list(1, sp = dat$sp, covar = Vphy)
#' 
#' # random slope with species independent
#' re.3 <- list(dat$X, sp = dat$sp, covar = diag(nspp))
#' 
#' # random slope with species showing phylogenetic covariances
#' re.4 <- list(dat$X, sp = dat$sp, covar = Vphy)
#' 
#' # random effect for site
#' re.site <- list(1, site = dat$site, covar = diag(nsite))
#'
#' simple <- communityPGLMM(Y ~ X, data = dat, family = "binomial", sp
#' = dat$sp, site = dat$site, random.effects = list(re.site),
#' REML=TRUE, verbose=FALSE)
#'
#' # The rest of these tests are not run to save CRAN server time;
#' # - please take a look at them because they're *very* useful!
#' \dontrun{ 
#' z.binary <- communityPGLMM(Y ~ X, data = dat, family = "binomial",
#' sp = dat$sp, site = dat$site, random.effects = list(re.1, re.2,
#' re.3, re.4), REML = TRUE, verbose = FALSE)
#' 
#' # output results
#' z.binary
#' plot(z.binary)
#'
#' # test statistical significance of the phylogenetic random effect
#' # on species slopes using a likelihood ratio test
#' communityPGLMM.binary.LRT(z.binary, re.number = 4)$Pr
#' 
#' # extract the predicted values of Y
#' communityPGLMM.predicted.values(z.binary, show.plot = TRUE)
#' 
#' # examine the structure of the overall covariance matrix
#' communityPGLMM.matrix.structure(Y ~ X, data = dat, family =
#' "binomial", sp = dat$sp, site = dat$site, random.effects =
#' list(re.1, re.2, re.3, re.4))
#' 
#' # look at the structure of re.1
#' communityPGLMM.matrix.structure(Y ~ X, data = dat, family =
#' "binomial", sp = dat$sp, site = dat$site, random.effects =
#' list(re.1))
#'
#' # compare results to glmer() when the model contains no
#' # phylogenetic covariance among species; the results should be
#' # similar.
#' communityPGLMM(Y ~ X, data = dat, family = "binomial", sp = dat$sp,
#' site = dat$site, random.effects = list(re.1, re.3), REML = FALSE,
#' verbose = FALSE)
#' 
#' # lmer
#' if(require(lme4)){
#' summary(glmer(Y ~ X + (1 | sp) + (0 + X | sp), data=dat, family =
#' "binomial"))
#' 
#' # compare results to lmer() when the model contains no phylogenetic
#' # covariance among species; the results should be similar.
#' communityPGLMM(Y ~ X, data = dat, family = "gaussian", sp = dat$sp,
#' site = dat$site, random.effects = list(re.1, re.3), REML = FALSE,
#' verbose = FALSE)
#' 
#' # lmer
#' summary(lmer(Y ~ X + (1 | sp) + (0 + X | sp), data=dat, REML = FALSE))
#' }
#' }


######Main PGLMM Function
                                                                                                                  
communityPGLMM <- function(formula, data = list(), family = "gaussian", sp = NULL, site = NULL, 
	random.effects = list(), REML = TRUE, s2.init = NULL, B.init = NULL, reltol = 10^-6, maxit = 500, 
	tol.pql = 10^-6, maxit.pql = 200, verbose = FALSE) {

    #Error checking
    if(!is.factor(sp))
        stop("'sp' must be a factor")
    if(!is.factor(site))
        stop("'site' must be a factor")

	if (family == "gaussian") 
		z <- communityPGLMM.gaussian(formula = formula, data = data, sp = sp, site = site, random.effects = random.effects, 
			REML = REML, s2.init = s2.init, B.init = B.init, reltol = reltol, maxit = maxit, 
			verbose = verbose)
	if (family == "binomial") {

		if (is.null(s2.init)) 
			s2.init <- 0.25
		z <- communityPGLMM.binary(formula = formula, data = data, sp = sp, site = site, random.effects = random.effects, 
			REML = REML, s2.init = s2.init, B.init = B.init, reltol = reltol, maxit = maxit, 
			tol.pql = tol.pql, maxit.pql = maxit.pql, verbose = verbose)
	}
	if (!is.element(family, c("gaussian", "binomial"))) 
		cat("\nSorry, but only binomial (binary) and gaussian options exist at this time")
	return(z)
}

######################################################
######################################################
# communityPLMM.gaussian
######################################################
######################################################
#' @rdname pglmm
#' @importClassesFrom Matrix RsparseMatrix dsCMatrix
#' @importMethodsFrom Matrix t solve %*% determinant diag
#' @importFrom stats pchisq model.frame model.matrix model.response lm var optim pnorm
#' @export
communityPGLMM.gaussian <- function(formula, data = list(), family = "gaussian", sp = NULL, 
	site = NULL, random.effects = list(), REML = TRUE, s2.init = NULL, B.init = NULL, reltol = 10^-8, 
	maxit = 500, verbose = FALSE) {

	# Begin pglmm.LL
	plmm.LL <- function(par, X, Y, Zt, St, nested = NULL, REML, verbose) {
		n <- dim(X)[1]
		p <- dim(X)[2]

		if (!is.null(St)) {
			q.nonNested <- dim(St)[1]
			sr <- Re(par[1:q.nonNested])
			iC <- sr[1] * St[1, ]
			if (length(sr) > 1) 
				for (i in 2:q.nonNested) {
					iC <- iC + sr[i] * St[i, ]
				}
			iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
		} else {
			q.nonNested <- 0
			sr <- NULL
		}
		if (is.null(nested[[1]])) {
			q.Nested <- 0
		} else {
			q.Nested <- length(nested)
		}

		if (q.Nested == 0) {
			sn <- NULL
		} else {
			sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
		}

		if (q.Nested == 0) {
			iA <- as(diag(n), "dsCMatrix")
			Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
			Ut.iA.U <- Ut %*% U
			# Woodbury identity
			iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
		} else {
			A <- as(diag(n), "dsCMatrix")
			for (j in 1:q.Nested) {
				A <- A + sn[j]^2 * nested[[j]]
			}
			iA <- solve(A)
			if (q.nonNested > 0) {
				Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
				Ut.iA.U <- Ut %*% iA %*% U
				iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
			} else {
				iV <- iA
			}
		}

		denom <- t(X) %*% iV %*% X
		num <- t(X) %*% iV %*% Y
		B <- solve(denom, num)
		B <- as.matrix(B)
		H <- Y - X %*% B

		if (q.Nested == 0) {
			# Sylvester identity
			logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
			if (is.infinite(logdetV)) 
				logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
		} else {
			logdetV <- -determinant(iV)$modulus[1]
			if (is.infinite(logdetV)) 
				logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
			if (is.infinite(logdetV)) 
				return(10^10)
		}


		if (REML == TRUE) {
			# concentrated REML likelihood function
			s2.conc <- t(H) %*% iV %*% H/(n - p)
			LL <- 0.5 * ((n - p) * log(s2.conc) + logdetV + (n - p) + determinant(t(X) %*% iV %*% 
				X)$modulus[1])
		} else {
			# concentrated ML likelihood function
			s2.conc <- t(H) %*% iV %*% H/n
			LL <- 0.5 * (n * log(s2.conc) + logdetV + n)
		}

		if (verbose == T) 
			show(c(as.numeric(LL), par))
		return(as.numeric(LL))
	}
	# End plmm.LL
	
	# Main program
	if (is.null(sp) | is.null(site)) 
		stop("Categorical variables for 'sp' and 'site' must be specified")
	nspp <- nlevels(sp)
	nsite <- nlevels(site)

	# order data first by site, second by species
	# sp.order <- order(sp)
# data <- data[sp.order, ]
# sp <- sp[sp.order]
# site <- site[sp.order]

	# site.order <- order(site)
	# data <- data[site.order, ]
# sp <- sp[site.order]
# site <- site[site.order]

	mf <- model.frame(formula = formula, data = data)
	X <- model.matrix(attr(mf, "terms"), data = mf)
	Y <- model.response(mf)

	re <- random.effects
	q <- length(re)

	Ztt <- list(NULL)
	St.lengths <- array(0, q)
	nested <- list(NULL)
	ii <- 0
	jj <- 0

	for (i in 1:q) {
		re.i <- re[[i]]
		# non-nested terms
		if (length(re.i) == 3) {
			counter <- 0
			Z.i <- matrix(0, nrow=nspp * nsite, ncol = nlevels(re.i[[2]]))
			for(i.levels in levels(re.i[[2]])) {
				counter <- counter + 1
				Z.i[,counter] <- re.i[[1]] * as.numeric(i.levels == re.i[[2]])
			}
			Zt.i <- chol(re.i[[3]]) %*% t(Z.i)
			ii <- ii + 1
			Ztt[[ii]] <- Zt.i
			St.lengths[ii] <- nlevels(re.i[[2]])
		}

		# nested terms
		if (length(re.i) == 4) {
			if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
				if (length(re.i[[1]]) > 1) 
					stop("Nested terms can only be for intercepts")
				nestedsp.j <- re.i[[3]]
				nestedsite.j <- diag(nsite)
				nested.j <- as(kronecker(nestedsite.j, nestedsp.j), "dgCMatrix")
			}
			if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
				if (length(re.i[[1]]) > 1) 
					stop("Nested terms can only be for intercepts")
				nestedsp.j <- diag(nspp)
				nestedsite.j <- re.i[[3]]
				nested.j <- as(kronecker(nestedsite.j, nestedsp.j), "dgCMatrix")
			}
			jj <- jj + 1
			nested[[jj]] <- nested.j
		}
	}
	q.nonNested <- ii
	q.Nested <- jj

	if (q.nonNested > 0) {
		St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
		Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * nsite)
		count <- 1
		for (i in 1:q.nonNested) {
			St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
			Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
			count <- count + St.lengths[i]
		}
		Zt <- as(Zt, "dgTMatrix")
		St <- as(St, "dgTMatrix")
	} else {
		Zt <- NULL
		St <- NULL
	}

	p <- ncol(X)
	n <- nrow(X)

	# Compute initial estimates assuming no phylogeny if not provided
	if (!is.null(B.init) & length(B.init) != p) {
		warning("B.init not correct length, so computed B.init using glm()")
	}
	if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & !is.null(s2.init)) {
		B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, ncol = p))
	}
	if (!is.null(B.init) & is.null(s2.init)) {
		s2.init <- var(lm(formula = formula, data = data)$residuals)/q
	}
	if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & is.null(s2.init)) {
		B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, ncol = p))
		s2.init <- var(lm(formula = formula, data = data)$residuals)/q
	}
	B <- B.init
	s <- as.vector(array(s2.init^0.5, dim = c(1, q)))

	if (q > 1) {
		opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, St = St, nested = nested, 
			REML = REML, verbose = verbose, method = "Nelder-Mead", control = list(maxit = maxit, 
				reltol = reltol))
	} else {
		opt <- optim(fn = plmm.LL, par = s, X = X, Y = Y, Zt = Zt, St = St, nested = nested, 
			REML = REML, verbose = verbose, method = "L-BFGS-B", control = list(maxit = maxit))

	}
	# Extract parameters
	par <- abs(Re(opt$par))
	LL <- opt$value
	if (!is.null(St)) {
		q.nonNested <- dim(St)[1]
		sr <- Re(par[1:q.nonNested])
		iC <- sr[1] * St[1, ]
		if (length(sr) > 1) 
			for (i in 2:q.nonNested) {
				iC <- iC + sr[i] * St[i, ]
			}
		iC <- as(diag(iC), "dsCMatrix")
		Ut <- iC %*% Zt
		U <- t(Ut)
	} else {
		q.nonNested <- 0
		sr <- NULL
	}
	if (is.null(nested[[1]])) {
		q.Nested <- 0
	} else {
		q.Nested <- length(nested)
	}
	if (q.Nested == 0) {
		sn <- NULL
	} else {
		sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
	}
	if (q.Nested == 0) {
		iA <- as(diag(n), "dsCMatrix")
		Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
		Ut.iA.U <- Ut %*% U
		# Woodbury identity
		iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
	} else {
		A <- as(diag(n), "dsCMatrix")
		for (j in 1:q.Nested) {
			A <- A + sn[j]^2 * nested[[j]]
		}
		iA <- solve(A)
		if (q.nonNested > 0) {
			Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
			Ut.iA.U <- Ut %*% iA %*% U
			iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
		} else {
			iV <- iA
		}
	}

	denom <- t(X) %*% iV %*% X
	num <- t(X) %*% iV %*% Y
	B <- solve(denom, num)
	B <- as.matrix(B)
	H <- Y - X %*% B

	if (q.Nested == 0) {
		# Sylvester identity
		logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
		if (is.infinite(logdetV)) 
			logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
	} else {
		logdetV <- -determinant(iV)$modulus[1]
		if (is.infinite(logdetV)) 
			logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
		if (is.infinite(logdetV)) 
			return(10^10)
	}

	if (REML == TRUE) {
		s2resid <- as.numeric(t(H) %*% iV %*% H/(n - p))
	} else {
		s2resid <- as.numeric(t(H) %*% iV %*% H/n)
	}

	s2r <- s2resid * sr^2
	s2n <- s2resid * sn^2
	ss <- c(sr, sn, s2resid^0.5)

	iV <- iV/s2resid

	B.cov <- solve(t(X) %*% iV %*% X)
	B.se <- as.matrix(diag(B.cov))^0.5
	B.zscore <- B/B.se
	B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)

	if (REML == TRUE) {
		logLik <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% X)$modulus[1] - 
			LL
	} else {
		logLik <- -0.5 * n * log(2 * pi) - LL
	}
	k <- p + q + 1
	AIC <- -2 * logLik + 2 * k
	BIC <- -2 * logLik + k * (log(n) - log(pi))

        #NOTE: the nested/non-nested returns are 'swapped' here, which looks dodgy but is better than them being the wrong way round!
	results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
		B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, ss = ss, 
		s2n = s2n, s2r = s2r, s2resid = s2resid, logLik = logLik, AIC = AIC, BIC = BIC, REML = REML, 
		s2.init = s2.init, B.init = B.init, Y = Y, X = X, H = H, iV = iV, mu = NULL, nested = nested, 
		sp = sp, site = site, Zt = Zt, St = St, convcode = opt$convergence, niter = opt$counts)

	class(results) <- "communityPGLMM"
	results
}


######################################################
######################################################
# communityPGLMM.binary
######################################################
######################################################
#' \code{communityPGLMM.binary} calculates the statistical
#' significance of the random effects in the generalized linear mixed
#' model from the marginal profile likelihood.
#' @rdname pglmm
#' @importClassesFrom Matrix dsCMatrix RsparseMatrix
#' @importMethodsFrom Matrix t solve %*% determinant diag
#' @importFrom methods as show
#' @importFrom stats model.frame model.matrix model.response glm binomial optim pchisq
#' @export
communityPGLMM.binary <- function(formula, data = list(), family = "binomial", sp = NULL, site = NULL, 
	random.effects = list(), REML = TRUE, s2.init = 0.25, B.init = NULL, reltol = 10^-5, maxit = 40, 
	tol.pql = 10^-6, maxit.pql = 200, verbose = FALSE) {
	plmm.binary.V <- function(par, Zt, St, mu, nested) {

		if (!is.null(St)) {
			q.nonNested <- dim(St)[1]
			sr <- Re(par[1:q.nonNested])
			iC <- sr[1] * St[1, ]
			if (length(sr) > 1) 
				for (i in 2:q.nonNested) {
					iC <- iC + sr[i] * St[i, ]
				}
			iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
		} else {
			q.nonNested <- 0
			sr <- NULL
		}

		if (q.Nested == 0) {
			sn <- NULL
		} else {
			sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
		}

		iW <- diag(as.vector((mu * (1 - mu))^-1))
		if (q.Nested == 0) {
			A <- iW
		} else {
			A <- iW
			for (j in 1:q.Nested) {
				A <- A + sn[j]^2 * nested[[j]]
			}
		}
		if (q.nonNested > 0) {
			V <- A + U %*% Ut
		} else {
			V <- A
		}
		return(V)
	}
	# End plmm.binary.V
	
	plmm.binary.iV <- function(par, Zt, St, mu, nested) {

		if (!is.null(St)) {
			q.nonNested <- dim(St)[1]
			sr <- Re(par[1:q.nonNested])
			iC <- sr[1] * St[1, ]
			if (length(sr) > 1) 
				for (i in 2:q.nonNested) {
					iC <- iC + sr[i] * St[i, ]
				}
			iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
		} else {
			q.nonNested <- 0
			sr <- NULL
		}

		if (is.null(nested[[1]])) {
			q.Nested <- 0
		} else {
			q.Nested <- length(nested)
		}

		if (q.Nested == 0) {
			sn <- NULL
		} else {
			sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
		}

		if (q.Nested == 0) {
			iA <- diag(as.vector((mu * (1 - mu))))
			Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
			Ut.iA.U <- Ut %*% iA %*% U
			# Woodbury identity
			iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
		} else {
			A <- as(diag(as.vector((mu * (1 - mu))^-1)), "dgCMatrix")
			for (j in 1:q.Nested) {
				A <- A + sn[j]^2 * nested[[j]]
			}
			iA <- solve(A)

			if (q.nonNested > 0) {
				Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
				Ut.iA.U <- Ut %*% iA %*% U
				iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
			} else {
				iV <- iA
			}
		}
		return(iV)
	}
	# End plmm.binary.iV
	
	plmm.binary.logdetV <- function(par, Zt, St, mu, nested) {

		if (!is.null(St)) {
			q.nonNested <- dim(St)[1]
			sr <- Re(par[1:q.nonNested])
			iC <- sr[1] * St[1, ]
			if (length(sr) > 1) 
				for (i in 2:q.nonNested) {
					iC <- iC + sr[i] * St[i, ]
				}
			iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
		} else {
			q.nonNested <- 0
			sr <- NULL
		}

		if (is.null(nested[[1]])) {
			q.Nested <- 0
		} else {
			q.Nested <- length(nested)
		}

		if (q.Nested == 0) {
			sn <- NULL
		} else {
			sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
		}

		if (q.Nested == 0) {
			iA <- diag(as.vector((mu * (1 - mu))))
			Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
			Ut.iA.U <- Ut %*% iA %*% U
			logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
			if (is.infinite(logdetV)) 
				logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
		} else {
			A <- as(diag(as.vector((mu * (1 - mu))^-1)), "dsCMatrix")
			for (j in 1:q.Nested) {
				A <- A + sn[j]^2 * nested[[j]]
			}
			iA <- solve(A)
			if (q.nonNested > 0) {
				Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
				Ut.iA.U <- Ut %*% iA %*% U
				iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
			} else {
				iV <- iA
			}
			logdetV <- -determinant(iV)$modulus[1]
			if (is.infinite(logdetV)) 
				logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
			if (is.infinite(logdetV)) 
				return(list(iV = NULL, logdetV = NULL))
		}

		return(logdetV)
	}
	# End plmm.binary.logdetV
	
	# Begin pglmm.binary.LL
	plmm.binary.LL <- function(par, H, X, Zt, St, mu, nested, REML = TRUE, verbose = FALSE) {
		par <- abs(par)
		n <- dim(H)[1]
		p <- dim(H)[2]

		iV <- plmm.binary.iV(par = par, Zt = Zt, St = St, mu = mu, nested = nested)
		logdetV <- plmm.binary.logdetV(par = par, Zt = Zt, St = St, mu = mu, nested = nested)
		if (REML == TRUE) {
			# REML likelihood function
			LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + determinant(t(X) %*% iV %*% X)$modulus[1])
		} else {
			# ML likelihood function
			LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
		}
		if (verbose == T) 
			show(c(as.numeric(LL), par))

		return(as.numeric(LL))
	}
	# End plmm.binary.LL
	
	############################################################
	# Begin main program
if (is.null(sp) | is.null(site)) 
		stop("Categorical variables for 'sp' and 'site' must be specified")
	nspp <- nlevels(sp)
	nsite <- nlevels(site)

	# order data first by site, second by species
	# sp.order <- order(sp)
# data <- data[sp.order, ]
# sp <- sp[sp.order]
# site <- site[sp.order]

	# site.order <- order(site)
	# data <- data[site.order, ]
# sp <- sp[site.order]
# site <- site[site.order]

	mf <- model.frame(formula = formula, data = data)
	X <- model.matrix(attr(mf, "terms"), data = mf)
	Y <- model.response(mf)

	re <- random.effects
	q <- length(re)

	Ztt <- list(NULL)
	St.lengths <- array(0, q)
	nested <- list(NULL)
	ii <- 0
	jj <- 0
	for (i in 1:q) {
		re.i <- re[[i]]
		# non-nested terms
		if (length(re.i) == 3) {
			counter <- 0
			Z.i <- matrix(0, nrow=nspp * nsite, ncol = nlevels(re.i[[2]]))
			for(i.levels in levels(re.i[[2]])) {
				counter <- counter + 1
				Z.i[,counter] <- re.i[[1]] * as.numeric(i.levels == re.i[[2]])
			}
			Zt.i <- chol(re.i[[3]]) %*% t(Z.i)
			ii <- ii + 1
			Ztt[[ii]] <- Zt.i
			St.lengths[ii] <- nlevels(re.i[[2]])
		}

		# nested terms
		if (length(re.i) == 4) {
			if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
				if (length(re.i[[1]]) > 1) 
					stop("Nested terms can only be for intercepts")
				nestedsp.j <- re.i[[3]]
				nestedsite.j <- diag(nsite)
				nested.j <- as(kronecker(nestedsite.j, nestedsp.j), "dgCMatrix")
			}
			if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
				if (length(re.i[[1]]) > 1) 
					stop("Nested terms can only be for intercepts")
				nestedsp.j <- diag(nspp)
				nestedsite.j <- re.i[[3]]
				nested.j <- as(kronecker(nestedsite.j, nestedsp.j), "dgCMatrix")
			}
			jj <- jj + 1
			nested[[jj]] <- nested.j
		}
	}
	q.nonNested <- ii
	q.Nested <- jj

	if (q.nonNested > 0) {
		St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
		Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * nsite)
		count <- 1
		for (i in 1:q.nonNested) {
			St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
			Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
			count <- count + St.lengths[i]
		}
		Zt <- as(Zt, "dgTMatrix")
		St <- as(St, "dgTMatrix")
	} else {
		Zt <- NULL
		St <- NULL
	}

	p <- ncol(X)
	n <- nrow(X)

	# Compute initial estimates assuming no phylogeny if not provided
	if (!is.null(B.init) & length(B.init) != p) {
		warning("B.init not correct length, so computed B.init using glm()")
	}
	if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p))) {
		B.init <- t(matrix(glm(formula = formula, data = data, family = binomial)$coefficients, 
			ncol = p))
	} else {
		B.init <- matrix(B.init, ncol = 1)
	}
	B <- B.init
	ss <- as.vector(array(s2.init^0.5, dim = c(1, q)))

	b <- matrix(0, nrow = n)
	beta <- rbind(B, b)
	mu <- exp(X %*% B)/(1 + exp(X %*% B))
	XX <- cbind(X, diag(1, nrow = n, ncol = n))

	est.ss <- ss
	est.B <- B
	oldest.ss <- 10^6
	oldest.B <- matrix(10^6, nrow = length(est.B))

	iteration <- 0
	exitflag <- 0
	rcondflag <- 0
	while (((t(est.ss - oldest.ss) %*% (est.ss - oldest.ss) > tol.pql^2) | (t(est.B - oldest.B) %*% 
		(est.B - oldest.B) > tol.pql^2)) & (iteration <= maxit.pql)) {

		iteration <- iteration + 1
		oldest.ss <- est.ss
		oldest.B <- est.B

		est.B.m <- B
		oldest.B.m <- matrix(10^6, nrow = length(est.B))
		iteration.m <- 0

		# mean component
		while ((t(est.B.m - oldest.B.m) %*% (est.B.m - oldest.B.m) > tol.pql^2) & (iteration.m <= 
			maxit.pql)) {
			iteration.m <- iteration.m + 1

			oldest.B.m <- est.B.m

			iV <- plmm.binary.iV(par = ss, Zt = Zt, St = St, mu = mu, nested = nested)

			Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
			denom <- t(X) %*% iV %*% X
			num <- t(X) %*% iV %*% Z
			B <- solve(denom, num)
			B <- as.matrix(B)

			V <- plmm.binary.V(par = ss, Zt = Zt, St = St, mu = mu, nested = nested)
			iW <- diag(as.vector((mu * (1 - mu))^-1))
			C <- V - iW
			b <- C %*% iV %*% (Z - X %*% B)
			beta <- rbind(B, matrix(b))
			mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))

			est.B.m <- B
			if (verbose == TRUE) 
				show(c(iteration, B))
		  
      if(any(is.nan(B))){
				stop("Estimation of B failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help.")
			}
    }

		# variance component
		Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
		H <- Z - X %*% B
		if (q > 1) {
			opt <- optim(fn = plmm.binary.LL, par = ss, H = H, X = X, Zt = Zt, St = St, mu = mu, 
				nested = nested, REML = REML, verbose = verbose, method = "Nelder-Mead", control = list(maxit = maxit, 
					reltol = reltol))
		} else {
			opt <- optim(fn = plmm.binary.LL, par = ss, H = H, X = X, Zt = Zt, St = St, mu = mu, 
				nested = nested, REML = REML, verbose = verbose, method = "L-BFGS-B", control = list(maxit = maxit))
		}
		ss <- abs(opt$par)
		LL <- opt$value

		est.ss <- ss
		est.B <- B
	}

	# Extract parameters
	if (q.nonNested > 0) {
		sr <- ss[1:q.nonNested]
	} else {
		sr <- NULL
	}
	if (q.Nested > 0) {
		sn <- ss[(q.nonNested + 1):(q.nonNested + q.Nested)]
	} else {
		sn <- NULL
	}

	s2r <- sr^2
	s2n <- sn^2

	B.cov <- solve(t(X) %*% iV %*% X)
	B.se <- as.matrix(diag(B.cov))^0.5
	B.zscore <- B/B.se
	B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)

	results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
		B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, ss = ss, 
		s2n = s2n, s2r = s2r, s2resid = NULL, logLik = NULL, AIC = NULL, BIC = NULL, REML = REML, 
		s2.init = s2.init, B.init = B.init, Y = Y, X = X, H = H, iV = iV, mu = mu, nested = nested, sp = sp, 
		site = site, Zt = Zt, St = St, convcode = opt$convergence, niter = opt$counts)
	class(results) <- "communityPGLMM"
	results
	return(results)
}

######################################################
######################################################
# communityPGLMM.binary.LRT
######################################################
######################################################
#' \code{communityPGLMM.binary.LRT} tests statistical significance of
#' the phylogenetic random effect on species slopes using a likelihood
#' ratio test
#' @param re.number which \code{random.effect} in \code{x} to be
#' tested
#' @param x \code{communityPGLMM} object
#' @rdname pglmm
#' @importClassesFrom Matrix dsCMatrix
#' @importMethodsFrom Matrix t solve %*% determinant diag
#' @importFrom methods as show
#' @export
communityPGLMM.binary.LRT <- function(x, re.number = 0, ...) {

	plmm.binary.iV <- function(par, Zt, St, mu, nested = NULL) {

		if (!is.null(St)) {
			q.nonNested <- dim(St)[1]
			sr <- Re(par[1:q.nonNested])
			iC <- sr[1] * St[1, ]
			if (length(sr) > 1) 
				for (i in 2:q.nonNested) {
					iC <- iC + sr[i] * St[i, ]
				}
			iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
		} else {
			q.nonNested <- 0
			sr <- NULL
		}
		if (is.null(nested[[1]])) {
			q.Nested <- 0
		} else {
			q.Nested <- length(nested)
		}

		if (q.Nested == 0) {
			sn <- NULL
		} else {
			sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
		}

		if (q.Nested == 0) {
			iA <- diag(as.vector((mu * (1 - mu))))
			Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
			Ut.iA.U <- Ut %*% iA %*% U
			# Woodbury identity
			iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
		} else {
			A <- as(diag(as.vector((mu * (1 - mu))^-1)), "dgCMatrix")
			for (j in 1:q.Nested) {
				A <- A + sn[j]^2 * nested[[j]]
			}
			iA <- solve(A)
			if (q.nonNested > 0) {
				Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
				Ut.iA.U <- Ut %*% iA %*% U
				iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
			} else {
				iV <- iA
			}
		}
		return(iV)
	}
	# End plmm.binary.iV
	
	plmm.binary.logdetV <- function(par, Zt, St, mu, nested = NULL) {

		if (!is.null(St)) {
			q.nonNested <- dim(St)[1]
			sr <- Re(par[1:q.nonNested])
			iC <- sr[1] * St[1, ]
			if (length(sr) > 1) 
				for (i in 2:q.nonNested) {
					iC <- iC + sr[i] * St[i, ]
				}
			iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
		} else {
			q.nonNested <- 0
			sr <- NULL
		}
		if (is.null(nested[[1]])) {
			q.Nested <- 0
		} else {
			q.Nested <- length(nested)
		}

		if (q.Nested == 0) {
			sn <- NULL
		} else {
			sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
		}

		if (q.Nested == 0) {
			iA <- diag(as.vector((mu * (1 - mu))))
			Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
			Ut.iA.U <- Ut %*% iA %*% U
			logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
			if (is.infinite(logdetV)) 
				logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
		} else {
			A <- as(diag(as.vector((mu * (1 - mu))^-1)), "dgCMatrix")
			for (j in 1:q.Nested) {
				A <- A + sn[j]^2 * nested[[j]]
			}
			iA <- solve(A)
			if (q.nonNested > 0) {
				Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
				Ut.iA.U <- Ut %*% iA %*% U
				iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
			} else {
				iV <- iA
			}
			logdetV <- -determinant(iV)$modulus[1]
			if (is.infinite(logdetV)) 
				logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
			if (is.infinite(logdetV)) 
				return(list(iV = NULL, logdetV = NULL))
		}

		return(logdetV)
	}
	# End plmm.binary.logdetV
	
	# Begin pglmm.binary.LL
	plmm.binary.LL <- function(par, H, X, Zt, St, mu, nested = NULL, REML = TRUE, verbose = FALSE) {
		par <- abs(par)
		n <- dim(H)[1]
		p <- dim(H)[2]

		iV <- plmm.binary.iV(par = par, Zt = Zt, St = St, mu = mu, nested = nested)
		logdetV <- plmm.binary.logdetV(par = par, Zt = Zt, St = St, mu = mu, nested = nested)
		if (REML == TRUE) {
			# REML likelihood function
			LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + determinant(t(X) %*% iV %*% X)$modulus[1])
		} else {
			# ML likelihood function
			s2.conc <- t(H) %*% iV %*% H/(n - p)
			LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
		}
		if (verbose == T) 
			show(c(as.numeric(LL), par))

		return(as.numeric(LL))
	}
	# End plmm.binary.LL
	
	############################################################
	# Begin main program

	n <- dim(x$X)[1]
	p <- dim(x$X)[2]
	df <- length(re.number)

	
	q.nonNested <- length(x$ss) - length(x$nested)
	ss0 <- x$ss
	ss0[re.number] <- 0

	LL <- plmm.binary.LL(par = x$ss, H = x$H, X = x$X, Zt = x$Zt, St = x$St, mu = x$mu, nested = x$nested, 
		REML = x$REML)
	if (x$REML == TRUE) {
		logLik <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - 
			LL
	} else {
		logLik <- -0.5 * n * log(2 * pi) - LL
	}

	LL0 <- plmm.binary.LL(par = ss0, H = x$H, X = x$X, Zt = x$Zt, St = x$St, mu = x$mu, nested = x$nested, 
		REML = x$REML)
	if (x$REML == TRUE) {
		logLik0 <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - 
			LL0
	} else {
		logLik0 <- -0.5 * n * log(2 * pi) - LL0
	}

	P.H0.s2 <- pchisq(2 * (logLik - logLik0), df = df, lower.tail = F)/2
	return(list(LR = logLik - logLik0, df = df, Pr = P.H0.s2))
}

######################################################
######################################################
# communityPGLMM.matrix.structure
######################################################
######################################################
#' \code{communityPGLMM.matrix.structure} produces the entire
#' covariance matrix structure (V) when you specify random effects.
#' @importClassesFrom Matrix dsCMatrix RsparseMatrix
#' @importClassesFrom Matrix dsCMatrix RsparseMatrix
#' @importFrom methods as
#' @importFrom stats model.frame model.matrix model.response
#' @param ss which of the \code{random.effects} to produce
#' @rdname pglmm
#' @export
communityPGLMM.matrix.structure <- function(formula, data = list(), family = "binomial", sp = NULL, 
	site = NULL, random.effects = list(), ss = 1) {
	plmm.binary.V.test <- function(par, Zt, St, X, nestedsp = NULL, nestedsite = NULL) {
		n <- nrow(X)

		if (!is.null(St)) {
			q.nonNested <- dim(St)[1]
			sr <- Re(par[1:q.nonNested])
			iC <- sr[1] * St[1, ]
			if (length(sr) > 1) 
				for (i in 2:q.nonNested) {
					iC <- iC + sr[i] * St[i, ]
				}
			iC <- as(diag(iC), "dsCMatrix")
			Ut <- iC %*% Zt
			U <- t(Ut)
		} else {
			q.nonNested <- 0
			sr <- NULL
		}
		if (is.null(nestedsp[[1]])) {
			q.Nested <- 0
		} else {
			q.Nested <- length(nestedsp)
		}

		if (q.Nested == 0) {
			sn <- NULL
		} else {
			sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
		}

		iW <- 0 * diag(n)
		if (q.Nested == 0) {
			A <- iW
		} else {
			A <- iW
			for (j in 1:q.Nested) {
				A <- A + sn[j]^2 * kronecker(nestedsite[[j]], nestedsp[[j]])
			}
		}
		if (q.nonNested > 0) {
			V <- A + U %*% Ut
		} else {
			V <- A
		}
		return(V)
	}
	# End plmm.binary.V
	
	############################################################
	# Begin main program
if (is.null(sp) | is.null(site)) 
		stop("Categorical variables for 'sp' and 'site' must be specified")
	nspp <- nlevels(sp)
	nsite <- nlevels(site)

	# order data first by site, second by species
	sp.order <- order(sp)
	data <- data[sp.order, ]
	sp <- sp[sp.order]
	site <- site[sp.order]

	site.order <- order(site)
	data <- data[site.order, ]
	sp <- sp[site.order]
	site <- site[site.order]

	mf <- model.frame(formula = formula, data = data)
	X <- model.matrix(attr(mf, "terms"), data = mf)
	Y <- model.response(mf)

	re <- random.effects
	q <- length(re)

	Ztt <- list(NULL)
	St.lengths <- array(0, q)
	nestedsp <- list(NULL)
	nestedsite <- list(NULL)
	ii <- 0
	jj <- 0

	for (i in 1:q) {
		re.i <- re[[i]]
		# non-nested terms
		if (length(re.i) == 3) {
			counter <- 0
			Z.i <- matrix(0, nrow=nspp * nsite, ncol = nlevels(re.i[[2]]))
			for(i.levels in levels(re.i[[2]])) {
				counter <- counter + 1
				Z.i[,counter] <- re.i[[1]] * as.numeric(i.levels == re.i[[2]])
			}
			Zt.i <- chol(re.i[[3]]) %*% t(Z.i)
			ii <- ii + 1
			Ztt[[ii]] <- Zt.i
			St.lengths[ii] <- nlevels(re.i[[2]])
		}

		# nested terms
		if (length(re.i) == 4) {
			if (setequal(levels(re.i[[2]]), levels(sp)) && all(re.i[[2]] == sp)) {
				if (length(re.i[[1]]) > 1) 
					stop("Nested terms can only be for intercepts")
				nestedsp.j <- re.i[[3]]
				nestedsite.j <- diag(nsite)
			}
			if (setequal(levels(re.i[[2]]), levels(site)) && all(re.i[[2]] == site)) {
				if (length(re.i[[1]]) > 1) 
					stop("Nested terms can only be for intercepts")
				nestedsp.j <- diag(nspp)
				nestedsite.j <- re.i[[3]]
			}
			jj <- jj + 1
			nestedsp[[jj]] <- nestedsp.j
			nestedsite[[jj]] <- nestedsite.j
		}
	}
	q.nonNested <- ii
	q.Nested <- jj

	if (q.nonNested > 0) {
		St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
		Zt <- matrix(0, nrow = sum(St.lengths), ncol = nspp * nsite)
		count <- 1
		for (i in 1:q.nonNested) {
			St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
			Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
			count <- count + St.lengths[i]
		}
		Zt <- as(Zt, "dgTMatrix")
		St <- as(St, "dgTMatrix")
	} else {
		Zt <- NULL
		St <- NULL
	}

	V <- plmm.binary.V.test(par = array(ss, c(1, q)), Zt = Zt, St = St, X = X, nestedsp = nestedsp, 
		nestedsite = nestedsite)
	return(V)
}



######################################################
######################################################
# summary.communityPGLMM
######################################################
######################################################
#' @rdname pglmm
#' @method summary communityPGLMM
#' @param digits minimal number of significant digits for printing, as
#' in \code{\link{print.default}}
#' @param object communityPGLMM object to be summarised
#' @importFrom stats printCoefmat
#' @export
summary.communityPGLMM <- function(object, digits = max(3, getOption("digits") - 3), ...) {
	if (object$family == "gaussian") {
		if (object$REML == TRUE) 
			cat("Linear mixed model fit by restricted maximum likelihood")
		else cat("Linear mixed model fit by maximum likelihood")
	}
	if (object$family == "binomial") {
		if (object$REML == TRUE) 
			cat("Generalized linear mixed model for binary data fit by restricted maximum likelihood")
		else cat("Generalized linear mixed model for binary data fit by maximum likelihood")
	}

	cat("\n\nCall:")
	print(object$formula)
	cat("\n")

	if (object$family == "gaussian") {

		logLik = object$logLik
		AIC = object$AIC
		BIC = object$BIC

		names(logLik) = "logLik"
		names(AIC) = "AIC"
		names(BIC) = "BIC"
		print(c(logLik, AIC, BIC), digits = digits)
	}
	cat("\nRandom effects:\n")
	w <- data.frame(Variance = matrix(c(object$s2r, object$s2n, object$s2resid), ncol = 1), Std.Dev = matrix(c(object$s2r^0.5, 
		object$s2n^0.5, object$s2resid^0.5), ncol = 1))

	re.names <- NULL
	if (length(object$s2r) > 0) 
		for (i in 1:length(object$s2r)) re.names <- c(re.names, paste("non-nested ", i, sep = ""))

	if (length(object$s2n) > 0) 
		for (i in 1:length(object$s2n)) re.names <- c(re.names, paste("nested ", i, sep = ""))

	if (object$family == "gaussian") 
		re.names <- c(re.names, "residual")

	row.names(w) <- re.names
	print(w, digits = digits)

	cat("\nFixed effects:\n")
	coef <- data.frame(Value = object$B, Std.Error = object$B.se, Zscore = object$B.zscore, Pvalue = object$B.pvalue)
	printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
	cat("\n")
}

print.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
	summary.communityPGLMM(x, digits = digits)
}

######################################################
######################################################
# plot.communityPGLMM
######################################################
######################################################
#' @rdname pglmm
#' @method plot communityPGLMM
#' @importFrom graphics par image
#' @importFrom stats reshape
#' @export
plot.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    #Wrangle data
    W <- data.frame(Y = x$Y, sp = x$sp, site = x$site)
    Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
    Y <- Y[, 2:dim(Y)[2]]
    
    #Wrangle for image
    Y <- t(Y)
    image(x=seq(1,nrow(Y)), y=seq(1,ncol(Y)), z=Y, ylab = "species", xlab = "sites", main = "Observed values", col=c("black","white"))
}


######################################################
######################################################
# communityPGLMM.predicted.values
######################################################
######################################################
#' \code{communityPGLMM.predicted.values} calculates the predicted
#' values of Y; for the generalized linear mixed model (family =
#' "binomial"), these values are in the logit-1 transformed space.
#' @rdname pglmm
#' @param show.plot if \code{TRUE} (default), display plot
#' @importFrom graphics par
#' @importFrom stats reshape
#' @export
communityPGLMM.predicted.values <- function(x, show.plot = TRUE, ...) {

	if (x$family == "gaussian") {
		V <- solve(x$iV)
		h <- matrix(0, nrow = length(x$Y), ncol = 1)
		for (i in 1:length(x$Y)) {
			h[i] <- as.numeric(V[i, -i] %*% solve(V[-i, -i]) %*% matrix(x$H[-i]))
		}
		predicted.values <- h
	}

	if (x$family == "binomial") {
		h <- x$H + x$X %*% x$B
		predicted.values <- as.numeric(h)
	}

	if (show.plot == TRUE) {
		W <- data.frame(Y = predicted.values, sp = x$sp, site = x$site)
		Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
		Y <- Y[, 2:dim(Y)[2]]
		par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2) - 0.1)
                Y <- t(Y)
                image(x=seq(1,nrow(Y)), y=seq(1,ncol(Y)), z=Y, ylab = "species", xlab = "sites", main = "Predicted values", col=c("black","white"))
	}
	return(predicted.values)
}
