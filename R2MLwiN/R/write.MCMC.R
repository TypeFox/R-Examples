#' Writes MLwiN macros to fit models using Markov chain Monte Carlo (MCMC)
#' methods
#'
#' write.MCMC is an internal function which creates an MLwiN macro file to
#' fit models using MCMC.
#'
#' @param indata A data.frame object containing the data to be modelled.
#' @param dtafile The file name of the dataset to be imported into MLwiN, which
#' is in Stata format (i.e. with extension .dta).
#' @param oldsyntax Specified as \code{FALSE} if new syntax has been used in
#' \code{Formula} object, else specified as \code{TRUE} (to enable
#' backcompatibility).
#' @param resp A character string (vector) of response variable(s).
#' @param levID A character string (vector) of the specified level ID(s). The
#' ID(s) should be sorted in the descending order of levels (e.g.
#' \code{levID = c('level2', 'level1')} where \code{'level2'} is the higher
#' level).
#' @param expl A character string (vector) of explanatory (predictor)
#' variable(s).
#' @param rp A character string (vector) of random part of random variable(s).
#' @param D A character string/vector specifying the type of distribution to be modelled, which
#' can include \code{'Normal'} (the default), \code{'Binomial'}, \code{'Poisson'},
#' \code{'Unordered Multinomial'}, \code{'Ordered Multinomial'},
#' \code{'Multivariate Normal'}, or \code{'Mixed'}. In the case of the latter,
#' \code{'Mixed'} precedes the response types which also need to be listed in
#' \code{D}, e.g. \code{c('Mixed', 'Normal', 'Binomial')}; these need to be
#' be listed in the same order to which they are referred to in the
#' \code{Formula} object (see \code{\link{runMLwiN}}, \code{\link{Formula.translate}},
#' \code{\link{Formula.translate.compat}}. \code{'Mixed'} combinations can only consist of
#' \code{'Normal'} and \code{'Binomial'} for MCMC estimation.
#' @param nonlinear A character vector specifying linearisation method for IGLS
#' starting values for discrete
#' response models (see Chapter 9 of Rasbash et al 2012, and Goldstein 2011).
#' \code{N = 0} specifies marginal quasi-likelihood
#' linearization (MQL), whilst \code{N = 1} specifies penalised quasi-
#' likelihood linearization (PQL); \code{M = 1} specifies first order
#' approximation, whilst \code{M = 2} specifies second order approximation.
#' \code{nonlinear = c(N = 0, M = 1)} by default. First order marginal
#' quasi-likelihood (MQL1) only option for single-level discrete response
#' models.
#' @param categ Specifies categorical variable(s) as a matrix. Each column
#' corresponds to a categorical variable; the first row specifies the name(s)
#' of variable(s); the second row specifies the name(s) of reference group(s),
#' \code{NA}(s) if no reference group; the third row states the number of
#' categories for each variable.
#' @param notation Specifies the model subscript notation to be used in the
#' MLwiN equations window. \code{'class'} means no multiple subscripts, whereas
#' \code{'level'} has multiple subscripts.
#' @param nonfp Removes the fixed part of random variable(s). \code{NA} if no
#' variable is removed.
#' @param clre A matrix used to estimate some, but not all, of the variances
#' and covariances for a set of coefficients at a particular level. Remove from
#' the random part at level <first row> the covariance matrix element(s)
#' defined by the pair(s) of rows <second row> <third row>. Each row
#' corresponds to a removed entry of the covariance matrix.
#' @param Meth Specifies the maximum likelihood estimation method to be used
#' when generating starting values via (R)IGLS.
#' If \code{Meth = 0} estimation method is set to RIGLS. If \code{Meth = 1}
#' estimation method is set to IGLS (the default setting). If \code{Meth} is
#' absent, alternate between IGLS and RIGLS.
#' @param merr A vector which sets-up measurement errors on predictor
#' variables. The first element \code{N} defines the number of variables that
#' have measurement errors. Then, for each variable with measurement error, a
#' pair of inputs is required: the first of which is the explanatory variable
#' name as a character string, and the second of which is the variance of
#' the measurement error for this variable. See \code{demo(MCMCGuide14)} for an
#' example.
#' @param carcentre If CAR model (i.e. if \code{car} is non-\code{NULL}),
#' \code{carcentre = TRUE} mean-centres all random effects at that level.
#' @param maxiter When generating starting values via (R)IGLS, a numeric
#' value specifying the total number of iterations, from
#' the start, before IGLS estimation halts (if \code{startval = NULL}).
#' @param convtol When generating starting values via (R)IGLS, a numeric
#' value specifying the IGLS convergence criterion, as
#' specified in the \code{tol} option within \code{estoptions}, where
#' \code{startval = NULL}) (see \code{\link{runMLwiN}}). If value of
#' \code{convtol} is m, estimation will be deemed to have converged when the
#' relative change in the estimate for any parameter from one iteration to the
#' next is less than 10(-m). Defaults to value of \code{2} for m if not
#' otherwise specified.
#' @param seed An integer specifying the random seed in MLwiN.
#' @param iterations An integer specifying the number of iterations after
#' burn-in.
#' @param burnin An integer specifying length of the burn-in.
#' @param scale An integer specifying the scale factor of proposed variances;
#' this number will be multiplied by the estimated
#' parameter variance (from IGLS/RIGLS) to give the proposal distribution variance.
#' @param thinning An integer specifying the frequency with which successive
#' values in the Markov chain are stored. By default \code{thinning = 1}.
#' @param priorParam A vector specifying the informative priors used, as output
#' from \code{\link{prior2macro}}.
#' @param refresh An integer specifying how frequently the parameter estimates
#' are refreshed on the screen during iterations; only applies if
#' \code{debugmode = TRUE} in \code{estoptions}:
#' see \code{\link[R2MLwiN]{runMLwiN}}.
#' @param fixM Specifies the fixed effect method: \code{1} for Gibbs Sampling,
#' \code{2} for univariate MH Sampling and \code{3} for multivariate MH
#' Sampling.
#' @param residM Specifies the residual method: \code{1} for Gibbs Sampling,
#' \code{2} for univariate MH Sampling and \code{3} for multivariate MH
#' Sampling.
#' @param Lev1VarM Specifies the level 1 variance method: \code{1} for Gibbs
#' Sampling, \code{2} for univariate MH Sampling and \code{3} for multivariate
#' MH Sampling.
#' @param OtherVarM Specifies the variance method for other levels: \code{1}
#' for Gibbs Sampling, \code{2} for univariate MH Sampling and \code{3} for
#' multivariate MH Sampling.
#' @param adaption \code{adaption = 1} indicates adaptation is to be used;
#' \code{0} otherwise.
#' @param priorcode A vector indicating which default priors are to be used
#' for the variance parameters. It defaults to \code{c(gamma = 1)} in which case
#' Gamma priors are used with MLwiN's defaults of Gamma a value (shape) = 0.001
#' and Gamma b value (scale) = 0.001, although alternative values for shape and
#' scale can be specified in subsequent elements of the vector,
#' e.g. \code{c(gamma = 1, shape = 0.5, scale = 0.2)}). Alternatively
#' \code{c(uniform = 1)} specifies Uniform priors on the variance scale. To allow
#' for back-compatibility with deprecated syntax used in versions of
#' \pkg{R2MLwiN} prior to 0.8-2, if \code{priorcode} is instead specified as
#' an integer, then \code{1} indicates that Gamma priors are used, whereas
#' \code{0} indicates that Uniform priors are used. See the section on 'Priors' in the
#' MLwiN help system for more details on the meaning of these priors.
#' @param rate An integer specifying the acceptance rate (as a percentage);
#' this command is ignored if \code{adaption = 0}.
#' @param tol An integer specifying tolerance (as a percentage) for the acceptance rate.
#' @param lclo This command toggles on/off the possible forms of complex level
#' 1 variation when using MCMC. \code{lclo = 0} expresses the level
#' 1 variation as a function of the predictors, whereas \code{lclo = 1} expresses the
#' log of the level 1 precision (1/variance) as a function of the predictors.
#' @param mcmcOptions A list of other MCMC options used. See `Details' below.
#' @param fact A list of objects specified for factor analysis. See `Details'
#' below.
#' @param xc Indicates whether model is cross-classified (\code{TRUE}) or
#' nested (\code{FALSE}). \code{xc = NULL} by default (corresponding to
#' \code{FALSE}), unless either \code{mm} or \code{car} are not null, in
#' which case \code{xc = TRUE}.
#' @param mm Specifies the structure of a multiple membership model.
#' Can be a list of variable names, a list of vectors, or a matrix (e.g. see
#' \code{\link{df2matrix}}). In the case of the former, each
#' element of the list corresponds to a level (classification) of the model,
#' in descending order. If a level is not a multiple membership classification,
#' then \code{NA} is specified. Otherwise, lists need to be assigned to
#' \code{mmvar} and \code{weights}, with the former containing columns
#' specifying the classification units, and the latter containing columns
#' specifying the weights. Ignored if \code{EstM = 0}, i.e. only applicable to models estimated via
#' MCMC. \code{mm = NULL} by default. Supersedes deprecated \code{xclass}.
#' E.g. (from \code{demo(MCMCGuide16)}) for
#' \code{logearn ~ 1 + age_40 + sex + parttime + (company|1) + (id|1)}, if
#' \code{company} is a multiple membership classification with the variables
#' indicating the classifications in \code{company}, \code{company2},
#' \code{company3}, \code{company4} and their weights in \code{weight1}, \code{weight2},
#' \code{weight3} and \code{weight4} then
#' \code{mm = list(list(mmvar = list('company', 'company2', 'company3', 'company4'),}
#' \code{weights = list('weight1', 'weight2', 'weight3', 'weight4')), NA)}
#' with the \code{NA}, listed last, corresponding to the level 1 identifier (\code{id}).
#' @param car A list specifying structure of a conditional autoregressive (CAR)
#' model. Each element of the list corresponds to a level (classification) of
#' the model, in descending order. If a level is not a spatial classification,
#' then \code{NA} is specified. Otherwise, lists need to be assigned to
#' \code{carvar} and \code{weights}, with the former containing columns
#' specifying the spatial classification units, and the latter containing
#' columns specifying the weights. See \code{demo(MCMCGuide17)} for examples.
#' \code{car = NULL} by default.
#' @param BUGO If non-\code{NULL} uses BUGS for MCMC estimation using files
#' specified in \code{modelfile}, \code{initfile} and \code{datafile}.
#' @param mem.init A vector which sets and displays worksheet capacities for
#' the current MLwiN session according to the value(s) specified. By default,
#' the number of levels is \code{nlev}+1; worksheet size in thousands of cells
#' is 6000; the number of columns is 2500; the number of explanatory variables
#' is \code{num_vars}+10; the number of group labels is 20. \code{nlev} is the
#' number of levels specified by \code{levID}, and \code{num_vars} is
#' approximately the number of explanatory variables calculated initially.
#' @param optimat This option instructs MLwiN to limit the maximum matrix size
#' that can be allocated by the (R)IGLS algorithm. Specify \code{optimat = TRUE}
#' if MLwiN gives the following error message 'Overflow allocating smatrix'.
#' This error message arises if one more higher-level units is extremely large
#' (contains more than 800 lower-level units).  In this situation runmlwin's
#' default behaviour is to instruct MLwiN to allocate a larger matrix size to
#' the (R)IGLS algorithm than is currently possible. Specifying
#' \code{optimat = TRUE} caps the maximum matrix size at 800 lower-level units,
#' circumventing the MLwiN error message, and allowing most MLwiN
#' functionality.
#' @param modelfile A file name where the BUGS model will be saved in .txt
#' format.
#' @param initfile A file name where the BUGS initial values will be saved
#' in .txt format.
#' @param datafile A file name where the BUGS data will be saved in .txt
#' format.
#' @param macrofile A file name where the MLwiN macro file will be saved.
#' @param IGLSfile A file name where the IGLS estimates will be saved.
#' @param MCMCfile A file name where the MCMC estimates will be saved.
#' @param chainfile A file name where the MCMC chains will be saved.
#' @param MIfile A file name where the missing values will be saved.
#' @param resifile A file name where the residual estimates will be saved.
#' @param resi.store A logical value to indicate if residuals are to be stored
#' (\code{TRUE}) or not (\code{FALSE}).
#' @param resioptions A string vector to specify the various residual options.
#' The \code{'variance'} option calculates the posterior variances instead of
#' the posterior standard errors; the \code{'standardised'} option calculates standardised
#' residuals.
#' @param resichains A file name where the residual chains will be saved.
#' @param FACTchainfile A file name where the factor chains will be saved.
#' @param resi.store.levs An integer vector indicating the levels at which the
#' residual chains are to be stored.
#' @param debugmode A logical value determining whether MLwiN is run in the
#' background or not. The default value is \code{FALSE}: i.e. MLwiN is run in
#' the background. If \code{TRUE} the MLwiN GUI is opened, and then pauses after the model
#' has been set-up, allowing user to check starting values; pressing 'Resume macro'
#' will then fit the model. Once fit, pressing 'Resume macro' once more will save
#' the outputs to the \code{workdir} ready to be read by \pkg{R2MLwiN}. Users can
#' instead opt to 'Abort macro' in which case the outputs are not saved to the
#' \code{workdir}. This option currently
#' works for 32 bit version of MLwiN only (automatically switches unless
#' \code{MLwiNPath} or \code{options(MLwiNPath)}
#' has been set directly to the executable).
#' @param startval A list of numeric vectors specifying the starting values
#' when using MCMC. \code{FP.b} corresponds to the estimates for the fixed
#' part; \code{FP.v} specifies the variance/covariance estimates for the fixed
#' part; \code{RP.b} specifies the variance estimates for the random part;
#' \code{RP.v} corresponds to the variance/covariance matrix of the variance
#' estimates for the random part.
#' @param dami This command outputs a complete (i.e. including non-missing
#' responses) response variable y. If \code{dami = c(0, <iter1>, <iter2>,...)} then
#' the response variables returned will be the value of y at the iterations
#' quoted (as integers \code{<iter1>, <iter2>}, etc.); these can be used for
#' multiple imputation. If \code{dami = 1} the value of y will be the mean
#' estimate from the iterations produced. \code{dami = 2} is as for \code{dami = 1}
#' but with the standard errors of the estimate additionally being stored.
#' @param namemap A mapping of column names to DTA friendly shorter names
#' @param saveworksheet A list of file names (one for each chain) used to store the
#' MLwiN worksheet after the model has been estimated.
#'
#' @details A list of other MCMC options as used in the argument
#' \code{mcmcOptions}:
#' \itemize{
#' \item \code{orth}: If \code{orth = 1}, orthogonal fixed effect
#' vectors are used; zero otherwise.
#' \item \code{hcen}: An integer specifying the
#' level where we use hierarchical centering.
#' \item \code{smcm}: If \code{smcm = 1},
#' structured MCMC is used; zero otherwise.
#' \item \code{smvn}: If \code{smvn = 1}, the
#' structured MVN framework is used; zero otherwise.
#' \item \code{paex}: A matrix of Nx2; in each row, if the second digit is \code{1}, parameter expansion
#' is used at level <the first digit>.
#' \item \code{mcco}: This
#' command allows the user to have constrained settings for the lowest level
#' variance matrix in a multivariate Normal model. If value is \code{0},
#' it estimates distinct variances for each residual error and distinct covariances
#' for each residual error pair. Four other
#' settings are currently available:\cr
#' \tabular{ll}{\code{1} \tab fits stuctured errors with a common correlation paramater and a common variance parameter;\cr
#' \code{2} \tab fits AR1 errors with a common variance parameter;\cr \code{3} \tab fits structured errors with a common
#' correlation parameter and independent variance parameters;\cr \code{4} \tab fits AR1 errors with independent variance
#' parameters.\cr }
#' }
#'
#' A list of objects specified for cross-classified and/or multiple membership
#' models, as used in the argument \code{xclass}:
#' \itemize{
#' \item \code{class}: An integer
#' (vector) of the specified class(es).
#' \item \code{N1}: This defines a multiple
#' membership across \code{N1} units at level \code{class}. \code{N1}>1 if
#' there is multiple membership.
#' \item \code{weight}: If there is multiple
#' membership then the column number \code{weight}, which is the length of the
#' dataset, will contain the first set of weights for the multiple membership.
#' Note that there should be \code{N1} weight columns and they should be
#' sequential in the worksheet starting from \code{weight}.
#' \item \code{id}: If the
#' response is multivariate then the column number \code{id} must be input and
#' this contains the first set of identifiers for the classification. Note that
#' for a p-variate model each lowest level unit contains p records and the
#' identifiers (sequence numbers) for each response variate need to be
#' extracted into \code{id} and following columns. There should be \code{N1} of
#' these identifier columns and they should be sequential starting from
#' \code{id} in the multivariate case.
#' \item \code{car}: \code{car = TRUE} indicates
#' the spatial CAR model; \code{FALSE} otherwise. \code{car = FALSE} if ignored.
#' }
#'
#' A list of objects specified for factor analysis, as used in the argument
#' \code{fact}:
#' \itemize{
#' \item \code{nfact}: Specifies the number of factors
#' \item \code{lev.fact}: Specifies the level/classification for the random part of
#' the factor for each factor.
#' \item \code{nfactcor}: Specifies the number of
#' correlated factors
#' \item \code{factcor}: a vector specifying the correlated
#' factors: the first element corresponds to the first factor number, the
#' second to the second factor number, the third element corresponds to the
#' starting value for the covariance and the fourth element to whether this
#' covariance is constrained
#' (\code{1}) or not (\code{0}). If more than one pair of factors is correlated,
#' then repeat this sequence for each pair.
#' \item \code{loading}: A matrix specifying the
#' starting values for the factor loadings and the starting value of the factor
#' variance. Each row corresponds to a factor.
#' \item \code{constr}: A matrix
#' specifying indicators of whether the factor loadings and the factor variance
#' are constrained (\code{1}) or not (\code{0}).
#' }
#'
#' @return Outputs a modified version of namemap containing newly generated
#' short names.
#' @note Note that for \code{FixM}, \code{residM}, \code{Lev1VarM} and
#'
#' \code{OtherVarM}, not all combinations of methods are available for all sets
#' of parameters and all models.
#'
#' @references
#' Goldstein, H. (2011) Multilevel Statistical Models. 4th Edition. London: John Wiley and Sons.
#'
#' Rasbash, J., Steele, F., Browne, W.J. and Goldstein, H. (2012)
#' A User's Guide to MLwiN Version 2.26. Centre for Multilevel Modelling,
#' University of Bristol.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso \code{\link{write.IGLS}}
#'
write.MCMC <- function(indata, dtafile, oldsyntax = FALSE, resp, levID, expl, rp, D = "Normal", nonlinear = c(0, 1), categ = NULL,
                         notation = NULL, nonfp = NULL, clre = NULL, Meth = 1, merr = NULL, carcentre = FALSE, maxiter = 20,
                         convtol = 2, seed = 1, iterations = 5000, burnin = 500, scale = 5.8, thinning = 1, priorParam = "default", refresh = 50,
                         fixM = 1, residM = 1, Lev1VarM = 1, OtherVarM = 1, adaption = 1, priorcode = c(gamma=1), rate = 50, tol = 10, lclo = 0,
                         mcmcOptions, fact = NULL, xc = NULL, mm = NULL, car = NULL, BUGO = NULL, mem.init = "default", optimat = FALSE,
                         modelfile, initfile, datafile, macrofile, IGLSfile, MCMCfile, chainfile, MIfile, resifile, resi.store = FALSE,
                         resioptions, resichains, FACTchainfile, resi.store.levs = NULL, debugmode = FALSE, startval = NULL, dami = NULL,
                         namemap = sapply(colnames(indata), as.character), saveworksheet = NULL) {
  
  shortname <- function(...) {
    name <- paste0(...)
    if (!name %in% names(namemap)) {
      sname <- paste0("v", digest(name, algo = "xxhash64", serialize = FALSE))
      names(sname) <- name
      namemap <<- c(namemap, sname)
    }
    return(namemap[[name]])
  }
  
  nlev <- length(levID)
  
  nrp <- length(rp)
  if (nrp > 0) {
    rp.names <- names(rp)
    if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
      for (i in 1:nrp) {
        temp <- rp.names[i]
        rp.names[i] <- paste("rp", as.integer(sub("rp", "", temp)) + 1, sep = "")
      }
    }
    names(rp) <- rp.names
  }
  num.expl.init <- function(p, nonfp, categ) {
    num_vars <- 0
    if (is.na(nonfp[1])) {
      if (is.null(categ) || sum(p == categ["var", ]) == 0) {
        num_vars <- num_vars + 1
      } else {
        if (is.na(categ["ref", which(p == categ["var", ])])) {
          num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])])
        } else {
          num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1
        }
      }
    } else {
      if (sum(p == nonfp) == 0) {
        if (is.null(categ) || sum(p == categ["var", ]) == 0) {
          num_vars <- num_vars + 1
        } else {
          if (is.na(categ["ref", which(p == categ["var", ])])) {
            num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])])
          } else {
            num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1
          }
        }
      }
    }
    num_vars
  }
  
  if (is.list(expl)) {
    sep.coeff <- expl$sep.coeff
    if (is.list(nonfp)) {
      nonfp.sep <- nonfp$nonfp.sep
      nonfp.common <- nonfp$nonfp.common
    }
    if (!is.na(sep.coeff[1])) {
      num_vars <- sum(sapply(sep.coeff, function(x) num.expl.init(x, nonfp.sep, categ)))
    } else {
      num_vars <- 0
    }
    
    if (D[1] == "Multinomial") {
      nresp <- length(levels(indata[, resp])) - 1
      ## it is true if adding each expl variables separately
      num_vars <- num_vars * nresp
    }
    if (D[1] == "Multivariate Normal") {
      nresp <- length(resp)
      num_vars <- num_vars * nresp
    }
    
    common.coeff <- expl$common.coeff
    num_vars <- num_vars + sum(sapply(common.coeff, function(x) num.expl.init(x, nonfp$nonfp.common, categ)))
    common.coeff.id <- expl$common.coeff.id
  } else {
    num_vars <- sum(sapply(expl, function(x) num.expl.init(x, nonfp, categ)))
    if (D[1] == "Multinomial") {
      nresp <- length(levels(indata[, resp])) - 1
      num_vars <- num_vars * nresp
    }
    if (D[1] == "Multivariate Normal") {
      nresp <- length(resp)
      num_vars <- num_vars * nresp
    }
  }
  
  if (nrp > 0) {
    for (ii in 1:nrp) {
      rpx <- rp[[rp.names[ii]]]
      nrpx <- length(rpx)
      if (nrpx == 1)
        num_vars <- num_vars + 1
      if (nrpx >= 2)
        num_vars <- num_vars + nrpx * (nrpx - 1)/2 + nrpx
    }
  }
  
  ## Write into macro file
  wrt <- function(x) write(x, macrofile, append = TRUE)
  
  cat(file = macrofile)
  wrt("ECHO    0")
  wrt("NOTE    *****************************************************************")
  wrt("NOTE      MLwiN macro created by rtomlwin command")
  wrt("NOTE    *****************************************************************")
  wrt("")
  
  wrt("NOTE    Initialise MLwiN storage")
  if (mem.init[1] == "default") {
    wrt(paste("INIT    ", nlev + 1, " 6000 2500 ", num_vars + 10, " 30", sep = ""))
  } else {
    wrt(paste("INIT    ", mem.init[1], mem.init[2], mem.init[3], mem.init[4], mem.init[5]))
  }
  wrt("NOTE     Limit the maximum matrix size")
  if (optimat) {
    wrt("OPTS   1")
  } else {
    wrt("OPTS   0")
  }
  
  wrt("MONI    0")
  wrt("MARK    0")
  wrt("NOTE    Import the R data into MLwiN")
  wrt(paste("RSTA    '", dtafile, "'", sep = ""))
  
  wrt("NOTE    Correct column names")
  for (name in names(namemap)) {
    wrt(paste0("COLN '", namemap[[name]], "' b50"))
    wrt(paste0("DESC cb50 '", name, "'"))
    wrt(paste0("NAME cb50 '", name, "'"))
  }
  
  if ("class" %in% notation) {
    wrt("INDE 1")
  }
  if (!(D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed")) {
    wrt("NOTE   Specify the response variable")
    for (p in resp) wrt(paste("RESP    '", p, "'", sep = ""))
    wrt("")
  }
  
  if (D[1] == "Binomial") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii]))
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
      
    }
    wrt("")
    
    wrt("RDISt 1 0")
    wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    if (D[2] == "logit") {
      wrt("LFUN 0")
      DD2 <- 0
    }
    if (D[2] == "probit") {
      wrt("LFUN 1")
      DD2 <- 1
    }
    if (D[2] == "cloglog") {
      wrt("LFUN 2")
      DD2 <- 2
    }
    
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  
  if (D[1] == "Mixed") {
    nresp <- length(resp)
    for (ii in 1:nresp) wrt(paste("MVAR 1   '", resp[ii], "'", sep = ""))
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:2
      if (!is.na(levID[ii]))
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("IDEN 1 'resp_indicator'")
    jj <- 1
    for (ii in 2:length(D)) {
      if (D[[ii]][1] == "Binomial") {
        wrt(paste("RDISt ", jj, " 0", sep = ""))
        wrt(paste("DOFFs ", jj, " '", D[[ii]][3], "'", sep = ""))
        if (D[[ii]][2] == "logit") {
          wrt("LFUN 0")
          DD2 <- 0
        }
        if (D[[ii]][2] == "probit") {
          wrt("LFUN 1")
          DD2 <- 1
        }
        if (D[[ii]][2] == "cloglog") {
          wrt("LFUN 2")
          DD2 <- 2
        }
      }
      if (D[[ii]][1] == "Poisson") {
        wrt(paste("RDISt ", jj, " 1", sep = ""))
        wrt("LFUN 3")
        DD2 <- 3
        DD2 <- 3
        if (!is.na(D[[jj]][3])) {
          wrt(paste("DOFFs ", jj, " '", D[[jj]][3], "'", sep = ""))
        }
      }
      jj <- jj + 1
    }
    wrt("")
    if (is.list(expl)) {
      if (!is.na(sep.coeff[1])) {
        interpos1 <- grep("\\:", sep.coeff)
        if (!isTRUE(oldsyntax) || length(interpos1) == 0) {
          for (x in 1:length(sep.coeff)) {
            p <- sep.coeff[x]
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        } else {
          exply <- sep.coeff[interpos1]
          explx <- sep.coeff[-interpos1]
          if (length(explx) > 0) {
            for (p in explx) {
              if (is.null(categ)) {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              } else {
                if (sum(p == categ["var", ]) != 0) {
                  if (is.na(categ["ref", which(p == categ["var", ])])) {
                    wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                  } else {
                    wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                               ])]), sep = ""))
                  }
                } else {
                  wrt(paste("ADDT    '", p, "'", sep = ""))
                }
              }
            }
          }
          for (i in 1:length(exply)) {
            TT <- ""
            interx <- unlist(strsplit(exply[i], "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          sep.coeff <- c(explx, exply)
        }
      }
      interpos2 <- grep("\\:", common.coeff)
      if (!isTRUE(oldsyntax) || length(interpos2) == 0) {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          
          p <- common.coeff[y]
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          p <- common.coeff[y]
          if (!(y %in% interpos2)) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          } else {
            TT <- ""
            interx <- unlist(strsplit(p, "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          wrt("RPAT")
          common.coeff[y] <- paste(p, ".", bb, sep = "")
        }
      }
      
    } else {
      interpos <- grep("\\:", expl)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in expl) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        exply <- expl[interpos]
        explx <- expl[-interpos]
        if (length(explx) > 0) {
          for (p in explx) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        }
        for (i in 1:length(exply)) {
          TT <- ""
          interx <- unlist(strsplit(exply[i], "\\:"))
          for (j in 1:length(interx)) {
            TT <- paste(TT, "'", interx[j], "' ", sep = "")
          }
          wrt(paste("ADDT    ", TT, sep = ""))
        }
        expl <- c(explx, exply)
      }
    }
    wrt("")
    
  }
  
  if (D[1] == "Multivariate Normal") {
    nresp <- length(resp)
    for (ii in 1:nresp) wrt(paste("MVAR 1   '", resp[ii], "'", sep = ""))
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:2
      if (!is.na(levID[ii]))
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("IDEN 1 'resp_indicator'")
    wrt("LFUN 0")
    wrt("")
    if (is.list(expl)) {
      if (!is.na(sep.coeff[1])) {
        interpos1 <- grep("\\:", sep.coeff)
        if (!isTRUE(oldsyntax) || length(interpos1) == 0) {
          for (x in 1:length(sep.coeff)) {
            p <- sep.coeff[x]
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        } else {
          exply <- sep.coeff[interpos1]
          explx <- sep.coeff[-interpos1]
          if (length(explx) > 0) {
            for (p in explx) {
              if (is.null(categ)) {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              } else {
                if (sum(p == categ["var", ]) != 0) {
                  if (is.na(categ["ref", which(p == categ["var", ])])) {
                    wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                  } else {
                    wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                               ])]), sep = ""))
                  }
                } else {
                  wrt(paste("ADDT    '", p, "'", sep = ""))
                }
              }
            }
          }
          for (i in 1:length(exply)) {
            TT <- ""
            interx <- unlist(strsplit(exply[i], "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          sep.coeff <- c(explx, exply)
        }
      }
      interpos2 <- grep("\\:", common.coeff)
      if (!isTRUE(oldsyntax) || length(interpos2) == 0) {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          
          p <- common.coeff[y]
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          p <- common.coeff[y]
          if (!(y %in% interpos2)) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          } else {
            TT <- ""
            interx <- unlist(strsplit(p, "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          wrt("RPAT")
          common.coeff[y] <- paste(p, ".", bb, sep = "")
        }
      }
      
    } else {
      interpos <- grep("\\:", expl)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in expl) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        exply <- expl[interpos]
        explx <- expl[-interpos]
        if (length(explx) > 0) {
          for (p in explx) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        }
        for (i in 1:length(exply)) {
          TT <- ""
          interx <- unlist(strsplit(exply[i], "\\:"))
          for (j in 1:length(interx)) {
            TT <- paste(TT, "'", interx[j], "' ", sep = "")
          }
          wrt(paste("ADDT    ", TT, sep = ""))
        }
        expl <- c(explx, exply)
      }
    }
    wrt("")
  }
  
  
  
  if (D[1] == "Multinomial") {
    wrt("LINK 2 G21")
    wrt("NAME G21[1] 'resp' G21[2] 'resp_indicator'")
    wrt("LINK 0 G21")
    wrt(paste("MNOM ", as.numeric(D[4]), " '", resp, "' ", "'resp' 'resp_indicator' ", as.numeric(D[5]), sep = ""))
    wrt("RESP   'resp'")
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:c(nlev - 1)) {
      aa <- nlev:1
      if (!is.na(levID[ii]))
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("IDEN 1 'resp_indicator'")
    wrt("")
    
    if (as.numeric(D[4]) == 0)
      wrt("RDISt 1 4") else wrt("RDISt 1 5")
    wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    if (D[2] == "logit") {
      wrt("LFUN 0")
      DD2 <- 0
    }
    if (D[2] == "probit") {
      wrt("LFUN 1")
      DD2 <- 1
    }
    if (D[2] == "cloglog") {
      wrt("LFUN 2")
      DD2 <- 2
    }
    
    wrt("")
    if (is.list(expl)) {
      if (!is.na(sep.coeff[1])) {
        interpos1 <- grep("\\:", sep.coeff)
        if (!isTRUE(oldsyntax) || length(interpos1) == 0) {
          for (x in 1:length(sep.coeff)) {
            p <- sep.coeff[x]
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        } else {
          exply <- sep.coeff[interpos1]
          explx <- sep.coeff[-interpos1]
          if (length(explx) > 0) {
            for (p in explx) {
              if (is.null(categ)) {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              } else {
                if (sum(p == categ["var", ]) != 0) {
                  if (is.na(categ["ref", which(p == categ["var", ])])) {
                    wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                  } else {
                    wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                               ])]), sep = ""))
                  }
                } else {
                  wrt(paste("ADDT    '", p, "'", sep = ""))
                }
              }
            }
          }
          for (i in 1:length(exply)) {
            TT <- ""
            interx <- unlist(strsplit(exply[i], "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          sep.coeff <- c(explx, exply)
        }
      }
      interpos2 <- grep("\\:", common.coeff)
      if (!isTRUE(oldsyntax) || length(interpos2) == 0) {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          
          p <- common.coeff[y]
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          p <- common.coeff[y]
          if (!(y %in% interpos2)) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          } else {
            TT <- ""
            interx <- unlist(strsplit(p, "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          wrt("RPAT")
          common.coeff[y] <- paste(p, ".", bb, sep = "")
        }
      }
    } else {
      interpos <- grep("\\:", expl)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in expl) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        exply <- expl[interpos]
        explx <- expl[-interpos]
        if (length(explx) > 0) {
          for (p in explx) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        }
        for (i in 1:length(exply)) {
          TT <- ""
          interx <- unlist(strsplit(exply[i], "\\:"))
          for (j in 1:length(interx)) {
            TT <- paste(TT, "'", interx[j], "' ", sep = "")
          }
          wrt(paste("ADDT    ", TT, sep = ""))
        }
        expl <- c(explx, exply)
      }
    }
    wrt("")
  }
  
  if (D[1] == "Normal") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii]))
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("")
    
    wrt("NOTE   Specify covariate(s) used anywhere in the model")
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  if (D[1] == "Poisson") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii]))
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("")
    
    wrt("RDISt 1 1")
    wrt("LFUN 3")
    DD2 <- 3
    if (!is.na(D[3])) {
      wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    }
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  
  
  if (D[1] == "Negbinom") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii]))
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("")
    
    wrt("RDISt 1 2")
    wrt("LFUN 3")
    DD2 <- 3
    if (!is.na(D[3])) {
      wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    }
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var",
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  
  if (is.list(nonfp)) {
    wrt("NOTE Turn off the fixed part of the explanatary variable(s)")
    nonfp.sep <- nonfp$nonfp.sep
    nonfp.common <- nonfp$nonfp.common
    if (!is.na(nonfp.sep[1])) {
      interpos <- grep("\\:", nonfp.sep)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in nonfp.sep) wrt(paste("FPAR 0  '", p, "'", sep = ""))
      } else {
        for (i in 1:length(nonfp.sep)) {
          if (i %in% interpos && oldsyntax) {
            wrt(paste("FPAR 0  '", gsub("\\:", "\\.", nonfp.sep[i]), "'", sep = ""))
          } else {
            wrt(paste("FPAR 0  '", nonfp.sep[i], "'", sep = ""))
          }
        }
      }
    }
    if (!is.na(nonfp.common[1])) {
      interpos <- grep("\\:", nonfp.common)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in nonfp.common) wrt(paste("FPAR 0  '", p, "'", sep = ""))
      } else {
        for (i in 1:length(nonfp.common)) {
          if (i %in% interpos && oldsyntax) {
            wrt(paste("FPAR 0  '", gsub("\\:", "\\.", nonfp.common[i]), "'", sep = ""))
          } else {
            wrt(paste("FPAR 0  '", nonfp.common[i], "'", sep = ""))
          }
        }
      }
    }
  } else {
    if (!is.na(nonfp[1])) {
      wrt("NOTE Turn off the fixed part of the explotary varible(s)")
      for (p in nonfp) {
        if(oldsyntax){
          p <- gsub("\\:", "\\.", p)
        }
        wrt(paste("FPAR 0  '", p, "'", sep = ""))
      }
    }
  }
  
  wrt("")
  wrt("NOTE   Specify random part covariate(s)")
  if (nrp > 0) {
    for (ii in 1:nrp) {
      for (p in rp[[ii]]) {
        if (oldsyntax) {
          p <- gsub("\\:", "\\.", p)
        }
        wrt(paste("SETV  ", as.numeric(sub("rp", "", rp.names[ii])), "   '", p, "'", sep = ""))
      }
    }
  }
  if (!is.null(clre)) {
    nclre <- ncol(clre)
    for (ii in 1:nclre) {
      if (oldsyntax){
        wrt(paste("CLRE  ", as.numeric(clre[1, ii]), " '", gsub("\\:", "\\.", clre[2, ii]), "' '",
                  gsub("\\:", "\\.", clre[3, ii]), "'", sep = ""))
      } else {
        wrt(paste("CLRE  ", as.numeric(clre[1, ii]), " '", clre[2, ii], "' '", clre[3, ii], "'", sep = ""))
      }
    }
  }
  
  nexpl <- length(expl)
  wrt("")
  
  wrt("NOTE   Set estimation method")
  if (Meth != 2) {
    wrt(paste("METH", Meth))
  }
  wrt(paste("LINE ", nonlinear[1], nonlinear[2]))
  wrt("")
  
  if (!is.null(fact)) {
    TT <- NULL
    for (i in 1:fact$nfact) {
      TT <- c(TT, fact$lev.fact[i] + 1, matrix(rbind(fact$loading[i, ], fact$constr[i, ]), nrow = 1))
    }
    if (fact$nfactcor > 0)
      TT <- c(TT, fact$factcor)
    FACT <- as.vector(c(length(resp), fact$nfact, fact$nfactcor, TT))
    rm(TT)
    TT <- ""
    for (i in 1:length(FACT)) {
      TT <- paste(TT, FACT[i])
    }
    wrt(paste("FACT ", TT))
    wrt("LINK 2 G21")
    wrt("SMFA 1 G21[1]")
    wrt("SMFA 2 G21[2]")
    
  }
  if (D[1] == "Normal") {
    wrt("PREF   0")
    wrt("POST   0")
  }

  if ("simple" %in% notation) {
    wrt("EXISt 'cons' b1000")
    wrt("SWITch b1000")
    wrt("CASE 0:")
    wrt("EXISt 'Intercept' b1000")
    wrt("SWITch b1000")
    wrt("CASE 1:")
    wrt("COLN 'Intercept' b1000")
    wrt("NAME cb1000 'cons'")
    wrt("ENDSWITch")
    wrt("ENDSWITch")
    wrt("NOTA 1")
  }

  wrt("NAME   c1098 '_FP_b'")
  wrt("NAME   c1099 '_FP_v'")
  wrt("NAME   c1096 '_RP_b'")
  wrt("NAME   c1097 '_RP_v'")
  
  wrt("NOTE   Fit the model")
  wrt("ECHO 1")
  wrt("BATC 1")
  wrt("MAXI 2")
  wrt("STAR")
  
  if (!is.null(startval)) {
    if (!is.null(startval$FP.b)) {
      wrt(paste("JOIN ", paste(startval$FP.b, collapse = " "), " '_FP_b'", sep = ""))
    }
    if (!is.null(startval$FP.v)) {
      wrt(paste("JOIN ", paste(startval$FP.v[!upper.tri(startval$FP.v)], collapse = " "), " '_FP_v'", sep = ""))
    }
    if (!is.null(startval$RP.b)) {
      wrt(paste("JOIN ", paste(startval$RP.b, collapse = " "), " '_RP_b'", sep = ""))
      if (D[1] == "Multinomial" && D[4] == 0) {
        wrt("JOIN '_RP_b' 1 '_RP_b'")
      }
    }
    if (!is.null(startval$RP.v)) {
      wrt(paste("JOIN ", paste(startval$RP.v[!upper.tri(startval$RP.v)], collapse = " "), " '_RP_v'", sep = ""))
      if (D[1] == "Multinomial" && D[4] == 0) {
        wrt(paste("JOIN '_RP_v' ", paste(rep(0, length(startval$RP.b)+1), collapse = " "), " '_RP_v'"))
      }
    }
  } else {
    wrt(paste("TOLE", convtol))
    wrt(paste("MAXI", maxiter))
    wrt("BATC 1")
    wrt("NEXT")
  }
  wrt("ECHO 0")
  wrt("MONI 1")
  wrt("ITNU 0 b21")
  wrt("CONV b22")
  wrt("")
  
  wrt("NOTE    *****************************************************************")
  wrt("")
  
  wrt("NOTE    *****************************************************************")
  wrt("NOTE       Export the model results to R")
  wrt("NOTE    *****************************************************************")
  
  wrt("LINK 1 G30")
  wrt("NAME   G30[1] '_Stats'")
  if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
    wrt("NOBS 2 b31 b32")
  } else {
    wrt("NOBS 1 b31 b32")
  }
  wrt("EDIT 1 '_Stats' b31")
  wrt("EDIT 2 '_Stats' b32")
  if (D[1] == "Normal" || D[1] == "Multivariate Normal") {
    wrt("LIKE   b100")
  }
  wrt("EDIT 3 '_Stats' b100")
  wrt("EDIT 7 '_Stats' b21")
  wrt("EDIT 8 '_Stats' b22")
  wrt("NAME   c1094 '_esample'")
  wrt("SUM '_esample' b1")
  wrt("EDIT 9 '_Stats' b1")
  wrt(paste("PSTA '", IGLSfile, "' ", "'_FP_b' ", "'_FP_v' ", "'_RP_b' ", "'_RP_v' ", "'_Stats'", sep = ""))
  wrt("LINK 0 G30")
  
  wrt("NOTE    *****************************************************************")
  wrt("NOTE Set estimation method to MCMC")
  wrt("NOTE    *****************************************************************")
  wrt("EMODe  3")
  
  if (!is.null(resi.store.levs)) {
    wrt(paste0("LINK ", length(resi.store.levs), " G22"))
    for (i in 1:length(resi.store.levs)) {
      resilev <- resi.store.levs[i]
      if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
        resilev <- resilev + 1
      }
      wrt(paste("SMRE ", resilev, " G22[", i, "]", sep = ""))
      wrt(paste("NAME ", " G22[", i, "] 'resi_lev", resi.store.levs[i], "'", sep = ""))
    }
    
  }
  if (!is.null(merr)) {
    nmerr <- as.numeric(merr[1])
    tt <- paste("MERR  ", nmerr)
    j <- 1
    for (ii in 1:nmerr) {
      tt <- paste(tt, " '", merr[j + 1], "' ", as.numeric(merr[j + 2]), sep = "")
      j <- j + 2
    }
    wrt(tt)
  }
  wrt(paste("ORTH", mcmcOptions$orth))
  if (mcmcOptions$hcen == 0) {
    wrt(paste("HCEN", mcmcOptions$hcen))
  } else {
    wrt(paste("HCEN", 1, mcmcOptions$hcen))
  }
  wrt(paste("SMCM", mcmcOptions$smcm))
  wrt(paste("SMVN", mcmcOptions$smvn))
  if (is.matrix(mcmcOptions$paex)) {
    apply(mcmcOptions$paex, 1, function(x) wrt(paste("PAEX", x[1], x[2])))
  } else {
    wrt(paste("PAEX", mcmcOptions$paex[1], mcmcOptions$paex[2]))
  }
  
  if (D[1] == "Multivariate Normal")
    wrt(paste("MCCO ", mcmcOptions$mcco))
  
  if (!is.null(xc) && isTRUE(xc)) {
    wrt("XCLA 1")
  }
  
  if (!is.null(mm)) {
    for (i in 1:length(mm)) {
      if (!any(is.na(mm[[i]]))) {
        mmlev <- (length(mm) - i) + 1
        if (D[1] != "Multivariate Normal") {
          wrt(paste0("MULM ", mmlev, " ", length(mm[[i]]$mmvar), " '", mm[[i]]$weights[[1]], "'"))
        } else {
          wrt(paste0("MULM ", mmlev, " ", length(mm[[i]]$mmvar), " '", mm[[i]]$weights[[1]], "' '", mm[[i]]$mmvars[[1]],
                     "'"))
        }
      }
    }
  }
  
  if (!is.null(car)) {
    for (i in 1:length(car)) {
      if (!any(is.na(car[[i]]))) {
        carlev <- (length(car) - i) + 1
        wrt(paste0("MULM ", carlev, " ", length(car[[i]]$carvar), " '", car[[i]]$weights[[1]], "' '", car[[i]]$carvar[[1]],
                   "'"))
        wrt(paste("CARP", carlev, "1"))
      }
    }
    if (!is.null(carcentre) && isTRUE(carcentre)) {
      wrt("CARC 1")
    }
  }
  wrt("")
  
  wrt("NOTE Set MCMC seed")
  wrt(paste("MCRS ", seed, sep = ""))
  wrt("")
  
  wrt("NOTE Set prior distribution parameters")
  if (priorParam[1] != "default") {
    wrt("PRIOR  c1092")
    lenpp <- length(priorParam)
    tempt <- " "
    for (i in 1:lenpp) tempt <- paste(tempt, priorParam[i])
    wrt(paste("JOIN", tempt, " c1092", sep = ""))
    wrt("")
  } else {
    wrt("")
  }
  
  if (D[1] == "Normal") {
    wrt("PREF   0")
    wrt("POST   0")
  }
  
  if (nlev > 1 && !isTRUE(xc)) {
    wrt("MISR   0")
    wrt("LINK 2 G30")
    if (D[1] == "Multinomial" && as.numeric(D[4]) == 0) {
      len.rpx <- 2
      wrt(paste0("LINK ", len.rpx, " G27"))
      wrt("LINK 1 G28")
      wrt("NOTE Calculate MCMC starting values for level 2 residuals")
      wrt("RLEV 2")
      wrt("RFUN")
      wrt("RCOV 2")
      wrt("ROUT G27 G28")
      wrt("RESI")
      wrt("JOIN G30[1] G27 G30[1]")
      wrt("JOIN G30[2] G28 G30[2]")
      wrt("ERAS G27")
      wrt("ERAS G28")
      wrt("LINK 0 G27")
      wrt("LINK 0 G28")
    }
    if (nrp > 0) {
      for (j in nrp:1) {
        if (as.numeric(sub("rp", "", rp.names[j])) != 1) {
          rpx <- rp[[j]]
          len.rpx <- length(rpx)
          wrt(paste0("LINK ", len.rpx, " G27"))
          wrt("LINK 1 G28")
          wrt(paste("NOTE Calculate MCMC starting values for level ", as.numeric(sub("rp", "", rp.names[j])),
                    " residuals", sep = ""))
          wrt(paste("RLEV   ", as.numeric(sub("rp", "", rp.names[j])), sep = ""))
          wrt("RFUN")
          wrt("RCOV   2")
          wrt("ROUT G27 G28")
          wrt("RESI")
          wrt("JOIN G30[1] G27 G30[1]")
          wrt("JOIN G30[2] G28 G30[2]")
          wrt("ERAS G27")
          wrt("ERAS G28")
          wrt("LINK 0 G27")
          wrt("LINK 0 G28")
        }
      }
    }
    wrt("MISR   1")
  }
  
  if (D[1] == "Normal")
    DD <- 1
  if (D[1] == "Binomial")
    DD <- 2
  if (D[1] == "Mixed")
    DD <- 5
  if (D[1] == "Poisson")
    DD <- 3
  if (D[1] == "Multivariate Normal")
    DD <- 4
  if (D[1] == "Multinomial") {
    if (as.numeric(D[4]) == 0)
      DD <- 6 else DD <- 7
    wrt("CLRV 2")
  }
  if (D[1] == "Negbinom")
    DD <- 8
  
  wrt(paste("LCLO   ", lclo, sep = ""))

  if (priorcode["gamma"] == 0) {
    wrt("UNIP")
  } else {
    if (length(priorcode) > 1) {
      wrt(paste("GAMP", priorcode["shape"], priorcode["scale"]))
    }
  }
  
  if (debugmode) {
    wrt("NOTE   Open the equations window")
    wrt("WSET 15 1")
    wrt("EXPA 3")
    wrt("ESTM 2")
    wrt("PAUS")
  }
  
  priorcol <- ""
  if (priorParam[1] != "default") {
    priorcol <- "c1092"
  }
  
  if ((!is.null(BUGO)) && !(D[1] == "Mixed") && nrp > 0) {
    if (D[1] == "Normal" || D[1] == "Multivariate Normal")
      DD2 <- 0
    if (nlev > 1 && !isTRUE(xc)) {
      wrt(paste("BUGO 6 ", DD, " ", DD2, " G30[1] ", priorcol, " '", modelfile, "' ", "'", initfile, "' ", "'",
                datafile, "'", sep = ""))
      wrt("ERAS   G30")
      wrt("LINK 0 G30")
    } else {
      wrt(paste("BUGO 6 ", DD, " ", DD2, " ", "'", modelfile, "' ", "'", initfile, "' ", "'", datafile, "'",
                sep = ""))
    }
  } else {
    wrt("NOTE   fit the model in MCMC")
    wrt(paste("MTOT   ", iterations, sep = ""))
    wrt("ECHO 1")
    
    residcols <- ""
    if (nlev > 1 && !isTRUE(xc)) {
      residcols <- "G30[1] G30[2]"
    }
    wrt(paste("MCMC   0", burnin, adaption, scale, rate, tol, residcols, priorcol, fixM, residM, Lev1VarM, OtherVarM,
              priorcode[1], DD))
    
    if (nlev > 1 && !isTRUE(xc)) {
      wrt("ERAS G30")
      wrt("LINK 0 G30")
    }
    wrt("ERAS  c1090 c1091")
    wrt("")
    if (!is.null(dami) && dami[1] == 0 && length(dami) > 1) {
      ndami <- length(dami)
      mvnames <- rep(NA, ndami - 1)
      wrt(paste("LINK", ndami, "G23"))
      for (i in 2:ndami) {
        wrt(paste("MCMC 1 ", dami[i] - dami[i - 1], " ", thinning, " c1090 c1091 c1003 c1004 1 ", DD, sep = ""))
        wrt("PUPN c1003 c1004")
        wrt("AVER c1091 b99 b100")
        wrt(paste0("DAMI 0 G23[", i - 1, "]"))
        mvnames[i - 1] <- paste0("'_est_", dami[i], "'")
        wrt(paste0("NAME G23[", i - 1, "] ", mvnames[i - 1]))
        wrt("PAUS 1")
      }
      wrt("LINK 0 G23")
      if (dami[ndami] < iterations) {
        wrt(paste("MCMC 1 ", (iterations/thinning) - dami[ndami], " ", thinning, " c1090 c1091 c1003 c1004 1 ",
                  DD, sep = ""))
        wrt("PUPN c1003 c1004")
        wrt("AVER c1091 b99 b100")
        wrt("PAUS 1")
      }
      wrt(paste("PSTA", paste0("'", MIfile, "'"), paste(mvnames, collapse = " "), "'resp_indicator'"))
      wrt(paste("ERAS  ", paste(mvnames, collapse = " "), sep = ""))
    } else {
      if (debugmode) {
        for (i in 1:floor((iterations/thinning)/refresh)) {
          wrt(paste("MCMC 1 ", refresh, " ", thinning, " c1090 c1091 c1003 c1004 1 ", DD, sep = ""))
          wrt("PUPN c1003 c1004")
          wrt("AVER c1091 b99 b100")
          wrt("PAUS 1")
        }
        is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
        if (!is.wholenumber(iterations/refresh)) {
          wrt(paste("MCMC 1 ", (iterations/thinning)%%refresh, " ", thinning, " c1090 c1091 c1003 c1004 1 ",
                    DD, sep = ""))
          wrt("PUPN c1003 c1004")
          wrt("AVER c1091 b99 b100")
          wrt("PAUS 1")
        }
      } else {
        wrt(paste("MCMC 1 ", iterations/thinning, " ", thinning, " c1090 c1091 c1003 c1004 1 ", DD, sep = ""))
        wrt("PUPN c1003 c1004")
        wrt("AVER c1091 b99 b100")
      }
    }
    wrt("ECHO 0")

    if (!is.null(saveworksheet)) {
      wrt(paste0("STOR ", saveworksheet))
    }
    
    if (debugmode) {
      wrt("WSET 15 1")
      wrt("EXPA 3")
      wrt("ESTM 2")
      wrt("PAUS")
    }
    
    wrt("NOTE    *****************************************************************")
    wrt("NOTE       Export the model results to R")
    wrt("NOTE    *****************************************************************")
    if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
      wrt("NOBS 2 b31 b32")
    } else {
      wrt("NOBS 1 b31 b32")
    }
    wrt("EDIT 1 '_Stats' b31")
    wrt("EDIT 2 '_Stats' b32")
    if (!(D[1] == "Mixed") && is.null(merr) && is.null(fact)) {
      wrt("BDIC b1 b2 b3 b4")
      wrt("EDIT 3 '_Stats' b1")
      wrt("EDIT 4 '_Stats' b2")
      wrt("EDIT 5 '_Stats' b3")
      wrt("EDIT 6 '_Stats' b4")
    }
    wrt("EDIT 7 '_Stats' b21")
    wrt("EDIT 8 '_Stats' b22")
    wrt("NAME   c1098 '_FP_b'")
    wrt("NAME   c1099 '_FP_v'")
    wrt("NAME   c1096 '_RP_b'")
    wrt("NAME   c1097 '_RP_v'")
    wrt("NAME   c1094 '_esample'")
    wrt("SUM '_esample' b1")
    wrt("EDIT 9 '_Stats' b1")
    if (is.null(fact)) {
      wrt(paste("PSTA '", MCMCfile, "' ", "'_FP_b' ", "'_FP_v' ", "'_RP_b' ", "'_RP_v' ", "'_Stats'", sep = ""))
    } else {
      wrt("LINK 6 G25")
      wrt("DAFA G25[1] G25[2]")
      wrt("DAFL G25[3] G25[4]")
      wrt("DAFV G25[5] G25[6]")
      wrt("NAME G21[1] '_FACT_load_b_chain'")
      wrt("NAME G21[2] '_FACT_load_v_chain'")
      wrt("NAME G25[1] '_FACT_value_b'")
      wrt("NAME G25[2] '_FACT_value_v'")
      wrt("NAME   G25[3]  '_FACT_load_b'")
      wrt("NAME   G25[4]  '_FACT_load_v'")
      wrt("NAME   G25[5]  '_FACT_var_b'")
      wrt("NAME   G25[6]  '_FACT_var_v'")
      wrt(paste("PSTA '", FACTchainfile, "' ", "'_FACT_load_b_chain' ", "'_FACT_load_v_chain' ", "'_FACT_value_b' ",
                "'_FACT_value_v' ", sep = ""))
      wrt(paste("PSTA '", MCMCfile, "' ", "'_FP_b' ", "'_FP_v' ", "'_RP_b' ", "'_RP_v' ", "'_FACT_load_b' ",
                "'_FACT_load_v' ", "'_FACT_var_b' ", "'_FACT_var_v' ", "'_Stats'", sep = ""))
      wrt("ERAS G21")
      wrt("LINK 0 G21")
      wrt("ERAS G25")
      wrt("LINK 0 G25")
    }
    wrt("ERAS '_Stats'")
    wrt("")
    
    if (!is.null(dami) && length(dami) == 1) {
      wrt("NOTE save imputed values if there are missing values")
      wrt("SWIT b1")
      wrt("CASE 0:")
      wrt("LEAVE")
      wrt("CASE:")
      wrt("LINK 3 G26")
      wrt("NAME G26[1] '_MissingInd'")
      wrt("CALC   '_MissingInd'=abso('_esample'-1)")
      if (dami == 1) {
        wrt("DAMI 1 G26[2]")
        wrt("NAME G26[2] '_est'")
        wrt(paste("PSTA '", MIfile, "' '_est' '_MissingInd' ", sep = ""))
        wrt("ERAS  '_est'")
      }
      if (dami == 2) {
        wrt("DAMI 2 G26[2] G26[3]")
        wrt("NAME G26[2] '_est'")
        wrt("NAME G26[3] '_SDs'")
        wrt(paste("PSTA '", MIfile, "' '_est' '_SDs' '_MissingInd' ", sep = ""))
        wrt("ERAS  '_est' '_SDs'")
      }
      wrt("LINK 0 G26")
      wrt("ENDS")
      wrt("")
    }
    
    wrt("NOTE export parameter chain")
    wrt("NAME   c1091 'deviance'")
    wrt("NAME   c1090 'mcmcchains'")
    
    wrt("LINK 0 G25")
    wrt("LINK 0 G26")
    
    if (D[1] == "Multinomial") {
      nresp <- length(levels(indata[, resp])) - 1
      resp.names <- levels(indata[, resp])[-as.numeric(D[5])]
      
      if (is.list(expl)) {
        nonfp.s <- nonfp.sep
        for (i in 1:length(resp.names)) {
          if (D["mode"] == 0) {
            nonfp.s <- gsub(paste(".", resp.names[i], sep = ""), "", nonfp.s)
          }
          if (D["mode"] == 1) {
            nonfp.s <- gsub(paste(".(>=", resp.names[i], ")", sep = ""), "", nonfp.s)
          }
        }
        nonfp.s <- unique(nonfp.s)
        if (is.na(sep.coeff[1]))
          sep.coeff <- character(0)
        for (p in sep.coeff) {
          if (is.na(nonfp.sep[1]) || sum(p == nonfp.s) == 0) {
            if (is.null(categ) || sum(p == categ["var", ]) == 0) {
              for (j in 1:nresp) {
                wrt("LINK 1 G25")
                wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", p), "_", resp.names[j]), "'"))
                wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", p), "_", resp.names[j], "'"))
                wrt("GSET 2 G26 G25 G26")
                wrt("LINK 0 G25")
              }
            } else {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                categ.names <- levels(indata[[p]])
                for (j in 1:nresp) {
                  for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                    wrt("LINK 1 G25")
                    wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                    wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                    wrt("GSET 2 G26 G25 G26")
                    wrt("LINK 0 G25")
                  }
                }
              } else {
                categ.names <- levels(indata[[p]])
                refx <- categ["ref", which(p == categ["var", ])]
                categ.names <- categ.names[-which(refx == categ.names)]
                for (j in 1:nresp) {
                  for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                    wrt("LINK 1 G25")
                    wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                    wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                    wrt("GSET 2 G26 G25 G26")
                    wrt("LINK 0 G25")
                  }
                }
              }
            }
          }
        }
        # svec.common <- character(0)
        kk <- 1
        tempid <- 1:(nresp + 1)
        tempid <- tempid[-as.numeric(D["ref.cat"])]
        for (p in common.coeff) {
          newp <- paste(p, paste(tempid[as.logical(common.coeff.id[kk, ])], collapse = ""), sep = ".")
          kk <- kk + 1
          nonfp.c <- nonfp.common
          if (is.na(nonfp.common[1]) || sum(newp == nonfp.c) == 0) {
            if (is.null(categ) || sum(p == categ["var", ]) == 0) {
              wrt("LINK 1 G25")
              wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", newp)), "'"))
              wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", newp), "'"))
              wrt("GSET 2 G26 G25 G26")
              wrt("LINK 0 G25")
            } else {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                categ.names <- levels(indata[[p]])
                for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                  wrt("LINK 1 G25")
                  wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i])), "'"))
                  wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "'"))
                  wrt("GSET 2 G26 G25 G26")
                  wrt("LINK 0 G25")
                }
              } else {
                categ.names <- levels(indata[[p]])
                refx <- categ["ref", which(p == categ["var", ])]
                categ.names <- categ.names[-which(refx == categ.names)]
                for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                  wrt("LINK 1 G25")
                  wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i])), "'"))
                  wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "'"))
                  wrt("GSET 2 G26 G25 G26")
                  wrt("LINK 0 G25")
                }
              }
            }
          }
        }
      } else {
        nonfp.s <- nonfp
        for (i in 1:length(resp.names)) {
          if (D["mode"] == 0) {
            nonfp.s <- gsub(paste(".", resp.names[i], sep = ""), "", nonfp.s)
          }
          if (D["mode"] == 1) {
            nonfp.s <- gsub(paste(".(>=", resp.names[i], ")", sep = ""), "", nonfp.s)
          }
        }
        nonfp.s <- unique(nonfp.s)
        expla <- expl
        for (p in expla) {
          if (is.na(nonfp[1]) || sum(p == nonfp.s) == 0) {
            if (is.null(categ) || sum(p == categ["var", ]) == 0) {
              for (j in 1:nresp) {
                wrt("LINK 1 G25")
                wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", p), "_", resp.names[j]), "'"))
                wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", p), "_", resp.names[j], "'"))
                wrt("GSET 2 G26 G25 G26")
                wrt("LINK 0 G25")
              }
              
            } else {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                categ.names <- levels(indata[[p]])
                for (j in 1:nresp) {
                  for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                    wrt("LINK 1 G25")
                    wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                    wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                    wrt("GSET 2 G26 G25 G26")
                    wrt("LINK 0 G25")
                  }
                }
              } else {
                categ.names <- levels(indata[[p]])
                refx <- categ["ref", which(p == categ["var", ])]
                categ.names <- categ.names[-which(refx == categ.names)]
                for (j in 1:nresp) {
                  for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                    wrt("LINK 1 G25")
                    wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                    wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                    wrt("GSET 2 G26 G25 G26")
                    wrt("LINK 0 G25")
                  }
                }
              }
            }
          }
        }
      }
    } else {
      if (D[1] == "Multivariate Normal" || D[1] == "Mixed") {
        nresp <- length(resp)
        
        if (is.list(expl)) {
          nonfp.s <- nonfp.sep
          for (i in 1:length(resp)) {
            nonfp.s <- gsub(paste(".", resp[i], sep = ""), "", nonfp.s)
          }
          nonfp.s <- unique(nonfp.s)
          if (is.na(sep.coeff[1]))
            sep.coeff <- character(0)
          for (p in sep.coeff) {
            if (is.na(nonfp.sep[1]) || sum(p == nonfp.s) == 0) {
              if (is.null(categ) || sum(p == categ["var", ]) == 0) {
                for (j in 1:nresp) {
                  wrt("LINK 1 G25")
                  wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", p), "_", resp[j]), "'"))
                  wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", p), "_", resp[j], "'"))
                  wrt("GSET 2 G26 G25 G26")
                  wrt("LINK 0 G25")
                }
              } else {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  categ.names <- levels(indata[[p]])
                  for (j in 1:nresp) {
                    for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                      wrt("LINK 1 G25")
                      wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                      wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                      wrt("GSET 2 G26 G25 G26")
                      wrt("LINK 0 G25")
                    }
                  }
                } else {
                  categ.names <- levels(indata[[p]])
                  refx <- categ["ref", which(p == categ["var", ])]
                  categ.names <- categ.names[-which(refx == categ.names)]
                  for (j in 1:nresp) {
                    for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                      wrt("LINK 1 G25")
                      wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                      wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                      wrt("GSET 2 G26 G25 G26")
                      wrt("LINK 0 G25")
                    }
                  }
                }
              }
            }
          }
          kk <- 1
          for (p in common.coeff) {
            newp <- paste(p, paste(which(as.logical(common.coeff.id[kk, ])), collapse = ""), sep = ".")
            kk <- kk + 1
            nonfp.c <- nonfp.common
            if (is.na(nonfp.common[1]) || sum(newp == nonfp.c) == 0) {
              if (is.null(categ) || sum(p == categ["var", ]) == 0) {
                wrt("LINK 1 G25")
                wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", newp)), "'"))
                wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", newp), "'"))
                wrt("GSET 2 G26 G25 G26")
                wrt("LINK 0 G25")
              } else {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  categ.names <- levels(indata[[p]])
                  for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                    wrt("LINK 1 G25")
                    wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i])), "'"))
                    wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "'"))
                    wrt("GSET 2 G26 G25 G26")
                    wrt("LINK 0 G25")
                  }
                } else {
                  categ.names <- levels(indata[[p]])
                  refx <- categ["ref", which(p == categ["var", ])]
                  categ.names <- categ.names[-which(refx == categ.names)]
                  for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                    wrt("LINK 1 G25")
                    wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i])), "'"))
                    wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "'"))
                    wrt("GSET 2 G26 G25 G26")
                    wrt("LINK 0 G25")
                  }
                }
              }
            }
          }
        } else {
          nonfp.s <- nonfp
          for (i in 1:length(resp)) {
            nonfp.s <- gsub(paste(".", resp[i], sep = ""), "", nonfp.s)
          }
          nonfp.s <- unique(nonfp.s)
          expla <- expl
          for (p in expla) {
            if (is.na(nonfp[1]) || sum(p == nonfp.s) == 0) {
              if (is.null(categ) || sum(p == categ["var", ]) == 0) {
                for (j in 1:nresp) {
                  wrt("LINK 1 G25")
                  wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", p), "_", resp[j]), "'"))
                  wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", p), "_", resp[j], "'"))
                  wrt("GSET 2 G26 G25 G26")
                  wrt("LINK 0 G25")
                }
              } else {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  categ.names <- levels(indata[[p]])
                  for (j in 1:nresp) {
                    for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                      wrt("LINK 1 G25")
                      wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                      wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                      wrt("GSET 2 G26 G25 G26")
                      wrt("LINK 0 G25")
                    }
                  }
                } else {
                  categ.names <- levels(indata[[p]])
                  refx <- categ["ref", which(p == categ["var", ])]
                  categ.names <- categ.names[-which(refx == categ.names)]
                  for (j in 1:nresp) {
                    for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                      wrt("LINK 1 G25")
                      wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j]), "'"))
                      wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "_", resp.names[j], "'"))
                      wrt("GSET 2 G26 G25 G26")
                      wrt("LINK 0 G25")
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        expla <- expl
        for (p in expla) {
          if (is.na(nonfp[1]) || sum(p == nonfp) == 0) {
            if (is.null(categ) || sum(p == categ["var", ]) == 0) {
              wrt("LINK 1 G25")
              wrt(paste0("NAME G25[1] '", shortname("FP_", p), "'"))
              wrt(paste0("DESC G25[1] 'FP:", p, "'"))
              wrt("GSET 2 G26 G25 G26")
              wrt("LINK 0 G25")
            } else {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                categ.names <- levels(indata[[p]])
                for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                  wrt("LINK 1 G25")
                  wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i])), "'"))
                  wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "'"))
                  wrt("GSET 2 G26 G25 G26")
                  wrt("LINK 0 G25")
                }
              } else {
                categ.names <- levels(indata[[p]])
                refx <- categ["ref", which(p == categ["var", ])]
                categ.names <- categ.names[-which(refx == categ.names)]
                for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                  wrt("LINK 1 G25")
                  wrt(paste0("NAME G25[1] '", shortname("FP_", chartr(".", "_", categ.names[i])), "'"))
                  wrt(paste0("DESC G25[1] 'FP:", chartr(".", "_", categ.names[i]), "'"))
                  wrt("GSET 2 G26 G25 G26")
                  wrt("LINK 0 G25")
                }
              }
            }
          }
        }
      }
    }
    
    wrt.resid <- function(rpx, resid.lev) {
      nrpx <- length(rpx)
      for (j in 1:nrpx) {
        for (i in 1:j) {
          if (i == j) {
            wrt("LINK 1 G25")
            wrt(paste0("NAME G25[1] '", shortname("RP", resid.lev, "_var_", chartr(".", "_", rpx[i])), "'"))
            wrt(paste0("DESC G25[1] 'RP", resid.lev, ":var(", chartr(".", "_", rpx[i]), ")'"))
            wrt("GSET 2 G26 G25 G26")
            wrt("LINK 0 G25")
          } else {
            wrt("LINK 1 G25")
            wrt(paste0("NAME G25[1] '", shortname("RP", resid.lev, "_cov_", chartr(".", "_", rpx[i]), "_", chartr(".", "_", rpx[j])), "'"))
            wrt(paste0("DESC G25[1] 'RP", resid.lev, ":cov(", chartr(".", "_", rpx[i]), ",", chartr(".", "_", rpx[j]), ")'"))
            wrt("GSET 2 G26 G25 G26")
            wrt("LINK 0 G25")
          }
        }
      }
    }
    
    wrt.resid2 <- function(rpx, resid.lev, clre) {
      nrpx <- length(rpx)
      nclre <- ncol(clre)
      k <- 1
      for (j in 1:nrpx) {
        for (i in 1:j) {
          if (i == j) {
            if (resid.lev == as.numeric(clre[1, k]) && rpx[i] == clre[2, k] && rpx[i] == clre[3, k]) {
              if (k < ncol(clre))
                k <- k + 1
            } else {
              wrt("LINK 1 G25")
              wrt(paste0("NAME G25[1] '", shortname("RP", resid.lev, "_var_", chartr(".", "_", rpx[i])), "'"))
              wrt(paste0("DESC G25[1] 'RP", resid.lev, ":var(", chartr(".", "_", rpx[i]), ")'"))
              wrt("GSET 2 G26 G25 G26")
              wrt("LINK 0 G25")
            }
          } else {
            if ((resid.lev == as.numeric(clre[1, k]) && rpx[i] == clre[2, k] && rpx[j] == clre[3, k]) || (resid.lev ==
                                                                                                            as.numeric(clre[1, k]) && rpx[j] == clre[2, k] && rpx[i] == clre[3, k])) {
              if (k < ncol(clre))
                k <- k + 1
            } else {
              wrt("LINK 1 G25")
              wrt(paste0("NAME G25[1] '", shortname("RP", resid.lev, "_cov_", chartr(".", "_", rpx[i]), "_", chartr(".", "_", rpx[j])), "'"))
              wrt(paste0("DESC G25[1] 'RP", resid.lev, ":cov(", chartr(".", "_", rpx[i]), ",", chartr(".", "_", rpx[j]), ")'"))
              wrt("GSET 2 G26 G25 G26")
              wrt("LINK 0 G25")
            }
          }
        }
      }
    }
    
    wrt.resid3 <- function(rpx, resid.lev) {
      nrpx <- length(rpx)
      for (j in 1:nrpx) {
        for (i in 1:j) {
          if (i == j) {
            wrt("LINK 1 G25")
            wrt(paste0("NAME G25[1] '", shortname("RP", resid.lev, "_var_", chartr(".", "_", rpx[i])), "'"))
            wrt(paste0("DESC G25[1] 'RP", resid.lev, ":var(", chartr(".", "_", rpx[i]), ")'"))
            wrt("GSET 2 G26 G25 G26")
            wrt("LINK 0 G25")
          }
        }
      }
    }
    
    if (nrp > 0) {
      for (ii in 1:nrp) {
        if (!is.null(fact)) {
          wrt.resid3(rp[[ii]], as.numeric(sub("rp", "", rp.names[ii])))
        } else {
          if (is.null(clre)) {
            wrt.resid(rp[[ii]], as.numeric(sub("rp", "", rp.names[ii])))
          } else {
            wrt.resid2(rp[[ii]], as.numeric(sub("rp", "", rp.names[ii])), clre)
          }
        }
      }
    }
    # Add in extra parameters ect.
    if (D[1] == "Multinomial" && as.numeric(D["mode"]) == 0) {
      wrt("LINK 1 G25")
      wrt(paste0("NAME G25[1] 'RP1_bcons_1'"))
      wrt(paste0("DESC G25[1] 'RP1:bcons_1'"))
      wrt("GSET 2 G26 G25 G26")
      wrt("LINK 0 G25")
    }
    
    wrt("GSIZ G26 b1000")
    wrt("LINK 3 G27")
    wrt(paste("CODE ", iterations/thinning, "b1000", 1, "G27[1]"))
    wrt(paste0("CALC G27[1] = G27[1] * ", thinning))
    wrt("NAME   G27[1] 'itnum'")
    wrt(paste("CODE b1000", 1, iterations/thinning, "G27[2]"))
    wrt("NAME   G27[2] 'parnum'")
    wrt("NAME   G27[3] 'iteration'")
    wrt("DESC   G27[3] '\\Iteration'")
    
    wrt("UNVE b1000 'parnum' 'itnum' 'mcmcchains' 'iteration' G26")
    if (D[1] == "Multinomial" && as.numeric(D["mode"]) == 1) {
      wrt("LINK 1 G25")
      wrt(paste("PUT ", iterations/thinning, 1, "G25[1]"))
      wrt("NAME G25[1] 'RP1_bcons_1'")
      wrt("DESC G25[1] 'RP1:bcons_1'")
      wrt("GSET 2 G26 G25 G26")
      wrt("LINK 0 G25")
    }
    wrt(paste0("PSTA '", chainfile, "' ", "'iteration' 'deviance' G26"))
    wrt("ERAS 'itnum' 'parnum' 'iteration' G26")
    wrt("LINK 0 G27")
    wrt("LINK 0 G26")
    
    calcresiduals <- function(level, displevel, rpx, resioptions, clre = clre) {
      wrt("")
      
      len.rpx <- length(rpx)
      # For MCMC there is always only one residual column at level 1
      if (level == 1) {
        len.rpx = 1
      }
      
      ii <- 1
      wrt(paste("LINK", len.rpx, "G26"))
      for (k in 1:len.rpx) {
        if (!is.null(clre)) {
          if (!(level == as.numeric(clre[1, ii]) && rpx[k] == clre[2, ii] && rpx[k] == clre[3, ii])) {
            wrt(paste0("NAME G26[", k, "] '", shortname("lev_", displevel, "_resi_est_", rpx[k]), "'"))
            wrt(paste0("DESC G26[", k, "] ", "'residual estimates'"))
          } else {
            if (ii < ncol(clre))
              ii <- ii + 1
          }
        } else {
          wrt(paste0("NAME G26[", k, "] '", shortname("lev_", displevel, "_resi_est_", rpx[k]), "'"))
          wrt(paste0("DESC G26[", k, "] ", "'residual estimates'"))
        }
      }
      
      wrt(paste("LINK", len.rpx, "G27"))
      if ("variance" %in% resioptions) {
        ii <- 1
        for (k in 1:len.rpx) {
          if (!is.null(clre)) {
            if (!(level == as.numeric(clre[1, ii]) && rpx[k] == clre[2, ii] && rpx[k] == clre[3, ii])) {
              wrt(paste0("NAME G27[", k, "] '", shortname("lev_", displevel, "_resi_variance_", rpx[k]), "'"))
              wrt(paste0("DESC G27[", k, "] ", "'residual variance'"))
            } else {
              if (ii < ncol(clre))
                ii <- ii + 1
            }
          } else {
            wrt(paste0("NAME G27[", k, "] '", shortname("lev_", displevel, "_resi_variance_", rpx[k]), "'"))
            wrt(paste0("DESC G27[", k, "] ", "'residual variance'"))
          }
        }
      } else {
        residual_se <- NULL
        ii <- 1
        for (k in 1:len.rpx) {
          if (!is.null(clre)) {
            if (!(level == as.numeric(clre[1, ii]) && rpx[k] == clre[2, ii] && rpx[k] == clre[3, ii])) {
              wrt(paste0("NAME G27[", k, "] '", shortname("lev_", displevel, "_resi_se_", rpx[k]), "'"))
              wrt(paste0("DESC G27[", k, "] ", "'residual standard error'"))
            } else {
              if (ii < ncol(clre))
                ii <- ii + 1
            }
          } else {
            wrt(paste0("NAME G27[", k, "] '", shortname("lev_", displevel, "_resi_se_", rpx[k]), "'"))
            wrt(paste0("DESC G27[", k, "] ", "'residual standard error'"))
          }
        }
      }
      wrt("RFUN")
      wrt("ROUT G26 G27")
      wrt("")
      
      wrt(paste("RLEV   ", level, sep = ""))
      wrt("RCOV   1")
      
      outgroups <- c("G26", "G27")
      
      if ("standardised" %in% resioptions) {
        wrt(paste("LINK", len.rpx, "G28"))
        ii <- 1
        for (k in 1:len.rpx) {
          if (!is.null(clre)) {
            if (!(level == as.numeric(clre[1, ii]) && rpx[k] == clre[2, ii] && rpx[k] == clre[3, ii])) {
              wrt(paste0("NAME G28[", k, "] '", shortname("lev_", displevel, "_std_resi_est_", rpx[k]), "'"))
              wrt(paste0("DESC G28[", k, "] ", "'std standardised residual'"))
            } else {
              if (ii < ncol(clre))
                ii <- ii + 1
            }
          } else {
            wrt(paste0("NAME G28[", k, "] '", shortname("lev_", displevel, "_std_resi_est_", rpx[k]), "'"))
            wrt(paste0("DESC G28[", k, "] ", "'std standardised residual'"))
          }
        }
        
        wrt("RTYP   0")
        wrt("MCRE")
        
        ccount <- 1
        for (k in 1:len.rpx) {
          wrt(paste0("CALC G28[", k, "]=G26[", k, "]/sqrt(G27[", k, "])", sep = ""))
        }
        outgroups <- c(outgroups, "G28")
      }
      
      wrt("RTYP   1")  # Compute comparative variances
      wrt("MCRE")
      
      if (!("variance" %in% resioptions)) {
        for (k in 1:len.rpx) {
          wrt(paste0("CALC G27[", k, "]=sqrt(G27[", k, "])"))  # Convert the variances to standard errors
        }
      }
      
      wrt("")
      wrt(paste0("NOBS ", level, " b30 b31"))
      wrt("LINK 1 G29")
      wrt("GENE 1 b30 1 G29[1]")
      wrt(paste0("NAME G29[1] 'lev_", displevel, "_residualid'"))
      outgroups <- c(outgroups, "G29")
      wrt("")
      wrt("LINK 0 G30")
      for (i in 1:length(outgroups)) {
        wrt(paste("GSET 2 G30", outgroups[i], "G30"))
        wrt(paste("LINK 0", outgroups[i]))
      }
    }
    
    if (resi.store && nrp > 0) {
      for (j in nrp:1) {
        rpx <- rp[[j]]
        len.rpx <- length(rp[[j]])
        wrt(paste("NOTE Calculate level ", as.numeric(sub("rp", "", rp.names[j])), " residuals", sep = ""))
        levtt <- as.numeric(sub("rp", "", rp.names[j]))
        displevel <- levtt
        if (D[1] == "Multivariate Normal" || D[1] == "Mixed" || D[1] == "Multinomial") {
          displevel <- displevel - 1
        }
        calcresiduals(levtt, displevel, rpx, resioptions, clre = clre)
        wrt(paste("PSTA '", resifile[j], "' G30", sep = ""))
      }
      wrt("ERAS G30")
      wrt("LINK 0 G30")
      
      if (!is.null(resi.store.levs)) {
        resiname <- rep(NA, length(resi.store.levs))
        for (i in 1:length(resi.store.levs)) {
          resiname[i] <- paste("'resi_lev", resi.store.levs[i], "'", sep = "")
        }
        wrt(paste("PSTA '", resichains, "' ", paste(resiname, collapse = " "), sep = ""))
      }
    }
  }
  if (debugmode) {
    wrt("WPMT 'Do you want to close MLwiN?' b50")
    wrt("SWITch b50")
    wrt("  CASE 1:")
    wrt("    EXIT ")
  } else {
    wrt("EXIT")
  }
  return(namemap)
}
