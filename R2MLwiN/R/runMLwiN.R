#' Calls MLwiN from R.
#'
#' This function executes MLwiN and then brings results back to R.
#'
#' @param Formula A \code{\link[stats]{formula}} object specifying the model
#' formula. See \code{\link{Formula.translate}} (\code{\link{Formula.translate.compat}}
#' details back-compatible functionality for deprecated syntax used in
#' versions of \pkg{R2MLwiN} prior to 0.8-0) and also `Details' below.
#' @param levID A character vector specifying the level ID(s). Deprecated
#' syntax: by default this is \code{NULL} and level ID(s) are specified
#' in the \code{Formula} object.
#' @param D A character string/vector specifying the type of distribution to be modelled, which
#' can include \code{'Normal'} (the default), \code{'Binomial'}, \code{'Poisson'},
#' \code{'Negbinom'}, \code{'Unordered Multinomial'}, \code{'Ordered Multinomial'},
#' \code{'Multivariate Normal'}, or \code{'Mixed'}. In the case of the latter,
#' \code{'Mixed'} precedes the response types which also need to be listed in
#' \code{D}, e.g. \code{c('Mixed', 'Normal', 'Binomial')}; these need to be
#' be listed in the same order to which they are referred to in the
#' \code{Formula} object (see \code{\link{Formula.translate}},
#' \code{\link{Formula.translate.compat}}). For (R)IGLS estimation (i.e. \code{EstM = 0}
#' in \code{estoptions}) \code{'Mixed'} combinations can consist of
#' \code{'Normal'} and \code{'Binomial'} or \code{'Normal'} and \code{'Poisson'};
#' for MCMC estimation (i.e. \code{EstM = 0}), on the other hand, only a combination
#' of \code{'Normal'} and \code{'Binomial'} is available.
#'
#' @param data A data.frame object containing the data to be modelled.
#' Optional (but recommended): if empty, data taken from environment of
#' \code{formula}.
#' @param estoptions A list of options used for estimating the model. See
#' `Details' below.
#' @param BUGO A vector specifying BUGS options. If non-null, then
#' WinBUGS/OpenBUGS, in conjunction with MLwiN, are used for modelling. Non-null
#' only applicable if \code{EstM = 1}. See `Details', below.
#' @param MLwiNPath A path to the MLwiN folder. By default, \code{MLwiNPath = NULL}
#' and path set by \code{options('MLwiN_path')}, the default for which can be
#' changed via \code{options(MLwiN_path = 'path/to/MLwiN vX.XX/')}).
#' @param stdout See \code{\link[base]{system2}}; \code{''} by default (i.e.
#' output to \code{stdout} sent to R console).
#' @param stderr See \code{\link[base]{system2}}; \code{''} by default (i.e.
#' output to \code{stderr} sent to R console).
#' @param workdir A path to the folder where the outputted files are to be saved.
#' If the folder specified does not exist, a new folder of that name is
#' created; \code{workdir = tempdir()} by default.
#' @param checkversion If \code{TRUE} (default), returns version number unless
#' (a) version detected is unknown or newer than MLwiN version available
#' when current version of R2MLwiN was released, in which case returns text
#' to this effect, or (b) version detected > 1 year older than MLwiN version
#' available when current version of R2MLwiN was released, in which case
#' function call stopped and user invited to update via usual channels. Can
#' disable via \code{FALSE} e.g. if slowing execution time down (for example
#' in a simulation).
#' @param indata A \code{data.frame} object containing the data to be modelled.
#' Deprecated syntax: by default this is \code{NULL} and the \code{data.frame}
#' is instead referenced via \code{data}.
#' @param saveworksheet A file name (or list of file names if more than one chain
#' is specified) used to store the MLwiN worksheet after the model has been estimated.
#'
#' @details
#' With regard to \code{runMLwiN}'s \code{Formula} object, see \code{\link[stats]{formula}}
#' for notes on general usage, noting the following differences:
#'
#' \itemize{
#' \item{The intercept is not included by default (this is keeping with the manner
#' in which models are specified in MLwiN). To include an intercept, then, one
#' can specify e.g. \code{normexam ~ 1 + standlrt + (1 | student)} or, assuming \code{cons}
#' is a constant of ones, \code{normexam ~ cons + standlrt + (cons | student)}. (Note also,
#' as further detailed below, for normal response models the level 1 ID (\code{student} in this example)
#' needs to be explicitly included in the random part of the model formula; this is not the
#' case for discrete response models.}
#' \item{The link function and denominator are included in the \code{Formula} object, e.g.
#' fitting a logistic model in which the variable \code{denom} is specified as the denominator:
#' \code{logit(resp, denom) ~ 1 + age + (1 | region)}.}
#' }
#'
#' Further details are as follows.
#'
#' The random part of the model is specified in sets of parentheses arranged in
#' descending order with respect to their hierarchy. E.g. in the case of a 3-level
#' model, the variable containing the level 3 ID is specified first, then
#' the variable containing the level 2 ID, etc. Note that the variable containing
#' the level 1 ID also needs to be explicitly specified unless
#' it is a discrete response model (in which case you should not specify it).
#'
#' The table below summarises the options for the \code{Formula} argument in
#' \pkg{R2MLwiN}. They assume an intercept is added (via \code{~ 1}; for alternative
#' specifications see \code{\link[stats]{formula}}). \code{<link>} denotes the link function,
#' \code{<y1>}, \code{<y2>}, etc. represent response variables, \code{<denom>} denotes
#' the denominator, \code{<offs>} the offset (optional), \code{<L2>}, \code{<L1>}, etc. the
#' variables containing the level 2 and level 1 identifying codes, and \code{<ref_cat>}
#' represents the reference category of a categorical response variable (optional:
#' if unspecified the lowest level of the factor is used as the reference category).
#' Explanatory variables are specified as e.g. \code{<x1> + <x2>}. For \code{'Ordered Multinomial'},
#' \code{'Multivariate Normal'} and \code{'Mixed'} responses, \code{[<common>]} indicates
#' a common coefficient (i.e. the same for each category) is to be fitted; here \code{<common>}
#' takes the form of a numeric identifier indicating the responses for which a common
#' coefficient is to be added (e.g. \code{[1:5]} to fit a common coefficient for
#' categories \code{1} to \code{5} of a 6-point ordered variable, \code{[1]} to fit a common
#' coefficient for the response variable specified first in the \code{Formula} object
#' for a \code{'Mixed'} response model, etc.) Otherwise a separate coefficient
#' (i.e. one for each category) is added. For \code{'Mixed'} response models, the
#' \code{Formula} arguments need to be grouped in the order the distributions
#' are listed in \code{D}.
#'
#' * denotes IGLS only in the table below.
#'
#' \tabular{lll}{
#' \strong{Distribution} \tab \strong{Format of \code{Formula} object} \tab \strong{Where \code{<link>} can equal...}\cr
#' \code{'Normal'} \tab \code{<y1> ~ 1 + <x1> + (1|<L2>) + (1|<L1>) + ...} \tab (identity link assumed)\cr
#' \code{'Poisson'} \tab \code{<link>(<y1>) ~ 1 + offset(<offs>) + <x1> + (1|<L2>) + ...} \tab \code{log}\cr
#' \code{'Negbinom'}* \tab \code{<link>(<y1>) ~ 1 + offset(<offs>) + (1|<L2>) + ...} \tab \code{log}\cr
#' \code{'Binomial'} \tab \code{<link>(<y1>, <denom>) ~ 1 + <x1> + (1|<L2>) + ...} \tab \code{logit},\code{probit},\code{cloglog}\cr
#' \code{'Unordered Multinomial'} \tab \code{<link>(<y1>, <denom>, <ref_cat>) ~ 1 + <x1> + (1|<L2>) + ...} \tab \code{logit}\cr
#' \code{'Ordered Multinomial'} \tab \code{<link>(<y1>, <denom>, <ref_cat>) ~ 1 + <x1> + <x2>[<common>] + (1[<common>]|<L3>) + (1|<L2>) + ...} \tab \code{logit},\code{probit},\code{cloglog}\cr
#' \code{'Multivariate Normal'} \tab \code{c(<y1>, <y2>, ...) ~ 1 + <x1> + <x2>[<common>] + (1[<common>]|<L3>) + (1|<L2>) + (1|<L1>) + ...} \tab (identity link assumed)\cr
#' \code{c('Mixed', 'Normal', 'Binomial')} \tab \code{c(<y1>, ..., <link> (<y2>, <denom>), ...) ~ 1 + <x1> + <x2>[<common>] + (1[<common>]|<L3>) + (1|<L2>) + (1|<L1>) + ...} \tab \code{logit}*,\code{probit},\code{cloglog}*\cr
#' \code{c('Mixed', 'Normal', 'Poisson')}* \tab \code{c(<y1>, ..., <link>(<y2>, <offset>), ...) ~ 1 + <x1> + <x2>[<common>] + (1[<common>]|<L3>) + (1|<L2>) + (1|<L1>) + ...} \tab \code{log}\cr
#' }
#'
#' The argument \code{estoptions} is a list which can contain the
#' following options used for estimating the model:
#'
#' \itemize{
#' \item \code{EstM}: specifies estimation method. When \code{EstM = 0} (default), estimation
#' method is (R)IGLS, otherwise \code{EstM = 1} specifies MCMC estimation.
#'
#' \item \code{resi.store}: a logical value indicating whether residuals are to be
#' stored or not. Defaults to \code{FALSE}.
#'
#' \item \code{resioptions}: a string vector to specify the various residual options.
#' The \code{'variance'} option calculates the posterior variances instead of
#' the posterior standard errors; the \code{'standardised'}, \code{'leverage'}, \code{'influence'}
#' and \code{'deletion'} options calculate standardised,
#' leverage, influence and deletion residuals respectively; the
#' \code{'sampling'} option calculates the sampling variance covariance matrix
#' for the residuals; the \code{'norecode'} option prevents residuals with values exceedingly close or
#' equal to zero from being recoded to missing.
#' When \code{EstM = 1} (i.e. MCMC estimation) \code{'variance'}
#' is default value, and the only other permissible value is \code{'standardised'}
#' (else function call stopped with appropriate error message).
#' When \code{EstM = 0} (i.e. (R)IGLS estimation), \code{'variance'}
#' cannot be specified together with \code{'standardised'}, \code{'leverage'} or
#' \code{'deletion'} (function call stopped with appropriate error message).
#' Default is \code{resioptions = c('variance')}.
#'
#' \item \code{resi.store.levs}: an integer vector indicating the levels at which the
#' residual chains are to be stored (\code{NULL} by default). Non-\code{NULL} values
#' not valid when \code{EstM = 0} (i.e. (R)IGLS estimation), else if \code{EstM = 0}
#' and \code{resi.store.levs} non-\code{NULL}, residual chains at specified levels
#' are returned.
#'
#' \item \code{debugmode}: a logical value determining whether MLwiN is run in the
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
#'
#' \item \code{x64}: a logical value indicating
#' whether the 64 bit version of MLwiN is used (unless \code{MLwiNPath} or \code{options(MLwiNPath)}
#' has been set directly to the executable). The default is determined by the characteristics
#' of the operating system on which the script is executed. If \code{FALSE},
#' the 32 bit version is called, if \code{TRUE} 64 bit version is called.
#'
#' \item \code{clean.files}: specifies whether the generated files are removed from
#' the \code{workdir} (\code{TRUE}, the default) or not (\code{FALSE}).
#'
#' \item \code{show.file}: a logical value indicating whether the output files (e.g.
#' MLwiN macro file) are shown on the screen. Defaults to \code{FALSE}.
#'
#' \item \code{clre}: a matrix used to define which elements of the random effects matrix
#' to remove (i.e. hold constant at zero). Removes
#' from the random part at level <first row> the covariance matrix element(s)
#' defined by the pair(s) of rows <second row> <third row>. Each column
#' corresponds to a removed entry of the covariance matrix. See e.g. \code{demo(UserGuide07)}
#' for an example.
#'
#' \item \code{notation}: specifies the model subscript notation
#' to be used in the MLwiN equations window. \code{'class'} means no multiple
#' subscripts, whereas \code{'level'} has multiple subscripts. If
#' \code{notation = NULL}, defaults to \code{'level'} if \code{'xc = NULL'} else
#' defaults to \code{'class'}.
#'
#' \item \code{mem.init}: sets and displays worksheet capacities for
#' the current MLwiN session. A vector of length 5 corresponding to
#' the following order: number of levels (defaults to 1 + the number of
#' levels specified in the function call); worksheet size in thousands of cells
#' (default is 6000); the number of columns (default is 2500); the number of
#' explanatory variables (default it 10 + number of explanatory variables
#' calculated initially); the number of group labels (default is 20).
#'
#' \item \code{optimat}: instructs MLwiN to limit the maximum matrix size
#' that can be allocated by the (R)IGLS algorithm. Specify \code{optimat = TRUE}
#' if MLwiN gives the following error message 'Overflow allocating smatrix'.
#' This error message arises if one or more higher-level units is/are extremely
#' large (containing more than 800 lower-level units). In this situation \code{runMLwiN}'s
#' default behaviour is to instruct MLwiN to allocate a larger matrix size to
#' the (R)IGLS algorithm than is currently possible. Specifying
#' \code{optimat = TRUE} caps the maximum matrix size at 800 lower-level units,
#' circumventing the MLwiN error message, and allowing most MLwiN
#' functionality.
#'
#' \item \code{nonlinear}: a character vector specifying linearisation method for discrete
#' response models estimated via IGLS (see Chapter 9 of Rasbash et al 2012,
#' and Goldstein 2011). \code{N = 0} specifies marginal quasi-likelihood
#' linearization (MQL), whilst \code{N = 1} specifies penalised quasi-
#' likelihood linearization (PQL); \code{M = 1} specifies first order
#' approximation, whilst \code{M = 2} specifies second order approximation.
#' \code{nonlinear = c(N = 0, M = 1)} by default. First order marginal
#' quasi-likelihood (MQL1) only option for single-level discrete response
#' models. Pertains to discrete response models estimated via IGLS: i.e. when
#' \code{EstM = 0} in \code{estoptions}, and for starting values when estimated via IGLS
#' for MCMC (\code{EstM = 1}).
#'
#' \item \code{Meth}: specifies which maximum likelihood estimation method is to be
#' used. If \code{Meth = 0} estimation method is set to RIGLS. If \code{Meth = 1}
#' estimation method is set to IGLS (the default setting).  Pertains to models
#' estimated via (R)IGLS: i.e. when \code{EstM = 0} in \code{estoptions}, and for starting
#' values when estimated via (R)IGLS for MCMC (\code{EstM = 1}).
#'
#' \item \code{merr}: a vector which sets-up measurement errors on predictor
#' variables. The first element \code{N} defines the number of variables that
#' have measurement errors. Then, for each variable with measurement error, a
#' pair of inputs are required: the first of these is the explanatory variable
#' name as a character string, and the second is the variance of
#' the measurement error for this variable. See \code{demo(MCMCGuide14)} for an
#' example.
#'
#' \item \code{fact}: a list of objects specified for factor analysis,
#' including:
#' \itemize{
#' \item \code{nfact}: Specifies the number of factors
#' \item \code{lev.fact}: Specifies the level/classification for the random part of
#' the factor for each factor.
#' \item \code{nfactcor}: Specifies the number of
#' correlated factors
#' \item \code{factcor}: A vector specifying the correlated
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
#' \item \code{weighting}: a deprecated option for specifying weights in IGLS estimation:
#' see \code{fpsandwich} and \code{rpsandwich} for new method of doing so.
#' \code{weighting} is a list of objects including \code{levels}, \code{weights},
#' \code{mode}, \code{FSDE} and \code{RSDE}; see \code{\link{write.IGLS}} for details.
#'
#' \item \code{centring}: deprecated method (only applicable when using old syntax
#' pre-\pkg{R2MLwiN} v.0.8-0) specifying function by
#' which explanatory variables are to be centred (users can instead transform
#' variables prior to \code{runMLwiN} call).
#' If non-\code{NULL}, centring is used for the selected explanatory
#' variables (\code{centring = NULL} by default). \code{centring} is a list of
#' objects specifying the methods to be used to centre specific explanatory
#' variables. E.g. \code{list(age = 1, ...)} specifies that the explanatory
#' variable \code{age} is to be centred around its grand mean;
#' \code{list(age = c(2, 'district'), ...)} specifies that \code{age} is to be
#' centred around its group mean, where group defined by the variable \code{district};
#' and \code{list(age = c(3, 18), ...)} specifies that \code{age} is to
#' be centred around the value \code{18}.
#'
#' \item \code{xclass}: a deprecated option for specifying cross-classified and/or
#' multiple membership models; see \code{xc} and \code{mm} for new method of
#' doing so. \code{xclass} is a list of objects including \code{class},
#' \code{N1}, \code{weight}, \code{id} and \code{car}; see \code{\link{write.MCMC}} for details.
#'
#' \item \code{mcmcOptions}: a list of objects specifying MCMC options, including the
#' following:
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
#' \item \code{drop.data}: If \code{TRUE} (default) only the data involved in the model
#' is passed to MLwiN, otherwise the entire dataset in \code{data} is passed.
#'
#' \item \code{drop.levels}: If \code{TRUE} (default) any unused levels are dropped from factors, otherwise the dataset
#' is left unchanged.
#'
#' \item \code{fpsandwich}: specifies standard error type for fixed parameters. If
#' \code{fpsandwich = TRUE}, robust or `sandwich' standard errors based on raw
#' residuals are used, if \code{fpsandwich = FALSE} (default) then standard,
#' uncorrected, IGLS or RIGLS computation used.
#'
#' \item \code{rpsandwich}: specifies standard error type for random parameters. If
#' \code{rpsandwich = TRUE}, robust or `sandwich' standard errors based on raw
#' residuals are used, if \code{rpsandwich = FALSE} (default) then standard,
#' uncorrected, IGLS or RIGLS `plug in' estimates used.
#'
#' \item \code{smat}: a matrix with two columns the levels at which a diagonal
#' matrix is to be specified. The first column specifies the level.
#' If the value of the second column is \code{1} then the random covariance matrix is
#' set to be diagonal.
#'
#' \item \code{maxiter}: a numeric value specifying the maximum number of iterations, from
#' the start, before (R)IGLS estimation halts. Pertains to models
#' estimated via (R)IGLS: i.e. when \code{EstM = 0} in \code{estoptions}, and for starting
#' values when estimated via (R)IGLS for MCMC (\code{EstM = 1}).
#'
#' \item \code{tol}: a numeric value specifying the convergence criterion.
#' If value is m, estimation will be
#' deemed to have converged when the relative change in the estimate for all
#' parameters from one iteration to the next is less than 10(-m). Defaults to
#' value of \code{2} for m if not otherwise specified.  Pertains to models
#' estimated via (R)IGLS: i.e. when \code{EstM = 0} in \code{estoptions}, and for starting
#' values when estimated via (R)IGLS for MCMC (\code{EstM = 1}).
#'
#' \item \code{extra}: if \code{TRUE}, extra binomial, extra negative binomial,
#' extra Poisson or extra multinomial distributions assumed, else \code{FALSE}.
#' can only be specified for discrete response models (i.e. \code{'Binomial'},
#' \code{'Negbinom'}, \code{'Poisson'}, \code{'Multinomial'})
#' estimated via (R)IGLS (i.e. \code{EstM = 0}).
#'
#' \item \code{reset}: a vector specifying the action to be
#' taken, at each level, if a variance parameter is estimated at a particular
#' iteration to be negative during estimation. Values specified in
#' ascending order of level hierarchy: if \code{0} a negative variance
#' estimate is reset to zero and so are any associated covariances; if \code{1}
#' a negative variance estimate is reset to zero but not the associated
#' covariances; if \code{2} no resetting takes place. E.g. \code{reset = c(0, 1)}
#' to assign value \code{0} to level 1 and value \code{1} to level 2 of
#' two-level model.
#'
#' \item \code{constraints}: \code{fixed.ui} and \code{fixed.ci} are used
#' to specify constraints on the fixed coefficients, and \code{random.ui}
#' and \code{random.ci} to specify constraints on the random parameters. The
#' syntax for specifying just fixed parameter constraints is
#' \code{constraints = list(fixed.ui = <fixed matrix>, fixed.ci = <fixed values>)},
#' where \code{<fixed matrix>} is a matrix where each row represents one fixed part
#' parameter, in the same order that they appear in the results table, each
#' column represents one constraint, and the values in the matrix are multipliers
#' for the parameters; and \code{<fixed values>} is a vector of values, one per
#' constraint, to which the parameters multiplied by the multipliers in the
#' corresponding column of \code{<fixed matrix>} should be equal. For example,
#' if we have a model with formula \code{y ~ 1 + x1 + x2 + x3 + x4 + (1|lev1ID)},
#' then \code{constraints = list(fixed.ui = matrix(c(0, 1, -1, 0, 0, 0, 0, 0, 1, 2), nrow = 5),
#' fixed.ci = c(0, 2))} specifies the constraints that the coefficient of \code{x1}
#' equals the coefficient of \code{x2} and that the coefficient of \code{x3} plus
#' twice the coefficient of \code{x4} equals \code{2}. Random constraints are
#' specified similarly, and fixed and random constraints may be applied
#' simultaneously. Applies to \code{EstM = 0} (i.e. estimation via (R)IGLS) only.
#'
#' \item \code{xc}: indicates whether model is cross-classified (\code{TRUE}) or
#' nested (\code{FALSE}). Ignored if \code{EstM = 0}, i.e. only applicable to
#' models estimated via MCMC. Defaults to \code{xc = FALSE}, unless either
#' \code{mm} or \code{car} are non-\code{NULL}, in which case \code{xc = TRUE}. Supersedes
#' deprecated \code{xclass}.
#'
#' \item \code{mm}: specifies the structure of a multiple membership model.
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
#' \code{logearn ~ 1 + age_40 + sex + parttime + (1 | company) + (1 | id)}, if
#' \code{company} is a multiple membership classification with the variables
#' indicating the classifications in \code{company}, \code{company2},
#' \code{company3}, \code{company4} and their weights in \code{weight1}, \code{weight2},
#' \code{weight3} and \code{weight4} then
#' \code{mm = list(list(mmvar = list('company', 'company2', 'company3', 'company4'),}
#' \code{weights = list('weight1', 'weight2', 'weight3', 'weight4')), NA)}
#' with the \code{NA}, listed last, corresponding to the level 1 identifier (\code{id}).
#'
#' \item \code{car}: specifies the structure of a conditional autoregressive (CAR)
#' model. Can be a list of variable names, a list of vectors, or a matrix (e.g. see
#' \code{\link{df2matrix}}). In the case of the former, each element of the list
#' corresponds to a level (classification) of
#' the model, in descending order. If a level is not a spatial classification,
#' then \code{NA} is specified. Otherwise, lists need to be assigned to
#' \code{carvar} and \code{weights}, with the former containing columns
#' specifying the spatial classification units, and the latter containing
#' columns specifying the weights. See \code{demo(MCMCGuide17)} for examples.
#' Ignored if \code{EstM = 0}, i.e. only applicable
#' to models estimated via MCMC. \code{car = NULL} by default. Supersedes
#' deprecated \code{xclass}. See \code{demo(MCMCGuide17)} for examples.
#'
#' \item \code{carcentre}: if CAR model (i.e. if \code{car} is non-\code{NULL}),
#' \code{carcentre = TRUE} mean-centres all random effects at that level.
#' \item \code{startval}: a list of numeric vectors specifying the starting values.
#' \code{FP.b} corresponds to the estimates for the fixed
#' part; \code{FP.v} specifies the variance/covariance estimates for the fixed
#' part; \code{RP.b} specifies the variance estimates for the random part;
#' \code{RP.v} corresponds to the variance/covariance matrix of the variance
#' estimates for the random part. \code{startval = NULL} by default: i.e. when
#' \code{EstM = 0} the OLS estimates are used, else if \code{EstM = 1} the
#' estimates obtained from IGLS are used as the starting values for MCMC.
#'
#' \item \code{sort.force}: If \code{TRUE} will sort data based on hierarchy as
#' determined by model formula; defaults to \code{FALSE}.
#'
#' \item \code{sort.ignore}: If \code{FALSE} will check data is sorted in a manner in
#' keeping with the hierarchy implied by the model formula, and will return a warning
#' if that is not the case.
#'
#' \item \code{mcmcMeth}: list of objects specifying MCMC methodology and prior
#' options, including the following (see \code{\link{write.MCMC}} for further details):
#' \itemize{
#' \item \code{iterations}: Number of main iterations post-burnin (i.e. monitoring chain length), defaults to 5000.
#' \item \code{burnin}: Length of burnin, defaults to 500.
#' \item \code{nchains}: Number of MCMC chains to run, defaults to 1.
#' \item \code{thinning}: Thinning factor, defaults to 1.
#' \item \code{seed}: MCMC random number seed, defaults to 1.
#' \item \code{priorParam}: A list specifying informative priors. This includes:
#' \code{fixe} -- for the fixed
#' parameters, if proper normal priors are used for some parameters, a list of
#' vectors of length two is provided, each of which specifies the mean and the
#' standard deviation. If not given, default ('flat' or 'diffuse') priors are
#' used for the parameters; \code{fixe.common} -- for multivariate normal,
#' multinomial and mixed response models, if common coefficients are added, use
#' \code{fixe.common} rather than \code{fixe}; \code{fixe.sep} -- if the common
#' coefficients are added, use \code{fixe.sep} for the separate coefficients;
#' \code{rp1} -- a list object specifying the Wishart or gamma prior for the
#' covariance matrix or scalar variance at level 1 (this consists of: (1)
#' \code{estimate} -- an estimate for the true value of the inverse of the
#' covariance matrix; (2) \code{size} -- the number of rows in the covariance
#' matrix. Note that this is a weakly-informative prior and the default prior
#' is used if missing); \code{rp2} -- a list object specifying the Wishart or
#' gamma prior for the covariance matrix or scalar variance at level 2 (this
#' consists of: (1) \code{estimate} -- an estimate for the true value of the
#' inverse of the covariance matrix; (2) \code{size} -- the number of rows in
#' the covariance matrix. Note that this is a weakly-informative prior and the
#' default prior is used if missing).
#' \item \code{scale}: Scale factor for proposal variances: this number will be
#' multiplied by the estimated parameter variance (from IGLS/RIGLS) to give the
#' proposal distribution variance. Defaults to 5.8.
#' \item \code{refresh}: Number of iterations after which screen (in MLwiN GUI) is
#' to be refreshed. Defaults to 50.
#' \item \code{fixM}: Specifies the estimation method for the fixed effects:
#' \code{1} for Gibbs sampling, \code{2} for univariate Metropolis-Hastings (MH)
#' sampling and \code{3} for multivariate MH sampling. Defaults to \code{2} if
#' Poisson, Multinomial, Binomial or Mixed model, else defaults to \code{1}.
#' \item \code{residM}: Specifies the estimation method for the random effects
#' (residuals): \code{1} for Gibbs sampling, \code{2} for univariate
#' Metropolis-Hastings (MH) sampling and \code{3} for multivariate MH sampling.
#' Defaults to \code{2} if Poisson, Multinomial, Binomial or Mixed model,
#' else defaults to \code{1}.
#' \item \code{Lev1VarM}: Specifies the estimation method for the level 1 variance:
#' \code{1} for Gibbs sampling, \code{2} for univariate
#' Metropolis-Hastings (MH) sampling and \code{3} for multivariate MH sampling.
#' Defaults to \code{2} if Poisson, Multinomial, Binomial or Mixed model,
#' else defaults to \code{1}.
#' \item \code{OtherVarM}: Specifies the estimation method for the higher level
#' variance matrices: \code{1} for Gibbs sampling, \code{2} for univariate
#' Metropolis-Hastings (MH) sampling and \code{3} for multivariate MH sampling.
#' Defaults to \code{1}.
#' \item \code{adaption}: \code{adaption = 1} (the default) indicates adaptation is to be used,
#' \code{adaption = 0} indicates it is not.
#' \item \code{tol}: An integer specifying tolerance (as a percentage; defaults to 10) when
#' \code{adaption = 1} (ignored if \code{adaption = 0}).
#' \item \code{rate}: An integer specifying the acceptance rate (as a percentage; defaults
#' to 50) when \code{adaption = 1} (ignored if \code{adaption = 0}).
#' \item \code{priorcode}: A vector indicating which default priors are to be used
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
#' \item \code{startval}: Deprecated: starting values are now specified directly
#' within \code{estoptions}.
#' \item \code{lclo}: Toggles on/off the possible forms of complex level
#' 1 variation when using MCMC. By default (\code{lclo = 0}) the level
#' 1 variation is expressed as a function of the predictors. Else
#' (\code{lclo = 1}) the log of the level 1 precision (1/variance) is expressed as
#' a function of the predictors. Defaults to \code{lclo = 0}.
#' \item \code{dami}: Outputs a complete (i.e. including non-missing
#' responses) response variable y. If \code{dami = c(0, <iter1>, <iter2>, ...)} then
#' the response variables returned will be the value of y at the iterations
#' quoted (as integers \code{<iter1>, <iter2>}, etc.); these can be used for
#' multiple imputation. If \code{dami = 1} the value of y will be the mean
#' estimate from the iterations produced. \code{dami = 2} is as for \code{dami = 1}
#' but with the standard errors of the estimate additionally being stored.
#' \code{dami = NULL} by default.
#' }
#' }
#' The argument \code{BUGO} is a vector specifying BUGS options as follows:
#' \itemize{
#' \item \code{n.chains}: specifies the
#' number of chains used by BUGS.
#' \item \code{debug}: determines
#' whether BUGS stays open following completion of the model run;
#' \code{debug = FALSE} by default.
#' \item \code{seed}: sets the random number
#' generator in BUGS.
#' \item \code{bugs}: specifies the path of the BUGS
#' executable.
#' \item \code{OpenBugs}: if \code{OpenBugs = TRUE}, OpenBUGS is used.
#' Otherwise (i.e. \code{OpenBugs = FALSE}, the default) WinBUGS is used.
#' }
#'
#' @return
#' If \code{BUGO} is non-NULL then the output is an \code{\link{mcmc.list}}
#' object.
#'
#' If the IGLS algorithm is used (i.e., \code{EstM = 0}), then returns \code{\link{mlwinfitIGLS-class}} object;
#' if MCMC estimation used (i.e., \code{EstM = 1}), then returns \code{\link{mlwinfitMCMC-class}} object.
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
#' @seealso
#' \code{\link[stats]{formula}}, \code{\link{Formula.translate}}, \code{\link{Formula.translate.compat}}, \code{\link{write.IGLS}}, \code{\link{write.MCMC}}
#'
#' @import doParallel foreach parallel
#' @importFrom stats acf as.formula cov density end getCall get_all_vars model.frame model.matrix model.offset na.omit pacf pnorm qnorm quantile sd start terms terms.formula update.formula var window complete.cases reshape
#' @importFrom grDevices dev.new
#' @importFrom graphics close.screen lines par plot points screen split.screen text
#' @importFrom utils read.delim stack
#' @importFrom methods is new validObject
#' @examples
#'
#' ## The R2MLwiN package includes scripts to replicate all the analyses in
#' ## Rasbash et al (2012) A User's Guide to MLwiN Version 2.26 and
#' ## Browne, W.J. (2012) MCMC estimation in MLwiN Version 2.26.
#' ## The MLwiN manuals are available online, see:
#' ## http://www.bristol.ac.uk/cmm/software/mlwin/download/manuals.html
#'
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' ## For a list of demo titles
#' demo(package = 'R2MLwiN')
#'
#' ## Take MCMCGuide03 as an example
#' ## To view file
#' file.show(system.file('demo', 'MCMCGuide03.R', package='R2MLwiN'))
#'
#' ## To run the demo
#' demo(MCMCGuide03)
#' }
#'
#' @export
runMLwiN <- function(Formula, levID = NULL, D = "Normal", data = NULL, estoptions = list(EstM = 0), BUGO = NULL, MLwiNPath = NULL,
                     stdout = "", stderr = "", workdir = tempdir(), checkversion = TRUE, indata = NULL, saveworksheet = NULL) {
  if (!is.null(indata) && !is.null(data)) {
    stop("Only one of data and indata can be specified")
  }
  if (!is.null(data)) {
    indata <- data
  }

  if (is.null(levID)) {
    oldsyntax <- FALSE
  } else {
    oldsyntax <- TRUE
    warning("This syntax has been superseded, see help for guidance on converting it.")
  }

  drop.data <- estoptions$drop.data
  if (is.null(drop.data)) {
    if (oldsyntax) {
      drop.data <- FALSE
    } else {
      drop.data <- TRUE
    }
  }

  drop.levels <- estoptions$drop.levels
  if (is.null(drop.levels)) {
    drop.levels <- TRUE
  }

  if (oldsyntax) {
    if (is.character(Formula)) {
      Formula <- gsub("\\{", "\\(", Formula)
      Formula <- gsub("\\}", "\\)", Formula)
      Formula <- gsub("[[:space:]]", "", Formula)
      cc <- c(0:length(levID))
      if (sum(grepl("\\({1}[[:digit:]]+\\|{2}", Formula)) > 0) {
        for (i in cc) {
          Formula <- sub(paste(i, "\\|{2}", sep = ""), paste("\\`", i, "c`\\|", sep = ""), Formula)
          Formula <- sub(paste(i, "\\|", sep = ""), paste("\\`", i, "s`\\|", sep = ""), Formula)
        }
      }
      if (sum(grepl("\\({1}[[:digit:]]+[[:alpha:]]{1}\\|", Formula)) > 0) {
        for (i in cc) {
          Formula <- sub(paste(i, "s\\|", sep = ""), paste("\\`", i, "s`\\|", sep = ""), Formula)
          Formula <- sub(paste(i, "c\\|", sep = ""), paste("\\`", i, "c`\\|", sep = ""), Formula)
        }
      }
      Formula <- as.formula(Formula)
    }
  } else {
    tmpvarnames <- unique(unlist(strsplit(all.vars(Formula), "\\.")))
    tForm <- as.formula(paste0("~", paste(tmpvarnames, collapse = "+")))
    if (drop.data) {
      indata <- get_all_vars(tForm, indata)
    } else {
      newdata <- get_all_vars(tForm, indata)
      newvars <- setdiff(colnames(newdata), colnames(indata))
      for (var in newvars) {
        indata[[var]] <- newdata[[var]]
      }
    }
  }

  if (drop.levels) {
    for (var in colnames(indata)) {
      if (is.factor(indata[[var]])) {
        if (length(setdiff(levels(indata[[var]]), levels(factor(indata[[var]])))) > 0) {
          indata[[var]] <- droplevels(indata[[var]])
          warning(paste0(var, " has unused factor levels defined. These were dropped from this model run, but we recommend removing them prior to calling runMLwiN."))
        }
      }
    }
  }

  EstM <- estoptions$EstM
  if (is.null(EstM))
    EstM <- 0
  if (EstM != 0 && EstM != 1) {
    stop("Invalid EstM option (can be zero or one)")
  }

  if (length(D) == 1) {
    if (!(D %in% c("Normal", "Binomial", "Poisson", "Negbinom", "Multivariate Normal", "Ordered Multinomial",
                   "Unordered Multinomial"))) {
      stop("Invalid distribution specified")
    }
  } else {
    if (D[1] != "Mixed") {
      stop("Invalid distribution specified")
    } else {
      for (i in 2:length(D)) {
        if (EstM == 0) {
          if (!(D[i] %in% c("Normal", "Binomial", "Poisson"))) {
            stop("Invalid distribution specified")
          }
        } else {
          if (!(D[i] %in% c("Normal", "Binomial"))) {
            stop("Invalid distribution specified")
          }
        }
      }
    }
  }

  # Check MLwiNPath is usable and set command/args
  debugmode <- estoptions$debugmode
  if (is.null(debugmode))
    debugmode <- FALSE

  x64 <- estoptions$x64
  if (is.null(x64)) {
    if (.Machine$sizeof.pointer == 8) {
      x64 <- TRUE
    } else {
      x64 <- FALSE
    }
  }

  if (is.null(MLwiNPath)) {
    MLwiNPath <- getOption("MLwiN_path")
  }

  pathinfo <- file.info(MLwiNPath)
  if (is.na(pathinfo$isdir)) {
    stop(paste0(MLwiNPath, " does not exist"))
  }

  if (!isTRUE(pathinfo$isdir)) {
    if (file.access(MLwiNPath, mode = 1) == 0) {
      cmd <- MLwiNPath
    } else {
      stop(paste0(MLwiNPath, " is not executable"))
    }
  }

  if (isTRUE(pathinfo$isdir)) {
    if (debugmode) {
      cmd <- paste0(MLwiNPath, "/i386/mlwin.exe")
      if (file.access(cmd, mode = 1) != 0) {
        cmd <- paste0(MLwiNPath, "/mlwin.exe")
      } else {
        if (file.access(cmd, mode = 1) != 0) {
          cmd <- paste0(MLwiNPath, "/i386/mlnscript.exe")
          if (file.access(cmd, mode = 1) != 0) {
            stop("Cannot find valid MLwiN executable")
          }
        }
      }
    } else {
      if (x64) {
        cmd <- paste0(MLwiNPath, "/x64/mlnscript.exe")
        if (file.access(cmd, mode = 1) != 0) {
          cmd <- paste0(MLwiNPath, "/i386/mlnscript.exe")
          if (file.access(cmd, mode = 1) != 0) {
            cmd <- paste0(MLwiNPath, "/mlwin.exe")
            if (file.access(cmd, mode = 1) != 0) {
              stop("Cannot find valid MLwiN executable")
            }
          }
        }
      } else {
        cmd <- paste0(MLwiNPath, "/i386/mlnscript.exe")
        if (file.access(cmd, mode = 1) != 0) {
          cmd <- paste0(MLwiNPath, "/mlwin.exe")
          if (file.access(cmd, mode = 1) != 0) {
            stop("Cannot find valid MLwiN executable")
          }
        }
      }
    }
  }


  versioninfostr <- '
version:date:md5:filename:x64:trial:platform
1.10:Jun 2001:44e840796f3c43113c45ac8fe7e0633a:mlwin.exe:FALSE:FALSE:win
1.10:Jul 2001:21e3a2d85e6f9c4bb8d926658981a020:mlwin.exe:FALSE:FALSE:win
2.00:May 2004:4ddb67c6426112bafc70bddca38cd63a:mlwin.exe:FALSE:FALSE:win
2.00:Jun 2004:bf9aff9fa66d8eadc8b3a170e616ab58:mlwin.exe:FALSE:FALSE:win
2.00:Jul 2004:f78c7d6d3a0dc82e7ea7e391a70ebb02:mlwin.exe:FALSE:FALSE:win
2.00:Nov 2004:e9fbafdc5715921dcadec601e3cec593:mlwin.exe:FALSE:FALSE:win
2.01:Dec 2004:359dbe8d728b2841b948504b5c272392:mlwin.exe:FALSE:FALSE:win
2.01:Feb 2004:f4a164b199b37b33158b194402034756:mlwin.exe:FALSE:FALSE:win
2.02:Jun 2005:cfb2ba2aea080ad69189709e87613ef7:mlwin.exe:FALSE:FALSE:win
2.10:Feb 2009:d96fdc4d9876206837d5c720bf37c8e1:mlwin.exe:FALSE:FALSE:win
2.11:Apr 2009:1da1348d7a65a3a7ae1f310d63520429:mlwin.exe:FALSE:FALSE:win
2.12:Jul 2009:7f44b98b0ca60ea6b34ee56f962869c7:mlwin.exe:FALSE:FALSE:win
2.13:Aug 2009:b435c8137676da09412ee6a57d7426cc:mlwin.exe:FALSE:FALSE:win
2.14:Sep 2009:2e493aa7cdf221caed82c0cdc4facb17:mlwin.exe:FALSE:FALSE:win
2.15:Oct 2009:3ca55fe4c04f546040fc4937f0ac1a9f:mlwin.exe:FALSE:FALSE:win
2.16:Nov 2009:afd80cecbe7e1164957f4530b07ea5ec:mlwin.exe:FALSE:FALSE:win
2.16:Nov 2009:a91958a92ee4e44bf80a58cb8c5a319a:mlwin.exe:FALSE:TRUE:win
2.17:Jan 2010:28e92f9aba0431d2a53ab4ea0c1471e6:mlwin.exe:FALSE:FALSE:win
2.18:Mar 2010:d0e7c52a33024ffe3bb176fac6fdd724:mlwin.exe:FALSE:FALSE:win
2.19:May 2010:46e40433d3f22f947ffc539bfddff58a:mlwin.exe:FALSE:FALSE:win
2.19:May 2010:2cd0159a1580452c43358511f763fb78:mlwin.exe:FALSE:TRUE:win
2.20:Jun 2010:cf9ba18ef770d1d5e761b99bd74cfb48:mlwin.exe:FALSE:FALSE:win
2.21:Oct 2010:71f36aecbbef624f70251d273547bae5:mlwin.exe:FALSE:FALSE:win
2.22:Dec 2010:d372e2ea4d3dd8202bddc0fc3e3be445:mlwin.exe:FALSE:FALSE:win
2.23:Apr 2011:0a150498818a6e519e1fb5f4c96863df:mlwin.exe:FALSE:FALSE:win
2.23:Apr 2011:dd6304ed39b0fd769e1e40e4b85e5b8f:mlwin.exe:FALSE:TRUE:win
2.24:Sep 2011:8b8a5d06d4440de87fa97359d06da8d6:mlwin.exe:FALSE:FALSE:win
2.24:Sep 2011:005a73f0c8af424520147151d504fffb:mlwin.exe:FALSE:TRUE:win
2.25:Feb 2012:a5a8e56a2da1faa75a0bf5a9c260f79b:mlwin.exe:FALSE:FALSE:win
2.25:Feb 2012:ce92f10b5146c3d4fd10c1cad78d01d3:mlwin.exe:FALSE:TRUE:win
2.26:Sep 2012:f915c285f8409fe66bf8ac0a90256fe7:mlwin.exe:FALSE:FALSE:win
2.26:Sep 2012:8bdcb5ef1b4a1c10b7a5f2e1c359fae4:mlwin.exe:FALSE:TRUE:win
2.26:Sep 2012:d4b3b6a97e0d413bf185debd18a7c388:mlnscript.exe:FALSE:FALSE:win
2.26:Sep 2012:6d06f90db77a9f4d3bf973cfd4be6aad:mlnscript.exe:TRUE:FALSE:win
2.27:Feb 2013:e25a7fb9431c024e2f09222434d9fc55:mlwin.exe:FALSE:FALSE:win
2.27:Feb 2013:d0330c49c3234474b0e7d79fcd83117d:mlwin.exe:FALSE:TRUE:win
2.27:Feb 2013:5bc7e8fade28bd8fb9f9ef110ec56afc:mlnscript.exe:FALSE:FALSE:win
2.27:Feb 2013:63fa77b06439f295231dd4795e4ed99e:mlnscript.exe:TRUE:FALSE:win
2.28:Jul 2013:6bdadad3615c49ca418cc63cb952d37f:mlwin.exe:FALSE:FALSE:win
2.28:Jul 2013:6a1f0d366ffa622e4a052695a6013e2c:mlwin.exe:FALSE:TRUE:win
2.28:Jul 2013:5fcbc7d0dd0a900c99ec411018ccdaa5:mlnscript.exe:FALSE:FALSE:win
2.28:Jul 2013:d299e1156e5f7ef909a182abf637bb0d:mlnscript.exe:TRUE:FALSE:win
2.29:Dec 2013:5f0a87e6cb7198d796f9664a05d5031a:mlwin.exe:FALSE:FALSE:win
2.29:Dec 2013:5afdf13c0406202aaf308b569052dd20:mlwin.exe:FALSE:TRUE:win
2.29:Dec 2013:47fbc35bf375d56d2291a3f85d2d838c:mlnscript.exe:FALSE:FALSE:win
2.29:Dec 2013:4d39f330c201e7614df17150f8aab74f:mlnscript.exe:TRUE:FALSE:win
2.30:Feb 2014:869c73b95daf1ec92c2b22277bd94724:mlwin.exe:FALSE:FALSE:win
2.30:Feb 2014:022ba981c2bf8751dad35c041f5f7db3:mlwin.exe:FALSE:TRUE:win
2.30:Feb 2014:b0f739262853e594242a6d4dad296eb6:mlnscript.exe:FALSE:FALSE:win
2.30:Feb 2014:c964df5ff4011eae94419c2f815a9450:mlnscript.exe:TRUE:FALSE:win
2.31:Sep 2014:befc087bb0e2b13ed01a57afa2d85bbe:mlwin.exe:FALSE:FALSE:win
2.31:Sep 2014:6038ba228ddde891b4673cae4b7aaa0c:mlwin.exe:FALSE:TRUE:win
2.31:Sep 2014:bfa10218aa4635ea2e5a4197faef98e7:mlnscript.exe:FALSE:FALSE:win
2.31:Sep 2014:3a4c5904a21788262ef8244958eb5302:mlnscript.exe:TRUE:FALSE:win
2.32:Jan 2015:eb320148ff952f5016c2aa4de7d8f363:mlwin.exe:FALSE:FALSE:win
2.32:Jan 2015:c3dbed5d07e14bd73fa0491a22749101:mlwin.exe:FALSE:TRUE:win
2.32:Jan 2015:7b86ba340c85d18d7e8ebd26fc30f09a:mlnscript.exe:FALSE:FALSE:win
2.32:Jan 2015:68541f384170dfe5f8a24adf187b9902:mlnscript.exe:TRUE:FALSE:win
2.32:Jan 2015:df5fabbed6204dd0277cad89b86c508a:mlnscript:TRUE:FALSE:lin
2.32:Jan 2015:fa29135101278f91c52bdb3f29773f8e:mlnscript:TRUE:FALSE:lin
2.32:Jan 2015:00e65cca7554bc5e82d21e0be564fe08:mlnscript:TRUE:FALSE:lin
2.32:Jan 2015:83d5b13cd8160a43e34603f35ef56ea3:mlnscript:TRUE:FALSE:lin
2.32:Jan 2015:e639bed035c855710d12004149f63465:mlnscript:TRUE:FALSE:lin
2.32:Jan 2015:47e627c655c52f7dc6ff1e76a670afe3:mlnscript:FALSE:FALSE:lin
2.32:Jan 2015:ab00a1cb783cf00ed6d95bbd5cebfb1a:mlnscript:TRUE:FALSE:mac
2.32:Jan 2015:9dc79007c3ac07cc04d7c34d8d936be6:mlnscript:TRUE:FALSE:bsd
2.33:May 2015:7bc55103dd0e093cedb3a61f1d297058:mlwin.exe:FALSE:FALSE:win
2.33:May 2015:1fa938cccf35f73669e55c2f1f022ff9:mlwin.exe:FALSE:TRUE:win
2.33:May 2015:c2c9953bbde950f11896ac1f0263c946:mlnscript.exe:FALSE:FALSE:win
2.33:May 2015:c561e82df447f972bb9ae564e8354d14:mlnscript.exe:TRUE:FALSE:win
2.33:May 2015:40b64d11a663b33516b1a1cb55a01c2b:mlnscript:TRUE:FALSE:lin
2.33:May 2015:97b235581ebc0b8af41747a3752323b1:mlnscript:TRUE:FALSE:lin
2.33:May 2015:8e47087a6ba730a4229b51c5bf7eacad:mlnscript:TRUE:FALSE:lin
2.33:May 2015:30d3a6d8e2ff0a9e21e601e72cecbf04:mlnscript:TRUE:FALSE:lin
2.33:May 2015:414c50ed7ccf6b68f50045092bcb5ac6:mlnscript:TRUE:FALSE:lin
2.33:May 2015:c3a7fbdb2ab067455443cfed4b620e87:mlnscript:TRUE:FALSE:lin
2.33:May 2015:bdb5089f1e25075824ba78d18006368a:mlnscript:TRUE:FALSE:lin
2.33:May 2015:7ceb282cf4e2c30cf071c485208800df:mlnscript:TRUE:FALSE:lin
2.33:May 2015:92e635a28f302f22470d552caeaa8cdc:mlnscript:TRUE:FALSE:lin
2.33:May 2015:5d8707133bad55b228909361b1a78b77:mlnscript:FALSE:FALSE:lin
2.33:May 2015:3f9e44e33daa6d63d8f6490382c85edc:mlnscript:TRUE:FALSE:mac
2.33:May 2015:c9c76531c4b01de9739b266eb1fa61c0:mlnscript:TRUE:FALSE:bsd
2.34:Jul 2015:513a13ad9ab8af09ffbc5995ff4dcffc:mlwin.exe:FALSE:FALSE:win
2.34:Jul 2015:2a891ff5bf102670774d55fdcc6fede6:mlwin.exe:FALSE:TRUE:win
2.34:Jul 2015:70b7ab62acbf003ae479a9235d20320b:mlnscript.exe:FALSE:FALSE:win
2.34:Jul 2015:5c1c2c78b965ba8c9278f83dba25e84c:mlnscript.exe:TRUE:FALSE:win
2.34:Jul 2015:a732fe7fdaed42b8c0d0f0428fb6768f:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:528298a0f57bad1d4fccd0daacf73912:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:120b475d744d8ee57186a91e87f77ccd:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:923cca1fd480143b1b3103a2c53e1b7f:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:d45cfa5489824870ffbc3bfa09140fd3:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:3c4fb8a62e24f27f8144c526382077cd:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:c197a4ec89768d9187929cc226447c05:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:fc3ad3ea9e44301a2d19f367348f803a:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:70749bcfcd9b4991a87f5b100f189e57:mlnscript:TRUE:FALSE:lin
2.34:Jul 2015:174b02d8bd5c78f23705738d49055984:mlnscript:FALSE:FALSE:lin
2.34:Jul 2015:30b35a64f45ebeff24bc98bdbdd5b149:mlnscript:TRUE:FALSE:mac
2.34:Jul 2015:072c3ac63fd07dff5f7396f6f336b995:mlnscript:TRUE:FALSE:bsd
2.35:Sep 2015:73fefcab85d5be673a6ba343dba5e49c:mlwin.exe:FALSE:FALSE:win
2.35:Sep 2015:eddefa11f05571290e6ffa247ffa2a8e:mlwin.exe:FALSE:TRUE:win
2.35:Sep 2015:1210ccb08f1d447109144830abb2b340:mlnscript.exe:FALSE:FALSE:win
2.35:Sep 2015:e51eb5d08247c12d6beef3c96608f767:mlnscript.exe:TRUE:FALSE:win
2.35:Sep 2015:6a7a4666f37fd3a707c5dea50aa8ebd0:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:b7d7c6fcc43c96d21accf8a86a86ae72:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:b095b8f9c205ab958552385f48cb2122:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:bc41e864ab1b9663bd6e0531317b0500:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:5670b2b3326d3da3d6be2cf84ba64f3c:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:2cc7e0dc180dd877706a7621fac1e755:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:83394f6ecbaf3842411fb4f29290b3a2:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:9c66d8b886b7a1f99a31e90a0eabfc98:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:f2f1ec1127cc88c583eab87d4f0b862a:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:896b3ccd73160abc2399d232dfd5a9a7:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:06aa5551c073948dff37a153c5404e24:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:a4a58452bef03243cfa873f0f2c5c6eb:mlnscript:TRUE:FALSE:lin
2.35:Sep 2015:521cd33275927f3dae57bedabf01124f:mlnscript:FALSE:FALSE:lin
2.35:Sep 2015:0d40e91fd5cb9bd6361d77e5316bf79f:mlnscript:FALSE:FALSE:lin
2.35:Sep 2015:807dccfc85dccb2955c38c953646c8d3:mlnscript:TRUE:FALSE:mac
2.35:Sep 2015:9b48a5dc3d5cb675a2c253d5009aa184:mlnscript:TRUE:FALSE:bsd
2.35:Sep 2015:b062e8c9116a001a0b73e541465e27c7:mlnscript:TRUE:FALSE:bsd
2.36:Mar 2016:f52da589b411a6636a6bede5914ce952:mlwin.exe:FALSE:FALSE:win
2.36:Mar 2016:f17aa345e767edc781f23e877d2e11ad:mlwin.exe:TRUE:FALSE:win
2.36:Mar 2016:fb25c32b4db90480789ec9208abf721a:mlwin.exe:FALSE:TRUE:win
2.36:Mar 2016:57fce6c83539daea0219e8916b0e1e40:mlnscript.exe:FALSE:FALSE:win
2.36:Mar 2016:ed98168c942104ad7c664315c834dd2b:mlnscript.exe:TRUE:FALSE:win
2.36:Mar 2016:3d73d31264a51f2fa0a2399727aa015a:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:3c4626d808d6b35c1f28909e9cec2034:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:e117a03eed7a297da7c255dce912b8be:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:30ea9ee4dad9ca7c5386a443c7802489:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:1a0abf5b2705dbb6ed66928ce9e689d5:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:bdbe97803b90f107b53b15ee538221c5:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:f03d43516a22f85295090d4244cdd62a:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:f19187c8f2c921b549e3162e92dc2962:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:7983bd105f45456f99cea2b0428bc2c2:mlnscript:TRUE:FALSE:lin
2.36:Mar 2016:df7f78276f22ee722ffa371c2fdf4321:mlnscript:FALSE:FALSE:lin
2.36:Mar 2016:8c33adfb5add5402a2df4c80c2d64183:mlnscript:TRUE:FALSE:mac
2.36:Mar 2016:88c5113d82d7013506c949c761689b65:mlnscript:TRUE:FALSE:bsd
2.36:Mar 2016:4b401e7a333ca3500959b72a6ed23afb:mlnscript:TRUE:FALSE:bsd
'
  versioninfo <- read.delim(textConnection(versioninfostr), header = TRUE, sep = ":", strip.white = TRUE)
  if (isTRUE(checkversion)) {
    # Allow disabling the version check if it is slowing things down (e.g. in a simulation study)
    currentver <- versioninfo[versioninfo$md5 == digest(cmd, algo = "md5", file = TRUE), ]
    if (nrow(currentver) == 0) {
      versiontext <- "MLwiN (version: unknown or >2.35)"
    } else {
      if (currentver$version < 2.28) {
        # Block versions >year older than current release
        stop("The current version of MLwiN is too old, please update it from http://www.bris.ac.uk/cmm/software/mlwin/download/upgrades.html")
      }
      versiontext <- paste0("MLwiN (version: ", currentver$version, ")")
    }
  } else {
    versiontext <- "MLwiN (version: unchecked)"
  }

  # the current function call
  cl <- match.call()

  centring <- estoptions$centring

  if (oldsyntax) {
    if (!is.null(centring)) {
      for (p in names(centring)) {
        if (as.integer(centring[[p]][1]) == 1) {
          indata[[p]] <- indata[[p]] - mean(indata[[p]])
        }
        if (as.integer(centring[[p]][1]) == 2) {
          grp <- as.factor(indata[[centring[[p]][2]]])
          indata[[p]] <- indata[[p]] - tapply(indata[[p]], grp, mean, na.rm = TRUE)[grp]  #mean(indata[[p]][as.logical(indata[[centring[[p]][2]]])])
        }
        if (as.integer(centring[[p]][1]) == 3) {
          indata[[p]] <- indata[[p]] - as.integer(centring[[p]][2])
        }
      }
    }
    invars <- Formula.translate.compat(Formula, levID, D, indata)
  } else {
    if (!is.null(centring)) {
      stop("Centring is not supported for new syntax")
    }
    invars <- Formula.translate(Formula, D, indata)
    newdata <- invars$indata
    invars$indata <- NULL
    if (!is.null(newdata)) {
      newvars <- setdiff(colnames(newdata), colnames(indata))
      for (var in newvars) {
        indata[[var]] <- newdata[[var]]
      }
    }
    rm(newdata)
  }

  resp <- invars$resp

  if (D[1] == "Normal") {
    if (is.factor(indata[[resp]])) {
      warning("You have specified a factor variable as a continuous response, you may wish to check its numeric values")
    }
  }

  if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
    if (!is.factor(indata[[resp]])) {
      indata[[resp]] <- as.factor(indata[[resp]])
    }
  }

  if (D[1] == "Binomial") {
    if (is.numeric(indata[[resp]])) {
      if (!all(is.na(indata[[resp]]) | (indata[[resp]] >= 0 & indata[[resp]] <= 1))) {
        stop("Binomial response variable must have values from zero to one")
      }
    }
    if (is.logical(indata[[resp]])) {
      indata[[resp]] <- as.integer(indata[[resp]])
    }
    if (is.factor(indata[[resp]])) {
      if (length(levels(indata[[resp]])) == 2) {
        indata[[resp]] <- as.integer(indata[[resp]] == max(levels(indata[[resp]])))
      } else {
        stop("Binomial responses must have two unique values")
      }
    }
  }

  if (D[1] == "Poisson" || D[1] == "Negbinom") {
    if (!all(is.na(indata[[resp]]) | (indata[[resp]] >= 0 & indata[[resp]] %% 1 == 0))) {
      stop("Poisson and Negative-binomial responses must be positive integers")
    }
  }

  if (D[1] == "Mixed") {
    for (i in 2:length(D)) {
      if (D[[i]][[1]] == "Normal") {
        if (is.factor(indata[[resp[i - 1]]])) {
          warning("You have specified a factor variable as a continuous response, you may wish to check its numeric values")
        }
      }
      if (D[[i]][[1]] == "Binomial") {
        if (is.numeric(indata[[resp[i - 1]]])) {
          if (!all(is.na(indata[[resp[i - 1]]]) | (indata[[resp[i - 1]]] >= 0 & indata[[resp[i - 1]]] <= 1))) {
            stop("Binomial response variable must have values from zero to one")
          }
        }
        if (is.logical(indata[[resp[i - 1]]])) {
          indata[[resp[i - 1]]] <- as.integer(indata[[resp[i - 1]]])
        }
        if (is.factor(indata[[resp[i - 1]]])) {
          if (length(levels(indata[[resp[i - 1]]])) == 2) {
            indata[[resp[i - 1]]] <- as.integer(indata[[resp[i - 1]]] == max(levels(indata[[resp[i - 1]]])))
          } else {
            stop("Binomial responses must have two unique values")
          }
        }
      }
      if (D[[i]][[1]] == "Poisson" || D[[i]][[1]] == "Negbinom") {
        if (!all(is.na(indata[[resp[i - 1]]]) | (indata[[resp[i - 1]]] >= 0 & indata[[resp[i - 1]]] %% 1 == 0))) {
          stop("Poisson and Negative-binomial responses must be positive integers")
        }
      }
    }
  }

  expl <- invars$expl

  D <- invars$D
  if (EstM == 0) {
    if (!is.element(D[1], c("Normal", "Binomial", "Poisson", "Negbinom", "Multivariate Normal", "Mixed", "Multinomial"))) {
      stop(cat("Invalid distribution specified:", D[1], "\n"))
    }
    if (D[1] == "Binomial") {
      if (!is.element(D[2], c("logit", "probit", "cloglog"))) {
        stop(cat("Invalid link function specified:", D[2], "\n"))
      }
      if (is.na(D[3])) {
        stop("A denominator must be specified for a Binomial response")
      }
    }
    if (D[1] == "Poisson") {
      if (!is.element(D[2], c("log"))) {
        stop(cat("Invalid link function specified:", D[2], "\n"))
      }
    }
    if (D[1] == "Mixed") {
      for (i in 2:length(D)) {
        if (!is.element(D[[i]][[1]], c("Normal", "Binomial", "Poisson"))) {
          stop(cat("Invalid distribution specified:", D[[i]][[1]], "\n"))
        }
        if (D[[i]][[1]] == "Binomial") {
          if (!is.element(D[[i]][[2]], c("logit", "probit", "cloglog"))) {
            stop(cat("Invalid link function specified:", D[[i]][[2]], "\n"))
          }
          if (is.na(D[[i]][[3]])) {
            stop("A denominator must be specified for Binomial responses")
          }
        }
        if (D[[i]][[1]] == "Poisson") {
          if (!is.element(D[[i]][[2]], c("log"))) {
            stop(cat("Invalid link function specified:", D[[i]][[2]], "\n"))
          }
        }
      }
    }
    if (D[1] == "Multinomial") {
      if (D[4] == 0) {
        # Unordered
        if (D[2] == "log") {
          warning("You specified as log link, but a logit link has been fitted")
          D[2] <- "logit"  # Backward compatibilty fix
        }
        if (D[2] != "logit") {
          stop("Invalid link function specified:", D[2])
        }
      }
      if (D[4] == 1) {
        # Ordered
        if (!is.element(D[2], c("logit", "probit", "cloglog"))) {
          stop(cat("Invalid link function specified:", D[2], "\n"))
        }
      }
    }
  } else {
    if (!is.element(D[1], c("Normal", "Binomial", "Poisson", "Negbinom", "Multivariate Normal", "Mixed", "Multinomial"))) {
      stop(cat("Invalid distribution specified:", D[1], "\n"))
    }
    if (D[1] == "Binomial") {
      if (!is.element(D[2], c("logit", "probit", "cloglog"))) {
        stop(cat("Invalid link function specified", D[2], "\n"))
      }
      if (is.na(D[3])) {
        stop("A denominator must be specified for a Binomial response")
      }
    }
    if (D[1] == "Poisson") {
      if (!is.element(D[2], c("log"))) {
        stop(cat("Invalid link function specified:", D[2], "\n"))
      }
    }
    if (D[1] == "Mixed") {
      for (i in 2:length(D)) {
        if (!is.element(D[[i]][[1]], c("Normal", "Binomial"))) {
          stop(cat("Invalid distribution specified:", D[[i]][[1]], "\n"))
        }
        if (D[[i]][[1]] == "Binomial") {
          if (!D[[i]][[2]] == "probit") {
            stop(cat("Invalid link function specified:", D[[i]][[2]], "\n"))
          }
          if (is.na(D[[i]][[3]])) {
            stop("A denominator must be specified for Binomial responses")
          }
        }
      }
    }
    if (D[1] == "Multinomial") {
      if (D[4] == 0) {
        # Unordered
        if (D[2] == "log")
          D[2] <- "logit"  # Backward compatibilty fix
        if (D[2] != "logit") {
          stop("Invalid link function specified", D[2], "\n")
        }
      }
      if (D[4] == 1) {
        # Ordered
        if (!is.element(D[2], c("logit", "probit", "cloglog"))) {
          stop(cat("Invalid link function specified", D[2], "\n"))
        }
      }
    }
  }

  if (is.null(levID)) {
    levID <- invars$levID
    isdiscrete <- D[1] %in% c("Binomial", "Poisson", "Negbinom", "Multinomial")
    if (D[1] == "Mixed") {
      for (i in 2:length(D)) {
        if (D[[i]][[1]] %in% c("Binomial", "Poisson")) {
          isdiscrete <- TRUE
        }
      }
    }
    if (isdiscrete) {
      indata[["l1id"]] <- 1:nrow(indata)
    }
  }

  rp <- invars$rp
  if (D[1] != "Normal" || D[1] != "Mixed") {
    if (D[1] == "Binomial") {
      if (any(rp$rp1 != "bcons.1"))
        stop("Variables cannot be made random at level one in Binomial models")
    }
    if (D[1] == "Poisson") {
      if (any(rp$rp1 != "bcons.1"))
        stop("Variables cannot be made random at level one in Poisson models")
    }
    if (D[1] == "Negbinom") {
      if (suppressWarnings(any(rp$rp1 != c("bcons.1", "bcons2.1"))))
        stop("Variables cannot be made random at level one in Negative-binomial models")
    }
    if (D[1] == "Mixed") {
      for (i in 2:length(D)) {
        if (D[[i]][[1]] == "Binomial") {
          rp1 <- rp$rp1[names(rp$rp1) == resp[i - 1]]
          if (length(rp1) > 0)
            stop("Variables cannot be made random at level one for Binomial responses")
        }
        if (D[[i]][[1]] == "Poisson") {
          rp1 <- rp$rp1[names(rp$rp1) == resp[i - 1]]
          if (length(rp1) > 0)
            stop("Variables cannot be made random at level one for Poisson responses")
        }
      }
    }
    if (D[1] == "Multinomial") {
      if (!is.null(rp$rp1)) {
        stop("Variables cannot be made random at level one in multinomial models")
      }
    }
  }
  nonfp <- invars$nonfp
  if (is.null(nonfp)) {
    if (is.list(expl)) {
      nonfp <- list(nonfp.sep = NA, nonfp.common = NA)
    } else {
      nonfp <- NA
    }
  }
  categ <- invars$categ
  if (is.null(categ)) {
    categ <- NULL
  } else {
    if (!isTRUE(oldsyntax)) {
      stop("categ not supported in new syntax")
    }
    if (is.null(rownames(categ))) {
      rownames(categ) <- c("var", "ref", "ncateg")
    }
  }

  if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
    levID <- c(levID, NA)
  }

  needsint <- FALSE
  if (is.list(expl)) {
    if ("1" %in% expl$sep.coeff || "1" %in% expl$common.coeff) {
      needsint <- TRUE
      expl$sep.coeff[expl$sep.coeff == "1"] <- "Intercept"
      expl$common.coeff[expl$common.coeff == "1"] <- "Intercept"
      rownames(expl$common.coeff.id) <- expl$common.coeff
    }
  } else {
    if ("1" %in% expl) {
      needsint <- TRUE
      expl[expl == "1"] <- "Intercept"
    }
  }

  if (is.list(nonfp)) {
    if ("1" %in% nonfp$nonfp.sep) {
      needsint <- TRUE
      nonfp$nonfp.sep[nonfp$nonfp.sep == "1"] <- "Intercept"
    }
    pos.cvar.one <- grepl("^1\\.{1}", nonfp$nonfp.common)
    if (any(pos.cvar.one)) {
      needsint <- TRUE
      nonfp$nonfp.common[pos.cvar.one] <- gsub("^1\\.", "Intercept.", nonfp$nonfp.common[pos.cvar.one])
    }
  } else {
    if ("1" %in% nonfp) {
      needsint <- TRUE
      nonfp[nonfp == "1"] <- "Intercept"
    }
  }

  for (ii in 1:length(rp)) {
    pos.cvar.one <- grepl("^1\\.{1}", rp[[ii]])
    if (any(pos.cvar.one)) {
      needsint <- TRUE
      rp[[ii]][pos.cvar.one] <- gsub("^1\\.", "Intercept.", rp[[ii]][pos.cvar.one])
    }
    if ("1" %in% rp[[ii]]) {
      needsint <- TRUE
      rp[[ii]][rp[[ii]] == "1"] <- "Intercept"
    }
  }

  if (isTRUE(needsint)) {
    indata[["Intercept"]] <- rep(1, nrow(indata))
  }

  show.file <- estoptions$show.file
  if (is.null(show.file))
    show.file <- FALSE

  resi.store <- estoptions$resi.store
  if (is.null(resi.store))
    resi.store <- FALSE

  resioptions <- estoptions$resioptions
  if (is.null(resioptions)) {
    resioptions <- c("variance")
  }
  if ("variance" %in% resioptions && ("standardised" %in% resioptions || "deletion" %in% resioptions || "leverage" %in%
                                        resioptions)) {
    stop("variance will not be calculated together with standardised or deletion or leverage. Please remove variance in resioptions, and then standard error will be calculated instead.")
  }

  if (EstM == 1 && ("sampling" %in% resioptions || "leverage" %in% resioptions || "influence" %in% resioptions ||
                      "deletion" %in% resioptions || "norecode" %in% resioptions)) {
    stop("Invalid residual option specified for MCMC estimation")
  }

  resi.store.levs <- estoptions$resi.store.levs
  if (EstM == 0 && !is.null(resi.store.levs)) {
    stop("resi.store.levs option is not valid for (R)IGLS estimation")
  }

  fpsandwich <- estoptions$fpsandwich
  rpsandwich <- estoptions$rpsandwich

  weighting <- estoptions$weighting
  if (EstM == 1 && !is.null(weighting)) {
    stop("weighting option is not valid for MCMC estimation")
  }
  if (!is.null(weighting) && (!is.element(D[1], c("Normal", "Poisson", "Binomial", "Negbinom")))) {
    stop("weighting can only be used in univariate models")
  }
  # Convert old weight syntax
  if (!is.null(weighting) && (!is.null(weighting$levels) || !is.null(weighting$weights) || !is.null(weighting$mode) ||
                                !is.null(weighting$FSDE) || !is.null(weighting$RSDE))) {
    warning("Old IGLS weighting options specified, see the help file for new syntax")
    if (weighting$FSDE == 2)
      fpsandwich <- TRUE
    if (weighting$RSDE == 2)
      rpsandwich <- TRUE
    if (weighting$mode == 2)
      weighting$standardised <- TRUE
    if (!is.null(weighting$levels) && !is.null(weighting$weights)) {
      if (length(weighting$levels) != length(weighting$weights)) {
        stop("The length of levels does not match with the length of weights.")
      }
      numlev <- length(levID)
      # weighting$weightvar = rep(NA, numlev)
      for (i in 1:length(weighting$levels)) {
        weighting$weightvar[[(numlev - weighting$levels[i]) + 1]] <- weighting$weights[i]
      }
    }
  }
  # Extract weights and add to output data
  if (!is.null(weighting)) {
    for (i in 1:length(weighting$weightvar)) {
      if (!is.na(weighting$weightvar[i])) {
        if (is.character(weighting$weightvar[[i]])) {
          wtvar <- model.frame(as.formula(paste0("~", weighting$weightvar[[i]])), data = data, na.action = NULL)
          indata <- cbind(indata, wtvar)
        } else {
          if (is.vector(weighting$weightvar[[i]])) {
            indata <- cbind(indata, weighting$weightvar[[i]])
            wtname <- paste0("_WEIGHT", (length(weighting$weightvar) - i) + 1)
            vnames <- colnames(indata)
            vnames[length(vnames)] <- wtname
            colnames(indata) <- vnames
            weighting$weightvar[[i]] <- wtname
          } else {
            stop("Invalid weights specification")
          }
        }
      }
    }
    if (is.null(fpsandwich))
      fpsandwich <- TRUE
    if (is.null(rpsandwich))
      rpsandwich <- TRUE
    if (is.null(weighting$standardised))
      weighting$standardised <- TRUE
  }

  if (is.null(fpsandwich))
    fpsandwich <- FALSE
  if (is.null(rpsandwich))
    rpsandwich <- FALSE

  clean.files <- estoptions$clean.files
  if (is.null(clean.files))
    clean.files <- TRUE

  clre <- estoptions$clre
  clre[2, ] <- gsub("^1$", "Intercept", clre[2, ])
  clre[3, ] <- gsub("^1$", "Intercept", clre[3, ])

  # There is no covariance between the negative-binomial bcons terms
  if (D[1] == "Negbinom") {
    clre <- cbind(clre, c(1, "bcons.1", "bcons2.1"))
  }

  smat <- estoptions$smat

  if (!is.null(smat)) {
    if (!is.matrix(smat)) {
      smat <- as.matrix(smat)
    }
    for (i in 1:ncol(smat)) {
      if (smat[2, i] == 1) {
        lev <- smat[1, i]
        for (j in 1:length(rp[[paste0("rp", lev)]])) {
          for (k in 1:j) {
            if (j != k)
              clre <- cbind(clre, c(lev, rp[[paste0("rp", lev)]][j], rp[[paste0("rp", lev)]][k]))
          }
        }
      }
    }
    smat <- NULL
  }


  mem.init <- estoptions$mem.init
  if (is.null(mem.init))
    mem.init <- "default"

  optimat <- estoptions$optimat
  if (is.null(optimat))
    optimat <- FALSE

  maxiter <- estoptions$maxiter
  if (is.null(maxiter))
    maxiter <- 20

  convtol <- estoptions$tol
  if (is.null(convtol))
    convtol <- 2

  extra <- estoptions$extra
  if (is.null(extra))
    extra <- FALSE
  if (extra == TRUE) {
    if (EstM != 0) {
      stop("extra can only be specified for (R)IGLS models")
    }
    if (D[1] != "Poisson" && D[1] != "Binomial" && D[1] != "Negbinom" && D[1] != "Multinomial") {
      stop("extra can only be specified for discrete outcome models")
    }
  }

  nonlinear <- estoptions$nonlinear
  if (!is.null(nonlinear)) {
    if (D[1] == "Normal" || D[1] == "Multivariate Normal") {
      stop("Nonlinear options are not valid for normal models")
    }
  }
  if (is.null(nonlinear))
    nonlinear <- c(0, 1)
  if (length(na.omit(levID)) == 1 && any(nonlinear != c(0, 1))) {
    stop("Only MQL1 is valid for one-level discrete models")
  }

  if (nonlinear[1] < 0 || nonlinear[1] > 1) {
    stop("Invalid nonlinear option")
  }

  if (nonlinear[2] < 1 || nonlinear[2] > 2) {
    stop("Invalid nonlinear option")
  }

  Meth <- estoptions$Meth
  if (is.null(Meth))
    Meth <- 1

  reset <- estoptions$reset
  if (is.null(reset)) {
    if (EstM == 0) {
      reset <- rep(0, length(na.omit(levID)))
      reset[1] <- 2
    }
  }
  if (!is.null(reset)) {
    if (EstM == 1) {
      stop("reset is only available for (R)IGLS models")
    } else {
      if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
        reset = c(0, reset)
      }
      if (length(reset) != length(levID)) {
        stop("reset vector is wrong length")
      }
      if (any(reset < 0 || reset > 2)) {
        stop("Invalid reset value")
      }
    }
  }

  fcon <- NULL
  rcon <- NULL
  if (!is.null(estoptions$constraints)) {
    if (EstM == 1) {
      stop("constraints are not available for MCMC estimation")
    } else {
      if (!is.null(estoptions$constraints$fixed.ui) || !is.null(estoptions$constraints$fixed.ci)) {
        if (is.null(estoptions$constraints$fixed.ui) || is.null(estoptions$constraints$fixed.ci)) {
          stop("fixed.ui and fixed.ci must be specified for fixed constaints")
        }
        if (ncol(estoptions$constraints$fixed.ui) != length(estoptions$constraints$fixed.ci)) {
          stop("number of columns in fixed.ui and fixed.ci must be equal for fixed constaints")
        }

        if (!is.list(expl)) {
          numfp <- length(setdiff(expl, nonfp))
        } else {
          numfp <- length(setdiff(expl$sep.coeff, nonfp$sep)) * length(resp) + setdiff(expl$common.coeff,
                                                                                       nonfp.common)
        }

        if (!is.vector(estoptions$constraints$fixed.ci)) {
          stop("fixed.ci should be a vector")
        }

        if (nrow(estoptions$constraints$fixed.ui) != numfp) {
          stop(paste("number of rows for fixed.ci must equal the number of fixed parameters:", numfp))
        }

        fcon <- rbind(estoptions$constraints$fixed.ui, estoptions$constraints$fixed.ci)
      }
      if (!is.null(estoptions$constraints$random.ui) || !is.null(estoptions$constraints$random.ci)) {
        if (is.null(estoptions$constraints$random.ui) || is.null(estoptions$constraints$random.ci)) {
          stop("random.ui and random.ci must be specified for fixed constaints")
        }
        if (ncol(estoptions$constraints$random.ui) != length(estoptions$constraints$random.ci)) {
          stop("number of columns in random.ui and random.ci must be equal for random constaints")
        }

        numrp <- 0
        for (ii in 1:length(rp)) {
          numrp <- numrp + (length(rp[[ii]]) * (length(rp[[ii]]) + 1)/2)
        }

        if (!is.null(clre)) {
          numrp <- numrp - ncol(clre)
        }

        if (!is.vector(estoptions$constraints$random.ci)) {
          stop("random.ci should be a vector")
        }

        if (nrow(estoptions$constraints$random.ui) != numrp) {
          stop(paste("number of rows for random.ui must equal the number of random parameters", numrp))
        }

        rcon <- rbind(estoptions$constraints$random.ui, estoptions$constraints$random.ci)
      }
    }
  }

  fact <- estoptions$fact
  if (EstM == 0 && !is.null(fact)) {
    stop("Factor models not available with (R)IGLS estimation")
  }
  if (!is.null(fact) && D[1] != "Multivariate Normal") {
    stop("Factor models can only be defined for multivariate models")
  }

  if (!is.null(fact)) {
    for (i in 1:fact$nfact) {
      for (j in 1:length(rp[[paste0("rp", fact$lev.fact[i])]])) {
        for (k in 1:j) {
          if (j != k)
            clre <- cbind(clre, c(fact$lev.fact[i], rp[[paste0("rp", fact$lev.fact[i])]][j], rp[[paste0("rp",
                                                                                                        fact$lev.fact[i])]][k]))
        }
      }
    }
  }

  xclass <- estoptions$xclass

  if (!is.null(xclass)) {
    warning("Old cross-classification option specified, see the help file for new syntax")
  }

  if (EstM == 0 && !is.null(xclass)) {
    stop("Cross-classification is not available with (R)IGLS estimation")
  }

  xc <- estoptions$xc
  if (is.null(xc)) {
    xc <- FALSE
  }

  mm <- estoptions$mm
  if (!is.null(mm)) {
    xc <- TRUE
    for (i in 1:length(mm)) {
      if (!any(is.na(mm[[i]]))) {
        if (is.matrix(mm[[i]]) || is(mm[[i]], "sparseMatrix")) {
          idstub <- paste0("id_", length(mm) - i, "_")
          weightstub <- paste0("weight_", length(mm) - i, "_")
          mmdf <- matrix2df(mm[[i]], standardise = TRUE, idstub = idstub, weightstub = weightstub)
          indata <- cbind(indata, mmdf)
          numids <- ncol(mmdf)/2
          mm[[i]] <- list(mmvar = paste0(idstub, 1:numids), weights = paste0(weightstub, 1:numids))
          levID[i] <- paste0(idstub[1], 1)
        } else {
          if (length(mm[[i]]$mmvar) != length(mm[[i]]$weights)) {
            stop("Number of multiple membership weights does not match identifiers")
          }
          mmlev <- (length(mm) - i) + 1
          for (j in 1:length(mm[[i]]$mmvar)) {
            var <- mm[[i]]$mmvar[[j]]
            if (is.character(var)) {
              if (var %in% colnames(indata))
                indata[[var]] <- NULL
              mmvar <- model.frame(as.formula(paste0("~", var)), data = data, na.action = NULL)
              indata <- cbind(indata, mmvar)
            } else {
              if (is.vector(var)) {
                indata <- cbind(indata, var)
                mmname <- paste0("_MMID", (length(mm) - i) + 1, "_", j)
                vnames <- colnames(indata)
                vnames[length(vnames)] <- mmname
                colnames(indata) <- vnames
                mm[[i]]$mmvar[[j]] <- mmname
              } else {
                stop("Invalid weights specification")
              }
            }
          }
          for (j in 1:length(mm[[i]]$weights)) {
            var <- mm[[i]]$weights[[j]]
            if (is.character(var)) {
              if (var %in% colnames(indata))
                indata[[var]] <- NULL
              mmweight <- model.frame(as.formula(paste0("~", var)), data = data, na.action = NULL)
              indata <- cbind(indata, mmweight)
            } else {
              if (is.vector(var)) {
                indata <- cbind(indata, var)
                mmwname <- paste0("_MMWEIGHT", (length(mm) - i) + 1, "_", j)
                vnames <- colnames(indata)
                vnames[length(vnames)] <- mmwname
                colnames(indata) <- vnames
                mms[[i]]$weight[[j]] <- mmwname
              } else {
                stop("Invalid weights specification")
              }
            }
          }
        }
      }
    }
  }

  car <- estoptions$car
  if (!is.null(car)) {
    xc <- TRUE
    for (i in 1:length(car)) {
      if (!any(is.na(car[[i]]))) {
        if (is.matrix(car[[i]]) || is(car[[i]], "sparseMatrix")) {
          idstub <- paste0("id_", length(car) - i, "_")
          weightstub <- paste0("weight_", length(car) - i, "_")
          cardf <- matrix2df(car[[i]], standardise = TRUE, idstub = idstub, weightstub = weightstub)
          indata <- cbind(indata, cardf)
          numids <- ncol(cardf)/2
          car[[i]] <- list(carvar = paste0(idstub, 1:numids), weights = paste0(weightstub, 1:numids))
          levID[i] <- paste0(idstub[1], 1)
        } else {
          if (length(car[[i]]$carvar) != length(car[[i]]$weights)) {
            stop("Number of CAR weights does not match identifiers")
          }
          carlev <- (length(car) - i) + 1
          for (j in 1:length(car[[i]]$carvar)) {
            var <- car[[i]]$carvar[[j]]
            if (is.character(var)) {
              if (var %in% colnames(indata)) {
                indata[[var]] <- NULL
              }
              carvar <- model.frame(as.formula(paste0("~", var)), data = data, na.action = NULL)
              indata <- cbind(indata, carvar)
            } else {
              if (is.vector(var)) {
                indata <- cbind(indata, var)
                carname <- paste0("_CARID", (length(car) - i) + 1, "_", j)
                vnames <- colnames(indata)
                vnames[length(vnames)] <- carname
                colnames(indata) <- vnames
                car[[i]]$carvar[[j]] <- carname
              } else {
                stop("Invalid weights specification")
              }
            }
          }
          for (j in 1:length(car[[i]]$weights)) {
            var <- car[[i]]$weights[[j]]
            if (is.character(var)) {
              if (var %in% colnames(indata))
                indata[[var]] <- NULL
              carweight <- model.frame(as.formula(paste0("~", var)), data = data, na.action = NULL)
              indata <- cbind(indata, carweight)
            } else {
              if (is.vector(var)) {
                indata <- cbind(indata, var)
                carwname <- paste0("_CARWEIGHT", (length(car) - i) + 1, "_", j)
                vnames <- colnames(indata)
                vnames[length(vnames)] <- carwname
                colnames(indata) <- vnames
                car[[i]]$weights[[j]] <- carwname
              } else {
                stop("Invalid weights specification")
              }
            }
          }
        }
      }
    }
  }

  carcentre <- estoptions$carcentre
  if (is.null(carcentre)) {
    carcentre <- FALSE
  }

  # Convert old version of syntax
  if (!is.null(xclass)) {
    xc <- TRUE
    for (i in 1:length(as.numeric(xclass$class))) {
      num <- as.numeric(xclass$N1[i])
      if (num > 1) {
        lev <- as.numeric(xclass$class[i])

        weightcol <- xclass$weight[i]
        idcol <- xclass$id[i]
        if (is.null(idcol) || is.na(idcol)) {
          idcol <- rev(na.omit(levID))[lev]
        }
        idstart <- which(colnames(indata) == idcol)
        idend <- idstart + (num - 1)
        idcols <- colnames(indata)[idstart:idend]

        weightstart <- which(colnames(indata) == weightcol)
        weightend <- weightstart + (num - 1)
        weightcols <- colnames(indata)[weightstart:weightend]

        if (!isTRUE(xclass$car)) {
          if (is.null(car))
            car <- list(carvar = list(), weights = list())
          car[[lev]]$carvar <- as.list(idcols)
          car[[lev]]$weights <- as.list(weightcols)
        } else {
          if (is.null(mm))
            mm <- list(mmvar = list(), weights = list())
          mm[[lev]]$mmvar <- as.list(idcols)
          mm[[lev]]$weights <- as.list(weightcols)
        }
      }
    }
  }


  if (!is.null(car) && !is.null(mm)) {
    mmlev <- max(length(car), length(mm))
    for (i in 1:mmlev) {
      if (any(!is.na(mm)) && any(!is.na(car))) {
        stop("Multiple membership and CAR cannot be defined at the same level")
      }
    }
  }

  if (!is.null(car)) {
    carcount <- 0
    for (i in 1:length(car)) {
      if (!any(is.na(car[[i]]))) {
        carcount <- carcount + 1
      }
    }
    if (carcount > 1 && !isTRUE(carcentre)) {
      stop("CAR can only apply to one level unless CAR centring is turned on")
    }
  }


  if (!is.null(mm)) {
    for (k in 1:length(mm)) {
      if (!any(is.na(mm[[k]]))) {
        varnames <- unlist(mm[[k]]$mmvar)
        weightnames <- unlist(mm[[k]]$weights)
        idmat <- indata[, varnames]
        weightmat <- indata[, weightnames]
        if (!all(idmat[, 1] == indata[, levID[k]])) {
          stop("The first multiple membership column should match the ID column")
        }
        levID[k] <- mm[[k]]$mmvar[[1]]

        # NOTE: These checks could probably be vectorised
        for (i in 1:nrow(idmat)) {
          for (j in 1:ncol(idmat)) {
            if (idmat[i, j] != 0 && weightmat[i, j] == 0) {
              stop(paste("The MM ID variable", j, "for observation", i, "is present but a zero MM weight has been specified for it"))
            }
            if (idmat[i, j] == 0 && weightmat[i, j] != 0) {
              stop(paste("The MM ID variable", j, "for observation", i, "is absent but a positive MM weight has been specified for it"))
            }
            for (k in 1:j) {
              if (k != j) {
                if (idmat[i, j] == idmat[i, k] && idmat[i, j] != 0) {
                  stop(paste("The MM ID variable", j, "for observation", i, "is duplicated in ID variable", k))
                }
              }
            }
          }
        }
      }
    }
  }

  if (!is.null(car)) {
    for (i in 1:length(car)) {
      if (!any(is.na(car[[i]]))) {
        varnames <- unlist(car[[i]]$carvar)
        weightnames <- unlist(car[[i]]$weights)
        idmat <- indata[, varnames]
        weightmat <- indata[, weightnames]
        # NOTE: These checks could probably be vectorised
        for (i in 1:nrow(idmat)) {
          for (j in 1:ncol(idmat)) {
            if (idmat[i, j] != 0 && weightmat[i, j] == 0) {
              stop(paste("The CAR ID variable", j, "for observation", i, "is present but a zero MM weight has been specified for it"))
            }
            if (idmat[i, j] == 0 && weightmat[i, j] != 0) {
              stop(paste("The CAR ID variable", j, "for observation", i, "is absent but a positive MM weight has been specified for it"))
            }
            for (k in 1:j) {
              if (k != j) {
                if (idmat[i, j] == idmat[i, k] && idmat[i, j] != 0) {
                  stop(paste("The CAR ID variable", j, "for observation", i, "is duplicated in ID variable", k))
                }
              }
            }
          }
        }
      }
    }
  }

  notation <- estoptions$notation
  if (!("level" %in% notation && "class" %in% notation)) {
    if (!isTRUE(xc)) {
      notation <- c(notation, "level")
    } else {
      notation <- c(notation, "class")
    }
  }

  merr <- estoptions$merr
  if (EstM == 0 && !is.null(merr)) {
    stop("Measurement error is not available with (R)IGLS estimation")
  }
  if (!is.null(merr) && (!is.element(D[1], c("Normal", "Poisson", "Binomial", "Negbinom")))) {
    stop("measurement error can only be used in univariate models")
  }
  mcmcMeth <- estoptions$mcmcMeth
  if (EstM == 0 && !is.null(merr)) {
    stop("MCMC method options cannot be specified with (R)IGLS estimation")
  }

  FP.names <- NULL

  if (D[1] == "Multinomial") {
    nresp <- length(levels(indata[, resp])) - 1
    resp.names <- levels(indata[, resp])[-as.numeric(D[5])]

    if (is.list(expl)) {
      nonfp.sep <- nonfp$nonfp.sep
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
      if (!is.na(expl$sep.coeff[1])) {
        for (p in expl$sep.coeff) {
          if (is.na(nonfp.sep[1]) || sum(p == nonfp.s) == 0) {
            if (is.null(categ) || sum(p == categ["var", ]) == 0) {
              for (j in 1:nresp) {
                FP.names <- c(FP.names, paste("FP_", chartr(".", "_", p), "_", resp.names[j], sep = ""))
              }
            } else {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                categ.names <- levels(indata[[p]])
                for (j in 1:nresp) {
                  for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                    FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j], sep = ""))
                  }
                }
              } else {
                categ.names <- levels(indata[[p]])
                refx <- categ["ref", which(p == categ["var", ])]
                categ.names <- categ.names[-which(refx == categ.names)]
                for (j in 1:nresp) {
                  for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                    FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j], sep = ""))
                  }
                }
              }
            }
          }
        }
      }
      kk <- 1
      tempid <- 1:(nresp + 1)
      tempid <- tempid[-as.numeric(D["ref.cat"])]
      for (p in expl$common.coeff) {
        newp <- paste(p, paste(tempid[as.logical(expl$common.coeff.id[kk, ])], collapse = ""), sep = ".")
        kk <- kk + 1
        nonfp.common <- nonfp$nonfp.common
        nonfp.c <- nonfp.common
        if (is.na(nonfp.common[1]) || sum(newp == nonfp.c) == 0) {
          if (is.null(categ) || sum(p == categ["var", ]) == 0) {
            FP.names <- c(FP.names, paste("FP_", chartr(".", "_", newp), sep = ""))
          } else {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              categ.names <- levels(indata[[p]])
              for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), sep = ""))
              }
            } else {
              categ.names <- levels(indata[[p]])
              refx <- categ["ref", which(p == categ["var", ])]
              categ.names <- categ.names[-which(refx == categ.names)]
              for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), sep = ""))
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
              FP.names <- c(FP.names, paste("FP_", chartr(".", "_", p), "_", resp.names[j], sep = ""))
            }
          } else {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              categ.names <- levels(indata[[p]])
              for (j in 1:nresp) {
                for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                  FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j], sep = ""))
                }
              }
            } else {
              categ.names <- levels(indata[[p]])
              refx <- categ["ref", which(p == categ["var", ])]
              categ.names <- categ.names[-which(refx == categ.names)]
              for (j in 1:nresp) {
                for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                  FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j],
                                                sep = ""))
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
        nonfp.sep <- nonfp$nonfp.sep
        nonfp.s <- nonfp.sep
        for (i in 1:length(resp)) {
          nonfp.s <- gsub(paste(".", resp[i], sep = ""), "", nonfp.s)
        }
        nonfp.s <- unique(nonfp.s)
        if (!is.na(expl$sep.coeff[1])) {
          for (p in expl$sep.coeff) {
            if (is.na(nonfp.sep[1]) || sum(p == nonfp.s) == 0) {
              if (is.null(categ) || sum(p == categ["var", ]) == 0) {
                for (j in 1:nresp) {
                  FP.names <- c(FP.names, paste("FP_", chartr(".", "_", p), "_", resp[j], sep = ""))
                }
              } else {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  categ.names <- levels(indata[[p]])
                  for (j in 1:nresp) {
                    for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                      FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j],
                                                    sep = ""))
                    }
                  }

                } else {
                  categ.names <- levels(indata[[p]])
                  refx <- categ["ref", which(p == categ["var", ])]
                  categ.names <- categ.names[-which(refx == categ.names)]
                  for (j in 1:nresp) {
                    for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                      FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j],
                                                    sep = ""))
                    }
                  }
                }
              }
            }
          }
        }
        kk <- 1
        for (p in expl$common.coeff) {
          newp <- paste(p, paste(which(as.logical(expl$common.coeff.id[kk, ])), collapse = ""), sep = ".")
          kk <- kk + 1
          nonfp.common <- nonfp$nonfp.common
          nonfp.c <- nonfp.common
          # for (i in 1:length(nonfp.c)){ nonfp.c[i]=gsub('\\.[[:digit:]]+$','',nonfp.c[i]) }
          if (is.na(nonfp.common[1]) || sum(newp == nonfp.c) == 0) {
            if (is.null(categ) || sum(p == categ["var", ]) == 0) {
              FP.names <- c(FP.names, paste("FP_", chartr(".", "_", newp), sep = ""))
            } else {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                categ.names <- levels(indata[[p]])
                for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                  FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), sep = ""))
                }
              } else {
                categ.names <- levels(indata[[p]])
                refx <- categ["ref", which(p == categ["var", ])]
                categ.names <- categ.names[-which(refx == categ.names)]
                for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                  FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), sep = ""))
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
                FP.names <- c(FP.names, paste("FP_", chartr(".", "_", p), "_", resp[j], sep = ""))
              }

            } else {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                categ.names <- levels(indata[[p]])
                for (j in 1:nresp) {
                  for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                    FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j],
                                                  sep = ""))
                  }
                }
              } else {
                categ.names <- levels(indata[[p]])
                refx <- categ["ref", which(p == categ["var", ])]
                categ.names <- categ.names[-which(refx == categ.names)]
                for (j in 1:nresp) {
                  for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                    FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), "_", resp.names[j],
                                                  sep = ""))
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
            FP.names <- c(FP.names, paste("FP_", chartr(".", "_", p), sep = ""))
          } else {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              categ.names <- levels(indata[[p]])
              for (i in 1:as.numeric(categ["ncateg", which(p == categ["var", ])])) {
                FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), sep = ""))
              }
            } else {
              categ.names <- levels(indata[[p]])
              refx <- categ["ref", which(p == categ["var", ])]
              categ.names <- categ.names[-which(refx == categ.names)]
              for (i in 1:(as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1)) {
                FP.names <- c(FP.names, paste("FP_", chartr(".", "_", categ.names[i]), sep = ""))
              }
            }
          }
        }
      }
    }
  }

  resid.names <- function(rpx, resid.lev, RP) {
    nrpx <- length(rpx)
    for (j in 1:nrpx) {
      for (i in 1:j) {
        if (i == j) {
          RP <- c(RP, paste("RP", resid.lev, "_var_", chartr(".", "_", rpx[i]), sep = ""))
        } else {
          RP <- c(RP, paste("RP", resid.lev, "_cov_", chartr(".", "_", rpx[i]), "_", chartr(".", "_", rpx[j]), sep = ""))
        }
      }
    }
    RP
  }

  resid2.names <- function(rpx, resid.lev, clre, RP) {
    nrpx <- length(rpx)
    nclre <- ncol(clre)
    k <- 1
    for (j in 1:nrpx) {
      for (i in 1:j) {
        if (!any(as.numeric(clre[1, ]) == resid.lev & ((clre[2, ] == rpx[i] & clre[3, ] == rpx[j]) | (clre[2, ] == rpx[j] & clre[3, ] == rpx[i])))) {
          if (i == j) {
            RP <- c(RP, paste("RP", resid.lev, "_var_", chartr(".", "_", rpx[i]), sep = ""))
          } else {
            RP <- c(RP, paste("RP", resid.lev, "_cov_", chartr(".", "_", rpx[i]), "_", chartr(".", "_", rpx[j]),
                              sep = ""))
          }
        }
      }
    }
    RP
  }

  RP.names <- NULL
  if (length(rp) > 0) {
    for (ii in 1:length(rp)) {
      if (is.null(clre)) {
        RP.names <- resid.names(rp[[ii]], as.numeric(sub("rp", "", names(rp)[ii])), RP.names)
      } else {
        RP.names <- resid2.names(rp[[ii]], as.numeric(sub("rp", "", names(rp)[ii])), clre, RP.names)
      }
    }
  }


  if (D[1] == "Multinomial") {
    RP.names <- c(RP.names, "RP1_bcons_1")
    if (EstM == 0 && as.numeric(D[4]) == 0) {
      RP.names <- c(RP.names, "RP1_bcons_2")
    }
  }

  FP <- rep(0, length(FP.names))
  names(FP) <- FP.names

  FP.cov <- matrix(0, length(FP.names), length(FP.names))
  colnames(FP.cov) <- FP.names
  rownames(FP.cov) <- FP.names

  RP <- rep(0, length(RP.names))
  names(RP) <- RP.names

  RP.cov <- matrix(0, length(RP.names), length(RP.names))
  colnames(RP.cov) <- RP.names
  rownames(RP.cov) <- RP.names


  if (EstM == 1) {
    iterations <- mcmcMeth$iterations
    if (is.null(iterations))
      iterations <- 5000
    burnin <- mcmcMeth$burnin
    if (is.null(burnin))
      burnin <- 500
    thinning <- mcmcMeth$thinning
    if (is.null(thinning))
      thinning <- 1
    nchains <- mcmcMeth$nchains
    if (is.null(nchains))
      nchains <- 1
    seed <- mcmcMeth$seed
    if (is.null(seed))
      seed <- 1:nchains
    if (length(seed) != nchains)
      seed <- rep(seed, nchains)
    priorParam <- mcmcMeth$priorParam
    if (is.list(priorParam))
      priorParam <- prior2macro(priorParam, Formula, levID, D, indata)
    if (is.null(priorParam))
      priorParam <- "default"
    scale <- mcmcMeth$scale
    if (is.null(scale))
      scale <- 5.8
    refresh <- mcmcMeth$refresh
    if (is.null(refresh))
      refresh <- 50
    # Ensure that refresh is not greater than total iterations
    refresh <- min(refresh, iterations)
    fixM <- mcmcMeth$fixM
    if (is.null(fixM)) {
      if (D[1] == "Poisson" || D[1] == "Multinomial" || D[1] == "Binomial" || D[1] == "Mixed") {
        fixM <- 2
      } else {
        fixM <- 1
      }
    }
    residM <- mcmcMeth$residM
    if (is.null(residM)) {
      if (D[1] == "Poisson" || D[1] == "Multinomial" || D[1] == "Binomial" || D[1] == "Mixed") {
        residM <- 2
      } else {
        residM <- 1
      }
    }
    Lev1VarM <- mcmcMeth$Lev1VarM
    if (is.null(Lev1VarM)) {
      if (D[1] == "Poisson" || D[1] == "Multinomial" || D[1] == "Binomial" || D[1] == "Mixed") {
        Lev1VarM <- 2
      } else {
        Lev1VarM <- 1
      }
    }
    OtherVarM <- mcmcMeth$OtherVarM
    if (is.null(OtherVarM))
      OtherVarM <- 1
    adaption <- mcmcMeth$adaption
    if (is.null(adaption))
      adaption <- 1
    priorcode <- mcmcMeth$priorcode
    if (is.null(priorcode))
      priorcode <- c(gamma=1)
    if (is.null(names(priorcode))) {
      if (length(priorcode) == 1) {
        names(priorcode) <- "gamma"
      } else if (length(priorcode == 3)) { 
          names(priorcode) <- c("gamma", "shape", "scale")
      }
    }
    if (names(priorcode)[1] == "uniform") {
      priorcode[1] <- !priorcode[1]
      names(priorcode)[1] <- "gamma"
    }
    if (names(priorcode)[1] != "gamma") {
      stop("Invalid prior option")
    }
    if (length(priorcode) > 1) {
      if (length(setdiff(names(priorcode), c("gamma", "shape", "scale"))) != 0) {
        stop("Invalid prior options")
      }
    }
    priorcode[1] <- as.integer(priorcode[1])
    rate <- mcmcMeth$rate
    if (is.null(rate))
      rate <- 50
    tol <- mcmcMeth$tol
    if (is.null(tol))
      tol <- 10
    lclo <- mcmcMeth$lclo
    if (is.null(lclo))
      lclo <- 0
    dami <- mcmcMeth$dami
    if (!is.null(dami)) {
      if (dami[1] %in% c(0, 1, 2)) {
        if (dami[1] == 0) {
          if (length(dami) == 1) {
            stop("No iterations specified for imputation")
          }
        } else {
          if (length(dami) > 1) {
            stop("Iterations cannot be specified when requesting imputation means")
          }
        }
      } else {
        stop("Invalid imputation option")
      }
    }
  }

  mcmcOptions <- estoptions$mcmcOptions
  if (EstM == 0 && !is.null(merr)) {
    stop("MCMC options cannot be specified with (R)IGLS estimation")
  }

  if (EstM == 1) {
    if (!is.null(mcmcOptions)) {
      if (D[1] == "Multivariate Normal") {
        mcmcOptions2 <- list(orth = 0, hcen = 0, smcm = 0, smvn = 0, paex = c(2, 0), mcco = 0)
        for (ii in names(mcmcOptions)) mcmcOptions2[[ii]] <- mcmcOptions[[ii]]
        mcmcOptions <- mcmcOptions2
      } else {
        mcmcOptions2 <- list(orth = 0, hcen = 0, smcm = 0, smvn = 0, paex = c(2, 0))
        for (ii in names(mcmcOptions)) mcmcOptions2[[ii]] <- mcmcOptions[[ii]]
        mcmcOptions <- mcmcOptions2
      }
    } else {
      if (D[1] == "Multivariate Normal") {
        mcmcOptions <- list(orth = 0, hcen = 0, smcm = 0, smvn = 0, paex = c(2, 0), mcco = 0)
      } else {
        mcmcOptions <- list(orth = 0, hcen = 0, smcm = 0, smvn = 0, paex = c(2, 0))
      }
    }
    if (mcmcOptions$hcen > 0) {
      if (mcmcOptions$hcen < 2 || mcmcOptions$hcen > length(na.omit(levID))) {
        stop("Invalid level for hierarchical centring")
      }
      if (D[1] == "Multivariate Normal" || D[1] == "Mixed" || D[1] == "Multinomial") {
        mcmcOptions$hcen = mcmcOptions$hcen + 1
      }
    }
    if (is.matrix(mcmcOptions$paex)) {
      for (i in 1:nrow(mcmcOptions$paex)) {
        if (mcmcOptions$paex[i, 2] == 1) {
          pelev <- mcmcOptions$paex[i, 1]
          if (pelev < 2 || pelev > length(na.omit(levID))) {
            stop("Invalid level for parameter expansion")
          }
          if (D[1] == "Multivariate Normal" || D[1] == "Mixed" || D[1] == "Multinomial") {
            mcmcOptions$paex[i, 1] = mcmcOptions$paex[i, 1] + 1
          }
        }
      }
    } else {
      if (mcmcOptions$paex[2] == 1) {
        pelev <- mcmcOptions$paex[1]
        if (pelev < 2 || pelev > length(na.omit(levID))) {
          stop("Invalid level for parameter expansion")
        }
        if (D[1] == "Multivariate Normal" || D[1] == "Mixed" || D[1] == "Multinomial") {
          mcmcOptions$paex[1] = mcmcOptions$paex[1] + 1
        }
      }
    }
  }

  sval <- estoptions$startval
  if (EstM == 1) {
    if (!is.null(mcmcMeth$startval)) {
      warning("startval is now specified directly within estoptions")
      sval <- mcmcMeth$startval
    }
  }

  if (!is.null(sval)) {
    # Check if we have a single set of starting values, if so expand to a list
    if (any(c("FP.b", "FP.v", "RP.b", "RP.v") %in% names(sval))) {
      if (EstM == 0) {
        svals <- list(sval)
      }
      if (EstM == 1) {
        svals <- list()
        for (i in 1:nchains) {
          svals[[i]] <- sval
        }
      }
    } else {
      if (length(sval) != nchains) {
        stop("Start values list is wrong size")
      }
      svals <- sval
    }

    #startval <- rep(list(list()), length(svals))
    startval <- replicate(length(svals), list())
    for (i in 1:length(svals)) {
      if (!is.null(svals[[i]]$FP.b) && is.null(names(svals[[i]]$FP.b))) {
        names(svals[[i]]$FP.b) <- FP.names
      }
      if (!is.null(svals[[i]]$FP.v) && (is.null(rownames(svals[[i]]$FP.v)) || is.null(colnames(svals[[i]]$FP.v)))) {
        rownames(svals[[i]]$FP.v) <- FP.names
        colnames(svals[[i]]$FP.v) <- FP.names
      }
      if (!is.null(svals[[i]]$RP.b) && is.null(names(svals[[i]]$RP.b))) {
        names(svals[[i]]$RP.b) <- RP.names
      }
      if (!is.null(svals[[i]]$RP.v) && (is.null(rownames(svals[[i]]$RP.v)) || is.null(colnames(svals[[i]]$RP.v)))) {
        rownames(svals[[i]]$RP.v) <- RP.names
        colnames(svals[[i]]$RP.v) <- RP.names
      }
      sharedFP <- intersect(FP.names, names(svals[[i]]$FP.b))
      if (!is.null(svals[[i]]$FP.b) && !is.null(sharedFP)) {
        startval[[i]]$FP.b <- FP
        startval[[i]]$FP.b[sharedFP] <- svals[[i]]$FP.b[sharedFP]
      }
      if (!is.null(svals[[i]]$FP.v) && !is.null(sharedFP)) {
        startval[[i]]$FP.v <- FP.cov
        startval[[i]]$FP.v[sharedFP, sharedFP] <- svals[[i]]$FP.v[sharedFP, sharedFP]
      }
      sharedRP <- intersect(RP.names, names(svals[[i]]$RP.b))
      if (!is.null(svals[[i]]$RP.b) && !is.null(sharedRP)) {
        startval[[i]]$RP.b <- RP
        startval[[i]]$RP.b[sharedRP] <- svals[[i]]$RP.b[sharedRP]
      }
      if (!is.null(svals[[i]]$RP.v) && !is.null(sharedRP)) {
        startval[[i]]$RP.v <- RP.cov
        startval[[i]]$RP.v[sharedRP, sharedRP] <- svals[[i]]$RP.v[sharedRP, sharedRP]
      }
    }
  } else {
    startval <- NULL
  }

  if (EstM == 0 && !is.null(BUGO)) {
    stop("BUGO requires MCMC estimation to be selected")
  }

  sort.force <- estoptions$sort.force
  if (is.null(sort.force))
    sort.force <- FALSE

  sort.ignore <- estoptions$sort.ignore
  if (is.null(sort.ignore))
    sort.ignore <- FALSE

  if (D[1] == "Binomial") {
    if (!all(is.na(indata[[resp]]) | (indata[[resp]] >= 0 & indata[[resp]] <= 1))) {
      stop("All values for a binomial response must lie between zero and one")
    }
  }

  if (D[1] == "Poisson") {
    if (!all(is.na(indata[[resp]]) | (indata[[resp]] >= 0 & indata[[resp]] %% 1 == 0))) {
      stop("All values for a Poisson response must be positive integers")
    }
  }

  if (D[1] == "Mixed") {
    mixlink <- NULL
    discreteresp <- NULL
    for (i in 2:length(D)) {
      if (D[[i]][1] == "Binomial") {
        if (!all(is.na(indata[[resp[i - 1]]]) | (indata[[resp[i - 1]]] >= 0 & indata[[resp[i - 1]]] <= 1))) {
          stop("All values for a binomial response must lie between zero and one")
        }
        discreteresp <- union(discreteresp, D[[i]][1])
        mixlink <- union(mixlink, D[[i]][2])
      }

      if (D[[i]][1] == "Poisson") {
        if (!all(is.na(indata[[resp[i - 1]]]) | (indata[[resp[i - 1]]] >= 0 & indata[[resp[i - 1]]] %% 1 == 0))) {
          stop("All values for a Poisson response must be positive integers")
        }
        discreteresp <- union(discreteresp, D[[i]][1])
        mixlink <- union(mixlink, D[[i]][2])
      }
    }
    if (length(discreteresp) > 1) {
      stop("Mixed response models cannot contain both Binomial and Poisson responses")
    }
    if (length(mixlink) > 1) {
      stop("Only one link type can be specified for mixed models")
    }
  }

  if (D[1] == "Multinomial") {
    if (length(unique(indata[[resp]])) < 2) {
      stop("Responses must have at least two categories for multinomial models")
    }
    if (is.na(D[5])) {
      stop("Invalid reference category")
    }
    if (D[4] == 1) {
      # Ordered multinomial
      if (is.na(as.integer(D[5]))) {
        D[5] = which(indata[[resp]] == D[5])[1] # Find corresponding numeric factor code
      }
      if (as.integer(D[5]) != min(as.integer(indata[[resp]])) && as.integer(D[5]) != max(as.integer(indata[[resp]]))) {
        stop(paste("Invalid reference category:", D[5]))
      }
    }
  }

  if (drop.data) {
    outvars <- c(resp)

    # Multiple membership IDs if applicable
    if (!is.null(mm)) {
      for (i in 1:length(mm)) {
        if (!any(is.na(mm[[i]]))) {
          varnames <- unlist(mm[[i]]$mmvar)
          weightnames <- unlist(mm[[i]]$weights)
          outvars <- union(outvars, varnames)
          outvars <- union(outvars, weightnames)
        }
      }
    }

    # CAR membership IDs if applicable
    if (!is.null(car)) {
      for (i in 1:length(car)) {
        if (!any(is.na(car[[i]]))) {
          varnames <- unlist(car[[i]]$carvar)
          weightnames <- unlist(car[[i]]$weights)
          outvars <- union(outvars, varnames)
          outvars <- union(outvars, weightnames)
        }
      }
    }

    outvars <- union(outvars, na.omit(levID))
    if (is.list(expl)) {
      if (!is.na(expl$sep.coeff[1])) {
        tsep.coeff <- expl$sep.coeff
      } else {
        tsep.coeff <- NULL
      }
      xvars <- c(tsep.coeff, expl$common.coeff)
    } else {
      xvars <- expl
    }
    outvars <- union(outvars, xvars)
    interpos <- grep("\\:", xvars)
    if (length(interpos) != 0) {
      explx <- xvars[-interpos]
      outvars <- c(outvars, explx)
      # This intersect is a bit of a hack, but avoids variables generated by MLwiN being referenced
      interx <- intersect(unlist(mapply(strsplit, as.character(expl[interpos]), "\\:"), use.names = FALSE), colnames(indata))
      outvars <- union(outvars, interx)
    }

    # Denominators/Offsets if applicable
    if (D[1] == "Binomial" || D[1] == "Poisson" || D[1] == "Multinomial" || D[1] == "Negbinom") {
      if (!is.na(D[[3]])) {
        outvars <- union(outvars, D[[3]])
      }
    }
    if (D[1] == "Mixed") {
      for (i in 2:length(D)) {
        if (D[[i]][1] == "Binomial" || D[[i]][1] == "Poisson")
          if (!is.na(D[[i]][[3]])) {
            outvars <- union(outvars, D[[i]][[3]])
          }
      }
    }

    # (R)IGLS Weights if applicable
    if (!is.null(weighting)) {
      for (w in weighting$weightvar) {
        if (!is.na(w)) {
          outvars <- union(outvars, w)
        }
      }
    }

    outdata <- indata[, outvars]
  } else {
    outdata <- indata
  }

  if (!isTRUE(sort.ignore)) {
    # Don't enforce sorting on level-1 in cases where it isn't used
    if (D[1] == "Normal" || D[1] == "Binomial" || D[1] == "Poisson" || D[1] == "Negbinom") {
      outdata[["_sortindex"]] <- seq(1, nrow(outdata))  # replace with sequence to keep sorting stable
      l1id <- levID[length(levID)]
      levID[length(levID)] <- "_sortindex"
    }

    # Check/sort data as approriate
    if (isTRUE(sort.force)) {
      outdata <- outdata[do.call(order, outdata[na.omit(levID)]), ]
    } else {
      if (is.null(xc) && !isTRUE(all(do.call(order, outdata[na.omit(levID)]) == seq(1, nrow(outdata))))) {
        stop("The input data are not sorted according to the model hierarchy")
      }
    }

    # Restore original level ID and drop temporary variable
    if (D[1] == "Normal" || D[1] == "Binomial" || D[1] == "Poisson" || D[1] == "Negbinom") {
      levID[length(levID)] <- l1id
      outdata[["_sortindex"]] <- NULL
    }
  }

  xcolumns <- NULL
  if (is.list(expl)) {
    if (!is.na(expl$sep.coeff[1])) {
      xcolumns <- c(expl$sep.coeff, expl$common.coeff)
    } else {
      xcolumns <- expl$common.coeff
    }
  } else {
    xcolumns <- expl
  }

  interpos <- grep("\\:", xcolumns)
  if (length(interpos) != 0) {
    explx <- xcolumns[-interpos]
    xcolumns <- union(xcolumns, explx)
    interx <- intersect(unlist(mapply(strsplit, as.character(expl[interpos]), "\\:"), use.names = FALSE), colnames(indata))
    xcolumns <- union(xcolumns, interx)
  }

  # Indicator that at least one response in the row is non-missing
  ymiss <- as.logical(apply(!is.na(outdata[, resp, drop=FALSE]), 1, max))

  # Exclude rows where any X or all responses are missing
  completerows <- complete.cases(outdata[, xcolumns]) & ymiss
  hierarchy <- NULL
  shortID <- na.omit(rev(levID))
  if (length(shortID) > 1) {
    for (lev in length(shortID):2) {
      if (!is.null(xc)) {
        groupsize <- by(outdata, outdata[, shortID[lev]], nrow)
        compgroupsize <- by(outdata[completerows, ], outdata[completerows, shortID[lev]], nrow)
      } else {
        test <- requireNamespace("reshape", quietly = TRUE)
        if (isTRUE(test)) {
          # If the level identifiers are factors with string labels then the following can produce the warning 'coercing
          # argument of type 'list' to logical' from within cbind2 in the reshape package.  This is due to the call:
          # 'all(lapply(list(...), is.numeric))' as the lapply returns a list which all doesn't like.  As the result is
          # still correct a suppressWarnings() call is added below to prevent this being passed onto the user
          groupsize <- as.vector(suppressWarnings(reshape::sparseby(outdata, outdata[, shortID[lev:length(shortID)]],
                                                                    nrow, GROUPNAMES = FALSE)))
          compgroupsize <- as.vector(suppressWarnings(reshape::sparseby(outdata[completerows, ], outdata[completerows, shortID[lev:length(shortID)]],
                                                                    nrow, GROUPNAMES = FALSE)))
        } else {
          groupsize <- na.omit(as.vector(by(outdata, outdata[, shortID[lev:length(shortID)]], nrow)))
          compgroupsize <- na.omit(as.vector(by(outdata[completerows, ], outdata[completerows, shortID[lev:length(shortID)]], nrow)))

        }
      }
      groupinfo <- cbind(length(groupsize), min(groupsize), mean(groupsize), max(groupsize), length(compgroupsize), min(compgroupsize), mean(compgroupsize), max(compgroupsize))
      colnames(groupinfo) <- c("N", "min", "mean", "max", "N_complete", "min_complete", "mean_complete", "max_complete")
      rownames(groupinfo) <- shortID[lev]
      hierarchy <- rbind(hierarchy, groupinfo)
    }
  }

  if (!file.access(workdir) == 0)
    dir.create(workdir)

  dtafile <- normalizePath(tempfile("dtafile_", tmpdir = workdir, fileext = ".dta"), winslash = "/", mustWork = FALSE)

  if (EstM == 0) {
    macrofile <- normalizePath(tempfile("macrofile_", tmpdir = workdir, fileext = ".txt"), winslash = "/", mustWork = FALSE)
    IGLSfile <- normalizePath(tempfile("IGLSfile_", tmpdir = workdir, fileext = ".dta"), winslash = "/", mustWork = FALSE)
  }

  gettempfile <- function(stub, ext) normalizePath(tempfile(stub, tmpdir = workdir, fileext = ext), winslash = "/", mustWork = FALSE)

  if (EstM == 1) {
    macrofile <- replicate(nchains, gettempfile("macrofile_", ".txt"))
    IGLSfile <- replicate(nchains, gettempfile("IGLSfile_", ".dta"))
    MCMCfile <- replicate(nchains, gettempfile("MCMCfile_", ".dta"))
    chainfile  <- replicate(nchains, gettempfile("chainfile_", ".dta"))
    if (!is.null(dami)) {
      MIfile <- replicate(nchains, gettempfile("MIfile_", ".dta"))
    } else {
      dami <- MIfile <- NULL
    }
    if (!is.null(fact)) {
      FACTchainfile <- replicate(nchains, gettempfile("factchainfile_", ".dta"))
    } else {
      FACTchainfile <- NULL
    }
  }
  if (!is.null(BUGO)) {
    modelfile <- replicate(nchains, gettempfile("modelfile_", ".txt"))
    initfile <- replicate(nchains, gettempfile("initfile_", ".txt"))
    datafile <- replicate(nchains, gettempfile("datafile_", ".txt"))
    scriptfile <- normalizePath(tempfile("scriptfile_", tmpdir = workdir, fileext = ".txt"), winslash = "/", mustWork = FALSE)
    bugEst <- normalizePath(tempfile("bugEst_", tmpdir = workdir, fileext = ".txt"), winslash = "/", mustWork = FALSE)
  }
  if (resi.store) {
    resifile <- NULL
    if (EstM == 0) {
      for (i in 1:length(rp)) {
        resifile <- c(resifile, normalizePath(tempfile(paste0("resifile_", i), tmpdir = workdir, fileext = ".dta"),
                                     winslash = "/", mustWork = FALSE))
      }
    } else {
      for (i in 1:length(rp)) {
        resifile <- rbind(resifile, replicate(nchains, gettempfile(paste0("resifile_", i), ".dta")))
      }
    }
  }
  if (!is.null(resi.store.levs))
    resichains  <- replicate(nchains, gettempfile("resichains_", ".dta"))
  if ((D[1] == "Multivariate Normal" || D[1] == "Mixed" || D[1] == "Multinomial") && !is.null(clre)) {
    clre[1, ] <- as.numeric(clre[1, ]) + 1
  }

  dups <- duplicated(tolower(colnames(outdata)))
  if (any(dups)) {
    stop(paste("variables name(s)", paste(colnames(outdata)[dups], collapse=","), "are duplicates when ignoring case"))
  }

  long2shortname <- sapply(colnames(outdata), digest, algo="xxhash64", serialize = FALSE)
  long2shortname[] <- paste0("v", long2shortname)
  colnames(outdata) <- long2shortname
  write.dta(outdata, dtafile, version = 10)
  colnames(outdata) <- names(long2shortname)

  finalClean <- function(clean.files) {
    if (clean.files) {
      file.remove(dtafile)
      file.remove(macrofile)
      file.remove(IGLSfile)
      if (EstM == 1 && is.null(BUGO))
        file.remove(MCMCfile)
      if (EstM == 1 && is.null(BUGO))
        file.remove(chainfile)
      if (!is.null(BUGO))
        file.remove(modelfile)
      if (resi.store && is.null(BUGO)) {
        file.remove(resifile)
      }
      if (EstM == 1 && is.null(BUGO)) {
        if (!is.null(resi.store.levs))
          file.remove(resichains)
        if (!is.null(fact))
          file.remove(FACTchainfile)
        if (!is.null(dami))
          file.remove(MIfile)
      }
    }
  }
  if (EstM == 0) {
      if (length(startval[[1]]) == 0) {
        svals <- NULL
      } else {
        svals <- startval[[1]]
      }
    long2shortnamemap <- write.IGLS(outdata, dtafile, oldsyntax, resp, levID, expl, rp, D, nonlinear, categ, notation, nonfp, clre,
                 Meth, extra, reset, rcon, fcon, maxiter, convtol, mem.init, optimat, weighting, fpsandwich, rpsandwich,
                 macrofile = macrofile, IGLSfile = IGLSfile, resifile = resifile, resi.store = resi.store, resioptions = resioptions,
                 debugmode = debugmode, startval = svals, namemap = long2shortname, saveworksheet = saveworksheet)
    iterations <- estoptions$mcmcMeth$iterations
    if (is.null(iterations))
      iterations <- 5000
    burnin <- estoptions$mcmcMeth$burnin
    if (is.null(burnin))
      burnin <- 500
    thinning <- estoptions$mcmcMeth$thinning
    if (is.null(thinning))
      thinning <- 1


    args <- paste0("/run ", "\"", macrofile, "\"")
    if (!debugmode) {
      args <- paste0("/nogui ", args)
    }
    time1 <- proc.time()
    system2(cmd, args = args, stdout = stdout, stderr = stderr)
    cat("\n")
    time2 <- proc.time() - time1

    estIGLS <- read.dta(IGLSfile)

    FP[] <- na.omit(estIGLS[, 1])

    estIGLS2 <- na.omit(estIGLS[, 2])
    k <- 1
    for (i in 1:length(FP)) {
      for (j in 1:i) {
        FP.cov[i, j] <- estIGLS2[k]
        FP.cov[j, i] <- FP.cov[i, j]
        k <- k + 1
      }
    }

    RP[] <- na.omit(estIGLS[, 3])

    estIGLS4 <- na.omit(estIGLS[, 4])
    k <- 1
    for (i in 1:length(RP)) {
      for (j in 1:i) {
        RP.cov[i, j] <- estIGLS4[k]
        RP.cov[j, i] <- RP.cov[i, j]
        k <- k + 1
      }
    }

    LIKE <- estIGLS[, dim(estIGLS)[2]][3]
    if (!is.na(LIKE)) {
      if (LIKE == 1)
        LIKE <- NA
    }

    NTotal <- estIGLS[, dim(estIGLS)[2]][1]
    NUsed <- estIGLS[, dim(estIGLS)[2]][2]

    if (estIGLS[, dim(estIGLS)[2]][8] == 1) {
      Converged <- TRUE
    } else {
      Converged <- FALSE
    }

    Iterations <- estIGLS[, dim(estIGLS)[2]][7]

    if (resi.store) {
      resiraw <- list()
      for (i in 1:length(rp)) {
        tmp <- as.list(read.dta(resifile[i]))
        for (name in names(long2shortnamemap)) {
          names(tmp) <- gsub(long2shortnamemap[[name]], name, names(tmp))
        }
        for (j in names(tmp)) {
          resiraw[[j]] <- tmp[[j]]
        }
      }
    }

  }

  # MCMC algorithm (using the starting values obtain from IGLS algorithm)
  if (EstM == 1 && is.null(BUGO)) {
    for (i in 1:nchains) {
      if (length(startval[[i]]) == 0) {
        svals <- NULL
      } else {
        svals <- startval[[i]]
      }
      long2shortnamemap <- write.MCMC(outdata, dtafile, oldsyntax, resp, levID, expl, rp, D, nonlinear, categ, notation, nonfp, clre,
                   Meth, merr, carcentre, maxiter, convtol, seed[i], iterations, burnin, scale, thinning, priorParam, refresh,
                   fixM, residM, Lev1VarM, OtherVarM, adaption, priorcode, rate, tol, lclo, mcmcOptions, fact, xc, mm, car,
                   BUGO, mem.init, optimat, modelfile = NULL, initfile = NULL, datafile = NULL, macrofile = macrofile[i],
                   IGLSfile = IGLSfile[i], MCMCfile = MCMCfile[i], chainfile = chainfile[i], MIfile = MIfile[i], resifile = resifile[,i],
                   resi.store = resi.store, resioptions = resioptions, resichains = resichains[i], FACTchainfile = FACTchainfile[i],
                   resi.store.levs = resi.store.levs, debugmode = debugmode, startval = svals, dami = dami,
                   namemap = long2shortname, saveworksheet = saveworksheet[i])

    }
    cat("MLwiN is running, please wait......\n")
    time1 <- proc.time()
    clust <- NULL
    regpar <- FALSE
    # Check whether the user has registered their own backend, if not use "doParallel"
    if (!getDoParRegistered()) {
      if (nchains > 1) {
        clust <- makeCluster(min(nchains, detectCores(logical = FALSE)), outfile=stdout)
        registerDoParallel(clust)
      } else {
        registerDoSEQ()
      }
      regpar <- TRUE
    }
    foreach(i=1:nchains) %dopar% {
      args <- paste0("/run ", "\"", macrofile[i], "\"")
      if (!debugmode) {
        args <- paste0("/nogui ", args)
      }
      system2(cmd, args = args, stdout = stdout, stderr = stderr)
    }
    if (isTRUE(regpar)) {
      if (!is.null(clust)) {
        stopCluster(clust)
      }
      # unregister doParallel as detailed on http://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster
      env <- get(".foreachGlobals", asNamespace("foreach"), inherits = FALSE)
      rm(list=ls(name=env), pos=env)
    }

    cat("\n")
    time2 <- proc.time() - time1

    chainslist <- list()
    factchainslist <- list()
    resichainslist <- list()
    factloadchainslist <- list()
    factcovchainslist <- list()
    factscorelist <- list()
    factscorevarlist <- list()
    resilist <- list()
    MIlist <- list()
    BDIC <- rep(0, 4)
    BDIC.names <- c("Dbar", "D(thetabar)", "pD", "DIC")
    names(BDIC) <- BDIC.names
    LIKE <- 0

    for (i in 1:nchains) {
      nlev <- length(levID)
      chains <- read.dta(chainfile[i])
      for (name in names(long2shortnamemap)) {
        colnames(chains) <- gsub(long2shortnamemap[[name]], name, colnames(chains))
      }

      chains <- coda::mcmc(data = chains[, -1], thin = thinning)
      chain.names <- colnames(chains)
      chain.names[grep("RP", chain.names)] <- RP.names
      colnames(chains) <- chain.names

      if (sum(grepl("bcons", colnames(chains))) > 0) {
        bcons.pos <- grep("bcons", colnames(chains))
        chains[1, bcons.pos] <- chains[1, bcons.pos] - 0.001
      }

      chainslist[[i]] <- chains

      estMCMC <- read.dta(MCMCfile[i])

      if (!(D[1] == "Mixed") && is.null(merr) && is.null(fact)) {
        BDIC <- BDIC + estMCMC[, dim(estMCMC)[2]][c(5, 6, 4, 3)]
      } else {
        LIKE <- LIKE + estMCMC[, dim(estMCMC)[2]][3]
      }

      NTotal <- estMCMC[, dim(estMCMC)[2]][1]
      NUsed <- estMCMC[, dim(estMCMC)[2]][2]
      if (!is.null(fact)) {
        load.names <- rep("", length(loadings))
        k <- 1
        for (j1 in 1:fact$nfact) {
          for (j2 in resp) {
            load.names[k] <- paste("load", j1, "_", j2, sep = "")
            k <- k + 1
          }
        }

        fact.cov.names <- rep("", (fact$nfact*(fact$nfact + 1))/2)
        k <- 1
        for (j1 in 1:fact$nfact) {
          for (j2 in 1:j1) {
            if (j1 == j2) {
              fact.cov.names[k] <- paste("var_fact", j1, sep = "")
            } else {
              fact.cov.names[k] <- paste("cov_fact", j1, "_fact", j2, sep = "")
            }
            k <- k + 1
          }
        }

        factchains <- read.dta(FACTchainfile[i])
        factscores <- matrix(na.omit(factchains[, "_FACT_value_b"]), ncol = fact$nfact, byrow = FALSE)
        factscores_v <- matrix(na.omit(factchains[, "_FACT_value_v"]), ncol = fact$nfact, byrow = FALSE)
        factloads <- matrix(na.omit(factchains[, "_FACT_load_b_chain"]), nrow = iterations/thinning, byrow = TRUE)
        factcovs <- matrix(na.omit(factchains[, "_FACT_load_v_chain"]), nrow = iterations/thinning, byrow = TRUE)
        nameloads <- NULL
        namefacts <- NULL
        namefacts_v <- NULL
        namecovs <- NULL
        for (j in 1:fact$nfact) {
          if (fact$lev.fact[j] > 1) {
            nunit <- nrow(unique(indata[rev(na.omit(levID))[fact$lev.fact[j]]]))
            if (length(factscores) > nunit) {
              factscores[(nunit + 1):nrow(factscores), j] <- NA
            }
          }
          namefacts <- c(namefacts, paste0("factorscores", j))
          namefacts_v <- c(namefacts_v, paste0("factorscores_var", j))
        }
        colnames(factscores) <- namefacts
        colnames(factscores_v) <- namefacts_v
        colnames(factloads) <- load.names
        colnames(factcovs) <- fact.cov.names
        factChains <- list(scores = factscores, scores_v = factscores_v, loadings = coda::mcmc(data = factloads,
                                                                                         thin = thinning), cov = coda::mcmc(data = factcovs, thin = thinning))
        factloadchainslist[[i]] <- factChains$loadings
        factcovchainslist[[i]] <- factChains$cov
        factscorelist[[i]] <- factChains$scores
        factscorevarlist[[i]] <- factChains$scores_v
      }

      if (!is.null(dami)) {
        MIdata <- read.dta(MIfile[i])
        MIlist[[i]] <- MIdata
      }

      if (!is.null(resi.store.levs)) {
        residata <- read.dta(resichains[i])
        for (name in names(long2shortnamemap)) {
          colnames(residata) <- gsub(long2shortnamemap[[name]], name, colnames(residata))
        }
        resiChains <- list()
        for (name in colnames(residata)) {
          lev <- as.integer(gsub("resi_lev", "", name))
          ucount <- length(rp[[paste0("rp", lev)]])
          if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
            lev = lev + 1
          }
          nunit <- nrow(unique(indata[rev(levID)[lev]]))
          pnames <- paste("u", (1:ucount)-1, rep(1:nunit, each=ucount), sep="_")
          resiChains[[name]] <- coda::mcmc(data = matrix(na.omit(residata[, name]), nrow = iterations/thinning, byrow = TRUE,
                                           dimnames = list(1:(iterations/thinning), pnames)), thin = thinning)
        }
        resichainslist[[i]] <- resiChains
      }
      if (resi.store) {
        resiraw <- list()
        for (j in 1:length(rp)) {
          tmp <- as.list(read.dta(resifile[j, i]))
          for (name in names(long2shortnamemap)) {
            names(tmp) <- gsub(long2shortnamemap[[name]], name, names(tmp))
          }
          for (k in names(tmp)) {
            resiraw[[k]] <- tmp[[k]]
          }
          resilist[[i]] <- resiraw
        }
      }
    }

    if (nchains != 1) {
      chains <- mcmc.list(chainslist)
      if (!is.null(resi.store.levs)) {
        resiChains <- mcmc.list(resichainslist)
      }
      if (!is.null(fact)) {
        factChains$loadings <- mcmc.list(factloadchainslist)
        factChains$cov <- mcmc.list(factcovchainslist)
      }
    }

    if (!is.null(dami)) {
      if (dami[1] == 0) {
        imputations <- list()
        for (i in 1:nchains) {
          for (j in 2:length(dami)) {
            impdata <- MIlist[[i]][, c("resp_indicator", paste0("_est_", dami[j]))]
            Nresp <- length(unique(impdata$resp_indicator))
            Nrecs = nrow(impdata) / Nresp
            impdata$id <- rep(1:Nrecs, each=Nresp)
            impdata <- reshape(impdata, timevar="resp_indicator", idvar="id", direction="wide")
            impdata$id <- NULL
            colnames(impdata) <- gsub(paste0("_est_", dami[j], "."), "", colnames(impdata))
            impout <- outdata
            impout[, colnames(impdata)] <- impdata
            imputations[[paste0("chain", i, "_iteration_", dami[j])]] <- impout
          }
        }
      } else {
        MIdata <- list()
        if ("_MissingInd" %in% colnames(MIlist[[1]])) {
          MIdata[["_MissingInd"]] <- MIlist[[1]][, "_MissingInd"]
        }
        for (i in 1:nchains) {
          impcols <- c(grep("_est_", colnames(MIlist[[i]])), grep("_SDs_", colnames(MIlist[[i]])))
          colnames(MIlist[[i]])[impcols] <- paste0("chain", i, colnames(MIlist[[i]])[impcols])
          MIdata <- c(MIdata, MIlist[[i]][, impcols])
        }
        MIdata <- data.frame(MIdata)
      }
    }

    combchains <- as.matrix(chains)
    FP[FP.names] <- colMeans(combchains[, FP.names, drop=FALSE])
    FP.cov[FP.names, FP.names] <- cov(combchains[, FP.names, drop=FALSE])
    RP[RP.names] <- colMeans(combchains[, RP.names, drop=FALSE])
    RP.cov[RP.names, RP.names] <- cov(combchains[, RP.names, drop=FALSE])
    ESS <- effectiveSize(chains)
    BDIC <- BDIC / nchains
    LIKE <- LIKE / nchains
    if (!is.na(LIKE)) {
      if (LIKE == 1) {
        LIKE <- NA
      }
    }

    if (resi.store) {
      resinames <- names(resilist[[1]])
      resiraw <- list()
      resi.means <- resinames[grep("_resi_est", resinames)]
      resi.std <- resinames[grep("_resi_est", resinames)]
      resi.var <- resinames[grep("_resi_variance", resinames)]
      resi.se <- resinames[grep("_resi_se", resinames)]
      for (name in resi.means) {
        for (i in 1:nchains) {
          if (is.null(resiraw[[name]])) {
            resiraw[[name]] <- resilist[[i]][[name]]
          } else {
            resiraw[[name]] <- resiraw[[name]] + resilist[[i]][[name]]
          }
        }
        resiraw[[name]] <- resiraw[[name]] / nchains
      }
      for (name in resi.std) {
        for (i in 1:nchains) {
          if (is.null(resiraw[[name]])) {
            resiraw[[name]] <- resilist[[i]][[name]]
          } else {
            resiraw[[name]] <- resiraw[[name]] + resilist[[i]][[name]]
          }
        }
        resiraw[[name]] <- resiraw[[name]] / nchains
      }
      for (name in resi.var) {
        meanname <- sub("_resi_variance_", "_resi_est_", name)
        for (i in 1:nchains) {
          if (is.null(resiraw[[name]])) {
            resiraw[[name]] <- resilist[[i]][[name]]
          } else {
            resiraw[[name]] <- resiraw[[name]] + resilist[[i]][[name]]
          }
        }
        resiraw[[name]] <- resiraw[[name]] / nchains
        if (nchains > 1) {
          B <- rep(0, length(resiraw[[name]]))
          for (i in 1:nchains) {
            B <- B + ((resilist[[i]][[meanname]] - resiraw[[meanname]])^2 / (nchains - 1))
          }
          resiraw[[name]] <- resiraw[[name]] + (1.0 + (1.0 / nchains)) * B
        }
      }
      for (name in resi.se) {
        meanname <- sub("_resi_se_", "_resi_est_", name)
        for (i in 1:nchains) {
          if (is.null(resiraw[[name]])) {
            resiraw[[name]] <- resilist[[i]][[name]]^2
          } else {
            resiraw[[name]] <- resiraw[[name]] + resilist[[i]][[name]]^2
          }
        }
        resiraw[[name]] <- resiraw[[name]] / nchains
        if (nchains > 1) {
          B <- rep(0, length(resiraw[[name]]))
          for (i in 1:nchains) {
            B <- B + ((resilist[[i]][[meanname]] - resiraw[[meanname]])^2 / (nchains - 1))
          }
          resiraw[[name]] <- resiraw[[name]] + (1.0 + (1.0 / nchains)) * B
        }
        resiraw[[name]] <- sqrt(resiraw[[name]])
      }
    }

    if (!is.null(fact)) {
      for (j in 1:fact$nfact) {
        factname <- paste0("factorscores", j)
        factvarname <- paste0("factorscores_var", j)
        factChains$scores[,factname] <- factscorelist[[1]][,factname]
        factChains$scores_v[,factvarname] <- factscorevarlist[[1]][,factvarname]
        if (nchains > 1) {
          for (i in 2:nchains) {
            factChains$scores[,factname] <- factChains$scores[,factname] + factscorelist[[i]][,factname]
            factChains$scores_v[,factvarname] <- factChains$scores_v[,factvarname] + factscorevarlist[[i]][,factvarname]
          }
          factChains$scores[,factname] <- factChains$scores[,factname] / nchains
          factChains$scores_v[,factvarname] <- factChains$scores_v[,factvarname] / nchains
          if (nchains > 1) {
            B <- rep(0, length(factChains$scores[,factname]))
            for (i in 1:nchains) {
              B <- B + ((factscorelist[[i]][,factname] - factChains$scores[,factname])^2 / (nchains - 1))
            }
            factChains$scores_v[,factvarname] <- factChains$scores_v[,factvarname] + (1.0 + (1.0 / nchains)) * B
          }
        }
      }
      loadings <- colMeans(as.matrix(factChains$loadings))
      loadings.sd <- sqrt(diag(cov(as.matrix(factChains$loadings))))
      fact.cov <- colMeans(as.matrix(factChains$cov))
      fact.cov.sd <- sqrt(diag(cov(as.matrix(factChains$cov))))
    }
  }

  if (EstM == 1 && !is.null(BUGO)) {
    for (i in 1:nchains) {
      if (length(startval[[i]]) == 0) {
        svals <- NULL
      } else {
        svals <- startval[[i]]
      }
      long2shortname <- write.MCMC(outdata, dtafile, oldsyntax, resp, levID, expl, rp, D, nonlinear, categ, notation, nonfp, clre,
                   Meth, merr, carcentre, maxiter, convtol, seed, iterations, burnin, scale, thinning, priorParam, refresh,
                   fixM, residM, Lev1VarM, OtherVarM, adaption, priorcode, rate, tol, lclo, mcmcOptions, fact, xc, mm, car,
                   BUGO, mem.init, optimat, modelfile = modelfile[i], initfile = initfile[i], datafile = datafile[i], macrofile = macrofile[i],
                   IGLSfile = IGLSfile[i], MCMCfile = MCMCfile[i], chainfile = chainfile[i], MIfile = MIfile[i], resifile = resifile[i],
                   resi.store = resi.store, resioptions = resioptions, resichains = resichains, FACTchainfile = FACTchainfile[i],
                   resi.store.levs = resi.store.levs, debugmode = debugmode, startval = svals, dami = dami,
                   namemap = long2shortname)
    }

    cat("MLwiN is running, please wait......\n")
    time1 <- proc.time()
    for (i in 1:nchains) {
      args <- paste0("/run ", "\"", macrofile[i], "\"")
      if (!debugmode) {
        args <- paste0("/nogui ", args)
      }
      system2(cmd, args = args, stdout = stdout, stderr = stderr)
    }
    cat("\n")
    time2 <- proc.time() - time1
    nlev <- length(levID)
  }

  if (show.file)
    file.show(macrofile)
  if ((!is.null(BUGO)) && !(D[1] == "Mixed")) {
    if (show.file)
      file.show(modelfile[1])
    n.iter <- iterations + burnin
    addmore <- NULL
    if (EstM == 1) {
      if (!is.null(car)) {
        addmore <- c(addmore, "carmean")
      }
      if (is.matrix(mcmcOptions$paex)) {
        if (sum(mcmcOptions$paex[, 2]) > 0) {
          vx <- sapply(2:nlev, function(x) paste("v", x, sep = ""))
          sigma2v <- sapply(2:nlev, function(x) paste("sigma2.v", x, sep = ""))
          addmore <- c(addmore, vx, sigma2v)
        }
      } else {
        if (is.vector(mcmcOptions$paex) && length(mcmcOptions$paex) == 2) {
          if (mcmcOptions$paex[2] > 0) {
            vx <- sapply(2:nlev, function(x) paste("v", x, sep = ""))
            sigma2v <- sapply(2:nlev, function(x) paste("sigma2.v", x, sep = ""))
            addmore <- c(addmore, vx, sigma2v)
          }
        }
      }
    }

    debug <- as.logical(BUGO["debug"])
    if (is.na(debug)) {
      debug <- FALSE
    }
    OpenBugs <- as.logical(BUGO["OpenBugs"])
    if (is.na(OpenBugs)) {
      OpenBugs <- FALSE
    }
    n.chains <- as.integer(BUGO["n.chains"])
    if (is.na(n.chains)) {
      n.chains <- 1
    }
    bugs.seed <- BUGO["seed"]
    bugs <- BUGO["bugs"]
    if (is.na(bugs.seed)) {
      bugs.seed <- NULL
    }
    if (is.na(bugs)) {
      stop("Need to specify path to the BUGS executable.")
    }
    chains.bugs.mcmc <- mlwin2bugs(D, levID, datafile[1], initfile, modelfile[1], bugEst, fact, addmore, n.chains = n.chains,
                                   n.iter = n.iter, n.burnin = burnin, n.thin = thinning, debug = debug, bugs = bugs, bugsWorkingDir = workdir,
                                   OpenBugs = OpenBugs, cleanBugsWorkingDir = clean.files, seed = bugs.seed)
    time2 <- proc.time() - time1
  } else {
    if (D[1] == "Mixed" && (!is.null(BUGO)))
      warning("The Mixed response model is currently not implemented in WinBUGS/OpenBUGS.")
  }

  if (EstM == 0) {
    if (is.null(BUGO)) {
      outIGLS <- new("mlwinfitIGLS")
      outIGLS["version"] <- versiontext
      outIGLS["Nobs"] <- NUsed
      outIGLS["DataLength"] <- NTotal
      outIGLS["Hierarchy"] <- hierarchy
      outIGLS["D"] <- D
      outIGLS["Formula"] <- Formula
      outIGLS["levID"] <- levID
      outIGLS["FP"] <- FP
      outIGLS["RP"] <- RP
      outIGLS["FP.cov"] <- FP.cov
      outIGLS["RP.cov"] <- RP.cov
      outIGLS["LIKE"] <- LIKE
      outIGLS["Meth"] <- Meth
      outIGLS["nonlinear"] <- nonlinear
      outIGLS["Converged"] <- Converged
      outIGLS["Iterations"] <- Iterations
      outIGLS["elapsed.time"] <- time2[3]
      outIGLS["call"] <- cl
      outIGLS["data"] <- as.data.frame(outdata)

      if (resi.store) {
        outIGLS["residual"] <- resiraw
      }

      finalClean(clean.files)
      return(outIGLS)
    } else {
      finalClean(clean.files)
      print(summary(chains.bugs.mcmc))
      return(chains.bugs.mcmc)
    }
  }
  if (EstM == 1) {
    if (is.null(BUGO)) {
      outMCMC <- new("mlwinfitMCMC")
      outMCMC["version"] <- versiontext
      outMCMC["Nobs"] <- NUsed
      outMCMC["DataLength"] <- NTotal
      outMCMC["Hierarchy"] <- hierarchy
      outMCMC["burnin"] <- burnin
      outMCMC["iterations"] <- iterations
      outMCMC["nchains"] <- nchains
      outMCMC["D"] <- D
      outMCMC["Formula"] <- Formula
      outMCMC["levID"] <- levID
      outMCMC["merr"] <- merr
      outMCMC["fact"] <- fact
      outMCMC["xc"] <- xc
      outMCMC["FP"] <- FP
      outMCMC["RP"] <- RP
      outMCMC["FP.cov"] <- FP.cov
      outMCMC["RP.cov"] <- RP.cov
      outMCMC["chains"] <- chains
      outMCMC["elapsed.time"] <- time2[3]
      outMCMC["call"] <- cl
      outMCMC["data"] <- as.data.frame(outdata)
      if (!(D[1] == "Mixed") && is.null(merr) && is.null(fact)) {
        outMCMC["BDIC"] <- BDIC
      } else {
        outMCMC["LIKE"] <- LIKE
      }
      if (!is.null(fact)) {
        outMCMC["fact.loadings"] <- loadings
        outMCMC["fact.loadings.sd"] <- loadings.sd
        outMCMC["fact.cov"] <- fact.cov
        outMCMC["fact.cov.sd"] <- fact.cov.sd
        outMCMC["fact.chains"] <- factChains
      }

      if (!is.null(resi.store.levs)) {
        outMCMC["resi.chains"] <- resiChains
      }
      if (!is.null(dami)) {
        if (dami[1] == 0) {
          outMCMC["imputations"] <- imputations
        } else {
          outMCMC["MIdata"] <- MIdata
        }
      }

      if (resi.store) {
        outMCMC["residual"] <- resiraw
      }
      finalClean(clean.files)
      return(outMCMC)
    } else {
      finalClean(clean.files)
      return(chains.bugs.mcmc)
    }
  }
}

##' @S3method summary mlwinfitIGLS
summary.mlwinfitIGLS <- function(object, ...) {
    summary(object)
}

##' @S3method print mlwinfitIGLS
print.mlwinfitIGLS <- function(x, ...) {
    print(x)
}

##' @S3method show mlwinfitIGLS
show.mlwinfitIGLS <- function(object, ...) {
    show(object)
}

##' @importFrom stats update
##' @S3method update mlwinfitIGLS
update.mlwinfitIGLS <- function(object, ...) {
    update(object)
}

##' @importFrom stats coef
##' @S3method coef mlwinfitIGLS
coef.mlwinfitIGLS <- function(object, ...) {
    coef(object)
}

##' @importFrom stats coefficients
##' @S3method coefficients mlwinfitIGLS
coefficients.mlwinfitIGLS <- function(object, ...) {
    coefficients(object)
}

##' @importFrom stats vcov
##' @S3method vcov mlwinfitIGLS
vcov.mlwinfitIGLS <- function(object, ...) {
    vcov(object)
}

##' @importFrom stats df.residual
##' @S3method df.residual mlwinfitIGLS
df.residual.mlwinfitIGLS <- function(object, ...) {
    df.residual(object)
}

##' @importFrom stats fitted
##' @S3method fitted mlwinfitIGLS
fitted.mlwinfitIGLS <- function(object, ...) {
    fitted(object)
}

##' @importFrom stats fitted.values
##' @S3method fitted.values mlwinfitIGLS
fitted.values.mlwinfitIGLS <- function(object, ...) {
    fitted.values(object)
}

##' @importFrom stats residuals
##' @S3method residuals mlwinfitIGLS
residuals.mlwinfitIGLS <- function(object, ...) {
    residuals(object)
}

##' @importFrom stats resid
##' @S3method resid mlwinfitIGLS
resid.mlwinfitIGLS <- function(object, ...) {
    resid(object)
}

##' @importFrom stats predict
##' @S3method predict mlwinfitIGLS
predict.mlwinfitIGLS <- function(object, ...) {
    predict(object)
}

##' @importFrom stats logLik
##' @S3method logLik mlwinfitIGLS
logLik.mlwinfitIGLS <- function(object, ...) {
    logLik(object)
}

##' @importFrom stats deviance
##' @S3method deviance mlwinfitIGLS
deviance.mlwinfitIGLS <- function(object, ...) {
    deviance(object)
}

##' @importFrom stats nobs
##' @S3method nobs mlwinfitIGLS
nobs.mlwinfitIGLS <- function(object, ...) {
    object@Nobs
}

##' @S3method summary mlwinfitMCMC
summary.mlwinfitMCMC <- function(object, ...) {
    summary(object)
}

##' @S3method print mlwinfitMCMC
print.mlwinfitMCMC <- function(x, ...) {
    print(x)
}

##' @S3method show mlwinfitMCMC
show.mlwinfitMCMC <- function(object, ...) {
    show(object)
}

##' @importFrom stats update
##' @S3method update mlwinfitMCMC
update.mlwinfitMCMC <- function(object, ...) {
    update(object)
}

##' @importFrom stats coef
##' @S3method coef mlwinfitMCMC
coef.mlwinfitMCMC <- function(object, ...) {
    coef(object)
}

##' @importFrom stats coefficients
##' @S3method coefficients mlwinfitMCMC
coefficients.mlwinfitMCMC <- function(object, ...) {
    coefficients(object)
}

##' @importFrom stats vcov
##' @S3method vcov mlwinfitMCMC
vcov.mlwinfitMCMC <- function(object, ...) {
    vcov(object)
}

##' @importFrom stats fitted
##' @S3method fitted mlwinfitMCMC
fitted.mlwinfitMCMC <- function(object, ...) {
    fitted(object)
}

##' @importFrom stats fitted.values
##' @S3method fitted.values mlwinfitMCMC
fitted.values.mlwinfitMCMC <- function(object, ...) {
    fitted.values(object)
}

##' @importFrom stats residuals
##' @S3method residuals mlwinfitMCMC
residuals.mlwinfitMCMC <- function(object, ...) {
    residuals(object)
}

##' @importFrom stats resid
##' @S3method resid mlwinfitMCMC
resid.mlwinfitMCMC <- function(object, ...) {
    resid(object)
}

##' @importFrom stats predict
##' @S3method predict mlwinfitMCMC
predict.mlwinfitMCMC <- function(object, ...) {
    predict(object)
}

##' @importFrom stats nobs
##' @S3method nobs mlwinfitMCMC
nobs.mlwinfitMCMC <- function(object, ...) {
    object@Nobs
}
