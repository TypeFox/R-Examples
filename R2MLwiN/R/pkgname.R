#' Running MLwiN from within R
#' 
#' \pkg{R2MLwiN} is an R command interface to the MLwiN multilevel modelling
#' software package, allowing users to fit multilevel models using
#' MLwiN (and also WinBUGS / OpenBUGS) from within the R environment.
#'
#' @importFrom rbugs runBugs rbugs2coda genBugsScript getBugsOutput rbugs
#' @importFrom coda mcmc mcmc.list as.mcmc.list read.coda effectiveSize raftery.diag thin is.mcmc
#' @importFrom lattice histogram densityplot xyplot qqmath Rows trellis.par.get trellis.par.set panel.xyplot panel.grid panel.abline panel.segments
#' @importFrom foreign read.dta write.dta
#' @importFrom digest digest
#' @importFrom methods show
#' @importFrom stats4 coef logLik summary vcov update
#' @importFrom stats formula
#' @importFrom Matrix nnzero sparseMatrix
#' 
#' @section Important differences between version 0.8-0 and earlier versions:
#' A number of wide-ranging changes, including a new model-fitting syntax more
#' in keeping with that conventionally used in R, were introduced in
#' \pkg{R2MLwiN} version 0.8-0.
#' 
#' The demos, which replicate both the User's Guide to MLwiN (Rasbash et al, 2012) and
#' MCMC Estimation in MLwiN (Browne, 2012) manuals, provide practical demonstrations of many
#' of these changes. See \code{demo(package = "R2MLwiN")} for a list of demo titles; to run one
#' type e.g. \code{demo(UserGuide03)} or view a demo's script via
#' \code{file.show(system.file("demo", "UserGuide03", package = "R2MLwiN"))}.
#' 
#' \itemize{
#'
#'  \item{The Formula is now specified via a \code{\link[stats]{formula}} object (with some differences
#'  in specification: see \code{\link{runMLwiN}}). So, for example,
#'  previously a 2-level model random intercept model would be specified by e.g.
#'  \code{normexam ~ (0|cons + standlrt) + (2|cons) + (1|cons), levID = c('school', 'student')},
#'  with \code{normexam} the response variable, \code{cons} a constant of ones forming the intercept,
#'  which is allowed to vary at level 1 (\code{student}) and level 2 (\code{school}), and
#'  \code{standlrt} included as a predictor in the fixed part of the model. Whilst back-compatibility
#'  is preserved (i.e. this specification will currently still work) the same model can now be more
#'  parsimoniously specified via \code{normexam ~ 1 + standlrt + (1 | school) + (1 | student)}.
#'  As well examples in the demos, see \code{\link{runMLwiN}} and \code{\link{Formula.translate}} for further info.}
#'  
#'  \item{As a means of specifying cross-classified, multiple membership or CAR models, \code{xclass} is now deprecated.
#'  Instead, cross-classified models are specified via \code{xc = TRUE}, multiple
#'  membership models are specified via \code{mm}, and CAR models are specified via \code{car},
#'  in the list of \code{estoptions}. \code{mm} and \code{car} can be a list of
#'  variable names, a list of vectors, or a matrix. See \code{\link{runMLwiN}} for further details.}
#'  
#'  \item{Multiple membership/CAR information can now be specified using matrices.
#'  \code{\link{df2matrix}} and \code{\link{matrix2df}} functions have also been added to convert such
#'  information between \code{\link[base]{data.frame}} and \code{\link[base]{matrix}} formats.}
#'  
#'  \item{As a means of specifying common (i.e. the same for each category) or separate (i.e. one for each
#'  category) coefficients in ordered multinomial
#'  and multivariate response models, \code{c} (for common) and \code{s} (for separate) have been
#'  replaced by the employment of square brackets after the relevant variable to indicate a common
#'  coefficient is to be fitted (a separate coefficient will be fitted otherwise). Within these square
#'  brackets needs to be placed a numeric identifier indicating the responses for which a common coefficient
#'  is to be added (see \code{\link{runMLwiN}} for further details). E.g. what would have been previously
#'  specified, within the \code{Formula} object, as \code{... (0s|cons + ravens) + (0c|fluent{1, 0}) ...} would now be specified by
#'  \code{... 1 + ravens + fluent[1] ...}.}
#'  
#'  \item{When added as a predictor, a variable encoded as a \code{\link[base]{factor}} is automatically
#'  handled as categorical, replacing the previous use of square brackets after the variable name.}
#'  
#'  \item{A number of generic s4 methods have been added to improve compatibility with statistical functions
#'  which use them (e.g. see \code{\link[stats4]{stats4-package}}). So, for example, the addition of 
#'  a \code{\link{logLik}} means a likelihood ratio test can now be conducted
#'  on two \code{\link{mlwinfitIGLS-class}} objects using the \code{\link[lmtest]{lrtest}} function,
#'  e.g. \code{lrtest(mymodel1, mymodel2)}. See \code{help(package = "R2MLwiN")} for
#'  the index listing these various methods.}
#'  }
#'
#' @examples
#' \dontrun{
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' data(tutorial, package = "R2MLwiN")
#' 
#' (mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student),
#'                      estoptions = list(EstM = 1), data = tutorial))
#' 
#' ## The R2MLwiN package includes scripts to replicate all the analyses in
#' ## Rasbash et al (2012) A User's Guide to MLwiN Version 2.26 and
#' ## Browne, W.J. (2012) MCMC estimation in MLwiN Version 2.26.
#' ## The MLwiN manuals are available online, see:
#' ## http://www.bristol.ac.uk/cmm/software/mlwin/download/manuals.html
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
#' @section References:
#' 
#' \subsection{MLwiN software and manuals}{
#' Browne, W.J. (2012) MCMC Estimation in MLwiN, v2.26.
#' Centre for Multilevel Modelling, University of Bristol.
#'
#' Rasbash, J., Charlton, C., Browne, W.J., Healy, M. and Cameron, B. (2009)
#' MLwiN Version 2.1. Centre for Multilevel Modelling, University of Bristol.
#' 
#' Rasbash, J., Charlton, C. and Pillinger, R. (2012) Manual Supplement to
#' MLwiN v2.26. Centre for Multilevel Modelling, University of Bristol.
#' 
#' Rasbash, J., Steele, F., Browne, W.J. and Goldstein, H. (2012)
#' A User's Guide to MLwiN Version 2.26. Centre for Multilevel Modelling,
#' University of Bristol.
#' }
#'
#' \subsection{OpenBUGS}{
#' Thomas, A., O'Hara, B., Ligges, U. and Sturtz, S. (2006) Making BUGS Open.
#' R News, 6, 12:17.
#' }
#'
#' \subsection{WinBUGS}{
#' Spiegelhalter, D.J., Thomas, A. and Best, N.G. (1999) WinBUGS Version 1.2
#' User Manual. MRC Biostatistics Unit.
#' }
#'
#' @section Maintainer:
#' Zhengzheng Zhang \email{zhengzheng236@@gmail.com}
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#' 
#' @docType package
#' @name R2MLwiN
NULL
 
