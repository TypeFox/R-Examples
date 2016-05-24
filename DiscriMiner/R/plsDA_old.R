#' PLS Discriminant Analysis
#' 
#' Performs a Partial Least Squares (PLS) Discriminant Analysis
#' 
#' When \code{validation=NULL} leave-one-out (loo) cross-validation is
#' performed. \cr When \code{validation="learntest"} validation is performed by
#' providing a learn-set and a test-set of observations. \cr
#' 
#' @param variables matrix or data frame with explanatory variables
#' @param group vector or factor with group memberships
#' @param autosel logical indicating automatic selection of PLS components by
#' cross-validation. Default \code{autosel=TRUE}
#' @param comps integer greater than one indicating the number of PLS
#' components to retain. Used only when \code{autosel=FALSE}
#' @param validation type of validation, either \code{NULL} or
#' \code{"learntest"}. Default \code{NULL}
#' @param learn optional vector of indices for a learn-set. Only used when
#' \code{validation="learntest"}. Default \code{NULL}
#' @param test optional vector of indices for a test-set. Only used when
#' \code{validation="learntest"}. Default \code{NULL}
#' @param scaled logical indicating whether to scale the data (default
#' \code{TRUE})
#' @return An object of class \code{"plsda"}, basically a list with the
#' following elements:
#' @return \item{functions}{table with discriminant functions}
#' @return \item{confusion}{confusion matrix}
#' @return \item{scores}{discriminant scores for each observation}
#' @return \item{classification}{assigned class}
#' @return \item{error_rate}{misclassification error rate}
#' @return \item{components}{PLS components}
#' @return \item{Q2}{quality of loo cross-validation}
#' @return \item{R2}{R-squared coefficients}
#' @return \item{VIP}{Variable Importance for Projection}
#' @return \item{comp_vars}{correlations between components and variables}
#' @return \item{comp_group}{correlations between components and groups}
#' @author Gaston Sanchez
#' @seealso \code{\link{classify}}, \code{\link{geoDA}}, \code{\link{linDA}},
#' \code{\link{quaDA}}
#' @references Tenenhaus M. (1998) \emph{La Regression PLS}. Editions Technip,
#' Paris.
#' 
#' Perez-Enciso M., Tenenhaus M. (2003) \emph{Prediction of clinical outcome
#' with microarray data: a partial least squares discriminant analysis (PLS-DA)
#' approach}. Human Genetics 112: 581-592.
#' @keywords internal
#' @note This is a previous version of plsDA. Not used anymore.
#' @examples
#' 
#'   \dontrun{
#'   # load iris dataset
#'   data(iris)
#' 
#'   # PLS discriminant analysis specifying number of components = 2
#'   my_pls1 = plsDA(iris[,1:4], iris$Species, autosel=FALSE, comps=2)
#'   my_pls1$confusion
#'   my_pls1$error_rate
#'   # plot circle of correlations
#'   plot(my_pls1)
#' 
#'   # PLS discriminant analysis with automatic selection of components
#'   my_pls2 = plsDA(iris[,1:4], iris$Species, autosel=TRUE)
#'   my_pls2$confusion
#'   my_pls2$error_rate
#'   
#'   # linear discriminant analysis with learn-test validation
#'   learning = c(1:40, 51:90, 101:140)
#'   testing = c(41:50, 91:100, 141:150)
#'   my_pls3 = plsDA(iris[,1:4], iris$Species, validation="learntest", learn=learning, test=testing)
#'   my_pls3$confusion
#'   my_pls3$error_rate
#'   }
#' 
plsDA_old <- 
function(variables, group, autosel = TRUE, comps = 2,
         validation = NULL, learn = NULL, test = NULL, scaled=TRUE)
{
  # Perform a PLS discriminant analysis
  # variables: matrix or data.frame with explanatory variables
  # group: vector or factor with group membership
  # autosel: logical indicating automatic selection of PLS comps
  # comps: number of PLS components (only when autosel=FALSE)
  # validation: NULL, "crossval", "learntest"
  # learn: vector of learn-set
  # test: vector of test-set
  # scaled: logical indicating whether to scale the data
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  X = verify_Xy$X
  y = verify_Xy$y
  # autosel
  if (!is.logical(autosel)) autosel = TRUE
  # number of components
  if (!autosel) {
    if (mode(comps)!="numeric" || length(comps)!=1 || comps<=1 || (comps%%1)!=0)
      stop("\nInvalid argument 'comps' (number of components)")
  }
  # type of validation
  if (is.null(validation)) {
    validation = "none"
  } else {
    vali = validation %in% c("crossval", "learntest")
    if (!vali)
      stop("\nIncorrect type of validation")
  }
  
  # how many observations and variables
  n = nrow(X)
  p = ncol(X)
  # how many groups
  ng = nlevels(y)
  # how many obs in each group
  nobs_group = as.vector(table(y))
  # group levels
  glevs = levels(y)
  
  ## plsDA with no validation
  if (validation %in% c("none","crossval")) {
    get_plsda = my_plsDA(X, y, 1:n, 1:n, autosel, comps)
    err = 1 - sum(diag(get_plsda$conf)) / n
  }
  
  ## plsDA with learn-test sets validation
  if (validation == "learntest")
  {
    if (any(learn) <= 0 || any(learn) > n)
      stop("\nsubscript out of bounds in 'learn' set")
    if (any(test) <= 0 || any(test) > n)
      stop("\nsubscript out of bounds in 'test' set")
    # apply plsDA
    get_plsda = my_plsDA(X, y, learn, test, autosel, comps, scaled)
    # misclassification error rate
    err = 1 - sum(diag(get_plsda$conf))/length(test)
  }
  
  ## specifications
  specs = list(n=n, p=p, ng=ng, glevs=glevs, 
               nobs_group=nobs_group, validation=validation)
  ## results
  structure(list(functions = get_plsda$coeffs, 
                 confusion = get_plsda$conf,
                 scores = get_plsda$Disc, 
                 classification = get_plsda$pred_class,
                 error_rate = err,
                 components = get_plsda$components,
                 Q2 = get_plsda$Q2T,
                 R2 = get_plsda$R2,
                 VIP = get_plsda$VIP,
                 comp_vars = get_plsda$cor_tx,
                 comp_group = get_plsda$cor_ty,
                 specs = specs),
            class = "plsda")
}
