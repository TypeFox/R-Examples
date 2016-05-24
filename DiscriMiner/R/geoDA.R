#' Geometric Predictive Discriminant Analysis
#' 
#' Performs a Geometric Predictive Discriminant Analysis
#' 
#' When \code{validation=NULL} there is no validation \cr When
#' \code{validation="crossval"} cross-validation is performed by randomly
#' separating the observations in ten groups. \cr When
#' \code{validation="learntest"} validationi is performed by providing a
#' learn-set and a test-set of observations. \cr
#' 
#' @param variables matrix or data frame with explanatory variables
#' @param group vector or factor with group memberships
#' @param validation type of validation, either \code{"crossval"} or
#' \code{"learntest"}. Default \code{NULL}
#' @param learn optional vector of indices for a learn-set. Only used when
#' \code{validation="learntest"}. Default \code{NULL}
#' @param test optional vector of indices for a test-set. Only used when
#' \code{validation="learntest"}. Default \code{NULL}
#' @return An object of class \code{"geoda"}, basically a list with the
#' following elements:
#' @return \item{functions}{table with discriminant functions}
#' @return \item{confusion}{confusion matrix}
#' @return \item{scores}{discriminant scores for each observation}
#' @return \item{classification}{assigned class}
#' @return \item{error_rate}{misclassification error rate}
#' @author Gaston Sanchez
#' @seealso \code{\link{classify}}, \code{\link{desDA}}, \code{\link{linDA}},
#' \code{\link{quaDA}}, \code{\link{plsDA}}
#' @references Lebart L., Piron M., Morineau A. (2006) \emph{Statistique
#' Exploratoire Multidimensionnelle}. Dunod, Paris.
#' 
#' Saporta G. (2006) \emph{Probabilites, analyse des donnees et statistique}.
#' Editions Technip, Paris.
#' 
#' Tuffery S. (2011) \emph{Data Mining and Statistics for Decision Making}.
#' Wiley, Chichester.
#' @export
#' @examples
#' 
#'   \dontrun{
#'   # load bordeaux wines dataset
#'   data(iris)
#' 
#'   # geometric predictive discriminant analysis with no validation
#'   my_geo1 = geoDA(iris[,1:4], iris$Species)
#'   my_geo1$confusion
#'   my_geo1$error_rate
#' 
#'   # geometric predictive discriminant analysis with cross-validation
#'   my_geo2 = geoDA(iris[,1:4], iris$Species, validation="crossval")
#'   my_geo2$confusion
#'   my_geo2$error_rate
#'   }
#' 
geoDA <- 
function(variables, group, validation = NULL, learn = NULL, test = NULL)
{
  # Perform a geometric predictive discriminant analysis
  # variables: matrix or data.frame with explanatory variables
  # group: vector or factor with group membership
  # validation: NULL, "crossval", "learntest"
  # learn: vector of learn-set
  # test: vector of test-set
  
  # check inputs
  verify_Xy = my_verify(variables, group, na.rm=FALSE)
  X = verify_Xy$X
  y = verify_Xy$y
  # type of validation
  if (is.null(validation)) {
    validation = "none"
  } else {
    vali = validation %in% c("crossval", "learntest")
    if (!vali)
      stop("\nIncorrect type of validation")
  }
  
  # how many observations
  n = nrow(X)
  # how many variables
  p = ncol(X)
  # how many groups
  ng = nlevels(y)
  glevs = levels(y)
  # how many obs in each group
  nobs_group = as.vector(table(y))
  # proportions
  props = nobs_group / n
  
  ## geoDA with no validation
  if (validation == "none") {
    get_geoda = my_geoDA(X, y, 1:n, 1:n)
    err = 1 - sum(diag(get_geoda$conf)) / n
  }
  
  ## geoDA with learn-test sets validation
  if (validation == "learntest")
  {
    if (any(learn) <= 0 || any(learn) > n)
      stop("\nsubscript out of bounds in 'learn' set")
    if (any(test) <= 0 || any(test) > n)
      stop("\nsubscript out of bounds in 'test' set")
    # apply linDA
    get_geoda = my_geoDA(X, y, learn, test)
    # misclassification error rate
    err = 1 - sum(diag(get_geoda$conf))/length(test)
  }
  
  ## geoDA with crossvalidation
  if (validation == "crossval")
  {
    # geoDA for all observations
    get_geoda = my_geoDA(X, y, 1:n, 1:n)
    # elements in each group 
    elems_group = vector("list", ng)
    for (k in 1:ng) {
      elems_group[[k]] = which(group == glevs[k])
    }
    # misclassification error rate
    mer = 0
    # 10 crossvalidation samples
    for (r in 1:10)
    {
      test = vector("list", ng)
      test_sizes = floor(n * props / 10)
      for (k in 1:ng) {
        test[[k]] = sample(elems_group[[k]], test_sizes[k])
      }
      test = unlist(test)
      learn = (1:n)[-test]
      # apply DA
      geoda_cv = my_geoDA(X, y, learn, test)
      # misclassification error rate
      mer = mer + sum(diag(geoda_cv$conf))/n
    }
    # total misclassification error rate
    err = 1 - mer
  }
  
  ## specifications
  specs = list(n=n, p=p, ng=ng, glevs=glevs, 
               nobs_group=nobs_group, validation=validation)
  ## results
  structure(list(functions = get_geoda$FDF, 
                 confusion = get_geoda$conf,
                 scores = get_geoda$Disc, 
                 classification = get_geoda$pred_class,
                 error_rate = err,               
                 specs = specs),
            class = "geoda")
}
