#' Calculate the Neyman-Pearson ROC
#'
#' \code{nproc} calculate the Neyman-Pearson ROC
#' curve for a given sequence of type I error values.
#' @export
#' @param x n * p observation matrix. n observations, p covariates.
#' @param y n 0/1 observatons.
#' @param methods classification method(s).
#' \itemize{
#' \item logistic: \link{glm} function with family = 'binomial'
#' \item penlog: \code{\link[glmnet]{glmnet}} in \code{glmnet} package
#' \item svm: \code{\link[e1071]{svm}} in \code{e1071} package
#' \item randomforest: \code{\link[randomForest]{randomForest}} in \code{randomForest} package
#' \item lda: \code{\link[MASS]{lda}} in \code{MASS} package
#' \item nb: \code{\link[e1071]{naiveBayes}} in \code{e1071} package
#' \item ada: \code{\link[ada]{ada}} in \code{ada} package
#' \item custom: a custom classifier. score vector needed.
#' }
#' @param kernel kernel used in the svm method. Default = 'radial'.
#' @param score score vector corresponding to y. Required when method  = 'custom'.
#' @param conf whether to generate two np roc curves representing a confidence band with prob 1-delta.
#' @param alphalist the sequence of type I error values. Default = seq(from=0,to=1,by=0.01).
#' @param delta the violation rate of the type I error. Default = 0.05.
#' @param split the number of splits for the class 0 sample. Default = 1. For ensemble
#' version, choose split > 1.  When method = 'custom',  split = 0 always.
#' @param cv whether cross-validation is performed for calculating the roc curve.
#' @param fold number of folds for the cross-validation. Default = 5.
#' @param loc.prob.lo the precalculated lower threshold locations in probability. Default = NULL.
#' @param loc.prob.up the precalculated upper threshold locations in probability. Default = NULL.
#' @param n.cores number of cores used for parallel computing. Default = 1.
#' @param randSeed the random seed used in the algorithm.
#' @seealso \code{\link{npc}}
#' @examples
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' #fit = nproc(x, y, method = 'svm')
#' #fit2 = nproc(x, y, method = 'svm', cv = TRUE)
#' fit3 = nproc(x, y, method = 'penlog')
#'
#' ##Plot the nproc curve
#' plot(fit3)
#' #fit3 = nproc(x, y, method = 'penlog',  n.cores = 2)
#' #In practice, replace 2 by the number of cores available 'detectCores()'
#' #fit4 = nproc(x, y, method = 'penlog', n.cores = detectCores())
#'
#' #Testing the custom method for nproc.
#' #fit = npc(x, y, method = 'lda', split = 0,  n.cores = 2) #use npc to get score list.
#' #obj = nproc(x = NULL, y = fit$y, method = 'custom', split = 0,
#' #score = fit$score,  n.cores = 2)
#'
#' #Compared different classifiers via nproc
#' #fit5 = nproc(x, y, method = c('lda','ada','randomforest'),  n.cores = detectCores())
#'
#' #Confidence nproc curves
#' #fit6 = nproc(x, y, method = 'lda', conf = TRUE)
#'
#' #nproc ensembled version
#' #fit7 = nproc(x, y, method = 'lda', split = 21)




nproc <- function(x = NULL, y, methods = c("logistic", "penlog", "svm", "randomforest",
                 "lda", "nb", "ada", "custom"), kernel = "radial", score = NULL, conf = FALSE,
                 alphalist = seq(from = 0, to = 1,by = 0.01), delta = 0.05,
                 split = 1, cv = FALSE, fold = 5, loc.prob.lo = NULL, loc.prob.up = NULL, n.cores = 1, randSeed = 0) {
  roc.lo = NULL
  roc.up = NULL
  methods = match.arg(methods, several.ok = TRUE)
  set.seed(randSeed)
  for(method in methods){
  alphalist = alphalist[alphalist >= 0 & alphalist <= 1]
  if (!conf) {
  v1 = nproc.core(x = x, y = y, method = method, score = score,
                 alphalist = alphalist, delta = delta,
                 split = split, cv = cv, fold = fold, loc.prob = loc.prob.lo, n.cores = n.cores)
  roc.lo = cbind(roc.lo,v1$v)
  loc.prob.lo = v1$loc.prob
  } else{
    v1 = nproc.core(x = x, y = y, method = method, score = score,
                   alphalist = alphalist, delta = delta/2,
                   split = split, cv = cv, fold = fold, loc.prob = loc.prob.lo, n.cores = n.cores)
    v2 = nproc.core(x = x, y = y, method = method, score = score,
                   alphalist = alphalist, delta = 1-delta/2,
                   split = split, cv = cv, fold = fold, loc.prob = loc.prob.up, n.cores = n.cores)
    roc.lo = cbind(roc.lo,v1$v)
    roc.up = cbind(roc.up,v2$v)
    loc.prob.lo = v1$loc.prob
    loc.prob.up = v2$loc.prob
  }
  }

  object = list(roc.lo = roc.lo, roc.up = roc.up, conf = conf, methods = methods, loc.prob.lo = loc.prob.lo, loc.prob.up = loc.prob.up, alphalist = alphalist, delta = delta)

  class(object) = "nproc"
  return(object)
}
