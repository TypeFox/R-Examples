##' 'kmlCov' re-launch the algorithm implemented in \link{glmClust}, for
##' clustering longitudinal data (trajectories), several
##' times with different starting conditions and various number of clusters.
##'
##' 
##' The purpose of \code{kmlCov} is clustering longitudinal data, as well as
##' \link{glmClust}, and automate the procedure of re-launching the algorithm
##' from different starting conditions by specifying \code{nRedraw}.\cr
##'
##' 
##' The algorithm depends greatly of the starting conditions (initial affectation on the
##'  trajectories/individuals), so it is recommanded
##' to run the algorithm multiple times in order to explore
##' the space of the solutions. \cr
##'
##'
##' 'kmlCov' return a list of list of \code{GlmCLuster}, the partitions
##' are compared using as criterion the \bold{classification log-likelihood},
##' the higher are the best partitions. 
##'
##'
###################################################################################' 
##' @title Clustering longitudinal data from different starting conditions
##' @usage kmlCov(formula, data, ident, timeVar, nClust = 2:6, nRedraw = 20, 
##' family = 'gaussian', effectVar = '', weights = rep(1,nrow(data)) ,
##' timeParametric = TRUE, separateSampling = TRUE, max_itr = 100, verbose = TRUE)
##' @param formula A symbolic description of the model. In the parametric case we write for
##' example 'y ~ clust(time+time2) + pop(sex)', here 'time' and 'time2' will have a different
##' effect according to the cluster, the 'sex' effect is the same for all the clusters. In the
##' non-parametric case only one covariate is allowed.
##' @param data  A [data.frame] in long format (no missing values) which means that each line
##' corresponds to one measure of the observed phenomenon, and one individual may have multiple
##' measures (lines) identified by an identity column. In the non-parametric case the totality
##' of patients must have all the measurements at all fixed times.
##' @param nClust The number of clusters, at leas 2 an at most 26.
##' @param nRedraw The number of time the algorithm is re-run with different starting conditions.
##' @param ident The name of the column identity.
##' @param timeVar Specify the column name of the time variable.
##' @param family  A description of the error distribution and link function to be used in the model,
##' by default 'gaussian'. This can be a character string naming a family function, a family
##' function or the result of a call to a family function.  (See 'family' for details of
##' family functions).
##' @param effectVar An effect, can be a level cluster effect or not.
##' @param weights Vector of 'prior weights' to be used in the fitting process, by default the weights are equal to one.
##' @param timeParametric By default [TRUE] thus parametric on the time. If [FALSE] then only one covariate is allowed in the formula and the algorithm used is the k-means.
##' @param separateSampling By default [TRUE] it means that the proportions of the clusters are supposed equal in the classification step, the log-likelihood maximised at each step of the algorithm is \eqn{\sum_{k=1}^{K}\sum_{y_i \in P_k} \log(f(y_i, \theta_k))}, otherwise the proportions of clusters are taken into account and the log-likelihood is \eqn{\sum_{k=1}^{K}\sum_{y_i \in P_k} \log(\lambda_{k}f(y_i, \theta_k))}.
##' @param max_itr The maximum number of iterations fixed at 100.
##' @param verbose Print the output in the console.
##' @return A an object of class \code{KmlCovList}.
##' @export
##' @seealso \link{glmClust} \cr \link{which_best}
##' @examples
##' data(artifdata)
##' res <- kmlCov(formula = Y ~ clust(time + time2), data = artifdata, ident = 'id',
##' timeVar = 'time', effectVar = 'treatment', nClust = 2:3, nRedraw = 2) #run 2 times for each cluster
##' 

kmlCov <- function(formula, data, ident, timeVar, nClust = 2:6, nRedraw = 20, family = 'gaussian',
                   effectVar = '', weights = rep(1,nrow(data)) , timeParametric = TRUE,
                   separateSampling = TRUE, max_itr = 100, verbose = TRUE) {
  part_clust <- list()                 # partitions par cluster 1,2...
  part_all <- list()                   # list de part_clust
  itr <- 1
  for (itrClust in nClust) {
    cat(itrClust,'Clusters : Running ')
    for (itrDraw in 1:nRedraw) {
      part_clust[[itrDraw]] <- glmClust(formula = formula, data = data, nClust = itrClust,
                                        ident = ident, timeVar = timeVar, family = family,
                                        effectVar = effectVar, weights = weights,
                                        timeParametric = timeParametric, separateSampling = separateSampling,
                                        max_itr = max_itr, verbose = FALSE)
      cat('.')
    }
    cat('End\n')
    part_all[[itr]] <- part_clust # stock la list des partitions ('GlmCluster')
    itr <- itr + 1
  }
  if (length(nClust) == 1) {
    return(new(Class = 'KmlCovList', list_part = part_clust))
  } else {
    return(new(Class = 'KmlCovList', list_part = part_all))
  }
}


##' Seek the best partitions in an object of class \code{KmlCovList} and
##' return the best one of each fixed number of cluster.
##' 
##' @title Seek the best partitions
##' @usage which_best(kmlcovar, crit = "log-class-likelihood")
##' @param kmlcovar An object of class \code{KmlCovList}.
##' @param crit Name of the criterion which have to be optimised, CLL for classification
##' log-likelihood AIC for Akaike information criterion and BIC for Bayesian information criterion.
##' @return An object of class \code{GlmCluster} or \code{KmlCovList}. 
##' @export
##' @seealso \link{kmlCov} 
##' @examples
##' data(artifdata)
##' res <- kmlCov(formula = Y ~ clust(time + time2), data = artifdata, ident = 'id',
##' timeVar = 'time', effectVar = 'treatment', nClust = 2:3, nRedraw = 2) # run 2 times the algorithm
##' best <- which_best(res) # return the best partition of each cluster
##' plot(best)

which_best <- function(kmlcovar, crit = 'log-class-likelihood') {
  if (class(kmlcovar)[1] == 'KmlCovList' & crit %in% c('log-class-likelihood', 'AIC', 'BIC')) {
    if (crit == 'log-class-likelihood') {
      min_max <- which.max
    } else {
      min_max <- which.min
    }
      
    if (class(kmlcovar@list_part[[1]])[1] == 'GlmCluster') {
      vec <- rep(NA, length(kmlcovar@list_part))

      for (itr in 1:length(kmlcovar@list_part)) {
        vec[itr] <- kmlcovar@list_part[[itr]]@criteria[, crit]
      }
      ind <- min_max(vec)
      return(kmlcovar@list_part[[ ind ]]) # renvoie GlmCluster
      
    } else {
      best_part <- list()

      for (itr_clust in 1:length(kmlcovar@list_part)) {
        vec <- rep(NA, length(kmlcovar@list_part[[itr_clust]]))
        for (itr in 1:length(kmlcovar@list_part[[ itr_clust ]])) {
          vec[itr] <- kmlcovar@list_part[[ itr_clust ]][[ itr ]]@criteria[, crit]
        }
        ind <- min_max(vec)
        best_part[[ itr_clust ]] <- kmlcovar@list_part[[ itr_clust ]][[ ind ]]
      }
      return(new(Class ='KmlCovList', list_part = best_part))
    }

  } else {
    cat('Wrong object and/or wrong name criterion')
    return(NULL)
  }
}
  

