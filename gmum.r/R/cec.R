#' @useDynLib gmum.r

#' @name CEC
#' @title Cross-Entropy Clustering
#' 
#' @description Create CEC model object
#'
#' @param x Numeric matrix of data.
#' @param k Initial number of clusters.
#' @param method.type Type of clustering (Gaussian family).
#' \enumerate{
#' \item 'diagonal' Gaussians with diagonal covariance. The clustering will try to divide the data into ellipsoid with radiuses parallel to coordinate axes
#' \item 'fixed_spherical' Spherical (radial) Gaussian densities (additional parameter - radius)
#' \item 'fixed_covariance' The clustering will have the tendency to divide the data into clusters resembling the unit circles in the Mahalanobis distance (additional parameter - covaraince matrix required)
#' \item 'func' Own function dependent on m and sigma (additional parameter)
#' \item 'standard' We divide dataset into ellipsoid-like clusters without any preferences (default)
#' \item 'spherical' The clustering will try to divide the data into circles of arbitrary sizes}
#' @param method.init Method to initialize clusters.
#' \enumerate{
#' \item 'centroids'
#' \item 'kmeans++'
#' \item 'random'}
#' @param params.r Radius for spherical family.
#' @param params.cov Covariance matrix for covariance family.
#' @param params.centroids List of centroids.
#' @param params.mix List of cluster with mixed Gaussian types.
#'  
#' @param params.function User energy function
#'
#' @param control.nstart How many times to perform algorithm.
#' @param control.eps What change of value should terminate algorithm.
#' @param control.itmax Maximum number of iterations at each start.
#' @param log.energy Records collected energy of all clusters in each iteration.
#' @param log.ncluster Records number of clusters in each iteration.
#' @param seed User seed
#' 
#' @examples
#' \dontrun{
#' CEC(k=3, x=dataset)
#' 
#' CEC(k=3, x=dataset, control.nstart=10, method.type='spherical', control.eps=0.05)
#' 
#' CEC(k=2, x=dataset, method.type='spherical', method.init='centroids', 
#'    params.centroids=list(c(-0.5,0.5),c(0,0)))
#' 
#' CEC(k=5, x=dataset, method.type='fixed_spherical', params.r=0.01, 
#'    control.nstart=10, control.eps=0.07)
#' 
#' CEC(k=5, x=dataset, method.type='fixed_covariance', 
#'    params.cov=matrix(c(0.03,0,0,0.01),2), control.nstart=10, control.eps=0.06)
#' 
#' CEC(k=1, x=dataset, method.type='func', 
#'    params.function='name_of_my_own_function')
#'    
#' fixed_spherical_cluster_param = list(method.type = 'fixed_spherical', params.r = 0.001),
#' covariance_cluster_param = list(method.type = 'fixed_covariance', 
#'    params.cov=matrix(c(0.05, 0, 0, 0.001), 2))
#' CEC(x = dataset, params.mix = list(covariance_cluster_param, 
#'    fixed_spherical_cluster_param, fixed_spherical_cluster_param, 
#'    fixed_spherical_cluster_param, fixed_spherical_cluster_param), control.nstart = 10)
#' 
#' p1 = list(method.type='spherical', k=3)
#' p2 = list(method.type='diagonal', k=2)
#' CEC(x=dataset, params.mix=list(p1, p2))
#' }
#' @import Rcpp
#' @export
CEC <- NULL


#' @name runAll
#' @title runAll
#' 
#' @aliases runAll,Rcpp_CecModel-method
#'
#' @description Starts whole algorithm again with same parameters
#'
#' @param c object Trained CEC model object.
#'
#' @examples
#' \dontrun{
#' runAll(c) 
#' }
#' @export
runAll <- NULL

#' @name runOneIteration
#' @title runOneIteration
#'
#' @aliases runOneIteration,Rcpp_CecModel-method
#' 
#' @description runs one iteration of algorithm
#'
#' @param c object Trained CEC model object.
#'
#' @examples
#' \dontrun{
#' runOneIteration(c) 
#' }
#' @export
runOneIteration <- NULL

#' @name energy
#' @title energy
#'
#' @aliases energy,Rcpp_CecModel-method
#'
#' @description Print result energy of clustering found
#'
#' @param c object Trained CEC model object.
#'
#' @examples
#' \dontrun{
#' energy(c) 
#' }
#' @export
energy <- NULL

#' @name clustering
#' @title clustering
#' @rdname clustering-methods
#' 
#' @description Print labels assigned
#'
#' @param c Object with clusters
#' @examples
#' \dontrun{
#' clustering(c) 
#' }
#' @export 
clustering <- function(c) UseMethod("clustering", c)

#' @name clustering.Rcpp_CecModel
#' @method clustering Rcpp_CecModel
#' @title clustering
#' @rdname clustering-methods
#' 
#' @description Print labels assigned
#' @examples
#' \dontrun{
#' clustering(c) 
#' }
#' @export 
clustering.Rcpp_CecModel <- NULL

#' @name getDataset 
#' @title getDataset 
#' 
#' @aliases getDataset,Rcpp_CecModel-method
#' 
#' @description Print input dataset 
#'
#' @param c object Trained CEC model object.
#'
#' @examples
#' \dontrun{
#' getDataset(c) 
#' }
#' @export
getDataset <- NULL

#' @name centers
#' @title centers
#' 
#' @aliases centers,Rcpp_CecModel-method
#'
#' @description Print centers of clusters
#'
#' @param c object Trained CEC model object.
#'
#' @examples
#' \dontrun{
#' centers(c) 
#' }
#' @export
centers <- NULL

#' @name covMatrix
#' @title covMatrix
#' 
#' @aliases covMatrix,Rcpp_CecModel-method
#' 
#' @description Print covariances of clusters
#'
#' @param c object Trained CEC model object.
#'
#' @examples
#' \dontrun{
#' covMatrix(c) 
#' }
#' @export
covMatrix <- NULL

#' @name predict
#' @rdname predict.cec
#' @method predict Rcpp_CecModel
#' @title predict
#' 
#' @aliases predict,Rcpp_CecModel-method
#'
#' @description Classify a new point according to the model (returns index of cluster where given point belong to)
#' 
#' @param object Trained CEC model object.
#' @param x Given point.
#' @param ... other arguments not used by this method.
#' @export 
predict.Rcpp_CecModel <- NULL

#' @name logClusters
#' @title logClusters
#' 
#' @aliases logClusters,Rcpp_CecModel-method
#'
#' @description Print number of clusters that has been recorded at each stage of learning.
#' Data is recorded only if you have chosen to when you created CEC model object.
#'
#' @param c object Trained CEC model object.
#' 
#' @usage logClusters(c)
#' 
#' @examples
#' \dontrun{
#' logClusters(c) 
#' }
#' @export
logClusters <- NULL

#' @name logEnergy
#' @title logEnergy
#' 
#' @aliases logEnergy,Rcpp_CecModel-method
#'
#' @description Print energy that has been recorded at each stage of learning.
#' Data is recorded only if you have chosen to when you created CEC model object.
#'  
#' @param c object Trained CEC model object.
#'
#' @examples
#' \dontrun{
#' logEnergy(c) 
#' }
#' @export
logEnergy <- NULL

#' @name iterations
#' @title iterations
#' 
#' @aliases iterations,Rcpp_CecModel-method
#'
#' @title iterations
#' 
#' @description Print how many iterations it took to learn CEC model
#'
#' @param c object Trained CEC model object.
#' @examples
#' \dontrun{
#' iterations(c) 
#' }
#' @export
iterations <- NULL

loadModule('cec', TRUE)

CEC <- function(x = NULL,
                 k = 0,
                 method.type = "standard",
                 method.init = "kmeans++",
                 params.r = 0,
                 params.cov = matrix(0),
                 params.centroids = NULL,
                 params.mix = NULL,
                 params.function = NULL,
                 control.nstart = 10,
                 control.eps = 0.05,
                 control.itmax = 25,
                 log.energy = TRUE,
                   log.ncluster= TRUE,
                   seed = NULL){
  
  # check for errors
  call <- match.call(expand.dots = TRUE)
  
  if (is.null(x))
    stop("Dataset is required!");
  
  if(is(x, "data.frame")){
    x = data.matrix(x);
  }
  
  if (is.null(params.mix) && k <= 0)
    stop("Number of clusters should be a positive integer!");
  
  if (control.nstart <= 0)
    stop("Number of starts should be a positive integer!");
  
  npoints <- dim(x)[1]
  if ( (control.eps < 0) || (control.eps > ((npoints - 1) / npoints)) )
    stop("control.eps = ", control.eps, " should be in range [0, (N-1)/N]!");  
  
  if (control.itmax < 0)
    stop("Maximum number of iterations should be a natural number!");
  
  if(is(params.cov, "data.frame")){
    params.cov = data.matrix(params.cov);
  }
  
  config <- new(CecConfiguration)

    if(is.null(seed) == FALSE) {
        config$.setSeed(seed)
    }
  config$.setDataSet(x)
  config$.setEps(control.eps)      
  config$.setNrOfClusters(k)
  
  if(is.null(params.mix) == FALSE) {
    config$.setMix(params.mix) 
  }
  
  if(is.null(params.function) == FALSE) {
    config$.setFunction(params.function)
  }

  config$.setLogEnergy(log.energy)
  config$.setLogCluster(log.ncluster)      
  config$.setNstart(control.nstart)
  config$.setCentroids(params.centroids)
  config$.setMethodType(method.type)             
  config$.setCov(params.cov)
  config$.setR(params.r)
  config$.setMethodInit(method.init) 
  config$.setItmax(control.itmax)
  config$.setAlgorithm('hartigan')
  
  model <- new(CecModel, config)
  
  
  
  assign("call", call, model)
  assign("energy", model$.getEnergy(), model)
  assign("clustering", model$.getClustering(), model)
  assign("centers", model$.getCenters(), model)
  assign("covMatrix", model$.getCovMatrix(), model)
  assign("logEnergy", model$.getLogEnergy(), model)
  assign("logNumberOfClusters", model$.getLogNumberOfClusters(), model)
  
  assign("iterations", model$.getIterations(), model)

  assign(".staticFields", c("call", "energy", "clustering", "centers", "covMatrix", 
                            "logNumberOfClusters", "logEnergy", "iterations"), model)
  
  model
}

runAll <- function(c) {
  c$.runAll()
}

runOneIteration <- function(c) {
  c$.runOneIteration()
}

energy <- function(c) {
  c$.getEnergy()
}

clustering.Rcpp_CecModel <- function(c) {
  c$.getClustering()
}

getDataset <- function(c) {
  c$.getDataset()
}

centers <- function(c) {
  c$.getCenters()
}

covMatrix <- function(c) {
  c$.getCovMatrix()
}

logEnergy <- function(c) {
  c$.getLogEnergy()
}

logNumberOfClusters <- function(c) {
  c$.getLogNumberOfClusters()
}


iterations <- function(c) {
  c$.getIterations()
}

predict.Rcpp_CecModel <- function(object, x, ...) {
  if ( !is(x, "data.frame") && !is(x, "matrix") && !is(x,"numeric")  ) {
    stop("Wrong target class, please provide data.frame, matrix or numeric vector")
  }
  
  if(is(x, "vector")){
    x = matrix(x, nrow=1, byrow=TRUE)
  }
  else if (!is(x, "matrix")) {
    x = data.matrix(x)
  }
  
  if(dim(object$.getDataset())[2] != dim(x)[2]){
    stop("Incompatible dimension!")
  }
  
  apply(x, 1, function(row) {
    object$.predict(row)
  })
}

#' Class Rcpp_CecModel.
#'
#' Class \code{Rcpp_CecModel} defines a CEC model class. 
#'
#' @name Rcpp_CecModel-class
#' @exportClass Rcpp_CecModel
setClass(Class = "Rcpp_CecModel")

#' Class Rcpp_CecConfiguration.
#'
#' Class \code{Rcpp_CecConfiguration} defines a CEC model configuration class. 
#'
#' @name Rcpp_CecConfiguration-class
#' @exportClass Rcpp_CecConfiguration
setClass(Class = "Rcpp_CecConfiguration")
