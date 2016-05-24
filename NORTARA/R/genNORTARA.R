#'Generates samples with specified input correlation matrix and marginal
#'distributions.
#'
#'The function simulates a data set with specified input correlation matrix
#'\code{cor_matrix} and pre-specified marignals \code{invcdfnames} using
#'bounding RA and NORTA methods.
#'
#'The function simulates a date set with varibles from arbitrary(continuous or
#'discrete) marginal distributions which have a correlation matrix
#'\code{cor_matrix}. The pre-specified marignals are described by
#'\code{invcdfnames},\code{paramslists},\code{defaultindex}, the later two
#'arguments will be combined into a full paramslists which has the same length
#'as \code{invcdfnames}. The function uses result of the function
#'\code{BoundingRA} which is an implementation of a specific RA(Retrospective
#'Approximation) algorithm called bounding RA. With the result, the function
#'uses the NORTA(NORmal To Anything) approach which generates a standard normal
#'random vector and then transforms it into a random vector with specified
#'marginal to generates the wanted samples.
#'
#'@param n Number of observations.
#'@param cor_matrix specified input correlation matrix.
#'@param invcdfnames A character sequence of the marginals' inverse
#'  cdf(cumulative distribution function) names.
#'@param paramslists A list contains lists of params of the marginals excluded
#'  the index(es) in \code{defaultindex} meanwhile as the same order as
#'  invcdfnames, the names of the arguments of the inner lists should keep the
#'  same with the function arguments matching rules with the arguments of
#'  invcdfnames functions.
#'@param defaultindex The index number sequence which indicates the
#'  corresponding inverse cdfs use the default argument values.
#'@param m1 The initial sample size.
#'@param c1 The sample-size multiplier(c1>1).
#'@param c2 The step-size multiplier(c2>0).
#'@param delta1 The initial step size(detla1>0).
#'@param sigma0 The standard error tolerance.
#'@param epsilon The initial error tolerance.
#'@param maxit The maximum number of numerical searches.
#'
#'@return A matrix of size n * (\code{ncol(cor_matrix)}) from pre-specified
#'  marignals which also have an asymptotically correlation matrix to specified
#'  input correlation matrix \code{cor_matrix}.
#'@author Po Su
#'@seealso \code{\link{BoundingRA}}, \code{\link{valid_input_cormat}},
#'  \code{\link{check_input_cormat}}
#'@references Huifen Chen, (2001) Initialization for NORTA: Generation of Random
#'  Vectors with Specified Marginals and Correlations. INFORMS Journal on
#'  Computing \bold{13(4):312-331.}'
#'@note \enumerate{ \item  The inverse cdf functions should have its first
#'  argument as a vector. \item The inverse cdf functions like \code{qweibull}
#'  should not use default values because the \code{shape} argument does not has
#'  a default value.So we should not put its index in \code{invcdfnames} into
#'  \code{defaultindex}. \item The function won't accept the invcdfnames all be
#'  'qnorm", you can generate multi-normal varibles by other packages.The
#'  function will give you a error message, but if you use a = qnorm, then using
#'  new name, the function won't find it, just be careful.\item You may get a
#'  warning message indicate the Nearest positive definite matrix is used. It
#'  happens when your inputs won't generate a positive definite intermediate
#'  normal correlation matrix. And in this case, the cor(res) may not very close
#'  to cor_matrix. }
#'@examples
#'\dontrun{
#'invcdfnames <- c("qt","qpois","qnorm","qweibull","qunif")
#'# The following usage :
#'# a <- qt; b <- qnorm; f <- stats::qweibull (It is also the way you can use functions
#'# from other packages)
#'# invcdfnames <- c("a","qpois","b","f","qunif") will also be ok!
#'paramslists <- list(
#'                m1 = list(df = 3),
#'                m2 = list(lambda = 5),
#'                m4 = list(shape = 1)
#'                  )
#'defaultindex <- c(3,5)
#'#It means the 3rd and 5th invcdf should use its default arguments.
#'#That means qnorm using mean = 0, sd = 1, qunif using min =0 ,max =1 and so on.
#'cor_matrix <- matrix(c(1.0,-0.4,0.1,0.7,-0.2,-0.4,
#'                       1.0,0.4,0.4,0.9,0.1,0.4,1.0,
#'                       0.5,0.5,0.7, 0.4,0.5,1.0,
#'                       0.7,-0.2,0.9,0.5,0.7,1.0),5,5)
#'
#'res <- genNORTARA(10000,cor_matrix,invcdfnames,paramslists,defaultindex)
#'#May get warning message indicating nearest positive definite is used,It's
#'#normal but the cor(res) may not very close to cor_matrix.
#'cor(res)
#'invcdfnames <- c("qt","qpois","qnorm")
#'paramslists <- list(
#'                m1 = list(df = 3),
#'                m2 = list(lambda = 5)
#'                  )
#'defaultindex <- 3
#'cor_matrix <- matrix(c(1,0.5,-0.3,0.5,1,0.4,-0.3,0.4,1), 3)
#'res <- genNORTARA(10000,cor_matrix,invcdfnames,paramslists,defaultindex)
#'cor(res) #This time cor(res) may very close to cor_matrix.
#'}
#'@export
genNORTARA <- function(n, cor_matrix, invcdfnames,
                      paramslists = NULL, defaultindex = NULL,
                      m1 = 60, c1 = 2, c2 = 1, delta1 = 0.0001,
                       sigma0 = 0.01, epsilon = 1e+50, maxit = 1000) {
  if (length(invcdfnames) != length(paramslists)+length(defaultindex)) {
    stop("inversecdfs should have the same length paramslists,check paramslists and defaultindex !")
  }
  if (all(invcdfnames == "qnorm")) {
    stop("You should not use genNORTARA to generate multi-normal varibles!\n
          Please choose anthor package to do so if you still want!")
  }
#----------------------------------------------------------
  # Create inversecdfs's paramslists from raw input
  tmp_paramslists <- list()
  cdfscount <- length(invcdfnames)
  length(tmp_paramslists) <- cdfscount   # create list of NULLs
  if (length(defaultindex) != 0) {
    count <- 0
    for (i in 1:(cdfscount)) {
      if (!i %in% defaultindex) {
        count <- count + 1
        tmp_paramslists[[i]] <- paramslists[[count]]
      }
    }
    paramslists <- tmp_paramslists
    rm(tmp_paramslists)
  }
  #------------------------------------------------------------
  #Compute the intermediate normal correlation matrix
  if (check_input_cormat(invcdfnames, paramslists, cor_matrix)) {
   Sigma <- BoundingRA(cor_matrix, invcdfnames,
                         paramslists,  m1,
                         c1, c2, delta1,
                         sigma0, epsilon, maxit)
  #Generate virables with NoRTA by using Sigma
  ndim <- ncol(cor_matrix)
  transform_mat <- NULL
  n_correlated_standard_normal <- matrix(rnorm(n * ndim), nrow = n) %*% chol(Sigma)

  for (i in 1:ndim) {

    funcall <- as.call(c(as.name(invcdfnames[i]),
                         list(pnorm(n_correlated_standard_normal)[ ,i]), paramslists[[i]]))
    transform_mat <- cbind(transform_mat, eval(funcall))
  }
  return(transform_mat)
  }
}
