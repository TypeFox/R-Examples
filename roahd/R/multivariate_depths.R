
#' (Modified) Band Depth for multivariate functional data
#'
#' These functions compute the Band Depth (BD) and Modified Band Depth (MBD) of
#' elements of a multivariate functional dataset.
#'
#' Given a multivariate functional dataset composed of \eqn{N} elements with
#' \eqn{L} components each, \eqn{\mathbf{X_1} =( X^1_1(t),} \eqn{X^2_1(t),
#' \ldots, X^L_1(t))}, and a set of \eqn{L} non-negative weights,
#'
#' \deqn{ w_1, w_2, \ldots, w_L, \qquad \sum_{i=1}^L w_i = 1,}
#'
#' these functions compute the BD and MBD of each element of the functional
#' dataset, namely:
#'
#' \deqn{ BD( \mathbf{X_j} ) = \sum_{i=1}^{L} w_i BD( X^i_j ), \quad \forall
#' j = 1, \ldots N.}
#'
#' \deqn{ MBD( \mathbf{X_j} ) = \sum_{i=1}^{L} w_i MBD( X^i_j ), \quad \forall
#' j = 1, \ldots N.}
#'
#'
#' @param Data specifies the the multivariate functional dataset.
#' It is either an object of class \code{mfData} or a list of 2-dimensional
#' matrices having as rows the elements of that component and as columns the
#' measurements of the functional data over the grid.
#' @param weights either a set of weights (of the same length of \code{listData}
#' ) or the string \code{"uniform"} specifying that a set of uniform weights
#' (of value \eqn{1 / L}, where \eqn{L} is the number of dimensions of the
#' functional dataset and thus the length of \code{listData}) is to be used.
#' @param manage_ties a logical flag specifying whether the check for ties and
#' the relative treatment is to be carried out while computing the MBDs in each
#' dimension. It is directly passed to \code{MBD}.
#'
#' @return The function returns a vector containing the depths of each element
#' of the multivariate functional dataset.
#'
#' @aliases multiBD
#'
#' @references
#'
#' Ieva, F. and Paganoni, A. M. (2013). Depth measures for multivariate
#' functional data, \emph{Communications in Statistics: Theory and Methods},
#' 41, 1265-1276.
#'
#' Tarabelloni, N., Ieva, F., Biasi, R. and Paganoni, A. M. (2015). Use of
#' Depth Measure for Multivariate Functional Data in Disease Prediction: An
#' Application to Electrocardiograph Signals, \emph{International Journal of
#' Biostatistics}, 11.2, 189-201.
#'
#' @seealso \code{\link{MBD}}, \code{\link{BD}}, \code{\link{toListOfValues}},
#' \code{\link{mfData}}
#'
#' @examples
#'
#' N = 20
#' P = 1e3
#'
#' grid = seq( 0, 10, length.out = P )
#'
#' # Generating an exponential covariance function to be used to simulate gaussian
#' # functional data
#' Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.8 )
#'
#' # First component of the multivariate guassian functional dataset
#' Data_1 = generate_gauss_fdata( N, centerline = rep( 0, P ), Cov = Cov )
#'
#' # First component of the multivariate guassian functional dataset
#' Data_2 = generate_gauss_fdata( N, centerline = rep( 0, P ), Cov = Cov )
#'
#' mfD = mfData( grid, list( Data_1, Data_2 ) )
#'
#' multiBD( mfD, weights = 'uniform' )
#' multiMBD( mfD, weights = 'uniform', manage_ties = TRUE )
#'
#' multiBD( mfD, weights = c( 1/3, 2/3 ))
#' multiMBD( mfD, weights = c( 1/3, 2/3 ), manage_ties = FALSE )
#'
#' multiBD( list( Data_1, Data_2 ), weights = 'uniform')
#' multiMBD( list( Data_1, Data_2 ), weights = 'uniform', manage_ties = TRUE )
#'
#' multiBD( list( Data_1, Data_2 ), weights = c( 1/3, 2/3 ))
#' multiMBD( list( Data_1, Data_2 ), weights = c( 1/3, 2/3 ), manage_ties = FALSE )
#'
#' @export
#'
multiMBD = function( Data, weights = 'uniform', manage_ties = FALSE )
{
  UseMethod( 'multiMBD', Data )
}

#' @rdname multiMBD
#'
#' @aliases multiMBD
#'
#' @export
multiMBD.mfData = function( Data, weights = 'uniform', manage_ties = FALSE )
{
  Data = toListOfValues( Data )
  NextMethod()
}

#' @rdname multiMBD
#'
#' @aliases multiMBD
#'
#' @export
multiMBD.default = function( Data, weights = 'uniform', manage_ties = FALSE )
{
  L = length( Data )

  N = nrow( Data[[ 1 ]] )
  P = ncol( Data[[ 2 ]] )

  if( ! all( sapply( Data, nrow ) - N == 0 ) |
      ! any( sapply( Data, ncol ) - P == 0 ) )
  {
    stop( ' Error in multiMBD: you provided a list with mismatching univariate
          functional datasets' )
  }

  if( is.character( weights ) )
  {
    if( weights == 'uniform' )
    {
      weights = rep( 1 / L, L )
    } else {
      stop( ' Error in multiMBD: misspecified weights characetr information')
    }
  } else {
    if( length( weights ) != L )
    {
      stop( ' Error in multiMBD: you provided more/fewer weights than
            required' )
    } else if( sum( weights ) != 1 ){
      stop( ' Error in multiMBD: sum of weights is not 1')
    }
    weights = as.numeric( weights )
  }

  return( as.numeric( sapply( Data, MBD, manage_ties ) %*% weights ) )
}

#' @rdname multiMBD
#'
#' @aliases multiMBD
#'
#' @export
multiBD = function( Data, weights = 'uniform' )
{
  UseMethod( 'multiBD', Data )
}

#' @rdname multiMBD
#'
#' @aliases multiMBD
#'
#' @export
multiBD.mfData =  function( Data, weights = 'uniform' )
{
  Data = toListOfValues( Data )
  NextMethod()
}

#' @rdname multiMBD
#'
#' @aliases multiMBD
#'
#' @export
multiBD.default = function( Data, weights = 'uniform' )
{
  L = length( Data )

  N = nrow( Data[[ 1 ]] )
  P = ncol( Data[[ 2 ]] )

  if( ! all( sapply( Data, nrow ) - N == 0 ) |
      ! any( sapply( Data, ncol ) - P == 0 ) )
  {
    stop( ' Error in multiBD: you provided a list with mismatching univariate
          functional datasets' )
  }

  if( is.character( weights ) )
  {
    if( weights == 'uniform' )
    {
      weights = rep( 1 / L, L )
    } else {
      stop( ' Error in multiBD: misspecified weights characetr information')
    }
  } else {
    if( length( weights ) != L )
    {
      stop( ' Error in multiBD: you provided more/fewer weights than
            required' )
    } else if( sum( weights ) != 1 ){
      stop( ' Error in multiMBD: sum of weights is not 1')
    }
    weights = as.numeric( weights )
    }

  return( as.numeric( sapply( Data, BD ) %*% weights ) )
  }
