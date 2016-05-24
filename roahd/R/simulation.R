#' Exponential covariance function over a grid
#'
#' This function computes the discretisation of an exponential
#' covariance function of the form:
#'  \deqn{C( s, t ) = \alpha e^{ - \beta | s - t | }}
#' over a 1D grid \eqn{[t_0, t_1, \ldots, t_{P-1}]}, thus obtaining the
#' \eqn{P \times P} matrix
#' of values:
#'  \deqn{ C_{i,j} = C( t_i, t_j ) = \alpha e^{ - \beta | t_i - t_j | } .}
#'
#' @param grid a vector of time points.
#' @param alpha the alpha parameter in the exponential covariance formula.
#' @param beta the beta parameter in the exponential covariance formula.
#'
#' @seealso \code{\link{generate_gauss_fdata}},
#' \code{\link{generate_gauss_mfdata}}
#'
#' @examples
#'
#' grid = seq( 0, 1, length.out = 5e2 )
#'
#' alpha = 0.2
#' beta = 0.3
#'
#' dev.new()
#' image( exp_cov_function( grid, alpha, beta ),
#'        main = 'Exponential covariance function',
#'        xlab = 'grid', ylab = 'grid')
#'
#'
#' @export
exp_cov_function = function( grid, alpha, beta )
{
  return(  outer( grid, grid,
                  function( s, t )( alpha * exp( - beta * abs( s - t ) ) ) ) )
}


#' Generation of gaussian univariate functional data
#'
#' \code{generate_gauss_fdata} generates a dataset of univariate functional data
#' with a desired mean and covariance function.
#'
#' In particular, the following model is considered for the generation of data:
#'
#'  \deqn{X(t) = m( t ) + \epsilon( t ), \quad t \in I = [a, b]}{ X(t) =
#'  m( t ) + \epsilon(t), for all  t in I = [a, b] }
#'
#' where \eqn{m(t)} is the center and \eqn{\epsilon(t)} is a centered gaussian
#' process with covariance function \eqn{C_i}.
#' That is to say:
#'
#'  \deqn{Cov( \epsilon(s), \epsilon(t) ) = C( s, t ), \quad \forall s, t \in
#'   I}{Cov( \epsilon(s), \epsilon(t) ) = C( s, t ), with s, t in I}
#'
#' All the functions are supposed to be observed on an evenly-spaced, one-
#' dimensional grid of P points: \eqn{[a = t_0, t_1, \ldots, t_{P-1} = b]
#' \subset I }.
#'
#'
#' @param N the number of distinct functional observations to generate.
#' @param centerline the centerline of the distribution, represented as a one-
#' dimensional data structure  of length \eqn{P} containing the measurement of
#' the centerline on grid points.
#' @param Cov the covariance operator (provided in form of a \eqn{P \times P}{
#' P x P} matrix) that has to be used in the generation of \eqn{\epsilon(t)}. At
#' least one argument between \code{Cov} and \code{CholCov} should be different
#' from \code{NULL}.
#' @param CholCov the Cholesky factor of the covariance operator (provided in
#' form of a \eqn{P \times P}{P x P} matrix) that has to be used in the
#' generation of observations from the process \eqn{\epsilon(t)}. At least one
#' argument between \code{Cov} and \code{CholCov} should be different from
#' \code{NULL}.
#'
#'
#' @return The function returns a matrix containing the discretized
#' values of the generated observations (in form of an \eqn{N \times P}{N x P}
#' matrix).
#'
#' @seealso \code{\link{exp_cov_function}}, \code{\link{fData}},
#' \code{\link{generate_gauss_mfdata}}
#'
#'
#' @examples
#'
#' N = 30
#' P = 1e2
#'
#' t0 = 0
#' tP = 1
#'
#' time_grid = seq( t0, tP, length.out = P )
#'
#' C = exp_cov_function( time_grid, alpha = 0.1, beta = 0.2 )
#'
#' CholC = chol( C )
#'
#' centerline = sin( 2 * pi * time_grid )
#'
#' generate_gauss_fdata( N, centerline, Cov = C )
#'
#' generate_gauss_fdata( N, centerline, CholCov = CholC )
#'
#' @importFrom stats rnorm
#'
#' @export
generate_gauss_fdata = function( N, centerline,
                                 Cov = NULL, CholCov = NULL )
{

  centerline = as.numeric( centerline )

  if( is.null( Cov ) & is.null( CholCov ) ){
    stop( 'Error: You have to provide at least either covariance matrix or
          its cholesky factor to .generate_gauss_fdata\n')
  } else if( ! is.null( CholCov ) ) {

    P = ncol( CholCov )

    if( length( centerline ) != nrow( CholCov ) | nrow( CholCov ) != P  ){
      stop( 'Error: You provided mismatching centerline and covaraince matrix
            Cholesky factor to generate_gauss_fdata\n')
    }
    } else if( ! is.null( Cov ) ){

      P = ncol( Cov )

      if( length( centerline ) != nrow( Cov ) | nrow( Cov ) != P  )
      {
        stop( 'Error: You provided mismatching centerline and covariance matrix
to generate_gauss_fdata\n')
      }

      CholCov = chol( Cov )
      }

  return( t( t( matrix( rnorm( N * P ),
                  nrow = N,
                  ncol = P ) %*% CholCov ) + centerline ) )
}

#' Generation of gaussian multivariate functional data
#'
#' \code{generate_gauss_mfdata} generates a dataset of multivariate functional
#' data with a desired mean and covariance function in each dimension and a
#' desired correlation structure among components.
#'
#' In particular, the following model is considered for the generation of data:
#'
#'  \deqn{X(t) = ( m_1( t ) + \epsilon_1( t ), \ldots, m_L(t) +
#'  \epsilon_L(t)), \quad t \in I = [a, b]}{ X(t) = ( m_1( t ) +
#'  \epsilon_1(t), ..., m_L(t) + \epsilon_L(t) ), for all  t in I = [a, b] }
#'
#' where \eqn{L} is the number of components of the multivariate functional
#' random variable, \eqn{m_i(t)} is the \eqn{i-}th component of the center and
#' \eqn{\epsilon_i(t)} is a centered gaussian process with covariance function
#' \eqn{C_i}. That is to say:
#'
#'  \deqn{Cov( \epsilon_{i}(s), \epsilon_{i}(t) ) = C( s, t ), \quad \forall i =
#'  1, \ldots, L, \quad \forall s, t \in I}{Cov( \epsilon_i(s),
#'   \epsilon_i(t) ) = C( s, t ),  with  i =
#'  1, \ldots, L, and with s, t in I}
#'
#' A correlation structure among \eqn{\epsilon_1(t),\ldots,\epsilon_L(t)} is
#' allowed in the following way:
#'
#' \deqn{ Cor( \epsilon_i(t), \epsilon_j(t)) = \rho_{i,j}, \quad \forall
#' i \neq j, \quad \forall t \in I.}{ Cor( \epsilon_i(t), \epsilon_j(t)  ) =
#' \rho_ij, for all i != j and for all t in I.}
#'
#' All the functions are supposed to be observed on an evenly-spaced, one-
#' dimensional grid of P points: \eqn{[ a = t_0, t_1, \ldots, t_{P-1} = b]
#' \subset I }.
#'
#'
#' @param N the number of distinct functional observations to generate.
#' @param L the number of components of the multivariate functional data.
#' @param centerline the centerline of the distribution, represented as a
#' 2-dimensional data structure with L rows (one for each dimension) having the
#' measurements along the grid as columns.
#' @param correlations is the vector containing the \eqn{1/2 L (L-1)}
#' correlation coefficients \eqn{\rho_{ij}} in the model generating data.
#' They have to be provided in the following order:
#' \deqn{(\rho_{1,2},\ldots,\rho_{1,L},\rho_{2,3},\ldots,\rho_{2,L},\ldots,
#' \rho_{L,L-1}),}
#' that is to say, the row-wise, upper triangular part of the correlation matrix
#' without the diagonal.
#' @param listCov a list containing the \eqn{L} covariance operators (provided
#' in form of a \eqn{P \times P}{P x P} matrix), one for each component of the
#' multivariate functional random vairable, that have to be used in the
#' generation of the processes \eqn{\epsilon_1(t), \ldots, \epsilon_L(t)}.
#' At least one argument between \code{listCov} and \code{listCholCov} must be
#' different from \code{NULL}.
#' @param listCholCov the Cholesky factor of the \eqn{L} covariance operators
#' (in \eqn{P \times P}{P x P} matrix form), one for each component of the
#' multivariate functional random vairable, that have to be used in the
#' generation of the processes \eqn{\epsilon_1(t), \ldots, \epsilon_L(t)}.
#' At least one argument between \code{listCov} and \code{listCholCov} must be
#' different from \code{NULL}.
#'
#' @return The function returns a list of L matrices, one for each component of
#' the multivariate functional random variable, containing the discretized
#' values of the generated observations (in form of \eqn{N \times P}{N x P}
#' matrices).
#'
#' @seealso \code{\link{exp_cov_function}}, \code{\link{mfData}},
#' \code{\link{generate_gauss_fdata}}
#'
#' @examples
#'
#' N = 30
#' P = 1e2
#' L = 3
#'
#' time_grid = seq( 0, 1, length.out = P )
#'
#' C1 = exp_cov_function( time_grid, alpha = 0.1, beta = 0.2 )
#' C2 = exp_cov_function( time_grid, alpha = 0.2, beta = 0.5 )
#' C3 = exp_cov_function( time_grid, alpha = 0.3, beta = 1 )
#'
#'
#' centerline = matrix( c( sin( 2 * pi * time_grid ),
#'                         sqrt( time_grid ),
#'                         10 * ( time_grid - 0.5 ) * time_grid ),
#'                      nrow = 3, byrow = TRUE )
#'
#' generate_gauss_mfdata( N, L, centerline,
#'                        correlations = c( 0.5, 0.5, 0.5 ),
#'                        listCov = list( C1, C2, C3 ) )
#'
#' CholC1 = chol( C1 )
#' CholC2 = chol( C2 )
#' CholC3 = chol( C3 )
#'
#' generate_gauss_mfdata( N, L, centerline,
#'                        correlations = c( 0.5, 0.5, 0.5 ),
#'                       listCholCov = list( CholC1, CholC2, CholC3 ) )
#'
#' @export
generate_gauss_mfdata = function( N, L, centerline, correlations,
                                  listCov = NULL, listCholCov = NULL )
{
  if( length( correlations ) != 0.5 * ( L ) * ( L - 1 ) )
  {
    stop( 'Error in generate_gauss_mfdata: you have to provide all the
          correlations among functional components' )
  }

  if( nrow( centerline ) != L )
  {
    stop( 'Error in generate_gauss_mfdata: you have to provide a centerline for
each dimension' )
  }

  if( is.null( listCov ) & is.null( listCholCov ) ){
    stop( 'Error: You have to provide at least either covariance matrices or
          their cholesky factors to generate_gauss_mfdata')
  } else if( ! is.null( listCholCov ) )
  {
      if( length( listCholCov ) != L )
      {
        stop( 'Error: You have to provide a covariance Cholesky factor for each
              dimension' )
      }

      P = ncol( listCholCov[[ 1 ]] )

      if( ncol( centerline ) != P | any( sapply( listCholCov, nrow ) != P ) |
          any( sapply( listCholCov, ncol ) != P ) )
      {
        stop( 'Error: You provided mismatching centerline and covariance
matrices Cholesky factors to generate_gauss_mfdata')
      }
    } else if( ! is.null( listCov ) )
    {

      P = ncol( listCov[[ 1 ]] )

        if( ncol( centerline ) != P | any( sapply( listCov, nrow ) != P ) |
            any( sapply( listCov, ncol ) != P ) )
        {
          stop( 'Error: You provided mismatching centerline and covariance
matrices to generate_gauss_mfdata')
        }
      listCholCov = lapply( listCov, chol )
    }

  # Generating the matrix of correlations among dimensions
  R = matrix( 1, ncol = L, nrow = L )

  R[ upper.tri( R ) ] = as.numeric( correlations )

  R[ lower.tri( R ) ] = as.numeric( correlations )

  R_chol = chol( R )

  # Generating gaussian data with correlations among dimensions
  Data = matrix( rnorm( N * L * P ), ncol = L, nrow = N * P )

  Data = Data %*% R_chol

  return( values = eval(
    parse( text =
             paste( 'list( ',
                    paste( 't( t( matrix( Data[ , ', 1 : L,
                           ' ], nrow = N, ncol = P ) %*% listCholCov[[ ',
                           1 : L, ' ]] ) + as.numeric( centerline[ ',
                           1 : L, ', ] ) )',
                           sep = '', collapse  = ', ' ),
                    ' )', sep = '' ) ) ) )
}
