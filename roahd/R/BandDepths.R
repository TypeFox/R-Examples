#' Band Depth for univariate functional data
#'
#' This function computes the Band Depth (BD) of elements of a functional
#' dataset.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' this function computes the sample BD of each element with respect to the
#' other elements of the dataset, i.e.:
#'
#' \deqn{ BD( X( t ) ) = {N \choose 2 }^{-1} \sum_{1 \leq i_1 < i_2 \leq N} I(
#' G(X) \subset B( X_{i_1}, X_{i_2} ) ),}
#'
#' where \eqn{G(X)} is the graphic of \eqn{X(t)}, \eqn{B(X_{i_1},X_{i_2})} is
#' the envelope of \eqn{X_{i_1}(t)} and \eqn{X_{i_2}(t)}, and \eqn{X \in
#' \left\{X_1, \ldots, X_N\right\}}.
#'
#' See the References section for more details.
#'
#'
#' @param Data either an object of class \code{fData} or a matrix-like dataset
#' of functional data (e.g. \code{fData$values}),
#' with observations as rows and measurements over grid points as columns.
#'
#' @return The function returns a vector containing the values of BD for the
#' given dataset.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2009). On the Concept of Depth for Functional
#' Data, \emph{Journal of the American Statistical Association}, 104, 718-734.
#'
#' Lopez-Pintado, S. and Romo. J. (2007). Depth-based inference for functional
#' data, \emph{Computational Statistics & Data Analysis} 51, 4957-4968.
#'
#' @seealso \code{\link{MBD}}, \code{\link{BD_relative}},
#' \code{\link{MBD_relative}}, \code{\link{fData}}
#'
#' @examples
#'
#' grid = seq( 0, 1, length.out = 1e2 )
#'
#'
#' D = matrix( c( 1 + sin( 2 * pi * grid ),
#'                0 + sin( 4 * pi * grid ),
#'                1 - sin( pi * ( grid - 0.2 ) ),
#'                0.1 + cos( 2 * pi * grid ),
#'                0.5 + sin( 3 * pi + grid ),
#'                -2 + sin( pi * grid ) ),
#'             nrow = 6, ncol = length( grid ), byrow = TRUE )
#'
#' fD = fData( grid, D )
#'
#' BD( fD )
#'
#' BD( D )
#'
#' @export
#'
BD = function( Data )
{
  UseMethod( 'BD', Data )
}

#' @rdname BD
#'
#' @aliases BD
#'
#' @export
BD.fData = function( Data )
{
  Data = Data$values
  NextMethod()
}

#' @rdname BD
#'
#' @aliases BD
#'
#' @export
BD.default = function( Data )
{
  # Number of rows
  N = nrow( Data )

  # Compute ranks of matrix-like representation of data with `min' tie-breaking rule
  rk = apply( Data, 2, function( v )( rank( v, ties.method = 'average' ) ) )

  # Actual number of function values strictly above
  N_a = N - apply( rk, 1, max )

  # Actual number of function values strictly below
  N_b = apply( rk, 1, min ) - 1

  Depths = ( N_a * N_b + ( N - 1 ) ) / ( N * ( N - 1 ) / 2 )

  return( Depths )
}

#' Relative Band Depth of functions in a univariate functional dataset
#'
#' This function computes Band Depth (BD) of elements of a univariate functional
#' dataset with respect to another univariate functional dataset.
#'
#' Given a univariate functional dataset of elements \eqn{X_1(t), X_2(t),
#' \ldots, X_N(t)}, and another univariate functional dataset of elements
#' \eqn{Y_1(t), Y_2(t) \ldots, Y_M(t)}, this function computes the BD of
#' elements of the former with respect to elements of the latter, i.e.:
#'
#'  \deqn{ BD( X_i( t ) ) = {M \choose 2 }^{-1} \sum_{1 \leq i_1 < i_2 \leq M} I(
#' G(X_i) \subset B( Y_{i_1}, Y_{i_2} ) ), \quad \forall i = 1, \ldots, N,}
#'
#' where \eqn{G(X_i)} is the graphic of \eqn{X_i(t)} and \eqn{B(Y_{i_1},Y_{i_2})} is
#' the envelope of \eqn{Y_{i_1}(t)} and \eqn{Y_{i_2}(t)}.
#'
#' @param Data_target is the univariate functional dataset, provided either as
#' a \code{fData} object or in matrix
#' form (N observations as rows and P measurements as columns), whose BD
#' have to be computed with respect to the reference dataset.
#' @param Data_reference is the dataset, provided either as a \code{fData}
#' object or in matrix form (N observations
#' as rows and P measurements as columns), containing the reference
#' univariate functional data that must be used to compute the band depths of
#' elements in \code{Data_target}. If \code{Data_target} is \code{fData}, it
#' must be of class \code{fData}
#'
#'
#' @return The function returns a vector containing the BD of elements in
#' \code{Data_target} with respect to elements in \code{Data_reference}.
#'
#' @seealso \code{\link{BD}}, \code{\link{MBD}}, \code{\link{MBD_relative}},
#' \code{\link{fData}}
#'
#' @examples
#'
#' grid = seq( 0, 1, length.out = 1e2 )
#'
#' Data_ref = matrix( c( 0  + sin( 2 * pi * grid ),
#'                       1  + sin( 2 * pi * grid ),
#'                       -1 + sin( 2 * pi * grid ) ),
#'                    nrow = 3, ncol = length( grid ), byrow = TRUE )
#'
#' Data_test_1 = matrix( c( 0.6 + sin( 2 * pi * grid ) ),
#'                       nrow = 1, ncol = length( grid ), byrow = TRUE )
#'
#' Data_test_2 = matrix( c( 0.6 + sin( 2 * pi * grid ) ),
#'                       nrow = length( grid ), ncol = 1, byrow = TRUE )
#'
#' Data_test_3 = 0.6 + sin( 2 * pi * grid )
#'
#' Data_test_4 = array( 0.6 + sin( 2 * pi * grid ), dim = length( grid ) )
#'
#' Data_test_5 = array( 0.6 + sin( 2 * pi * grid ), dim = c( 1, length( grid ) ) )
#'
#' Data_test_6 = array( 0.6 + sin( 2 * pi * grid ), dim = c( length( grid ), 1 ) )
#'
#' Data_test_7 = matrix( c( 0.5  + sin( 2 * pi * grid ),
#'                          -0.5 + sin( 2 * pi * grid ),
#'                          1.1 + sin( 2 * pi * grid ) ),
#'                       nrow = 3, ncol = length( grid ), byrow = TRUE )
#'
#' fD_ref = fData( grid, Data_ref )
#' fD_test_1 = fData( grid, Data_test_1 )
#' fD_test_2 = fData( grid, Data_test_2 )
#' fD_test_3 = fData( grid, Data_test_3 )
#' fD_test_4 = fData( grid, Data_test_4 )
#' fD_test_5 = fData( grid, Data_test_5 )
#' fD_test_6 = fData( grid, Data_test_6 )
#' fD_test_7 = fData( grid, Data_test_7 )
#'
#' BD_relative( fD_test_1, fD_ref )
#' BD_relative( Data_test_1, Data_ref )
#'
#' BD_relative( fD_test_2, fD_ref )
#' BD_relative( Data_test_2, Data_ref )
#'
#' BD_relative( fD_test_3, fD_ref )
#' BD_relative( Data_test_3, Data_ref )
#'
#' BD_relative( fD_test_4, fD_ref )
#' BD_relative( Data_test_4, Data_ref )
#'
#' BD_relative( fD_test_5, fD_ref )
#' BD_relative( Data_test_5, Data_ref )
#'
#' BD_relative( fD_test_6, fD_ref )
#' BD_relative( Data_test_6, Data_ref )
#'
#' BD_relative( fD_test_7, fD_ref )
#' BD_relative( Data_test_7, Data_ref )
#'
#' @export
#'
BD_relative = function( Data_target, Data_reference )
{
  if( class( Data_target ) != class( Data_reference ) )
  {
    if( ! ( class( Data_target ) %in% c( 'numeric', 'array', 'matrix' ) &
            class( Data_reference ) %in% c( 'numeric', 'array', 'matrix' ) ) )
    {
      stop( 'Error in BD_relative: you have to provide target and reference data
          with the same class')
    }
  }

  UseMethod( 'BD_relative', Data_target )
}

#' @rdname BD_relative
#'
#' @aliases BD_relative
#'
#' @export
BD_relative.fData = function( Data_target, Data_reference )
{
  Data_target = Data_target$values
  Data_reference = Data_reference$values
  NextMethod()
}

#' @rdname BD_relative
#'
#' @aliases BD_relative
#'
#' @export
BD_relative.default = function( Data_target, Data_reference )
{
  # Observations
  N = nrow( Data_reference )

  # Number of columns
  P = ncol( Data_reference )

  Data_target = toRowMatrixForm( Data_target )

  if( ncol( Data_target ) != P ) stop( 'Error: you provided a Data_target with
not compliant dimensions to BD_relative')

  # Number of target_functions
  N_target = nrow( Data_target )

  Depths = rep( 0, N_target );

  for( iObs in 1 : N_target )
  {
    rk = apply( rbind( Data_target[ iObs, ], Data_reference ),
                2, rank, ties.method = 'average' )

    N_a = N + 1 - max( rk[ 1, ] )

    N_b = min( rk[ 1, ] ) - 1

    Depths[ iObs ] = ( N_a * N_b ) / ( N * ( N - 1 ) / 2 )
  }

  return( Depths )
}

#' Modified Band Depth for univariate functional data
#'
#' This function computes the Modified Band Depth (MBD) of elements of a
#' functional dataset.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I= [a,b]},
#' this function computes the sample MBD of each element with respect to the
#' other elements of the dataset, i.e.:
#'
#' \deqn{ MBD( X( t ) ) = {N \choose 2 }^{-1} \sum_{1 \leq i_1 < i_2 \leq N}
#' \tilde{\lambda}\big( {t : \min( X_{i_1}(t), X_{i_2}(t) ) \leq X(t) \leq
#' \max( X_{i_1}(t), X_{i_2}(t) ) } \big), }
#'
#' where \eqn{\tilde{\lambda}(\cdot)} is the normalised Lebesgue measure over
#' \eqn{I=[a,b]}, that is \eqn{\tilde{\lambda(A)} = \lambda( A ) / ( b - a )}.
#'
#' See the References section for more details.
#'
#'
#' @param Data either a \code{fData} object or a matrix-like dataset of functional
#' data (e.g. \code{fData$values}), with observations as rows and measurements
#' over grid points as columns.
#' @param manage_ties a logical flag specifying whether a check for ties and
#' relative treatment must be carried out or not (default is \code{FALSE}).
#'
#' @return The function returns a vector containing the values of MBD for the
#' given dataset.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2009). On the Concept of Depth for Functional
#' Data, \emph{Journal of the American Statistical Association}, 104, 718-734.
#'
#' Lopez-Pintado, S. and Romo. J. (2007). Depth-based inference for functional
#' data, \emph{Computational Statistics & Data Analysis} 51, 4957-4968.
#'
#' @seealso \code{\link{BD}}, \code{\link{MBD_relative}},
#' \code{\link{BD_relative}}, \code{\link{fData}}
#'
#' @examples
#' grid = seq( 0, 1, length.out = 1e2 )
#'
#'
#' D = matrix( c( 1 + sin( 2 * pi * grid ),
#'                0 + sin( 4 * pi * grid ),
#'                1 - sin( pi * ( grid - 0.2 ) ),
#'                0.1 + cos( 2 * pi * grid ),
#'                0.5 + sin( 3 * pi + grid ),
#'                -2 + sin( pi * grid ) ),
#'             nrow = 6, ncol = length( grid ), byrow = TRUE )
#'
#' fD = fData( grid, D )
#'
#' MBD( fD )
#'
#' MBD( D )
#'
#' @export
MBD = function( Data, manage_ties = FALSE )
{
  UseMethod( 'MBD', Data )
}

#' @rdname MBD
#'
#' @aliases MBD
#'
#' @export
MBD.fData = function( Data, manage_ties = FALSE )
{
  Data = Data$values
  NextMethod()
}


#' @rdname MBD
#'
#' @aliases MBD
#'
#' @export
MBD.default = function( Data, manage_ties = FALSE )
{
  # Number of rows
  N = nrow( Data )

  # Number of time points in the discrete grid
  P = ncol( Data )

  if( manage_ties )
  {
    # Compute ranks of matrix-like representation of data with `min'
    # tie-breaking rule
    rk_min = apply( Data, 2, function( v )( rank( v, ties.method = 'min' ) ) )

    # Compute ranks of matrix-like representation of data with `max'
    # tie-breaking rule
    rk_max = apply( Data, 2, function( v )( rank( v, ties.method = 'max' ) ) )

    # Times each function value is repeated in the dataset, for each time point
    # ( matrix is N x P)
    Repetitions = rk_max - rk_min + 1

    # Actual number of function values strictly above
    N_a = N - rk_max

    # Actual number of function values strictly below
    N_b = rk_min - 1

    # if( any( Repetitions > 1 ) )
    # {
      # Now, owing to repetitions we gain more bands containing each functional
      # observation.
      # In standard cases (no coincident values) we have N - 1, in general we
      # have N - 1 + N - 2 + ... + N - K, where K is the number of repetitions
      # of the function value in the dataset, time by time.
      # If you split the sum and use Gauss formula, you can get the following
      # expression

      added_bands = N * Repetitions - 0.5 * ( Repetitions * Repetitions +
                                                Repetitions )

      Depths = rowSums( N_a * N_b  + added_bands ) / ( P * ( N - 1 ) * N / 2 )

    } else {

      rk =  apply( Data, 2, function( v )( rank( v ) )  )

      Depths = ( rowSums( ( N - rk ) * ( rk - 1 ) ) / P +
                   N - 1 ) / ( N * ( N - 1 ) / 2 )

    }

  return( Depths )
}

#' Relative Modified Band Depth of functions in a univariate functional dataset
#'
#' This function computes Modified Band Depth (BD) of elements of a univariate
#' functional dataset with respect to another univariate functional dataset.
#'
#' Given a univariate functional dataset of elements \eqn{X_1(t), X_2(t),
#' \ldots, X_N(t)}, and another univariate functional dataset of elements
#' \eqn{Y_1(t), Y_2(t) \ldots, Y_M(t)}, defined over the same compact interval
#' \eqn{I=[a,b]}, this function computes the MBD of
#' elements of the former with respect to elements of the latter, i.e.:
#'
#'\deqn{ MBD( X_i( t ) ) = {M \choose 2 }^{-1} \sum_{1 \leq i_1 < i_2 \leq M}
#' \tilde{\lambda}\big( {t : \min( Y_{i_1}(t), Y_{i_2}(t) ) \leq X_i(t) \leq
#' \max( Y_{i_1}(t), Y_{i_2}(t) ) } \big),}
#'
#' \eqn{\forall i = 1, \ldots, N}, where \eqn{\tilde{\lambda}(\cdot)} is the
#' normalised Lebesgue measure over \eqn{I=[a,b]}, that is
#' \eqn{\tilde{\lambda(A)} = \lambda( A ) / ( b - a )}.
#'
#' @param Data_target is the univariate functional dataset, provided either as
#' an \code{fData} object or in matrix
#' form (N observations as rows and P measurements as columns), whose MBD
#' have to be computed with respect to the reference dataset.
#' @param Data_reference is the dataset, provided either as an \code{fData}
#' object or in matrix form (N observations
#' as rows and P measurements as columns), containing the reference
#' univariate functional data that must be used to compute the MBD of
#' elements in \code{Data_target}. If \code{Data_target} is \code{fData}, it
#' must be of class \code{fData}.
#'
#'
#' @return The function returns a vector containing the MBD of elements in
#' \code{Data_target} with respect to elements in \code{Data_reference}.
#'
#' @seealso \code{\link{MBD}}, \code{\link{BD}}, \code{\link{BD_relative}},
#' \code{\link{fData}}
#'
#' @examples
#'
#' grid = seq( 0, 1, length.out = 1e2 )
#'
#' Data_ref = matrix( c( 0  + sin( 2 * pi * grid ),
#'                       1  + sin( 2 * pi * grid ),
#'                       -1 + sin( 2 * pi * grid ) ),
#'                    nrow = 3, ncol = length( grid ), byrow = TRUE )
#'
#' Data_test_1 = matrix( c( 0.6 + sin( 2 * pi * grid ) ),
#'                       nrow = 1, ncol = length( grid ), byrow = TRUE )
#'
#' Data_test_2 = matrix( c( 0.6 + sin( 2 * pi * grid ) ),
#'                       nrow = length( grid ), ncol = 1, byrow = TRUE )
#'
#' Data_test_3 = 0.6 + sin( 2 * pi * grid )
#'
#' Data_test_4 = array( 0.6 + sin( 2 * pi * grid ), dim = length( grid ) )
#'
#' Data_test_5 = array( 0.6 + sin( 2 * pi * grid ), dim = c( 1, length( grid ) ) )
#'
#' Data_test_6 = array( 0.6 + sin( 2 * pi * grid ), dim = c( length( grid ), 1 ) )
#'
#' Data_test_7 = matrix( c( 0.5  + sin( 2 * pi * grid ),
#'                          -0.5 + sin( 2 * pi * grid ),
#'                          1.1 + sin( 2 * pi * grid ) ),
#'                       nrow = 3, ncol = length( grid ), byrow = TRUE )
#'
#' fD_ref = fData( grid, Data_ref )
#' fD_test_1 = fData( grid, Data_test_1 )
#' fD_test_2 = fData( grid, Data_test_2 )
#' fD_test_3 = fData( grid, Data_test_3 )
#' fD_test_4 = fData( grid, Data_test_4 )
#' fD_test_5 = fData( grid, Data_test_5 )
#' fD_test_6 = fData( grid, Data_test_6 )
#' fD_test_7 = fData( grid, Data_test_7 )
#'
#' MBD_relative( fD_test_1, fD_ref )
#' MBD_relative( Data_test_1, Data_ref )
#'
#' MBD_relative( fD_test_2, fD_ref )
#' MBD_relative( Data_test_2, Data_ref )
#'
#' MBD_relative( fD_test_3, fD_ref )
#' MBD_relative( Data_test_3, Data_ref )
#'
#' MBD_relative( fD_test_4, fD_ref )
#' MBD_relative( Data_test_4, Data_ref )
#'
#' MBD_relative( fD_test_5, fD_ref )
#' MBD_relative( Data_test_5, Data_ref )
#'
#' MBD_relative( fD_test_6, fD_ref )
#' MBD_relative( Data_test_6, Data_ref )
#'
#' MBD_relative( fD_test_7, fD_ref )
#' MBD_relative( Data_test_7, Data_ref )
#'
#'
#' @export
#'
MBD_relative = function( Data_target, Data_reference )
{
  if( class( Data_target ) != class( Data_reference ) )
  {
    if( ! ( class( Data_target ) %in% c( 'numeric', 'array', 'matrix' ) &
            class( Data_reference ) %in% c( 'numeric', 'array', 'matrix' ) ) )
    {
      stop( 'Error in MBD_relative: you have to provide target and reference data
            with the same class')
    }
    }


  UseMethod( 'MBD_relative', Data_target )
  }

#' @rdname MBD_relative
#'
#' @aliases MBD_relative
#'
#' @export
MBD_relative.fData = function( Data_target, Data_reference )
{
  Data_target = Data_target$values
  Data_reference = Data_reference$values
  NextMethod()
}

#' @rdname MBD_relative
#'
#' @aliases MBD_relative
#'
#' @export
MBD_relative.default = function( Data_target, Data_reference )
{
  # Observations
  N = nrow( Data_reference )

  # Number of columns
  P = ncol( Data_reference )

  Data_target = toRowMatrixForm( Data_target )

  if( ncol( Data_target ) != P ) stop( 'Error: you provided a Data_target with
not compliant dimensions to MBD_relative')

  # Number of target_functions
  N_target = nrow( Data_target )

  Depths = rep( 0, N_target );

  for( iObs in 1 : N_target )
  {
    rk = apply( rbind( Data_target[ iObs, ], Data_reference ),
                2, rank, ties.method = 'average' );

    N_a = N + 1 - rk[ 1, ]

    N_b = rk[ 1, ] - 1

    Depths[ iObs ] = sum( N_a * N_b ) / (  P *  N * ( N - 1 ) / 2 )
  }

  return( Depths )
}
