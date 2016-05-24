#' Epigraph Index of univariate functional dataset
#'
#' This function computes the Epigraphic Index (EI) of elements of a univariate
#' functional dataste.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I=[a,b]}, this function computes the
#' EI, i.e.:
#'
#' \deqn{EI( X(t) ) = \frac{1}{N} \sum_{i=1}^N I( G( X(t) ) \subset
#' epi( X_i(t) ) ) = \frac{1}{N} \sum_{i=1}^N I( X(t) \geq X_i(t), \ \
#' \forall t \in I), }
#'
#' where \eqn{G(X(t))} indicates the graph of \eqn{X(t)}, \eqn{epi( X_i(t))}
#' indicates the epigraph of \eqn{X_i(t)}.
#'
#' @param Data either an \code{fData} object or a matrix-like dataset of
#' functional data (e.g. \code{fData$values}), with observations as rows and
#' measurements over grid points as columns.
#'
#' @return The function returns a vector containing the values of EI for each
#' element of the functional dataset provided in \code{Data}.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2012). A half-region depth for functional
#' data, \emph{Computational Statistics and Data Analysis}, 55, 1679-1695.
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and
#' visualization for functional data: the outliergram, \emph{Biostatistics},
#' 15(4), 603-619.
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 2 * pi * grid ),
#'                              C )
#' fD = fData( grid, Data )
#'
#' EI( fD )
#'
#' EI( Data )
#'
#' @seealso \code{\link{MEI}}, \code{\link{HI}}, \code{\link{MHI}},
#' \code{\link{fData}}
#'
#' @export
#'
EI = function( Data )
{
  UseMethod( 'EI', Data )
}

#' @rdname EI
#'
#' @aliases EI
#'
#' @export
EI.fData = function( Data )
{
  Data = Data$values
  NextMethod()
}

#' @rdname EI
#'
#' @aliases EI
#'
#' @export
#'
EI.default = function( Data )
{
  # Number of observations
  N = nrow( Data )

  rk = apply( Data, 2, function( v )( rank( v, ties.method = 'min' ) ) )

  N_a = N - apply( rk, 1, max ) + 1

  EI = N_a / N

  return( EI )
}


#' Modified Epigraph Index of univariate functional dataset
#'
#' This function computes the Modified Epigraphic Index (MEI) of elements of a
#' univariate functional dataste.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I=[a,b]}, this function computes the
#' MEI, i.e.:
#'
#' \deqn{MEI( X(t) ) = \frac{1}{N} \sum_{i=1}^N \tilde{\lambda}( X(t) \geq
#' X_i(t) ), }
#'
#' where \eqn{\tilde{\lambda}(\cdot)} is the normalised Lebesgue measure over
#' \eqn{I=[a,b]}, that is \eqn{\tilde{\lambda(A)} = \lambda( A ) / ( b - a )}.
#'
#' @param Data either an \code{fData} object or a matrix-like dataset of
#' functional data (e.g. \code{fData$values}),
#' with observations as rows and measurements over grid points as columns.
#'
#' @return The function returns a vector containing the values of MEI for each
#' element of the functional dataset provided in \code{Data}.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2012). A half-region depth for functional
#' data, \emph{Computational Statistics and Data Analysis}, 55, 1679-1695.
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and
#' visualization for functional data: the outliergram, \emph{Biostatistics},
#' 15(4), 603-619.
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 2 * pi * grid ),
#'                              C )
#'
#' fD = fData( grid, Data )
#'
#' MEI( fD )
#'
#' MEI( Data )
#'
#' @seealso \code{\link{EI}}, \code{\link{MHI}}, \code{\link{HI}},
#' \code{\link{fData}}
#'
#' @export
#'
MEI = function( Data )
{
  UseMethod( 'MEI', Data )
}

#' @rdname MEI
#'
#' @aliases MEI
#'
#' @export
MEI.fData = function( Data )
{
  Data = Data$values
  NextMethod()
}


#' @rdname MEI
#'
#' @aliases MEI
#'
#' @export
MEI.default = function( Data )
{
  # Number of observations
  N = nrow( Data )

  # Number of time points
  P = ncol( Data )

  # Matrix of ranks, with `min' policy for breaking ties
  rk = apply( Data, 2, function( v ) ( rank( v, ties.method = 'min' ) ) )

  # Number of curves equal or above, time by time
  N_a = N - rk + 1

  MEI = rowSums( N_a ) / ( N * P )

  return( MEI )
}



#' Hypograph Index of univariate functional dataset
#'
#' This function computes the Hypograph Index (HI) of elements of a univariate
#' functional dataste.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I=[a,b]}, this function computes the
#' HI, i.e.:
#'
#' \deqn{EI( X(t) ) = \frac{1}{N} \sum_{i=1}^N I( G( X(t) ) \subset
#' hypo( X_i(t) ) ) = \frac{1}{N} \sum_{i=1}^N I( X(t) \leq X_i(t), \ \
#' \forall t \in I), }
#'
#' where \eqn{G(X(t))} indicates the graph of \eqn{X(t)}, \eqn{epi( X_i(t))}
#' indicates the hypograph of \eqn{X_i(t)}.
#'
#' @param Data either an \code{fData} object or a matrix-like dataset of
#' functional data (e.g. \code{fData$values}),
#' with observations as rows and measurements over grid points as columns.
#'
#' @return The function returns a vector containing the values of HI for each
#' element of the functional dataset provided in \code{Data}.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2012). A half-region depth for functional
#' data, \emph{Computational Statistics and Data Analysis}, 55, 1679-1695.
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and
#' visualization for functional data: the outliergram, \emph{Biostatistics},
#' 15(4), 603-619.
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 2 * pi * grid ),
#'                              C )
#' fD = fData( grid, Data )
#'
#' HI( fD )
#'
#' HI( Data )
#'
#' @seealso \code{\link{MHI}}, \code{\link{EI}}, \code{\link{MEI}},
#' \code{\link{fData}}
#'
#' @export
#'
HI = function( Data )
{
  UseMethod( 'HI', Data )
}


#' @rdname HI
#'
#' @aliases HI
#'
#' @export
HI.fData = function( Data )
{
  Data = Data$values
  NextMethod()
}


#' @rdname HI
#'
#' @aliases HI
#'
#' @export
HI.default = function( Data )
{
  # Number of observations
  N = nrow( Data )

  rk = apply( Data, 2, function( v )( rank( v, ties.method = 'max' ) ) )

  N_b = apply( rk, 1, min )

  HI = N_b / N

  return( HI )
}


#' Modified Hypograph Index of univariate functional dataset
#'
#' This function computes the Modified Hypograph Index (MEI) of elements of a
#' univariate functional dataste.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I=[a,b]}, this function computes the
#' MHI, i.e.:
#'
#' \deqn{MHI( X(t) ) = \frac{1}{N} \sum_{i=1}^N \tilde{\lambda}( X(t) \geq
#' X_i(t) ), }
#'
#' where \eqn{\tilde{\lambda}(\cdot)} is the normalised Lebesgue measure over
#' \eqn{I=[a,b]}, that is \eqn{\tilde{\lambda(A)} = \lambda( A ) / ( b - a )}.
#'
#' @param Data either an \code{fData} object or a matrix-like dataset of
#' functional data (e.g. \code{fData$values}),
#' with observations as rows and measurements over grid points as columns.
#'
#' @return The function returns a vector containing the values of MHI for each
#' element of the functional dataset provided in \code{Data}.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2012). A half-region depth for functional
#' data, \emph{Computational Statistics and Data Analysis}, 55, 1679-1695.
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and
#' visualization for functional data: the outliergram, \emph{Biostatistics},
#' 15(4), 603-619.
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 2 * pi * grid ),
#'                              C )
#' fD = fData( grid, Data )
#'
#' MHI( fD )
#'
#' MHI( Data )
#'
#' @seealso \code{\link{HI}}, \code{\link{MEI}}, \code{\link{EI}},
#' \code{\link{fData}}
#'
#' @export
#'
MHI = function( Data )
{
  UseMethod( 'MHI', Data )
}

#' @rdname MHI
#'
#' @aliases MHI
#'
#' @export
MHI.fData = function( Data )
{
  Data = Data$values
  NextMethod()
}

#' @rdname MHI
#'
#' @aliases MHI
#'
#' @export
MHI.default = function( Data )
{
  # Number of observations
  N = nrow( Data )

  # Number of time points
  P = ncol( Data )

  # Matrix of ranks, with `max' policy for breaking ties
  rk = apply( Data, 2, function( v ) ( rank( v, ties.method = 'max' ) ) )

  # Number of curves equal or below, time by time
  # N_b = rk

  MHI = rowSums( rk ) / ( N * P )

  return( MHI )
}


#' Half-Region Depth for univariate functional data
#'
#' This function computes the Half-Region Depth (HRD) of elements of a univariate
#' functional dataset.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I=[a,b]}, this function computes the HRD
#' of its elements, i.e.:
#'
#' \deqn{HRD(X(t)) = \min( EI( X(t) ), HI(X(t)) ),}
#'
#' where \eqn{EI(X(t))} indicates the Epigraph Index (EI) of \eqn{X(t)} with
#' respect to the dataset, and \eqn{HI(X(t))} indicates the Hypograph Index of
#' \eqn{X(t)} with respect to the dataset.
#'
#' @param Data either an \code{fData} object or a matrix-like dataset of
#' functional data (e.g. \code{fData$values}),
#' with observations as rows and measurements over grid points as columns.
#'
#' @return The function returns a vector containing the values of HRD for each
#' element of the functional dataset provided in \code{Data}.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2012). A half-region depth for functional
#' data, \emph{Computational Statistics and Data Analysis}, 55, 1679-1695.
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and
#' visualization for functional data: the outliergram, \emph{Biostatistics},
#' 15(4), 603-619.
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 2 * pi * grid ),
#'                              C )
#'
#' fD = fData( grid, Data )
#'
#' HRD( fD )
#'
#' HRD( Data )
#'
#' @seealso \code{\link{MHRD}}, \code{\link{EI}}, \code{\link{HI}},
#' \code{\link{fData}}
#'
#' @export
#'
HRD = function( Data )
{
  UseMethod( 'HRD', Data )
}

#' @rdname HRD
#'
#' @aliases HRD
#'
#' @export
HRD.fData = function( Data )
{
  Data = Data$values
  NextMethod()
}


#' @rdname HRD
#'
#' @aliases HRD
#'
#' @export
HRD.default = function( Data )
{
  ei = EI( Data )

  hi = HI( Data )

  return( mapply( min, ei, hi ) )
}


#' Modified Half-Region Depth for univariate functional data
#'
#' This function computes the Modified Half-Region Depth (MHRD) of elements of
#' a univariate functional dataset.
#'
#' Given a univariate functional dataset, \eqn{X_1(t), X_2(t), \ldots, X_N(t)},
#' defined over a compact interval \eqn{I=[a,b]}, this function computes the MHRD
#' of its elements, i.e.:
#'
#' \deqn{MHRD(X(t)) = \min( MEI( X(t) ), MHI(X(t)) ),}
#'
#' where \eqn{MEI(X(t))} indicates the Modified Epigraph Index (MEI) of
#' \eqn{X(t)} with respect to the dataset, and \eqn{MHI(X(t))} indicates the
#' Modified Hypograph Index of \eqn{X(t)} with respect to the dataset.
#'
#' @param Data either an \code{fData} object or a matrix-like dataset of
#' functional data (e.g. \code{fData$values}),
#' with observations as rows and measurements over grid points as columns.
#'
#' @return The function returns a vector containing the values of MHRD for each
#' element of the functional dataset provided in \code{Data}.
#'
#' @references
#'
#' Lopez-Pintado, S. and Romo, J. (2012). A half-region depth for functional
#' data, \emph{Computational Statistics and Data Analysis}, 55, 1679-1695.
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and
#' visualization for functional data: the outliergram, \emph{Biostatistics},
#' 15(4), 603-619.
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 2 * pi * grid ),
#'                              C )
#' fD = fData( grid, Data )
#'
#' MHRD( fD )
#'
#' MHRD( Data )
#'
#' @seealso \code{\link{HRD}}, \code{\link{MEI}}, \code{\link{MHI}}
#'
#' @export
#'
MHRD = function( Data )
{
  UseMethod( 'MHRD', Data )
}


#' @rdname MHRD
#'
#' @aliases MHRD
#'
#' @export
MHRD.fData = function( Data )
{
  Data = Data$values
  NextMethod()
}


#' @rdname MHRD
#'
#' @aliases MHRD
#'
#' @export
MHRD.default = function( Data )
{
  mei = MEI( Data )

  mhi = MHI( Data )

  return( mapply( min, mei, mhi ) )
}
