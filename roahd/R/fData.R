
#' \code{S3} Class for univariate functional datasets.
#'
#'  This function implements a constructor for elements of \code{S3} class
#'  \code{fData}, aimed at implementing a representation of a functional
#'  dataset.
#'
#'  The functional dataset is represented as a collection of measurement of the
#'  observations on an evenly spaced, 1D grid of discrete points (representing,
#'  e.g. time), namely, for functional data defined over a grid \eqn{[t_0,
#'  t_1, \ldots, t_{P-1}]}:
#'
#'  \deqn{ f_{i,j} = f_i( t_0 + j h ), \quad h =  \frac{t_P - t_0}{N},
#'  \quad \forall j = 1, \ldots, P, \quad \forall i = 1, \ldots
#'  N.}
#'
#' @param grid the evenly spaced grid over which the functional observations are
#' measured. It must be a numeric vector of length \code{P}.
#' @param values the values of the observations in the functional dataset,
#' prodived in form of a 2D data structure (e.g. matrix or array) having as
#' rows the observations and as columns their measurements over the 1D grid of
#' length \code{P} specified in \code{grid}.
#'
#' @return The function returns a \code{S3} object of class \code{fData}, containing
#' the following elements:
#' \itemize{
#'  \item{"\code{N}"}{: the number of elements in the dataset;}
#'  \item{"\code{P}"}{: the number of points in the 1D grid over which elements
#'  are measured;}
#'  \item{"\code{t0}"}{: the starting point of the 1D grid;}
#'  \item{"\code{tP}"}{: the ending point of the 1D grid;}
#'  \item{"\code{values}"}{: the matrix of measurements of the functional
#'  observations on the 1D grid provided with \code{grid}.}
#' }
#'
#' @seealso \code{\link{generate_gauss_fdata}}, \code{\link{sub-.fData}}
#'
#' @examples
#' # Defining parameters
#' N = 20
#' P = 1e2
#'
#' # One dimensional grid
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating an exponential covariance function (see related help for more
#' # information )
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Generating a synthetic dataset with a gaussian distribution and
#' # required mean and covariance function:
#' values = generate_gauss_fdata( N,
#'                                centerline = sin( 2 * pi * grid ),
#'                                Cov = C )
#'
#' fD = fData( grid, values )
#'
#' @export
#'
fData = function( grid, values )
{
  all( abs( diff( unique( diff( grid ) ) ) ) < 1e-14 ) ||
    stop( ' Error in fData: you provided an unevenly spaced grid')

  P = length( grid )

  values = toRowMatrixForm( values )

  if( P != ncol( values ) )
  {
    stop( ' Error in fData: you provided mismatching grid and data structure
          of values\n')
  }

  h = grid[ 2 ] - grid[ 1 ]

  return( structure( list( t0 = grid[ 1 ],
                           tP = grid[ P ],
                           h = h,
                           P = P,
                           N = nrow( values ),
                           values = values ),
                     class = c( 'fData' ) ) )
}

#' Specialised method to plot \code{fData} objects
#'
#' This function performs the plot of a functional univariate dataset stored in
#' an object of class \code{fData}. It is able to accept all the usual
#' customisable graphical parameters, otherwise it will use the default ones.
#'
#' @param x the univariate functional dataset in form of \code{fData} object.
#' @param ... additional graphical parameters to be used in plotting functions
#'
#' @seealso \code{\link{fData}}
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' # One dimensional grid
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating an exponential covariance function (see related help for more
#' # information )
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Generating a synthetic dataset with a gaussian distribution and
#' # required mean and covariance function:
#' values = generate_gauss_fdata( N,
#'                                centerline = sin( 2 * pi * grid ),
#'                                Cov = C )
#'
#' fD = fData( grid, values )
#'
#' plot( fD )
#'
#'
#' @export
#'
plot.fData = function( x, ... )
{
  .plot_fData( x, ...  )
}

.plot_fData = function( x,
                        type = 'l', lty = 1,
                        col = fDColorPalette( min( c( x$N,
                                                      30 + x$N %% 30 ) ) ),
                        xlab = '', ylab = '', main = '',
                        ... )
{
  matplot( seq( x$t0, x$tP, length.out = x$P ),
           t( x$values ), type = type, lty = lty,
           col = col, xlab = xlab, ylab = ylab, main = main, ... )
}



#' \code{S3} class for multivariate functional datasets
#'
#' This function implements a constructor for elements of \code{S3} class
#' \code{mfData}, aimed at implementing a representation of a multivariate
#' functional dataset.
#'
#' The functional dataset is represented as a collection of \code{L} components,
#' each one an object of class \code{fData}. Each component must contain elements
#' defined on the same grid as the others, and must contain the same number of
#' elements (\code{N}).
#'
#' @param grid the (evenly spaced) grid over which the functional dataset is
#' defined.
#' @param Data_list a \code{list} containing the \code{L} components of the
#' multivariate functional dataset, defined as 2D data structures (e.g. matrix
#' or array) having as rows the \code{N} observations and as columns the
#' \code{P} measurements on the grid provided by \code{grid}.
#'
#' @return
#' The function returns a \code{S3} object of class \code{mfData}, containing
#' the following elements:
#' \itemize{
#'  \item{"\code{N}"}{: the number of elements in the dataset;}
#'  \item{"\code{L}"}{: the number of components of the functional dataset;}
#'  \item{"\code{P}"}{: the number of points in the 1D grid over which elements
#'  are measured;}
#'  \item{"\code{t0}"}{: the starting point of the 1D grid;}
#'  \item{"\code{tP}"}{: the ending point of the 1D grid;}
#'  \item{"\code{fDList}"}{: the list of \code{fData} objects representing the
#'  \code{L} components as corresponding unviariate functional datasets.}
#' }
#'
#' @seealso \code{\link{fData}}, \code{\link{generate_gauss_fdata}},
#' \code{\link{generate_gauss_mfdata}}
#'
#' @examples
#' # Defining parameters
#' N = 1e2
#'
#' P = 1e3
#'
#' t0 = 0
#' t1 = 1
#'
#' # Defining the measurement grid
#' grid = seq( t0, t1, length.out = P )
#'
#' # Generating an exponential covariance matrix to be used in the simulation of
#' # the functional datasets (see the related help for details)
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Simulating the measurements of two univariate functional datasets with
#' # required center and covariance function
#' Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#' Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#'
#' # Building the mfData object
#' mfData( grid, list( Data_1, Data_2 ) )
#'
#' @export
#'
mfData = function( grid, Data_list )
{
  all( abs( diff( unique( diff( grid ) ) ) ) < 1e-14 ) ||
    stop( ' Error in mfData: you provided an unevenly spaced grid')

  dimMatrix = sapply( Data_list, dim )

  if( any( sapply( dimMatrix, is.null ) ) )
  {
    if( ! all( sapply( dimMatrix, is.null ) ) )
    {
      stop( ' Error in mfData: you provided mismatching datasets as Data_list')
    }

    Data_list = lapply( Data_list, toRowMatrixForm )
  } else if( unique( apply( dimMatrix, 1,
                     function( x )( length( unique( x ) ) ) ) ) != 1 )
  {
    stop( ' Error in mfData: you provided mismatching datasets as Data_list')
  }

  L = length( Data_list )

  if( is.null( names( Data_list ) ) )
  {
    nms = 1 : L
  } else {
    nms = names( Data_list )
  }

  fDList = NULL

  for( nmL in nms )
  {
    fDList[[ nmL ]] = fData( grid, toRowMatrixForm( Data_list[[ nmL ]] ) )
  }

  return( structure( list( N = fDList[[ 1 ]]$N,
                           L = L,
                           P = fDList[[ 1 ]]$P,
                           t0 = fDList[[ 1 ]]$t0,
                           tP = fDList[[ 1 ]]$tP,
                           fDList = fDList ),
                     class = c( 'mfData' ) ) )
}

#' Specialised method to plot \code{mfData} objects
#'
#' This function performs the plot of a functional multivariate dataset stored
#' in an object of class \code{mfData}. It is able to accept all the usual
#' customisable graphical parameters, otherwise it will use the default ones.
#'
#' The current active graphical device is split into a number of sub-figures,
#' each one meant to contain the plot of the corresponding dimension of the
#' \code{mfData} object. In particular, they are arranged in a rectangular
#' lattice with a number of rows equal to \eqn{ \lfloor \sqrt{ L } \rfloor }
#' and a number of columns equal to \eqn{ \lceil L / \lfloor \sqrt{L} \rfloor
#' \rceil }.
#'
#' A special use of the graphical parameters allows to set up y-labels and
#' titles for all the sub-figures in the graphical window. In particular,
#' parameters \code{ylab} and \code{main} can take as argument either a single
#' string, that are repeatedly used for all the sub-graphics, or a list of
#' different strings (one for each of the \code{L} dimensions) that have to be
#' used in the corresponding graphic.
#'
#' @param x the multivariate functional dataset in form of \code{mfData} object.
#' @param ... additional graphical parameters to be used in plotting functions
#' (see \code{Details} for the use of \code{ylab} and \code{main}).
#'
#' @seealso \code{\link{mfData}}, \code{\link{fData}}, \code{\link{plot.fData}}
#'
#' @examples
# Defining parameters
#' N = 1e2
#'
#' P = 1e3
#'
#' t0 = 0
#' t1 = 1
#'
#' # Defining the measurement grid
#' grid = seq( t0, t1, length.out = P )
#'
#' # Generating an exponential covariance matrix to be used in the simulation of
#' # the functional datasets (see the related help for details)
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Simulating the measurements of two univariate functional datasets with
#' # required center and covariance function
#' Data_1 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#' Data_2 = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = C )
#'
#' # Building the mfData object and plotting tt
#' plot( mfData( grid, list( Data_1, Data_2 ) ),
#'       xlab = 'time', ylab = list( '1st dim.', '2nd dim.' ),
#'       main = list( 'An important plot here', 'And another one here' ) )
#'
#' @export
#'
plot.mfData = function( x, ... )
{
  .plot_mfData( x, ... )

}

.plot_mfData = function( x,
                         type = 'l',lty = 1,
                         col = fDColorPalette( min( c( x$N,
                                                       30 + x$N %% 30 ) ) ),
                         xlab = NULL, ylab = NULL, main = NULL,
                         add = FALSE, ... )
{

  if( add == FALSE )
  {
    mfrow_rows = floor( sqrt( x$L ) )
    mfrow_cols = ceiling( x$L / floor( sqrt( x$L ) ) )

    par( mfrow = c( mfrow_rows, mfrow_cols ) )
  }

  if( ! is.null( ylab ) )
  {
    if( length( ylab ) == 1 )
    {
      ylab = rep( ylab, x$L )
    } else if( length( ylab ) != x$L )
    {
      stop( 'Error in plot_mfData_default: you specified a wrong number of y
            labels' )
    }
  } else {
    ylab = rep( list( '' ), x$L )
  }

  if( ! is.null( main ) )
  {
    if( length( main ) == 1 )
    {
      main = rep( main, x$L )
    } else if( length( main ) != x$L )
    {
      stop( 'Error in plot_mfData_default: you specified a wrong number of
            subtitles' )
    }
  } else {
    main = rep( list( '' ), x$L )
  }

  # for( iL in 1 : x$L )
  # {
  #   plot.fData( x$fDList[[ iL ]],
  #               ylab = ylab[[ iL ]],
  #               main = main[[ iL ]],
  #               col = col,
  #               lty = lty,
  #               add = add,
  #               ... )
  #
  #   if( iL < x$L )
  #   {
  #     dev.set( dev.next() )
  #   }
  # }

  plot_aux = function( i, ... )( plot.fData( x$fDList[[ i ]],
                                             ylab = ylab[[ i ]],
                                             main = main[[ i ]],
                                             col = col,
                                             lty = lty,
                                             add = add,
                                             ... ) )

  invisible( sapply( 1 : x$L, plot_aux, ... ) )
}


#' Operator \code{+} and \code{-} for \code{fData} objects
#'
#' These methods provide operators \code{+} and \code{-} to perform sums
#' or differences between an \code{fData} object and either another
#' \code{fData} object or other compliant data structures, like matrices or
#' vectors or arrays, representing the pointwise measurements of the second
#' term of the  sum.
#'
#' If the second term of the operation is an \code{fData} object, it must be
#' defined over the same grid as the first.
#'
#' @param fD the univariate functional dataset in form of \code{fData} object.
#' @param A either an \code{fData} object, defined on the very same grid of
#' \code{fD}, or a 1D data structure (such as 1D array or raw
#' numeric vector), or a 2D data structure (such as 2D array or raw numeric
#' matrix ), that specifies the second term of the sum.
#' In case of a 1D data structure, the sum is performed element-wise between
#' each element of  \code{fD} and \code{A}, and \code{A} must have length
#' \code{P}, size of \code{fD}'s grid.
#' In case of a 2D data structure, the sum is performed element-wise between
#' corresponding elements of \code{fD} and \code{A}'s rows. In this case,
#' \code{A} must have \code{P} columns, as the size of \code{fD}'s grid.
#'
#' @name plus-.fData
#'
#' @return The function returns an \code{fData} object, whose function values
#' have undergone the sum/difference.
#'
#'
NULL


#' @rdname plus-.fData
#'
#' @examples
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = 21, ncol = 10, byrow = TRUE ) )
#' fD + 1 : 10
#'
#' fD + array( 1, dim = c( 1, 10 ) )
#'
#' fD + fD
#'
#' @export
"+.fData" = function( fD, A )
{
  if( class( A ) == 'fData' )
  {
    if( fD$t0 != A$t0 || fD$tP != A$tP || fD$h != A$h || fD$P != A$P )
    {
      stop( 'Error in +.fData: functional data defined over
            mismatching intervals' )
    }

    if( A$N == 1 )
    {
      fD$values = t( t( fD$values ) + as.vector( A$values ) )
    } else {
      fD$values = fD$values + A$values
    }

  } else if( is.null( dim( A ) ) || nrow( A ) == 1 ){

    A = as.vector( A )

    if( length( A ) != fD$P )
    {
      stop( 'Error in +.fData: mismatching arguments')
    }

    fD$values = t( t( fD$values ) +
                     as.vector( rep( A, fD$P / length( A ) ) ) )

  } else if( is.matrix( A ) ) {

    if( ncol( A ) != fD$P || nrow( A ) != fD$N )
    {
      stop('Error in +.fData: mismatching arguments' )
    }

    fD$values = fD$values + A
  }

  return( fD )
}

#' @rdname plus-.fData
#'
#' @examples
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = 21, ncol = 10, byrow = TRUE ) )
#' fD - 2 : 11
#'
#' fD - array( 1, dim = c( 1, 10 ) )
#'
#' fD - fD
#'
#' @export
#'
"-.fData" = function( fD, A )
{
  if( class( A ) == 'fData' )
  {
    if( fD$t0 != A$t0 || fD$tP != A$tP || fD$h != A$h || fD$P != A$P )
    {
      stop( 'Error in -.fData: functional data defined over
            mismatching intervals' )
    }

    if( A$N == 1 )
    {
      fD$values = t( t( fD$values ) - as.vector( A$values ) )
    } else {
      fD$values = fD$values - A$values
    }

  } else if( is.null( dim( A ) ) || nrow( A ) == 1 ){

    A = as.vector( A )

    if( length( A ) != fD$P )
    {
      stop( 'Error in -.fData: mismatching arguments')
    }

    fD$values = t( t( fD$values ) -
                     as.vector( rep( A, fD$P / length( A ) ) ) )

  } else if( is.matrix( A ) ) {

    if( ncol( A ) != fD$P || nrow( A ) != fD$N )
    {
      stop('Error in -.fData: mismatching arguments' )
    }

    fD$values = fD$values - A
  }

  return( fD )
}

#' Operator \code{*} and \code{/} for \code{fData} objects
#'
#' These methods provide operators \code{*} and \code{/} to perform products
#' or divisions between an \code{fData} object and either a number or a
#' compliant 1D data structure, like numeric vector, array or
#' matrix. The operation is computed by performing the element-wise product
#' or division between \code{fD}'s observations and the provided value(s).
#'
#' If the second argument is a 1D data structure, it must have length \code{N}
#' equal to the number of observations in \code{fD}.
#'
#'
#' @param fD the univariate functional dataset in form of \code{fData} object.
#' @param a either a single number or a 1D data structure (such as numeric
#' raw vector, matrix or array) specifying the factor(s) to use in the
#' multiplication/division of \code{fD} elements' values.
#' In the latter case, each factor is used with the corresponding element in
#' \code{fD}, hence a must have length \code{N}, number of observations in
#' \code{fD}.
#'
#' @name times-.fData
#'
#' @return The function returns an \code{fData} object, whose function values
#' have undergone the product/division.
#'
NULL

#' @rdname times-.fData
#'
#' @examples
#'
#' N = 11
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = N, ncol = 10, byrow = TRUE ) )
#' fD * 2
#'
#' fD * seq( 1, N )
#'
#' @export
#'
"*.fData" = function( fD, a )
{
  if( ! 1 %in% dim( as.matrix( a ) ) )
  {
      stop( 'Error in *.fData: dimensions are not compliant' )
  }

  fD$values = fD$values * as.numeric( a )

  return( fD )
}

#' @rdname times-.fData
#'
#' @examples
#'
#' N = 11
#' fD = fData( seq( 0, 1, length.out = 10 ),
#'             values = matrix( seq( 1, 10 ),
#'                              nrow = N, ncol = 10, byrow = TRUE ) )
#' fD / 2
#'
#' fD / rep( 10, N )
#'
#' @export
"/.fData" = function( fD, a )
{
  if( ! 1 %in% dim( as.matrix( a ) ) )
  {
    stop( 'Error in *.fData: dimensions are not compliant' )
  }

  fD$values = fD$values / as.numeric( a )


  return( fD )
}

#' Cross-sectional mean of of a fData object.
#'
#' This \code{S3} method implements the \bold{cross-sectional} mean of a
#' univariate functional dataset stored in a \code{fData} object, i.e. the
#' mean computed point-by-point along the grid over which the dataset is
#' defined.
#'
#' @param x the univariate functional dataset whose cross-sectional mean must be
#' computed, in form of \code{fData} object.
#' @param ... possible additional parameters. This argument is kept for
#' compatibility with the \code{S3} definition of \code{mean}, but it is not
#' actually used.
#'
#' @return The function returns a \code{fData} object with one observation
#' defined on the same grid as the argument \code{x}'s representing the
#' desired cross-sectional mean.
#'
#' @seealso \code{\link{fData}}
#'
#' @examples
#'
#' N = 1e2
#' P = 1e2
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating a gaussian functional sample with desired mean
#' target_mean = sin( 2 * pi * grid )
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.2 )
#' fD = fData( grid, generate_gauss_fdata( N,
#'                                       centerline = target_mean,
#'                                        Cov = C ) )
#'
#' # Graphical representation of the mean
#' plot( fD )
#' plot( mean( fD ), col = 'black', lwd = 2, lty = 2, add = TRUE )
#'
#'
#' @export
mean.fData = function( x, ... )
{
  return( fData( seq( x$t0, x$tP, length.out = x$P ),
                colMeans( x$values ) ) )
}

#' Cross-sectional mean of of a mfData object.
#'
#' This \code{S3} method implements the \bold{cross-sectional} mean of a
#' multivariate functional dataset stored in a \code{mfData} object, i.e. the
#' mean computed point-by-point along the grid over which the dataset is
#' defined.
#'
#' @param x the multivariate functional dataset whose cross-sectional mean must
#' be computed, in form of \code{mfData} object.
#' @param ... possible additional parameters. This argument is kept for
#' compatibility with the \code{S3} definition of \code{mean}, but it is not
#' actually used.
#'
#' @return The function returns a \code{mfData} object with one observation
#' defined on the same grid as the argument \code{x}'s representing the
#' desired cross-sectional mean.
#'
#' @seealso \code{\link{mfData}}
#'
#' @examples
#'
#' N = 1e2
#' L = 3
#' P = 1e2
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating a gaussian functional sample with desired mean
#' target_mean = sin( 2 * pi * grid )
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.2 )
#' # Independent components
#' correlations = c( 0, 0, 0 )
#' mfD = mfData( grid,
#'               generate_gauss_mfdata( N, L,
#'                                      correlations = correlations,
#'                                      centerline = matrix( target_mean,
#'                                                           nrow = 3,
#'                                                           ncol = P,
#'                                                           byrow = TRUE ),
#'                                      listCov = list( C, C, C ) )
#' )
#'
#' # Graphical representation of the mean
#' par( mfrow = c( 1, 3 ) )
#'
#' for( iL in 1 : L )
#' {
#'   plot( mfD$fDList[[ 1 ]] )
#'   plot( mean( mfD )$fDList[[ 1 ]], col = 'black',
#'         lwd = 2, lty = 2, add = TRUE )
#' }
#' @export
mean.mfData = function( x, ... )
{
  return( mfData( seq( x$t0, x$tP, length.out = x$P ),
                  lapply( x$fDList, function( y )( mean( y )$values ) ) ) )
}

#'
#' Median of a univariate functional dataset
#'
#' This method computes the sample median of a univariate functional dataset
#' based on a definition of depth for univariate functional data.
#'
#' Provided a definition of functional depth for univariate data,
#' the corresponding median (i.e. the deepest element of the sample) is returned as the desired median.
#' This method does \bold{not} coincide with the computation of the
#' cross-sectional median of the sample of the point-by-point measurements on
#' the grid. Hence, the sample median is a member of the dataset provided.
#'
#' @param fData the univariate functional dataset whose
#' median is required, in form of \code{fData} object.
#' @param type a string specifying the name of the function defining the depth
#' for univariate data to be used. It must be a valid name of a function defined
#' in the current environment, default is \code{MBD}.
#' @param ... additional parameters to be used in the function specified by
#' argument \code{type}.
#'
#' @return The function returns a \code{fData} object containing the desired
#' sample median.
#'
#' @seealso \code{\link{fData}}, \code{\link{mean.fData}},
#' \code{\link{median_mfData}}
#'
#' @examples
#'
#' N = 1e2
#' P = 1e2
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating a gaussian functional sample with desired mean
#' # Being the distribution symmetric, the sample mean and median are coincident
#' target_median = sin( 2 * pi * grid )
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.2 )
#' fD = fData( grid, generate_gauss_fdata( N,
#'                                       centerline = target_median,
#'                                        Cov = C ) )
#'
#' # Graphical representation of the mean
#' plot( fD )
#' plot( median_fData( fD ), col = 'black', lwd = 2, lty = 2, add = TRUE )
#'
#' @export
median_fData = function( fData, type = 'MBD', ... )
{
  Depths = eval( parse( text = paste( type, '( fData$values, ... )', sep = '' ) ) )

  return( fData( seq( fData$t0, fData$tP, length.out = fData$P ),
                 fData$values[ which.max( Depths ), ] ) )
}


#'
#' Median of a multivariate functional dataset
#'
#' This method computes the sample median of a multivariate functional dataset
#' based on a definition of depth for multivariate functional data.
#'
#' Provided a definition of functional depth for multivariate data,
#' the corresponding median (i.e. the deepest element of the sample) is returned
#' as the desired median.
#' This method does \bold{not} coincide with the computation of the
#' cross-sectional median of the sample of the point-by-point measurements on
#' the grid. Hence, the sample median is a member of the dataset provided.
#'
#' @param mfData the multivariate functional dataset whose
#' median is required, in form of \code{mfData} object.
#' @param type a string specifying the name of the function defining the depth
#' for multivariate data to be used. It must be a valid name of a function
#' defined in the current environment, default is \code{multiMBD}.
#' @param ... additional parameters to be used in the function specified by
#' argument \code{type}.
#'
#' @return The function returns a \code{mfData} object containing the desired
#' sample median.
#'
#' @seealso \code{\link{mfData}}, \code{\link{mean.mfData}},
#' \code{\link{median_fData}}
#'
#' @examples
#'
#' N = 1e2
#' L = 3
#' P = 1e2
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating a gaussian functional sample with desired mean
#' # Being the distribution symmetric, the sample mean and median are coincident
#' target_median = sin( 2 * pi * grid )
#' C = exp_cov_function( grid, alpha = 0.2, beta = 0.2 )
#'
#' # Strongly dependent components
#' correlations = c( 0.9, 0.9, 0.9 )
#' mfD = mfData( grid,
#'               generate_gauss_mfdata( N, L,
#'                                      correlations = correlations,
#'                                      centerline = matrix( target_median,
#'                                                           nrow = 3,
#'                                                           ncol = P,
#'                                                           byrow = TRUE ),
#'                                      listCov = list( C, C, C ) )
#' )
#'
#' med_mfD = median_mfData( mfD, type = 'multiMBD', weights = 'uniform' )
#'
#' # Graphical representation of the mean
#' par( mfrow = c( 1, 3 ) )
#'
#' for( iL in 1 : L )
#' {
#'   plot( mfD$fDList[[ 1 ]] )
#'   plot( med_mfD$fDList[[ 1 ]], col = 'black',
#'         lwd = 2, lty = 2, add = TRUE )
#' }
#' @export
median_mfData = function( mfData, type = 'multiMBD', ... )
{
  Depths = eval( parse( text = paste( type, '( toListOfValues( mfData ), ... )',
                                      sep = '' ) ) )

  return( mfData( seq( mfData$t0, mfData$tP, length.out = mfData$P ),
                  lapply( toListOfValues( mfData ), `[`,
                          which.max( Depths ), ) ) )
}

#' Operator \code{sub-.fData} to subset \code{fData} obejcts
#'
#' This method provides an easy and natural way to subset a functional dataset
#' stored in a \code{fData} object, without having to deal with the inner
#' representation of \code{fData} class.
#'
#' @param fD the univariate functional dataset in form of \code{fData} object.
#' @param i a valid expression to subset rows ( observations ) of the univariate
#' functional dataset
#' @param j a valid expression to subset columns ( measurements over the grid )
#' of the univariate functional dataset.
#' @param as_fData logical flag to specify whether the output should be returned
#' as an \code{fData} object containing the required subset or as a matrix of
#' values, default is \code{TRUE}.
#'
#' @return The method returns either an \code{fData} object ( if \code{as_fData
#' = TRUE } ) or a \code{matrix} ( if \code{as_fData = FALSE } ) containing the
#' required subset ( both in terms  of observations and measurement points ) of
#' the univariate functional dataset.
#'
#' @name sub-.fData
#'
#' @seealso \code{\link{fData}}
#'
#' @examples
#'
#' N = 20
#' P = 1e2
#'
#' # One dimensional grid
#' grid = seq( 0, 1, length.out = P )
#'
#' # Generating an exponential covariance function (see related help for more
#' # information )
#' C = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Generating a synthetic dataset with a gaussian distribution and
#' # required mean and covariance function:
#' fD = fData( grid,
#'             generate_gauss_fdata( N,
#'                                   centerline = sin( 2 * pi * grid ),
#'                                   Cov = C ) )
#'
#' dev.new()
#' par( mfrow = c( 2, 2 ) )
#'
#' # Original data
#' plot( fD )
#'
#' # Subsetting observations
#' plot( fD[ c(1,2,3), , as_fData = TRUE ] )
#'
#' # Subsetting measurements
#' plot( fD[ , 1 : 30 ] )
#'
#' # Subsetting both observations and measurements
#' plot( fD[ 1 : 10, 50 : P ] )
#'
#' # Subsetting both observations and measurements but returning a matrix
#' fD[ 1 : 10, 50 : P, as_fData = FALSE ]
#'
#' @export
"[.fData" = function( fD, i, j, as_fData = TRUE )
{
  if( as_fData == TRUE )
  {
    if( missing( j ) )
    {
      return( structure( list( t0 = fD$t0,
                               tP = fD$tP,
                               h = fD$h,
                               P = fD$P,
                               N = ifelse( missing( i ), fD$N, length( i ) ),
                               values = toRowMatrixForm( fD$values[ i, ] ) ),
                         class = c( 'fData' ) ) )
    } else {
      return( structure( list( t0 = fD$t0 + ( min( j ) - 1 ) * fD$h,
                               tP = fD$t0 + ( max( j ) - 1 ) * fD$h,
                               h = fD$h,
                               P = length( j ),
                               N = ifelse( missing( i ), fD$N, length( i ) ),
                               values = toRowMatrixForm( fD$values[ i, j ] ) ),
                         class = c( 'fData' ) ) )
    }
  } else {
    return( fD$values[ i, j ] )
  }
}


#' Manipulation of \code{mfData} list of values
#'
#' This utility function manipulates a \code{mfData} object in order to extract
#' from the list of its \code{fData} objects ( namely, \code{mfData$fDList} )
#' the measurement values of each component and stores them into a list.
#'
#' Given a \code{mfData} of \code{L} components, the function is equivalent to
#' \code{ list( mfData$fDList[[ 1 ]]$values,} \code{..., }
#' \code{ mfData$fDList[[ L ]]$values ) }.
#'
#' @param mfData the multivariate functional dataset in form of \code{mfData}
#' object.
#'
#' @return  The function returns the list of values of each \code{fData} object
#' representing the components of \code{mfData}.
#'
#' @seealso \code{\link{mfData}}
#'
#' @examples
#'
#' grid = seq( 0, 1, length.out = 5 )
#'
#' D_1 = matrix( 1 : 5, nrow = 10, ncol = 5, byrow = TRUE )
#' D_2 = 2 * D_1
#' D_3 = 3 * D_1
#'
#' mfD = mfData( grid, list( D_1, D_2, D_3 ) )
#' mfD
#'
#' toListOfValues( mfD )
#'
#' @export
#'
toListOfValues = function( mfData )
{
  eval( parse( text = paste( 'list(', paste( 'mfData$fDList[[ ',
                                             1 : mfData$L, ' ]]$values',
                         sep = '', collapse = ', ' ), ')' ) ) )
}

#' Unfolding a univariate functional dataset
#'
#' This function operates on a univariate functional dataset and transforms its
#' observations unfolding their values and turning them into monotone functions.
#'
#' Each function of the \code{fData} object is transformed into a nonmonotone
#' function into a monotone function by ``unfolding'' it at any of its maxima.
#' For more details about the definition of the transform, see the reference.
#'
#' @param fData the unvariate functional dataset in form of \code{fData} object.
#'
#' @return The function returns an \code{fData} object whose observations are
#' the unfolded version of the corresponding observations in the argument
#' \code{fData}.
#'
#' @references Arribas-Gil, A. and Romo, J. (2012). Robust depth-based estimation
#' in the time warping model, \emph{Biostatistics}, 13 (3), 398--414.
#'
#' @seealso \code{\link{fData}}, \code{\link{warp}}
#'
#' @examples
#' P = 1e3
#'
#' time_grid = seq( 0, 1, length.out = P )
#'
#' D = matrix( c( sin( 2 * pi * time_grid ),
#'                cos( 2 * pi * time_grid ),
#'                sin( 10 * pi * time_grid ) * time_grid + 2 ),
#'             ncol = P, nrow = 3, byrow = TRUE )
#'
#' # Functional dataset
#' fD = fData( time_grid, D )
#'
#' # Unfolded version
#' fD_unfold = unfold( fD )
#'
#' dev.new()
#' par( mfrow = c( 1, 2 ) )
#' plot( fD, main = 'Original data' )
#' plot( fD_unfold, main = 'Unfolded data' )
#'
#' @export
unfold = function( fData )
{
  return( fData( seq( fData$t0, fData$tP, length.out = fData$P ),
                 t( apply( fData$values,
                           1,
                           function( x )( c( x[ 1 ],
                                             x[ 1 ] +
                                               cumsum( abs( diff( x ) ) ) ) ) )
                 ) ) )
}

#' Warp elements of a univariate functional dataset
#'
#' This function carries out the warping of elements of a univariate functional
#' dataset by using a set of pre-computed warping functions.
#'
#' Given a univariate functional dataset \eqn{X_1(t), \ldots, X_N(t)} and a set
#' of warping functions \eqn{H_1(t), \ldots, H_N(t)}, such that:
#' \deqn{ H_i : s \longrightarrow t = H_i(s), \quad \forall i = 1, \ldots, N,}
#' where \eqn{s} spans the warped (or registered) grid and \eqn{t} spans the
#' original grid, the function computes the warping given by the following
#' composition:
#' \deqn{ X_1 \circ H_1(t), \ldots, X_N \circ H_N(t).}
#'
#' @param fData the functional dataset whose observations must be warped in
#' form of \code{fData} object.
#' @param warpings the warping functions \eqn{H_1, \ldots, H_N}, in form of
#' \code{fData} object, defined over the registered/warped grid.
#'
#' @return The function returns the univariate functional dataset of warped
#' functions, in form of \code{fData} object.
#'
#' @seealso \code{\link{fData}}
#'
#' @examples
#'
#' set.seed( 1618033 )
#'
#' N = 30
#'
#' t0 = 0
#' t1 = 1
#' P = 1e3 + 1
#'
#' time_grid = seq( t0, t1, length.out = P )
#'
#' means = round( runif( N,
#'                       t0 + (t1 - t0) / 8,
#'                       t1 - (t1 - t0) / 8  ), 3 )
#'
#' Data = matrix( sapply( means,
#'                        function( m )( dnorm( time_grid, mean = m, sd = 0.05 ) ) ),
#'                ncol = P, nrow = N, byrow = TRUE )
#'
#' fD = fData( time_grid, Data )
#'
#' # Piecewise linear warpings
#' template_warping = function( m )( c( time_grid[ time_grid <= 0.5 ] * m / 0.5,
#'                                      ( time_grid[ time_grid > 0.5 ]
#'                                        - 0.5 ) * (1 - m ) / 0.5 + m ) )
#'
#'
#' warpings = matrix( sapply( means, template_warping ),
#'                    ncol = P,
#'                    nrow = N, byrow = TRUE )
#'
#' wfD = fData( time_grid, warpings )
#'
#' fD_warped = warp( fD, wfD )
#'
#' dev.new()
#' par( mfrow = c( 1, 3 ) )
#' plot( fD,
#'      main = 'Unregistered functions', xlab = 'actual grid', ylab = 'values'  )
#' plot( wfD,
#'      main = 'Warping functions', xlab = 'registered grid',
#'      ylab = 'actual grid' )
#' plot( fD_warped,
#'      main = 'Warped functions', xlab = 'registered grid',
#'      ylab = 'values' )
#'
#' @importFrom stats approx
#'
#' @export
warp = function( fData, warpings )
{
  if(
    fData$N != warpings$N
    )
  {
    stop( ' Error in warp: you have to provide a warping for each functional
observations' )
  }

  if( any( diff( warpings$values[ , 1 ] ) > .Machine$double.eps ) |
      any( diff( warpings$values[ , fData$P ] ) > .Machine$double.eps ) )
  {
    stop( ' Error in warp: you have to prescribe warpings with same starting
          point and ending point (at least up to .Machine$double.eps = ',
          .Machine$double.eps, ')' )
  }

  time_grid = seq( fData$t0,
                   fData$tP,
                   length.out = fData$P )

  return( fData( seq( warpings$t0,
                      warpings$tP,
                      length.out = warpings$P ),
                 t( sapply( 1 : fData$N,
                            function( i )(
                              approx( time_grid,
                                      fData$values[ i, ],
                                      xout = warpings$values[ i, ],
                                      yright = fData$values[ i, fData$P ],
                                      yleft = fData$values[ i, 1 ] )$y ) ) ) ) )
}



#' Converting object to \code{mfData} class
#'
#' This S3 method provides a way to convert some objects to the class
#' \code{mfData}, thus obtaining a multivariate functional dataset.
#'
#' @param ... additional parameters.
#'
#' @return The function returns a \code{mfData} object, obtained starting from
#' argument \code{x}.
#'
#' @export
#'
as.mfData = function( x, ... )
{
  UseMethod( 'as.mfData', x )
}

#'@rdname as.mfData
#'
#' @param x a list of univariate functional datasets, provided in form of
#' \code{fData} objects.
#'
#' @examples
#'
#' grid = seq( 0, 1, length.out = 100 )
#'
#' fD_1 = fData( grid, sin( 2 * pi * grid ) )
#' fD_2 = fData( grid, cos( 2 * pi * grid ) )
#'
#' plot( as.mfData( list( fD_1, fD_2 ) ) )
#'
#' @export
as.mfData.list = function( x, ... )
{
  if( any( sapply( x, class ) != 'fData' ) )
  {
    stop( ' Error in as.mfData.list: you have to provide a list of fData
          objects')
  }

  if( any( diff( sapply( x,
                         function( y ) ( y$t0 ) ) ) >= .Machine$double.eps ) |
      any( diff( sapply( x,
                         function( y ) ( y$tP ) ) ) >= .Machine$double.eps ) |
      any( diff( sapply( x,
                         function( y ) ( y$P ) ) ) >= .Machine$double.eps )  |
      any( diff( sapply( x,
                         function( y ) ( y$N ) ) ) >= .Machine$double.eps )  )
  {
    stop( ' Error in as.mfData.list: you have provided fData list elements with
          different parameters')
  }

  return( mfData( seq( x[[ 1 ]]$t0,
                       x[[ 1 ]]$tP,
                       length.out = x[[ 1 ]]$P ),
                  eval( parse( text = paste( 'list(',
                                             paste( 'x[[ ',
                                                    1 : length( x ),
                                                    ' ]]$values',
                                                    sep = '', collapse = ', ' ),
                                             ')' ) ) ) ) )
}
