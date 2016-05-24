#' Functional boxplot of univariate and multivariate functional data
#'
#' This function can be used to perform the functional boxplot of univariate or
#' multivariate functional data.
#'
#' @section Adjustment:
#'
#' In the \bold{univariate functional case}, when the adjustment option is selected,
#' the value of \eqn{F} is optimised for the univariate functional dataset
#' provided with \code{Data}.
#'
#' In practice, a number \code{adjust$N_trials} of times a synthetic population
#' (of size \code{adjust$tiral_size} with the same covariance (robustly
#' estimated from data) and centerline as \code{fData} is simulated without
#' outliers and each time an optimised value \eqn{F_i} is computed so that a
#' given proportion (\code{adjust$TPR}) of observations is flagged as outliers.
#' The final value of \code{F} for the functional boxplot is determined as an
#' average of \eqn{F_1, F_2, \ldots, F_{N_{trials}}}.
#' At each time step the optimisation problem is solved using
#' \code{stats::uniroot} (Brent's method.
#'
#'
#' @param Data the univariate or multivariate functional dataset whose functional
#' boxplot must be determined, in form of \code{fData} or \code{mfData} object.
#' @param Depths either a vector containing the depths for each element of the
#' dataset, or:
#' \itemize{
#' \item{"\emph{univariate case}"}{: a string containing the name of the method
#' you want to use to compute it. The default is \code{'MBD'}};
#' \item{"\emph{multivariate case}"}{: a list with elements \code{def},
#' containing the name of the depth notion to be used to compute depths
#' (\code{BD} or \code{MBD}), and \code{weights}, containing the value
#' of parameter \code{weights} to be passed to the depth function. Default is
#' \code{list( def = 'MBD', weights = 'uniform' ) }. }
#' }
#' In both cases the name of the functions to compute depths must be available
#' in the caller's environment.
#' @param Fvalue the value of the inflation factor \eqn{F}, default is
#' \code{F = 1.5}.
#' @param adjust either \code{FALSE} if you would like the default value for the
#' inflation factor, \eqn{F = 1.5}, to be used, or (for now \bold{only in the
#' univariate functional case}) a list specifying the parameters required by
#' the adjustment:
#'  \itemize{
#'  \item{"\code{N_trials}"}{: the number of repetitions of the adujustment
#'  procedure based on the simulation of a gaussisan population of functional
#'  data, each one producing an adjusted value of \eqn{F}, which will lead
#'  to the averaged adjusted value \eqn{\bar{F}}. Default is 20;}
#'  \item{"\code{trial_size}"}{: the number of elements in the gaussian
#'  population of functional data that will be simulated at each repetition of
#'  the adjustment procedure. Default is \code{Data$N};}
#'  \item{"\code{TPR}"}{: the True Positive Rate of outliers, i.e. the proportion
#'  of observations in a dataset without amplitude outliers that have to be
#'  considered outliers. Default is \code{2 * pnorm( 4 * qnorm( 0.25 ) )};}
#'  \item{"\code{F_min}"}{: the minimum value of \eqn{F}, defining the left
#'  boundary for the optimisation problem aimed at finding, for a given dataset
#'  of simulated gaussian data associated to \code{Data}, the optimal value of
#'  \eqn{F}. Default is 0.5;}
#'  \item{"\code{F_max}"}{: the maximum value of \eqn{F}, defining the right
#'  boundary for the optimisation problem aimed at finding, for a given dataset
#'  of simulated gaussian data associated to \code{Data}, the optimal value of
#'  \eqn{F}. Default is 5;}
#'  \item{"\code{tol}"}{: the tolerance to be used in the optimisation problem
#'  aimed at finding, for a given dataset of simulated gaussian data associated
#'  to \code{Data}, the optimal value of \eqn{F}. Default is \code{1e-3};}
#'  \item{"\code{maxiter}"}{: the maximum number of iterations to solve the
#'  optimisation problem aimed at finding, for a given dataset of simulated
#'  gaussian data associated to \code{Data}, the optimal value of \eqn{F}.
#'  Default is \code{100};}
#'  \item{"\code{VERBOSE}"}{: a parameter controlling the verbosity of the
#'  adjustment process;}
#'  }
#' @param display either a logical value indicating wether you want the
#' outliergram to be displayed, or the number of the graphical device
#' where you want the outliergram to be displayed.
#' @param xlab the label to use on the x axis when displaying the functional
#' boxplot.
#' @param ylab the label (or list of labels for the multivariate functional case)
#' to use on the y axis when displaying the functional boxplot.
#' @param main the main title (or list of titles for the multivariate functional
#' case) to be used when displaying the functional boxplot.
#' @param ... additional graphical parameters to be used in plotting functions.
#'
#' @return
#' Even when used in graphical way to plot the functional boxplot, the function
#' returns a list of three elements: the first, \code{Depths}, contains the depths
#' of each element of the functional dataset; the second, \code{Fvalue}, is the
#' value of F used to obtain the outliers, and the third, \code{ID_out}, contains
#' the vector of indices of dataset's elements flagged as outliers (if any).
#'
#' @references
#'
#'
#'	Sun, Y., & Genton, M. G. (2012). Functional boxplots. Journal of
#'	Computational and Graphical Statistics.
#'
#'	Sun, Y., & Genton, M. G. (2012). Adjusted functional boxplots for spatio-
#'	temporal data visualization and outlier detection. Environmetrics, 23(1), 54-64.
#'
#' @examples
#'
#' # UNIVARIATE FUNCTIONAL BOXPLOT - NO ADJUSTMENT
#'
#' N = 2 * 10 + 1
#' P = 2e2
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' D = matrix( sin( 2 * pi * grid ), nrow = N, ncol = P, byrow = TRUE )
#'
#' D = D + c( 0, 1 : (( N - 1 )/2), -( ( ( N - 1 ) / 2 ) : 1 ) )
#'
#' fD = fData( grid, D )
#'
#' dev.new()
#' par( mfrow = c(1,2) )
#' plot( fD, lwd = 2, main = 'Functional dataset',
#'       xlab = 'time', ylab = 'values' )
#'
#' fbplot( fD, main = 'Functional boxplot', xlab = 'time', ylab = 'values' )
#'
#' # UNIVARIATE FUNCTIONAL BOXPLOT - WITH ADJUSTMENT
#'
#'
#' set.seed( 161803 )
#'
#' P = 2e2
#' grid = seq( 0, 1, length.out = P )
#'
#' N = 1e2
#'
#' # Generating a univariate synthetic gaussian dataset
#' Data = generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ),
#'                              Cov = exp_cov_function( grid,
#'                                                      alpha = 0.3,
#'                                                      beta  = 0.4 ) )
#' fD = fData( grid, Data )
#'
#' dev.new()
#' \dontrun{
#' fbplot( fD, adjust = list( N_trials = 10,
#'                            trial_size = 5 * N,
#'                            VERBOSE = TRUE ),
#'                      xlab = 'time', ylab = 'Values',
#'                      main = 'My adjusted functional boxplot' )
#' }
#'
#' # MULTIVARIATE FUNCTIONAL BOXPLOT - NO ADJUSTMENT
#'
#' set.seed( 1618033 )
#'
#' P = 1e2
#' N = 1e2
#' L = 2
#'
#' grid = seq( 0, 1, length.out = 1e2 )
#'
#' C1 = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#' C2 = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' # Generating a bivariate functional dataset of gaussian data with partially
#' # correlated components
#' Data = generate_gauss_mfdata( N, L,
#'                               centerline = matrix( sin( 2 * pi * grid ),
#'                                                    nrow = 2, ncol = P,
#'                                                    byrow = TRUE ),
#'                               correlations = rep( 0.5, 1 ),
#'                               listCov = list( C1, C2 ) )
#'
#' mfD = mfData( grid, Data )
#'
#' dev.new()
#' fbplot( mfD, Fvalue = 2.5, xlab = 'time', ylab = list( 'Values 1',
#'                                                        'Values 2' ),
#'         main = list( 'First component', 'Second component' ) )
#'
#'
#'
#' @seealso \code{\link{fData}}, \code{\link{MBD}}, \code{\link{BD}},
#' \code{\link{mfData}}, \code{\link{multiMBD}}, \code{\link{multiBD}}
#'
#' @export
#'
fbplot = function( Data,
                   Depths = 'MBD',
                   Fvalue = 1.5,
                   adjust = FALSE,
                   display = TRUE,
                   xlab = NULL,
                   ylab = NULL,
                   main = NULL,
                   ... )
{
  UseMethod( 'fbplot', Data )
}

#' @rdname fbplot
#' @aliases fbplot
#'
#' @importFrom stats quantile
#'
#' @export
#'
fbplot.fData = function( Data,
                         Depths = 'MBD',
                         Fvalue = 1.5,
                         adjust = FALSE,
                         display = TRUE,
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         ... )
{
  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths

    if( Depths_spec == 'MBD' )
    {
      Depths = MBD( Data$values, manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( Data$values )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == Data$N )
  }

  if( ! is.list( adjust ) )
  {
    # Plain functional boxplot with default F value: F = 1.5

    out = .fbplot_fData( time_grid, Data$values, Depths, Fvalue )

  } else {

    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )

    trial_size = ifelse( is.null( adjust$trial_size ),
                         Data$N,
                         adjust$trial_size )

    TPR = ifelse( is.null( adjust$TPR ),
                  2 * pnorm( 4 * qnorm( 0.25 ) ),
                  adjust$TPR )

    F_min = ifelse( is.null( adjust$F_min ),
                    0.5,
                    adjust$F_min )

    F_max= ifelse( is.null( adjust$F_max ),
                   5,
                   adjust$F_max )

    tol = ifelse( is.null( adjust$tol ),
                  1e-3,
                  adjust$tol )

    maxiter = ifelse( is.null( adjust$maxiter ),
                      100,
                      adjust$maxiter )

    VERBOSE = ifelse( is.null( adjust$VERBOSE ),
                      FALSE,
                      adjust$VERBOSE )

    # Estimation of robust covaraince matrix
    Cov = robustbase::covOGK( Data$values, sigmamu = robustbase::s_Qn )$cov

    # Cholesky factor
    CholCov <- chol( Cov )

    # Centerline of the dataset
    centerline = Data$values[ which.max( Depths ), ]

    Fvalues = rep( 0, N_trials )

    cost_functional = function( F_curr )( length(
      .fbplot_fData( time_grid,
                     Data_gauss,
                     Depths = Depths_spec,
                     Fvalue = F_curr )$ID_out ) /
        trial_size - TPR )

    for( iTrial in 1 : N_trials )
    {
      if( VERBOSE > 0 )
      {
        cat( ' * * * Iteration ', iTrial, ' / ', N_trials, '\n' )
      }

      Data_gauss = generate_gauss_fdata( trial_size, centerline, CholCov = CholCov )

      if( VERBOSE > 0 )
      {
        cat( ' * * * * beginning optimisation\n' )
      }

      opt = uniroot( cost_functional,
                     interval = c( F_min, F_max ),
                     tol = tol,
                     maxiter = maxiter )
      if( VERBOSE > 0 )
      {
        cat( ' * * * * optimisation finished.\n')
      }

      Fvalues[ iTrial ] = opt$root
    }

    Fvalue = mean( Fvalues )

    out = .fbplot_fData( time_grid, Data$values, Depths, Fvalue = Fvalue  )
  }

  ID_out = out$ID_out

  # Plotting part
  if( is.numeric( display ) )
  {
    dev.set( display )
  }

  if( ! display == FALSE )
  {
    # Creating color palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( Data$N - length( ID_out ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                    c = 150 )( length( ID_out ) )
    col_envelope = set_alpha( 'blue', alpha = 0.4 )
    col_center = set_alpha( 'blue', alpha = 1 )
    col_fence_structure = 'darkblue'

    time_grid = seq( Data$t0, Data$tP, length.out = Data$P )

    xlab = ifelse( is.null( xlab ), '', xlab )
    ylab = ifelse( is.null( ylab ), '', ylab )
    main = ifelse( is.null( main ), '', main )

    if( length( ID_out ) > 0 )
    {
      # Plotting non-outlying data
      matplot( time_grid,
               t( Data$values[ - ID_out, ] ), lty = 1, type = 'l',
               col = col_non_outlying,
               ylim = range( Data$values ),
               xlab = xlab, ylab = ylab, main = main, ... )

      # Computing maximum and minimum envelope
      max_envelope_limit = apply( Data$values[ - ID_out, ], 2, max )
      min_envelope_limit = apply( Data$values[ - ID_out, ], 2, min )
    } else {
      # Plotting all data
      matplot( time_grid,
               t( Data$values ), lty = 1, type = 'l',
               col = col_non_outlying,
               ylim = range( Data$values ),
               xlab = xlab, ylab = ylab, main = main, ... )

      # Computing maximum and minimum envelope
      max_envelope_limit = apply( Data$values, 2, max )
      min_envelope_limit = apply( Data$values, 2, min )
    }


    # Filling in the central envelope

    polygon( c(time_grid, rev( time_grid) ),
             c( out$min_envelope_central, rev( out$max_envelope_central ) ),
             col = col_envelope, border = NA)
    lines( time_grid, out$max_envelope_central, lty = 1, col = col_envelope, lwd = 3 )
    lines( time_grid, out$min_envelope_central, lty = 1, col = col_envelope, lwd = 3 )

    # Plotting the sample median
    lines( time_grid, Data$values[ which.max( Depths ), ], lty = 1, type = 'l',
           col = col_center, lwd = 3)

    lines( time_grid, max_envelope_limit, lty = 1,
    col = col_fence_structure, lwd = 3 )
    lines( time_grid, min_envelope_limit, lty = 1,
    col = col_fence_structure, lwd = 3 )

    # Plotting vertical whiskers
    half.time_grid = which.min( abs( time_grid - 0.5 ) )
    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$max_envelope_central[ half.time_grid ],
              max_envelope_limit[ half.time_grid ] ),
           lty = 1, col = col_fence_structure, lwd = 3 )

    lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
           c( out$min_envelope_central[ half.time_grid ],
              min_envelope_limit[ half.time_grid ] ),
           lty = 1, col = col_fence_structure, lwd = 3 )

    # Plotting outlying data
    if( length( ID_out ) > 0 )
    {
      matplot( time_grid, t( toRowMatrixForm( Data$values[ ID_out, ] ) ),
               lty = 1, type = 'l', col = col_outlying, lwd = 3, add = T )
    }
  }

  return( list( Depth = Depths,
                Fvalue = Fvalue,
                ID_outliers = ID_out ) )
}


.fbplot_fData = function( time_grid, Data, Depths = 'MBD', Fvalue = 1.5 )
{
  # Number of observations
  N = nrow( Data )

  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths = eval( parse( text = paste( Depths, '( Data )', sep = '' ) ) )
  } else {
    stopifnot( length( Depths ) == N )
  }

  Data_center = Data[ which.max( Depths ), ]

  id_central_region = which( Depths >= quantile( Depths, prob = 0.5 ) )

  max_envelope_central = apply( Data[ id_central_region, ], 2, max )
  min_envelope_central = apply( Data[ id_central_region, ], 2, min )

  fence_upper = ( max_envelope_central - Data_center ) * Fvalue + Data_center
  fence_lower = ( min_envelope_central - Data_center ) * Fvalue + Data_center

  ID_outlying = which( apply( Data, 1, function(x)( any( x > fence_upper  ) |
                                                      any ( x < fence_lower ) ) ) )

  return( list( ID_out = ID_outlying,
                min_envelope_central = min_envelope_central,
                max_envelope_central = max_envelope_central,
                fence_lower = fence_lower,
                fence_upper = fence_upper ) )
}

#' @rdname fbplot
#' @aliases fbplot
#' @export
fbplot.mfData = function( Data,
                          Depths = list( def = 'MBD',
                                         weights = 'uniform' ),
                          Fvalue = 1.5,
                          adjust = FALSE,
                          display = TRUE,
                          xlab = NULL,
                          ylab = NULL,
                          main = NULL,
                          ... )
{
  if( is.list( adjust ) )
  {
    stop( ' Error in fbplot.mfData: for now the adjustment support is not
provided in the multivariate version of the functional boxplot' )
  }

  listOfValues = toListOfValues( Data )

  # Checking if depths have already been provided or must be computed
  if( is.list( Depths ) )
  {
    if( length( Depths ) != 2 & ! identical( names( Depths ),
                                             c( 'def', 'weights' ) ) )
    {
      stop( " Error in fbplot.mfData: you have to provide both a specification
            for the depth definition and one for the set of weights.")
    }

    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths

    if( Depths_spec[[ 'def' ]] == 'MBD' )
    {
      Depths = multiMBD( listOfValues,
                         weights = Depths_spec[[ 'weights' ]],
                         manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( listOfValues )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == Data$N )
  }

  out = .fbplot_mfData( time_grid, listOfValues, Depths, Fvalue  )

  ID_out = out$ID_out

  if( is.numeric( display ) )
  {
    dev.set( display )
  }

  if( ! display == FALSE )
  {
    if( ! is.null( ylab ) )
    {
      if( length( ylab ) == 1 )
      {
        ylab = rep( ylab, Data$L )
      } else if( length( ylab ) != Data$L )
      {
        stop( 'Error in fbplot.mfData: you specified a wrong number of y
              labels' )
      }
      } else {
        ylab = rep( list( '' ), Data$L )
    }

    if( ! is.null( main ) )
    {
      if( length( main ) == 1 )
      {
        main = rep( main, Data$L )
      } else if( length( main ) != Data$L )
      {
        stop( 'Error in fbplot.mfData: you specified a wrong number of
              subtitles' )
      }
      } else {
        main = rep( list( '' ), Data$L )
    }

    # Subdividing the graphical window
    mfrow_rows = ceiling( Data$L / 2 )
    mfrow_cols = 2

    par( mfrow = c( mfrow_rows, mfrow_cols ) )

    # Creating color palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( Data$N - length( ID_out ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                    c = 150 )( length( ID_out ) )
    col_envelope = set_alpha( 'blue', alpha = 0.4 )
    col_center = set_alpha( 'blue', alpha = 1 )
    col_fence_structure = 'darkblue'

    time_grid = seq( Data$t0, Data$tP, length.out = Data$P )

    xlab = ifelse( is.null( xlab ), '', xlab )

    for( iL in 1 : Data$L )
    {
      Data_curr = Data$fDList[[ iL ]]$values

      if( length( ID_out ) > 0 )
      {
        # Plotting non-outlying data
        matplot( time_grid,
                 t( Data_curr[ - ID_out, ] ), lty = 1, type = 'l',
                 col = col_non_outlying,
                 ylim = range( Data_curr ),
                 xlab = xlab,
                 ylab = ylab[[ iL ]],
                 main = main[[ iL ]], ... )

        # Computing maximum and minimum envelope
        max_envelope_limit = apply( Data_curr[ - ID_out, ], 2, max )
        min_envelope_limit = apply( Data_curr[ - ID_out, ], 2, min )
      } else {
        # Plotting all data
        matplot( time_grid,
                 t( Data_curr ), lty = 1, type = 'l',
                 col = col_non_outlying,
                 ylim = range( Data_curr ),
                 xlab = xlab, ylab = ylab, main = main, ... )

        # Computing maximum and minimum envelope
        max_envelope_limit = apply( Data_curr, 2, max )
        min_envelope_limit = apply( Data_curr, 2, min )
      }

      # Filling in the central envelope

      polygon( c(time_grid, rev( time_grid) ),
               c( as.numeric( out$min_envelope_central[ iL, ] ),
                  rev( as.numeric( out$max_envelope_central[ iL, ] ) ) ),
               col = col_envelope, border = NA)
      lines( time_grid, as.numeric( out$max_envelope_central[ iL, ] ),
             lty = 1, col = col_envelope, lwd = 3 )
      lines( time_grid, as.numeric( out$min_envelope_central[ iL, ] ),
             lty = 1, col = col_envelope, lwd = 3 )

      # Plotting the sample median
      lines( time_grid, Data_curr[ which.max( Depths ), ], lty = 1, type = 'l',
             col = col_center, lwd = 3)

      lines( time_grid, max_envelope_limit, lty = 1,
             col = col_fence_structure, lwd = 3 )
      lines( time_grid, min_envelope_limit, lty = 1,
             col = col_fence_structure, lwd = 3 )

      # Plotting vertical whiskers
      half.time_grid = which.min( abs( time_grid - 0.5 ) )
      lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
             c( out$max_envelope_central[ iL, half.time_grid ],
                max_envelope_limit[ half.time_grid ] ),
             lty = 1, col = col_fence_structure, lwd = 3 )

      lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
             c( out$min_envelope_central[ iL, half.time_grid ],
                min_envelope_limit[ half.time_grid ] ),
             lty = 1, col = col_fence_structure, lwd = 3 )


      # Plotting outlying data
      if( length( ID_out ) > 0 )
      {
        matplot( time_grid, t( toRowMatrixForm( Data_curr[ ID_out, ] ) ),
                 lty = 1, type = 'l', col = col_outlying, lwd = 3, add = T )
      }
    }
  }

  return( list( Depth = Depths,
                Fvalue = Fvalue,
                ID_outliers = ID_out ) )
}


.fbplot_mfData = function( time_grid,
                           listOfValues,
                           Depths = list( def = 'MBD',
                                          weights = 'uniform' ),
                           Fvalue = 1.5 )
{

  L = length( listOfValues )
  N = nrow( listOfValues[[ 1 ]] )
  P = ncol( listOfValues[[ 2 ]] )

  # Checking if depths have already been provided or must be computed
  if( is.list( Depths ) )
  {
    if( length( Depths ) != 2 & ! identical( names( Depths ),
                                             c( 'def', 'weights' ) ) )
    {
      stop( " Error in .fbplot_mfData: you have to provide both a specification
            for the depth definition and one for the set of weights.")
    }

    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths

    if( Depths_spec[[ 'def' ]] == 'MBD' )
    {
      Depths = multiMBD( listOfValues,
                         weights = Depths_spec[[ 'weights' ]],
                         manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( listOfValues )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == N )
  }

  # Nice hack to extract selected row from each matrix of the list and
  # concatenate the result into a row-major matrix
  Data_center = t( sapply( listOfValues, `[`,  which.max( Depths ), 1 : P ) )

  id_central_region = which( Depths >= quantile( Depths, prob = 0.5 ) )

  max_envelope_central = t( sapply( 1 : L, function( i ) (
    apply( listOfValues[[ i ]][ id_central_region, ], 2, max ) ) ) )

  min_envelope_central = t( sapply( 1 : L, function( i ) (
    apply( listOfValues[[ i ]][ id_central_region, ], 2, min ) ) ) )

  fence_upper = ( max_envelope_central - Data_center ) * Fvalue + Data_center
  fence_lower = ( min_envelope_central - Data_center ) * Fvalue + Data_center

  ID_outlying = unique( unlist( sapply( 1 : L, function( iL ) ( which(
    apply( listOfValues[[ iL ]], 1,
           function( x ) ( any( x > as.numeric( fence_upper[ iL, ] ) ) |
                             any( x < as.numeric( fence_lower[ iL, ] ) ) ) ) )
  ) ) ) )

  return( list( ID_out = ID_outlying,
                min_envelope_central = min_envelope_central,
                max_envelope_central = max_envelope_central,
                fence_lower = fence_lower,
                fence_upper = fence_upper ) )
}
