
#' Outliergram for univariate functional datasets
#'
#' This function performs the outliergram of a univariate functional dataset,
#' possibly with an adjustment of the true positive rate of outliers discovered
#' under assumption of gaussianity.
#'
#' @section Adjustment:
#'
#' When the adjustment option is selected, the value of \eqn{F} is optimised for
#' the univariate functional dataset provided with \code{fData}. In practice,
#' a number \code{adjust$N_trials} of times a synthetic population
#' (of size \code{adjust$trial_size} with the same covariance (robustly
#' estimated from data) and centerline as \code{fData} is simulated without
#' outliers and each time an optimised value \eqn{F_i} is computed so that a
#' given proportion (\code{adjust$TPR}) of observations is flagged as outliers.
#' The final value of \code{F} for the outliergram is determined as an average
#' of \eqn{F_1, F_2, \ldots, F_{N_{trials}}}. At each time step the optimisation
#' problem is solved using \code{stats::uniroot} (Brent's method).
#'
#' @param fData the univariate functional dataset whose outliergram has to be
#' determined.
#' @param MBD_data a vector containing the MBD for each element of the dataset.
#' If missing, MBDs are computed.
#' @param MEI_data a vector containing the MEI for each element of the dataset.
#' If not not provided, MEIs are computed.
#' @param q_low parameter used in the part where data are shifted toward the
#' center of the dataset. It indicates the quantile to be used to compute the
#' target to compare functions in the secondary check for outliers. Defult is 0,
#' i.e. High MEI functions (lying at the bottom of the dataset) are compared
#' to the minimum of all the remaining functions.
#' @param q_high parameter used in the part where data are shifted toward the
#' center of the dataset. It indicates the quantile to be used to compute the
#' target to compare functions in the secondary check for outliers.
#' Defult is 1, i.e. Low MEI functions (lying at the top of the dataset) are
#' compared to the maximum of all the remaining functions.
#' @param p_check percentage of observations with either low or high MEI to be
#' checked for outliers in the secondary step (shift towards the center of the
#' dataset).
#' @param Fvalue the \eqn{F} value to be used in the procedure that finds the
#' shape outliers by looking at the lower parabolic limit in the outliergram.
#' Default is \code{1.5}. You can also leave the default value and, by providing
#' the parameter \code{adjust}, specify that you want \code{Fvalue} to be
#' adjusted for the dataset provided in \code{fData}.
#' @param adjust either \code{FALSE} if you would like the default value for the
#' inflation factor, \eqn{F = 1.5}, to be used, or a list specifying the
#' parameters required by the adjustment.
#'  \itemize{
#'  \item{"\code{N_trials}"}{: the number of repetitions of the adujustment
#'  procedure based on the simulation of a gaussisan population of functional
#'  data, each one producing an adjusted value of \eqn{F}, which will lead
#'  to the averaged adjusted value \eqn{\bar{F}}. Default is 20;}
#'  \item{"\code{trial_size}"}{: the number of elements in the gaussian
#'  population of functional data that will be simulated at each repetition of
#'  the adjustment procedure. Default is \code{5 * fData$N};}
#'  \item{"\code{TPR}"}{: the True Positive Rate of outleirs, i.e. the proportion
#'  of observations in a dataset without shape outliers that have to be considered
#'  outliers. Default is \code{2 * pnorm( 4 * qnorm( 0.25 ) )};}
#'  \item{"\code{F_min}"}{: the minimum value of \eqn{F}, defining the left
#'  boundary for the optimisation problem aimed at finding, for a given dataset
#'  of simulated gaussian data associated to \code{fData}, the optimal value of
#'  \eqn{F}. Default is 0.5;}
#'  \item{"\code{F_max}"}{: the maximum value of \eqn{F}, defining the right
#'  boundary for the optimisation problem aimed at finding, for a given dataset
#'  of simulated gaussian data associated to \code{fData}, the optimal value of
#'  \eqn{F}. Default is 20;}
#'  \item{"\code{tol}"}{: the tolerance to be used in the optimisation problem
#'  aimed at finding, for a given dataset of simulated gaussian data associated
#'  to \code{fData}, the optimal value of \eqn{F}. Default is \code{1e-3};}
#'  \item{"\code{maxiter}"}{: the maximum number of iterations to solve the
#'  optimisation problem aimed at finding, for a given dataset of simulated
#'  gaussian data associated to \code{fData}, the optimal value of \eqn{F}.
#'  Default is \code{100};}
#'  \item{"\code{VERBOSE}"}{: a parameter controlling the verbosity of the
#'  adjustment process;}
#'  }
#' @param display either a logical value indicating wether you want the
#' outliergram to be displayed, or the number of the graphical device
#' where you want the outliergram to be displayed.
#' @param xlab a list of two labels to use on the x axis when displaying the
#' functional dataset and the outliergram
#' @param ylab a list of two labels to use on the y axis when displaying the
#' functional dataset and the outliergram;
#' @param main a list of two titles to be used on the plot of the functional
#' dataset and the outliergram;
#' @param ... additional graphical parameters to be used \emph{only} in the plot
#' of the functional dataset
#'
#' @return
#'
#' Even when used graphically to plot the outliergram, the function returns a
#' list containing a numeric vector with the IDs of observations in
#' \code{fData} that are considered as shape outliers and the value of
#' \code{Fvalue} that has been used in determining them.
#'
#' @references
#'
#' Arribas-Gil, A., and Romo, J. (2014). Shape outlier detection and visualization
#' for functional data: the outliergram, \emph{Biostatistics}, 15(4), 603-619.
#'
#' @seealso \code{\link{fData}}, \code{\link{MEI}}, \code{\link{MBD}},
#' \code{\link{fbplot}}
#'
#' @examples
#'
#'
#' set.seed( 1618 )
#'
#' N = 200
#' P = 200
#' N_extra = 4
#'
#' grid = seq( 0, 1, length.out = P )
#'
#' Cov = exp_cov_function( grid, alpha = 0.2, beta = 0.8 )
#'
#' Data = generate_gauss_fdata( N,
#'                              centerline = sin( 4 * pi * grid ),
#'                              Cov = Cov )
#'
#' Data_extra = array( 0, dim = c( N_extra, P ) )
#'
#' Data_extra[ 1, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid + pi / 2 ),
#'                                          Cov = Cov )
#'
#' Data_extra[ 2, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid - pi / 2 ),
#'                                          Cov = Cov )
#'
#' Data_extra[ 3, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid + pi/ 3 ),
#'                                          Cov = Cov )
#'
#' Data_extra[ 4, ] = generate_gauss_fdata( 1,
#'                                          sin( 4 * pi * grid - pi / 3),
#'                                          Cov = Cov )
#' Data = rbind( Data, Data_extra )
#'
#' fD = fData( grid, Data )
#'
#' outliergram( fD, display = TRUE )
#'
#' outliergram( fD, Fvalue = 10, display = TRUE )
#' \dontrun{
#' outliergram( fD,
#'              adjust = list( N_trials = 10,
#'                             trial_size = 5 * nrow( Data ),
#'                             TPR = 0.01,
#'                             VERBOSE = FALSE ),
#'              display = TRUE )
#' }
#'
#' @importFrom grDevices dev.set dev.cur
#' @importFrom stats cor pnorm rnorm qnorm uniroot
#' @importFrom graphics text lines polygon plot points matplot par
#'
#' @export
outliergram = function( fData, MBD_data = NULL, MEI_data = NULL,
                        q_low = 0, q_high = 1, p_check = 0.05,
                        Fvalue = 1.5,
                        adjust = FALSE, display = TRUE,
                        xlab = NULL, ylab = NULL, main = NULL, ... )
{
  N = fData$N

  grid = seq( fData$t0,
              fData$tP,
              length.out = fData$P )

  if( ! is.list( adjust ) )
  {
    # Plain outliergram with default F value: F = 1.5

    if( Fvalue == 1.5 )
    {
      out = .outliergram( fData,
                          MBD_data = MBD_data, MEI_data = MEI_data,
                          p_check = p_check, q_low = q_low, q_high = q_high,
                          shift = TRUE )

    } else {
      out = .outliergram( fData,
                          MBD_data = MBD_data, MEI_data = MEI_data,
                          p_check = p_check, q_low = q_low, q_high = q_high,
                          Fvalue = Fvalue,
                          shift = TRUE )
    }
  } else {

    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )

    trial_size = ifelse( is.null( adjust$trial_size ),
                         5 * fData$N,
                         adjust$trial_size )

    TPR = ifelse( is.null( adjust$TPR ),
                  2 * pnorm( 4 * qnorm( 0.25 ) ),
                  adjust$TPR )

    F_min = ifelse( is.null( adjust$F_min ),
                    0.5,
                    adjust$F_min )

    F_max= ifelse( is.null( adjust$F_max ),
                    20,
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

    Cov = robustbase::covOGK( fData$values, sigmamu = robustbase::s_Qn )$cov

    CholCov <- chol( Cov )

    if( is.null ( MBD_data ) )
    {
      MBD_data = MBD( fData$values, manage_ties = TRUE )
    }

    centerline = fData$values[ which.max( MBD_data ), ]

    Fvalues = rep( 0, N_trials )

    obj_function = function( F_curr )( length(
      .outliergram( fData_gauss,
                    MBD_data = NULL,
                    MEI_data = NULL,
                    q_low = q_low,
                    q_high = q_high,
                    Fvalue = F_curr,
                    shift = FALSE )$ID_SO ) / trial_size - TPR )

    for( iTrial in 1 : N_trials )
    {
      if( VERBOSE > 0 )
      {
        cat( ' * * * Iteration ', iTrial, ' / ', N_trials, '\n' )
      }

      fData_gauss = fData( grid,
                           generate_gauss_fdata( N = trial_size,
                                                 centerline = centerline,
                                                 CholCov = CholCov ) )

      if( VERBOSE > 0 )
      {
        cat( ' * * * * beginning optimisation\n' )
      }

      opt = uniroot( obj_function,
                     interval = c( F_min, F_max ),
                     tol = tol,
                     maxiter = maxiter )

      Fvalues[ iTrial ] = opt$root
    }

    Fvalue = mean( Fvalues )

    out = .outliergram( fData, MBD_data, MEI_data, p_check, q_low,
                        q_high, Fvalue = Fvalue, shift = TRUE  )
  }

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Plot
  if( display )
  {

    if( is.null( xlab ) )
    {
      xlab = list( '', 'MEI' )
    }

    if( is.null( ylab ) )
    {
      ylab = list( '', 'MBD' )
    }

    if( is.null( main ) )
    {
      main = list( '', 'Outliergram' )
    }

    # Setting up palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( length( out$ID_NO ) )

    col_non_outlying = set_alpha( col_non_outlying, 0.5 )

    col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                    c = 150 )( length( out$ID_SO ) )

    dev.cur()
    par( mfrow = c( 1, 2 ) )

    # Plotting functional data
    if( length( out$ID_SO ) > 0 )
    {
      matplot( grid, t( fData$values[ - out$ID_SO, ] ), type = 'l', lty = 1,
               ylim = range( fData$values ),
               col = col_non_outlying,
               xlab = xlab[[1]],
               ylab = ylab[[1]],
               main = main[[1]],
               ... )
      matplot( grid, t( toRowMatrixForm( fData$values[ out$ID_SO, ] ) ),
               type = 'l', lty = 1, lwd = 3, ylim = range( fData$values ),
               col = col_outlying, add = TRUE )
    } else {
      matplot( grid, t( fData$values ), type = 'l', lty = 1,
               ylim = range( fData$values ),
               col = col_non_outlying,
               xlab = xlab[[1]],
               ylab = ylab[[1]],
               main = main[[1]],
               ... )
    }


    # Adding text labels with curve ID
    w_spacing = diff( range( grid ) ) / ( 2 * length( out$ID_SO ) )

    for( iOut in seq_along( out$ID_SO ) )
    {
      text( grid[ 1 ] + ( 2 * iOut - 1 ) * w_spacing,
            fData$values[ out$ID_SO[ iOut ],
                  which.min( abs( grid - grid[ 1 ] -
                                    ( 2 * iOut - 1 ) * w_spacing ) ) ] +
              diff( range( fData$values[ out$ID_SO[ iOut ]  ] ) ) / 30,
            out$ID_SO[ iOut ],
            col = col_outlying[ iOut ] )
    }

    # Plotting outliergram

    # Upper parabolic limit
    grid_1D = seq( 0, 1, length.out = 100 )

    plot( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2,
          lty = 2, type = 'l', col = 'darkblue', lwd = 2,
          ylim = c( 0, a_0_2 + a_1 / 2 + a_0_2 * N^2/4 ),
          xlab = xlab[[2]],
          ylab = ylab[[2]],
          main = main[[2]] )

    if( length( out$ID_SO ) > 0 )
    {
      points( out$MEI_data[ - out$ID_SO ], out$MBD_data[ - out$ID_SO ],
              pch = 16, col = col_non_outlying )
      points( out$MEI_data[ out$ID_SO ], out$MBD_data[ out$ID_SO ],
              pch = 16, cex = 1.5, col = col_outlying )
      for( idOut in out$ID_SO )
      {
        text( out$MEI_data[ idOut ],
              out$MBD_data[ idOut ] + 0.5 / 30,
              idOut,
              col = col_outlying[ match( idOut, out$ID_SO ) ] )
      }
    } else {
      points( out$MEI_data, out$MBD_data,
              pch = 16, col = col_non_outlying )
    }

    # lower parabolic limit
    if( Fvalue == 1.5 )
    {
      lines( grid_1D, a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 -
               out$Q_d3 - 1.5 * out$IQR_d,
             lty = 2, lwd = 2, col = 'lightblue' )
    }
    else
    {
      lines( grid_1D,  a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 -
               Fvalue * out$Q_d1,
             lty = 2, lwd = 2, col = 'lightblue' )
    }

  }

  return( list( Fvalue = Fvalue,
                ID_outliers = out$ID_SO ) )
}

.outliergram = function( fData, MBD_data = NULL, MEI_data = NULL,
                         p_check = 0.05, q_low = 0, q_high = 1,
                         Fvalue = NULL, shift = TRUE )
{
  N = fData$N

  a_0_2 = -2 / ( N * ( N - 1 ) )

  a_1 = 2 * ( N + 1 ) / ( N - 1 )

  # Computing MBD
  if( is.null( MBD_data ) ){

    MBD_data = MBD( fData$values )
  }

  # Computing MEI
  if( is.null( MEI_data ) )
  {
    MEI_data = MEI( fData$values )
  }

  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data

  Q = quantile( d )

  Q_d3 = Q[ 4 ]

  Q_d1 = Q[ 2 ]

  IQR_d = Q[ 4 ] - Q[ 2 ]

  # Computing surely outlying curves

  if( is.null( Fvalue ) )
  {
    ID_shape_outlier = which( d >= Q_d3 + 1.5 * IQR_d )
  } else {

    ID_shape_outlier = which( d >= Fvalue * Q_d1 )
  }

  # Computing non outlying curves ids
  ID_non_outlying = setdiff( 1 : nrow( fData$values ), ID_shape_outlier )

  if( shift )
  {
    # # Low MEI curves will be checked for upward shift
    # ID_non_outlying_Low_MEI = ID_non_outlying[
    #   which( MEI_data[ - ID_shape_outlier ] >=
    #            quantile( MEI_data,
    #                      probs = 1 - p_check ) ) ]
    #
    # # High MEI curves will be checked for downward shift
    # ID_non_outlying_High_MEI = ID_non_outlying[
    #   which( MEI_data[ - ID_shape_outlier ] <=
    #            quantile( MEI_data, probs = p_check ) ) ]

    # Low MEI curves will be checked for upward shift
    ID_non_outlying_Low_MEI = ID_non_outlying[
      which( MEI_data[ - ID_shape_outlier ] <=
               quantile( MEI_data,
                         probs = p_check ) ) ]

    # High MEI curves will be checked for downward shift
    ID_non_outlying_High_MEI = ID_non_outlying[
      which( MEI_data[ - ID_shape_outlier ] >=
               quantile( MEI_data, probs = 1 - p_check ) ) ]


    aux_function = function( ID )( min( fData$values[ ID, ] -
                                          apply( fData$values[ - ID, ], 2,
                                                 quantile,
                                                 probs = q_low ) ) )

    # Managing High MEI data
    min_diff_min = sapply( ID_non_outlying_High_MEI, aux_function )

    ID_to_check = ID_non_outlying_High_MEI[ min_diff_min < 0 ]

    aux_function_MBD = function( ID )(
      MBD( rbind( fData$values[ - ID, ],
                  Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )

    aux_function_MEI = function( ID )(
      MEI( rbind( fData$values[ - ID, ],
                  Data_tilde[ grep( ID, ID_to_check ), ] ) )[ N ] )

    if( length( ID_to_check ) > 0 )
    {
      Data_tilde = toRowMatrixForm( fData$values[ ID_to_check, ] -
                                      min_diff_min[ min_diff_min < 0 ] )

      MBD_curr = sapply( ID_to_check, aux_function_MBD )

      MEI_curr = sapply( ID_to_check, aux_function_MEI )

      d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

      if( is.null( Fvalue ) )
      {
        ID_out_extra = ID_to_check[ which( d_curr >= Q_d3 + 1.5 * IQR_d ) ]

      } else {
        ID_out_extra = ID_to_check[ which( d_curr >= Q_d1 * Fvalue ) ]
      }

      ID_shape_outlier = c( ID_shape_outlier, ID_out_extra )
      ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra )
    }

    # Managing Low MEI data

    aux_function = function( ID )( max( fData$values[ ID, ] -
                                          apply( fData$values[ - ID, ],
                                                 2,
                                                 quantile,
                                                 probs = q_high ) ) )

    max_diff_max = sapply( ID_non_outlying_Low_MEI, aux_function )

    ID_to_check = ID_non_outlying_Low_MEI[ max_diff_max > 0 ]

    if( length( ID_to_check ) > 0 )
    {
      Data_tilde = toRowMatrixForm( fData$values[ ID_to_check, ] -
                                      max_diff_max[ max_diff_max > 0 ] )

      MBD_curr = sapply( ID_to_check, aux_function_MBD )

      MEI_curr = sapply( ID_to_check, aux_function_MEI )

      d_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2 - MBD_curr

      if( is.null( Fvalue ) )
      {
        ID_out_extra = ID_to_check[ which( d_curr >= Q_d3 + 1.5 * IQR_d ) ]
      } else {
        ID_out_extra = ID_to_check[ which( d_curr >= Q_d1 * Fvalue ) ]
      }


      ID_shape_outlier = c( ID_shape_outlier, ID_out_extra )
      ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra )
    }
  }

  return( list( ID_SO = ID_shape_outlier,
                ID_NO = ID_non_outlying,
                MEI_data = MEI_data,
                MBD_data = MBD_data,
                Q_d3 = Q_d3,
                Q_d1 = Q_d1,
                IQR_d = IQR_d ) )

}


# check_outlier_Low_MEI = function( Data, idObs, q_high = 1  )
# {
#   diff_max_vector = Data[ idObs, ] - apply( Data[ - idObs, ],
#                                             2,
#                                             ifelse( q_high == 1,
#                                                     max,
#                                                     function( x )( quantile( x, probs = q_high ) ) ) )
#
#   max_diff_max = max( diff_max_vector )
#
#   Data_tilde = Data[ idObs, ] - any ( diff_max_vector > 0 ) * max_diff_max
#
#   MBD_curr = MBD( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   MEI_curr = MEI( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   P_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2
#
#   return( as.logical( MBD_curr <= P_curr - Q_d3 - 1.5 * IQR ) )
# }
#
# check_outlier_High_MEI = function( Data, idObs, q_low = 0 )
# {
#   diff_min_vector = Data[ idObs, ] - apply( Data[ - idObs, ],
#                                             2,
#                                             ifelse( q_low == 0,
#                                                     min,
#                                                     function( x ) ( quantile( x, probs = q_low ) ) ) )
#
#   min_diff_min = min( diff_min_vector )
#
#   Data_tilde = Data[ idObs, ] - any ( diff_min_vector < 0 ) * min_diff_min
#
#   MBD_curr = MBD( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   MEI_curr = MEI( rbind( Data[ - idObs, ], Data_tilde ) )[ N ]
#
#   P_curr = a_0_2 + a_1 * MEI_curr + N^2 * a_0_2 * MEI_curr^2
#
#   return( as.logical( MBD_curr <= P_curr - Q_d3 - 1.5 * IQR ) )
# }


# par_outliergram = function( grid, Data, MBD_data = NULL, MEI_data = NULL, q_low = 0, q_high = 1, Fvalue = NULL )
# {
#   N = nrow( Data )
#
#   a_0_2 = -2 / ( N * ( N - 1 ) )
#
#   a_1 = 2 * ( N + 1 ) / ( N - 1 )
#
#   # Computing MBD
#   if( is.null( MBD_data ) ){
#
#     MBD_data = MBD( Data )
#   }
#
#   # Computing MEI
#   if( is.null( MEI_data ) )
#   {
#     MEI_data = MEI( Data )
#   }
#
#   d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data
#
#   Q = quantile( d )
#
#   Q_d3 = Q[ 4 ]
#
#   Q_d1 = Q[ 2 ]
#
#   IQR_d = Q[ 4 ] - Q[ 2 ]
#
#   # Computing surely outlying curves
#   ID_shape_outlier = ifelse( is.null( Fvalue ),
#                              which( d >= Q_d3 + 1.5 * IQR_d ),
#                              which( d >= Fvalue * Q_d1 ) )
#
#   # Computing non outlying curves ids
#   ID_non_outlying = setdiff( 1 : nrow( Data ), ID_shape_outlier )
#
#   # Low MEI curves will be checked for upward shift
#   ID_non_outlying_Low_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] < 0.5 ) ]
#
#   # High MEI curves will be checked for downward shift
#   ID_non_outlying_High_MEI = ID_non_outlying[ which( MEI_data[ ID_non_outlying ] >= 0.5 ) ]
#
#   registerDoParallel( 8 )
#   min_diff_min = foreach( ID = ID_non_outlying_High_MEI, .combine = 'c' ) %dopar% {
#     min( Data[ ID, ] -
#            apply( Data[ - ID, ],
#                   2, min ) )
#   }
#   max_diff_max = foreach( ID = ID_non_outlying_Low_MEI, .combine = 'c' ) %dopar% {
#     max( Data[ ID, ] -
#            apply( Data[ - ID, ],
#                   2, max ) )
#   }
#
#   ID_to_check_high = ID_non_outlying_High_MEI[ min_diff_min < 0 ]
#   ID_to_check_low = ID_non_outlying_Low_MEI[ max_diff_max < 0 ]
#
#   Data_tilde_high = t( t( Data[ ID_to_check_high, ] ) - min_diff_min[ ID_to_check_high ] )
#   Data_tilde_low = t( t( Data[ ID_to_check_low, ] ) - max_diff_max[ ID_to_check_low ] )
#
#   MBD_curr_high = foreach( ID = ID_to_check_high, .combine = 'c' ) %dopar% {
#     MBD( rbind( Data[ - grep( ID, ID_non_outlying_High_MEI ), ],
#                 Data_tilde_high[ grep( ID, ID_to_check_high ), ] ) )[ N ]
#   }
#   MBD_curr_low = foreach( ID = ID_to_check_low, .combine = 'c' ) %dopar% {
#     MBD( rbind( Data[ - grep( ID, ID_non_outlying_Low_MEI ), ],
#                 Data_tilde_low[ grep( ID, ID_to_check_low ), ] ) )[ N ]
#   }
#
#   MEI_curr_high = foreach( ID = ID_to_check_high, .combine = 'c'  ) %dopar% {
#     MEI( rbind( Data[ - grep( ID, ID_non_outlying_High_MEI ), ],
#                 Data_tilde_high[ grep( ID, ID_to_check_high ), ] ) )[ N ]
#   }
#   MEI_curr_low = foreach( ID = ID_to_check_low, .combine = 'c' ) %dopar% {
#     MEI( rbind( Data[ - grep( ID, ID_non_outlying_Low_MEI ), ],
#                 Data_tilde_low[ grep( ID, ID_to_check_low ), ] ) )[ N ]
#   }
#
#   d_curr_high = a_0_2 + a_1 * MEI_curr_high + N^2 * a_0_2 * MEI_curr_high^2 - MBD_curr_high
#   d_curr_low = a_0_2 + a_1 * MEI_curr_low + N^2 * a_0_2 * MEI_curr_low^2 - MBD_curr_low
#
#   ID_out_extra_high = ID_to_check_high[ which( ifelse( is.null( Fvalue ),
#                                              d_curr_high >= Q_d3 + 1.5 * IQR_d,
#                                              d_curr_high >= Q_d1 * Fvalue ) ) ]
#   ID_out_extra_low = ID_to_check_low[ which( ifelse( is.null( Fvalue ),
#                                                      d_curr_low >= Q_d3 + 1.5 * IQR_d,
#                                                      d_curr_low >= Q_d1 * Fvalue ) ) ]
#
#   ID_shape_outlier = c( ID_shape_outlier, ID_out_extra_high )
#   ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra_high )
#   ID_shape_outlier = c( ID_shape_outlier, ID_out_extra_low )
#   ID_non_outlying = setdiff( ID_non_outlying, ID_out_extra_low )
#
#
#   return( list( ID_SO = ID_shape_outlier,
#                 ID_NO = ID_non_outlying,
#                 MEI_data = MEI_data,
#                 MBD_data = MBD_data ) )
#
# }

