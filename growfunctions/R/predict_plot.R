#' Plot estimated functions both at estimated and predicted time points with 95\% credible intervals.
#'
#' Uses as input the output object from the gpdpgrow.predict() and 
#' gmrfdpgrow.predict() methods.
#'
#' @param object A \code{gpdpgrow.predict} or \code{gmrfdpgrow.predict} object, obtained from 
#'        \code{predict_functions(object,...)}. 
#' @param units_label A vector of labels to apply to experimental units with length equal to the number of
#'        unique units.  Defaults to sequential numeric values as input with data, \code{y}.
#' @param type_label A character vector assigning a "fitted" or "predicted
#'        label for the \code{time_points} input. Defaults to \code{type_label = c("fitted","predicted")}.
#' @param time_points A list input of length 2 with each entry containing a numeric vector
#'        of times -  one for the observed times for the set of "fitted" functions and the other denotes 
#'        time values at which "predicted" values were rendered for the functions.  This input variable
#'        only applies to \code{gpdpgrow} objects and not \code{gmrfdpgrow} objects because the latter 
#'        covariance structure is based on adjacency for equally-spaced time points.
#'        Defaults to \code{1:T_train} for the list entry pointed to "fitted" and 
#'        {(T_train+1):(T_train + T_test)} for the list entry pointed to "predicted".
#' @param date_label A vector of \code{Date} values for labeling the x-axis tick marks.
#'        Defaults to \code{1:T}  .
#' @param x.axis.label Text label for x-axis. Defaults to \code{"time"}.
#' @param y.axis.label Text label for y-axis. Defaults to \code{"function values"}.
#' @param single_unit A scalar boolean indicating whether to plot the fitted vs data curve for
#'        only a single experimental units (versus a random sample of 6). 
#'        Defaults to \code{FALSE}. 
#' @param credible A scalar boolean indicating whether to plot 95 percent credible intervals for
#'        estimated functions, \code{bb}, when plotting fitted functions versus data.  Defaults to
#'        \code{credible = TRUE}  .              
#' @return A list object containing the plot of estimated functions, faceted by cluster,
#'     	and the associated \code{data.frame} object.
#'     \item{p.cluster}{A \code{ggplot2} plot object}
#'     \item{dat.cluster}{A \code{data.frame} object used to generate \code{p.cluster}.}
#' @seealso \code{\link{gpdpgrow}}, \code{\link{gmrfdpgrow}}
#' @examples 
#' \dontrun{
#' library(growfunctions)
#' data(cps)
#' y_short             <- cps$y[,(cps$yr_label %in% c(2008:2013))]
#' t_train             <- ncol(y_short)
#' N                   <- nrow(y_short)
#' t_test              <- 4
#'  
#' ## Model Runs
#'
#' res_gmrf            <- gmrfdpgrow(y = y_short, 
#'                                 q_order = c(2,4), 
#'                                 q_type = c("tr","sn"), 
#'                                 n.iter = 40, 
#'                                 n.burn = 20, 
#'                                 n.thin = 1) 
#'
#' res_gp              <- gpdpgrow(y = y_short
#'                               n.iter = 10, 
#'                               n.burn = 4, 
#'                               n.thin = 1, 
#'                               n.tune = 0) 
#'
#' ## Prediction Model Runs
#' T_test             <- 4
#'
#' pred_gmrf          <- predict_functions( object = res_gmrf,
#'                                      J = 1000, 
#'                                      T_test = T_test )
#'
#' T_yshort           <- ncol(y_short)
#' pred_gp            <- predict_functions( object = res_gp, 
#'                      test_times = (T_yshort+1):(T_yshort+T_test) )
#'
#' ## plot estimated and predicted functions
#' plot_gmrf       <- predict_plot(object = pred_gmrf, 
#'                                units_label = cps$st, 
#'                                single_unit = TRUE, 
#'                                credible = FALSE)
#'
#' plot_gp         <- predict_plot(object = pred_gp, 
#'                                units_label = cps$st, 
#'                                single_unit = FALSE, 
#'                                credible = TRUE)  
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases predict_plot
#' @export
predict_plot <- function(object = NULL, units_label = NULL, 
                         type_label = c("fitted","predicted"), time_points = NULL, 
                         date_label = NULL, x.axis.label = NULL, y.axis.label = NULL, 
                         single_unit = FALSE, credible = TRUE)
{
     ## read in the data
     pred                <- object$E_bb_pred ## N x T_test matrix of mean predicted function values
     bb_pred             <- object$bbpred_draws ## a J x T_test x N array, where J denotes number of 
     ## samples drawn from the posterior predictive distribution
     bb                  <- object$bb ## nkeep x (N*T) fitted functions to training data. N is fast-moving
     N                   <- nrow(pred) ## number of experimental units
     T_test              <- ncol(pred) ## number of time-indexed test points
     T                   <- ncol(bb) / N ## number of training time points
     J                   <- dim(bb_pred)[1] ## number of MCMC draws
     L                   <- length(type_label) ## c("fitted","predicted")
     time_count          <- vector("numeric",L) ## capture time counts for training and test sets
     if( length(time_points) != L ) ## time labels not input for training and test set
     {
          time_points    <- vector("list",L) ## capture the time point labels
          for(l in 1:L)
          {
               if(type_label[l] == "fitted")
               {
                    
                    time_count[l]       <- T 
                    time_points[[l]]    <- as.vector(1:T)
               }else{
                    time_count[l]       <- T_test
                    time_points[[l]]    <- as.vector((T+1):(T+T_test))
               }
          }
     } ## end conditional statement on whether user inputs time_count
     
     ## create link of units_label to replace unit identifying numbers with labels
     units_numeric                     <- 1:N 
     if( length(units_label) > 0 )
     {
          units_label                       <- sort(unique(units_label)) ## a vector of entries
     }else{ ## just set units_label = to numeric values from bigSmin is units_label == NULL
          units_label                       <- units_numeric
     }
     
     stopifnot(length(units_label) == length(units_numeric))
     tmp                                <- data.frame(units_label,units_numeric)
     names(tmp)                         <- c("unit","units_numeric")
     
     
     ## create plot data.frame
     dat_gcp       <- lapply(1:L,function(l){
          if(type_label[l] == "fitted")
          {    
               hat_l      <- as.vector(colMeans(bb)) ## 1 x (N*T) function samples. N is fast-moving
               lo_l       <- apply(bb,2,function(x){quantile(x,probs = 0.025)}) # 1 x N*T - N is fast
               hi_l       <- apply(bb,2,function(x){quantile(x,probs = 0.975)}) # 1 x N*T - N is fast
          }else{ ## type_label[l] == "predicted"
               hat_l     <- as.vector(pred) ## takes N x T_test matrix, by column, so N is fast-moving
               ## bb_pred is a K x T_test x N array.  
               ## Transform to an N*T_test x K matrix, where N is fast-moving
               bbp_mat   <- t(matrix(aperm(bb_pred,c(3,2,1)),nrow=(N*T_test),ncol=J))
               lo_l      <- apply(bbp_mat,2,function(x){quantile(x,probs = 0.025)}) # 1 x N*T_test - N is fast
               hi_l      <- apply(bbp_mat,2,function(x){quantile(x,probs = 0.975)}) # 1 x N*T_test - N is fast       
          }
          ## create data.frame with unit, time point, fit quantiles where nrow = N*T
          x              <- data.frame(rep(tmp[,"unit"],times=time_count[l]),
                                       rep(time_points[[l]],each=N),eval(type_label[l]),hat_l,lo_l,hi_l)
          names(x)       <- c("unit","time","form","fit","lo","hi")
          x$est_form     <- interaction(x$form,x[,"unit"])
          if( length(date_label) == 2  )
          {
               dat_date_l          <- data.frame(sort(unique(time_points[l])),sort(unique(date_label[l])))
               names(dat_date)     <- c("time","date")
               x                   <- merge(x,dat_date_l, by = "time", all.x=TRUE)
               x$time              <- as.Date(x$date)
               x$date              <- NULL
          }
          x})
     ## concatenate list of data.frames indexed by model type
     dat_gcp             <- do.call("rbind",dat_gcp)
     
     ## randomly select records to plot
     if(single_unit == FALSE) ## plot the fit for a single experimental unit?
     {
          num_st_plot         <- 6
     }else{
          num_st_plot         <- 1
     }
     
     ## select from among observational units to plot
     plot_st             <- sample(tmp[,"unit"],num_st_plot,replace = FALSE)
     
     plot_gcp            <- subset(dat_gcp, dat_gcp[,"unit"] %in% plot_st)    
     ## compose plot
     p.t               	= ggplot(data=plot_gcp,aes(x=time, y = fit) )
     l				= geom_line(aes(group = "unit", colour = form), size = 1.2,
                                   alpha = 0.6)
     if( length(y.axis.label) > 0 )
     {
          if( length(x.axis.label) > 0 )
          {
               axis <- labs(x = eval(x.axis.label), y = eval(y.axis.label))
          }else{ ## no x label
               axis <- labs(x = "time", y = eval(y.axis.label) )
          }
     }else{ 
          if( length(x.axis.label) > 0 )
          {
               axis <- labs(x = eval(x.axis.label), y = "normalized y units")
          }else{ ## no labels
               axis <- labs(x = "time", y = "normalized y units" )
          }
     }
     
     f				<- facet_wrap(as.formula(paste("~", "unit")), scales="fixed")
     if( !credible )
     {
          p.t          	<- p.t + l + f + axis + theme_bw() 
     }else{ ## plot credible intervals
          ci             <- geom_ribbon(aes(ymin = lo,ymax = hi, colour = form), alpha=0.2)
          p.t            <- p.t + l + ci + f + axis + theme_bw() 
     }
     
     print(p.t)
     
     value <- fit <- lo <- hi <- form <- NULL
     return(invisible(list(p.fit = p.t, dat.fit = dat_gcp)))
     
} ## end function predict_plot    