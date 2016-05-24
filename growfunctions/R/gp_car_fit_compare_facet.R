#' Side-by-side plot panels that compare latent function 
#' values to data for different estimation models
#'
#' Uses as input the output object from the gpdpgrow() and 
#' gmrfdpgrow() functions.
#'
#' @param objects A list input where each element is a returned object
#'             from estimation with either of \code{gpdpgrow} or 
#'             \code{gmrfdpgrow} or an object that contains true 
#'             \code{N x T} matrix of true latent
#'             function values, \code{f}.  This latter input 
#'             is only needed if want to compare 
#'             estimated to true latent function values.
#' @param H An \code{N x 1} with entries in \code{1,...,M} of cluster 
#'        assignments for the \code{N}
#'        units of \code{y} under a known clustering.
#' @param label.object A character vector of length equal to \code{objects} 
#'        that contains labels for each element of \code{objects} 
#'        to be used in rendering comparison plots.
#'        Defaults to \code{label.object = c("gp_rq","gmrf_rw2")}.
#' @param units_name A character input that provides a label 
#'        for the set of \code{N} observation units.
#'        Defaults to \code{units_name = "Observation_Unit"}.
#' @param units_label A vector of labels to apply to the observation units 
#'        with length equal to the 
#'        number of unique units.  Defaults to sequential 
#'        numeric values as input with data, \code{y}.
#' @param date_field A vector of \code{Date} values for labeling the x-axis tick marks.
#'        Defaults to \code{1:T}  .
#' @param x.axis.label Text label for x-axis. Defaults to \code{"time"}.
#' @param y.axis.label Text label for y-axis. Defaults to \code{"function values"}.              
#' @return A list object containing the plot of estimated functions, faceted by cluster,
#'		and the associated \code{data.frame} object.
#'     \item{p.t}{A \code{ggplot2} plot object}
#'     \item{map}{A \code{data.frame} object that contains clustering structure of
#'        observation units.}
#' @seealso \code{\link{gpdpgrow}}, \code{\link{gmrfdpgrow}}
#' @examples 
#' {
#' library(growfunctions)
#' 
#' ## load the monthly employment count data 
#' ## for a collection of 
#' ## U.S. states from the Current 
#' ## Population Survey (cps)
#' data(cps)
#' ## subselect the columns of N x T, y, 
#' ## associated with 
#' ## the years 2009 - 2013
#' ## to examine the state level 
#' ## employment levels 
#' ## during the "great recession"
#' y_short <- cps$y[,(cps$yr_label %in% 
#'                  c(2010:2013))]
#'
#' ## run DP mixture of GP's to 
#' ## estimate posterior distributions 
#' ## for model parameters
#' ## uses default setting of a 
#' ## single "rational quadratic" 
#' ## covariance formula
#' res_gp         <- gpdpgrow(
#'                      y = y_short, 
#'                      n.iter = 3, 
#'                      n.burn = 1, 
#'                      n.thin = 1, 
#'                      n.tune = 0)  
#' ## 2 plots of estimated functions: 
#' ## 1. faceted by cluster and fit;
#' ## 2.  data for experimental units.
#' ## for a group of randomly-selected 
#' ## functions
#' fit_plots_gp   <- cluster_plot( 
#'  object = res_gp,  units_name = "state", 
#'  units_label = cps$st, single_unit = FALSE, 
#'  credible = TRUE )
#'                                    
#' ## Run the DP mixture of iGMRF's to 
#' ## estimate posterior 
#' ## distributions for model parameters
#' ## Under default 
#' ## RW2(kappa) = order 2 trend 
#' ## precision term
#' res_gmrf     <- gmrfdpgrow(y = y_short, 
#'                        n.iter = 13, 
#'                        n.burn = 4, 
#'                        n.thin = 1) 
#'                                      
#' ## 2 plots of estimated functions: 
#' ## 1. faceted by cluster and fit;
#' ## 2.  data for experimental units.
#' ## for a group of randomly-selected functions
#' fit_plots_gmrf   <- cluster_plot( object = res_gmrf, 
#'   units_name = "state", units_label = cps$st, 
#'   single_unit = FALSE, 
#'   credible = TRUE )                                    
#'                                      
#' ## visual comparison of fit performance 
#' ## between gpdpgrow() and gmrfdpgrow()
#' ## or any two objects returned from any
#' ## combination of these estimation
#' ## functions
#' objects        <- vector("list",2)
#' objects[[1]]   <- res_gmrf
#' objects[[2]]   <- res_gp
#' label.object   <- c("gmrf_tr2","gp_rq")
#' ## the map data.frame object 
#' ## from fit_plots gp 
#' ## includes a field that 
#' ## identifies cluster assignments
#' ## for each unit (or domain)
#' H        <- fit_plots_gp$map$cluster
#' fit_plot_compare_facet <- 
#' fit_compare( objects = objects, 
#'  H = H, label.object = label.object,
#'  y.axis.label = "normalized y",
#'  units_name = "state", units_label = cps$st)                                  
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases fit_compare
#' @export
fit_compare <- function(objects, H = NULL, label.object = c("gp_rq","gmrf_rw2"), 
                        units_name = "Observation_Unit", 
                        units_label = NULL, 
                        date_field = NULL, x.axis.label = NULL, y.axis.label = NULL)
{
     ## read in the data
     y                   <- objects[[1]]$optpartition$y ## N x T data matrix used for modeling 
     ## for gmrfdpgrow object, y may have intermittant missing values (= NA).
     ## there is also returned an N x T, y_bar, that estimates the missing values, not used here.
     N                   <- nrow(y) ## number of experimental units
     T                   <- ncol(y) ## numer of time-indexed observations for each unit
     L                   <- length(objects) ## number of models to compare fit
     
     ###############################################################################
     ## compose plot of fitted functions vs. actual data for randomly-selected units
     ###############################################################################
     ## create link of units_label to replace numbers with labels
     units_numeric                          <- 1:N ##sort(unique(map$units_numeric))
     if( length(units_label) > 0 )
     {
          ## units_label                       <- sort(unique(units_label)) ## a vector of entries
          if( length(units_label) != nrow(y) )
          {
               stop("\nLength of units_label must be equal to nrow(y) since labels apply to
                    observed units.\n")
          }
     }else{ ## just set units_label = to numeric values from bigSmin is units_label == NULL
          units_label                       <- units_numeric
     }
     
     d3             <- data.frame(units_label,y) ## y is N x T
     names(d3)      <- c(eval(units_name),1:T) ## even though month is date/numeric, names will be 
                                               ## char
     dat_gca        <- melt(d3,measure.vars = as.character(1:T), variable.name = "time",
                            value.name = "fit")
     
     if( !is.null(date_field))
     {
          dat_date            <- data.frame(1:T,sort(unique(date_field)))
          names(dat_date)     <- c("time","date")
          dat_gca             <- merge(dat_gca,dat_date, by = "time", all.x=TRUE)
          dat_gca$time        <- as.Date(dat_gca$date)
          dat_gca$date        <- NULL
     }
        
     ##
     ## compose fit predictions
     ##
     
     ## bases
     ## create plot data.frame
     dat_gcp       <- lapply(1:L,function(l){
                              if(label.object[l] == "true")  ## can choose not to input truth
                              {    ## 1 x (N*T) function samples. N is fast-moving
                                   fit_l          <- objects[[l]]$bb ## already N x T
                              }else{
                                   hat_l          <- as.vector(colMeans(objects[[l]]$bb))
                                   ## assuming T is the slow-moving index
                                   fit_l          <- matrix(hat_l,N,T,byrow=FALSE)
                              }
                              
                              ## build data.frame
                              x              <- data.frame(d3[,eval(units_name)],
                                                           eval(label.object[l]),fit_l)
                              names(x)       <- c(eval(units_name),"type",1:T)
                              ## render predicted data in long form
                              x              <- melt(x,measure.vars = as.character(1:T),
                                                     variable.name="time",
                                                     value.name = "fit")
                              x$est_type    <- interaction(x$type,x[,eval(units_name)])
                              if( !is.null(date_field))
                              {
                                   x             <- merge(x,dat_date, by = "time", all.x=TRUE)
                                   x$time        <- as.Date(x$date)
                                   x$date        <- NULL
                              }
                         x})
     ## concatenate list of data.frames indexed by model type
     dat_gcp             <- do.call("rbind",dat_gcp)
     
     
     #####
     ## randomly select records to plot
     #####
     
     ## build data.frame that maps experimental units to clusters
     ## then subset the data.frame by cluster membership and randomly select a unit
     ## from each cluster for plotting
     if( is.null(H) & (length(grep("true",label.object) >= 1)) )
     {
          H              <- objects[[which(label.object == "true")]]$H 
     }else{ ## did not input H
          if( is.null(H) )
          {
               stop("need to input cluster assignment vector, H, or dat_sim object with true values")
          }     
     } ## end conditional statement on whether input H
     
     map            <- data.frame(1:length(H),H)
     names(map)       <- c("units_numeric","cluster")
     
     ## add in units_label to map in order to merge with plot sets to convey cluster labels
     stopifnot(length(units_label) == length(units_numeric))
     tmp                                <- data.frame(units_label,units_numeric)
     names(tmp)                         <- c(eval(units_name),"units_numeric")
     map                                <- merge(map,tmp,by="units_numeric",all.x=TRUE)
     map                                <- map[order(map$cluster,map[,eval(units_name)]),]
     ## merge in cluster labels
     dat_gcp                            <- merge(dat_gcp,map,all.x=TRUE,by=eval(units_name),
                                                 sort = FALSE)
     dat_gca                            <- merge(dat_gca,map,all.x=TRUE,by=eval(units_name),
                                                 sort = FALSE)
     
     ## select random unit to plot from within each cluster
     plot_st                  <- vector("list",length(map$cluster))
     for(m in 1:length(unique(map$cluster)))
     {
          map_clu             <- subset( map, cluster == m )
          plot_st[[m]]        <- sample(as.character(map_clu[,eval(units_name)]),1,replace = FALSE)
     }
     plot_st                  <- unlist(plot_st)
     plot_gca                 <- subset(dat_gca, dat_gca[,eval(units_name)] %in% plot_st)
     plot_gcp                 <- subset(dat_gcp, dat_gcp[,eval(units_name)] %in% plot_st)
     
     
     ##
     ## render growth curve plot
     ##
     p.t          		= ggplot(data=plot_gcp,aes(x=time, y = fit) )
     l				= geom_line(aes(group = est_type),
                              size = 1.2, alpha = 0.6)
     l.2     			= geom_point(data=plot_gca,aes_string(group = eval(units_name)),size=3,shape=1,
                             colour="black") ## Add the data
     l.3          		= geom_line(data=plot_gca,aes_string(group = eval(units_name)),lty="dotted")  
     ## Add the data  
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

     ## f				= facet_wrap(as.formula(paste(" ~", cluster)), scales="fixed")
     f                   = facet_grid(cluster~type,scales="free")
     p.t     		     <- p.t + l + f + axis + l.2 + l.3 + 
                              theme_bw() + theme(axis.text.x = element_blank(), 
                                                 axis.ticks.x = element_blank())  
     suppressWarnings(print(p.t))
     
     
     fit <- cluster <- est_type <- NULL
     return(invisible(list(p.t = p.t, map = map)))
     
} ## end function gp_car_fit_compare