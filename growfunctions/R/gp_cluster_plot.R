#' Plot estimated functions for experimental units faceted by cluster versus data to assess fit.
#'
#' Uses as input the output object from the gpdpgrow() and gmrfdpgrow() functions.
#'
#' @param object A \code{gpdpgrow} or \code{gmrfdpgrow} object. 
#' @param N_clusters Denotes the number of largest sized (in terms of membership) clusters to plot.
#'        Defaults to all clusters.
#' @param time_points Inputs a vector of common time points at which the collections of functions were
#'        observed (with the possibility of intermittent missingness).  The length of \code{time_points}
#'        should be equal to the number of columns in the data matrix, \code{y}.  Defaults to 
#'        \code{time_points = 1:ncol(y)}.
#' @param units_name The plot label for observation units.  Defaults to \code{units_name = "function"}.
#' @param units_label A vector of labels to apply to the observation units with length equal to the number of
#'        unique units.  Defaults to sequential numeric values as input with data, \code{y}.
#' @param date_field A vector of \code{Date} values for labeling the x-axis tick marks.
#'        Defaults to \code{1:T}  .
#' @param x.axis.label Text label for x-axis. Defaults to \code{"time"}.
#' @param y.axis.label Text label for y-axis. Defaults to \code{"function values"}.
#' @param smoother  A scalar boolean input indicating whether to co-plot a smoother line 
#'                  through the functions in each cluster.
#' @param sample_rate A numeric value in (0,1] indicating percent of functions to randomly sample within
#'                    each cluster to address over-plotting.  Defaults to 1.
#' @param single_unit A scalar boolean indicating whether to plot the fitted vs data curve for
#'        only a single experimental units (versus a random sample of 6). 
#'        Defaults to \code{single_unit = FALSE}. 
#' @param credible A scalar boolean indicating whether to plot 95 percent credible intervals for
#'        estimated functions, \code{bb}, when plotting fitted functions versus data.  Defaults to
#'        \code{credible = FALSE}
#' @param num_plot A scalar integer indicating how many randomly-selected functions to plot
#'        (each in it's own plot panel) in the plot of functions versus the observed time series
#'        in the case that \code{single_unit == TRUE}.
#'        Defaults to \code{num_plot = 6}.      
#' @return A list object containing the plot of estimated functions, faceted by cluster,
#'		and the associated \code{data.frame} object.
#'     \item{p.cluster}{A \code{ggplot2} plot object}
#'     \item{dat.cluster}{A \code{data.frame} object used to generate \code{p.cluster}.}
#' @seealso \code{\link{gpdpgrow}}, \code{\link{gmrfdpgrow}}
#' @examples 
#' {
#' library(growfunctions)
#' 
#' ## load the monthly employment count data for a collection of 
#' ## U.S. states from the Current 
#' ## Population Survey (cps)
#' data(cps)
#' ## subselect the columns of N x T, y, associated with 
#' ## the years 2008 - 2013
#' ## to examine the state level employment levels 
#' ## during the "great recession"
#' y_short             <- cps$y[,(cps$yr_label %in% c(2008:2013))]
#'
#' ## Run the DP mixture of iGMRF's to estimate posterior 
#' ## distributions for model parameters
#' ## Under default RW2(kappa) = order 2 trend 
#' ## precision term
#' res_gmrf            <- gmrfdpgrow(y = y_short, 
#'                                      n.iter = 40, 
#'                                      n.burn = 20, 
#'                                      n.thin = 1) 
#'                                      
#' ## 2 plots of estimated functions: 1. faceted by cluster and fit;
#' ## 2.  data for experimental units.
#' ## for a group of randomly-selected functions
#' fit_plots_gmrf      <- cluster_plot( object = res_gmrf, 
#'                                      units_name = "state", 
#'                                      units_label = cps$st, 
#'                                      single_unit = FALSE, 
#'                                      credible = TRUE )   
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases cluster_plot
#' @export
cluster_plot <- function(object, N_clusters = NULL, time_points = NULL, units_name = "unit", 
                         units_label = NULL, date_field = NULL, x.axis.label = NULL, 
                         y.axis.label = NULL, smoother = TRUE, 
                         sample_rate = 1.0, single_unit = FALSE, credible = FALSE,
                         num_plot = NULL)
{
     ## read in the data
     y                   <- object$optpartition$y ## N x T data matrix used for modeling 
     ## for gmrfdpgrow object, y may have intermittant missing values (= NA).
     ## there is also returned an N x T, y_bar, that estimates the missing values, not used here.
     N                   <- nrow(y) ## number of experimental units
     T                   <- ncol(y) ## numer of time-indexed observations for each unit
     ## capture time points, t_j, at which functions are observed
     if(is.null(time_points))
     {
          time_points    <- 1:T
     }
     ## build data.frame that maps experimental units to clusters
     cluster     		<- object$bigSmin
     c.sizes			<- sapply(cluster,length)
     if( (length(N_clusters) == 0) || (N_clusters > length(cluster)) )
     {    ## plot all the clusters
          clusterstoplot     	<- sort(c.sizes,decreasing = TRUE,
                                      index.return=TRUE)$ix[1:length(cluster)]
     }else{ ## plot the functions in the N_cluster largest clusters
          clusterstoplot          <- sort(c.sizes,decreasing = TRUE,
                                          index.return=TRUE)$ix[1:N_clusters]
     } ## end conditional statement on which clusters to plot
     
     map			     <- vector(mode="list",length = length(clusterstoplot))
     
     for(i in 1:length(clusterstoplot))
     {
          cluster.i			<- cluster[[clusterstoplot[i]]] 
          map[[i]]     		<- as.data.frame(cbind(cluster.i,i),stringsAsFactors = FALSE)
          names(map[[i]]) 	<- c("units_numeric","cluster")
     }
     map			                         <- do.call("rbind",map)
     ## experimental unit entries are numeric from bigSmin, with minimum value set to 1 (not 0 as in c++)
     map$units_numeric	     <- as.numeric(map$units_numeric)
     
     ## create link of units_label to replace numbers with labels
     units_numeric                     <- 1:N ##sort(unique(map$units_numeric))
     if( length(units_label) > 0 )
     {
          units_label                       <- sort(unique(units_label)) ## a vector of entries
     }else{ ## just set units_label = to numeric values from bigSmin is units_label == NULL
          units_label                       <- units_numeric
     }
     
     stopifnot(length(units_label) == length(units_numeric))
     tmp                                <- data.frame(units_label,units_numeric)
     names(tmp)                         <- c(eval(units_name),"units_numeric")
     map                                <- merge(map,tmp,by="units_numeric",all.x=TRUE)
     map                                <- map[order(map$cluster,map[,eval(units_name)]),]
     
     ###############################################################################
     ## compose plot of fitted functions grouped into membership clusters
     ###############################################################################
     
     ## create plot data.frame
     bb                            <- object$bb ## nkeep x (N*T) function samples. N is fast-moving
     bb.hat     			     <- colMeans(bb)
     units_label_dat               <- rep(eval(units_label),times=T)
     if( length(date_field) > 0 )
     {
          stopifnot(length(date_field) == T)
          month                    <- rep(sort(unique(date_field)), each = N)
     }else{ ## no date field input
          month                    <- rep(time_points, each = N)
     }
     
     dat.b				     <- data.frame(bb.hat,month,units_label_dat)
     names(dat.b)	          	<- c("value","time",eval(units_name))
     datb.clust			     <- merge(dat.b,map,all.x=TRUE,by=eval(units_name),sort = FALSE)
     
     ## sample records
     rate     		<- sample_rate
     tmp			<- split(datb.clust,list(datb.clust$cluster))
     tmp			<- unlist(sapply(tmp,function(x){
          tot_recs  	<- length(unique(x[,eval(units_name)]))
          u_recs	     <- sort(unique(x[,eval(units_name)]))
          inc_recs	     <- sample(u_recs,round(rate*tot_recs),replace = FALSE)
     }))
     datb_plot		<- subset(datb.clust, datb.clust[,eval(units_name)] %in% tmp)		
     
     ## bases faceted by cluster
     p.c     		<- ggplot(data=datb_plot,aes(x = time, y = value))
     l			<- geom_line(aes_string(group = eval(units_name)), alpha = 0.2)
     l.2     	     <- geom_smooth(aes(group=1),alpha = 1.0, 
                                  size = 1, linetype = 2, se = FALSE, colour = "brown",
                                  method = "loess")
     f              <- facet_wrap(~cluster, scales = "fixed")
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
               axis <- labs(x = eval(x.axis.label), y = "normalized y units" )
          }else{ ## no labels
               axis <- labs(x = "time", y = "normalized y units" )
          }
     }
     dev.new()
     if( smoother == TRUE)
     {
          p.c     		<- p.c + l + l.2 + f + theme_bw() + axis + 
                                   theme(axis.text.x=element_text(angle=90, hjust=0))
     }else{ ## no smoother
          p.c          	<- p.c + l + f + theme_bw() + axis + 
                              theme(axis.text.x=element_text(angle=90, hjust=0))
     }
     suppressWarnings(print(p.c))
     
     ###############################################################################
     ## compose plot of fitted functions vs. actual data for randomly-selected units
     ###############################################################################
     
     d3             <- data.frame(units_label,y) ## y is N x T
     names(d3)      <- c(eval(units_name),eval(time_points)) ## even though month is date/numeric, names will be char
     dat_gca        <- merge(d3,map,by=eval(units_name),all.x=TRUE)
     dat_gca        <- melt(dat_gca,measure.vars = as.character(time_points), variable.name = "time",
                            value.name = "fit")
     dat_gca$time   <- as.numeric(dat_gca$time)
     
     if( !is.null(date_field))
     {
          dat_date            <- data.frame(sort(unique(time_points)),sort(unique(date_field)))
          names(dat_date)     <- c("time","date")
          dat_gca             <- merge(dat_gca,dat_date, by = "time", all.x=TRUE)
          dat_gca$time        <- as.Date(dat_gca$date)
          dat_gca$date        <- NULL
     }
     
        
     ##
     ## compose fit predictions
     ##
     
     ## build data.frame with N*T rows (where N is fast-moving)
     lo                  <- apply(bb,2,function(x){quantile(x,probs = 0.025)}) # 1 x N*T - N is fast
     hi                  <- apply(bb,2,function(x){quantile(x,probs = 0.975)}) # 1 x N*T - N is fast
     dat_gcp             <- data.frame(rep(d3[,eval(units_name)],times=T),
                                  rep(1:T,each=N),bb.hat,lo,hi)
     names(dat_gcp)      <- c(eval(units_name),"time","fit","lo","hi")
     ## merge in cluster locations to facet fit v data plot panels by cluster membership
     dat_gcp             <- merge(dat_gcp,map,by=eval(units_name),all.x=TRUE)
     
     if( !is.null(date_field) )
     {
          dat_gcp             <- merge(dat_gcp,dat_date, by = "time", all.x=TRUE)
          dat_gcp$time        <- as.Date(dat_gcp$date)
          dat_gcp$date        <- NULL
     }
     
     ##
     ## render growth curve plot
     ##
     
     ## randomly select records to plot
     if(single_unit == FALSE) ## plot the fit for a single experimental unit?
     {
          if(!is.null(num_plot))
          { 
               num_plot       <- floor(num_plot) ## ensure num_plot is an integer
          }else{ ## is.null(num_plot) = TRUE
               num_plot       <- 6
          } 
          ## number of units per cluster
          cluster_counts                <- as.vector(table(map$cluster))
          M                             <- length(cluster_counts) ## number of clusters
          ## map[,cluser] is a vector of length N of cluster assignments 
          ## that is sorted by cluster.  We sample num_plot from this, without replacement,
          ## to ensure we don't sample more units per cluster than there are members.
          ## This is effectively sampling with probability proportional to cluster size.
          units_to_sample               <- map[,eval(units_name)]
          cluster_labels                <- map[,"cluster"]
          sample_ids                    <- sort(sample(cluster_labels,num_plot,
                                                       replace=FALSE))
          table_clust                   <- table(sample_ids) 
          ## labels of clusers sampled
          clusters_sampled              <- as.numeric(names(table_clust)) 
          ## number of each cluster label sampled
          num_per_clust                 <- as.vector(table_clust)
          ## generate num_per_clust randomly sampled units with clusters_sampled
          M_to_sample                   <- length(num_per_clust) ## won't sample all M clusters
          plot_st                       <- vector("list",M_to_sample)
          for( m in 1:M_to_sample )
          {
               units_to_sample_m        <- units_to_sample[cluster_labels == clusters_sampled[m]]
               plot_st[[m]]             <- sample(units_to_sample_m, num_per_clust[m],replace=FALSE)
          } ## end loop m sampling num_per_clust[m] units within cluster, cluster_sampled[m]
          plot_st         <- unlist(plot_st)
     }else{ ## sample only a single unit
          plot_st         <- sample(map[,eval(units_name)],1,replace=FALSE)
     
     } ## end conditional statement on whether single unit or group of units to be sampled
     ## outputs unit_name ids to be sampled in plot_st
     
     plot_gca            <- subset(dat_gca, dat_gca[,eval(units_name)] %in% plot_st)
     plot_gcp            <- subset(dat_gcp, dat_gcp[,eval(units_name)] %in% plot_st)    
     ## compose plot
     p.t          		= ggplot(data=plot_gcp,aes(x=time, y = fit) )
     l				= geom_line(aes_string(group = eval(units_name)),size = 1.2,
                              alpha = 0.6, colour="#FF9999")
     l.2				= geom_point(data=plot_gca,aes_string(group = eval(units_name)),size=3,shape=1,
                         colour="black") ## Add the data
     l.3          		= geom_line(data=plot_gca,aes_string(group = eval(units_name)),lty="dotted") ## Add the data
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

     f				<- facet_wrap(as.formula(paste("cluster ~", units_name)), scales="fixed")
     if( credible == FALSE )
     {
          p.t          	<- p.t + l + f + axis + l.2 + l.3 + theme_bw() 
     }else{ ## plot credible intervals
          ci             <- geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.2)
          p.t            <- p.t + l + ci + f + axis + l.2 + l.3 + theme_bw() 
     }
     
     suppressWarnings(print(p.t))
     
     value <- fit <- NULL
     return(invisible(list(p.cluster = p.c, p.fit = p.t, dat.fit = dat_gcp, 
                           dat.cluster = datb_plot, map = map)))
     
} ## end function cluster_plot