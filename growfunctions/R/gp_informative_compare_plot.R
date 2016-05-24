#' Plot credible intervals for parameters to compare ignoring with weighting an informative sample
#'
#' Uses as input the output object from the gpdpgrow() and gmrfdpgrow() functions.
#'
#' @param objects A list of objects, either all outputs from gpdpgrow(), or all from gmrfdpgrow().
#'        \code{objects} includes a model estimated under ignoring the informativeness of the 
#'        sampling design and another that employs weighting to account for the informativeness.
#'        An additional object may be added to represent a separate "iid" (or non-informative)
#'        sample from the same population, which will typically be available if the dataset
#'        was generated as synthetic data using function, "gen_informative_sample()".  
#' @param objects_labels A character vector of length equal to \code{objects} that provides
#'        labels for each entry in \code{objects}.  Allowed entries in \code{objects_labels} 
#'        are \code{c("ignore","weight","iid")}, where "ignore" denotes a model that ignores
#'        the informativeness, while "weight" denotes a model that employs sampling weights and
#'        "iid" denotes a model run on a non-informative, iid sample from the same population.
#'        Defaults to \code{objects_labels = c("ignore","weight")}.       
#' @param map A list matrices, where each entry is produced 
#'        from \code{cluster_plot(object)$map}. It is comprised of
#'        unit labels and cluster assignments for each object in \code{objects}.
#'        The length of \code{map} must be equal to the length of \code{objects}.
#' @param units_name The label in each "map" matrix for the observation units. Will be the same as the
#'        \code{units_name} entry for the previously run, \code{cluster_plot()} function.
#' @param model A scalar character input indicating the estimation model for \emph{all} of the entries
#'        in \code{objects}. Allowable values for \code{model} are c("gp","gmrf"). Defaults to
#'        \code{model = "gp"}.
#' @param true_star An optional, \emph{P x M}  matrix, of true parameter location values, where
#'        \code{P} denotes the number of parameters per cluster and \code{M} denotes the number of
#'        clusters.  For example, in \code{model = "gp"} with a single, rational quadratic covariance,
#'        \code{P = 3} and if there are 3 clusters, then \code{M = 3}.  For a \code{model = "gmrf"},
#'        with a single covariance, \code{P = 1}.  
#' @param map_true An optional \code{data.frame} object with \code{n} rows, the size of the
#'        informative sample used for \code{c("ignore","weight")} objects that maps the
#'        \code{units_name} to a true cluster.  \code{map_true} must have 2 columns (and the
#'        rest are ignored), one must be named the same value as input for \code{units_name}.
#'        The second column must be named, \code{cluster}.  If the true values derive
#'        from running \code{gen_informative_sample()} as the source of the true values,
#'        one may just input the \code{map_obs} \code{data.frame} that is listed in the object
#'        returned by \code{gen_informative_sample()}.
#' @return A list object containing the plot of estimated functions, faceted by cluster,
#'     	and the associated \code{data.frame} object.
#'     \item{p.compare}{A \code{ggplot2} plot object}
#'     \item{dat.compare}{A \code{data.frame} object used to generate \code{p.compare}.}
#' @seealso \code{\link{gpdpgrow}}, \code{\link{gmrfdpgrow}}
#' @examples 
#' \dontrun{
#' library(growfunctions)
#' ## use gen_informative_sample() to generate an 
#' ## N X T population drawn from a dependent GP
#' ## By default, 3 clusters are used to generate 
#' ## the population.
#' ## A single stage stratified random sample of size n 
#' ## is drawn from the population using I = 4 strata. 
#' ## The resulting sample is informative in that the 
#' ## distribution for this sample is
#' ## different from the population from which 
#' ## it was drawn because the strata inclusion
#' ## probabilities are proportional to a feature 
#' ## of the response, y (in the case, the variance.
#' ## The stratified random sample over-samples 
#' ## large variance strata).
#' ## (The user may also select a 2-stage 
#' ## sample with the first stage
#' ## sampling "blocks" of the population and 
#' ## the second stage sampling strata within blocks). 
#' dat_sim        <- gen_informative_sample(N= 10000, 
#'                                 n = 500, T = 5,
#'                                 noise_to_signal = 0.1)
#'
#' y_obs                       <- dat_sim$y_obs
#' T                           <- ncol(y_obs)
#'
## estimate parameters using substitution method under 
#' an informative sampling design that inputs inclusion
#' probabilities, ipr
#' res_gp_w            <- gpdpgrow(y = y_obs, 
#'                                ipr = dat_sim$map_obs$incl_prob, 
#'                                n.iter = 5, n.burn = 2,  
#'                                n.thin = 1, n.tune = 0)
## plots of estimated functions, faceted by cluster 
#' and fit vs. data for experimental units
#' fit_plots_w         <- cluster_plot( object = res_gp_w,  
#'                            units_name = "establishment", 
#'                            units_label = dat_sim$map_obs$establishment, 
#'                            single_unit = FALSE, credible = TRUE )
#'
#' ## estimate parameters ignoring sampling design
#' res_gp_i            <- gpdpgrow(y = y_obs, 
#'                                n.iter = 5, n.burn = 2, 
#'                                n.thin = 1, n.tune = 0)
#' ## plots of estimated functions, faceted by cluster and fit vs. 
#' ## data for experimental units
#' fit_plots_i         <- cluster_plot( object = res_gp_i,  
#'                                     units_name = "establishment", 
#'                                     units_label = dat_sim$map_obs$establishment, 
#'                                     single_unit = FALSE, credible = TRUE )
#'
#' ## We also draw an iid (non-informative, exchangeable) 
#' ## sample from the same population to 
#' ## compare estimation results to our weighted 
#' ## (w) and unweighted/ignoring (i) models
#'
#' ## estimate parameters under an iid sampling design
#' res_gp_iid          <- gpdpgrow(y = dat_sim$y_iid, 
#'                                n.iter = 5, n.burn = 2,   
#'                                n.thin = 1, n.tune = 0)
#' ## plots of estimated functions, faceted by cluster and 
#' ## fit vs. data for experimental units
#' fit_plots_iid       <- cluster_plot( object = res_gp_iid,  
#'                            units_name = "establishment", 
#'                            units_label = dat_sim$map_iid$establishment, 
#'                            single_unit = FALSE, credible = TRUE )
#'
#' ## compare estimations of covariance parameter credible 
#' ## intervals when ignoring informativeness vs.
#' ## weighting to account for informativeness
#' objects                  <- map <- vector("list",3)
#' objects[[1]]             <- res_gp_i
#' objects[[2]]             <- res_gp_iid
#' objects[[3]]             <- res_gp_w
#' map[[1]]                 <- fit_plots_i$map
#' map[[2]]                 <- fit_plots_iid$map
#' map[[3]]                 <- fit_plots_w$map
#' objects_labels           <- c("ignore","iid","weight")
#' 
#' parms_plots_compare      <- informative_plot( objects = objects, 
#'                                 objects_labels = objects_labels,
#'                                 map = map, units_name = "establishment", 
#'                                 model = "gp",
#'                                 true_star = dat_sim$theta_star, 
#'                                 map_true = dat_sim$map_obs)
#' 
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases informative_plot
#' @export
informative_plot <- function(objects = NULL, objects_labels = c("ignore","weight"), map = NULL, 
                              units_name = NULL, model = "gp",
                              true_star = NULL, map_true = NULL)
{ ## begin function informative_plot()
     ## test inputs to make have minimal set-up for comparison plot of covariance parameters
     if( length(objects) < 2 )
     {
          stop("\nMust have at least 2 models in 'object' in order to perform comparison plot
               with other models.\n")
     } ## end condition on having at least 2 models to compare
     
     if( length(objects) != length(map) )
     {
          stop("\nThe number of 'objects' and the number of entries in 'map' must be equal
               and each entry in map must correspond to the same in 'objects'.\n")
     } ## ensure length of objects and map are the same
     
     if( !("ignore" %in% objects_labels) | !("weight" %in% objects_labels) )
     {
          stop("\nThe comparison must include an 'ignore' model that ignores the sampling
               design and another, 'weight', that uses weights to account for an informative
               sampling design.\n")
     } ## end conditional statement on inclusion of "weight" model
     
     if( (model != "gp") & (model != "gmrf") )
     {
          stop("\nInput for 'model' must be either 'gp' or 'gmrf'.\n")
     } ## check to ensure model is either "gp" or "gmrf"
     
     ## some dimensions
     N          <- nrow(map[[which(objects_labels == "ignore")]])
     if( model == "gp" )
     {
          P   <- ncol(objects[[which(objects_labels == "ignore")]]$Theta)/N
     }else{ ## gmrf
          P   <- ncol(objects[[which(objects_labels == "ignore")]]$Kappa)/N
     } ## end generating dimensions N, number of observed units and P, number of parms per unit
     
     if( is.null(true_star) )
     {
          ## more likely to have the right number of clusters if use weighting
          M               <- length(unique(map[[which(objects_labels == "weight")]]$cluster))
     }else{
          if( is.null(map_true) )
          {
               stop("\nMust input 'map_true' that maps units to true clusters in the case
                    want to co-plot the true values of cluster parameters.\n")
          } ## end condition to ensure that if true_star is input, then so is map_true
          M               <- ncol(true_star)
     }
     
     ## units name must represent the user input labels from cluster_plot(), not the numeric labels
     ## since we will have to use common labels to match units from iid to ignore objects
     if( is.null(units_name) ) 
     {
          units_name <- names(map[[which(objects_labels == "ignore")]])[3] ## the third column is always the labels
          ## if exclude units_label option from cluster_plot the 3rd column duplicates the 1st.
     } ## end conditional statement on assigning units_name
     
     ## capture overlapping units to use for plotting covariance parameters, facetted by cluster
     if( "iid" %in% objects_labels )
     {
          ## first, need to be sure that the user assigned labels to sampled units in both iid and ignore
          ## objects, so the units can be compared.
          if( all( map[[which(objects_labels == "iid")]][,1] == map[[which(objects_labels == "iid")]][,3] ) )
          {
               stop("\nYou must input labels for the sampled units in the 'map' objects
                    extracted from the cluster_plot() functions for 'iid' and c('ignore','weight') 
                    so that the units can be matched when performing plotting.  Remember that
                    because the c(ignore, weight) and iid are distinct samples, the included units in 
                    each will not perfectly overlap, so that overlap must be found to ensure
                    we're plotting parameter comparisons for the same units.\n")
          } ## ensure labels for units are being units when include 'iid' and 'ignore'.
          
          ## set indicator in map(which[[objects_labels == "weight"]]) where overlaps in units sampled
          ## map(which[[objects_labels == "iid"]]).
          map_overlap         <- subset(map[[which(objects_labels == "weight")]], 
                                        map[[which(objects_labels == "weight")]][,eval(units_name)] %in% 
                                             map[[which(objects_labels == "iid")]][,eval(units_name)])
          if( nrow(map_overlap) == 0 )
          {
               stop("\nCannot render plot because have no sampled units in common between
                    'ignore' and 'weight' objects, on the one hand,  and the 'iid' object,
                    on the other hand. To render plot, must remove 'iid' object or re-generate
                    iid sample and estimation to overlap with units from 'ignore' or 'weight'\n")
          }
          
          map[[which(objects_labels == "weight")]]$overlap     <- 0
          map[[which(objects_labels == "weight")]]$overlap[map[[which(objects_labels == "weight")]][,eval(units_name)] 
                                                   %in%  map_overlap[,eval(units_name)]] = 1
     } ## end conditional statement whether including iid sampling
     
     
     ##
     ## generate data.frame of covariance parameters for objects_labels = "weight"
     ##
     
     ## determine row indices of units to sample from map objects for 'ignore' and 'weight' based
     ## on either overlap with units in iid sample or just choose the units randomly, if don't
     ## include iid.  
     unit            <- vector("integer",M)
     if( "iid" %in% objects_labels )
     {
          ## sample a row number, not a unit name, since not all names included in informative sample
          obs_cand        <- which(map[[which(objects_labels == "weight")]]$overlap == 1) ## row indices in 1:nrow(map_obs)
          if( !is.null(true_star) )
          {
               map_true                 <- data.frame(map_true[,eval(units_name)], map_true$cluster)
               names(map_true)          <- c(eval(units_name),"cluster_true")
               map[[which(objects_labels == "weight")]] <- merge(map[[which(objects_labels == "weight")]],
                                                                 map_true, by = eval(units_name),
                                                                 all.x = TRUE, sort = FALSE)
               obs_cand_clust  <- map[[which(objects_labels == "weight")]][obs_cand,"cluster_true"]
          }else{## use clustering from weighted estimation (since has same data as "ignore")
               obs_cand_clust  <- map[[which(objects_labels == "weight")]][obs_cand,"cluster"]
          } ## end creating map linking informative sampled ids to cluster
          
          for( m in 1:M ) 
          {
               ## randomly draw a row number (which we will later tie to a unit for iid)
               unit[m]                <- sample(obs_cand[obs_cand_clust == m],1) 
          }
     }else{ ## not including iid, so just pick a random unit from each cluster
          rows_weight     <- 1:nrow(map[[which(objects_labels == "weight")]])
          for( m in 1:M )
          {
               ## which rows have cluster m
               pos_m                  <- which(rows_weight[map[[which(objects_labels == "weight")]]$
                                                                cluster == m]) 
               unit[m]                <- sample(pos_m,1)
          }
          
     } ## end sampling of row positions in map for "ignore" and "weight" (which use the same sample)
     
     
     ## use the row position to select covariance parameters for inclusion in plots
     sel            <- vector("list",M)
     for( m in 1:M )
     {
          ## res_gp_i uses same same y_obs as does re-sampling, so row indices are the same
          sel[[m]]            <- (unit[m]-1)*P+1:P # a P x 1 vector for each cluster, m
     }
     sel            <- unlist(sel) ## an integer vector
     
     ## build data.frame of sampled Theta values for estimation type = "ignore"
     if( model == "gp" )
     {
          Theta_w   <- objects[[which(objects_labels == "weight")]]$Theta[,sel]
     }else{ ## model == "gmrf"
          Theta_w   <- objects[[which(objects_labels == "weight")]]$Kappa[,sel]
     } ## end conditional statement on building estimation matrix of selected parameters
     nsamps         <- nrow(Theta_w)
     colnames(Theta_w) <- NULL
     parameter      <- rep(1:P,times = M)
     cluster        <- rep(1:M,each=P) ## 1:P,1:P,1:P,
     unit_vec       <- rep(unit,each=P)
     type           <- rep("weight",each=(P*M)) 
     dat_alt          <- data.frame(t(Theta_w),parameter,cluster,unit_vec,type)
     names(dat_alt)     <- c(1:nsamps,"parameter","cluster","unit","type")
     dat_alt		<- melt(dat_alt, id.vars=c("parameter","cluster","unit","type"), 
                      value.name = "draw",variable.name = "number")
     
     ##
     ## generate data.frame of covariance parameters for objects_labels = "ignore"
     ##
     if( model == "gp" )
     {
          Theta_i        <- objects[[which(objects_labels == "ignore")]]$Theta[,sel]
     }else{ ## model == "gmrf"
          Theta_i        <- objects[[which(objects_labels == "ignore")]]$Kappa[,sel]
     } ## end conditional statement on building estimation matrix of selected parameters
     nsamps         <- nrow(Theta_i)
     parameter      <- rep(1:P,times = M)
     cluster        <- rep(1:M,each=P) ## 1:P,1:P,1:P,
     unit_vec       <- rep(unit,each=P)
     type           <- rep("ignore",each=(P*M)) 
     dat_i            <- data.frame(t(Theta_i),parameter,cluster,unit_vec,type)
     names(dat_i)     <- c(1:nsamps,"parameter","cluster","unit","type")
     dat_i     	<- melt(dat_i, id.vars=c("parameter","cluster","unit","type"), 
                          value.name = "draw",variable.name = "number")
     dat_alt        <- rbind(dat_alt,dat_i)
     
     if( "iid" %in% objects_labels )
     {
          ##
          ## generate data.frame of covariance parameters for objects_labels = "iid"
          ##
          
          ## place estimates from an "iid" sample in the same
          ## data.frame
          sel_iid                  <- vector("list",M)
          unit_iid                 <- vector("integer",M) 
          est_weight               <- map[[which(objects_labels == "weight")]][,eval(units_name)]
          for( m in 1:M )
          {
               est_m               <- est_weight[unit[m]]
               ## row index
               unit_iid[m]         <- which(map[[which(objects_labels == "iid")]][,eval(units_name)] == est_m) 
               sel_iid[[m]]        <- (unit_iid[m]-1)*P+1:P # a P x 1 vector for each cluster, m
          }
          sel_iid            <- unlist(sel_iid) ## an integer vector
          
          ## build data.frame of sampled Theta values for estimation type = "iid"
          if( model == "gp" )
          {
               Theta_iid        <- objects[[which(objects_labels == "iid")]]$Theta[,sel_iid]
          }else{ ## model == "gmrf"
               Theta_iid        <- objects[[which(objects_labels == "iid")]]$Kappa[,sel_iid]
          } ## end conditional statement on building estimation matrix of selected parameters
          nsamps              <- nrow(Theta_iid)
          unit                <- rep(unit_iid,each=P)
          parameter           <- rep(1:P,times = M)
          cluster             <- rep(1:M,each=P) ## 1:P,1:P,1:P,
          type                <- rep("iid",each=(P*M)) 
          dat_iid             <- data.frame(t(Theta_iid),parameter,cluster,unit,type)
          names(dat_iid)          <- c(1:nsamps,"parameter","cluster","unit","type")
          dat_iid		     <- melt(dat_iid, id.vars=c("parameter","cluster","unit","type"), 
                                value.name = "draw",variable.name = "number")
          
          dat_alt             <- rbind(dat_alt,dat_iid)
          
          ## concatenate the 2 datasets and control order of plotting models
          type_order          <- data.frame(c("iid","ignore","weight"),c(1,2,3))
          names(type_order)   <- c("type","sort_order")     
     }else{ ## "iid" is not in objects_labels for plotting
          ## concatenate the 2 datasets
          type_order          <- data.frame(c("ignore","weight"),c(1,2))
          names(type_order)   <- c("type","sort_order")        
     }## end conditional statement to build data.frame for iid sample if "iid" is in objects_labels
     df                  <- merge(dat_alt,type_order,by="type",all.x = TRUE)
     
     
     ## insert true values
     if(!is.null(true_star))
     {
               dat_true        <- data.frame(matrix(true_star,(M*P),1),rep(1:P,times=M),rep(1:M,each=P))
               names(dat_true) <- c("value","parameter","cluster")  
     } ## end conditional statement on whether true parameters values are included for plotting
     
     
     ##
     ## render plot for covariance parameters credible intervals
     ##
     
     ## box-plots of by-factor posterior distributions
     f <- function(x) {
          r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
          names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
          r
     }
     p.ci	<- ggplot(data=df,aes(x = reorder(type,sort_order), y = draw, fill = factor(reorder(type,sort_order))))
     l	<- stat_summary(fun.data = f, geom="boxplot",alpha=0.3)
     f	<- facet_wrap(cluster~parameter,scales = "free") 
     ## m	<- stat_summary(fun.y=mean, geom="line", aes(group=1)) 
     if( model == "gp" )
     {
          axis     <- labs(x = "Parameter", y = expression(paste("95% CI  ", theta[ip])), 
                           fill = "Estimation Type") 
     }else{ ## model == "gmrf"
          axis     <- labs(x = "Parameter", y = expression(paste("95% CI  ", kappa[ip])), 
                           fill = "Estimation Type")
     }
     options <- theme(axis.text.x=element_text(angle=90, hjust=0))##,legend.position="none") 
     if( !is.null(true_star) ) ## incorporating true location values
     {
          l.2  <- geom_hline(data = dat_true, aes(yintercept = value), linetype = "dotted")#, size = 1)
          p.ci <- p.ci + l + l.2 + f + theme_bw() + axis + options
     }else{ ## true values, true_star, not included in plot
          p.ci <- p.ci + l + f + theme_bw() + axis + options
     }
     
     suppressWarnings(print(p.ci))
     return(invisible(list(p.compare = p.ci, dat.compare = df)))
     
     sort_order <- draw <- value <- NULL

} ## end function informative_plot()

