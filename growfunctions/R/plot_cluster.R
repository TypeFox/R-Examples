#' Plot estimated functions, facetted by cluster numbers, for a known clustering
#'
#' An internal function of the \code{growfunctions} package.
#'
#' @param y An \code{N x T} matrix for a collection of functions.
#' @param H An \code{N x 1} with entries in \code{1,...,M} of cluster assignments for the \code{N}
#'        units of \code{y} under a known clustering.
#' @param sort An optional boolean input on whether to sort the cluster-indexed plot panels 
#'        of function by size of cluster.  Defaults \code{sort = FALSE}.
#' @param sample_rate An optional numeric value in (0,1] indicating percent of functions to 
#'        randomly sample within each cluster to address over-plotting.  Defaults to 1.
#' @param y.axis.label An optional text label for y-axis. Defaults to \code{"function values"}.
#' @param smoother An optional scalar boolean input indicating whether to co-plot a smoother line 
#'                  through the functions in each cluster.
#' @param fade An optional numeric input in \code{(0,1)} indicating the degree of fading to apply to the
#'        plots of functions in each cluster-indexed panel.  Defaults to \code{fade = 0.2}.
#' @param cluster_order An optional numeric vector of length \code{M}, the number of clusters, 
#'        indicating the order, from-to-right, for plotting the cluster-indexed panels.
#' @param plot_render An optional boolean input indicating whether to render the plot.
#'        Defaults to \code{plot_render = TRUE}.             
#' @return A list object containing the plot of estimated functions, faceted by cluster,
#'     	and the associated \code{data.frame} object.
#'     \item{p.basis}{A \code{ggplot2} plot object}
#'     \item{map}{A \code{data.frame} object listing the unit and associated cluster membership.}
#' @seealso \code{\link{gpdpgrow}}, \code{\link{gmrfdpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases plot_cluster
#' @export

plot_cluster     <- function(y,H,sort=FALSE,sample_rate = 0.05,
                             y.axis.label = NULL,
                             smoother = TRUE,
                             fade = 0.2, cluster_order = NULL, plot_render = TRUE)
{
     ## convert cluster assignment vector, H, to a list of length M, where each 
     ## element contains the member establishment id's for the cluster associated to that element
     N                             <- nrow(y)
     T                             <- ncol(y)
     M                             <- length(unique(H))
     cluster          		     <- vector("list", M)
     for(m in 1:M) 
     {
          if(is.null(cluster_order))
          {
               cluster[[m]]      		     <- which(H == m) 
          }else{
               cluster[[m]]                  <- which(H == cluster_order[m])
          }
          names(cluster[[m]]) 	     <- NULL
     }
     
     ## create subject-cluster map matrix to use for plotting (instead of list object, lastS)
     c.sizes               <- sapply(cluster,length)
     if(!sort)
     {
          clusterstoplot      <- 1:M
     }else{
          clusterstoplot     	<- sort(c.sizes,decreasing = TRUE,index.return=TRUE)$ix[1:length(cluster)]
     }
     
     map			     <- vector(mode="list",length = length(clusterstoplot))
     
     for(i in 1:length(clusterstoplot))
     {
          cluster.i			<- cluster[[clusterstoplot[i]]] 
          map[[i]]     		<- as.data.frame(cbind(cluster.i,i),stringsAsFactors = FALSE)
          names(map[[i]]) 	<- c("establishment","cluster")
     }
     
     map			          <- do.call("rbind",map)
     map$establishment	     <- as.numeric(map$establishment)
     
     ## create plot data.frame
     y.hat          		     <- matrix(y,(N*T),1,byrow=FALSE) ## load by column, N is fast index
     establishment                 <- rep(1:N,times=T)
     month                         <- rep(1:T,each=N)
     dat.b				     <- data.frame(y.hat,month,establishment)
     names(dat.b)	          	<- c("value","time","establishment")
     map                           <- map[order(map$establishment),]
     datb.clust			     <- merge(dat.b,map,all.x=TRUE,by="establishment",sort = FALSE)
     
     ## sample records
     rate          	<- sample_rate
     tmp			<- split(datb.clust,list(datb.clust$cluster))
     tmp			<- unlist(sapply(tmp,function(x){
          tot_recs  	<- length(unique(x$establishment))
          u_recs	     <- sort(unique(x$establishment))
          inc_recs	     <- sample(u_recs,round(rate*tot_recs),replace = FALSE)
     }))
     datb_plot		<- subset(datb.clust, establishment %in% tmp)
     
     ## bases faceted by cluster
     p          	<- ggplot(data=datb_plot,aes(x = time, y = value))
     l			<- geom_line(aes(group = establishment), alpha = fade)
     
     if(is.null(y.axis.label))
     {
          axis      <- labs( x = "time", 
                             y = substitute(expression(paste(Delta," (Simulated) Employment Counts"))) )
     }else{
          axis          	<- labs( x = "time", y = eval(y.axis.label) )
     }
     
     if(smoother == TRUE) ## plot emphasizes overall theme
     {
          l.2                 <- geom_smooth(aes(group=1),method = "loess", alpha = 0.5, 
                                           linetype = 2, se = FALSE, colour = "brown")
          f                   <- facet_wrap(~cluster, scales = "fixed")
          p.basis     		<- p + l + l.2 + f + theme_bw() + axis + 
                     theme(axis.text.x=element_text(angle=90, hjust=0))
     }else{ ## plot emphasizes each curve
          f              <- facet_grid(cluster~., scales = "fixed")
          p.basis          	<- p + l + f + theme_bw() + axis + 
                     theme(axis.text.x=element_text(angle=90, hjust=0))
     }
     
     if(plot_render)
     {
          suppressWarnings(print(p.basis))
     }
     
     
     ## dev.new()
     value <- NULL
     
     return(list(map = map, p.basis = p.basis))
}