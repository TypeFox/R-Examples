#' Dot-plots of community presence/absence or abundance
#' 
#' @param x \code{\link{comparative.comm}} object
#' @param sites names of sites to plot (default: all); see examples
#' @param abundance make size proportional to species abundance
#' (default: FALSE)
#' @param dot.cex function to determine point size; see examples, this
#' isn't as terrible-sounding as it seems.
#' @param site.col colours to use when plotting sites; if not same
#' length as number of sites, only the first element is used (no
#' recycling)
#' @param fraction fraction of plot window to be taken up with
#' phylogeny; e.g., 3 (default) means phylogeny is 1/3 of plot
#' @param pch plotting character to be used for sites (see
#' \code{\link{pch}})
#' @param x.increment specify exact spacing of points along plot; see
#' examples
#' @param show.tip.label whether to plot species names on phylogeney
#' (default: \code{FALSE})
#' @param ... additional arguments passed to plotting functions
#' @details Take a look at the examples: this is (hopefully!) a lot
#' more straightforward than it might seem. Getting the right spacing
#' of dots on the phylogeny may take some playing around with the
#' \code{fraction} and \code{x.increment} arguments. It may seem a
#' little strange to set point size using a function, however, this
#' gives you much more flexibility and the ability to (usefully)
#' transform your data.
#' @return List containing plot.phylo information, as well as the used
#' x.adj values (compare with your \code{x.increment})
#' @seealso \code{\link{comparative.comm}} \code{\link{traitgram.cc}}
#' @author Will Pearse, Matt Helmus
#' @method plot comparative.comm
#' @importFrom ape tiplabels plot.phylo
#' @importFrom graphics plot
#' @export
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' plot(data)
#' plot(data, sites=c("AT", "BP"), fraction=1.5)
#' settings <- plot(data, sites=c("AT", "BP"), site.col=rainbow(2), fraction=1.5)
#' plot(data, sites=c("AT", "BP"), site.col=rainbow(2),
#' fraction=1.2, x.increment=settings$x.increment/4)
#' #dot.cex isn't as scary as it sounds...
#' plot(data, site.col=rainbow(2), fraction=1.2, abundance=TRUE, dot.cex=sqrt)
#' #...or other trivial variants...
#' abund.sqrt <- function(x) ifelse(x>0, sqrt(x), 0)
#' plot(data, sites=c("AT", "BP"), site.col=rainbow(2), fraction=1.2,
#' x.increment=settings$x.increment/4, abundance=TRUE, dot.cex=abund.sqrt)
#' plot(data, sites=c("AT", "BP"), site.col=rainbow(2), fraction=1.2,
#' x.increment=settings$x.increment/4, abundance=TRUE, dot.cex=function(x) sqrt(x))
plot.comparative.comm <- function(x, sites=NULL, abundance=FALSE, pch=20, dot.cex=NULL, site.col="black", fraction=3, x.increment=NULL, show.tip.label=FALSE, ...){
  #Assertions and argument checking
  data <- x
  if (!inherits(data, "comparative.comm"))
    stop("ERROR:", deparse(substitute(data)), "is not of class 'comparative.comm'")
  if(is.null(sites)){
    sites <- rownames(data$comm)
  } else {
    if(any(!sites %in% rownames(data$comm))){
      missing <- paste(sites[!sites %in% rownames(data$comm)], collapse=" ")
      stop("ERROR:", missing, "not site(s) in", deparse(substitute(data)))
    }
  }
  if(!is.null(dot.cex) & !is.function(dot.cex))
    stop("ERROR:", deparse(substitute(dot.cex)), "is not a function!")
  
  #Process data
  sites <- data$comm[sites, ]
  if(!abundance)
    sites <- sites > 0
  if(is.function(dot.cex))
    sites <- dot.cex(sites) else sites
  if(length(site.col) != nrow(sites))
    site.col <- rep(site.col[1], nrow(sites))
  
  #Plot phylogeny and dots
  display <- plot(data$phy, show.tip.label=show.tip.label, plot=FALSE, no.margin=TRUE, ...)
  plot(data$phy, x.lim=c(0, display$x.lim[2] * fraction), no.margin=TRUE, show.tip.label=show.tip.label, ...)
  if(!is.null(x.increment)){
    if(length(x.increment) != nrow(sites))
      stop("ERROR:", deparse(substitute(x.increment)), "'s length does not match the number of sites to be plotted") else x.adj <- x.increment
  } else {
    x.adj <- (display$x.lim[2] * fraction) / (nrow(sites) + 2)
    x.adj <- seq(from=x.adj/2, by=x.adj, length.out=nrow(sites))
  }
  for(i in seq(nrow(sites))){
      if(is.function(dot.cex)){
          tiplabels(tip=match(colnames(sites), data$phy$tip.label), col=site.col[i], adj=x.adj[i], pch=pch, cex=dot.cex(sites[i,]), ...)
      } else tiplabels(tip=match(colnames(sites), data$phy$tip.label), col=site.col[i], adj=x.adj[i], pch=pch, ...)
   }
  
  #Invisibly return
  output <- display
  output$x.increment <- x.adj
  invisible(output)
}
