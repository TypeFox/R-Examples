#' Plot a gitter dat file 
#' 
#' This function will plot a heatmap or bubble plot of a data.frame produced by \code{\link{gitter}} or a dat file saved to disk.
#' 
#' @keywords plot visualize display bubble heatmap
#' 
#' @export 
#' 
#' @param x The data.frame produced by \code{\link{gitter}} or the path to a dat file saved by \code{\link{gitter}}.
#' @param title Title of plot. Default is blank.
#' @param type Type of plot. "heatmap" for a heatmap, "bubble" for a bubble plot. Default is "heatmap".
#' @param low Color for the lower bound of colony sizes. Default is "turquoise".
#' @param mid Color of the middle value of colony sizes. Default is "black".
#' @param high Color for the upper bound of colony sizes. Default is "yellow".
#' @param show.text Logical indicating if text representation of colony sizes should be overlaid on the plot. Default is \code{TRUE}.
#' @param text.color Color of text if show.text is \code{TRUE}. Default is "white".
#' @param norm Logical indicating if colony sizes should be normalized by dividing colony sizes the middle mean of values and capping them between 0-2. Default is \code{TRUE}.
#' @param show.flags Logical indicating if dots should be overlaid on the plot for flagged colonies. Default is \code{TRUE}.
#' @param flag.color Color of flag dot if show.flags is \code{TRUE}. Default is "white".
#' @param ... Additional arguments. Not used.
#' 
#' 
#' @return a ggplot heatmap or bubble plot.
#' 
#' @examples
#' f = system.file("extdata", "sample.jpg.dat", package="gitter")
#' # Read in path as a gitter data object
#' g = gitter.read(f)
#' # Plot a heatmap
#' plot(g, type="heatmap")
#' # Show a bubble plot
#' plot(g, type="bubble", low="black", high="red")
plot.gitter <- function(x, title='', type='heatmap', low='turquoise', mid='black', high='yellow', 
                        show.text=F, text.color='white', norm=T, 
                        show.flags=T, flag.color='white', ...){
  dat = x
  #   if(is.character(dat)) dat = read.table(dat, stringsAsFactors=F, header=T)
  if(!is.data.frame(dat) & ! 'gitter' %in% class(dat)) stop('Argument must be a gitter data object')
  if(!type %in% c('heatmap', 'bubble')) stop('Invalid plot type. Use "heatmap" or "bubble"')
  
  if(length(dat) > 5 | length(dat) < 3) stop('Invalid number of columns for dat file')
  
  names(dat) = c('r', 'c', 's')
  dat.cs = s = NULL
  
  r = max(dat$r)
  c = max(dat$c)
  
  t = r:1
  names(t) = 1:r
  dat$r = t[as.character(dat$r)]
  
  srt = sort(dat$s)
  pmm = mean( srt[ (0.4*length(srt)) : (0.6*length(srt)) ], na.rm=T)
  if(norm){
    dat$s = dat$s / pmm
    pmm = 1
    #dat$s = nrm(dat$s, 0.5, 1.5)
    dat$s[dat$s < 0.5] = 0.5
    dat$s[dat$s > 1.5] = 1.5
  }
  
  if(length(dat) == 5){
    names(dat)[5] = 'flags'
    dat.cs = dat[dat$flags != "",] 
  }
  
  gitter_theme <- theme( axis.ticks = element_blank(), axis.ticks.margin=grid::unit(0,'cm') )
  if(type == 'heatmap'){
    p <- ggplot(dat, aes(x = c, y = r, fill = s)) + 
      geom_tile() +
      scale_fill_gradient2(midpoint=pmm, low=low, high=high, mid=mid) + 
      scale_x_discrete(expand = c(0, 0), limits = as.character(1:c)) +
      scale_y_discrete(expand = c(0, 0), limits = as.character(r:1)) + 
      coord_equal() +
      theme_bw() +
      labs(list(title = title, x = "Column", y = "Row", fill = "Size")) + 
      theme(legend.position="right",title=element_text(size=14,face="bold")) +
      gitter_theme
  }
  if(type == 'bubble'){
    p = ggplot(dat, aes(x = c, y = r)) + 
      geom_point(aes(x = c, y = r, size = s, colour = s),shape=16, alpha=0.80) +
      scale_colour_gradient2(midpoint=pmm, low=low, mid=mid, high=high) +
      scale_x_discrete(limits = as.character(1:c)) +
      scale_y_discrete(limits = as.character(r:1)) +
      coord_equal() +
      labs(list(title = title, x = "Column", y = "Row", color = "Size", size=""))+
      theme_bw() + 
      gitter_theme
  }
  if(show.flags & !is.null(dat.cs) & nrow(dat.cs) != 0){
    p = p + geom_point(aes(x=c, y=r), data=dat.cs, color=flag.color, size=1)
  }
  
  if(show.text) p = p + geom_text(color = text.color, aes(label=s), size=3)
  return(p)
}