plotGapfillSeries <- function(
  ##title<< Plot an overview of a the results of a SSA gapfilling (from a ncdf file)
  file.orig                ##<< object to plot: file name or file.con object linking to a ncdf file
  , file.filled = sub('[.]nc', '_gapfill.nc', file.orig) ##<< character string: name of the
                           ##   filled file.
  , data.orig = c()        ##<< array: Unfilled data. Can be supplied to prevent loading
                           ##   the data from a ncdf file again. This is read from
                           ##   'file.filled' if no value is given.
  , data.filled = c()      ##<< array: Filled data. Can be supplied to prevent loading
                           ##   the data from a ncdf file again. This is read from
                           ##   'file.filled' if no value is given.  
  , ...
)
  ##description<<
  ## This function plots some overview statistics of the results of a gapfilling run
  ## in a netCDF file,  i.e. the results of a call to gapfillNcdf().
  ##seealso<<
  ##\code{\link{gapfillSSA}}, \code{\link{gapfillNcdf}}
  ##\if{html}{\out{<img src="../doc/visualize_ncdf_demo.png" alt="image ..visualize_ncdf_demo should be here"/>}}\ifelse{latex}{}{}
{
  ##TODO facilitate datacube input
  ##TODO include plotNLines capabilites
  set.seed(12233)
  
  ## preparation
  con.orig   <- open.nc(file.orig)
  con.filled <- open.nc(file.filled)
  vars.all <- infoNcdfVars(con.filled)[,'name']
  vars.filled <- sub('_', '', sub('flag.orig$', '', vars.all[grepl('flag.orig$', vars.all)]))


  for (var.t in vars.filled) {
    assign(paste(var.t, '.orig', sep =''), var.get.nc(con.orig, var.t))
    assign(paste(var.t, '.filled', sep =''), var.get.nc(con.filled, var.t))
  }
  
 

  for (var.t in vars.filled) {
    
    if(interactive() && names(dev.cur()) == 'X11')
      dev.new()
    layout(matrix(c(1,2,3,3,4,5,6,7,8,9), byrow=TRUE, ncol=2),
           heights = c(1, 1))
    
    par(tcl = 0.2, mgp = c(1, 0, 0), mar = c(2, 0, 0, 2), oma = c(0, 2, 0, 0))

    filled.t = get(paste(var.t, 'filled', sep='.'))
    orig.t   = get(paste(var.t, 'orig', sep='.'))
    if (sum(!is.na(orig.t)) == 0) {
      for (i in 1:9) 
        plot.new()
      text(userCoords(0.1,0.9), labels =  var.t, cex = 2)
      next
    }
    
    breaks = seq(min(c(orig.t, filled.t), na.rm = TRUE),
      max(c(orig.t, filled.t), na.rm = TRUE), length.out = 100)
    hst.orig    <- hist(get(paste(var.t, 'orig', sep='.')), plot = FALSE, breaks = breaks)
    hst.filled  <- hist(get(paste(var.t, 'filled', sep='.')), plot = FALSE, breaks = breaks)
    plotBG(rgb(.5,.5,.5))
    plot(hst.filled, xlim = range(c(hst.orig$mids, hst.filled$mids)),
         ylim = c(0, max( range(c(hst.orig$counts, hst.filled$counts)))),
         col = 'red', xlab = '', yaxt ='n', main = '')
    par(new=TRUE)  
    plot(hst.orig, xlim = range(c(hst.orig$mids, hst.filled$mids)),
         ylim = c(0, max( range(c(hst.orig$counts, hst.filled$counts)))), xlab = '',
         main = '', col = 'black', yaxt ='n')
    box()
    text(userCoords(c(0.8,0.9),c(0.9, 0.9)), labels =  c('filled', 'orig'),
         col = c('red', 'black'), cex = 2)
    text(userCoords(0.1,0.9), labels =  var.t, cex = 2)
    
    hst.filled  <- logHist(filled.t, breaks = breaks, col = 'red',
                           pch = 16, xlab = '', main = '')
    par(new = TRUE)
    logHist(orig.t, breaks = breaks, ylim = hst.filled$ylim,
            cex = 1.1, xlab = 'ratio of missing values per grid point', main = '')
    
    plot(filled.t, pch = 16, cex = 0.5, col = 'red')
    points(orig.t, pch = 16, cex = 0.5, col = 'black')
    
    n.points <- ceiling(length(filled.t)/50)
    index    <- rep(1:ceiling(length(filled.t)/n.points), each = n.points)[1:length(rep(1:length(filled.t), times = ))]
    na.pp    <- aggregate(orig.t, list(index), function(x)sum(is.na(x)))
    na.pp    <- na.pp[na.pp$x != n.points,]
    na.pp    <- na.pp[-dim(na.pp)[1],]
    na.pp    <- na.pp[na.pp$x < n.points/1.5,]
    n.plots  <- min(c(6,dim(na.pp)[1]))

    inds.plot<- na.pp$Group.1[sort(sample(order(na.pp$x, decreasing = TRUE)[1:min(c(n.plots, 30))], n.plots))]
    abline(v = (inds.plot - 1) * n.points, lty = 2, col ='green', lwd = 2)
    for (i in 1:n.plots) {
      index.t = 1:n.points + ((inds.plot[i]-1) * n.points)
      plot(filled.t[index.t], pch = 16, cex = 0.5, col = 'red')
      points(orig.t[index.t], pch = 16, cex = 0.5, col = 'black')
    }    
  }  
}
