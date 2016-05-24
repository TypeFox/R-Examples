plotDecompSeries <- function(
  ##title<< Visualize/plot an overview of a SSA decomposed ncdf file.
  file.orig                ##<< object to plot: file name or file.con object linking
                           ##   to the original ncdf file
  , file.decomp = sub('[.]nc', '_specdecomp.nc', file.orig) ##<< object to plot:
                           ##   file name or file.con object linking to the decomposed ncdf file.
  , file.plot = ''         ##<< character string: name of the file containing the
                           ##   plots. If not given, a plot window is opened.
  , ...
)
##description<<
## This function plots an visualisation of a SSA decomposed ncdf file.  
##\if{html}{\out{<img src="../doc/visualize_ncdf_demo.png" alt="image ..visualize_ncdf_demo should be here"/>}}\ifelse{latex}{}{}
{
  ##TODO facilitate datacube input
  ##TODO include plotNLines capabilites
  set.seed(12233)
  
  ## preparation
  con.orig    <- open.nc(file.orig)
  con.decomp  <- open.nc(file.decomp)
  vars.all    <- infoNcdfVars(con.decomp)[,'name']
  vars.decomp <- as.vector(na.omit(vars.all[infoNcdfVars(con.decomp)[,'2.dim'] == 'spectral_bands' & infoNcdfVars(con.decomp)[,'1.dim'] == 'time']))
  n.bands     <- infoNcdfDims(con.decomp)[infoNcdfDims(con.decomp)[, 'name']== 'spectral_bands', 'length']


  ## load data
  for (var.t in vars.decomp) {
    assign(paste(var.t, '.orig', sep =''), var.get.nc(con.orig, var.t))
    assign(paste(var.t, '.decomp', sep =''), var.get.nc(con.decomp, var.t))
  }
  
  if (is.element('time', infoNcdfVars(con.orig)[,'name'])) {
    x.vals <- convertDateNcdf2R(con.orig)
  } else {
    x.vals <- 1:length(var.get.nc(con.orig, vars.decomp[1]))
  }  

  for (var.t in vars.decomp) {
    printStatus(paste('Plotting variable ', var.t, sep = ''))
    n.plots = 5
    if (interactive() && is.element(names(dev.cur()), c('X11', 'null device')) && nchar(file.plot) == 0) {
      dev.new()
    } else if ( nchar(file.plot) > 0) {
      png(filename = paste(file.plot, '_', var.t, '.png', sep = ''), width = 3200, height = 1600, type = 'cairo',
          pointsize = 48)
    }
    
    mean.window = (var.get.nc(con.decomp, 'borders.up') + var.get.nc(con.decomp, 'borders.low'))/2
    borders.string = paste(var.get.nc(con.decomp, 'borders.low'),  var.get.nc(con.decomp, 'borders.up'), sep = '-')

    decomp.t = get(paste(var.t, 'decomp', sep='.'))
    
    orig.t   = get(paste(var.t, 'orig', sep='.'))
    plots.full <- mean.window > length(orig.t)/20
    layout.mat = matrix(c(rep(1:(sum(plots.full)+1),n.plots)), ncol = n.plots)
    layout.final <- rbind(layout.mat, matrix(1:(n.plots * sum(!plots.full)) + max(layout.mat), ncol = n.plots, byrow = TRUE))
    layout(layout.final, heights = c(2, rep(1, n.bands)))
    par(tcl = 0.2, mgp = c(1, 0, 0), mar = c(1, 0, 0, 2), oma = c(0, 2, 2, 0))

    
    if (sum(!is.na(decomp.t)) == 0) {
      for (i in 2:max(layout.final)) 
        plot.new()
      mtext(var.t, cex = 2, outer = TRUE, side = 3)
      next      
    }
    plot(x.vals, orig.t, pch = 16, cex = 0.5, col = 'black', xlab = '', ylim = range(c(orig.t, decomp.t), na.rm = TRUE))
    mtext(var.t, cex = 2, outer = TRUE, side = 3)

    plots.start = sort(sample(which(!is.na(decomp.t))[which(!is.na(decomp.t)) < (length(orig.t) - max(mean.window[!plots.full] * 10)) ], n.plots))
    abline(v = x.vals[plots.start], lty = 2, col = 'red', lwd = 2)

    
    for (i in 1:n.bands) {
      points(x.vals, decomp.t[,i], pch = 16, cex = 0.5, col = i + 1)
    }
    
    for (band.t in which(plots.full)) {
      plot(x.vals, orig.t, ylim = range(c(orig.t, decomp.t[, band.t]), na.rm = TRUE), xlab = '')
      text(userCoords(0.1, 0.8), borders.string[band.t], col = 'red', cex = 1.5)
      points(x.vals, decomp.t[, band.t], pch = 16, cex = 0.5, col = band.t + 1)
    }
    
      
    for (band.t in which(!plots.full)) {
      series.length = mean.window[band.t] * 10 
      for (j in 1:n.plots) {
        index = plots.start[j] + (0:series.length)
        plot(x.vals[index], orig.t[index], ylim = range(c(orig.t[index], decomp.t[index, band.t] ), na.rm = TRUE), xlab = '')
        if (j == 1)
          text(userCoords(0.1*n.plots, 0.8), borders.string[band.t], col = 'red', cex = 1.5)

        points(x.vals[index], decomp.t[index, band.t], pch = 16, cex = 0.5, col = band.t + 1)
      }
    }

    if ( nchar(file.plot) > 0)
      dev.off()

  }
  
  ##value<<
  ##nothing is returned.
}

