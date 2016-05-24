
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2015-04-07 16:09:44 +0200 (Tue, 07 Apr 2015) $
# $Rev: 339 $
# ----------------------------------------------------------------

##########################################################################
##----------------------------ScatterPlotWux----------------------------##
##########################################################################

plot.wux.df <- function(x,
                        var1.name = "delta.air_temperature",
                        var2.name = "perc.delta.precipitation_amount",
                        subreg.subset = NULL,
                           season.subset = NULL,
                           boxplots = TRUE,
                           label.only.these.models = NULL,
                           highlight.models = NULL,
                           no.text = FALSE, 
                           vert.box.col = "cyan",
                           horiz.box.col = "coral",
                           zero.line.col = "gray80",
                           median.line.col = "black",
                           draw.legend = TRUE,
                           draw.seperate.legend = FALSE,
                           draw.median.lines = TRUE,
                           use.rainbow.colors = TRUE,
                           xlim = NULL,
                           ylim = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           main = NULL,
                           out.file.directory = NULL,
                           out.file.name = NULL,
                           copyright = FALSE,
                           ...) {

  ## Draws a climate change scatterplot for two parameters of the WUX dataframe
  ##
  ## Args:
  ##   x: WUX dataframe (wux.df object) obtained from
  ##              'models2wux(user.input)'
  ##   var1.name: Character name of 1st parameter in the WUX dataframe
  ##              (default is temperature change)
  ##   var2.name: Character name of 2nd parameter in the WUX dataframe
  ##              (default is precipitation change)
  ##   season.subset: Vector of seasons to be plotted (e.g. c("MAM", "DJF"))
  ##   subreg.subset: Vector of subregions to be plotted (eg. c("ACQWA", "GAR"))
  ##   boxplots: Boolean. TRUE if marginal boxplots for the two input parameters
  ##     should be plotted
  ##   label.only.these.models: Character vector of model names (acronyms) to
  ##     be labeled in the scatterplot
  ##   no.text: Boolean. FALSE if no models should be labeled
  ##   vert.box.col: Color character for vertical boxplot
  ##   horiz.box.col: Color character for horizontal boxplot
  ##   zero.line.col = Color of zero line,
  ##   median.line.col = Color of median line
  ##   draw.median.lines:  Draw median lines for both parameters. Default is
  ##                       TRUE.
  ##   draw.legend: Boolean. TRUE if a legend indicating the GCMs and emission
  ##                scenarios should be plotted.
  ##   draw.seperate.legend: Boolean. Should legend with GCMs be plotted on a
  ##                         seperate screen?  Default is FALSE.
  ##                         Draws legend even if draw.legend is set FALSE}.
  ##   use.rainbow.colors: Boolean. If TRUE, use rainbow() colors, else
  ##                       a custom color palette with discrete colors is used.
  ##                       Default is TRUE.
  ##   xlim: Range vector for 1st parameter (e.g. c(0, 10))  
  ##   ylim: Range vector for 2nd parameter (e.g. c(-100, 100))
  ##   xlab: Label of x-axis (i.e. 1st parameter)
  ##   ylab: Label of y-axis (i.e. 2nd parameter)
  ##   main: Main title
  ##   out.file.directory: String of the directory where the plots are exported
  ##     (e.g. "/tmp/plots/"
  ##   out.file.name: Preceding string of the file name of the plots
  ##     (e.g. "scatterplot"). Scatterplot will be saved as
  ##     out.file.name_subreg_season.eps
  ##   copyright: Boolean. If a copyright message should be plotted.
  ##              Default is FALSE.
  ##   '...': Further arguments to be passed `plot`, such as graphical
  ##     parameters (see 'par').
  ##
  ## History:
  ##   2010-10-21 | original code (thm)
  ##   2011-07-14 | loop over subregions and seasons (geh)
  ##   2011-10-21 | now we can plot legend seperatly (thm)
  ##   2011-10-21 | copyright added (thm)
  ##   2011-12-02 | added optional custom color palette (msu)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)
  ##   2015-04-03 | !!! changed  from function 'ScatterplotWux' to method
  ##                plot.wux.df !!! (thm)

   datain.df <- x
  ## extract specified seasons and regions
  ## default subregions
  if (is.null(subreg.subset)) { 
    subreg.subset <- levels(datain.df$subreg)
  }
  ## default seasons
  if (is.null(season.subset)) { 
    season.subset <- levels(datain.df$season)
  } 
  ## extract from dataframe
  data.df <- datain.df[datain.df[["season"]] %in% season.subset &
                       datain.df[["subreg"]] %in% subreg.subset, ]
  data.df <- droplevels(data.df)

  ## do write plot to harddrive?
  write.to.drive.list <-
    DoWriteToDiskQuestionmark(out.file.directory, out.file.name)
  do.write.plot.to.disk <- write.to.drive.list$do.write.plot.to.disk
  out.file.name <- write.to.drive.list$out.file.name
  out.file.directory <- write.to.drive.list$out.file.directory

  ## loop over seasons and subregions
  season.levels <- levels(data.df$season)
  subreg.levels <- levels(data.df$subreg)
  for (ii in subreg.levels) {
    for (jj in season.levels) {

      ## plot setup
      if (do.write.plot.to.disk){
        file.path <- paste(out.file.directory, out.file.name,
                           "_", ii, "_", jj, ".eps", sep = "")
        postscript(file.path, width = 6, height = 6, horizontal = F,
                   paper="special")
      } else {
        dev.new()
      }

      ## subset season and region of dataframe
      subset.df <- data.df[data.df[["subreg"]] == ii &
                           data.df[["season"]] == jj, ]

      gcm.level <- levels(subset.df$gcm)
      emission.level <- levels(subset.df$em.scn)
      model.names <- as.character(subset.df$acronym)
      ## erase model names not to be labeled
      if (!is.null(label.only.these.models)) {
        model.names[!(model.names %in% label.only.these.models)] <- ""
      }

      ## get data vectors
      ## TODO(thm): not only reading in dataset, but vectors and lists as well
      ## (then acronym- and gcm vector have to be passed)
      x.data <- subset.df[[var1.name]]
      y.data <- subset.df[[var2.name]]

      ## set default x and y ranges
      if (is.null(xlim)) {
        tmp.xlim <-  c(min(x.data, na.rm = TRUE), max(x.data, na.rm = TRUE))
      } else {
        tmp.xlim <- xlim
      }

      if (is.null(ylim)) {
        tmp.ylim <-  c(min(y.data, na.rm = TRUE), max(y.data, na.rm = TRUE))
      } else {
        tmp.ylim <- ylim
      }

      ## set default xlab and ylab
      if (is.null(xlab)) {
        tmp.xlab <- var1.name  
      } else {
        tmp.xlab <- xlab
      }

      if (is.null(ylab)) {
        tmp.ylab <- var2.name
      } else {
        tmp.ylab <- ylab
      }
      
      ## split device in 3 and draw boxplots
      if (boxplots == TRUE) {
        ## defining outer margins and split device
        par(oma = c(1,1,4,2) + 0.1)
        nf <- layout(matrix(c(1, 0, 3, 2), 2, 2, byrow = TRUE),
                     widths = c(7, 1), heights = c(1, 7), TRUE)
        ## horizontal boxplot (var1.name)
        par(mar = c(0, 4, 0, 0) + 0.1)
        boxplot(x.data, axes = FALSE, ylim = xlim, space = 0, horizontal = TRUE,
                col = horiz.box.col)
        ## vertical boxplot (var2.name)
        par(mar = c(5, 0, 0, 0) + 0.1)
        boxplot(y.data, axes = FALSE, ylim = ylim, space = 0,
                col = vert.box.col)
        ## margins for scatterplot
        par(mar = c(5, 4, 0, 0) + 0.1)
        print.outer.title <- TRUE
      } else {
        print.outer.title <- FALSE
      }
      
      ## Scatterplot initialization
      plot(x.data, y.data, pch="", xlab = tmp.xlab, ylab = tmp.ylab,
           xlim = xlim, ylim = ylim, ...)
      Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = 0.5)
      abline(h = 0, v = 0, col = zero.line.col, lty = 2)
      if (draw.median.lines == TRUE) {
        abline(v = median(x.data), col = median.line.col, lty = 1)
        abline(h = median(y.data), col = median.line.col, lty = 1)
      }

      ## point, triangle, square, diamond
      emission.pch <- c(19, 17, 15, 18, 21, 22, 23, 24, 25)

      ## define a color palette
      if ((use.rainbow.colors == FALSE) && (length(gcm.level) > 17)) {
        cat("\nThere are too many GCMs in this dataset, not sufficient\n",
            "colors available - fall back on rainbow palette!\n")
        use.rainbow.colors <- TRUE
      }
      if (use.rainbow.colors == TRUE ) {
        colors <- rainbow(length(gcm.level))
      } else {
        colors.avail <- c(
                          "#800000","#004000","#000040","#46230A","#800080",
                          "#D90000","#008000","#000080","#68340E","#BFBF00",
                          "#FF8080","#00D900","#0000D9","#A8744E","#80FF80",
                          "#8080FF","#C6A38A"
                          )
        colors <- colors.avail[1:length(gcm.level)]
      }

      # draw points in scatterplot for each GCM in another color
      counter = 1
      for (kk in gcm.level) {
        gcm.df <- subset.df[subset.df[["gcm"]] == kk, ]
        for (ll in seq(along = emission.level)) {
          em.df <- gcm.df[gcm.df[["em.scn"]] ==  emission.level[ll], ]
          points(em.df[[var1.name]], em.df[[var2.name]], col = colors[counter],
                 pch = emission.pch[ll])
        }
        counter <- counter + 1
      }

      ## highlight models
      if (!is.null(highlight.models)){
        highlited.models.idx <- subset.df$acronym %in% highlight.models
        x <- subset.df[[var1.name]][highlited.models.idx]
        y <- subset.df[[var2.name]][highlited.models.idx]
        points(x,y, col = "red"
               , cex = 3, lwd = 3,...)
      }

      ## plot acronym text to plot points
      if (no.text == FALSE) {
        text(x.data, y.data, model.names, cex = 0.55, pos = 4)
      }
        
      ## print main title
      if (is.null(main)) {
        tmp.main <- paste("region: ", ii, ", season: ", jj, sep = "")
      } else {
        tmp.main <- paste(main, "\n", "region: ", ii,
                                ", season: ", jj, sep = "")
      }      
      title(main = tmp.main, cex.main = 1.2, font = 2,
            outer = print.outer.title, line = 0.05)

      ## plot copyright
      if (copyright == TRUE)
        plot.copyright()

      ## draw emission scenario legend if there are more than one
      if (length(emission.level) > 1)     
        legend("topright" , emission.level, title = "Scenarios",
               pch = emission.pch, cex = 0.7, inset = 0.01)
      
      ## plot GCM legend; here case seperation if we want to draw legend
      ## seperatly
      if (draw.legend || draw.seperate.legend) { #do we want to draw at all?
        if (draw.seperate.legend) {
          if (do.write.plot.to.disk){
            ## start new postscript device
            dev.off()
            file.path <- paste(out.file.directory, out.file.name, "_", ii,
                               "_LEGEND.eps", sep = "")
            postscript(file.path, width = 8, height = 6, horizontal = F)
          } else {
            dev.new()
          }
          ## make new plot for legend anyway
          plot.new()
        }
        legend("topleft" , gcm.level, title = "GCM", pch = 19,
               col = colors, cex = 0.7, inset = 0.01)
      }

      ## if writing to postscript file, then device has to be closed
      if (do.write.plot.to.disk)
        dev.off()
      
    }       
  }
}
##########################################################################


##########################################################################
##-----------------------------HistPlotWux------------------------------##
##########################################################################

hist.wux.df <- function(x,
                        datain2.df = NULL,
                        var.name = "delta.air_temperature",
                        subreg.subset = NULL,
                        season.subset = NULL,
                        plot.density = TRUE,
                        ## plot.lmm.density = FALSE,
                        hist1.col = "red",
                        hist2.col = "blue",
                        bw = "nrd0",
                        kernel = "gaussian",
                        mark.df = NULL,
                        plot.legend = FALSE, 
                        xlim = NULL,
                        ylim = NULL,
                        xtick.number = 10,
                        ytick.number = 10,
                        xminor.tick = FALSE,
                        yminor.tick = FALSE,
                        xlab = NULL,
                        ylab = "Probability Density",
                        main = NULL,
                        out.file.directory = NULL,
                        out.file.name = NULL,
                        copyright = FALSE,
                        ...) {

  ## Compares the histograms and kernel density estimates of two distributions
  ##   of the same parameter of the WUX dataframe 
  ##
  ## Args:
  ##   x: 1st WUX dataframe
  ##   datain2.df: 2nd WUX dataframe. Is optionally. If indicated, the
  ##     distributions of datain1.df$var.name and datain2.df$var.name
  ##     will be compared. If not, only datain1.df$var1.name will be plotted.
  ##   var.name: Character name of the parameter in the WUX dataframe
  ##   season.subset: Vector of seasons to be plotted (e.g. c("MAM", "DJF"))
  ##   subreg.subset: Vector of subregions to be plotted (eg. c("ACQWA", "GAR"))
  ##   plot.density: Boolean. TRUE if kernel density estimates should be plotted
  ##   hist1.col: Color name of the 1st histogram (e.g. "red", "blue") 
  ##   hist2.col: Color name of the 2nd histogram (e.g. "red", "blue")
  ##   bw: The smoothing bandwidth to be used in 'density'. 
  ##   kernel: A character string giving the smoothing kernel to be used in 'density'.
  ##     This must be one of "gaussian", "rectangular", "triangular", "epanechnikov",
  ##     "biweight", "cosine" or "optcosine" with default "gaussian"
  ##   mark.df: Subset of WUX dataframe indicating the models to be marked
  ##   plot.legend: Boolean. TRUE if a plot legend indicating the models of mark.df and
  ##     sample size should be plotted
  ##   xlim: Range vector for the x-axis (e.g. c(0, 10))  
  ##   ylim: Range vector for the y-axis (e.g. c(-100, 100))
  ##   xtick.number: Number of ticks for the x-axis
  ##   ytick.number: Number of ticks for the y-axis
  ##   xminor.tick: Boolean. TRUE if minor ticks for the x-axis should be plotted
  ##   yminor.tick: Boolean. TRUE if minor ticks for the y-axis should be plotted
  ##   xlab: Title for the x-axis
  ##   ylab: Title for the y-axis. Default is "Probability Density"
  ##   main: Main title
  ##   out.file.directory: String of the directory where the plots are exported
  ##     (e.g. "/tmp/plots/"
  ##   out.file.name: Preceding string of the file name of the plots
  ##     (e.g. "histogram"). Histogram will be saved as out.file.name_subreg.eps 
  ##   copyright: Boolean. If a copyright message should be plotted.
  ##              Default is FALSE.
  ##   '...': Further arguments passed to `hist`
  ##
  ## History:
  ##   2011-07-15 | original code (geh)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)
  ##   2015-04-03 | !!! changed  from function 'HistPlotWux' to method
  ##                hist.wux.df !!! (thm)

  datain1.df <- x
  ## extract specified seasons and regions
  ## default subregions
  if (is.null(subreg.subset)) { 
    subreg.subset <- levels(datain1.df$subreg)
  }
  ## default seasons
  if (is.null(season.subset)) { 
    season.subset <- levels(datain1.df$season)
  } 
  ## extract from 1st dataframe
  data1.df <- datain1.df[datain1.df[["season"]] %in% season.subset &
                         datain1.df[["subreg"]] %in% subreg.subset, ]
  data1.df <- droplevels(data1.df)
  ## extract from 2nd dataframe
  if (!is.null(datain2.df)) {
    data2.df <- datain2.df[datain2.df[["season"]] %in% season.subset &
                           datain2.df[["subreg"]] %in% subreg.subset, ]
    data2.df <- droplevels(data2.df)
  }

  ## do write plot to harddrive?
  write.to.drive.list <-
    DoWriteToDiskQuestionmark(out.file.directory, out.file.name)
  do.write.plot.to.disk <- write.to.drive.list$do.write.plot.to.disk
  out.file.name <- write.to.drive.list$out.file.name
  out.file.directory <- write.to.drive.list$out.file.directory

  ## loop over subregions
  season.levels <- levels(data1.df$season)
  subreg.levels <- levels(data1.df$subreg)
  for (ii in subreg.levels) {
    ## plot setup
    if (do.write.plot.to.disk){
      file.path <- paste(out.file.directory, out.file.name, "_", ii, ".eps", sep = "")
      postscript(file.path, width = 8, height = 6, horizontal = F, paper = "special") 
    } else {
      dev.new()
    }
    
    par(mfcol = c(2, 2), mar = c(1.5, 1.5, 3.5, 1), oma = c(1.3, 1.3, 0, 0))
    ## loop over seasons
    for (jj in season.levels) {
 
      ## subset season and region of 1st dataframe
      subset1.df <- data1.df[data1.df[["subreg"]]== ii &
                             data1.df[["season"]]== jj, ]
      subset1.hist <- hist(subset1.df[[var.name]], plot = FALSE, ...)

      ## set default xlim and ylim for 1 histogram
      if (is.null(xlim)) {
        tmp.xlim <- c(2 * min(subset1.hist$breaks),
                      2 * max(subset1.hist$breaks))
      } else {
        tmp.xlim <- xlim
      }

      if (is.null(ylim)) {
        tmp.ylim <- c(0, 1.5 * round(max(subset1.hist$density), digits = 2))
      } else {
        tmp.ylim <- ylim
      }

      ## set default xlim and ylim for 2 histograms
      if (!is.null(datain2.df)) {

        ## subset season and region of 2nd dataframe
        subset2.df <- data2.df[data2.df[["subreg"]] == ii &
                               data2.df[["season"]] == jj, ]
        subset2.hist <- hist(subset2.df[[var.name]], plot = FALSE, ...)
        
        if (is.null(xlim)) {
          tmp.xlim <- c(2 * min(min(subset1.hist$breaks),
                                min(subset2.hist$breaks)),
                        2 * max(max(subset1.hist$breaks),
                                max(max(subset2.hist$breaks))))
        } else {
          tmp.xlim <- xlim
        }

        if (is.null(ylim)) {
          tmp.ylim <- c(0, 1.5 * round(max(max(subset1.hist$density),
                                           max(subset2.hist$density)),
                                       digits = 2))
        } else {
          tmp.ylim <- ylim
        }

      }
      
      ## plot 1st histogram
      hist(subset1.df[[var.name]], col = hist1.col, density = 35, angle = 45,
           prob = T, axes = FALSE, main = "", xlab = "", ylab = "",
           xlim = tmp.xlim, ylim = tmp.ylim, ...)
      lines(rep(quantile(subset1.df[[var.name]], c(0.05)) ,2),
            c(tmp.ylim[1], tmp.ylim[2]), col = hist1.col, lty = 2)
      lines(rep(quantile(subset1.df[[var.name]], c(0.5)) ,2),
            c(tmp.ylim[1], tmp.ylim[2]), col = hist1.col, lty = 2)
      lines(rep(quantile(subset1.df[[var.name]], c(0.95)) ,2),
            c(tmp.ylim[1], tmp.ylim[2]), col = hist1.col, lty = 2)
      
      ## plot legend for sample size 1st histogram
      if (is.null(datain2.df) & plot.legend == TRUE) {
        legend("topright" ,
               paste("n = ", as.character(length(subset1.df[[var.name]])),  sep = ""),
               title = "Sample Size", pch = 20, col = hist1.col, cex = 0.6, inset = 0.08)   
    }
      
      ## plot 2nd histogram
      if (!is.null(datain2.df)) {
        hist(subset2.df[[var.name]], col = hist2.col, density = 35, angle = -45,
             prob = TRUE, add = TRUE, ...)
        lines(rep(quantile(subset2.df[[var.name]], c(0.05)) ,2),
              c(tmp.ylim[1], tmp.ylim[2]), col = hist2.col, lty = 2)
        lines(rep(quantile(subset2.df[[var.name]], c(0.5)) ,2),
              c(tmp.ylim[1], tmp.ylim[2]), col = hist2.col, lty = 2)
        lines(rep(quantile(subset2.df[[var.name]], c(0.95)) ,2),
              c(tmp.ylim[1], tmp.ylim[2]), col = hist2.col, lty = 2)

        if (plot.legend == TRUE) {
          ## plot legend for sample size 1st and 2nd histogram
          legend("topright" , paste(c("n = ", "n = "),
                                    c(as.character(length(subset1.df[[var.name]])),
                                      as.character(length(subset2.df[[var.name]]))),
                                    sep = ""),
                 title = "Sample Size", pch = 20,
                 col = c(hist1.col, hist2.col), cex = 0.6, inset = 0.08)
        }
      }

      ## plot kernel density estimates
      if (plot.density == TRUE) {
        ## 1st density estimate
        subset1.density <- density(subset1.df[[var.name]],
                                   bw = bw, kernel = kernel, na.rm = TRUE)
        lines(subset1.density, col = hist1.col, lwd = 1.5)

        ## 2nd density estimate
        if (!is.null(datain2.df)) {
          subset2.density <- density(subset2.df[[var.name]],
                                     bw = bw, kernel = kernel, na.rm = TRUE)
          lines(subset2.density, col = hist2.col, lwd = 1.5)
        }

        ## if (plot.lmm.density == TRUE) {
        ##   ## fit a linear mixed model for unbalanced data
        ##   form <- formula(paste(var.name, " ~ 1  +  (1|gcm) +  (1|rcm)", sep =""))
        ##   lmm <- lmer(form, data = subset1.df)
        ##   ## mean
        ##   avg <- lmm@fixef
        ##   ## sum over all VCs and claim thats our marginal var
        ##   sd <- sqrt(sum(as.numeric(summary(lmm)@REmat[,3])))
        ##   ## and plot the normal distr density
        ##   x   <-subset1.density$x
        ##     ## x   <- seq(xlim[1], xlim[2], 0.05)
        ##   y   <- dnorm(x,mean=avg, sd=sd)
        ##   ## subset1.density$x <- x
        ##   subset1.density$y <- y
        ##   lines(subset1.density, col = "red", lwd = 4)
        ## }
      }
   

      ## plot model marks and according legend
      if (!is.null(mark.df)) {
        mark.colors <- rainbow(length(levels(mark.df$acronym)))
        for (kk in 1:length(levels(mark.df$acronym))) {
          mark.num <- mark.df[mark.df[["acronym"]] == as.character(levels(mark.df$acronym)[kk]) &
                              mark.df[["subreg"]] == ii &
                              mark.df[["season"]] == jj, ][[var.name]]
          points(mark.num, 0, pch = 18, col = mark.colors[kk], cex = 2.5)
        }
        ## plot legend for the marks
        if (plot.legend == TRUE) {
            legend("topleft" , as.character(levels(mark.df$acronym)),
                   title = "RCM", pch = 18, col = mark.colors, cex = 0.6, inset = 0.08)
          }
      }
      
      ## plot axes
      ## lower x axis
      axis(1, at = seq(tmp.xlim[1], tmp.xlim[2],
                (tmp.xlim[2] - tmp.xlim[1]) / (xtick.number + 1)),
           labels = FALSE, pos = 0, tcl = 0.5)
      if (xminor.tick == TRUE) {
        axis(1, at = seq(tmp.xlim[1], tmp.xlim[2],
                  (tmp.xlim[2] - tmp.xlim[1]) / ((xtick.number + 1) * 2)),
             labels = FALSE, pos = 0, tcl = 0.25)
      }

      ## left y axis
      axis(2, at = seq(tmp.ylim[1], tmp.ylim[2],
                (tmp.ylim[2] - tmp.ylim[1]) / ((ytick.number + 1))),
           labels = FALSE, pos = tmp.xlim[1], tcl = 0.5)
      if (yminor.tick == TRUE) {
        axis(2, at = seq(tmp.ylim[1], tmp.ylim[2],
                  (tmp.ylim[2] - tmp.ylim[1]) / ((ytick.number + 1) * 2)),
             labels = FALSE, pos = tmp.xlim[1], tcl = 0.25)      
      }

      ## upper x axis
      axis(3, at = seq(tmp.xlim[1], tmp.xlim[2],
                (tmp.xlim[2] - tmp.xlim[1]) / (xtick.number + 1)),
           labels = FALSE, pos = tmp.ylim[2], tcl = 0.5)
      if (xminor.tick == TRUE) {   
        axis(3, at = seq(tmp.xlim[1], tmp.xlim[2],
                  (tmp.xlim[2] - tmp.xlim[1]) / ((xtick.number + 1) * 2)),
              labels = FALSE, pos = tmp.ylim[2], tcl = 0.25)
       }

      ## right y axis
      axis(4, at = seq(tmp.ylim[1], tmp.ylim[2],
                (tmp.ylim[2] - tmp.ylim[1]) / (ytick.number + 1)),
           labels = FALSE, pos = tmp.xlim[2], tcl = 0.5)
      if (yminor.tick == TRUE) {   
        axis(4, at = seq(tmp.ylim[1], tmp.ylim[2],
                  (tmp.ylim[2] - tmp.ylim[1]) / ((ytick.number + 1) * 2)),
              labels = FALSE, pos = tmp.xlim[2], tcl = 0.25)
       }

      ## label axes
      ## lower x axis
      axis(1, at = seq(tmp.xlim[1], tmp.xlim[2],
                (tmp.xlim[2] - tmp.xlim[1]) / (xtick.number + 1)),
           labels = as.character(seq(tmp.xlim[1], tmp.xlim[2],
             (tmp.xlim[2] - tmp.xlim[1]) / (xtick.number + 1))),
           cex.axis = 0.8, line = -1.2, tick = FALSE, outer = FALSE)      

      ## left y axis
      axis(2, at = seq(tmp.ylim[1], tmp.ylim[2],
                (tmp.ylim[2] - tmp.ylim[1]) / (ytick.number + 1)),
           labels = as.character(format(seq(tmp.ylim[1], tmp.ylim[2],
             (tmp.ylim[2] - tmp.ylim[1]) / (ytick.number + 1)), digits = 2)),
           cex.axis = 0.8, line = -1.4, tick = FALSE, outer = FALSE, las = 1)   

      ## print main title
      if (is.null(main)) {
        tmp.main <- paste("region: ", ii, ", season: ", jj, sep = "")
      } else {
        tmp.main <- paste(main, "\n", "region: ", ii,
                                ", season: ", jj, sep = "")
      }      
      title(main = tmp.main, cex.main = 1.0, font = 2, line = 0.2)

      ## print axes labels
      mtext(ylab, side = 2, line = 1.2, cex = 0.8)
      mtext(xlab, side = 1, line = 0.8, cex = 0.8)

    }

    ## print copyright
    if (copyright == TRUE)
      plot.copyright()

     ## if writing to postscript file, then device has to be closed
    if (do.write.plot.to.disk)
      dev.off()
  }
  
}
##########################################################################


##########################################################################
##----------------------------AnovaBarplotWux---------------------------##
##########################################################################

## method corresponding to AnovaBarplotWux for simpler user interface
plot.wux.aov <- function(x,
                            ss.relative = TRUE,
                            subreg.subset = NULL,
                            cex.names = 1.2,
                            cex.lab = 1.2,
                            legend.text = NULL,
                            sd.text = TRUE,
                            sd.unit = "",
                            ylim = NULL,
                            ylab = NULL,
                            main = NULL,
                            out.file.directory = NULL,
                            out.file.name = NULL,
                            copyright = FALSE,
                           ...) {

  ## Barplots of the 'WuxANOVA' results displaying the relative or absolute
  ## contribution of the factors to the overall variance
  ##
  ## Args:
  ##   x: List result obtained from `aovWux`
  ##   ss.relative: Boolean. TRUE if the relative contribution should be calculated
  ##   subreg.subset: Vector of subregions to be plotted (e.g. c("ACQWA", "GAR"))
  ##   ylim: Range vector for the y-axis (e.g. c(-100, 100))
  ##   cex.names: Expansion factor for numeric axis labels for `barplot`
  ##   cex.lab: Expansion factor for axis names (bar labels) for `barplot`
  ##   ylab: Title for the y-axis
  ##   main: Main title
  ##   legend.text: String vector of the factors
  ##   sd.text: Boolean. TRUE if the overall standard deviation should be displayed
  ##   sd.unit: Character string of the  standard deviation unit (e.g. "K")
  ##   out.file.directory: String of the directory where the plots are exported
  ##     (e.g. "/tmp/plots/"
  ##   out.file.name: Preceding string of the file name of the plots
  ##     (e.g. "anova"). Barplot will be saved like out.file.name_subreg.eps 
  ##   copyright: Boolean. If a copyright message should be plotted.
  ##              Default is FALSE.
  ##   '...': Further arguments to be passed to `barplot`
  ##
  ## History:
  ##   2011-07-15 | original code (geh)
  ##   2011-10-24 | copyright statment added (thm)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)

  datain.list <- x
  
  ## extract information on subregions and seasons from list
  ## default subregions
  if (is.null(subreg.subset)) { 
    subreg.names <- unique(unlist(lapply(strsplit(names(datain.list),";"), "[", 1)))
  } else {
    subreg.names <- paste("subreg = ", subreg.subset, sep = "")
  }

  names.list <- strsplit(names(datain.list), ";")
  season.names <- unique(unlist(lapply(names.list, "[", 2)))

  if (length(season.names) == 4){
    season.names <- c("season = DJF", "season = MAM",
                      "season = JJA", "season = SON")
  }
  
  ## do write plot to harddrive?
  write.to.drive.list <-
    DoWriteToDiskQuestionmark(out.file.directory, out.file.name)
  do.write.plot.to.disk <- write.to.drive.list$do.write.plot.to.disk
  out.file.name <- write.to.drive.list$out.file.name
  out.file.directory <- write.to.drive.list$out.file.directory

  ## loop over subregions
  for (ii in subreg.names) {

    ## plot setup
    if (do.write.plot.to.disk){
      file.path <- paste(out.file.directory, out.file.name, "_",
                         stringr::str_trim(unlist(strsplit(ii, "="))[2]),
                         ".eps", sep = "")
      postscript(file.path, width = 8, height = 6, horizontal = F)
    } else {
      dev.new()
    }
    
    par(mar = c(1.5, 5, 3.5, 1), oma = c(1, 1, 0, 0))

    ## create the plot matrix
    plot.matrix <- c()
    season.vector <- c()
    sd.vector <- c()
    ## loop over seasons
    for (jj in season.names) {
      subset.name <- paste(ii, ";", jj, sep = "")
      tmp.df <- data.frame(summary(datain.list[[subset.name]])[[1]])

      dummy.season <- unlist(strsplit(jj, "= "))[2]
      season.vector <- cbind(season.vector, dummy.season)

      ## calculate the standard deviations
      dummy.sd <- sd(datain.list[[subset.name]]$model[[1]])
      sd.vector <- cbind(sd.vector, dummy.sd)

      ## calculate the relative or absolute contribution to the
      ## overall variance
      if (ss.relative == TRUE) {
        plot.vector <- tmp.df$Sum.Sq / sum(tmp.df$Sum.Sq) * 100
      } else {
        plot.vector <- tmp.df$Sum.Sq / sum(tmp.df$Sum.Sq) * (dummy.sd)^2
      }

      plot.matrix <- cbind(plot.matrix, plot.vector)
    }

    tmp.df$factors <- stringr::str_trim(row.names(tmp.df))
    attr(plot.matrix, "dimnames") <- list(tmp.df$factors, season.vector)

    ## set default values for the bar plot
    ## default colors
    col <- rainbow(dim(plot.matrix)[1])

    ## default y-axis range
    if (is.null(ylim)) {
      if (ss.relative == TRUE) {
        tmp.ylim <- c(0, 100)  
      } else {
        tmp.ylim <- c(0, 1.5 * max(plot.matrix))
      }
    } else {
      tmp.ylim <- ylim
    }

    ## default y-axis title 
    if (is.null(ylab)) {
      if (ss.relative == TRUE) {
        tmp.ylab <- "Variance [%]"  
      } else {
        tmp.ylab <- "Variance"
      }
    } else {
      tmp.ylab <- ylab
    }

    ## default labels for the legend
    if (is.null(legend.text)) {
      legend.text <- tmp.df$factors} 

    ## bar plot
    barplot <- barplot(plot.matrix, beside = TRUE, col = col,
                       ylim = tmp.ylim, border = TRUE, ylab = tmp.ylab,
                       cex.names = cex.names, cex.lab = cex.lab, ...)
    Hmisc::minor.tick(nx = 0, ny = 2, tick.ratio = 0.5)
    
    ## plot legend
    legend("topleft" ,legend.text, title = "Factor", pch = 15, col = col,
           cex = 0.8, inset = 0.01)
    box()
    
    ## add text 
    if (sd.text == TRUE) {
      x.pos <- grconvertX(barplot, from = "user", to = "npc")
      y.pos <- grconvertY(max(plot.matrix), from = "user", to = "npc")
      usr <- par("usr")
      par(usr = c(0, 1, 0, 1))

      text(as.vector(x.pos), 0.05, format(round(as.vector(plot.matrix), digits = 1)),
           font = 2, cex = 0.9, srt = 90)

      text(colMeans(matrix(x.pos, nrow = dim(plot.matrix)[1], ncol = dim(plot.matrix)[2])),
           y.pos + 0.05, paste("SD:", format(round(sd.vector, digits = 1)), sd.unit, sep = " "),
           font = 4, cex = 0.9, adj = c(0.5, 1))
    }

    ## print main title
    if (is.null(main)) {
      tmp.main <- paste("region:", unlist(strsplit(ii, "="))[2], sep = "")
    } else {
      tmp.main <- paste(main, "\n" ,"region:", unlist(strsplit(ii, "="))[2], sep = "")
    }      
    title(main = tmp.main, cex.main = 1.2, font = 2, line = 0.8)

    ## print copyright
    if (copyright == TRUE)
      plot.copyright()

    ## if writing to postscript file, then device has to be closed
    if (do.write.plot.to.disk)
      dev.off()
    
  }

}

###########################################################################


##########################################################################
##---------------------------AnnualCycleplotWux-------------------------##
##########################################################################

plotAnnualCycle <- function(datain.df,
                               var.name = NULL,
                               subreg.subset = NULL,
                               season.subset = NULL,
                               plot.quantiles = NULL,
                               quantile.method = 7,
                               mark.df = NULL,
                               plot.legend = FALSE,
                               cex.names = 1.2,
                               cex.lab = 1.2,
                               ylab = NULL,
                               main = NULL,
                               out.file.directory = NULL,
                               out.file.name = NULL,
                               copyright = FALSE,
                               ...) {

  ## Plots the annual cycle of indicated models and the box-whisker
  ## plots of the underlying model distribution
  ##
  ## Args:
  ##   datain1.df: WUX dataframe
  ##   var.name: Character name of the parameter in the WUX dataframe
  ##   season.subset: Vector of seasons to be plotted (e.g. c("MAM", "DJF"))
  ##   subreg.subset: Vector of subregions to be plotted (e.g. c("ACQWA", "GAR"))
  ##   plot.quantiles: 5 element vector indictaing the quantiles to be plotted
  ##     (e.g. c(0.02,0.25,0.5,0.75,0.98))
  ##   quantile.method: An integer between 1 and 9 selecting one of the nine
  ##     quantile types in `quantiles` with default 7 
  ##   mark.df: Subset of WUX dataframe indicating the models to be marked
  ##   plot.legend: Boolean. TRUE if a plot legend indicating the models of mark.df and
  ##     sample size should be plotted
  ##   cex.names: Expansion factor for numeric axis labels in `bxp`
  ##   cex.lab: Expansion factor for axis names (bar labels) in `bxp`
  ##   ylab: Title for the y-axis. Default is "Probability Density"
  ##   main: Main title
  ##   out.file.directory: String of the directory where the plots are exported
  ##     (e.g. "/tmp/plots/"
  ##   out.file.name: Preceding string of the file name of the plots
  ##     (e.g. "annualcycle"). Histogram will be saved as out.file.name_subreg.eps 
  ##   copyright: Boolean. If a copyright message should be plotted.
  ##              Default is FALSE.
  ##   '...': Further arguments to be passed to `bxp` (e.g. notch = FALSE)
  ##
  ## History:
  ##   2011-07-15 | original code (geh)
  ##   2011-10-24 | copyright added (thm)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)

  ## extract specified seasons and regions
  ## default subregions
  if (is.null(subreg.subset)) { 
    subreg.subset <- levels(datain.df$subreg)
  }
  ## default seasons
  if (is.null(season.subset)) { 
    season.subset <- levels(datain.df$season)
  } 
  ## extract seasons and subregions from dataframe
  data.df <- datain.df[datain.df[["season"]] %in% season.subset &
                       datain.df[["subreg"]] %in% subreg.subset, ]
  data.df <- droplevels(data.df)

  season.levels <- levels(data.df$season)  

  if (length(season.levels) == 4){
    season.levels <- c("DJF", "MAM", "JJA", "SON")
  }
  if (length(season.levels) == 12){
    season.levels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                       "Aug", "Sep", "Oct", "Nov", "Dec")
  }

  subreg.levels <- levels(data.df$subreg)

  ## do write plot to harddrive?
  write.to.drive.list <-
    DoWriteToDiskQuestionmark(out.file.directory, out.file.name)
  do.write.plot.to.disk <- write.to.drive.list$do.write.plot.to.disk
  out.file.name <- write.to.drive.list$out.file.name
  out.file.directory <- write.to.drive.list$out.file.directory

  ## loop over subregions
  for (ii in subreg.levels) {
    ## plot setup
    if (do.write.plot.to.disk) {
      file.path <- paste(out.file.directory, out.file.name, "_", ii, ".eps", sep = "")
      postscript(file.path, width = 8, height = 6, horizontal = F)
    } else {
      dev.new()
    }
    
    par(mar = c(1.5, 4, 3.5, 1), oma = c(1, 1, 0, 0))  
    ## subregion subset of dataframe
    subset.df <- data.df[data.df[["subreg"]] == ii, ]

    ## check if quantiles are specified by user and initialize
    ## output variables
    if (!is.null(plot.quantiles)) {
      quantiles.matrix <- vector()
      outlier.values <- vector()
      outlier.group <- vector()
    }
      
    ## create the list required for bxp() and calculate user
    ## specified quantiles and according outliers
    season.list <- list()
    for (jj in c(1:length(season.levels))) {
      act.season <- season.levels[jj]
      season.list[[act.season]] <- subset.df[subset.df[["season"]] == act.season, ][[var.name]]
      if (!is.null(plot.quantiles)) {
        calc.quantiles <- quantile(season.list[[act.season]], plot.quantiles,
                                   type = quantile.method)

        quantiles.matrix <- cbind(quantiles.matrix, calc.quantiles)
        outlier.pos <- which((season.list[[act.season]] < calc.quantiles[1]) |
                             (season.list[[act.season]] > calc.quantiles[length(calc.quantiles)]))
        outlier.values <- c(outlier.values, season.list[[act.season]][outlier.pos])
        outlier.group <- c(outlier.group, rep(jj, length(outlier.pos)))
      }    
    }

    boxplot.list <- boxplot(season.list, plot = FALSE)
    if (!is.null(plot.quantiles)) {
      boxplot.list$stats <- quantiles.matrix
      boxplot.list$out <- outlier.values
      boxplot.list$group <- outlier.group
    }
      
    ## plot box-whiskers for the seasons
    bxp(boxplot.list, cex.names = cex.names, cex.lab = cex.lab, ...)
    Hmisc::minor.tick(ny = 2, tick.ratio = 0.5)

    ## add zero lines
    abline(h = 0, lty = 2)   

      ## add lines for specified models
      ## default colors
    if (!is.null(mark.df)) {
      ## mark.color <- rainbow(length(levels(mark.df$acronym)))
      n.models <- length(levels(mark.df$acronym))
      mark.color <- rainbow(n.models)

      for (kk in 1:length(levels(mark.df$acronym))) {
        mark.subset <- mark.df[mark.df[["acronym"]] == as.character(
                                        levels(mark.df$acronym)[kk]) & 
                               mark.df[["subreg"]] == ii, ]
        mark.subset$season <- gdata::reorder.factor(mark.subset$season, new.order =
                                     c("Jan", "Feb", "Mar", "Apr", "May",
                                       "Jun", "Jul", "Aug", "Sep", "Oct",
                                       "Nov", "Dec"))
        ## mark.subset$season <- gdata::reorder.factor(mark.subset$season,
        ##                                             new.order =
        ##                                             c("DJF", "MAM", "JJA", "SON"))
        
        lines(mark.subset[[var.name]][order(mark.subset$season)], lwd = 5,
              col = mark.color[kk])
        points(mark.subset[[var.name]][order(mark.subset$season)], pch = 23,
               cex = 1.2, col = mark.color[kk])       
      }

      ## plot legend for the marks
      if (plot.legend == TRUE) {
        legend("topleft" , as.character(levels(mark.df$acronym)),
               title = "RCM", pch = 23, col = mark.color, cex = 0.7, inset = 0.01)
      }
    }
#    stop()

    ## print main title
    if (is.null(main)) {
      tmp.main <- paste("region: ", ii, sep = "")
    } else {
      tmp.main <- paste(main, "\n" , "region: ", ii, sep = "")
    }      
    title(main = tmp.main, cex.main = 1.2, font = 2, line = 0.8)

    ## print y-axis title 
    if (is.null(ylab)) {
      tmp.ylab <- var.name
    } else {
      tmp.ylab <- ylab
    }
    mtext(tmp.ylab, side = 2, line = 2.5, cex = 1.2)

    ## print copyright
    if (copyright == TRUE)
      plot.copyright()
    
    ## if writing to postscript file, then device has to be closed
    if (do.write.plot.to.disk)
      dev.off()
    
  }
}
##########################################################################


## Helper Functions:


DoWriteToDiskQuestionmark <- function(out.file.directory, out.file.name) {
  ## Decide whether to write plots to harddisk or not.
  ## If any "out.file.name" is given, the plot will be stored either in
  ## location of "out.file.directory" or in the current working directory.
  ## If "out.file.directory" is given without "out.file.name", the plot will
  ## be stored and its filename will be very random.
  ## If both "out.file.directory" and "out.file.name" are NULL, this function
  ## decides not to write the plot to disk.
  ## 
  ## If this function decides to write the plot to the harddrive, proper
  ## "out.file.directory "and "out.file.name" will be returned.
  ##
  ## Args:
  ##   out.file.directory: Character. The directory in which the plot will
  ##                       be written to. If NULL passed, either no plot will
  ##                       be written to drive, or (if out.file.name is not
  ##                       NULL) the current working directory will be used as
  ##                       out.file.directory.
  ##   out.file.name: Character. The basename of the filename for the plot
  ##                  on disk. 
  ##
  ## Returns:
  ##   List of boolean variable (do.write.plot.to.disk = TRUE -> plot it
  ##   do.write.plot.to.disk = FALSE print to screen) and out.file.directory
  ##   and out.file.name.
  ##
  ## History:
  ##   2011-09-02: original code (thm)
  
  if (is.null(out.file.directory) & is.null(out.file.name)) {
    do.write.plot.to.disk <- FALSE
  } else if (is.null(out.file.directory) & !is.null(out.file.name)) {
    do.write.plot.to.disk <- TRUE
    out.file.directory <- paste(getwd(), "/", sep = "")
  } else if (!is.null(out.file.directory) & is.null(out.file.name)) {
    do.write.plot.to.disk <- TRUE
    out.file.name <- paste("plot_", paste(sample(c(rep(0:9,each=5),
                                                   LETTERS,letters),8,
                                                 replace=TRUE),collapse=''),
                           sep = "")
  } else if (!is.null(out.file.directory) & !is.null(out.file.name)) {
    do.write.plot.to.disk <- TRUE
  }
  return(list(do.write.plot.to.disk = do.write.plot.to.disk,
              out.file.directory = out.file.directory,
              out.file.name = out.file.name))
}



















