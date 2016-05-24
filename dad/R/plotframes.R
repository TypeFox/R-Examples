plotframes <-
function(x, y, font.size=12, layout = NULL)
{ 
  # x:          data frame of the variables which is displayed on the
  #             abscissa.
  # y:          data frame of the variables which are displayed on the
  #             ordinate.
  # font.size:  size of the characters in the strips of the graphs.
  # layout:     numeric vector of length 2 or 3 defining the layout of the
  #             graphical device, as in "xyplot" (package "lattice"). It is
  #             either c(i, j) or c(i, j, k):
  #               - i: number of lines
  #               - j: number of columns
  #               - k: number of pages (if omitted, it is set to as many as is
  #                    required to plot all the panels, see "xyplot")
  #             If omitted: layout=c(3, 3)
  
  # Convert x and y into data frames (if they are vectors or matrices)
  x <- as.data.frame(x)
  y <- as.data.frame(y)

  # Build a 4-column-data frame:
  #   1) factor: groups (levels of the factor)
  #   2) numeric: stacked x values (variables on the abscissa)
  #   3) factor: names of x variables
  #   4) numeric: will contain successively each one of y variables (ordinate)
  xy <- stack(x); names(xy) <- c("x", "variable")
  xy <- data.frame(nom = rep(rownames(x), ncol(x)), xy, y = 0)
  lev <- levels(xy$variable)
  
  # Layout of the graphical device
  if (is.null(layout))
    {layout = c(3, 3)
    }
    
  # Function defining the order of the graphs on the graphical device
  my.packet.panel = function(layout, condlevels, page, row, column, ...)
    { tlayout <- layout[c(2, 1, 3)] # switch row and column
      packet.panel.default(tlayout, condlevels, page=page, row=column, column=row, ...)
    }

  for (j in 1:ncol(y))  {
    # Graph of 'x[, k] ~ y[, j]' for each column x[, k] of x,
    # and computation of the correlation of each variable x[, k] with y[, j]
    # (these correlations will be displayed in the strips of each graph made by
    # 'xyplot()')
    xy$y = rep(y[,j], ncol(x))
    correl <- numeric(0)
    for (l in lev)
      {correl <- c(correl, round(cor(xy[xy$variable == l, "x"], xy[xy$variable == l, "y"]), 2))
      }
    levels(xy$variable) <- paste(lev, " (r=", correl, ")", sep = "")
    
    # Graphs with "xyplot()":

    pages <- (as.numeric(xy$variable) - 1) %/% prod(layout) + 1
    npages <- max(pages)

    for (n.p in 1:npages)
      {if (.Device %in% c("null device", "X11", "windows", "quartz"))
          {dev.new(title = colnames(y)[j])
          }
      trellis.par.set(list(fontsize = list(text=font.size)))
      plot(xyplot(x ~ y|variable, data = xy[pages == n.p, ], xlab = colnames(y)[j], ylab = "",
          main = "", layout = layout, as.table = TRUE, type = "n",
          panel = function(x, y, subscripts, ...) {
          panel.xyplot(x, y, ...); ltext(x, y, xy$nom[subscripts], cex = 0.7)
          }, scales = list(y = list(relation = "free"), tick.number=4)),
          packet.panel = my.packet.panel)
      }
    
    levels(xy$variable) <- lev
  }

  return(invisible(NULL))
}
