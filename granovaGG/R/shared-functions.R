LayoutFourPlotsPerPage <- function(list.of.plots) {
  four.plot.pages <- floor(length(list.of.plots) / 4)
  remainder.plots <- length(list.of.plots) %% 4 # a %% b is notation for a modulo b

  LayoutFourPlotPages(list.of.plots, four.plot.pages)
  LayoutRemainderPlots(list.of.plots, remainder.plots)
}

LayoutFourPlotPages <- function(plot.list, pages) {
  i <- 1
  while (i < 4 * pages) {
    args.list <- list(plot.list[[i]],
                      plot.list[[i+1]],
                      plot.list[[i+2]],
                      plot.list[[i+3]]
                 )
    args.list <- c(args.list, list(nrow = 2, ncol = 2))
    do.call(grid.arrange, args.list)
    DisplayEndOfPageMessage()
    i <- i + 4
  }
}

LayoutRemainderPlots <- function(plot.list, remainder.plots) {
  if (remainder.plots > 0) {
    remainder.start <- length(plot.list) - remainder.plots + 1
    remainder.end   <- length(plot.list)
    args.list       <- lapply(X   = remainder.start : remainder.end,
                              FUN = GetListElementByIndex,
                              x   = plot.list
                       )
             
    args.list <- c(args.list, list(nrow = 2, ncol = 2))
    do.call(grid.arrange, args.list)
  }
}

GetListElementByIndex <- function(x, index) {
  return(x[[index]])
}


DisplayEndOfPageMessage <- function() {
  if (interactive()) {
    readline("Examine the contrast plots and consider printing. When you're done, press <Return>")
  }
}

Theme <- function(plot.theme) {
    return(
      eval(call(plot.theme))
    )
}

GetContrastName <- function(contrast.data, index) {
  if (is.null(dimnames(contrast.data))) {
    contrast.name <- paste(index)
  }
  else {
    provided.contrast.name <- dimnames(contrast.data)[[2]][index]
    contrast.name <- paste(provided.contrast.name)
  }
  
  return(paste("Contrast ", contrast.name, sep=""))
}

# ReorderDataByColumn reorders a data structure according to the values in "column"
# For example, "sort the data on all these cars according to their highway mileage"
ReorderDataByColumn <- function(x, column) { 
  return(x[order(x[[column]]), ])
}

OverlapWarning <- function(data, tolerance) {
  first <- FirstElementOverlaps(data, tolerance)
  inner <- InnerElementsOverlap(data, tolerance)
  last  <- LastElementOverlaps(data, tolerance)
  
  return(c(first, inner, last))
}

FirstElementOverlaps <- function(data, tolerance) {
  return(
    abs(data[1] - data[2]) < tolerance
  )
}

InnerElementsOverlap <- function(data, tolerance) {
  return(sapply(
           X         = 2:(length(data) - 1),
           FUN       = NeighborOverlaps,
           data      = data,
           tolerance = tolerance
         )
  )
}

LastElementOverlaps <- function(data, tolerance) {
  last.index <- length(data)
  return(
    abs(data[last.index] - data[(last.index) - 1]) < tolerance
  )
}

NeighborOverlaps <- function(data, tolerance, index) {
  left.overlap  <- abs(data[index] - data[(index - 1)]) < tolerance
  right.overlap <- abs(data[index] - data[(index + 1)]) < tolerance
  
  return(left.overlap | right.overlap)
}
