nberShade.ggplot <- function(gg = ggplot2::last_plot(),
                             fill = "grey80",
                             color = NA,
                             size = 0.5, 
                             alpha = 0.5,
                             xrange = NULL,
                             openShade = TRUE, ...){
  Start <- End <- ymin <- ymax <-
    "Defined here to defeat 'no visible binding' R CMD check warnings"
  yrng <- extendrange(eval(gg$mapping$y, gg$data))
  nber.dates <- transform(data.frame(as.data.frame(nberDates()),
                                     ymin = yrng[1],
                                     ymax = yrng[2]), 
                          Start = as.Date(as.character(Start), format = "%Y%m%d"),
                          End = as.Date(as.character(End), format = "%Y%m%d"))
  
  if( is.null(xrange) ){
    if( openShade ){
      nr <- nrow(nber.dates)
      if(is.na(nber.dates[nr, "End"])) nber.dates[nr, "End"] <- Sys.Date()
    }  
  }
  else {
    xrange <- as.Date(xrange)
    nber.dates <- subset(nber.dates, (End > xrange[1]) | (is.na(End)))
    if(nber.dates[1, "Start"] < xrange[1]) 
      nber.dates[1, "Start"] <- xrange[1]    
    nber.dates <- subset(nber.dates, Start < xrange[2])
    shouldReplaceLast <- (nber.dates[dim(nber.dates)[1], "End"] > xrange[2])
    shouldReplaceLast <- ifelse(is.na(shouldReplaceLast), TRUE, shouldReplaceLast)
    if( openShade | shouldReplaceLast ){
      nber.dates[ dim(nber.dates)[1], "End"] <- xrange[2]    
    }
  }
  
  gg <- gg + ggplot2::geom_rect( ggplot2::aes(xmin = Start, xmax = End, x = NULL, y = NULL, 
                                                ymin = ymin, ymax = ymax), color = color, fill = fill, 
                                 size = 0.5, alpha = 0.5, data = nber.dates) + 
                                   ggplot2::scale_y_continuous(expand = c(0,0))
  
  if( ! openShade ){
    ## color = fill is the behavior in nberShade.default
    gg <- gg + ggplot2::geom_vline(xintercept = as.numeric(tail(nber.dates, 1)$Start), 
                                    color = fill, size = 1.5)
  }
  
  return(gg)
} 
