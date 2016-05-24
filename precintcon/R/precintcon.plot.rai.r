#' @export
precintcon.plot.rai <- function(
  ..., 
  granularity     = "m",
  xlab            = "Month",
  ylab            = "RAI", 
  ylim            = c(-3,3),
  legend          = NULL,
  fontsize        = 10, 
  axis.text.color = "black", 
  export          = FALSE, 
  export.name     = "rai_plot.png", 
  width           = 8.6, 
  height          = 7.5, 
  units           = "cm",
  args            = NA
) {
  
  l <- list(...)
  
  if (length(l) > 0) {
    
    if (length(l) > 1)
      par(ask = T)
    
    l <- lapply(l, FUN = precintcon.rai.analysis, granularity = granularity)
    
    varl <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)
    
    if (!is.null(legend) && length(varl) != length(legend))
      stop(paste("legend should has length equals to the number of input data. legend parameter length", 
                 length(legend), ": number of input data", length(varl)))
    
    else if (!is.null(legend))
      varl <- as.list(legend)
    
    plotl <- mapply(p.plot.rai, l, varl, 
                    MoreArgs = list(g = granularity, 
                                    xlab = xlab, ylab = ylab, ylim = ylim, fontsize = fontsize, 
                                    axis.text.color = axis.text.color),
                    SIMPLIFY = FALSE)
    
    for (i in 1:length(plotl)) {
      if (!export){
        print(plotl[[i]])
      }else{
        ggsave(paste(varl[[i]], export.name, sep = "_"), plotl[[i]], width = width, height = height, units = units)
      }
    }
    
    par(ask = F)
    
  } else {
    stop("empty input data in precintcon.plot.spi function.")
  }
}

p.plot.rai <- function(d, n, g,
                       xlab = ifelse(g == "m", "Month", "Year"),
                       ylab = "RAI", 
                       ylim = c(-3,3),
                       fontsize = 10, 
                       axis.text.color = "black"
) {
  
  if (is.element("precintcon.rai", class(d))) {
    
    data <- NA
    
    if (g == "m" || g == "a") {
      
      if (g == "m") {
        d[which(nchar(d[,2]) < 2), 2] <- paste("0", d[which(nchar(d[,2]) < 2), 2], sep = "")
      }
      
      data <- data.frame(
        x = as.Date(
          paste(d[ ,1], 
                if(g == "m"){
                  d[ ,2]
                }else{
                  "01"
                }
                , "01", sep = "/"), "%Y/%m/%d"), 
        y =  d$rai,
        dataset = paste(n, sep="")
      )
      
    } else {
      stop("invalid granularity value. It should be either 'a' for annual or 'm' for monthly")
    }
    
    p <- ggplot(data, aes_string(x = "x", y = "y")) + geom_bar(stat = "identity", position = "identity") +  
      xlab(xlab) + ylab(ylab) + 
      scale_x_date(expand            = c(1/48, 1/48), 
                   limits            = as.Date(c(data[1, 1], tail(data$x, n = 1))), 
                   date_breaks       = ifelse(g == "m", "20 months", "years"), 
                   labels            = date_format(ifelse(g == "m", "%b %y", "%Y")), 
                   date_minor_breaks = ifelse(g == "m", "1 month", "1 year")) +
      scale_y_continuous(breaks = seq(round(min(data$y), digits = 2) - .5, round(max(data$y), digits = 2) + .5, by = 1), 
                         limits = c(round(min(data$y), digits = 2) - .5, round(max(data$y), digits = 2) + .5)) +
      theme(text               = element_text(size = fontsize), 
            axis.text          = element_text(color = axis.text.color),
            axis.text.x        = element_text(angle = 25),
            axis.title.x       = element_text(vjust = .1),
            panel.grid.minor.x = element_blank()) +
      facet_grid(. ~ dataset)
    
    return(p)
    
  } else {
    stop("invalid input data type. It should be of type precintcon.rai")		
  }
}