#' Add cost benefit ratio axis to a decision curve plot.
#'
#'
#' @param xlim range of x-axis.
#' @param cost.benefits Character vector of the form c("c1:b1", "c2:b2", ..., "cn:bn") with integers ci, bi corresponding to specific cost:benefit ratios to print.
#' @param n.cost.benefits number of cost:benefit ratios to print if cost.benefit.axis = TRUE (default n.cost.benefit = 6).
#' @param line x-axis line to print the axis (default is 4).
#' @param ... other options sent to `axis`.
#'
#' @return List with components threshold, value and name.
#'
#' @seealso \code{\link{summary.decision_curve}},  \code{\link{decision_curve}}
#'
#' @export


Add_CostBenefit_Axis <- function(xlim,
                                 cost.benefits,
                                 n.cost.benefits = 6,
                                 line = 4, ...){

  if(missing(cost.benefits)){
    #no cost.benefits are given, so we must choose them (wisely) from available defaults
    #listed here
    tmp <- c( 1/2, 1/3, 1/4, 2/3, 3/4, 1/5, 1/10, 2/5, 3/5, 4/5, 7/10, 9/10, 5/6, 1/15, 1/20,  1/25, 1/30, 1/40, 1/50, 1/60, 1/75, 1/80, 1/100)
    dd <-  data.frame( "name" = c("1:1", "1:2", "1:3", "1:4", "2:3", "3:4", "1:5", "1:10", "2:5", "3:5", "4:5", "7:10", "9:10", "5:6", "1:15", "1:20", "1:25", "1:30", "1:40", "1:50", "1:60", "1:75","1:80",  "1:100",
                                  "2:1", "3:1", "4:1", "3:2", "4:3", "5:1", "10:1", "5:2", "5:3", "5:4", "10:7", "10:9", "6:5", "15:1", "20:1", "25:1", "30:1", "40:1", "50:1",  "60:1","75:1", "80:1",  "100:1"),
                       "value" = c(1/1, tmp, 1/tmp))
    dd$threshold <- dd$value/(1+dd$value)
    dd <- dd[order(dd$threshold),]
    threshold <- NULL
    dd.sub <- subset(dd, threshold >= min(xlim) & threshold <= max(xlim))
    #check to make sure there are default choices within xlim
    if(nrow(dd.sub)==0){
      warning("Range of risk thresholds printed do not span available default Cost:Benefit ratios. Cost:Benefit axis will not be printed.")
    }else{

      #evenly spaced increments of threshold from xlim
      thresh.quants <- quantile(seq(min(xlim), max(xlim), length.out = 100),
                                probs = seq(0, 1, length.out = n.cost.benefits),
                                type = 1)

      #now find C:B that are near thresh.quants
      match.quant <- numeric(length(thresh.quants))
      for(i in seq_along(thresh.quants)){
        match.quant[i] <- which.min(abs(thresh.quants[i] - dd.sub$threshold))
      }

      #just print the C:B that have been chosen
      out <- dd.sub[match.quant, ]

      axis(1, at   = out$threshold, labels = out$name, line = line, ...)
      invisible(out)
    }
  }else{
    #cost benefits are given, so we plot them.
    dd <- data.frame("name" = cost.benefits)

    #extract data from the
    CB.split <- strsplit(as.character(cost.benefits), split = ":")
    C <- sapply(CB.split, function(x) as.numeric(x[1]))
    B <- sapply(CB.split, function(x) as.numeric(x[2]))

    dd$value = C/B
    dd$threshold <- dd$value/(1+dd$value)
    dd <- dd[order(dd$threshold),]
    dd.sub <- subset(dd, threshold >= min(xlim) & threshold <= max(xlim))
    #check to make sure there are default choices within xlim
    if(nrow(dd.sub)==0){
      warning("Range of risk thresholds printed do not span provided Cost:Benefit ratios. Cost:Benefit axis will not be printed.")
    }else{

      out <- dd.sub
      axis(1, at   = out$threshold, labels = out$name, line = line, ...)
      invisible(out)
    }
  }
}


