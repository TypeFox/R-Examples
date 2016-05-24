evidence <-
  function(aic.table, model.high = "top", model.low = "second.ranked") {
    ##if multComp object, extract relevant table
    if(identical(class(aic.table)[1], "multComp")) {
      if(!is.data.frame(aic.table)) {
        aic.table <- aic.table$model.table
      }  
      ##coerce to aictab
      class(aic.table) <- c("aictab", "data.frame")
    }

    if(identical(class(aic.table)[1], "boot.wt")) {
      ##coerce to aictab
      class(aic.table) <- c("aictab", "data.frame")
    }
    
    if(!identical(class(aic.table)[1], "aictab")) {stop("\nThe input object must be of class 'aictab'\n")}

    ##sort model table in case it is not
    sort.tab <- aic.table[order(aic.table[, 4]), ]

    ##top model
    if(identical(model.high, "top")) {
      ##determine which is the highest ranking model
      top.name <- sort.tab[1, 1]
      top.wt <- sort.tab[1, 6]
    } else {
      top.name <- model.high
      top.wt <- sort.tab[which(sort.tab$Modnames == paste(model.high)), 6]
    }

    ##model compared
    if(identical(model.low, "second.ranked")) {
      sec.name <- sort.tab[2, 1]
      sec.wt <- sort.tab[2, 6]
    } else {
      sec.name <- model.low
      sec.wt <- sort.tab[which(sort.tab$Modnames == paste(model.low)), 6]
    }

    ##compute evidence ratio
    ev.ratio <- top.wt/sec.wt
    ev.ratio.list <- list("Model.high" = paste(top.name), "Model.low" = paste(sec.name), "Ev.ratio" = ev.ratio)
    class(ev.ratio.list) <- c("evidence", "list")
    return(ev.ratio.list)
  }


print.evidence <- function(x, digits = 2, ...) {
  cat("\nEvidence ratio between models '", x$Model.high,"' and '", x$Model.low, "':\n", sep = "")
  cat(round(x$Ev.ratio, digits = digits), "\n\n")
}
