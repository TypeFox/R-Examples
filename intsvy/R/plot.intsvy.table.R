plot.intsvy.table <- function(x, se=FALSE, stacked=FALSE, ...) {
  # it is assumed that last three columns are:  "Freq"       "Percentage" "Std.err."  
  vars <- setdiff(colnames(x), c("Freq", "Percentage", "Std.err."))
  colnames(x)[ncol(x)] = "se"
  nvar <- length(vars)
  x$PercentageL <- x$Percentage - x$se
  x$PercentageH <- x$Percentage + x$se
  pl <- NA
  if (stacked) {
    if (nvar == 1) {
      pl <- ggplot(data=x, aes_string(y = "Percentage", x=vars[1])) + 
        geom_bar(stat="identity", position="fill") + theme_bw() + coord_flip() + 
        theme(legend.position="top")
    } else {
      if (nvar == 2) {
        pl <- ggplot(data=x, aes_string(y = "Percentage", x=vars[1], fill=vars[2])) + 
          geom_bar(stat="identity", position="fill") + theme_bw() + coord_flip() + 
          theme(legend.position="top")
      } else {
        pl <- ggplot(data=x, aes_string(y = "Percentage", x=vars[1], fill=vars[3])) + 
          geom_bar(stat="identity", position="stack") + theme_bw() + coord_flip() + 
          facet_wrap(as.formula(paste0("~", vars[2])))+ 
          theme(legend.position="top")
      }
    }
  } else {
    if (nvar == 1) {
      pl <- ggplot(data=x, aes_string(y = "Percentage", x=vars[1])) + 
        geom_bar(stat="identity", position=position_dodge(width = 0.9)) + 
        theme_bw() + coord_flip() + 
        theme(legend.position="top")
    } else {
      if (nvar == 2) {
        pl <- ggplot(data=x, aes_string(y = "Percentage", x=vars[1], fill=vars[2])) + 
          geom_bar(stat="identity", position=position_dodge(width = 0.9)) + 
          theme_bw() + coord_flip() + 
          facet_wrap(as.formula(paste0("~", vars[2])))+ 
          theme(legend.position="top")
      } else {
        pl <- ggplot(data=x, aes_string(y = "Percentage", x=vars[1], fill=vars[2], col= vars[2])) + 
          geom_bar(stat="identity", position=position_dodge(width = 0.9)) + 
          theme_bw() + coord_flip() + 
          facet_wrap(as.formula(paste0("~", vars[3])))+ 
          theme(legend.position="top")
      }
    }
    
    if (se) {
      pl <- pl + geom_errorbar(aes_string(ymin="PercentageL", ymax="PercentageH"), color="black", 
                               position=position_dodge(width = 0.9), width=.35) 
    } 
  }
  pl
}

