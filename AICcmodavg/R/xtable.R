##methods for xtable

##aictab
xtable.aictab <- function(x, caption = NULL, label = NULL, align = NULL,
                          digits = NULL, display = NULL, nice.names = TRUE,
                          include.AICc = TRUE, include.LL = TRUE,
                          include.Cum.Wt = FALSE, ...) {

  ##change to nicer names
  if(nice.names) {
    new.delta <- names(x)[4]
    new.weight <- names(x)[6]
    names(x)[1] <- "Model"
    names(x)[2] <- "K"
    #names(x)[4] <- paste("$\\delta$", unlist(strsplit(new.delta, "_"))[2], collapse = " ") #requires sanitize.text.function( )
    names(x)[4] <- paste(unlist(strsplit(new.delta, "_")), collapse = " ")
    names(x)[6] <- paste(unlist(strsplit(new.weight, "Wt")), "weight", collapse = " ")
    names(x)[7] <- "log-Likelihood"
    names(x)[8] <- "Cumulative weight"
  }

  #format to data.frame
  x <- data.frame(x, check.names = FALSE)
  class(x) <- c("xtable","data.frame")
  
  ##with AICc and LL but not Cum.Wt
  if(include.AICc && include.LL && !include.Cum.Wt) {
    x <- x[, c(1:4, 6:7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f","f"))
  }

  ##with AICc, but not LL and Cum.Wt
  if(include.AICc && !include.LL && !include.Cum.Wt) {
    x <- x[, c(1:4, 6)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  ##without AICc, but with LL but not Cum.Wt
  if(!include.AICc && include.LL && !include.Cum.Wt) {
    x <- x[, c(1:2, 4, 6:7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  ##without AICc and LL and Cum.Wt
  if(!include.AICc && !include.LL && !include.Cum.Wt) {
    x <- x[, c(1:2, 4, 6)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f"))
  }

  ##with AICc and LL and Cum.Wt
  if(include.AICc && include.LL && include.Cum.Wt) {
    x <- x[, c(1:4, 6:8)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f","f","f"))
  }

  ##with AICc, but not LL but with Cum.Wt
  if(include.AICc && !include.LL && include.Cum.Wt) {
    x <- x[, c(1:4, 6, 8)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f","f"))
  }

  ##without AICc, but with LL and Cum.Wt
  if(!include.AICc && include.LL && include.Cum.Wt) {
    x <- x[, c(1:2, 4, 6:8)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f","f"))
  }

  ##without AICc and LL but with Cum.Wt
  if(!include.AICc && !include.LL && include.Cum.Wt) {
    x <- x[, c(1:2, 4, 6, 8)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  caption(x) <- caption
  label(x) <- label
  return(x)
}



##modavg
xtable.modavg <- function(x, caption = NULL, label = NULL, align = NULL,
                          digits = NULL, display = NULL,
                          nice.names = TRUE, print.table = FALSE, ...) {

  if(print.table) {
    ##extract model selection table
    modavg.table <- data.frame(x$Mod.avg.table[, c(1:4, 6, 8:9)], check.names = FALSE)

    ##change to nicer names
    if(nice.names) {
      new.delta <- names(modavg.table)[4]
      new.weight <- names(modavg.table)[5]
      names(modavg.table)[1] <- "Model"
      names(modavg.table)[2] <- "K"
      ##names(x)[4] <- paste("$\\delta$", unlist(strsplit(new.delta, "_"))[2], collapse = " ") #requires sanitize.text.function( )
      names(modavg.table)[4] <- paste(unlist(strsplit(new.delta, "_")), collapse = " ")
      names(modavg.table)[5] <- paste(unlist(strsplit(new.weight, "Wt")), "weight", collapse = " ")
      names(modavg.table)[6] <- paste("Beta(", x$Parameter, ")", sep = "")
      names(modavg.table)[7] <- paste("SE(", x$Parameter, ")", sep = "")
    }

    ##format to data.frame
    class(modavg.table) <- c("xtable","data.frame")

    align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r","r"))
    digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2,2,2,2))
    display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f","f","f","f"))
  }  
    

    ##print model-averaged estimate, unconditional SE, CI
    if(!print.table) {
      
      ##model-averaged estimate                          
      modavg.table <- data.frame(Mod.avg.beta = x$Mod.avg.beta, Uncond.SE = x$Uncond.SE,
                                 Lower.CL = x$Lower.CL, Upper.CL = x$Upper.CL, check.names = FALSE)
      rownames(modavg.table) <- x$Parameter
    
      ##change to nicer names
      if(nice.names) {
        names(modavg.table)[1] <- "Model-averaged beta estimate" 
        names(modavg.table)[2] <- "Unconditional SE"
        names(modavg.table)[3] <- paste(100*x$Conf.level, "%", " lower limit", sep = "")
        names(modavg.table)[4] <- paste(100*x$Conf.level, "%", " upper limit", sep = "")
      }
  
      ##format to data.frame
      class(modavg.table) <- c("xtable","data.frame")

      align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
      digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2))
      display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f"))
    }

  caption(modavg.table) <- caption
  label(modavg.table) <- label
  return(modavg.table)
}



##modavgShrink
xtable.modavgShrink <- function(x, caption = NULL, label = NULL, align = NULL,
                                digits = NULL, display = NULL, 
                                nice.names = TRUE, print.table = FALSE, ...) {

  if(print.table) {
    ##extract model selection table
    modavg.table <- data.frame(x$Mod.avg.table[, c(1:4, 6, 8:9)], check.names = FALSE)

    ##change to nicer names
    if(nice.names) {
      new.delta <- names(modavg.table)[4]
      new.weight <- names(modavg.table)[5]
      names(modavg.table)[1] <- "Model"
      names(modavg.table)[2] <- "K"
      ##names(x)[4] <- paste("$\\delta$", unlist(strsplit(new.delta, "_"))[2], collapse = " ") #requires sanitize.text.function( )
      names(modavg.table)[4] <- paste(unlist(strsplit(new.delta, "_")), collapse = " ")
      names(modavg.table)[5] <- paste(unlist(strsplit(new.weight, "Wt")), "weight", collapse = " ")
      names(modavg.table)[6] <- paste("Beta(", x$Parameter, ")", sep = "")
      names(modavg.table)[7] <- paste("SE(", x$Parameter, ")", sep = "")
    }

    ##format to data.frame
    class(modavg.table) <- c("xtable","data.frame")

    align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r","r"))
    digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2,2,2,2))
    display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f","f","f","f"))
  }  


  ##print model-averaged estimate, unconditional SE, CI
  if(!print.table) {
      
    ##model-averaged estimate                          
    modavg.table <- data.frame(Mod.avg.beta = x$Mod.avg.beta, Uncond.SE = x$Uncond.SE,
                               Lower.CL = x$Lower.CL, Upper.CL = x$Upper.CL, check.names = FALSE)
    rownames(modavg.table) <- x$Parameter
    
    ##change to nicer names
    if(nice.names) {
      names(modavg.table)[1] <- "Model-averaged beta estimate" 
      names(modavg.table)[2] <- "Unconditional SE"
      names(modavg.table)[3] <- paste(100*x$Conf.level, "%", " lower limit", sep = "")
      names(modavg.table)[4] <- paste(100*x$Conf.level, "%", " upper limit", sep = "")
    }
  
    ##format to data.frame
    class(modavg.table) <- c("xtable","data.frame")
    
    align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
    digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2))
    display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f"))
  }

  caption(modavg.table) <- caption
  label(modavg.table) <- label
  return(modavg.table)
}



##modavgPred
xtable.modavgPred <- function(x, caption = NULL, label = NULL, align = NULL,
                              digits = NULL, display = NULL, 
                              nice.names = TRUE, ...) {
  modavg.pred <- data.frame(x$matrix.output, check.names = FALSE)
  
  ##change to nicer names
  if(nice.names) {
    names(modavg.pred)[1] <- "Model-averaged predictions" 
    names(modavg.pred)[2] <- "Unconditional SE"
    names(modavg.pred)[3] <- paste(100*x$conf.level, "%", " lower limit", sep = "")
    names(modavg.pred)[4] <- paste(100*x$conf.level, "%", " upper limit", sep = "")
  }
  
  ##format to data.frame
  class(modavg.pred) <- c("xtable","data.frame")

  align(modavg.pred) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
  digits(modavg.pred) <- switch(1+is.null(digits), digits, c(0,2,2,2,2))
  display(modavg.pred) <- switch(1+is.null(display), display, c("s","f","f","f","f"))

  caption(modavg.pred) <- caption
  label(modavg.pred) <- label
  return(modavg.pred)
}



##dictab
xtable.dictab <- function(x, caption = NULL, label = NULL, align = NULL,
                          digits = NULL, display = NULL, nice.names = TRUE,
                          include.DIC = TRUE, include.Cum.Wt = FALSE, ...) {

  ##change to nicer names
  if(nice.names) {
    names(x)[1] <- "Model"
    names(x)[2] <- "pD"
    names(x)[3] <- "DIC"
    names(x)[4] <- "Delta DIC"
    names(x)[6] <- "DIC weight"
    names(x)[7] <- "Cumulative weight"
  }

  #format to data.frame
  x <- data.frame(x, check.names = FALSE)
  class(x) <- c("xtable","data.frame")
  
  ##with DIC and Cum.Wt
  if(include.DIC && include.Cum.Wt) {
    x <- x[, c(1:4, 6:7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f","f"))
  }

  ##with DIC but not Cum.Wt
  if(include.DIC && !include.Cum.Wt) {
    x <- x[, c(1:4, 6)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  ##without DIC but with Cum.Wt
  if(!include.DIC && include.Cum.Wt) {
    x <- x[, c(1:2, 4, 6:7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  ##without DIC and Cum.Wt
  if(!include.DIC && !include.Cum.Wt) {
    x <- x[, c(1:2, 4, 6)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f"))
  }

  caption(x) <- caption
  label(x) <- label
  return(x)
}



##modavgEffect
xtable.modavgEffect <- function(x, caption = NULL, label = NULL, align = NULL,
                                digits = NULL, display = NULL, 
                                nice.names = TRUE, print.table = FALSE, ...) {

  if(print.table) {
    ##extract model selection table
    modavg.table <- data.frame(x$Mod.avg.table[, c(1:4, 6, 8:9)], check.names = FALSE)

    ##change to nicer names
    if(nice.names) {
      new.delta <- names(modavg.table)[4]
      new.weight <- names(modavg.table)[5]
      names(modavg.table)[1] <- "Model"
      names(modavg.table)[2] <- "K"
      ##names(x)[4] <- paste("$\\delta$", unlist(strsplit(new.delta, "_"))[2], collapse = " ") #requires sanitize.text.function( )
      names(modavg.table)[4] <- paste(unlist(strsplit(new.delta, "_")), collapse = " ")
      names(modavg.table)[5] <- paste(unlist(strsplit(new.weight, "Wt")), "weight", collapse = " ")
      names(modavg.table)[6] <- paste("Effect(", x$Group1, " - ", x$Group2, ")", sep = "")
      names(modavg.table)[7] <- paste("SE(", x$Group1, " - ", x$Group2, ")", sep = "")
    }

    ##format to data.frame
    class(modavg.table) <- c("xtable","data.frame")

    align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r","r"))
    digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2,2,2,2))
    display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f","f","f","f"))
  }  


  ##print model-averaged estimate, unconditional SE, CI
  if(!print.table) {
      
    modavg.table <- data.frame(Mod.avg.beta = x$Mod.avg.eff, Uncond.SE = x$Uncond.se,
                               Lower.CL = x$Lower.CL, Upper.CL = x$Upper.CL, check.names = FALSE)
    rownames(modavg.table) <- x$Group.variable

    ##change to nicer names
    if(nice.names) {
      names(modavg.table)[1] <- "Model-averaged effect size" 
      names(modavg.table)[2] <- "Unconditional SE"
      names(modavg.table)[3] <- paste(100*x$Conf.level, "%", " lower limit", sep = "")
      names(modavg.table)[4] <- paste(100*x$Conf.level, "%", " upper limit", sep = "")
    }
  
    ##format to data.frame
    class(modavg.table) <- c("xtable","data.frame")

    align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
    digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2))
    display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f"))
  }

  caption(modavg.table) <- caption
  label(modavg.table) <- label
  return(modavg.table)
}



##multComp
xtable.multComp <- function(x, caption = NULL, label = NULL, align = NULL,
                            digits = NULL, display = NULL, 
                            nice.names = TRUE, print.table = FALSE, ...) {

  if(print.table) {
    ##extract model selection table
    modavg.table <- data.frame(x$model.table[, c(1:4, 6)], check.names = FALSE)

    ##change to nicer names
    if(nice.names) {
      new.delta <- names(modavg.table)[4]
      new.weight <- names(modavg.table)[5]
      names(modavg.table)[1] <- "Group structure"
      names(modavg.table)[2] <- "K"
      ##names(x)[4] <- paste("$\\delta$", unlist(strsplit(new.delta, "_"))[2], collapse = " ") #requires sanitize.text.function( )
      names(modavg.table)[4] <- paste(unlist(strsplit(new.delta, "_")), collapse = " ")
      names(modavg.table)[5] <- paste(unlist(strsplit(new.weight, "Wt")), "weight", collapse = " ")
    }

    ##format to data.frame
    class(modavg.table) <- c("xtable","data.frame")

    align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2,2))
    display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f","f"))
  }  

  
  ##print model-averaged estimate, unconditional SE, CI
  if(!print.table) {
    modavg.table <- data.frame(Group = x$ordered.levels,
                               Model.avg.est = x$model.avg.est[, "Mod.avg.est"],
                               Uncond.SE = x$model.avg.est[, "Uncond.SE"],
                               Lower.CL = x$model.avg.est[, "Lower.CL"],
                               Upper.CL = x$model.avg.est[, "Upper.CL"], check.names = FALSE)
  
    ##change to nicer names
    if(nice.names) {
      names(modavg.table)[2] <- paste("Model-averaged estimates (", x$factor.id, ")", sep = "") 
      names(modavg.table)[3] <- "Unconditional SE"
      names(modavg.table)[4] <- paste(100*x$conf.level, "%", " lower limit", sep = "")
      names(modavg.table)[5] <- paste(100*x$conf.level, "%", " upper limit", sep = "")
    }
  
    ##format to data.frame
    class(modavg.table) <- c("xtable","data.frame")

    align(modavg.table) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(modavg.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2,2))
    display(modavg.table) <- switch(1+is.null(display), display, c("s","f","f","f","f","f"))

  }

  caption(modavg.table) <- caption
  label(modavg.table) <- label
  return(modavg.table)
}



##boot.wt - class aictab - potentially create new class for boot.wt
xtable.boot.wt <- function(x, caption = NULL, label = NULL, align = NULL,
                           digits = NULL, display = NULL, nice.names = TRUE,
                           include.AICc = TRUE, include.AICcWt = FALSE, ...) {

  ##change to nicer names
  if(nice.names) {
    new.delta <- names(x)[4]
    new.weight <- names(x)[6]
    names(x)[1] <- "Model"
    names(x)[2] <- "K"
    names(x)[4] <- paste(unlist(strsplit(new.delta, "_")), collapse = " ")
    names(x)[6] <- paste(unlist(strsplit(new.weight, "Wt")), "weight", collapse = " ")
    names(x)[7] <- "Pi weight"
  }

  #format to data.frame
  x <- data.frame(x, check.names = FALSE)
  class(x) <- c("xtable","data.frame")
  
  ##with AICc but not AICc.Wt
  if(include.AICc && !include.AICcWt) {
    x <- x[, c(1:4, 7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  ##with AICc and AICc.Wt
  if(include.AICc && include.AICcWt) {
    x <- x[, c(1:4, 6:7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f","f"))
  }

  ##without AICc, but with AICc.Wt
  if(!include.AICc && include.AICcWt) {
    x <- x[, c(1:2, 4, 6:7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  ##without AICc and AICc.Wt
  if(!include.AICc && !include.AICcWt) {
    x <- x[, c(1:2, 4, 7)]
    align(x) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
    digits(x) <- switch(1+is.null(digits), digits, c(0,0,2,2,2))
    display(x) <- switch(1+is.null(display), display, c("s","s","f","f","f"))
  }

  caption(x) <- caption
  label(x) <- label
  return(x)
}



##mb.chisq - for single-season model
xtable.mb.chisq <- function(x, caption = NULL, label = NULL, align = NULL,
                            digits = NULL, display = NULL, nice.names = TRUE,
                            include.detection.histories = TRUE, ...) {


  ##stop if dynamic occupancy model
  if(identical(x$model.type, "dynamic")) {
    stop("\n'xtable' does not yet support dynamic occupancy models\n")
  }

  ##extract names
  x.names <- names(x)

  ##stop if no chi-square table is present
  if(!any(x.names == "chisq.table")) {
    stop("\nchi-square table must be included in the object\n")
  }

  
  ##extract table
  chisq.table <- data.frame("Detection.history" = rownames(x$chisq.table),
                            x$chisq.table, check.names = FALSE)

    
  ##change to nicer names
  if(nice.names) {
    colnames(chisq.table)[1] <- "Detection history"
    ##extract rownames
    rows <- rownames(x$chisq.table)
    ##replace NA by "-", here " -" to avoid creating endash in LaTeX
    ##also possible to use {} between consecutive dashes - this requires sanitization in print( )
    new.rows <- gsub(pattern = "NA", replacement = " -", rows)
    chisq.table[, "Detection history"] <- new.rows
    rownames(chisq.table) <- new.rows
  }

  ##format to data.frame
  class(chisq.table) <- c("xtable", "data.frame")

  ##do not include detection history as a column (only rownames)
  if(!include.detection.histories) {

    ##exclude column with detection histories
    chisq.table <- chisq.table[, 2:5]
    
    align(chisq.table) <- switch(1+is.null(align), align, c("l","r","r","r","r"))
    digits(chisq.table) <- switch(1+is.null(digits), digits, c(0,2,2,2,2))
    display(chisq.table) <- switch(1+is.null(display), display, c("s","f","f","f","f"))
  }

  ##include detection history as a column
  if(include.detection.histories) {

    ##exclude column with detection histories
    chisq.table <- chisq.table[, 1:5]
    
    align(chisq.table) <- switch(1+is.null(align), align, c("l","r","r","r","r","r"))
    digits(chisq.table) <- switch(1+is.null(digits), digits, c(0,0,2,2,2,2))
    display(chisq.table) <- switch(1+is.null(display), display, c("s","s","f","f","f","f"))
  }

  caption(chisq.table) <- caption
  label(chisq.table) <- label
  return(chisq.table)
}



##add method for detHist class
xtable.detHist <- function(x, caption = NULL, label = NULL, align = NULL,
                           digits = NULL, display = NULL, 
                           nice.names = TRUE, table.detHist = "freq",
                           ...) {
  ##for single season, display frequencies in a single matrix
  if(x$n.seasons == 1) {

    ##display detection histories
    if(identical(table.detHist, "hist")){
      det.hist <- x$hist.table.full
      det.mat <- matrix(det.hist, nrow = 1)
      
      if(nice.names) {
        det.names <- names(det.hist)
        new.names <- gsub(pattern = "NA", replacement = " -", det.names)
        colnames(det.mat) <- new.names
        rownames(det.mat) <- "Season-1"
      }
      
      det.frame <- as.data.frame(det.mat)
      n.cols <- ncol(det.frame)
    }

    ##display frequencies
    if(identical(table.detHist, "freq")) {
        det.frame <- as.data.frame(x$out.freqs)
        n.cols <- ncol(det.frame)
      }

    ##display proportions
    if(identical(table.detHist, "prop")){
      det.frame <- as.data.frame(x$out.props)
      n.cols <- ncol(det.frame)
    }
       
  } else {

      ##display entire detection histories
      if(identical(table.detHist, "hist")) {
        det.hist <- x$hist.table.full
        det.mat <- matrix(det.hist, nrow = 1)
      
        if(nice.names) {
          det.names <- names(det.hist)
          new.names <- gsub(pattern = "NA", replacement = " -", det.names)
          colnames(det.mat) <- new.names
          rownames(det.mat) <- "All seasons"
        }
        
        det.frame <- as.data.frame(det.mat)
        n.cols <- ncol(det.frame)
      }

      ##display frequencies
      if(identical(table.detHist, "freq")) {
        det.frame <- as.data.frame(x$out.freqs)
        #det.frame[1, 3:6] <- "."
        n.cols <- ncol(det.frame)
      }

      ##display proportions
      if(identical(table.detHist, "prop")) {
        det.frame <- as.data.frame(x$out.props)
#        det.frame[1, 2:4] <- "."
        n.cols <- ncol(det.frame)
      }
    }

    
  ##format to data.frame
  class(det.frame) <- c("xtable","data.frame")
  
  align(det.frame) <- switch(1+is.null(align), align, c("l", rep("r", n.cols)))

  if(identical(table.detHist, "prop")) {
    digits(det.frame) <- switch(1+is.null(digits), digits, c(0, rep(2, n.cols)))
  } else {
    digits(det.frame) <- switch(1+is.null(digits), digits, c(0, rep(0, n.cols)))
  }
  display(det.frame) <- switch(1+is.null(display), display, c("s", rep("f", n.cols)))
  
  caption(det.frame) <- caption
  label(det.frame) <- label
  return(det.frame)
}



##add method for countHist class
xtable.countHist <- function(x, caption = NULL, label = NULL, align = NULL,
                             digits = NULL, display = NULL, 
                             nice.names = TRUE, table.countHist = "count",
                             ...) {
  ##for single season, display frequencies in a single matrix
  if(x$n.seasons == 1) {

    ##display count histories
    if(identical(table.countHist, "hist")){
      det.hist <- x$hist.table.full
      det.mat <- matrix(det.hist, nrow = 1)
      
      if(nice.names) {
        det.names <- names(det.hist)
        new.names <- gsub(pattern = "NA", replacement = " -", det.names)
        colnames(det.mat) <- new.names
        rownames(det.mat) <- "Season-1"
      }
      
      det.frame <- as.data.frame(det.mat)
      n.cols <- ncol(det.frame)
    }

    ##display frequencies
    if(identical(table.countHist, "freq")) {
        det.frame <- as.data.frame(x$out.freqs)
        n.cols <- ncol(det.frame)
      }

    ##display counts
    if(identical(table.countHist, "count")) {
        det.mat <- matrix(x$count.table.full, nrow = 1)
        colnames(det.mat) <- names(x$count.table.full)

        if(nice.names) {
          rownames(det.mat) <- "Season-1"
        }
        
        n.cols <- ncol(det.mat)
        det.frame <- as.data.frame(det.mat)
      }

    ##display proportions
    if(identical(table.countHist, "prop")){
      det.frame <- as.data.frame(x$out.props)
      n.cols <- ncol(det.frame)
    }
       
  } else {

      ##display entire detection histories
      if(identical(table.countHist, "hist")) {
        det.hist <- x$hist.table.full
        det.mat <- matrix(det.hist, nrow = 1)
      
        if(nice.names) {
          det.names <- names(det.hist)
          new.names <- gsub(pattern = "NA", replacement = " -", det.names)
          colnames(det.mat) <- new.names
          rownames(det.mat) <- "All seasons"
        }
        
        det.frame <- as.data.frame(det.mat)
        n.cols <- ncol(det.frame)
      }

      ##display frequencies
      if(identical(table.countHist, "freq")) {
        det.frame <- as.data.frame(x$out.freqs)
        #det.frame[1, 3:6] <- "."
        n.cols <- ncol(det.frame)
      }

      ##display counts
      if(identical(table.countHist, "count")) {
        det.mat <- matrix(x$count.table.full, nrow = 1)
        colnames(det.mat) <- names(x$count.table.full)

        if(nice.names) {
          rownames(det.mat) <- "All seasons"
        }
        n.cols <- ncol(det.mat)
        det.frame <- as.data.frame(det.mat)
      }

      ##display proportions
      if(identical(table.countHist, "prop")) {
        det.frame <- as.data.frame(x$out.props)
#        det.frame[1, 2:4] <- "."
        n.cols <- ncol(det.frame)
      }

    }

    
  ##format to data.frame
  class(det.frame) <- c("xtable","data.frame")
  
  align(det.frame) <- switch(1+is.null(align), align, c("l", rep("r", n.cols)))

  if(identical(table.countHist, "prop")) {
    digits(det.frame) <- switch(1+is.null(digits), digits, c(0, rep(2, n.cols)))
  } else {
    digits(det.frame) <- switch(1+is.null(digits), digits, c(0, rep(0, n.cols)))
  }
  display(det.frame) <- switch(1+is.null(display), display, c("s", rep("f", n.cols)))
  
  caption(det.frame) <- caption
  label(det.frame) <- label
  return(det.frame)
}



##add method for countDist class
xtable.countDist <- function(x, caption = NULL, label = NULL, align = NULL,
                             digits = NULL, display = NULL, 
                             nice.names = TRUE, table.countDist = "distance",
                             ...) {
  ##for single season, display frequencies in a single matrix
  if(x$n.seasons == 1) {

    ##display counts across distance classes
    if(identical(table.countDist, "distance")){
      det.dist <- matrix(x$dist.sums.full, nrow = 1)
      colnames(det.dist) <- names(x$dist.sums.full)

      if(nice.names) {
        rownames(det.dist) <- "Season-1"
      }
      
      n.cols <- ncol(det.dist)
      det.frame <- as.data.frame(det.dist)
    }

    ##display frequencies
    if(identical(table.countDist, "freq")) {
        det.frame <- as.data.frame(x$out.freqs)
        n.cols <- ncol(det.frame)
      }

    ##display counts
    if(identical(table.countDist, "count")) {
        det.mat <- matrix(x$count.table.full, nrow = 1)
        colnames(det.mat) <- names(x$count.table.full)

        if(nice.names) {
          rownames(det.mat) <- "Season-1"
        }
        
        n.cols <- ncol(det.mat)
        det.frame <- as.data.frame(det.mat)
      }

    ##display proportions
    if(identical(table.countDist, "prop")){
      det.frame <- as.data.frame(x$out.props)
      n.cols <- ncol(det.frame)
    }
       
  } else {

      ##display counts across distance classes
      if(identical(table.countDist, "distance")) {
        ##assumes the same distance classes were used each year
        det.dist <- matrix(unlist(x$dist.sums.seasons),
                           nrow = x$n.seasons)
        colnames(det.dist) <- names(x$dist.sums.full)

        if(nice.names) {
          rownames(det.dist) <- paste("Season-", 1:x$n.seasons, sep = "")
        }
      
        n.cols <- ncol(det.dist)
        det.frame <- as.data.frame(det.dist)
      }

      ##display frequencies
      if(identical(table.countDist, "freq")) {
        det.frame <- as.data.frame(x$out.freqs)
        #det.frame[1, 3:6] <- "."
        n.cols <- ncol(det.frame)
      }

      ##display counts
      if(identical(table.countDist, "count")) {
        det.mat <- matrix(x$count.table.full, nrow = 1)
        colnames(det.mat) <- names(x$count.table.full)

        if(nice.names) {
          rownames(det.mat) <- "All seasons"
        }
        
        n.cols <- ncol(det.mat)
        det.frame <- as.data.frame(det.mat)
      }

      ##display proportions
      if(identical(table.countDist, "prop")) {
        det.frame <- as.data.frame(x$out.props)
#        det.frame[1, 2:4] <- "."
        n.cols <- ncol(det.frame)
      }

    }

    
  ##format to data.frame
  class(det.frame) <- c("xtable","data.frame")
  
  align(det.frame) <- switch(1+is.null(align), align, c("l", rep("r", n.cols)))

  if(identical(table.countDist, "prop")) {
    digits(det.frame) <- switch(1+is.null(digits), digits, c(0, rep(2, n.cols)))
  } else {
    digits(det.frame) <- switch(1+is.null(digits), digits, c(0, rep(0, n.cols)))
  }
  display(det.frame) <- switch(1+is.null(display), display, c("s", rep("f", n.cols)))
  
  caption(det.frame) <- caption
  label(det.frame) <- label
  return(det.frame)
}



