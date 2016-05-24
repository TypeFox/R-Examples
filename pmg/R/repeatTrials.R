## Function and GUI suggested by Daniel Kaple
## Daniel Kaplan 

## make these pmg functions. resample? too generic?
pmgRepeatTrials = function(expr, n = 10) {
  ## return n values
  if(!is.expression(expr))
    expr = as.expression(substitute(expr))
  out = try(sapply(1:n, function(...) eval(expr, envir=.GlobalEnv)),silent=TRUE)
  if(inherits(out,"try-error")) {
    cat(gettext("Error occurred."))
    return(NA)
  }
    
  if(is.matrix(out))
    return(t(out))
  else
    return(out)
}



repeatTrialsGUI = function(container = NULL) {
  
  ## main dialog
  actions = c("","print", "plot","hist","boxplot")
  
  g = gframe("Repeat trials", container=container)
  tbl <- glayout(container = g)
  
  tbl[1,1, anchor=c(1,0)] <- "expression"
  tbl[1,2] <- (.expr <- gedit("", container=tbl))
  tbl[2,1, anchor=c(1,0)] <- "No. of times"
  tbl[2,2] <- (.n <- gedit(10, container=tbl, coerce.with=as.numeric))
  tbl[3,1:2] <- gseparator(container=tbl)
  tbl[4,1, anchor=c(1,0)] <- "Quick action"
  tbl[4,2, anchor=c(-1,0)] <- (.actions <- gdroplist(actions, editable=TRUE, container=tbl))
  tbl[5,1, anchor=c(1,0)] <- "Save as"
  tbl[5,2] <- (.saveAs <- gedit("", container=tbl))
  tbl[6,1:2] <- gseparator(container=tbl)

  
  tbl[7,2, anchor=c(-1,0)] <- (bg  <- ggroup(container=tbl))

  .resample <- gbutton("Repeat trials", container = bg)
  addSpace(bg, 10)
  .help <- gbutton("help", container= bg, handler = function(h,...) {
    ghelp("repeatTrials",package="pmg",
          container = pmgWC$new("repeatTrials", width=400, height=250))
  })
  
  visible(tbl) <- TRUE
  
  ## add handler to resample button
  addHandlerClicked(.resample, handler = function(h,...) {
    theExpression = svalue(.expr)
    if(theExpression == "") {
      cat("Need an expression\n")
      return()
    }
    theExpression = try(parse(text = theExpression),silent=TRUE)
    if(inherits(theExpression,"try-error")) {
      cat(gettext("Error in the expression\n"))
      return(FALSE)
    }
    
    n <- svalue(.n)
    if(is.na(n) || (!is.numeric(n) && n < 0)) {
      cat("No. of times is a positive intefer\n")
      return()
    }
    
    res = pmgRepeatTrials(theExpression, n)
    
    if((action <- svalue(.actions)) != "") {
      l = list(res)
      if(action %in% c("plot","hist","boxplot"))
        l = c(l, list(xlab="Simulation",ylab="", main=""))
      try(do.call(action, l), silent=TRUE)
    }
    
    saveAs = svalue(.saveAs)
    if(saveAs != "")
        assign_global(saveAs, res)
    
  })

  return(g)
}
                    
