getDataSets = function(...) {
  dataSets = data()$results
  dataSets = dataSets[, c(3,1,4)]
  return(dataSets)
}

## uses data() to move data set into environment
pmg.viewDataSets = function(width=550, height=400) {
  
  win = pmgWC$new("Load data set", v=T)
  size(win) <-  c( width, height)
  group = ggroup(horizontal=FALSE, container=win, expand=TRUE)

  

  dataSetHandler = function(h,...) {
    dataSets = svalue(dataSetList, drop=FALSE)
    for(i in 1:nrow(dataSets)) {
      dataset = dataSets[i,1]
      package = dataSets[i,2]
      command = Paste("data(",dataset,",package=\"",package,"\")")
      cat(pmg.prompt,command,"\n")
      svalue(status) <- Paste("attach data set ",dataset)
      do.call("data",list(dataset, package=package))
      svalue(status)
    }
  }

  
  dataSetList = gtable(getDataSets(), multiple=TRUE, filter.column = 2,
    handler =  dataSetHandler)
  add(group, dataSetList, expand=TRUE)

  ## add buttons
  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
  gbutton("cancel",container=buttonGroup, handler = function(h,...) dispose(win))

  status = gstatusbar("Double click data set to load",container=group)

  ## return window if desired -- can use destroy then.
  invisible(win)
}
