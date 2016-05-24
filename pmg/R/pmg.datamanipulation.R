## Dialogs for data manipulation


## dialog for finding subsets
## return a group object
pmg.subset.dialog = function(container=NULL) {

  group  = ggroup(horizontal=FALSE, container=container)

  frame = gframe("<b>Data</b>",markup=TRUE, container=group)

  table = glayout()
  table[1,1] = glabel("x=")
  dataEntry = gedit("",width=30)
  table[1,2] = dataEntry
  ## subset
  table[2,1] = glabel("subset=")
  subsetEntry = gedit("",width=30)
  table[2,2] = subsetEntry
  subsetButton = gbutton("edit",handler = function(h,...) {
    editSubsetDialog(data=svalue(dataEntry),
                     widget=subsetEntry)})
  table[2,3] = subsetButton
  ## select
  table[3,1] = glabel("select=")
  selectEntry = gedit("", width=30)
  table[3,2] = selectEntry
  selectButton = gbutton("edit",handler = function(h,..) {
    editSelectDialog(data=svalue(dataEntry),
                     widget=selectEntry
                     )})
  table[3,3] = selectButton

  table[4,1] = glabel("drop=")
  dropEntry = gradio(c("TRUE","FALSE"),index=FALSE,selected=2)
  table[4,2] = dropEntry

  table[5,1] = glabel("assign to:")
  assignEntry = gedit("",width=30)
  table[5,2] = assignEntry
  
  submitButton = gbutton("submit",handler=function(h,...) {
    dataValue = svalue(dataEntry)
    subsetValue = svalue(subsetEntry)
    selectValue = svalue(selectEntry)
    dropValue = svalue(dropEntry)
    assignValue = svalue(assignEntry)
    if(assignValue == "")
      assignValue = NULL
    ## use pmg.cli to evaluate
    if(dataValue == "") {
      warning("No dataset chosen")
      return(NULL)
    }
    string = Paste("subset(", dataValue)
    if(nchar(subsetValue)>0) 
      string = Paste(string,", subset=",subsetValue)
    if(nchar(selectValue)>0)
      string = Paste(string, ", select=",selectValue)
    string = Paste(string,",drop=",dropValue,")")

    names(string) = assignValue
    svalue(pmg.cli) <- string

    ## close dialog?
    ##     if(!is.null(win))
    ##       dispose(win)
  })
  
  table[6,3] = submitButton

  add(frame, table,expand=TRUE)
  visible(table) <-  TRUE

  return(group)
}
  


##################################################

## edit data frame properties

pmg.edit.dataframe.properties.dialog = function(container=NULL) {

  ## in gWIdgetsRGtkw, but not exported?
  lsType = function(type, envir=.GlobalEnv) {
    x = with(.GlobalEnv, sapply(ls(), function(i) class(get(i))))
    objects = names(x)[sapply(x, function(i) any(i %in% type))]
  return(objects)
  }
  
  

  ## need means to select the data frame, popup this for editing

  g = ggroup(horizontal=FALSE, container=container)

  add(g, glabel("This dialog allows you to change the\n name, and data type for the columns of a data frame."))
  add(g, gseparator())
  
  tbl = glayout(container=g)

  allDFs = lsType("data.frame")
  selHandler = function(h,...) {
    newDFName = svalue(selectDF)
    .editDFProperties(newDFName)
  }
  selectDF = gdroplist(c("",allDFs), editable=TRUE, handler = selHandler)

  tbl[1,1] = glabel("Select a data frame:")
  tbl[2,1] = selectDF
  tbl[2,2] = gbutton("edit",handler=selHandler)
  
  visible(tbl) <- TRUE
  
  return(g)
  
}



.editDFProperties = function(dfname,envir=.GlobalEnv) {
  dlg = BasicGUI$new(sprintf("Edit properties for %s.",dfname))

  ## some defs
  dlg$allTypes = c("","numeric","integer","character","factor","logical")
  dlg$getType = function(.,i) head(class(i),n=1)



  ## validate name
  df = try(get(dfname,envir=envir), silent=TRUE)
  if(inherits(df,"try-error")) {
    cat("Need to have a data frame name.\n")
    return()
  }


  ## Store the data. We make changes to df as we update
  dlg$dfname <- dfname
  dlg$df <- df                          # make a copy
  
  dlg$colTypes = function(.) sapply(.$df,function(i) i$getType)

  
  ## Display dialog
  dlg$makeBody = function(., container) {
    g <- ggroup(horizontal=FALSE, container=container, expand=TRUE)
    glabel(gettext("Edit names and column types"),container=g)
    tbl <- glayout(container=g, expand=TRUE)

    tbl[1,1] <- "Which column:"         # no ?, : ala apple
    tbl[1,2] <- (.$columnDroplistGroup <- ggroup(container=tbl))
    .$columnDroplist <- gdroplist(names(.$df), container=.$columnDroplistGroup)
    
    tbl[2,1] <- "Column type:"
    tbl[2,2] <- (.$columnTypeDroplist <- gdroplist(.$allTypes, container=tbl))
    svalue(.$columnTypeDroplist) <- .$getType(.$df[,1]) ## initialize

    tbl[3,1] <- "Column name:"
    tbl[3,2] <- (.$columnNameDroplist <- gedit(names(.$df)[1], container=tbl))

    visible(tbl) <- TRUE

    gseparator(container=g)

    ## Show the current data frame
    .$dfGroup = ggroup(container=g, expand=TRUE)
    .$dfShow = gdf(head(.$df), container=.$dfGroup,expand=TRUE)
#    enabled(.$dfShow) <- FALSE
    
    bg <- ggroup(container=g)
    glabel("Save data frame as", container=bg)
    .$saveName <- gedit(.$dfname, container=bg)

    ## helper
    getIndex = function(.) {
      svalue(.$columnDroplist, index=TRUE)
    }      


    
    ## Now add handlers
    ## change colType
    addHandlerChanged(.$columnTypeDroplist,action=.,
                      handler = function(h,...) {
      ## get current var by index
      . = h$action
      ind = getIndex(.)
      coerceTo = svalue(h$obj)
      ## commit change
      .$df[,ind] <- do.call(paste("as.",coerceTo,sep="",collapse=""),
                            list(.$df[,ind]))
      updateDF(.)
    })

    addHandlerChanged(.$columnNameDroplist, action=.,
                      handler = function(h,...) {
                        ## handler for updating names
                        . = h$action
                        ## get variable index
                        ind = getIndex(.)

                        newName = make.names(svalue(h$obj))
                        ## validate
                        if(newName %in% names(.$df)) {
                          ## uniqueness not essential for a data frame
                          ## but is included here. Could leave out
                          cat(gettext("Specify a new, unique column name\n"))
                          return(FALSE)
                        }
                        ## aok
                        names(.$df)[ind] <- newName
                        ## update things
                        updateDF(.)
                        updateNames(.,ind)
                      })
    ## replace data frame in disply
    updateDF = function(.) {
      delete(.$dfGroup, .$dfShow)
      .$dfShow = gdf(head(.$df), container=.$dfGroup, expand=TRUE)
#      enabled(.$dfShow) <- FALSE
    }

    ## name droplist is trickier
    updateNames = function(.,ind) {
      ## now update names droplist, add in handler
      delete(.$columnDroplistGroup,.$columnDroplist)
      .$columnDroplist <- gdroplist(names(.$df),container=.$columnDroplistGroup)
      svalue(.$columnDroplist,index=TRUE) <- ind
      addHandlerChanged(.$columnDroplist, action=.,
                        handler = function(h,...) {
                          . = h$action
                          ind = svalue(h$obj, index=TRUE)
                          svalue(.$columnTypeDroplist) <- .$getType(.$df[,ind])
                          svalue(.$columnNameDroplist) <- names(.$df)[ind]
                        })
    }
    ## initialize -- to add handler, isn't added above
    updateNames(.,1)
  }
  
  dlg$clearButtonHandler = NULL
  dlg$okButtonHandler = function(.,h,...) {
    ## verify that name is okay
    outputName = svalue(.$saveName)
    if(outputName %in% ls(envir=.GlobalEnv)) {
      out = gconfirm(sprintf("There is already a varibable named %s. Overwrite?", outputName))
      if(!out) {
        return(FALSE)
      }
    }
    ## write
    assign_global(outputName, df)
    ## close up
    dispose(.$window)
  }


  ## now sohw the dialog
  dlg$show()
}

