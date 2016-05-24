#############################################################################################
## Project: PKgraph
## File: functionOK.R
## Author: Xiaoyong Sun
## Date: 11/02/2009
## Goal: note
##        - interface between function.R and proto.R
## Notes:

################################################################################
## global function
################################################################################
getSubHeight <- function() .pk$getSubHeight()
getSubWidth <- function() .pk$getSubWidth()

getTerm <- function() .pk$getTerm()
setTerm <- function(term.df) .pk$setTerm(term.df)

getValidateData <- function() .pk$getValidateData()
setValidateData <- function(vdata) .pk$setValidateData(vdata)

getNameDataSpecialPlot <- function() .pk$getNameDataSpecialPlot()


## for some ggobi data, it is not directly from main screen. data need some modification,
## e.g. absoluate value, etc.
getDataSpecialPlot <- function(i) .pk$getDataSpecialPlot(i)
setDataSpecialPlot <- function(tdata, tname) .pk$setDataSpecialPlot(tdata, tname)
cleanDataSpecialPlot <- function() .pk$cleanDataSpecialPlot()

getDataLayoutPlot <- function() .pk$getDataLayoutPlot()
setDataLayoutPlot <- function(tdata) .pk$setDataLayoutPlot(tdata)
cleanDataLayoutPlot <- function() .pk$cleanDataLayoutPlot()

getDatasets <- function() .pk$getDatasets()
setDatasets <- function(tmp.data, dataname) .pk$setDatasets(tmp.data, dataname)

getTotalDataLen <- function() .pk$getTotalDataLen()
getTotalDataName <- function() .pk$getTotalDataName()

getCurrentData <- function(currentNo) 
{
  if (missing(currentNo))
  {
     currentMain <- svalue(pmg.dialog.notebook)
    .pk$getCurrentData(currentMain)
  }
  else
  {
    .pk$getCurrentData(currentNo)
  }

}
getCurrentDataType <- function(currentNo) .pk$getCurrentDataType(currentNo)
setCurrentDataType <- function(thisDataType, dataname) .pk$setCurrentDataType(thisDataType, dataname)

getItDataName <- function() .pk$getItDataName()
setItDataName <- function(itname) .pk$setItDataName(itname)

getItMap <- function() .pk$getItMap()
setItMap <- function(key) .pk$setItMap(key)

getComDataName <- function() .pk$getComDataName()
setComDataName <- function(itname) .pk$setComDataName(itname)

getComMap <- function() .pk$getComMap()
setComMap <- function(key.df) .pk$setComMap(key.df)

getCall <- function(package, plotType)
{
    if (package == "lattice")
    {
        if (plotType=="hist") mycall <- "histogram"
        else if (plotType=="scatter") mycall <- "xyplot"
        else if (plotType=="smatrix") mycall <- "splom"
        else if (plotType=="bwplot")  mycall <- "bwplot"
    }
    else
    {
        if (plotType=="hist") mycall <- "qplot"
        else if (plotType=="scatter") mycall <- "qplot"
        else if (plotType=="smatrix") mycall <- "plotmatrix"
    }
    return(mycall)
}

getPKCode <- function(i) .pk$getPKCode(i)
setPKCode <- function(newlist) .pk$setPKCode(newlist)

cleanPKCode <- function() .pk$cleanPKCode()

getPKGGobi <- function(i) .pk$getPKGGobi(i)
setPKGGobi <- function(newx) .pk$setPKGGobi(newx)
cleanPKGGobi <- function() .pk$cleanPKGGobi()


getSaveFormat <- function() .pk$getSaveFormat()
setSaveFormat <- function(newformat) .pk$setSaveFormat(newformat)

getFigConfig <- function() .pk$getFigConfig()
setFigConfig <- function(newconfig) .pk$setFigConfig(newconfig)

getGGobiPlotType <- function(currentNo) .pk$getGGobiPlotType(currentNo)  
setGGobiPlotType <- function(typelist, dataname) .pk$setGGobiPlotType(typelist, dataname)   

ggobiPlotType <- function()
{
    currentMain <- svalue(pmg.dialog.notebook)
    ggobi.map = gwindow("Configure interactive graphics", horizontal=FALSE)

    gtgroup1 = ggroup(cont=ggobi.map, horizontal=FALSE)

    gf1 <- gframe(text = "", markup = FALSE, pos = 0, horizontal=TRUE, container = gtgroup1)
    tbl <- glayout(cont=gf1)
        
      cline <- 0
      tbl.list <- list()
      tbl.list[[1]] = gdroplist(items=colnames(getCurrentData(currentMain)))
      tbl.list[[2]] = gdroplist(items=colnames(getCurrentData(currentMain)))
      tbl.list[[3]] = gdroplist(items=colnames(getCurrentData(currentMain)))

      name.order <- colnames(getCurrentData(currentMain))
      if (length(getGGobiPlotType()) != 0)
      {
          current.list <- getGGobiPlotType(currentMain)
          if (length(current.list) > 0)
          {
              ID.index <- match(current.list$ID, name.order)
              tbl.list[[1]] <- gdroplist(items=c(name.order[c(ID.index)], name.order[-ID.index]))
              Time.index <- match(current.list$Time, name.order)
              tbl.list[[2]] <- gdroplist(items=c(name.order[c(Time.index)], name.order[-Time.index]))
              Conc.index <- match(current.list$Conc, name.order)
              tbl.list[[3]] <- gdroplist(items=c(name.order[c(Conc.index)], name.order[-Conc.index]))
          }

      }
      
      dataname <- c("ID variable:", "Time variable:", "Concentration variable(DV):")
      
      for (i in 1: length(dataname))
      {
          cline <- cline + 1

          tbl[cline, 1, anchor = c(-1,-1)] = dataname[i]

          #tbl.list[[i]] = gdroplist(items=colnames(getCurrentData()))
          tbl[cline, 2, anchor = c(-1,-1)] = tbl.list[[i]]

      }

    gb1 = gbutton(text="Set PK data type", horizontal=FALSE )
    addhandlerclicked(gb1, function(h,...)
                    {
                        key <- list( ID=svalue(tbl.list[[1]]) , Time=svalue(tbl.list[[2]]),
                                     Conc=svalue(tbl.list[[3]]))

                        setGGobiPlotType(key, as.character(currentMain))
                        
                        tmp.data <- getCurrentData(currentMain)
                        tmp.name <- getTotalDataName()[currentMain]
                        tmp.data[,key$ID] <- factor(tmp.data[,key$ID])
                        setDatasets(tmp.data, tmp.name)
                        
                        dispose(ggobi.map)
                        svalue(pmg.statusBar) <- "Data types are configured successfully."                        

                    })
    cline <- cline + 1
    tbl[cline, 2, anchor = c(-1,-1)] = gb1
}

cleanAll <- function() .pk$cleanAll()

## Goal: repeat code for runing ggobi time series plot
ggobiRun <- function()       
{
      currentMain <- svalue(pmg.dialog.notebook)  
      ggobi.data <- getCurrentData(currentMain)

      Time.name <- getGGobiPlotType(as.character(currentMain))$Time
      Conc.name <- getGGobiPlotType(as.character(currentMain))$Conc
      ID.name <- getGGobiPlotType(as.character(currentMain))$ID

      if (!is.null(Time.name))
      {
          ## rearrange time, conc as the first and second variable
          time.ind <- which(colnames(ggobi.data)== Time.name)
          conc.ind <- which(colnames(ggobi.data)== Conc.name)
          old.ind <- c(1:length(colnames(ggobi.data)))
          ggobi.data <- ggobi.data[c(time.ind, conc.ind, old.ind[-c(time.ind, conc.ind)])]

          ggobi.text <- paste("g <- ggobi_longitudinal(ggobi.data", Time.name, ID.name, ")", sep=",")
          eval(parse(text=ggobi.text))
      }
      else
      {
          g <- ggobi(ggobi.data)
      }
      
      return(g)
}

## check data exist
checkDataExist <- function()
{
   current.no <- svalue(pmg.dialog.notebook)
   if (current.no == 0) return(FALSE)
   
   tmp.data <- getCurrentData(current.no)
   if (is.null(tmp.data))
   {
      return(FALSE)
   }
   else return(TRUE)
}

## check data type
checkDataType <- function(requiredType)
{
    currentDataNo <- svalue(pmg.dialog.notebook)
    dataType <- getCurrentDataType(currentDataNo)
    if (dataType==requiredType) return(TRUE)
    else return(FALSE)
}

## check data config
checkDataConfig <- function()
{
    if (nrow(getTerm()) > 0) return(TRUE)
    else return(FALSE)
}

checkSaveFormat <- function()
{
    if (length(getSaveFormat()) > 0) return(TRUE)
    else return(FALSE)
}

checkforPKModel <- function()
{
## check data exist
    if(!checkDataExist())
    {
        gmessage("No data is available for configuration!",
                icon=c("warning"), title="Warning")
        return(FALSE)
    }
## check dataType is "Data"
  if(!checkDataType(requiredDataType.PKmodel))
  {
      gmessage("The current data set is NOT DATA type!",
              icon=c("warning"), title="Warning")
      return(FALSE)
  }

## check data configure
  if (!checkDataConfig())
  {
      gmessage("Data is Not configured. Please use Menu: PK Models, Configure model result!",
              icon=c("warning"), title="Warning")
      return(FALSE)
  }
  return(TRUE)
}

ErrorMessage <- function(mymessage)
{
    gmessage(mymessage, icon="warning", title="warning")
}


################################################################################
## OK button for function.R and proto.R
################################################################################

## "cond" mainly for lattice, NOT for ggplot2
addList1 <- function(mylist, x, package, cond)
{
   if (package == "lattice")
   {
        if (missing(cond)) mylist$x <- formula(paste("~", x, sep=""))
        else mylist$x <- formula(paste("~", x, " | ", cond, sep=""))

   }
   else
   {
        mylist$x <- as.name(x)

   }

   mylist$xlab <- x
   #mylist$main <- x

   return(mylist)
}

addList2 <- function(mylist, x, y, package, cond)
{
   if (package == "lattice")
   {
        if (missing(cond)) mylist$x <- formula(paste(y, "~", x, sep=""))
        else mylist$x <- formula(paste(y, "~", x, " | ", cond, sep=""))
        mylist$xlab <- x
        mylist$ylab <- y
        #mylist$main <- paste(y, "vs", "x", sep=" ")
   }
   else
   {
        mylist$x <- as.name(x)
        mylist$y <- as.name(y)
        mylist$xlab <- x
        mylist$ylab <- y
        #mylist$main <- paste(y, "vs", "x", sep=" ")

   }
        return(mylist)
}

## TODO later; now only for individual plots
translateList <- function(para.list, config.list)
{
    mylist <- list()
    
    ## check data x exist
    if (para.list[["x"]]=="")
    {
        gmessage("Please choose x value!", icon=c("warning"), title="Warning")
        return(mylist)
    }

    mypackage <- config.list[["graphics"]]
    
    # decide call & list
    if (mypackage == "lattice")
    {
        ## decide call
        if (para.list[["y"]]=="") mycall <- "histogram"
        else mycall <- "xyplot"

        ## decide list
        #   - choose x or y
        #   - choose condition

        if (para.list[["y"]]=="")
        {
            if (config.list[["cond"]]=="") mylist <- addList1(mylist, para.list[["x"]], mypackage)
            else mylist <- addList1(mylist, para.list[["x"]], mypackage, config.list[["cond"]])
        }
        else
        {
            if (config.list[["cond"]]=="") mylist <- addList2(mylist, para.list[["x"]], para.list[["y"]],mypackage)
            else mylist <- addList2(mylist, para.list[["x"]], para.list[["y"]], mypackage, config.list[["cond"]])
        }
        
        #   - choose type
        if (para.list[["type"]]!="")
        {
            mylist$type <- para.list[["type"]]
        }

        #   - choose layout
        mylist$layout <- as.numeric(c(config.list[["layout_x"]], config.list[["layout_y"]]))
        
    }
    else
    {
        mycall <- "qplot"

## REMEMBER:
        #   - choose condition is DONE in the last step, since do.call syntax NOT work for ggplot2
        #   - choose x or y
        if (para.list[["y"]]=="")
        {
            mylist <- addList1(mylist, para.list[["x"]], mypackage)
        }
        else
        {
            mylist <- addList2(mylist, para.list[["x"]], para.list[["y"]], mypackage)
        }
        
        #   - choose type
        if (para.list[["type"]]!="")
        {
            mylist$geom <- switch(para.list[["type"]],
                                    p = c("point"),
                                    l = c("line"),
                                    psmooth = c("point", "smooth"),
                                    lsmooth = c("line", "smooth"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram")
                                    )
        }
        
    }
    
    #   - choose xlab, ylab
    if (para.list[["xlab"]]!="") mylist$xlab <- para.list[["xlab"]]
    if (para.list[["ylab"]]!="") mylist$xlab <- para.list[["ylab"]]

    #   - choose main
    if (para.list[["main"]]!="") mylist$xlab <- para.list[["main"]]

    return(list(mycall=mycall, mylist=mylist))
}

extractSimData <- function(dir.path, folder.name, file.name, id.var, cond.var)
{
        final.df <- NULL
        
        old.pwd <- getwd()
        on.exit(setwd(old.pwd))
        try1 <- try(setwd(dir.path), silent=TRUE)
        if (inherits(try1, "try-error"))
        {
            ErrorMessage("Path is wrong! Please input right path.")
            return(invisible(NULL))
        }
        else
        {
            
            target.dir <- dir(pattern=folder.name)
            sapply(1:length(target.dir), function(i)
                  {
                      if (file.name %in% list.files(path=target.dir[i]))
                      {
                          fileName <- paste(target.dir[i], file.name, sep="/")
                          sim.data <- read.table(fileName, header=T, skip=1)
                          mydata <- unique(sim.data[,c(id.var, cond.var)])

                          ## if cond.var is non-subject-specific,
                          ## - abs.sum or abs.sum^2
                          if (nrow(mydata) != length(unique(mydata[[id.var]])))
                          {
                              mydata[[cond.var]] <- abs(mydata[[cond.var]])
                              mean.tmp <- tapply(mydata[[cond.var]], list(mydata[[id.var]]), mean, simplify=TRUE, na.rm=TRUE)
                              mydata <- data.frame(unique(mydata[[id.var]]), mean.tmp)
                              colnames(mydata) <- c(id.var, cond.var)
                          }

                      }
                      else
                      {
                          ErrorMessage("There is no such file in these folders!")
                          return(invisible(NULL))
                      }
                      rowLabel <- "resample"
                      thisName <- paste(rowLabel, i, sep="")
                      final.df[[thisName]] <<- mydata[[cond.var]]
                      return(invisible(NULL))

                  })
        }
        
        return(final.df)
}

extractCddData <- function(dir.path, folder.name, file.name, id.var, cond.var, total.id, rowLabel)
{
        final.df <- NULL
        v.delete.id <- NULL
        
        old.pwd <- getwd()
        on.exit(setwd(old.pwd))
        try1 <- try(setwd(dir.path), silent=TRUE)
        if (inherits(try1, "try-error"))
        {
            ErrorMessage("Path is wrong! Please input right path.")
            return(invisible(NULL))
        }
        else
        {

            target.dir <- dir(pattern=folder.name)
            if (length(target.dir)==0)      
            {
                ErrorMessage("In the Target directory, there is NO such folder pattern!")
                return(invisible(NULL))
            }
                        
            sapply(1:length(target.dir), function(i)
                  {
                      if (file.name %in% list.files(path=target.dir[i]))
                      {
                          fileName <- paste(target.dir[i], file.name, sep="/")
                          sim.data <- read.table(fileName, header=T, skip=1)
                          if (!all(c(id.var, cond.var) %in% colnames(sim.data)))
                          {
                              ErrorMessage(paste("In simulation folder - ", target.dir[i], ", Patiend ID or Plot variable
                                           does not match the NONMEM result file !", sep=""))
                              return(invisible(NULL))
                          }                          
                          
                          mydata <- unique(sim.data[,c(id.var, cond.var)])

                          ## if cond.var is non-subject-specific,
                          ## - abs.sum or abs.sum^2
                          if (nrow(mydata) != length(unique(mydata[[id.var]])))
                          {
                              mydata[[cond.var]] <- abs(mydata[[cond.var]])
                              mean.tmp <- tapply(mydata[[cond.var]], list(mydata[[id.var]]), mean, simplify=TRUE, na.rm=TRUE)
                              mydata <- data.frame(unique(mydata[[id.var]]), mean.tmp)
                              colnames(mydata) <- c(id.var, cond.var)
                          }

                          ## since it is cdd file, some ids are deleted;
                          ## to get same length for original data, we add NA to position of deleted IDs
                          delete.id <- which(!(unique(total.id) %in% unique(mydata[[id.var]])))
                          if (length(delete.id) > 0)
                          {
                              v.delete.id <<- c(v.delete.id, delete.id)
                              delete.id.df <- data.frame(delete.id, rep(NA, length(delete.id)))
                              colnames(delete.id.df) <- colnames(mydata)
                              mydata <- rbind(mydata, delete.id.df)
                          }

                          mydata <- mydata[order(mydata[[id.var]]),]
                      }
                      else
                      {
                          ErrorMessage(paste("In simulation folder - ", target.dir[i], ", there is NO such NONMEM result file!", sep=""))
                          return(invisible(NULL))
                      }

                      thisName <- paste(rowLabel, i, sep="")
                      final.df[[thisName]] <<- mydata[[cond.var]]
                      return(invisible(NULL))

                  })
        }

        return(list(data=final.df, deleteID=v.delete.id))
}

extractBootData <- function(dir.path, folder.name, file.name, id.var, cond.var, bootKey.table, totalID, missingIDValue)
{
        final.df <- NULL

        old.pwd <- getwd()
        on.exit(setwd(old.pwd))
        try1 <- try(setwd(dir.path), silent=TRUE)
        if (inherits(try1, "try-error"))
        {
            ErrorMessage("Path is wrong! Please input right path.")
            return(invisible(NULL))
        }
        else
        {

            real.dir <- dir(pattern=folder.name)
            if (length(real.dir) < 1)
            {
                ErrorMessage(paste("There is no folder named as",folder.name, sep=" "))
                return(invisible(NULL))
            }

            target.dir <- gsub(folder.name, "", real.dir)
            target.dir <- order(target.dir)
            
            sapply(1:length(target.dir), function(i)
                  {
                      real.dir <- paste(folder.name,  target.dir, sep="")
                      if (file.name %in% list.files(path=real.dir[i]))
                      {

                          fileName <- paste(real.dir[i], file.name, sep="/")
                          sim.data <- read.table(fileName, header=T, skip=1)
            
                          if (!all(c(id.var, cond.var) %in% colnames(sim.data)))
                          {
                              ErrorMessage(paste("In simulation folder - ", real.dir[i], ", Patiend ID or Plot variable
                                           does not match the NONMEM result file !", sep=""))
                              return(invisible(NULL))
                          }
                          
                          if (is.null(bootKey.table))
                          {
                              mydata0 <- unique(sim.data[,c(id.var, cond.var)])
                          }
                          else
                          {
                              mydata0 <- data.frame(rep(as.integer(bootKey.table[i,]), each=nrow(sim.data)/ncol(bootKey.table)), sim.data[,c(cond.var)])
                              colnames(mydata0) <- c(id.var, cond.var)
                          }

                          mydata <- unique(mydata0)

                          ## if cond.var is non-subject-specific: some subject ID matches more than value
                          ## - abs.sum or abs.sum^2

                          if (nrow(mydata) != length(unique(mydata[[id.var]])))
                          {
                              mydata[[cond.var]] <- abs(mydata[[cond.var]])
                              mean.tmp <- tapply(mydata[[cond.var]], list(mydata[[id.var]]), mean, simplify=TRUE, na.rm=TRUE)
                              mydata <- data.frame(unique(mydata[[id.var]]), mean.tmp)
                              colnames(mydata) <- c(id.var, cond.var)
                          }

                          ## FOR unbootstrap ID
                          delete.id <- which(! (c(totalID) %in% as.integer(bootKey.table[i,])))

                          if (missingIDValue == 0)  delete.mydata <- data.frame(delete.id, rep(0, length(delete.id)))
                          else delete.mydata <- data.frame(delete.id, rep(NA, length(delete.id)))

                          colnames(delete.mydata) <- c(id.var, cond.var)

                          mydata <- rbind(mydata, delete.mydata)

                          mydata <- mydata[order(mydata[[id.var]]),]

                      }
                      else
                      {
                          ErrorMessage(paste("In simulation folder - ", real.dir[i], ", there is NO such NONMEM result file!", sep=""))
                          return(invisible(NULL))
                      }
                      rowLabel <- "resample"
                      thisName <- paste(rowLabel, i, sep="")
                      final.df[[thisName]] <<- mydata[[cond.var]]

                      return(invisible(NULL))

                  })
        }

        return(final.df)
}

################################################################################
# Dialog, Table selection for dataset
################################################################################
# labelMessage1 <- "Please move data for diagnostics from left TABLE to right TABLE."
# labelMessage2 <- "After choosing data, click to interactive diagnostics."
# winTitle <- "Configure datasets"
# statusMessage <- "Data is ready for interactive diagnostics."
# menuOption: 1 for model comparison; 2 for interactive graphics

selectDataDialog <- function(winTitle, labelMessage1, labelMessage2, statusMessage, menuOption)
{
   # check data exist
    if(!checkDataExist())
    {
        gmessage("No data is available!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }

    alldata.name <- getTotalDataName()

    ggobi.gwin <- gwindow(title=winTitle)
    ggobi.group <- ggroup(cont=ggobi.gwin, horizontal=FALSE)
    g0 <- glabel(text=labelMessage1, cont=ggobi.group)
    g1 <- gframe(cont=ggobi.group)
        size(g1) <- c(getSubWidth()*0.5, getSubHeight()*0.5)
    #tbl1 = gtable(data.frame(AllData=alldata.name), cont = g1, expand=TRUE)
    tbl1 = gtable(alldata.name, cont = g1, expand=TRUE)
    size(tbl1) <- c(getSubWidth()*0.2, getSubHeight()*0.4)

    arrowButton = gbutton(">", cont = g1);
    enabled(arrowButton) <- FALSE

    #tbl2 = gtable(data.frame(TargetData=colnames(Theoph)), cont=g1, expand = TRUE)
    tbl2 = gtable(colnames(Theoph), cont=g1, expand = TRUE)
    tbl2[] <- c() ## clear out initialized values. Can't start gtable empty.
    size(tbl2) <- c(getSubWidth()*0.2, getSubHeight()*0.4)

    g3 <- glabel(text=labelMessage2, cont=ggobi.group)
    confirmButton = gbutton("Go to next step", cont= ggobi.group)


    addHandlerClicked(tbl1, function(h,...)
    {
      enabled(arrowButton) <- TRUE
      svalue(arrowButton) <- ">"
      ## no means to clear all selections
    })
    addHandlerClicked(tbl2, function(h,...)
    {
      enabled(arrowButton) <- TRUE
      svalue(arrowButton) <- "<"
    })

    addHandlerClicked(arrowButton, function(h,...)
    {
      if(svalue(arrowButton) == ">") {
        curVal = svalue(tbl1)
        tbl1[] <- setdiff(tbl1[], curVal)
        if(any(is.na(tbl2[])))
          tbl2[] <- curVal
        else
          tbl2[] <- sort(unique(c(tbl2[],curVal)))          # adjust for initial  NA
      } else {
        curVal = svalue(tbl2)
        tbl2[] <- setdiff(tbl2[], curVal)
        tbl1[] <- sort(unique(c(tbl1[],curVal)))
      }
      enabled(arrowButton) <- FALSE
    })

    addHandlerClicked(confirmButton, function(h,...)
    {
       it.data <- tbl2[]
       if(length(it.data)==0)
       {
          ErrorMessage("Please use button in the middle to choose data from left column")
          return(invisible(NULL))
       }
       
       if (menuOption==2)
       {
          # for interactive
          setItDataName(it.data)
       }
       else
       {
          # for model comparison
          if (length(it.data)!=2)
          {
              ErrorMessage("You need to pick TWO datasets for comparison!")
              return(invisible(NULL))
          }
          setComDataName(it.data)
       }
       svalue(pmg.statusBar) <- statusMessage
       dispose(ggobi.gwin)

    })


}
subset.f1 <- function(sid, ssign, svalue, tmp.data)
{
    check.value <- c(VarName=sid, Sign=ssign, Value=svalue)

    ## check data exist or not
    if (any(check.value==""))
    {

         ErrorMessage(paste(names(check.value[check.value==""]), " is Not chosen!", sep=""))
    }

    ## TODO: diff numerical value vs categorical value
    #check.value[3] <- as.numeric(check.value[3])

    if (check.value[2]=="==")
     {
          tmp1.data <- tmp.data[tmp.data[[check.value[1]]] == check.value[3],]
     }
     else if (check.value[2]==">")
     {
          tmp1.data <- tmp.data[tmp.data[[check.value[1]]] > as.numeric(check.value[3]),]
     }
     else if (check.value[2]==">=")
     {
          tmp1.data <- tmp.data[tmp.data[[check.value[1]]] >= as.numeric(check.value[3]),]
      }
     else if (check.value[2]=="<")
     {
          tmp1.data <- tmp.data[tmp.data[[check.value[1]]] < as.numeric(check.value[3]),]
     }
     else if (check.value[2]=="<=")
     {
          tmp1.data <- tmp.data[tmp.data[[check.value[1]]] <= as.numeric(check.value[3]),]
     }
     else if (check.value[2]=="!=")
     {
          tmp1.data <- tmp.data[tmp.data[[check.value[1]]] != check.value[3],]
     }

     
    return(tmp1.data)
}

subset.okButtonHandler = function(., h,...)
{
   subset.var <- svalue(.$table.widget)
   dispose(.$window)
   
   if (length(subset.var) < 1)
   {
      ErrorMessage("No data is chosen for subset!")
      return(invisible(NULL))
   }
   
   gwin.subset = gwindow("Subset")
   gframe.subset = gframe(text = "Subset", markup = FALSE, pos = 0, horizontal=FALSE, container = gwin.subset)
   tbl = glayout(cont=gframe.subset)
   tbl.value <- list()
   
   currentPage <- svalue(pmg.dialog.notebook)
   currentData <- getCurrentData(currentPage)
   
   for ( i in 1:length(subset.var))
   {
       tbl.value[[i*3-2]] = gdroplist(items=subset.var[i])
       tbl[i,1] = tbl.value[[i*3-2]]

       if (is.factor(currentData[[subset.var[i]]]) || is.character(currentData[[subset.var[i]]]))
          tbl.value[[i*3-1]] = gdroplist(items=c("==", "!="))
       else
          tbl.value[[i*3-1]] = gdroplist(items=c("==", ">", ">=", "<", "<="))

       tbl[i,2] = tbl.value[[i*3-1]]
       tbl.value[[i*3]] = gedit(text="")
       tbl[i,3] = tbl.value[[i*3]]
   }
   

   gbutton(text = "Subset", border=TRUE, cont=gframe.subset, handler = function(h,...)
           {
                  sapply(1:length(subset.var), function(i)
                        {
                            currentData <<- subset.f1(svalue(tbl.value[[i*3-2]]), svalue(tbl.value[[i*3-1]]),
                                                 svalue(tbl.value[[i*3]]), currentData)

                        })

                  thisDatano <- getTotalDataLen() + 1
                  datatype <- "Subset"
                  thisDataName <- paste(getTotalDataLen() + 1, "_", datatype , sep="")
                  
                  setDatasets(currentData, thisDatano) # use no as data name
                  setCurrentDataType(datatype, thisDataName)
                  
                  ptable=gtable(currentData, multiple=TRUE, expand=TRUE)
                  pkmain.add(ptable, as.character(thisDataName), override.closebutton = TRUE)

                  dispose(gwin.subset)
                  svalue(PKENV$pmg.statusBar) <- "Subset sucessfully!"
           })
}

factor.okButtonHandler = function(., h, ...)
{
    factor.var <- svalue(.$table.widget)
    if (length(factor.var) < 1)
    {
        ErrorMessage("No variable is chosen!")
        return(invisible(NULL))
    }
    
    currentPage <- svalue(pmg.dialog.notebook)
    currentData <- getCurrentData(currentPage)

    factor.call.value <- svalue(.$checkbox.widget)
    if (factor.call.value == "factor") factor.call <- "factor"
    else if (factor.call.value == "numeric") factor.call <- "as.numeric"
    else factor.call <- "as.character"
    ## factor
    sapply(1:length(factor.var), function(i)
          {
              currentData[[factor.var[i]]] <<- do.call(factor.call, list(currentData[[factor.var[i]]]))
          })

    thisDataName <- getTotalDataName()[currentPage]

    setDatasets(currentData, thisDataName) # use no as data name
    # No need to reset data type
    #setCurrentDataType(svalue(datatype), thisDataName)

    svalue(pmg.statusBar) <- "Factor/Unfactor sucessfully!"
    dispose(.$window)
}

saveImageHandler = function(.,h,...)
{
## check saving format is setup!
    if (!checkSaveFormat())
    {
        ErrorMessage("Please Menu configure to set save format first in main panel!")
        return(invisible(NULL))
    }

    currentMain <- svalue(pmg.dialog.notebook) 
    currentPage <- svalue(pk.dialog.notebook)
    pkcode <- getPKCode(currentPage)

    ## datagroup: data, boot, outlier
    # datagroup = 0, save single plot; otherwise, multiple plot
    # specialPlot, TRUE, - data, using others; otherwise, use getCurrentData()
    if (.$datagroup == 0)
    {
        # check data absolute or not
        abs.no <- getNameDataSpecialPlot()
        abs.check <- currentPage %in% abs.no
        if (currentPage %in% abs.no) pkcode$pklist$data <- getDataSpecialPlot(as.character(currentPage))
        else pkcode$pklist$data <- getCurrentData(currentMain)
    }
    else pkcode$pklist$data <- getValidateData()

    saveFormat <- getSaveFormat()

    gfile("Save file",type="save", handler = function(h,...)
      {
          mycommand <- saveFormat$command

          sapply(1:length(mycommand), function(i)
          {   
              # figureGroup: save one file or multiple file
              if (.$figureGroup == 0)
              {
                  layout.plot <- getDataLayoutPlot()
                  if (currentPage %in% layout.plot)  filename <- paste(h$file, "%3d.", mycommand[i], sep="")
                  else  filename <- paste(h$file, mycommand[i], sep=".")
              }
              else
                  filename <- paste(h$file, "%3d.", mycommand[i], sep="")

              if (!is.na(saveFormat$width) && !is.na(saveFormat$width))
                  mysavelist <- list(file=filename, width=saveFormat$width, height=saveFormat$height)
              else
                  mysavelist <- list(file=filename)

              do.call(mycommand[i], mysavelist)
              print(pkcode$pklist)
              dev.off()
          })
          
          
          focus(.$window) <- TRUE
          
      })
}

saveImageHandler.pkmodel = function(.,h,...)   
{
## check saving format is setup!
    if (!checkSaveFormat())
    {
        ErrorMessage("Please Menu configure to set save format first in main panel!")
        return(invisible(NULL))
    }
    currentMain <- svalue(pmg.dialog.notebook) 
    currentPage <- svalue(pk.dialog.notebook)
    pkcode <- getPKCode(currentPage)

    ## datagroup: data, boot, outlier
    # datagroup = 0, save single plot; otherwise, multiple plot
    # specialPlot, TRUE, - data, using others; otherwise, use getCurrentData()
    if (.$datagroup == 0)
    {
        mydata <- getPKGGobi(currentPage)
        x.name <- mydata$x        
        if ((length(x.name) > 1) && (svalue(.$savewd[["graphics"]]) != "lattice")) 
        {
            ## for scatterplot matrix, do nothing
        }
        else
        {   
            # check data absolute or not
            abs.no <- getNameDataSpecialPlot()
            abs.check <- currentPage %in% abs.no
            if (currentPage %in% abs.no) pkcode$pklist$data <- getDataSpecialPlot(as.character(currentPage))
            else pkcode$pklist$data <- getCurrentData(currentMain)
        }
    }
    else pkcode$pklist$data <- getValidateData()

    saveFormat <- getSaveFormat()

    gfile("Save file",type="save", handler = function(h,...)
      {
          mycommand <- saveFormat$command

          sapply(1:length(mycommand), function(i)
          {   
              # figureGroup: save one file or multiple file
              if (.$figureGroup == 0)
              {
                  layout.plot <- getDataLayoutPlot()
                  if (currentPage %in% layout.plot)  filename <- paste(h$file, "%3d.", mycommand[i], sep="")
                  else  filename <- paste(h$file, mycommand[i], sep=".")
              }
              else
                  filename <- paste(h$file, "%3d.", mycommand[i], sep="")

              if (!is.na(saveFormat$width) && !is.na(saveFormat$width))
                  mysavelist <- list(file=filename, width=saveFormat$width, height=saveFormat$height)
              else
                  mysavelist <- list(file=filename)

              do.call(mycommand[i], mysavelist)
              print(pkcode$pklist)
              dev.off()
          })
          
          
          focus(.$window) <- TRUE
          
      })
}

saveImageHandler.matrix = function(.,h,...)
{
## check saving format is setup!
    if (!checkSaveFormat())
    {
        ErrorMessage("Please Menu configure to set save format first in main panel!")
        return(invisible(NULL))
    }

    currentMain <- svalue(pmg.dialog.notebook) 
    currentPage <- svalue(pk.dialog.notebook)
    pkcode <- getPKCode(currentPage)

    ## datagroup: data, boot, outlier
    # datagroup = 0, save single plot; otherwise, multiple plot
    # specialPlot, TRUE, - data, using others; otherwise, use getCurrentData()

    if (.$datagroup == 0)
    {
        # THAT'S SPECIAL PART FOR ggplot2::plotmatrix. No need to replace data
        if (svalue(.$savewd[["graphics"]]) == "lattice")
        {
            abs.no <- getNameDataSpecialPlot()
            abs.check <- currentPage %in% abs.no
            if (currentPage %in% abs.no) pkcode$pklist$data <- getDataSpecialPlot(as.character(currentPage))
            else pkcode$pklist$data <- getCurrentData(currentMain)
        }

    }
    else pkcode$pklist$data <- getValidateData()

    saveFormat <- getSaveFormat()

    gfile("Save file",type="save", handler = function(h,...)
      {
          mycommand <- saveFormat$command

          sapply(1:length(mycommand), function(i)
          {
              # figureGroup: save one file or multiple file
              if (.$figureGroup == 0)
              {
                  layout.plot <- getDataLayoutPlot()
                  if (currentPage %in% layout.plot)  filename <- paste(h$file, "%3d.", mycommand[i], sep="")
                  else  filename <- paste(h$file, mycommand[i], sep=".")
              }
              else
                  filename <- paste(h$file, "%3d.", mycommand[i], sep="")

              if (!is.na(saveFormat$width) && !is.na(saveFormat$width))
                  mysavelist <- list(file=filename, width=saveFormat$width, height=saveFormat$height)
              else
                  mysavelist <- list(file=filename)

              do.call(mycommand[i], mysavelist)
              print(pkcode$pklist)
              dev.off()
          })


          focus(.$window) <- TRUE

      })
}

model.ggobiImageHandler = function(.,h,...)
{
      # get ggobi <-> pk.dialog page: variable list
      currentPage <- svalue(pk.dialog.notebook)
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }
      
      mydata <- getPKGGobi(currentPage)

        # check data absolute or not
      abs.no <- getNameDataSpecialPlot()
      currentMain <- svalue(pmg.dialog.notebook) 
      if (is.null(abs.no)) tmp.data <- getCurrentData(currentMain)
      else
      {
            abs.check <- currentPage %in% abs.no
            if (currentPage %in% abs.no) tmp.data <- getDataSpecialPlot(as.character(currentPage))
            else tmp.data <- getCurrentData()
      }

      x.name <- mydata$x      
      y.name <- mydata$y

      if (is.null(y.name))
      {
          g <- ggobi(tmp.data)
          if (length(x.name) == 1)
          {  
             display(g[1], pmode="Barchart", vars=list(X = mydata$x))
          }
          else # for scatterplot matrix
          {
             display(g[1], "Scatterplot Matrix", vars=list(X = mydata$x) )
          }
      }
      else
      {
	x.ind <- which(colnames(tmp.data)== x.name)
          y.ind <- which(colnames(tmp.data)== y.name)
          old.ind <- c(1:length(colnames(tmp.data)))
          old.ind <- old.ind[-c(x.ind, y.ind)]

          tmp.data <- tmp.data[c(x.ind, y.ind, old.ind)]

          g <- ggobi(tmp.data)
      }
}

validation.ggobiImageHandler = function(.,h,...)
{
      # get ggobi <-> pk.dialog page: variable list
      currentPage <- svalue(pk.dialog.notebook)
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }
      
      mydata <- getPKGGobi(currentPage)$x
      g <- ggobi(mydata)
      display(g[1], "Parallel Coordinates Display")
}

psn.outlier.ggobiImageHandler = function(.,h,...)
{
    ErrorMessage("No ggobi instance available for this data.")
}


vis.outlier.ggobiImageHandler = function(.,h,...)
{
      currentPage <- length(getNameDataSpecialPlot()) # new
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }
      
      id.var <- svalue(.$widgets[["Patient ID:"]])
      
      mydata <- getDataSpecialPlot(as.character(currentPage))
      # make it factor, so can link all together
      mydata[[id.var]] <- factor(mydata[[id.var]])
      
      currentMain <- svalue(pmg.dialog.notebook)           
      oridata <- getCurrentData(currentMain)
      oridata[[id.var]] <- factor(oridata[[id.var]])
      
      # replace name, so that all plots have same ID name to link
      oriname <- colnames(oridata)
      if (id.var %in% oriname)
      {
          oriname <- gsub(id.var, "resampleID", oriname)
          colnames(oridata) <- oriname
      }

      Time.name <- getGGobiPlotType(as.character(currentMain))$Time
      Conc.name <- getGGobiPlotType(as.character(currentMain))$Conc
      ID.name <- getGGobiPlotType(as.character(currentMain))$ID
      if (length(Time.name)==0 || length(Conc.name)==0 || length(ID.name)==0)
      {
          ErrorMessage("Please make sure the input file type is PK data!")
          return(invisible(NULL))
      }

      Time.ind <- which(colnames(oridata)==Time.name)
      Conc.ind <- which(colnames(oridata)==Conc.name)
      oridata <- cbind(oridata[,c(Time.ind, Conc.ind)], oridata[,-c(Time.ind, Conc.ind)])

      gtext <- paste("ggobi_longitudinal(oridata,", Time.name, ", resampleID)", sep="")
      g <- eval(parse(text=gtext))
            
      # not g[2], because the edges of longitudinal data is g[2] 
      g["moreplot"] <- mydata
      display(g[3], pmode="Scatterplot Display",
                 vars=list(X="cor1", Y="cor2"))
      display(g[3], pmode="Scatterplot Display",
                 vars=list(X="resampleID", Y="para"))
}

boot.vis.ggobiImageHandler = function(.,h,...)
{
      currentPage <- length(getNameDataSpecialPlot()) 
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }
      
      id.var <- svalue(.$widgets[["Patient ID:"]])

      mydata <- getDataSpecialPlot(as.character(currentPage))

      # make it factor, so can link all together
      mydata[[id.var]] <- factor(mydata[[id.var]])

      currentMain <- svalue(pmg.dialog.notebook)
      oridata <- getCurrentData(currentMain)
      Time.name <- getGGobiPlotType(as.character(currentMain))$Time
      Conc.name <- getGGobiPlotType(as.character(currentMain))$Conc
      ID.name <- getGGobiPlotType(as.character(currentMain))$ID

      if (length(Time.name)==0 || length(Conc.name)==0 || length(ID.name)==0)
      {
          ErrorMessage("Please make sure the input file type is PK data!")
          return(invisible(NULL))
      }

      Time.ind <- which(colnames(oridata)==Time.name)
      Conc.ind <- which(colnames(oridata)==Conc.name)
      oridata <- cbind(oridata[,c(Time.ind, Conc.ind)], oridata[,-c(Time.ind, Conc.ind)])      

      oridata[[id.var]] <- factor(oridata[[id.var]])
      gtext <- paste("ggobi_longitudinal(oridata,", Time.name, ",", ID.name, ")", sep="")
      g <- eval(parse(text=gtext))

      g["moreplot"] <- mydata
      # not g[2], because the edges of longitudinal data is g[2]
      display(g[3], pmode="Scatterplot Display",
                 vars=list(X="ID", Y="VAR"))

}


getHistCall <- function(hist.graphics, hist.x,
                hist.bin,
                hist.main,
                hist.xlab,
                hist.ylab,
                hist.type,
                hist.cond,
                hist.layout_x,
                hist.layout_y,
                hist.data
                )
{
    hist.col <- getFigConfig()$col

    currentMain <- svalue(pmg.dialog.notebook)  
    if (missing(hist.data)) hist.data <- getCurrentData(currentMain)  

    if (hist.graphics == "lattice")
    {
        x <- paste( "~", hist.x, sep="")

        if ( !is.null(hist.bin) && hist.bin != "" )
        {
            if ( hist.cond == "")
            {
                x <- as.formula(x)
                lattice.final <- histogram(x=x, xlab=hist.xlab, nint = as.numeric(hist.bin),
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                main=hist.main , data= hist.data )
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- histogram(x=x, xlab=hist.xlab, nint = as.numeric(hist.bin),
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }
        }
        else
        {
            if ( is.null(hist.cond) || hist.cond =="" )
            {
                x <- as.formula(x)
                lattice.final <- histogram(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type,  col=hist.col,
                                main=hist.main , data= hist.data)
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- histogram(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }

        }
        
        return(lattice.final)

     }
     else ## start ggplot2
     {

        mytype <- switch(hist.type,
                                    #p = c("point"),
                                    #l = c("line"),
                                   # b = c("point", "smooth"),
                                    #psmooth = c("point", "smooth"),
                                    #lsmooth = c("line", "smooth"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram"),
                                    c("histogram")
                                    )

        myx <- hist.data[[hist.x]]

        if ( !is.null(hist.bin) && hist.bin != "" )
        {
            checkbin <- 1
            f <- as.numeric(diff(range(myx)/as.numeric(hist.bin)))
        }
        else checkbin <- 0
        
        if ( !is.null(hist.cond) && hist.cond != "" )
        {
            checkcond <- 1
            if (hist.layout_x == "" || hist.layout_y == "")
            {
               ErrorMessage("layout_x or layout_y is NOT specified for conditional variable")
               return(invisible(NULL))
            }

            mycond <- hist.data[[hist.cond]]
            total.figno <- as.numeric(hist.layout_x) * as.numeric(hist.layout_y)

            if (length(unique(mycond)) > total.figno)
            {
                newno <- unique(mycond)[1:total.figno]
                part.data <- hist.data[which(mycond%in%newno),]
            }
            else  part.data <- hist.data

        }
        else checkcond <- 0
        
        if (checkcond == 1)
        {
            if (checkbin == 1)
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab,  #colour = hist.col,
                                    main=hist.main, data=hist.data)+ geom_histogram(binwidth = f)+ facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))
            }
            else
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = mytype, #colour = hist.col,
                                    main=hist.main, data=hist.data)+ facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))
            }

        }
        else
        {
            if (checkbin == 1)
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab,  #colour = hist.col,
                                    main=hist.main, data=hist.data) + geom_histogram(binwidth = f)
            }
            else
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = mytype, #colour = hist.col,
                                    main=hist.main, data=hist.data)
            }
        }


        return(ggplot.final)
     }
}

getScatterCall <- function(hist.graphics, hist.x, hist.y,
                hist.main,
                hist.xlab,
                hist.ylab,
                hist.type,
                hist.cond,
                hist.layout_x,
                hist.layout_y, hist.data
                )
{
    currentMain <- svalue(pmg.dialog.notebook)  
    if (missing(hist.data)) hist.data <- getCurrentData(currentMain)  

    hist.col <- getFigConfig()$col
    
    if ((missing(hist.type)) || hist.type=="") hist.type = "p"
    if (!is.null(getFigConfig()$loess) && getFigConfig()$loess == 1) hist.type = c(hist.type, "smooth")

    ## two graph packages
    if (hist.graphics == "lattice")
    {
        x <- paste(hist.y, "~", hist.x, sep="")

            if ( hist.cond == "")
            {
                x <- as.formula(x)
                lattice.final <- xyplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                main=hist.main , data= hist.data)
                ## type: ts
                if (hist.type == "ts" && (!is.null(getGGobiPlotType(currentMain)$ID)) )
                {
                    lattice.final <- xyplot(x=x, xlab=hist.xlab, ylab=hist.ylab, 
                                  type=c("l","p"), groups= hist.data[,getGGobiPlotType(currentMain)$ID],
                                  main=hist.main , data= hist.data)
                 }
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- xyplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

                ## type: ts
                if (hist.type == "ts" && (!is.null(getGGobiPlotType(currentMain)$ID)) )
                {
                    lattice.final <- xyplot(x=x, xlab=hist.xlab, ylab=hist.ylab, 
                                type=c("l","p"), groups= hist.data[,getGGobiPlotType(currentMain)$ID],
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)
                 }
            }

        return(lattice.final)

     }
     else  ## start ggplot
     {
        tmp.type <- NULL
        if (length(hist.type) > 1)
        {
            tmp.type <- hist.type
            hist.type <- hist.type[1]
        }
        hist.type <- switch(hist.type,
                                    p = c("point"),
                                    l = c("line"),
                                    b = c("point", "line"),
                                    ts = c("ts"),
                                    #psmooth = c("point", "smooth"),
                                    #lsmooth = c("line", "smooth"),
                                    a = c("point", "line"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram")
                                    )
        if (!is.null(tmp.type)) hist.type <- c(hist.type, tmp.type[-1])

        myx <- hist.data[[hist.x]]
        myy <- hist.data[[hist.y]]

        if ( !is.null(hist.cond) && hist.cond != "" )
        {
            # TODO
            if (hist.layout_x == "" || hist.layout_y == "")
            {
               ErrorMessage("layout_x or layout_y is NOT specified for conditional variable")
               return(invisible(NULL))
            }

            mycond <- hist.data[[hist.cond]]
            total.figno <- as.numeric(hist.layout_x) * as.numeric(hist.layout_y)

            if (length(unique(mycond)) > total.figno)
            {
                newno <- unique(mycond)[1:total.figno]
                part.data <- hist.data[which(mycond%in%newno),]
            }
            else part.data <- hist.data

             ## type: ts
             if (hist.type == "ts" && (!is.null(getGGobiPlotType(currentMain)$ID)) )
             {
                ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab, ylab=hist.ylab, 
                                    geom = c("point", "line"), colour = part.data[,getGGobiPlotType(currentMain)$ID],
                                    main=hist.main , data= part.data, se=F) + facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x)) + opts(legend.position="none")
             
             }
             else
             {
                ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = hist.type, #colour = hist.col,
                                    main=hist.main , data= part.data, se=F) + facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))
             
             }
        }
        else
        {

             if (hist.type == "ts" && (!is.null(getGGobiPlotType(currentMain)$ID)) )
             {
                ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab, ylab=hist.ylab, 
                                    geom = c("point", "line"), colour = hist.data[,getGGobiPlotType(currentMain)$ID],
                                    main=hist.main , data= hist.data, se=F) + opts(legend.position="none")            
             }
             else
             {
                ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = hist.type, #colour = hist.col,
                                    main=hist.main , data= hist.data, se=F)             
             }                                    
        }

        return(ggplot.final)
     }
}

getScatterCall.back <- function(hist.graphics, hist.x, hist.y,
                hist.main,
                hist.xlab,
                hist.ylab,
                hist.type,
                hist.cond,
                hist.layout_x,
                hist.layout_y, hist.data
                )
{
    currentMain <- svalue(pmg.dialog.notebook)  
    if (missing(hist.data)) hist.data <- getCurrentData(currentMain)  

    hist.col <- getFigConfig()$col
    if (hist.type == "")
    {

        if (is.null(getFigConfig()$loess)) hist.type = "p"
        else if (getFigConfig()$loess == 1) hist.type = c("p", "smooth")
        else
        {
            ErrorMessage("Wrong with loess options!")
            return(invisible(NULL))
        }
    }


    ## two graph packages
    if (hist.graphics == "lattice")
    {
        x <- paste(hist.y, "~", hist.x, sep="")

            if ( hist.cond == "")
            {
                x <- as.formula(x)
                lattice.final <- xyplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                main=hist.main , data= hist.data)
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- xyplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }

        return(lattice.final)

     }
     else  ## start ggplot
     {
        hist.type <- switch(hist.type,
                                    p = c("point"),
                                    l = c("line"),
                                    b = c("point", "line"),
                                    #psmooth = c("point", "smooth"),
                                    #lsmooth = c("line", "smooth"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram")
                                    )

        myx <- hist.data[[hist.x]]
        myy <- hist.data[[hist.y]]

        if ( !is.null(hist.cond) && hist.cond != "" )
        {
            # TODO
            if (hist.layout_x == "" || hist.layout_y == "")
            {
               ErrorMessage("layout_x or layout_y is NOT specified for conditional variable")
               return(invisible(NULL))
            }

            mycond <- hist.data[[hist.cond]]
            total.figno <- as.numeric(hist.layout_x) * as.numeric(hist.layout_y)

            if (length(unique(mycond)) > total.figno)
            {
                newno <- unique(mycond)[1:total.figno]
                part.data <- hist.data[which(mycond%in%newno),]
            }
            else part.data <- hist.data

            ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = hist.type, #colour = hist.col,
                                    main=hist.main , data= part.data) + facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))

        }
        else
        {
            ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = hist.type, #colour = hist.col,
                                    main=hist.main , data= hist.data)
        }

        return(ggplot.final)
     }
}

getPKHistCall <- function(hist.graphics, hist.x,
                hist.bin,
                hist.main,
                hist.xlab,
                hist.ylab,
                hist.type,
                hist.cond,
                hist.layout_x,
                hist.layout_y,
                hist.data
                )
{
    hist.col <- getFigConfig()$col
    if (missing(hist.data)) hist.data <- getCurrentData()
    if (missing(hist.type) || hist.type == "") hist.type <- "percent"

    if (hist.graphics == "lattice")
    {
        x <- paste( "~", hist.x, sep="")

        if ( !is.null(hist.bin) && hist.bin != "" )
        {
            if ( hist.cond == "")
            {
                x <- as.formula(x)
                lattice.final <- histogram(x=x, xlab=hist.xlab, nint = as.numeric(hist.bin),
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                main=hist.main , data= hist.data)
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- histogram(x=x, xlab=hist.xlab, nint = as.numeric(hist.bin),
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }
        }
        else
        {
            if ( is.null(hist.cond) || hist.cond =="" )
            {
                x <- as.formula(x)
                lattice.final <- histogram(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type,  col=hist.col,
                                main=hist.main , data= hist.data)
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- histogram(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }

        }

        return(lattice.final)

     }
     else ## start ggplot2
     {
        mytype <- switch(hist.type,
                                    #p = c("point"),
                                    #l = c("line"),
                                    #b = c("point", "smooth"),
                                    #psmooth = c("point", "smooth"),
                                    #lsmooth = c("line", "smooth"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram"),
                                    c("histogram")
                                    )

        myx <- hist.data[[hist.x]]

        if ( !is.null(hist.bin) && hist.bin != "" )
        {
            checkbin <- 1
            f <- as.numeric(diff(range(myx)/as.numeric(hist.bin)))
        }
        else checkbin <- 0

        if ( !is.null(hist.cond) && hist.cond != "" )
        {
            checkcond <- 1
            if (hist.layout_x == "" || hist.layout_y == "")
            {
               ErrorMessage("layout_x or layout_y is NOT specified for conditional variable")
               return(invisible(NULL))
            }

            mycond <- hist.data[[hist.cond]]
            total.figno <- as.numeric(hist.layout_x) * as.numeric(hist.layout_y)

            if (length(unique(mycond)) > total.figno)
            {
                newno <- unique(mycond)[1:total.figno]
                part.data <- hist.data[which(mycond%in%newno),]
            }
            else  part.data <- hist.data

        }
        else checkcond <- 0

        if (checkcond == 1)
        {
            if (checkbin == 1)
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = mytype, #colour = hist.col,
                                    main=hist.main , data= part.data)+ geom_histogram(binwidth = f)+ facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))
            }
            else
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = mytype, #colour = hist.col,
                                    main=hist.main , data= part.data)+ facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))
            }

        }
        else
        {
            if (checkbin == 1)
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = mytype, #colour = hist.col,
                                    main=hist.main , data= hist.data) + geom_histogram(binwidth = f)
            }
            else
            {
                ggplot.final <- qplot(x=myx, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = mytype, #colour = hist.col,
                                    main=hist.main , data= hist.data)
            }
        }

        return(ggplot.final)
     }
}


getPKScatterCall <- function(hist.graphics, hist.x, hist.y,
                hist.bin,
                hist.main,
                hist.xlab,
                hist.ylab,
                hist.type,
                hist.cond,
                hist.layout_x,
                hist.layout_y, hist.data
                )
{
    if (missing(hist.data)) hist.data <- getCurrentData()

    hist.col <- getFigConfig()$col

    ## two graph packages
    if (hist.graphics == "lattice")
    {
        if (hist.type == "")
        {
            if (is.null(getFigConfig()$loess)) hist.type = "p"
            else if (getFigConfig()$loess == 1) hist.type = c("p", "smooth")
            else
            {
                ErrorMessage("Wrong with loess options!")
                return(invisible(NULL))
            }
        }
        
        x <- paste(hist.y, "~", hist.x, sep="")

            if ( hist.cond == "")
            {
                x <- as.formula(x)
                lattice.final <- xyplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                main=hist.main , data= hist.data)
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- xyplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }

        return(lattice.final)

     }
     else  ## start ggplot
     {
        if (hist.type == "")
        {
            if (is.null(getFigConfig()$loess)) hist.type = "point"
            else if (getFigConfig()$loess == 1) hist.type = c("point", "smooth")
            else
            {
                ErrorMessage("Wrong with loess options!")
                return(invisible(NULL))
            }
        }

        myx <- hist.data[[hist.x]]
        myy <- hist.data[[hist.y]]

        if ( !is.null(hist.cond) && hist.cond != "" )
        {
            # TODO
            if (hist.layout_x == "" || hist.layout_y == "")
            {
               ErrorMessage("layout_x or layout_y is NOT specified for conditional variable")
               return(invisible(NULL))
            }

            mycond <- hist.data[[hist.cond]]
            total.figno <- as.numeric(hist.layout_x) * as.numeric(hist.layout_y)

            if (length(unique(mycond)) > total.figno)
            {
                newno <- unique(mycond)[1:total.figno]
                part.data <- hist.data[which(mycond%in%newno),]
            }
            else part.data <- hist.data

            ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = hist.type, #colour = hist.col,
                                    main=hist.main , data= part.data, se=F) + facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))

        }
        else
        {
            ggplot.final <- qplot(x=myx, y=myy, xlab=hist.xlab,
                                    ylab=hist.ylab, geom = hist.type, #colour = hist.col,
                                    main=hist.main , data= hist.data, se=F)
        }


        return(ggplot.final)
     }
}

getPKQqmathCall <- function(hist.graphics, hist.x,
                hist.main,
                hist.xlab,
                hist.ylab,
                hist.type,
                hist.cond,
                hist.layout_x,
                hist.layout_y,
                hist.data
                )
{
    hist.col <- getFigConfig()$col
    if (missing(hist.data)) hist.data <- getCurrentData()

    if (hist.graphics == "lattice")
    {
        x <- paste( "~", hist.x, sep="")

        if (hist.type == "")
        {
            if (is.null(getFigConfig()$loess)) hist.type = "p"
            else if (getFigConfig()$loess == 1) hist.type = c("p", "smooth")
            else
            {
                ErrorMessage("Wrong with loess options!")
                return(invisible(NULL))
            }
        }

            if ( is.null(hist.cond) || hist.cond =="" )
            {
                x <- as.formula(x)
                lattice.final <- qqmath(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type,  col=hist.col,
                                main=hist.main , data= hist.data)
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- qqmath(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }

        return(lattice.final)

     }
     else ## start ggplot2
     {
        ggplot.text <- paste("ggplot(hist.data) + geom_point(aes(sample =", hist.x, "), stat = \"qq\")", sep="")
        ggplot.final <- eval(parse(text=ggplot.text))
        return(ggplot.final)
     }
}
getPKBwplotCall <- function(hist.graphics, hist.x, hist.y,
                hist.bin,
                hist.main,
                hist.xlab,
                hist.ylab,
                hist.type,
                hist.cond,
                hist.layout_x,
                hist.layout_y, hist.data
                )
{
    if (missing(hist.data)) hist.data <- getCurrentData()

    hist.col <- getFigConfig()$col

    ## two graph packages
    if (hist.graphics == "lattice")
    {
            if (is.null(getFigConfig()$loess)) hist.type = "p"
            else if (getFigConfig()$loess == 1) hist.type = c("p", "smooth")
            else
            {
                ErrorMessage("Wrong with loess options!")
                return(invisible(NULL))
            }

        x <- paste(hist.y, "~", hist.x, sep="")

            if ( hist.cond == "")
            {
                x <- as.formula(x)
                lattice.final <- bwplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                main=hist.main , data= hist.data)
            }
            else
            {
                x <-  as.formula(paste(x, "|", hist.cond, sep=" "))
                lattice.layout <- as.numeric(c(hist.layout_x, hist.layout_y))
                lattice.final <- bwplot(x=x, xlab=hist.xlab,
                                ylab=hist.ylab, type= hist.type, col=hist.col,
                                layout = lattice.layout,
                                main=hist.main , data= hist.data)

            }

        return(lattice.final)

     }
     else  ## start ggplot
     {
        if (hist.type == "")
        {
            if (is.null(getFigConfig()$loess)) hist.type = "point"
            else if (getFigConfig()$loess == 1) hist.type = c("point", "smooth")
            else
            {
                ErrorMessage("Wrong with loess options!")
                return(invisible(NULL))
            }
        }

        myx <- hist.data[[hist.x]]
        myy <- hist.data[[hist.y]]

        tmp1 <- paste("ggplot(hist.data, aes(factor(", hist.x, "),", hist.y, sep="")     
        tmp2 <- paste(")) + geom_boxplot() +", "labs(x=hist.xlab, y=hist.ylab)", sep="")
        ggplot.final <- eval(parse(text=paste(tmp1, tmp2, sep="")))

        if ( !is.null(hist.cond) && hist.cond != "" )
        {
            # TODO
            if (hist.layout_x != "" )
                  ggplot.final <- ggplot.final + facet_wrap(hist.cond, ncol = as.numeric(hist.layout_x))
            else if (hist.layout_y != "")
                  ggplot.final <- ggplot.final + facet_wrap(hist.cond, nrow = as.numeric(hist.layout_y))
            else
            {
               ErrorMessage("No layout is specified for conditional variable")
               return(invisible(NULL))
            }

        }

        return(ggplot.final)
     }
}

getPKMatrixplotCall <- function(hist.graphics, hist.data)
{
    currentMain <- svalue(pmg.dialog.notebook)  
    if (missing(hist.data)) hist.data <- getCurrentData(currentMain)  

    hist.col <- getFigConfig()$col

    ## two graph packages
    if (hist.graphics == "lattice")
    {
        lattice.final <- splom(hist.data)
        return(lattice.final)

     }
     else  ## start ggplot
     {
        ggplot.final <- plotmatrix(hist.data)
        return(ggplot.final)
     }
}
################################################################################
################################################################################
#getSubHeight() = 500
#getSubWidth() = 800
#size(pmg.dialog.notebook) <- c(getSubWidth()*0.6,getSubHeight()*.67)
################################################################################
cleanFigureButtonHandler = function(.,h,...)  #1031
{
    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()
    cleanDataSpecialPlot()
    
    for(i in 1:length(pk.dialog.notebook))
      dispose(pk.dialog.notebook)
}

summary.uni.okButtonHandler = function(.,h,...)
{
    tmp.para <- NULL

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          override.closebutton = TRUE)

    ## lattice
    call.command <- "histogram"
    call.final <- getHistCall(hist.graphics = svalue(.$savewd[["graphics"]]), hist.x=svalue(.$widgets[["x"]]),
                              hist.bin = svalue(.$widgets[["number of bins"]]),
                              hist.main=svalue(.$widgets[["main"]]),
                              hist.xlab=svalue(.$widgets[["xlab"]]),
                              hist.ylab=svalue(.$widgets[["ylab"]]),
                              hist.type= svalue(.$widgets[["type"]]),
                              hist.cond = svalue(.$savewd[["cond"]]),
                              hist.layout_x = svalue(.$savewd[["layout_x"]]),
                              hist.layout_y = svalue(.$savewd[["layout_y"]])
                              )
                

        print(call.final)
        setPKCode(list(pkcall=call.command, pklist=call.final))
        setPKGGobi(list(x=svalue(.$widgets[["x"]])))



    #dispose(.$window)
}

summary.uni.ggobiImageHandler = function(.,h,...)   
{
      currentPage <- svalue(pk.dialog.notebook)
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }

      g <- ggobiRun()
      
      mydata <- getPKGGobi(currentPage)
      display(g[1], pmode="Barchart", vars=list(X = mydata$x))
}



################################################################################
summary.bi.okButtonHandler = function(.,h,...)
{
    #tmp.para <- NULL

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                              hist.x=svalue(.$widgets[["x"]]),
                              hist.y=svalue(.$widgets[["y"]]),
                              #hist.bin = svalue(.$widgets[["number of bins"]]),
                              hist.main=svalue(.$widgets[["main"]]),
                              hist.xlab=svalue(.$widgets[["xlab"]]),
                              hist.ylab=svalue(.$widgets[["ylab"]]),
                              hist.type= svalue(.$widgets[["type"]]),
                              hist.cond = svalue(.$savewd[["cond"]]),
                              hist.layout_x = svalue(.$savewd[["layout_x"]]),
                              hist.layout_y = svalue(.$savewd[["layout_y"]])
                              )


    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=svalue(.$widgets[["x"]]), y=svalue(.$widgets[["y"]])))


    #dispose(.$window)
}


summary.bi.ggobiImageHandler = function(.,h,...)  # 0603
{
      currentPage <- svalue(pk.dialog.notebook)
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }

      g <- ggobiRun()

      mydata <- getPKGGobi(currentPage)
      display(g[1], pmode="Scatterplot Display",
              vars=list(X=mydata$x, Y=mydata$y))
}

## 3d does Not have ggplot2 implementation
summary.tri.okButtonHandler = function(.,h,...)
{
    tmp.para <- NULL

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = FALSE)

###################################################

        lattice.call <- "cloud"
        x <- paste(svalue(.$widgets[["z"]]), "~", svalue(.$widgets[["x"]]), "*",
                          svalue(.$widgets[["y"]]), sep="")
        lattice.list <- list(x=x, xlab=svalue(.$widgets[["xlab"]]),
                            ylab=svalue(.$widgets[["ylab"]]), zlab=svalue(.$widgets[["zlab"]]),
                            type= svalue(.$widgets[["type"]]),
                            main=svalue(.$widgets[["main"]]) , data= getCurrentData())

        #if (svalue(.$widgets[["xlim"]])=="")
        if (is.null(.$savewd[["conditional var"]]) ||(svalue(.$savewd[["conditional var"]])==""))
        {
            lattice.list$x <- as.formula(lattice.list$x)
            print(do.call(lattice.call, lattice.list))
        }
        else
        {
            lattice.list$x <-  as.formula(paste(lattice.list$x, "|", svalue(.$savewd[["conditional var"]]), sep=" "))
            lattice.list$layout <- as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]])))
            print(do.call(lattice.call, lattice.list))
        }


}
summary.tri.ggobiImageHandler = function(.,h,...)
{
      x.name <- svalue(.$widgets[["x"]])
      x.ind <- which(colnames(tmp.data)== x.name)
      y.name <- svalue(.$widgets[["y"]])
      y.ind <- which(colnames(tmp.data)== y.name)
      old.ind <- c(1:length(colnames(tmp.data)))
      old.ind <- old.ind[-c(x.ind, y.ind)]

      tmp.data <- tmp.data[c(x.ind, y.ind, old.ind)]

      g <- ggobi(tmp.data)
      #display(g[1], pmode="Barchart", vars=list(X=svalue(.$widgets[[1]])))
      #display(g[1], pmode="Scatterplot Display",
              #vars=list(X=svalue(.$widgets[["x"]]), Y=svalue(.$widgets[["y"]])))
}


summary.para.okButtonHandler = function(.,h,...)
{

    currentMain <- svalue(pmg.dialog.notebook) 

    ## directly apply
    need.var <- svalue(.$widgets[["x"]])
    if (length(need.var) < 1)
    {
        ErrorMessage("Please choose x variables!")
        return(invisible(NULL))
    }
    all.data <- getCurrentData(currentMain)
    part.data <- getCurrentData(currentMain)[,need.var]

    call.command <- "parallel"
    if (is.null(.$savewd[["conditional var"]]) ||(svalue(.$savewd[["conditional var"]])==""))
    {
        x <- as.formula("~part.data")
        call.final <- parallel(x=x, data=all.data,
                            main=svalue(.$widgets[["main"]]) ,
                            horizontal.axis= as.logical(svalue(.$widgets[["horizontal"]])))
    }
    else
    {
         cond <- svalue(.$savewd[["conditional var"]])
         cond.logic <- colnames(Theoph) %in% cond
         if (any(cond.logic))
         {
            cond.ind <- which(cond.logic)
            part.data <- part.data[,-c(cond.ind)]
         }
         x <-  as.formula(paste("~part.data", "|", cond, sep=""))
         call.final <- parallel(x=x, data=all.data,
                            main=svalue(.$widgets[["main"]]) ,
                            horizontal.axis= as.logical(svalue(.$widgets[["horizontal"]])),
                            layout= as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]]))))
    }
    

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          override.closebutton = TRUE)
     print(call.final)
     setPKCode(list(pkcall=call.command, pklist=call.final))
     setPKGGobi(list(x=need.var ))
    
    #print(parallel(~ para.data | factor(cond) , all.data, horizontal.axis= as.logical(svalue(.$widgets[["horizontal"]]))))


}

summary.para.ggobiImageHandler = function(.,h,...)
{
      currentPage <- svalue(pk.dialog.notebook)
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }
      
    g <- ggobiRun()
   
    mydata <- getPKGGobi(currentPage)

    display(g[1], "Parallel Coordinates Display", vars=list(X=mydata$x))
     
}

summary.heat.okButtonHandler = function(.,h,...)
{
    #tmp.para <- NULL
    #for(i in names(.$widgetList))
    #{
      ## store vals in props of super
     # .$.super$props[[i]] <- svalue(.$widgets[[i]]) # pre 0.4-0
     #h$action$super$props[[i]] <- svalue(.$widgets[[i]])
    #}
    #size(.$widgets[["x"]]) <- c(200, 200)
    
    need.var <- svalue(.$widgets[["x"]])
    if (length(need.var) < 1)
    {
        ErrorMessage("Please choose x variables!")
        return(invisible(NULL))
    }
    part.data <- getCurrentData()[,need.var]

    #print(parallel(~ para.data | factor(cond) , all.data, horizontal.axis= as.logical(svalue(.$widgets[["horizontal"]]))))
    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = TRUE)
          
    call.command <- "heatmap"

    #old.par <- par(no.readonly = TRUE)
    t <- par()
    par(t)
    #on.exit(par(old.par))
    
    x <- as.matrix(part.data)
    if (svalue(.$widgets[["dendrogram for row"]]) == "no")
    {
       if (svalue(.$widgets[["dendrogram for column"]]) == "no")
       {
            print(call.final <- heatmap(x, Rowv=NA, Colv=NA, scale= svalue(.$widgets[["scale by"]]), main=svalue(.$widgets[["main"]])))
            setPKCode(list(pkcall=call.command, pklist= call.final))
       }
       else
       {
            print(call.final <- heatmap(x, Rowv=NA, scale= svalue(.$widgets[["scale by"]]), main=svalue(.$widgets[["main"]])))
            setPKCode(list(pkcall=call.command, pklist= call.final))
       }
    }
    else
    {
       if (svalue(.$widgets[["dendrogram for column"]]) == "no")
       {
            print(call.final <- heatmap(x, Colv=NA, scale= svalue(.$widgets[["scale by"]]), main=svalue(.$widgets[["main"]])))
            setPKCode(list(pkcall=call.command, pklist= call.final))
       }
       else
       {
            print(call.final <- heatmap(x, scale= svalue(.$widgets[["scale by"]]), main=svalue(.$widgets[["main"]])))
            setPKCode(list(pkcall=call.command, pklist= call.final))
       }
    }

}

summary.heat.ggobiImageHandler = function(.,h,...)
{
      x.name <- svalue(.$widgets[["x"]])
      x.ind <- which(colnames(tmp.data)== x.name)
      y.name <- svalue(.$widgets[["y"]])
      y.ind <- which(colnames(tmp.data)== y.name)
      old.ind <- c(1:length(colnames(tmp.data)))
      old.ind <- old.ind[-c(x.ind, y.ind)]

      tmp.data <- tmp.data[c(x.ind, y.ind, old.ind)]

      g <- ggobi(tmp.data)
      #display(g[1], pmode="Barchart", vars=list(X=svalue(.$widgets[[1]])))
      #display(g[1], pmode="Scatterplot Display",
              #vars=list(X=svalue(.$widgets[["x"]]), Y=svalue(.$widgets[["y"]])))
}

summary.matrix.okButtonHandler = function(.,h,...)
{

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          override.closebutton = TRUE)

    chose.var <- svalue(.$widgets[["x"]])
    if (length(chose.var) < 1)
    {
        ErrorMessage("Please choose x variables!")
        return(invisible(NULL))
    }

    currentMain <- svalue(pmg.dialog.notebook)
    tmp.data <- getCurrentData(currentMain)
    part.data <- tmp.data[, chose.var]
    call.final <- getPKMatrixplotCall(hist.graphics= svalue(.$savewd[["graphics"]]), part.data)

    call.command <- "splom"
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=chose.var ))
    
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(part.data, as.character(currentPage))

}


summary.matrix.ggobiImageHandler = function(.,h,...)
{
      currentPage <- svalue(pk.dialog.notebook)
      if (currentPage == 0)
      {
          ErrorMessage("Please draw figure first!")
          return(NULL)
      }

    g <- ggobiRun()
   
    mydata <- getPKGGobi(currentPage)
    display(g[1], "Scatterplot Matrix", vars=list(X=mydata$x))
}
################################################################################
################################################################################
model.ind.okButtonHandler = function(.,h,...)
{
    # ind plot
    match.term <- getTerm()
    hist.cond <- match.term[match.term$TermName == "ID",]$VarName

    if (svalue(.$widgets[["main"]])=="")
        newstr <- paste(svalue(.$widgets[["y"]]), "vs", svalue(.$widgets[["x"]]), sep="")
    else newstr <- svalue(.$widgets[["main"]])

    pgraph = ggraphics(ps=6)
    add(pk.dialog.notebook, pgraph, label = newstr,
              override.closebutton = TRUE)

    if (svalue(.$widgets[["x"]]) == "")
    {
        ErrorMessage("You have to choose one x value!")
        return(invisible(NULL))
    }
    
    if (svalue(.$widgets[["y"]]) != "")
    {

      call.command <- "xyplot"
      call.final <- getScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=svalue(.$widgets[["x"]]),
                                hist.y=svalue(.$widgets[["y"]]),
                                hist.main=svalue(.$widgets[["main"]]),
                                hist.xlab=svalue(.$widgets[["xlab"]]),
                                hist.ylab=svalue(.$widgets[["ylab"]]),
                                hist.type= svalue(.$widgets[["type"]]),
                                hist.cond = hist.cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]])
                                )
    }
    else
    {
      call.command <- "histogram"
      call.final <- getHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=svalue(.$widgets[["x"]]),
                                hist.bin = "",
                                hist.main=svalue(.$widgets[["main"]]),
                                hist.xlab=svalue(.$widgets[["xlab"]]),
                                hist.ylab=svalue(.$widgets[["ylab"]]),
                                hist.type= svalue(.$widgets[["type"]]),
                                hist.cond = hist.cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]])
                                )
    }


    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=svalue(.$widgets[["x"]]), y=svalue(.$widgets[["y"]]) ))

    
}


###############################################################################
model.gof.okButtonHandler = function(.,h,...)
{
    match.term <- getTerm()
    mywres <- match.term[match.term$TermName == "WRES",]$VarName
    mypred <- match.term[match.term$TermName == "PRED",]$VarName
    myipre <- match.term[match.term$TermName == "IPRE",]$VarName
    mydv <- match.term[match.term$TermName == "DV",]$VarName
    myidv <- match.term[match.term$TermName == "IDV",]$VarName

    pkdata <- getCurrentData()
    mypackage <- svalue(.$savewd[["graphics"]])

    # DV vs PRED
    mylist <- list()
    myx <- svalue(.$widgets[["PRED_1"]])
    myy <- svalue(.$widgets[["DV_1"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    #print(do.call(mycall, mylist))
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))

    # DV vs IPRED:
    mylist <- list()
    myx <- svalue(.$widgets[["IPRE"]])
    myy <- svalue(.$widgets[["DV_2"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))

    # WRES vs IDV
    mylist <- list()
    myx <- svalue(.$widgets[["IDV_3"]])
    myy <- svalue(.$widgets[["WRES_3"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))

    # PRED vs IDV:
    mylist <- list()
    myx <- svalue(.$widgets[["IDV_4"]])
    myy <- svalue(.$widgets[["PRED_4"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # IPRED vs IDV
    mylist <- list()
    myx <- svalue(.$widgets[["IDV_5"]])
    myy <- svalue(.$widgets[["IPRE_5"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
}

model.struct.okButtonHandler = function(.,h,...)
{
    match.term <- getTerm()
    mywres <- match.term[match.term$TermName == "WRES",]$VarName
    mypred <- match.term[match.term$TermName == "PRED",]$VarName
    mycov <- match.term[match.term$TermName == "COV",]$VarName
    myipre <- match.term[match.term$TermName == "IPRE",]$VarName

    pkdata <- getCurrentData()
    mypackage <- svalue(.$savewd[["graphics"]])

    # PRED vs DV|IDV
    mylist <- list()
    myx <- svalue(.$widgets[["DV_1"]])
    myy <- svalue(.$widgets[["PRED_1"]])
    cond <- svalue(.$widgets[["IDV_1"]])
    part.data <- pkdata[c(myx,myy,cond)]
    
    newstr <- paste(myy, "vs", myx, "|", cond, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main=newstr,
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = part.data
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # IPRED vs DV|IDV:
    mylist <- list()
    myx <- svalue(.$widgets[["DV_2"]])
    myy <- svalue(.$widgets[["IPRE"]])
    cond <- svalue(.$widgets[["IDV_2"]])
    part.data <- pkdata[c(myx,myy,cond)]

    newstr <- paste(myy, "vs", myx, "|", cond, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main=newstr,
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = part.data
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # WRES vs IDV:
    mylist <- list()
    myx <- svalue(.$widgets[["IDV_3"]])
    myy <- svalue(.$widgets[["WRES_3"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # WRES vs IDV (bw):
    mylist <- list()
    myx <- svalue(.$widgets[["IDV_4"]])
    myy <- svalue(.$widgets[["WRES_4"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "bwplot"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              override.closebutton = TRUE)

    pkdata[[myx]] <- factor(pkdata[[myx]]) 
    call.command <- "bwplot"
    call.final <- getPKBwplotCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # WRES vs PRED:
    mylist <- list()
    myx <- svalue(.$widgets[["PRED_5"]])
    myy <- svalue(.$widgets[["WRES_5"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # TODO: figure does NOT make sense
    # WRES vs PRED (bw):
    mylist <- list()
    myx <- svalue(.$widgets[["PRED_6"]])
    myy <- svalue(.$widgets[["WRES_6"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "bwplot"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    pkdata[[myx]] <- factor(pkdata[[myx]])
    call.command <- "bwplot"
    call.final <- getPKBwplotCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # PRED vs DV|Covariates
    mylist <- list()
    myx <- svalue(.$widgets[["DV_7"]])
    myy <- svalue(.$widgets[["PRED_7"]])
    cond <- svalue(.$widgets[["COV_7"]])

    newstr <- paste(myy, "vs", myx, "|", cond, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main=newstr,
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    # IPRED vs DV|Covariates
    mylist <- list()
    myx <- svalue(.$widgets[["DV_8"]])
    myy <- svalue(.$widgets[["IPRE_8"]])
    cond <- svalue(.$widgets[["COV_8"]])

    newstr <- paste(myy, "vs", myx, "|", cond, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
}


model.resid.okButtonHandler = function(.,h,...)
{

    match.term <- getTerm()
    mywres <- match.term[match.term$TermName == "WRES",]$VarName
    mypred <- match.term[match.term$TermName == "PRED",]$VarName
    mycov <- match.term[match.term$TermName == "COV",]$VarName
    myipre <- match.term[match.term$TermName == "IPRE",]$VarName

    pkdata <- getCurrentData()
    mypackage <- svalue(.$savewd[["graphics"]])

    # Distribution of WRES:
    mylist <- list()
    myx <- svalue(.$widgets[["Distribution of WRES:"]])

    newstr <- paste("Distribution of", myx, sep=" ")
    plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "histogram"
    call.final <- getPKHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.bin = "",
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL))
    
    # Distribution of WRES (QQ):
    mylist <- list()
    myx <- svalue(.$widgets[["Distribution of WRES(QQ):"]])

    newstr <- paste("Distribution of", myx, "(QQ)", sep=" ")
    #plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "qqmath"
    call.final <- getPKQqmathCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL ))

        
    ## TODO: only work for lattice
    # Individual distribution of WRES:
    mylist <- list()
    myx <- svalue(.$widgets[[" of WRES:"]])

    newstr <- paste("Individual distribution of", myx, sep=" ")
    plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    cond <- match.term[match.term$TermName == "ID",]$VarName
    
    pkdata[[cond]] <- factor(pkdata[[cond]])

    call.command <- "histogram"
    call.final <- getPKHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.bin = "",
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond =cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    
    currentPage <- svalue(pk.dialog.notebook)
    setDataLayoutPlot(as.character(currentPage))
    setPKGGobi(list(x=myx, y=NULL))
    
    # Individual distribution of WRES (QQ):
    mylist <- list()
    myx <- svalue(.$widgets[[" of WRES(QQ)"]])

    newstr <- paste("Distribution of", myx, "(QQ)", sep=" ")
    #plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    ## TODO: make ggplot for qq works too
        
        call.command <- "qqmath"
        call.final <- getPKQqmathCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

        print(call.final)
        setPKCode(list(pkcall=call.command, pklist=call.final))
        currentPage <- svalue(pk.dialog.notebook)
        setDataLayoutPlot(as.character(currentPage))
        setPKGGobi(list(x=myx, y=NULL ))
        
    ## |WRES| vs PRED
    mylist <- list()
    myx <- svalue(.$widgets[["PRED_1"]])
    myy <- svalue(.$widgets[["|WRES|_1"]])
    ## get |WRES|
    part.data <- pkdata[c(myx,myy)]
    part.data[,2] <- abs(part.data[,2])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab="|WRES|",
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = part.data
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(part.data, as.character(currentPage))
    setPKGGobi(list(x=myx, y=myy ))
    
    # Covariates vs |WRES| (bw):
    mylist <- list()
    myx <- svalue(.$widgets[["Covariates_2"]])
    myy <- svalue(.$widgets[["|WRES|_2"]])
    
    # get |WRES|
    part.data <- pkdata[c(myx,myy)]
    part.data[,1] <- abs(part.data[,1])
    #myx <- paste("|",myx,"|", sep="")
    colnames(part.data) <- c(myx, myy)

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "bwplot"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    part.data[[myx]] <- factor(part.data[[myx]])   
    call.command <- "bwplot"
    call.final <- getPKBwplotCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab="|WRES|",
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = part.data
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(part.data, as.character(currentPage))
    setPKGGobi(list(x=myx, y=myy ))
    
    # |WRES| vs PRED|Covariates
    mylist <- list()
    myx <- svalue(.$widgets[["PRED_3"]])
    myy <- svalue(.$widgets[["|WRES|_3"]])
    cond <- svalue(.$widgets[["Covariates_3"]])

    part.data <- pkdata[c(myx, myy, cond)]
    part.data[,2] <- abs(part.data[,2])

    newstr <- paste("|", myy, "|", "vs", myx, "|", cond, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main=newstr,
                                hist.xlab=myx,
                                hist.ylab="|WRES|",
                                hist.type= "",
                                hist.cond = cond,
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = part.data
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(part.data, as.character(currentPage))
    setDataLayoutPlot(as.character(currentPage))
    setPKGGobi(list(x=myx, y=myy ))
    
    #TODO:
    # |IWRES| vs IPRED|Covariates
    
    # Autocorrelation of WRES: Plot of WRESi against WRESi+1
    mylist <- list()
    mywres <- svalue(.$widgets[["Autocorrelation of WRES:"]])
    tmp <- pkdata[[mywres]]
    mydf <- data.frame(WRESi=tmp[-length(tmp)], WRESi_1=tmp[-1])

    myx <- "WRESi"
    myy <- "WRESi_1"
    
    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = mydf
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(mydf, as.character(currentPage))

    setPKGGobi(list(x=myx, y=myy ))
    
}

model.para.okButtonHandler = function(.,h,...)
{

    match.term <- getTerm()
    mypara <- match.term[match.term$TermName == "PARAMETERS",]$VarName

    pkdata <- getCurrentData()
    mypackage <- svalue(.$savewd[["graphics"]])

    # Distribution of parameters:
    mylist <- list()
    myx <- svalue(.$widgets[["Distribution of parameters:"]])

    newstr <- paste("Distribution of", myx, sep=" ")
    plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              override.closebutton = TRUE)

    call.command <- "histogram"
    call.final <- getPKHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.bin = "",
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL))
    
    # Distribution of parameters (QQ):
    mylist <- list()
    myx <- svalue(.$widgets[["Distribution of parameters (QQ):"]])
    
    newstr <- paste("Distribution of", myx, "(QQ)", sep=" ")
    #plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    ## TODO: make ggplot for qq works too

    call.command <- "qqmath"
    call.final <- getPKQqmathCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL))
        
    # Scatterplot matrix of parameters
    newstr <- "Scatterplot matrix of parameters"
    plotType <- "smatrix"

        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    tmp.data <- pkdata[c(mypara)]

    ## TODO
    if (ncol(tmp.data) > 1)
    {
        call.command <- "matrix"
        call.final <- getPKMatrixplotCall(hist.graphics= svalue(.$savewd[["graphics"]]), tmp.data)
        print(call.final)
        setPKCode(list(pkcall=call.command, pklist=call.final))
        setPKGGobi(list(x=colnames(tmp.data), y=NULL))
    }

    # Parameter vs parameter:
    mylist <- list()
    myx <- svalue(.$widgets[["Parameters_x:"]])
    myy <- svalue(.$widgets[["Parameters_y:"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
}

model.cov.okButtonHandler = function(.,h,...)
{
    mylist <- list()
    match.term <- getTerm()
    myeta <- match.term[match.term$TermName == "ETA",]$VarName
    mycov <- match.term[match.term$TermName == "COV",]$VarName
    mypara <- match.term[match.term$TermName == "PARAMETERS",]$VarName

    pkdata <- getCurrentData()
    mypackage <- svalue(.$savewd[["graphics"]])

    # Scatterplot matrix of covariates
    newstr <- "Scatterplot matrix of covariates"
    plotType <- "smatrix"
    
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)
            
    tmp.data <- pkdata[c(mycov)]

    ## TODO
    if (ncol(tmp.data) > 1)
    {
        call.command <- "matrix"
        call.final <- getPKMatrixplotCall(hist.graphics= svalue(.$savewd[["graphics"]]), tmp.data)
        print(call.final)
        setPKCode(list(pkcall=call.command, pklist=call.final))
        setPKGGobi(list(x=colnames(tmp.data), y=NULL)) 

    }
    
    #Parameters vs covariates:
    mylist <- list()
    myx <- svalue(.$widgets[["Cov_P:"]])
    myy <- svalue(.$widgets[["Parameters:"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    #ETAs vs covariates::
    mylist <- list()
    myx <- svalue(.$widgets[["Cov_E:"]])
    myy <- svalue(.$widgets[["ETAS:"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)


    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))
    
    #WRES vs covariates:
    mylist <- list()
    myx <- svalue(.$widgets[["Cov_W:"]])
    myy <- svalue(.$widgets[["WRES:"]])

    newstr <- paste(myy, "vs", myx, sep=" ")
    plotType <- "scatter"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

    call.command <- "xyplot"
    call.final <- getPKScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.y=myy,
                                hist.main="",
                                hist.xlab=myx,
                                hist.ylab=myy,
                                hist.type= "",
                                hist.cond = "",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=myy ))

}

###############################################################################
model.random.okButtonHandler = function(.,h,...)
{

    pkdata <- getCurrentData()
    mylist <- list()
    match.term <- getTerm()
    myterm <- match.term[match.term$TermName == "ETA",]$VarName
    mypackage <- svalue(.$savewd[["graphics"]])
    
    #Distribution of ETAS:
    newstr <- "Distribution of ETAS"
    plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)
    myx <- svalue(.$widgets[["Distribution of ETAS"]])
    
    call.command <- "histogram"
    call.final <- getPKHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.bin = "",
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL))
    
    #Distribution of ETAs (QQ):
    newstr <- "Distribution of ETAS(QQ)"
    #plotType <- "hist"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              override.closebutton = TRUE)

        mylist <- list()
        myx <- svalue(.$widgets[["Distribution of ETAs (QQ)"]])
        
    call.command <- "qqmath"
    call.final <- getPKQqmathCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.main=newstr,
                                hist.xlab="",
                                hist.ylab="",
                                hist.type= "",
                                hist.cond ="",
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = pkdata
                                )

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL))
     
    #Scatterplot matrix of ETAs:
    newstr <- "Scatterplot matrix of ETAs"
    plotType <- "smatrix"

    tmp.data <- pkdata[c(myterm)]

    # TODO
    currentPage <- svalue(pk.dialog.notebook)

    if (ncol(tmp.data) > 1)
    {
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = newstr,
              #pageno = 3,
              override.closebutton = TRUE)

        call.command <- "matrix"
        call.final <- getPKMatrixplotCall(hist.graphics= svalue(.$savewd[["graphics"]]), tmp.data)
        print(call.final)
        setPKCode(list(pkcall=call.command, pklist=call.final))
        setPKGGobi(list(x=colnames(tmp.data), y=NULL)) 
        
    }
    
}


################################################################################
## Menu "Model validation"
################################################################################
 
psn.outlier.okButtonHandler = function(.,h,...)
{
    file1 <- svalue(.$widgets[["Result file:"]])
    file2 <- svalue(.$widgets[["Deleted ID file:"]])

    if (file1!="" && file2!="")
    {
        checkdata <- psn.cdd(file1, file2)
        if (is.null(checkdata)) 
        {
            ErrorMessage("Two files are not right format for PsN!")
            return(invisible(NULL))
        }
        
        mydata <- data.frame(checkdata)    

        ## TODO: what does this function for?
        #setValidateData(mydata)

        message <- "PsN summary for influence analysis"
        pgraph = ggraphics(ps=6)
        add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)
        
            if (svalue(.$savewd[["graphics"]]) == "lattice")
            {
                call.final <- xyplot(cov.ratios~cook.scores, data= mydata,
                                        groups=ID,
                                        type="l",
                                        xlab="Cook.scores",
                                        ylab="Cov.ratios",
                                        main="",
                                        panel= function(x, y,groups, subscripts, ...)
                                        {
                                            panel.xyplot(x, y, groups = groups, subscripts = subscripts, ...)
                                            panel.text(x, y, groups[subscripts])
                                        }
                                      )
            }
            else
            {
                call.final <- qplot(cook.scores, cov.ratios, data = mydata, label = ID,
                                geom=c("point", "text")) + xlab("Cook.scores") + ylab("Cov.ratios")
            }



                print(call.final)
                call.command <- "xyplot"
                setPKCode(list(pkcall=call.command, pklist=call.final))
                currentPage <- svalue(pk.dialog.notebook)
                setDataSpecialPlot(mydata, as.character(currentPage))
    }
    else
        ErrorMessage("Please choose proper files first!")

}

vis.outlier.okButtonHandler = function(.,h,...)
{

    target.dirpath <- svalue(.$widgets[["Target directory path:"]])
    sim.pattern <- svalue(.$widgets[["Simulation folder pattern:"]])
    startFileName <- svalue(.$widgets[["NONMEM result file name:"]])
    cond.var <- svalue(.$widgets[["Plot variable:"]])
    id.var <- svalue(.$widgets[["Patient ID:"]])
    rowLabel <- "resampleID"

    if (target.dirpath=="" || sim.pattern=="" || startFileName=="")
    {
        ErrorMessage("Please input all parameters first!")
        return(invisible(NULL))
    }    

    # match the current data set
    currentPage <- svalue(pmg.dialog.notebook)
    total.id <- unique(getCurrentData(currentPage)[[id.var]])
    final.df <- data.frame(ID=total.id)

    if (cond.var==id.var)
    {
        ErrorMessage("You can NOT have same variable for both Patient ID and Plot variable!")
        return(invisible(NULL))
    }

    if (target.dirpath!="")
    {
    
            ## start process bar
            convertW = gwindow(title="Processing...", parent=c(150,200), height=getSubHeight()*0.05, width=getSubWidth()*0.5, horizontal=FALSE)
            convert.group =ggroup(horizontal=FALSE, spacing=0, expand=TRUE)
            convert.bar <- gtkProgressBar()
      
            add(convert.group, convert.bar)
            add(convertW, convert.group)
      
            # start with something
            gtkProgressBarSetFraction(convert.bar, 0.2)
      
            cdd.list <- extractCddData(target.dirpath, sim.pattern, startFileName, id.var, cond.var, total.id, rowLabel)
            if (is.null(cdd.list))
            {
                dispose(convertW)
                return(invisible(NULL))
            }
                        
            gtkProgressBarSetFraction(convert.bar, 0.7)
            
            extract.data <- cdd.list$data
            extract.data <- data.frame(extract.data)
            v.delete.id <- cdd.list$deleteID
            final.para <- extract.data

            ####################################
            ## plot
            # plot 1 ---------- Patient ID deletion plot
            if (length(v.delete.id) == length(total.id))
            {
                mydata1 <- data.frame(deleteID=v.delete.id, runID=total.id)

                ## parallel plot
                message <- "Patient ID deletion plot"
                pgraph = ggraphics(ps=6)
                add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)

                call.command <- "xyplot"
                myy <- "deleteID"
                myx <- "runID"

                if (svalue(.$savewd[["graphics"]]) == "lattice")
                {
                    call.final <- xyplot(deleteID~runID, main="Patient ID deletion plot for case deletion diagnostics",
                                        xlab= "Case deletion run ID", ylab="Deleted ID", type=c("p", "l"),
                                        data=mydata1)
                }
                else
                {
                    call.final <- qplot(runID, deleteID, main="Patient ID deletion plot for case deletion diagnostics",
                                  xlab= "Case deletion run ID", ylab="Deleted ID",
                                  geom=c("point", "line"),data=mydata1)
                }
                print(call.final)
                setPKCode(list(pkcall=call.command, pklist=call.final))
                currentPage <- svalue(pk.dialog.notebook)
                setDataSpecialPlot(mydata1, as.character(currentPage))
                # setGGobi
            }

            # plot 2 ---------- histogram
            boot.all <- NULL
            resampleID.all <- NULL
            id.all <- NULL

            sapply(1:ncol(extract.data), function(i)
                  {
                     #miss.id <- which(! id.unique %in% boot.data[,i])
                     boot.all <<- c(boot.all, extract.data[,i])
                     resampleID.all <<- c(resampleID.all, rep(i, nrow(extract.data)))
                     id.all <<- c(id.all, rownames(extract.data))
                     invisible(NULL)
                  })

            boot.df <- data.frame(ID=id.all, resampleID=resampleID.all, para=boot.all)

            message <- "Grouped by patient ID"
            pgraph = ggraphics(ps=6)
            add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)

            if (svalue(.$savewd[["graphics"]]) == "lattice")
            {
                call.final <- densityplot(~para, groups=ID,
                              data=boot.df, main="grouped by patient ID", xlab=cond.var)
            }
            else
            {
                call.final <- ggplot(boot.df, aes(para)) + stat_density(geom = "path",
                      position = "identity", aes(colour = factor(ID)))
            }
            print(call.final)
            call.command <- "densityplot"
            setPKCode(list(pkcall=call.command, pklist=call.final))
            currentPage <- svalue(pk.dialog.notebook)
            setDataSpecialPlot(boot.df, as.character(currentPage))

            message <- "Grouped by case deletion run ID"
            pgraph = ggraphics(ps=6)
            add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)

            if (svalue(.$savewd[["graphics"]]) == "lattice")
            {
                call.final <- densityplot(~para, groups=resampleID, data=boot.df,
                              main="grouped by case deletion run ID", xlab=cond.var)
            }
            else
            {
                call.final <- ggplot(boot.df, aes(para)) + stat_density(geom = "path",
                      position = "identity", aes(colour = factor(resampleID))) + xlab(cond.var)
            }
            print(call.final)
            call.command <- "densityplot"
            setPKCode(list(pkcall=call.command, pklist=call.final))
            currentPage <- svalue(pk.dialog.notebook)
            setDataSpecialPlot(boot.df, as.character(currentPage))

            gtkProgressBarSetFraction(convert.bar, 0.8)

           # plot 3 ---------- parallel coor plot
            message <- "Parallel coordinate plot"
            pgraph = ggraphics(ps=6)
            add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)
            scale.data <- rbind(extract.data, min(extract.data, na.rm=T), max(extract.data, na.rm=T))



            call.final <- parallel(~scale.data, main="bounded by global min and max", ylab="Case deletion ID")
            call.command <- "parallel"
            setPKCode(list(pkcall=call.command, pklist=call.final))
            currentPage <- svalue(pk.dialog.notebook)
            setDataSpecialPlot(scale.data, as.character(currentPage))

            if (svalue(.$savewd[["graphics"]]) != "lattice")
            {
                scale.data <- namerows(scale.data, col.name = "ID")
                df <- melt(scale.data[-10], id.var = c("ID"))
                dfm <- ddply(df, .(variable), transform, rng = rescaler(value,
                     type = "range"))

                call.final <- ggplot(dfm, aes(group = ID, colour = factor(ID))) +
                     geom_line(aes(variable, rng)) + xlab("Case deletion ID")+
                    coord_flip()

            }
            print(call.final)

            gtkProgressBarSetFraction(convert.bar, 0.9)
            
           # plot 4 ---------- MDS plot
            e2.data <- t(extract.data)
            #colnames(e2.data) <- 1:ncol(e2.data)
            #parallel(~e2.data)

            mydata <- e2.data
            d <- dist(mydata) # euclidean distances between the rows
            fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
            fit # view results

            # plot solution
            x <- fit$points[,1]
            y <- fit$points[,2]
            message <- "Metric MDS"
            pgraph = ggraphics(ps=6)
            add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)
            #call.final <- xyplot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
            #  main="Metric MDS", type="n")
            #print(call.final)
            #call.text <- text(x, y, labels = row.names(mydata), cex=.7)
            #print(call.text)

            #NOTE: NEED TO UPDATE -- "SIM" with proper one
            mydata2 <- data.frame(cor1=x, cor2=y, resampleID=gsub(rowLabel, "", names(x)))

            if (svalue(.$savewd[["graphics"]]) == "lattice")
            {
                call.final <- xyplot(cor2~cor1, data= mydata2,
                                        groups=resampleID,
                                        type="l",
                                        xlab="Coordinate 1",
                                        ylab="Coordinate 2",
                                        main="Metric MDS",
                                        panel= function(x, y,groups, subscripts, ...)
                                        {
                                            panel.xyplot(x, y, groups = groups, subscripts = subscripts, ...)
                                            panel.text(x, y, groups[subscripts])
                                        }
                                      )
            }
            else
            {
                call.final <- qplot(cor1, cor2, data = mydata2, label = resampleID,
                                geom=c("point", "text")) + xlab("Coordinate 1") + ylab("Coordinate 2")
            }


                print(call.final)
                call.command <- "xyplot"
                setPKCode(list(pkcall=call.command, pklist=call.final))
                currentPage <- svalue(pk.dialog.notebook)
                setDataSpecialPlot(mydata2, as.character(currentPage))

         ## interactive graphics data
         ig.data <- merge(boot.df, mydata2, by="resampleID")
         ig.data <- ig.data[,c("resampleID", "para", "ID", "cor1", "cor2")]
         ig.data$resampleID <- factor(ig.data$resampleID)

         currentPage <- length(getNameDataSpecialPlot()) + 1
         setDataSpecialPlot(ig.data, as.character(currentPage))
        
         gtkProgressBarSetFraction(convert.bar, 1.0)
         dispose(convertW)
      }

}


psn.bootstrap.vis.okButtonHandler = function(.,h,...)
{

    target.dirpath <- svalue(.$widgets[["Target directory path:"]])
    sim.pattern <- svalue(.$widgets[["Bootstrap folder pattern:"]])
    startFileName <- svalue(.$widgets[["NONMEM result file name:"]])
    
    boot.key.path <- svalue(.$widgets[["Bootstrap key table path:"]])
    boot.key.name <- svalue(.$widgets[["Bootstrap key table name:"]])
    
    cond.var <- svalue(.$widgets[["Plot variable:"]])
    id.var <- svalue(.$widgets[["Patient ID:"]])
    rowLabel <- "resampleID"

    if (target.dirpath=="" || sim.pattern=="" || startFileName=="" || boot.key.path=="" || boot.key.name=="")
    {
        ErrorMessage("Please input all parameters first!")
        return(invisible(NULL))
    }

    # match the current data set
    currentPage <- svalue(pmg.dialog.notebook)
    total.id <- unique(getCurrentData(currentPage)[[id.var]])
    final.df <- data.frame(ID=total.id)
    

    if (cond.var==id.var)
    {
        ErrorMessage("You can NOT have same variable for both Patient ID and Plot variable!")
        return(invisible(NULL))
    }

    if (Sys.info()[["sysname"]] == "Windows")
    {
    	check.bootkey <- try(bootKey.table <- read.csv(paste(boot.key.path, boot.key.name, sep="\\"), header=F))
    	
    }
    else
    {
    	check.bootkey <- try(bootKey.table <- read.csv(paste(boot.key.path, boot.key.name, sep="/"), header=F))
    }
    
    if (inherits(check.bootkey, "try-error"))
    {
       ErrorMessage("Bootstrap key table can NOT be read!")
       return(invisible(NULL))
    }    

## start process bar
      convertW = gwindow(title="Processing...", parent=c(150,200), height=getSubHeight()*0.05, width=getSubWidth()*0.5, horizontal=FALSE)
      convert.group =ggroup(horizontal=FALSE, spacing=0, expand=TRUE)
      convert.bar <- gtkProgressBar()

      add(convert.group, convert.bar)
      add(convertW, convert.group)

      # start with something
      gtkProgressBarSetFraction(convert.bar, 0.2)

    ori.data <- extractBootData(target.dirpath, sim.pattern, startFileName, id.var, cond.var, bootKey.table, total.id, 1)
    if (is.null(ori.data))
    {
        dispose(convertW)
        return(invisible(NULL))
    }
    
    extract.data <- data.frame(ori.data)
    gtkProgressBarSetFraction(convert.bar, 0.7)

    # dat prepare for bootstrap number
    boot.all <- NULL
    resampleID.all <- NULL
    id.all <- NULL

    sapply(1:ncol(extract.data), function(i)
          {
             #miss.id <- which(! id.unique %in% boot.data[,i])
             boot.all <<- c(boot.all, extract.data[,i])
             resampleID.all <<- c(resampleID.all, rep(i, nrow(extract.data)))
             id.all <<- c(id.all, rownames(extract.data))
             invisible(NULL)
          })

    boot.df <- data.frame(ID=id.all, resampleID=resampleID.all, para=boot.all)

    ####################################
    ## plot
    # plot 1 ---------- Patient ID deletion plot
    
    # ---------- Bootstrap randomization - patient ID deletion plot
    deleted.df <- unique(boot.df[is.na(boot.df$para),][,c(1,2)])
    choose.df <- unique(boot.df[!is.na(boot.df$para),][,c(1,2)])

    message <- "Resampling design"
    pgraph = ggraphics(ps=6)
    add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)

    if (svalue(.$savewd[["graphics"]]) == "lattice")
    {
        call.final <- xyplot(resampleID~ID, data=choose.df, type="p", xlab="Bootstrap ID",
                      main="Bootstrap randomization check")
    }
    else
    {
        call.final <- qplot(ID, resampleID, data=choose.df, geom=c("point"), xlab="Bootstrap ID",
                      main="Bootstrap randomization check")
    }
    
    print(call.final)
    call.command <- "xyplot"
    setPKCode(list(pkcall=call.command, pklist=call.final))
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(choose.df, as.character(currentPage))
    gtkProgressBarSetFraction(convert.bar, 0.8)
    
    message <- "grouped by bootstrap run ID"
    pgraph = ggraphics(ps=6)
    add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)
                      
    if (svalue(.$savewd[["graphics"]]) == "lattice")
    {
        call.final <- densityplot(~para, groups=resampleID, data=boot.df,
                      main="grouped by bootstrap run ID", xlab=cond.var )
    }
    else
    {
        call.final <- ggplot(boot.df, aes(para)) + stat_density(geom = "path",
                      position = "identity", aes(colour = factor(resampleID))) + xlab(cond.var)
    }
    print(call.final)
    call.command <- "densityplot"
    setPKCode(list(pkcall=call.command, pklist=call.final))
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(boot.df, as.character(currentPage))
    gtkProgressBarSetFraction(convert.bar, 0.9)
    
    ## plot n ------------- rank boostrap variability
    extract.var <- NULL
    sapply(1:nrow(extract.data), function(i)
          {
              extract.var <<- c(extract.var, var(unlist(extract.data[i,]), na.rm=T))
              invisible(NULL)
          })

    var.df <- data.frame(ID=order(total.id), VAR=extract.var)
    var.plot <- var.df[order(var.df$VAR),]
    var.plot$ID <- factor(var.plot$ID, levels= var.plot$ID)
    
    message <- "variability of parameters for ordered ID"
    pgraph = ggraphics(ps=6)
    add(pk.dialog.notebook, pgraph, label = message,
                      override.closebutton = TRUE)
    if (svalue(.$savewd[["graphics"]]) == "lattice")
    {
        call.final <- xyplot(VAR~ID, data=var.plot, xlab= "Ordered ID" ,
                             ylab= paste("Variance of ", cond.var, sep=""))
    }
    else
    {
        call.final <- qplot(ID, VAR, data=var.plot,xlab= "Ordered ID" ,
                             ylab= paste("Variance of ", cond.var, sep=""))
    }
    print(call.final)
    call.command <- "xyplot"
    setPKCode(list(pkcall=call.command, pklist=call.final))
    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(var.plot, as.character(currentPage))

    currentPage <- length(getNameDataSpecialPlot()) + 1
    setDataSpecialPlot(var.plot, as.character(currentPage))

    #currentPage <- as.numeric(svalue(pk.dialog.notebook))
    #sapply(1:currentPage, function(i)
          #{
          #    setDataSpecialPlot(var.plot, as.character(i))
          #})
    gtkProgressBarSetFraction(convert.bar, 1.0)
    dispose(convertW)    
}

## Code modified from PsN
psn.bootstrap.sum.okButtonHandler = function(.,h,...)
{
   while (svalue(pk.dialog.notebook) > 1)
   {
      dispose(pk.dialog.notebook)
   }

    file1 <- svalue(.$widgets[["PsN result file:"]])
    file2 <- svalue(.$widgets[["Bootstrap key file:"]])

    min.failed    <- FALSE      # do we want to omit minimization failed runs?
    cov.failed    <- FALSE      # do we want to omit covariance failed runs?
    cov.warnings  <- TRUE       # do we want to omit covariance failed runs?
    boundary      <- TRUE       # do we want to omit boundary runs?
    showoriginal  <- TRUE       # show line for original estimate
    showmean      <- TRUE       # show line for mean
    showmedian    <- FALSE      # show line for median
    show95CI      <- TRUE       # show line for 95 % confidence interval (percentile)
    showquart     <- FALSE      # show line for quartiles

    excl.id <- c()              # exclude samples that have this individual

    ## read files
    b.try1 <- try(bootstrap.data <- read.csv(file1, header=T)) # read.csv("raw_results1.csv", header=T)
    in.try2 <- try(incl.ids <- read.csv(file2, header=F)) # read.csv("included_individuals1.csv", header=F)

    if (inherits(b.try1, "try-error") || inherits(in.try2, "try-error"))
    {
        ErrorMessage("PsN result file or key file can NOT be read!")
        return(invisible(NULL))
    }

    check.col <- c("minimization_successful", "covariance_step_successful",
                   "covariance_step_warnings", "estimate_near_boundary")
    if (!all(check.col %in% names(bootstrap.data))) 
    {
        ErrorMessage("PsN result file or key file can NOT be read!")
        return(invisible(NULL))
    }
    
    ## replace underscores
    for (i in 1:length(names(bootstrap.data))) {
      names(bootstrap.data)[i] <- gsub("_", ".", names(bootstrap.data)[i])
    }

    ## find ofv column index
    index <- 0
    seen  <- FALSE

    for (i in names(bootstrap.data)) {
      if (!seen) {
        index <- index + 1
      }
      if (i == "ofv") {
        seen <- TRUE
      }
    }

    ## get number of parameters
    n       <- length(colnames(bootstrap.data)) - index
    nparams <- length(colnames(bootstrap.data))

    ## separate out original model fit
    p1 <- subset(bootstrap.data, bootstrap.data$model != 0)
    o1 <- subset(bootstrap.data, bootstrap.data$model == 0)

    incl.flag <- rep(0,length(rownames(p1)))
    for( i in excl.id ) {
      incl.flag <- incl.flag + rowSums( incl.ids == i )
    }

    p1 <- p1[(incl.flag==0),]

    #names(p1)[2] <- "minimization.successful"
    #names(p1)[3] <- "covariance.step.successful"
    #names(p1)[4] <- "covariance.step.warnings"
    #names(p1)[5] <- "estimate.near.boundary"

    #cat(nrow(p1))
    if (min.failed) {
      p1 <- subset(p1, minimization.successful == 1)
    }
    if (cov.failed) {
      p1 <- subset(p1, covariance.step.successful == 1)
    }
    if (cov.warnings) {
      p1 <- subset(p1, covariance.step.warnings == 0)
    }
    if (boundary) {
      p1 <- subset(p1, estimate.near.boundary == 0)
    }

    ## stats and plots for each- single

    for (i in index:nparams)
    {
      if (mode(p1[[i]]) == "numeric" &&
          sum(p1[[i]],na.rm=T))
      {
        sp <- summary(p1[[i]])
        # IQR <- diff(summary(p1[[i]])[c(5,2)])
        dp <- density(p1[[i]], na.rm=T)
        parmlabel <- names(p1)[i]

        #pdf(file=paste("bootstrap.", parmlabel, ".pdf", sep=""), paper="special",
        #  title=paste("Bootstrap results - ", parmlabel, sep=""),width=10,height=7 )

        qu <- quantile(p1[[i]], c(0.025, 0.975), na.rm=T)

        legend=paste("n = ", nrow(p1), sep="")
        if (showmean) {
          legend=paste(legend, "; Mean = ", sp[4], sep="")
        }
        if (showmedian) {
          legend=paste(legend, "; Median = ", sp[3], sep="")
        }
        if (showoriginal) {
          legend=paste(legend, "; Orig = ", o1[[i]], sep="")
        }
######################################
    mylabel <- paste("Bootstrap results - ", parmlabel, sep="")
    pgraph = ggraphics(ps=6)
    add(pk.dialog.notebook, pgraph, label = mylabel,
          override.closebutton = TRUE)
######################################
        
        part.data <- data.frame(p1[[i]])
        colnames(part.data) <- parmlabel

        call.command <- "histogram"
        call.final <- getPKHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                    hist.x=parmlabel,
                                    hist.bin = svalue(.$widgets[["number of bins"]]),
                                    hist.main="",
                                    hist.xlab=parmlabel,
                                    hist.ylab="",
                                    hist.type= "",
                                    hist.cond ="",
                                    hist.layout_x = "",
                                    hist.layout_y = "",
                                    hist.data = part.data
                                    )

        print(call.final)
        setPKCode(list(pkcall=call.command, pklist=call.final))
        #setPKGGobi(list(x=myx, y=NULL))
      }
    }

    ## stats and plots for each - 6 per sheet

    total  <- 0
    bspage <- 0

    for (i in index:nparams) {
      if (mode(p1[[i]]) == "numeric" &&
          sum(p1[[i]],na.rm=T)) {
        sp <- summary(p1[[i]])
        # IQR <- diff(summary(p1[[i]])[c(5,2)])
        dp <- density(p1[[i]], na.rm=T)
        parmlabel <- names(p1)[i]

        if (total == 0) {
          bspage <- bspage + 1
          #pdf(file=paste("bootstrap.page", bspage, ".pdf", sep=""), paper="special",
           # title="Bootstrap results",width=10,height=7)
          par(mfrow = c(3,3))
        }
        total <- total + 1

        qu <- quantile(p1[[i]], c(0.025, 0.975), na.rm=T)

        legend=paste("n = ", nrow(p1), sep="")
        if (showmean) {
          legend=paste(legend, "; Mean = ", sp[3], sep="")
        }
        if (showmedian) {
          legend=paste(legend, "; Median = ", sp[4], sep="")
        }
        if (showoriginal) {
          legend=paste(legend, "; Orig = ", o1[[i]], sep="")
        }

######################################
    mylabel <- paste("Bootstrap results - ", parmlabel, sep="")
    pgraph = ggraphics(ps=6)
    add(pk.dialog.notebook, pgraph, label = mylabel,
          override.closebutton = TRUE)
######################################

        part.data <- data.frame(p1[[i]])
        colnames(part.data) <- parmlabel

        call.command <- "histogram"
        call.final <- getPKHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                    hist.x=parmlabel,
                                    hist.bin = svalue(.$widgets[["number of bins"]]),
                                    hist.main=paste("Bootstrap results - ", parmlabel, sep=""),
                                    hist.xlab=parmlabel,
                                    hist.ylab="",
                                    hist.type= "",
                                    hist.cond ="",
                                    hist.layout_x = "",
                                    hist.layout_y = "",
                                    hist.data = part.data
                                    )

        print(call.final)
        setPKCode(list(pkcall=call.command, pklist=call.final))
        #setPKCode(list(pkcall=lattice.call, pklist=lattice.list))
        
        if (total == 9) {
          total <- 0
          #dev.off()
        }
      }
    }



}


psn.bootstrap.ggobiImageHandler = function(.,h,...)
{
    ErrorMessage("No ggobi instance available for this data.")
}


val.okButtonHandler = function(.,h,...)
{

    pgraph = ggraphics(ps=6)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = FALSE)
    x <- getCurrentData()[[svalue(.$widgets[[1]])]]
    y <- getCurrentData()[[svalue(.$widgets[[2]])]]

    print(xyplot(y~x, xlab=svalue(.$widgets[[1]]), ylab=svalue(.$widgets[[2]])))


    #dispose(.$window)
}

val.ggobiImageHandler = function(.,h,...)
{
      g <- ggobi(getCurrentData())
      display(g[1], pmode="Scatterplot Display", vars=list(X=svalue(.$widgets[[1]]), Y=svalue(.$widgets[[2]])))
}

psn.sim.okButtonHandler = function(.,h,...)
{
    model.1.path <- svalue(.$widgets[["Directory path_Model 1:"]])
    model.1.folder.pattern <- svalue(.$widgets[["Simulation folder pattern_Model 1:"]])
    model.1.fileName <- svalue(.$widgets[["NONMEM result file name_Model 1:"]])
    
    model.2.path <- svalue(.$widgets[["Directory path_Model 2:"]])
    model.2.folder.pattern <- svalue(.$widgets[["Simulation folder pattern_Model 2:"]])
    model.2.fileName <- svalue(.$widgets[["NONMEM result file name_Model 2:"]])

    cond.var <- svalue(.$widgets[["Plot variable:"]])
    id.var <- svalue(.$widgets[["Patient ID:"]])
    
    if (any(c(model.1.path,model.2.path)) == "")
    {
        ErrorMessage("Please input Directory path first!")
        return(invisible(NULL))
    }
    else
    {
        extractSimData(model.1.path, model.1.folder.pattern, model.1.fileName, id.var, cond.var)
    }

}

###############################################################################
com.map.okButtonHandler = function(.,h,...)
{
    data.name <- getComDataName()
    data.name1 <- data.name[1]
    name1 <- colnames(getCurrentData(data.name1))

    key.df <- data.frame(data1=rep(NA, length(name1)), data2=rep(NA, length(name1)))
    sapply(1:length(name1), function(i) key.df[i,] <<- c(name1[i], svalue(.$widgets[[name1[i]]])) )
    setComMap(key.df)

    match.ind <- which(key.df[,2] != "")
    if (length(match.ind) < 1)
    {
        ErrorMessage("Please choose at least One mapping colnames!")
        return(invisible(NULL))
    }

    dispose(.$window)
    
    data1 <- getCurrentData(data.name[1])
    data2 <- getCurrentData(data.name[2])
    
    match.term <- key.df[match.ind,]
    merge.data <- merge(data1, data2, by.x = match.term[,1], by.y=match.term[,2], all=FALSE)
    ## TODO: check merge completeness ?
    
    datatype <- "ModelComparison"
    thisDataName <- paste(getTotalDataLen() + 1, "_", datatype , sep="")
    setDatasets(merge.data, thisDataName) # use no as data name
    setCurrentDataType(datatype, thisDataName)

    ptable=gtable(merge.data, multiple=TRUE, expand=TRUE)
    pkmain.add(ptable, as.character(thisDataName), override.closebutton = TRUE)
    
    svalue(pmg.statusBar) <- "Mapping for model comparison is configured successfully."
}
com.hist.okButtonHandler = function(.,h,...)
{

    ## translate parameters
    currentPage <- svalue(pmg.dialog.notebook)
    currentData <- getCurrentData(currentPage)

    call.command <- "histogram"
    myx <- svalue(.$widgets[["x1"]])

    call.final <- getHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.bin = svalue(.$widgets[["number of bins"]]),
                                hist.main=svalue(.$widgets[["main"]]),
                                hist.xlab=myx,
                                hist.ylab=svalue(.$widgets[["ylab"]]),
                                hist.type= svalue(.$widgets[["type"]]),
                                hist.cond =svalue(.$savewd[["cond"]]),
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = currentData
                                )

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          override.closebutton = TRUE)
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL))

    call.command <- "histogram"
    myx <- svalue(.$widgets[["x2"]])

    call.final <- getHistCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                                hist.x=myx,
                                hist.bin = svalue(.$widgets[["number of bins"]]),
                                hist.main=svalue(.$widgets[["main"]]),
                                hist.xlab=myx,
                                hist.ylab=svalue(.$widgets[["ylab"]]),
                                hist.type= svalue(.$widgets[["type"]]),
                                hist.cond =svalue(.$savewd[["cond"]]),
                                hist.layout_x = svalue(.$savewd[["layout_x"]]),
                                hist.layout_y = svalue(.$savewd[["layout_y"]]),
                                hist.data = currentData
                                )

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          override.closebutton = TRUE)
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=myx, y=NULL))

    ## figure 3
    if (svalue(.$savewd[["graphics"]]) == "lattice")
    {
        lattice.call3 <- "densityplot"

        merge.name <- colnames(currentData)
        data3.index <- merge.name[-c(grep("x", merge.name), grep("y", merge.name))]
        a.index <- c(data3.index, svalue(.$widgets[["x1"]]))
        a <- currentData[,a.index]
        a$com.group <- getComDataName()[1]
        colnames(a) <- c(data3.index, svalue(.$widgets[["x1"]]), "com.group")
        b.index <- c(data3.index, svalue(.$widgets[["x2"]]))
        b <- currentData[,b.index]
        b$com.group <- getComDataName()[2]
        colnames(b) <- c(data3.index, svalue(.$widgets[["x1"]]), "com.group")

        data3 <- rbind(a, b)

        x3 <- paste( "~", svalue(.$widgets[["x1"]]), sep="")
        lattice.list.3 <- list(x=x3, xlab=svalue(.$widgets[["xlab"]]), group=quote(com.group),
                            ylab=svalue(.$widgets[["ylab"]]), type= svalue(.$widgets[["type"]]),
                            main=svalue(.$widgets[["main"]]) , data= data3, auto.key=TRUE)



        ## TODO: should para compare...require conditional var too?
        if ( is.null(.$savewd[["cond"]]) || (svalue(.$savewd[["cond"]])==""))
        {
            lattice.list.3$x <- as.formula(lattice.list.3$x)
        }
        else
        {
            lattice.list.3$x <-  as.formula(paste(lattice.list.3$x, "|", svalue(.$savewd[["cond"]]), sep=" "))
            lattice.list.3$layout <- as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]])))
        }
        pgraph = ggraphics(ps=6)
        size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
        add(pk.dialog.notebook, pgraph, label = message,
              override.closebutton = TRUE)

        print(call.final <- do.call(lattice.call3, lattice.list.3))
        setPKCode(list(pkcall=lattice.call3, pklist=call.final))
        setPKGGobi(list(x=svalue(.$widgets[["x1"]]), y=svalue(.$widgets[["x2"]])))
    }
    else
    {
        ggplot.call3 <- "ggplot2"

        data3 <- currentData[,c(svalue(.$widgets[["x1"]]), svalue(.$widgets[["x2"]]))]
        ggplot.txt1 <- "df <- data.frame(data3)"
        ggplot.txt2 <- "df.m <- melt(df)"

        ggplot.txt3 <- "ggplot(df.m) + geom_freqpoly(aes(x = value,
               y = ..density.., colour = variable))"
        ggplot.txt <- paste(ggplot.txt1, ggplot.txt2, ggplot.txt3, sep=";")

        pgraph = ggraphics(ps=6)
        size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
        add(pk.dialog.notebook, pgraph, label = message,
              override.closebutton = TRUE)

        print(call.final <- print(eval(parse(text=ggplot.txt))))
        setPKCode(list(pkcall=ggplot.call3, pklist=call.final))
        setPKGGobi(list(x=svalue(.$widgets[["x1"]]), y=svalue(.$widgets[["x2"]])))
    }    
    
}
com.hist.okButtonHandler.0125 = function(.,h,...)
{

    ## translate parameters

    
    if (svalue(.$savewd[["graphics"]]) == "lattice")
    {
        currentData <- getCurrentData()
        lattice.call <- "histogram"
        x1 <- paste( "~", svalue(.$widgets[["x1"]]), sep="")
        x2 <- paste( "~", svalue(.$widgets[["x2"]]), sep="")
        lattice.list.1 <- list(x=x1, xlab=svalue(.$widgets[["x1"]]),
                            ylab=svalue(.$widgets[["ylab"]]), type= svalue(.$widgets[["type"]]),
                            main=svalue(.$widgets[["main"]]) , data= currentData) 
        setPKGGobi(list(x=svalue(.$widgets[["x1"]]), y=NULL))
        lattice.list.2 <- list(x=x2, xlab=svalue(.$widgets[["x2"]]),
                            ylab=svalue(.$widgets[["ylab"]]), type= svalue(.$widgets[["type"]]),
                            main=svalue(.$widgets[["main"]]) , data= currentData)
        setPKGGobi(list(x=svalue(.$widgets[["x2"]]), y=NULL))

        lattice.call3 <- "densityplot"

        merge.name <- colnames(currentData)

        data3.index <- merge.name[-c(grep("x", merge.name), grep("y", merge.name))]
        #com.name <- unlist(strsplit(svalue(.$widgets[["x1"]]), "\\."))[1]
        a.index <- c(data3.index, svalue(.$widgets[["x1"]]))
        a <- currentData[,a.index]
        a$com.group <- getComDataName()[1]
        colnames(a) <- c(data3.index, svalue(.$widgets[["x1"]]), "com.group")
        b.index <- c(data3.index, svalue(.$widgets[["x2"]]))
        b <- currentData[,b.index]
        b$com.group <- getComDataName()[2]
        colnames(b) <- c(data3.index, svalue(.$widgets[["x1"]]), "com.group")
        
        data3 <- rbind(a, b)

        x3 <- paste( "~", svalue(.$widgets[["x1"]]), sep="")
        lattice.list.3 <- list(x=x3, xlab=svalue(.$widgets[["xlab"]]), group=quote(com.group),
                            ylab=svalue(.$widgets[["ylab"]]), type= svalue(.$widgets[["type"]]),
                            main=svalue(.$widgets[["main"]]) , data= data3, auto.key=TRUE)
        setPKGGobi(list(x=svalue(.$widgets[["x1"]]), y=svalue(.$widgets[["x2"]])))

        if ( svalue(.$widgets[["number of bins"]])!= "" )
        {
              lattice.list.1$nint = as.numeric(svalue(.$widgets[["number of bins"]]))
              lattice.list.2$nint = as.numeric(svalue(.$widgets[["number of bins"]]))
        }
        
        ## TODO: should para compare...require conditional var too?
        if ( is.null(.$savewd[["cond"]]) || (svalue(.$savewd[["cond"]])==""))
        {
            lattice.list.1$x <- as.formula(lattice.list.1$x)
            lattice.list.2$x <- as.formula(lattice.list.2$x)
            lattice.list.3$x <- as.formula(lattice.list.3$x)
        }
        else
        {
            lattice.list.1$x <-  as.formula(paste(lattice.list.1$x, "|", svalue(.$savewd[["cond"]]), sep=" "))
            lattice.list.1$layout <- as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]])))

            lattice.list.2$x <-  as.formula(paste(lattice.list.2$x, "|", svalue(.$savewd[["cond"]]), sep=" "))
            lattice.list.2$layout <- as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]])))

            lattice.list.3$x <-  as.formula(paste(lattice.list.3$x, "|", svalue(.$savewd[["cond"]]), sep=" "))
            lattice.list.3$layout <- as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]])))

        }

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = TRUE)
            print(do.call(lattice.call, lattice.list.1))

    setPKCode(list(pkcall=lattice.call, pklist=lattice.list.1))

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = TRUE)
            print(do.call(lattice.call, lattice.list.2))
    setPKCode(list(pkcall=lattice.call, pklist=lattice.list.2))

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = TRUE)
            print(do.call(lattice.call3, lattice.list.3))

    currentPage <- svalue(pk.dialog.notebook)
    setDataSpecialPlot(data3, as.character(currentPage))
    setPKCode(list(pkcall=lattice.call3, pklist=lattice.list.3))

     }
     else
     {
        ggplot.call <- "qplot"
        mytype <- switch(svalue(.$widgets[["type"]]),
                                    p = c("point"),
                                    l = c("line"),
                                    psmooth = c("point", "smooth"),
                                    lsmooth = c("line", "smooth"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram")
                                    )
    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = FALSE)
        ggplot.list.1 <- list(x=as.name(svalue(.$widgets[["x1"]])), xlab=svalue(.$widgets[["xlab"]]),
                            ylab=svalue(.$widgets[["ylab"]]), geom = mytype,
                            main=svalue(.$widgets[["main"]]) , data= getCurrentData())
         print(do.call(ggplot.call, ggplot.list.1))
                    
    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = FALSE)
        ggplot.list.2 <- list(x=as.name(svalue(.$widgets[["x2"]])), xlab=svalue(.$widgets[["xlab"]]),
                            ylab=svalue(.$widgets[["ylab"]]), geom = mytype,
                            main=svalue(.$widgets[["main"]]) , data= getCurrentData())

        print(do.call(ggplot.call, ggplot.list.2))

     }

    #dispose(.$window)
}
com.scatter.okButtonHandler = function(.,h,...)
{
    currentPage <- svalue(pmg.dialog.notebook)
    currentData <- getCurrentData(currentPage)

    call.command <- "xyplot"
    call.final <- getScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                              hist.x=svalue(.$widgets[["x"]]),
                              hist.y=svalue(.$widgets[["y"]]),
                              #hist.bin = svalue(.$widgets[["number of bins"]]),
                              hist.main=svalue(.$widgets[["main"]]),
                              hist.xlab=svalue(.$widgets[["xlab"]]),
                              hist.ylab=svalue(.$widgets[["ylab"]]),
                              hist.type= svalue(.$widgets[["type"]]),
                              hist.cond = svalue(.$savewd[["cond"]]),
                              hist.layout_x = svalue(.$savewd[["layout_x"]]),
                              hist.layout_y = svalue(.$savewd[["layout_y"]])
                              )

        pgraph = ggraphics(ps=6)
        size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
        add(pk.dialog.notebook, pgraph, label = message,
              override.closebutton = TRUE)
              
    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))
    setPKGGobi(list(x=svalue(.$widgets[["x"]]), y=svalue(.$widgets[["y"]])))
    
}
com.scatter.okButtonHandler.0125 = function(.,h,...)
{

    ## translate parameters
    if (svalue(.$savewd[["graphics"]]) == "lattice")
    {
        currentData <- getCurrentData()
        lattice.call <- "xyplot"
        x <- paste(svalue(.$widgets[["y"]]), "~", svalue(.$widgets[["x"]]), sep="")

        lattice.list<- list(x=x, xlab=svalue(.$widgets[["x"]]),
                            ylab=svalue(.$widgets[["y"]]), type= svalue(.$widgets[["type"]]),
                            main=svalue(.$widgets[["main"]]) , data= currentData)
        setPKGGobi(list(x=svalue(.$widgets[["x"]]), y=svalue(.$widgets[["y"]])))

        ## TODO: should para compare...require conditional var too?
        if ( is.null(.$savewd[["cond"]]) || (svalue(.$savewd[["cond"]])==""))
        {
            lattice.list$x <- as.formula(lattice.list$x)
        }
        else
        {
            lattice.list$x <-  as.formula(paste(lattice.list$x, "|", svalue(.$savewd[["cond"]]), sep=" "))
            lattice.list$layout <- as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]])))

        }

        pgraph = ggraphics(ps=6)
        size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
        add(pk.dialog.notebook, pgraph, label = message,
              override.closebutton = TRUE)
                print(do.call(lattice.call, lattice.list))

        setPKCode(list(pkcall=lattice.call, pklist=lattice.list))

     }
     else
     {
        ## TODO
        ggplot.call <- "qplot"
        mytype <- switch(svalue(.$widgets[["type"]]),
                                    p = c("point"),
                                    l = c("line"),
                                    psmooth = c("point", "smooth"),
                                    lsmooth = c("line", "smooth"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram")
                                    )
    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = FALSE)
        ggplot.list.1 <- list(x=as.name(svalue(.$widgets[["x1"]])), xlab=svalue(.$widgets[["xlab"]]),
                            ylab=svalue(.$widgets[["ylab"]]), geom = mytype,
                            main=svalue(.$widgets[["main"]]) , data= getCurrentData())
         print(do.call(ggplot.call, ggplot.list.1))

     }

    #dispose(.$window)
}

com.time.okButtonHandler = function(.,h,...)
{
    currentPage <- svalue(pmg.dialog.notebook)
    currentData <- getCurrentData(currentPage)

        ## TODO, devision:0, careful...try?
        y.ratio <- currentData[,svalue(.$widgets[["y1"]])] / currentData[,svalue(.$widgets[["y2"]])]

        if (svalue(.$savewd[["Transform"]]) == "y1/y2")
        {
            newdata <- data.frame(currentData, y.ratio=y.ratio)
            myy <- "y.ratio"

            ## note: you have NOT generated figure yet, so we need to + 1
            current.pk.page <- svalue(pk.dialog.notebook) + 1
            setDataSpecialPlot(newdata, as.character(current.pk.page))
            setPKGGobi(list(x=svalue(.$widgets[["x"]]), y="y.ratio"))
        }
        else
        {
            y.ratio <- log(y.ratio)
            newdata <- data.frame(currentData, log.ratio=y.ratio)
            myy <- "log.ratio"
            
            current.pk.page <- svalue(pk.dialog.notebook) + 1
            setDataSpecialPlot(newdata, as.character(current.pk.page))
            setPKGGobi(list(x=svalue(.$widgets[["x"]]), y="log.ratio"))
        }

    call.command <- "xyplot"
    call.final <- getScatterCall(hist.graphics = svalue(.$savewd[["graphics"]]),
                              hist.x=svalue(.$widgets[["x"]]),
                              hist.y= myy ,
                              #hist.bin = svalue(.$widgets[["number of bins"]]),
                              hist.main=svalue(.$widgets[["main"]]),
                              hist.xlab=svalue(.$widgets[["xlab"]]),
                              hist.ylab=svalue(.$widgets[["ylab"]]),
                              hist.type= svalue(.$widgets[["type"]]),
                              hist.cond = svalue(.$savewd[["cond"]]),
                              hist.layout_x = svalue(.$savewd[["layout_x"]]),
                              hist.layout_y = svalue(.$savewd[["layout_y"]]),
                              hist.data=newdata
                              )

    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
              override.closebutton = TRUE)

    print(call.final)
    setPKCode(list(pkcall=call.command, pklist=call.final))


}

com.time.okButtonHandler.back = function(.,h,...)
{

    ## translate parameters
    if (svalue(.$savewd[["graphics"]]) == "lattice")
    {
        currentData <- getCurrentData()
        ## TODO, devision:0, careful...try?
        y.ratio <- currentData[,svalue(.$widgets[["y1"]])] / currentData[,svalue(.$widgets[["y2"]])]

        if (svalue(.$savewd[["Transform"]]) == "y1/y2")
        {
            newdata <- data.frame(currentData, y.ratio=y.ratio)
            myx <- paste("y.ratio", "~", svalue(.$widgets[["x"]]), sep="")
            setPKGGobi(list(x=svalue(.$widgets[["x"]]), y="y.ratio"))
        }
        else
        {
            y.ratio <- log(y.ratio)
            newdata <- data.frame(currentData, log.ratio=y.ratio)
            myx <- paste("log.ratio", "~", svalue(.$widgets[["x"]]), sep="")
            setPKGGobi(list(x=svalue(.$widgets[["x"]]), y="log.ratio"))
        }
        
        lattice.call <- "xyplot"
        lattice.list<- list(x=myx, xlab=svalue(.$widgets[["x"]]),
                            ylab=svalue(.$savewd[["Transform"]]), type= svalue(.$widgets[["type"]]),
                            main=svalue(.$widgets[["main"]]) , data= newdata)

        ## TODO: should para compare...require conditional var too?
        if ( is.null(.$savewd[["cond"]]) || (svalue(.$savewd[["cond"]])==""))
        {
            lattice.list$x <- as.formula(lattice.list$x)
        }
        else
        {
            lattice.list$x <-  as.formula(paste(lattice.list$x, "|", svalue(.$savewd[["cond"]]), sep=" "))
            lattice.list$layout <- as.numeric(c(svalue(.$savewd[["layout_x"]]), svalue(.$savewd[["layout_y"]])))

        }

        pgraph = ggraphics(ps=6)
        size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
        add(pk.dialog.notebook, pgraph, label = message,
              override.closebutton = TRUE)
                print(do.call(lattice.call, lattice.list))


        currentPage <- svalue(pk.dialog.notebook)
        setDataSpecialPlot(newdata, as.character(currentPage))
        setPKCode(list(pkcall=lattice.call, pklist=lattice.list))
     }
     else
     {
        ## TODO
        ggplot.call <- "qplot"
        mytype <- switch(svalue(.$widgets[["type"]]),
                                    p = c("point"),
                                    l = c("line"),
                                    psmooth = c("point", "smooth"),
                                    lsmooth = c("line", "smooth"),
                                    percent = c("histogram"),
                                    count = c("histogram"),
                                    density = c("histogram")
                                    )
    pgraph = ggraphics(ps=6)
    size(pgraph) <- c(getSubHeight()*0.5, getSubWidth()*0.5)
    add(pk.dialog.notebook, pgraph, label = message,
          #pageno = 3,
          override.closebutton = FALSE)
        ggplot.list.1 <- list(x=as.name(svalue(.$widgets[["x1"]])), xlab=svalue(.$widgets[["xlab"]]),
                            ylab=svalue(.$widgets[["ylab"]]), geom = mytype,
                            main=svalue(.$widgets[["main"]]) , data= getCurrentData())
         print(do.call(ggplot.call, ggplot.list.1))

     }

    #dispose(.$window)
}

compare.pred.ggobiImageHandler = function(.,h,...)
{
    x <- getDatasets()[[svalue(.$widgets[[1]])]]
    y <- getDatasets()[[svalue(.$widgets[[2]])]]

    current.lab <- names(pk.dialog.notebook)[svalue(pk.dialog.notebook)]
    current.name <- unlist(strsplit(current.lab, "\\."))[1]

    mypred <- c("PRED","RES","IPRE","WRES")
    choice <- names(x) %in% mypred
    base.fit <- x[,!choice]
    fit1 <- x[,mypred]
    fit2 <- y[,mypred]

    names(fit1) <- paste(names(fit1),".fit1",sep="")
    names(fit2) <- paste(names(fit2),".fit2",sep="")

    ## calculate ind ratio and group ratio
    ## group ratio NEED TO RETHINK
    all.data <- data.frame(base.fit, fit1, fit2)

    all.data <- all.data[all.data$PRED.fit2!=0,]
    all.data <- all.data[all.data$RES.fit2!=0,]
    all.data <- all.data[all.data$IPRE.fit2!=0,]
    all.data <- all.data[all.data$WRES.fit2!=0,]
    all.data <- data.frame(all.data,
                  PRED.ratio = log(abs(all.data$PRED.fit1/all.data$PRED.fit2)),
                  RES.ratio = log(abs(all.data$RES.fit1/all.data$RES.fit2)),
                  IPRE.ratio = log(abs(all.data$IPRE.fit1/all.data$IPRE.fit2)),
                  WRES.ratio = log(abs(all.data$WRES.fit1/all.data$WRES.fit2)))
################
      x.name <- "TIME"
      x.ind <- which(colnames(all.data)== x.name)
      y.name <- current.name
      y.ind <- which(colnames(all.data)== y.name)
      old.ind <- c(1:length(colnames(all.data)))
      old.ind <- old.ind[-c(x.ind, y.ind)]

      all.data <- all.data[c(x.ind, y.ind, old.ind)]

################
    g <- ggobi(all.data)

}
