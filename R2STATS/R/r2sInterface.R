#--------------------------------------------------------------------------------------------------
#
#                      R2STATS: A Graphical User Interface for GLM and GLMM in R
#                       Yvonnick Noel, University of Brittany, Rennes 2, France
#                                            2006-2011
#
#--------------------------------------------------------------------------------------------------
#                                      INTERFACE PROTOTYPE
#--------------------------------------------------------------------------------------------------

r2stats = proto(

  # Create the GUI
  create = function(.) {

    # Main Window
   .$version = packageDescription("R2STATS")$Version
    mainWindowTitle = paste("R2STATS",Sys.info()["sysname"],.$version,sep="-")
   .$mainWindow = gwindow(mainWindowTitle,visible=FALSE)
    add(.$mainWindow,bigGroup <- ggroup(horizontal=FALSE),expand=TRUE)

    # Menus
    aLoad    = gaction(label=.$translate("Load models"),handler=.$loadModels)
    aSave    = gaction(label=.$translate("Save models"),handler=.$saveModels)
    aClose   = gaction(label=.$translate("Quit"),icon="quit",handler=function(h,...) dispose(.$mainWindow))
    aOptions = gaction(label=.$translate("Options"),handler=.$editOptions)
    aAbout   = gaction(label=.$translate("About..."),icon="about",handler=.$aboutR2stats)

    tmp = list(Session=list(load=aLoad,save=aSave,sep=list(separator=TRUE),quit=aClose),
               Help  =list(options=aOptions,about=aAbout))
    names(tmp) = .$translate(names(tmp))
   .$menu = gmenu(tmp,cont=.$mainWindow)

    # Tableau d'onglets
    add(bigGroup,.$mainNotebook <- gnotebook(closebuttons=TRUE,dontCloseThese=1:5),expand=TRUE)

    #-------------------------------------------------- File Tab ---------------------------------------------
    add(.$mainNotebook,   dataGroup <- ggroup(horizontal=FALSE),label=.$translate("Files"),override.closebutton=TRUE,expand=TRUE)
    add(dataGroup,dataNb    <- gnotebook(),expand=TRUE)

    # Load a file from a local or remote URL
    add(dataNb,        urlGroup     <- ggroup(horizontal=FALSE),label=.$translate("URL"),expand=TRUE)
    add(urlGroup,      urlFrame     <- gframe(.$translate("Load a local or remote file"),horizontal=FALSE),expand=TRUE)
    add(urlFrame,      tmp          <- ggroup())
    add(tmp,          .$dataUrl     <- gedit("http://yvonnick.noel.free.fr/data/",handler=.$dataLoad),expand=TRUE)
    addhandlerchanged(.$dataUrl,.$updateDFLists)
    add(tmp,                           gbutton(.$translate("Browse"),handler=.$fileLoad))
    add(urlFrame,     .$hasHeader   <- gcheckbox(.$translate("Column headers on the first row"),checked=TRUE))
    add(urlFrame,     .$hasRowNames <- gcheckbox(.$translate("Row labels in the first column"),checked=FALSE))
    addSpring(urlFrame)
    add(urlFrame,      tmp          <- ggroup())
    addSpring(tmp)
    add(tmp,          .$loadFromFileButton <- gbutton(.$translate("Load"),handler=.$dataLoad))

    # Load from a library
    add(dataNb,       libGroup      <- ggroup(horizontal=FALSE),label=.$translate("Package"),expand=TRUE)
    add(libGroup,     listlibFrame  <- gframe(.$translate("Load a data file from an R package"),horizontal=FALSE),expand=TRUE)
    add(listlibFrame, g1            <- ggroup(),expand=TRUE)
    add(g1,           tmp           <- ggroup(horizontal=FALSE))
    add(tmp,                           glabel(.$translate("Available packages")))
    add(tmp,         .$libList      <- gtable(.$getLibList()),expand=TRUE)
    addHandlerClicked(.$libList,      .$updateLibDataList)
    add(g1,           tmp           <- ggroup(horizontal=FALSE),expand=TRUE)
    add(tmp,                           glabel(.$translate("Data files")))
    add(tmp,         .$libDataList  <- gtable(.$getLibDataList(),handler=.$loadDataLib), expand=TRUE)
    add(listlibFrame, tmp           <- ggroup())
    addSpring(tmp)
    add(tmp,                           gbutton(.$translate("Details"),handler=.$getDataDescription))
    add(tmp,          .$loadFromLibButton  <- gbutton(.$translate("Load"),handler=.$loadDataLib))

    svalue(dataNb) = 1

    #---------------------------------------------- Data grids Tab --------------------------------------
    add(.$mainNotebook,  gridTab    <- ggroup(),label=.$translate("Data"),override.closebutton=TRUE,expand=TRUE)
    add(gridTab,         gridGroup  <- ggroup(horizontal=FALSE), expand=TRUE)

    # Recode utility
    # add(gridGroup,  transfrm        <- gexpandgroup(.$translate("Recode and transform")))
    transfrm        <- gexpandgroup(.$translate("Recode and transform"),cont=gridGroup)
    # add(transfrm,   tmp             <- ggroup(horizontal=FALSE))
    tmp             <- ggroup(horizontal=FALSE,cont=transfrm)
    add(tmp,        glabel(.$translate("With")))
    add(tmp,       .$currentFactor  <- gdroplist(.$translate("No factor"),handler=.$printCat))
    add(transfrm,   tmp             <- ggroup(horizontal=FALSE),expand=TRUE)
    add(tmp,        glabel(.$translate("Apply (ex. \"cat1,cat2=cat12;cat3,cat4=cat34\" or log(.))")))
    add(tmp,       .$toCompute      <- gedit("",width=40))
    add(transfrm,   tmp             <- ggroup(horizontal=FALSE))
    add(tmp,        glabel(.$translate("Store in")))
    add(tmp,       .$newVarName     <- gdroplist(.$translate("No variable"),editable=TRUE))
    add(transfrm,   tmp             <- ggroup(horizontal=FALSE))
    addSpring(tmp)
    add(tmp,                           gbutton(.$translate("Run"),handler=.$varTransform))

    # Data grids
    add(gridGroup, .$gridNotebook   <- gnotebook(closebuttons=TRUE), expand=TRUE)
   .$handler.ID['onTabChange'] = addhandlerchanged(.$gridNotebook,.$changeDataset,action="tabChange")
    addhandlerchanged(.$gridNotebook,.$updateTransf, action="tabChange")

    # Define a default data grid
    dfList = .$getAllTables()
   .$currentDataName = dfList[1]

    # Data selector
    add(gridGroup,tmp             <- ggroup())
    add(tmp,          .$openGrid  <- gdroplist(dfList))
   .$handler.ID['onOpenGrid'] = addhandlerchanged(.$openGrid, handler=.$changeDataset,action="openGrid")
    add(tmp,             myicon1  <- gimage("refresh",dir="stock",handler=.$updateDFLists))
    tooltip(myicon1) = .$translate("Click to update list of loaded data frames")

    addSpring(tmp)
    add(tmp,                           gbutton(.$translate("Save"),handler=.$saveGrid))
    

    #-------------------------------------------------- Model Tab --------------------------------------------
    add(.$mainNotebook,bigFrame     <- ggroup(),label=.$translate("Models"),override.closebutton=TRUE,expand=TRUE)
    add(bigFrame,     .$modelPanedGroup   <- gpanedgroup(),expand=TRUE)
    add(.$modelPanedGroup,  leftGroup    <- ggroup(horizontal=FALSE),expand=TRUE)
    add(leftGroup,     tabFrame     <- gframe(.$translate("Active data frame"),horizontal=FALSE),expand=TRUE)
    add(tabFrame,      tmp          <- ggroup())
    add(tmp,          .$currentData <- gdroplist(dfList,handler=.$updateVarList),expand=TRUE)
    addhandlerchanged(.$currentData, .$clearModelFields)
    addhandlerchanged(.$currentData, .$updateWeightList)
    addhandlerchanged(.$currentData, .$updateConstrFactor)
    addhandlerchanged(.$currentData, .$updateGraphNumVarList)
   .$handler.ID['onChangeModelData'] = addhandlerchanged(.$currentData, .$changeDataset,action="changeModelData")
    add(tmp,             myicon2 <-   gimage("refresh",dir="stock",handler=.$updateDFLists))
    tooltip(myicon2) = .$translate("Click to update list of loaded data frames")

    add(tabFrame,     .$varList     <- gtable(.$getVarList(.$getDataName()),multiple=TRUE),expand=TRUE)
    
    addHandlerDoubleclick(.$varList, handler = .$updateFIVField, action = "add")
    addhandlerclicked(    .$varList, handler = .$updateVarSummary)
    
    add(leftGroup,     varFrame <- gframe(.$translate("Variable summary")),expand=TRUE)
    add(varFrame,     .$varSummary  <- gtable(cbind(Attribute=.$translate("None"),Value=.$translate("None"))),expand=TRUE)
    names(.$varSummary) = .$translate(names(.$varSummary))

    add(.$modelPanedGroup,  rightFrame   <- ggroup(horizontal=FALSE),expand=TRUE)
    add(rightFrame,    modelFrame   <- gframe(.$translate("Model definition"),horizontal=FALSE),expand=TRUE)

    # .$distribList = gdroplist(c("Normale","Binomiale","Poisson","Gamma","Gaussienne inverse","Multinomiale","Multinomiale ordonnée"),handler=.$updateLink)
   .$distribList = gdroplist(.$translate(c("Gaussian","Binomial","Poisson","Gamma","Inverse gaussian")),handler=.$updateLink)
   .$linkLists = list(Gaussian     = .$translate(c("Identity","Log","Inverse")),
                       Binomial     = .$translate(c("Logit","Probit","Cauchit","Log","Cloglog")),
                       Poisson      = .$translate(c("Log","Square root","Identity")),
                       Gamma        = .$translate(c("Inverse","Log","Identity")),
                       "Inverse gaussian" = .$translate(c("Inverse","Log","Identity","1/mu2")))
         #              Multinomiale = c("Logit", "Probit", "Cloglog", "Loglog","Cauchit", "Aranda-Ordaz", "Log-gamma"),
         #   "Multinomiale ordonnée" = c("Logit", "Probit", "Cloglog", "Loglog","Cauchit", "Aranda-Ordaz", "Log-gamma"))
    names(.$linkLists) = .$translate(names(.$linkLists))
   .$currentLinkList = gdroplist(.$linkLists[[svalue(.$distribList)]])
   .$modelName = gdroplist(.$translate(c("New")),editable=TRUE,handler=.$retrieveModel)
   .$dvList  = gtext(height=64,width=120)
   .$fivList = gtext(height=64,width=120)

    add(modelFrame,tmp <- ggroup())
    add(tmp,glabel(.$translate("Model name")))
    add(tmp,.$modelName)
    addSpring(tmp)
    add(tmp,gbutton(.$translate("New"),handler=.$clearModelFields))

    add(modelFrame,frm <- gframe(.$translate("Dependent variables"),horizontal=FALSE),expand=TRUE)
    add(frm,.$dvList,expand=TRUE)
    add(frm,tmp <- ggroup())
    add(tmp,gbutton(.$translate("Add"),handler=.$updateDVField))
    addSpring(tmp)
    add(tmp,gbutton(.$translate("Clear"),handler=.$clearDVField))

    add(modelFrame,frm <- gframe(.$translate("Independent variables"),horizontal=FALSE),expand=TRUE)
    add(frm,.$fivList,expand=TRUE)
    add(frm,tmp <- ggroup())
    add(tmp,gbutton(.$translate("Add"),  handler=.$updateFIVField,action="add"))
    addSpring(tmp)
    add(tmp,gbutton(" + ",      handler=.$updateFIVField,action="+"))
    add(tmp,gbutton(" : ",      handler=.$updateFIVField,action=":"))
    add(tmp,gbutton(" * ",      handler=.$updateFIVField,action="*"))
    add(tmp,gbutton(" - ",      handler=.$updateFIVField,action="-"))
    add(tmp,gbutton(" () ",     handler=.$updateFIVField,action="()"))
    add(tmp,gbutton(.$translate("Fixed"),    handler=.$updateFIVField,action="offset"))
    add(tmp,gbutton(" +1 ",     handler=.$updateFIVField,action="1"))
    add(tmp,gbutton(" +0 ",     handler=.$updateFIVField,action="-1"))
    addSpring(tmp)
    add(tmp,gbutton(" | ",      handler=.$updateFIVField,action="|"))
    add(tmp,gbutton(" (1|.) ",  handler=.$updateFIVField,action="(1|)"))
    add(tmp,gbutton(" (.| ) ",  handler=.$updateFIVField,action="(.|)"))
    add(tmp,gbutton(" / ",      handler=.$updateFIVField,action="/"))
    addSpring(tmp)
    add(tmp,gbutton(.$translate("Clear"),handler=.$clearFIVField))

    addSpring(modelFrame)
    add(modelFrame,tmp <- ggroup(),expand=TRUE)
    addSpring(tmp)
    add(tmp,layout <- glayout(homogeneous=FALSE),expand=TRUE)
    layout[1,1] = glabel(.$translate("Distribution family"))
    layout[1,2] = glabel(.$translate("Link function"))
    layout[1,3] = glabel(.$translate("Weighting variable"))
    layout[1,4] = glabel(.$translate("Constraint factor"))
    layout[2,1] =  .$distribList
    layout[2,2] =  .$currentLinkList
    layout[2,3] <- .$weightList <- gdroplist(c(.$translate("No variable"),.$getNumVarList(.$getCurrentDataName())))
    layout[2,4] <- .$structList <- gdroplist(c(.$translate("No factor"),.$translate("Constant"),.$getCatVarList(.$getCurrentDataName())))
    layout[3,1] = glabel(.$translate("Observation selection"))
    layout[3,2:4] <- .$subsetVar <- gedit("")
    addSpring(tmp)

    add(rightFrame,tmp <- ggroup())
    addSpring(tmp)
    add(tmp,gbutton(.$translate("Estimate"),handler=.$run))

    # Open the first data grid (must appear after all droplists have been created)
    add(.$gridNotebook, tmp <- ggroup(),label = .$currentDataName,expand=TRUE)
    add(tmp,  gdfedit(eval(parse(text=.$currentDataName),envir=.GlobalEnv),name=.$currentDataName), expand=TRUE)

    #------------------------------------------------ Result Tab ------------------------------------------
    add(.$mainNotebook,resBigFrame <- ggroup(horizontal=FALSE),override.closebutton=TRUE,label=.$translate("Results"))
    add(resBigFrame,.$results <- gtext(wrap=FALSE),expand=TRUE)
    add(resBigFrame,tmp <- ggroup())
    add(tmp,gbutton(.$translate("Clear"),handler=.$clearResults))

    #------------------------------------------------- Plot Tab -----------------------------------------
    add(.$mainNotebook,  graphBigFrame   <- ggroup(),override.closebutton=TRUE,label=.$translate("Plots"))
    add(graphBigFrame,  .$graphPanedGroup <- gpanedgroup(),expand=TRUE)
    add(.$graphPanedGroup, graphLeftGroup  <- ggroup(horizontal=FALSE),expand=TRUE)
    add(graphLeftGroup,                     glabel(.$translate("Plot type")))
    add(graphLeftGroup, .$plotType       <- gdroplist(.$translate(c("Regression plot","Response distribution","Fitted and observed values","Quantile-quantile plot","Residuals distribution","Fitted values and residuals","Quantile residuals")),handler=.$plotCurrentModel,action="plot"))
    add(graphLeftGroup,  graphNb         <- gnotebook(),expand=TRUE)
    add(graphNb,         graphModelGroup <- ggroup(horizontal=FALSE),label=.$translate("Model"),expand=TRUE)
    add(graphModelGroup,.$graphModelList <- gtable(.$getModelList()),expand=TRUE)
    addhandlerclicked(  .$graphModelList,  .$plotCurrentModel,action="plot")
    add(graphNb,         graphParamGroup <- ggroup(horizontal=FALSE),label=.$translate("Options"),expand=TRUE)

    add(graphParamGroup,layout1 <- glayout(),fill="")
    
    layout1[1,1,anchor=c(-1,0)] = glabel(.$translate("Legend"))
    layout1[1,2] <- .$legendLoc <- gdroplist(.$translate(c("None","Right","Left","Top","Bottom")),handler=.$plotCurrentModel,action="plot")
    svalue(.$legendLoc,index=TRUE) = 2
    layout1[2,2] <- .$legendCols <- gdroplist(paste(1:10,"col."),handler=.$plotCurrentModel,action="plot")
    layout1[3,1,anchor=c(-1,0)] = glabel(.$translate("X-axis"))
    layout1[3,2] <- .$graphAxisX <- gdroplist(c(.$translate("Default"),.$getNumVarList(.$getCurrentDataName())),handler=.$plotCurrentModel,action="plot")
    layout1[4,1,anchor=c(-1,0)] = glabel(.$translate("X-limits"))
    layout1[4,2] <- .$graphLimitsX <- gedit("",handler=.$plotCurrentModel,action="plot")
    layout1[5,1,anchor=c(-1,0)] = glabel(.$translate("Y-limits"))
    layout1[5,2] <- .$graphLimitsY <- gedit("",handler=.$plotCurrentModel,action="plot")
    layout1[6,1,anchor=c(-1,0)] <- .$translate("Selection")
    layout1[6,2] <- .$groupList <- gdroplist(c(.$translate("All groups")),handler=.$plotCurrentModel,action="plot")
    
    add(graphParamGroup, addGroup <- gframe(.$translate("Add")), expand=TRUE)
    add(addGroup, layout2 <- glayout(),expand=TRUE)
    layout2[1,1] <- .$addData <- gcheckbox(.$translate("Data"),checked=TRUE,handler=.$plotCurrentModel,action="plot")
    layout2[2,1] <- .$addModel <- gcheckbox(.$translate("Model"),checked=TRUE,handler=.$plotCurrentModel,action="plot")
    layout2[3,1] <- .$addCondMeans <- gcheckbox(.$translate("Conditional\nmeans"),checked=FALSE,handler=.$plotCurrentModel,action="plot")
    layout2[4,1] <- .$addRandCurves <- gcheckbox(.$translate("Random\ncurves"),checked=TRUE,handler=.$plotCurrentModel,action="plot")
    layout2[1,2] <- .$addRefLine <- gcheckbox(.$translate("Reference\nline"),handler=.$plotCurrentModel,action="plot")
    layout2[2,2] <- .$addGrid <- gcheckbox(.$translate("Grid"),handler=.$plotCurrentModel,action="plot")
    layout2[3,2] <- .$addNoise <- gcheckbox(.$translate("Noise"),handler=.$plotCurrentModel,action="plot")
    layout2[4,2] <- .$addSmooth <- gcheckbox(.$translate("Smoothing"),handler=.$plotCurrentModel,action="plot")
    
    add(.$graphPanedGroup,graphRightGroup <- ggroup(horizontal=FALSE),expand=TRUE)
    add(graphRightGroup,.$plotArea <- ggraphics())
    svalue(graphNb) = 1

    addhandlerclicked(.$plotArea,.$getXY)
    add(graphRightGroup, tmp <- ggroup())
    addSpring(tmp)
    add(tmp,.$obsId <- glabel(""))
    addSpring(tmp)
    add(tmp,gbutton(.$translate("Save"), handler=.$savePlot,action="save"))
    add(tmp,gbutton(.$translate("Copy"), handler=.$copyPlot,action="save"))
    add(tmp,gbutton(.$translate("Replot"),handler=.$plotCurrentModel,action="plot"))

    #--------------------------------------------------- Model comparison tab -------------------------------------
    add(.$mainNotebook,compBigFrame <- ggroup(horizontal=FALSE),override.closebutton=TRUE,label=.$translate("Comparisons"))
    add(compBigFrame,modelGroup <- gframe(.$translate("Fitted models")),expand=TRUE,horizontal=FALSE)
    add(modelGroup,.$modelList <- gtable(.$getModelList(),multiple=TRUE),expand=TRUE)
    add(compBigFrame,tmp <- ggroup())
    add(tmp,gbutton(.$translate("Delete"),handler=.$deleteModels))
    addSpring(tmp)
    add(tmp,gbutton(.$translate("Select all"),handler=.$selectAllModels))
    add(tmp,gbutton(.$translate("Compare"),handler=.$compareModels))

    # Status bar
    add(bigGroup,.$status <- gstatusbar(.$translate("Status: Ready.")))

    # Global variables
    .$currentPlot.XY = NULL
  },
  ### Make the R2STATS main window visible
  show = function(.) {

    if(is.null(.$mainWindow)) {
      gmessage(.$translate("Error: GUI not created."))
      return()
    }
    
    # Seems necessary to have the graphics tab appear first for it to be the default plot
    # or the default cairoDevice would pop up
    svalue(.$mainNotebook) = 5

    # Popup main window
    visible(.$mainWindow) = TRUE
   .$setPlotParams()
    
    # Model tab appears first
    svalue(.$mainNotebook) = 3
    
    # Set model and graph panel sizes
    svalue(.$modelPanedGroup) = .25
    svalue(.$graphPanedGroup) = .40

  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS FILE MANAGEMENT METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Data file selector
  fileChoose = function(.,type) {
  
    # Under W32 (change backslashes into slashes beforehand)
    if(Sys.info()["sysname"] == "Windows") invisible(file.choose())
    
    # Under Linux
    else   {
      filter = list("All"=list(patterns = c("*")))
      names(filter) = .$translate("All")

      invisible(gfile(text=.$translate("File selector"),type=type,filter=filter))
    }
  },
  ### Data file selection (from local disk)
  fileLoad = function(.,h,...) {

    filename = .$fileChoose("open")

    if(is.na(filename)) return()
    if(filename == "")  return()

    # The working directory is implicitly set
    setwd(dirname(filename))

    # This change triggers dataLoad() and updateDFLists() handlers in chain
    svalue(.$dataUrl) = filename
  },
  ### Data file selection (from direct path or URL input)
  dataLoad = function(.,h,...) {
    
    filename = .$trim(svalue(.$dataUrl))
    filename = gsub("\\\\","/",filename)

    if(is.na(filename)) return()
    if(filename == "")  return()
    
    .$setStatus(.$translate("Status: Download in progress. Please wait..."))

    # Extract file extension
    name.ext = .$getFileExtension(filename)
    tabname =  .$removeSpaces(.$getBaseName(filename))
    r.names = NULL
    if(svalue(.$hasRowNames)) r.names = 1

    if(tolower(name.ext) == "csv") res = try(eval(parse(text=paste("assign('",tabname,"',read.csv2('",filename,"',header=",svalue(.$hasHeader),",row.names=",r.names,"),envir = .GlobalEnv)",sep=""))),silent=TRUE)
    else                           res = try(eval(parse(text=paste("assign('",tabname,"',read.table('",filename,"',header=",svalue(.$hasHeader),",row.names=",r.names,"),envir = .GlobalEnv)",sep=""))),silent=TRUE)
    
    if(inherits(res,"try-error")) {
      gmessage(.$translate("Download failure: The server might be down\nor the file name incorrect."))
     .$setStatus(.$translate("Status: Ready."))
      return()
    }
    
    # Only close the corresponding grid in case we are reloading an existing one
    blockhandler(.$gridNotebook,.$handler.ID['onTabChange'])
    displayed.grids = names(.$gridNotebook) 
    if(tabname %in% displayed.grids) {
      pos = which(displayed.grids == tabname)
      svalue(.$gridNotebook) = pos
      dispose(.$gridNotebook)
    }
    unblockhandler(.$gridNotebook,.$handler.ID['onTabChange'])

    # Update list of data frames
   .$updateDFLists(h,...)

    # Reset current data frame: Need the DF list to be refreshed first
   .$currentDataName = tabname
    svalue(.$currentData) = tabname
    
    # Switch to data grid tab
    svalue(.$mainNotebook) = 2

   .$setStatus(.$translate("Status: Ready."))

  },
  ### Get the list of installed libraries
  getLibList = function(.) {
    .packages(all.available = TRUE)
  },
  ### Get the list of datasets available in the selected library
  getLibDataList = function(.) {
  
    dl = data(package=svalue(.$libList))$results[,3:4]
    if(is.matrix(dl)) ndset = nrow(dl)
    if(is.vector(dl)) ndset = 1        # data() doesn't return a matrix when there is only one dataset!

    if(ndset==0)      { dl = data.frame(Table=.$translate("No table"),Description=.$translate("No description"),stringsAsFactors=FALSE) }
    else if(ndset==1) { dl = data.frame(Table=dl[1],  Description=dl[2],   stringsAsFactors=FALSE) }
    # else            { colnames(dl) <- .$translate(c("Tableau","Description")) }
    colnames(dl) <- .$translate(colnames(dl))
    dl
  },
  ### Update the list of available datasets upon library selection
  updateLibDataList = function(.,h,...) {
    .$libDataList[,] = .$getLibDataList()
  },
  ### Load a dataset from a selected library
  loadDataLib = function(.,h,...) {

    if(.$debug) cat("Function: LoadDataLib\n")
    
    # Load a data frame from a package
    tabname = svalue(.$libDataList)
    eval(parse(text=paste("data(list='",tabname,"',package='",svalue(.$libList),"')",sep="")),envir=.GlobalEnv)
    cl = class(eval(parse(text=tabname),envir=.GlobalEnv))

    if(cl != "data.frame")	{
      gmessage(.$translate("This table class cannot be edited in the GUI.\nIt has nevertheless been loaded in memory."))
      return()
    }
	
    # Only close the corresponding grid in case we are reloading an existing one
    blockhandler(.$gridNotebook,.$handler.ID['onTabChange'])
    displayed.grids = names(.$gridNotebook) 
    if(tabname %in% displayed.grids) {
      pos = which(displayed.grids == tabname)
      svalue(.$gridNotebook) = pos
      dispose(.$gridNotebook)
    }
    unblockhandler(.$gridNotebook,.$handler.ID['onTabChange'])

    # Refresh DF list
   .$updateDFLists(h,...)
   
    # Set to current data
   .$currentDataName = tabname
    svalue(.$currentData) = tabname

    # Switch to the datagrid tab
    svalue(.$mainNotebook) = 2
  },
  ### Get the vector of level names in a given categorical variable in a data table
  getFacLevelList = function(.,dataname,varname) {

     level.names = .$translate("None")
     datatab = eval(parse(text=dataname),envir=.GlobalEnv)
     if(is.factor(datatab[,varname]))    level.names = levels(datatab[,varname])
     if(is.character(datatab[,varname])) level.names = sort(unique(datatab[,varname]))
     
     l = cbind(Levels=level.names)
     names(l) = .$translate("Levels")
     l
  },
  ### Display the help file about this dataset
  getDataDescription = function(.,h,...) {
  
    if(Sys.info()["sysname"] == "Windows")
      eval(parse(text=paste("help('",svalue(.$libDataList),"',package='",svalue(.$libList),"')",sep="")),envir=.GlobalEnv)

    else
      eval(parse(text=paste("ghelp('",svalue(.$libDataList),"',package='",svalue(.$libList),"',container=TRUE)",sep="")),env=.GlobalEnv)

  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS DATA MANAGEMENT METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Get the selected data frame
  getDataName = function(.) {
    return(svalue(.$currentData))
  },
  ### Get the active dataset
  getCurrentDataName = function(.) {
    return(.$currentDataName)
  },
  ### Update the table list in the grid tab
  updateDFLists = function(.,h,...) {

    if(.$debug) cat("Function: UpdateDfList\n")
    
    available.tables = .$getAllTables()

    if(is.null(available.tables) || !length(available.tables)) {
      available.tables = .$translate("No table")
      nt = 0
      return()
    }
    
    # Block all other handlers
    blockhandler(.$currentData, .$handler.ID['onChangeModelData'])
    blockhandler(.$gridNotebook,.$handler.ID['onTabChange'])
    blockhandler(.$openGrid,    .$handler.ID['onOpenGrid'])
    
    # Add the new file names to droplists
   .$openGrid[]    = available.tables
    svalue(.$openGrid) = .$currentDataName
   .$currentData[] = available.tables
    svalue(.$currentData) = .$currentDataName

    # Unblock handlers
    unblockhandler(.$currentData, .$handler.ID['onChangeModelData'])
    unblockhandler(.$gridNotebook,.$handler.ID['onTabChange'])
    unblockhandler(.$openGrid,    .$handler.ID['onOpenGrid'])
    
    if(.$debug) cat("Active data table:",.$currentDataName,"\n")
    
    # Display number of tables in status bar
    nt = length(available.tables)
   .$setStatus(paste(.$translate("Status:"),nt,.$translate("table(s) in workspace.")))
  },
  ### Update the data grids upon grid selector change
  updateGrid = function(.,h,...) {

    if(.$debug) cat("Function: updateGrid\n")
    
    # dataName = .$currentDataName
    dataName = svalue(.$openGrid)

    # Warning: a droplist may be temporarily set to NULL when changed
    if(is.null(dataName)) return()
    
    # No available table in working memory
    if(dataName == .$translate("No table")) return()
    
    # Which tables are already displayed?
    displayed.grids = names(.$gridNotebook)
    
    # None: Open the table
    if(is.null(displayed.grids)) {
     .$showData(dataName)
      return()
    }
    
    # Which one is visible?
    which.visible = displayed.grids[svalue(.$gridNotebook)]

    # Already selected: Useless to trigger all other handlers
    if(which.visible == dataName) return()

    # Put the data grid to the foreground if already opened
    if(dataName %in% displayed.grids) {
      svalue(.$gridNotebook) = which(displayed.grids == dataName)
      return()
    }
    
    # ... or open it
    .$showData(dataName)
  },
  ### Variable transformations
  varTransform = function(.,h,...) {

    # Current data tab
    displayed.grids = names(.$gridNotebook)
    if(length(displayed.grids)==0) return()
    
    # Current data  table and variables
    dataname = displayed.grids[svalue(.$gridNotebook)]
    if(.$debug) cat("Function: varTransform, dataname: ",dataname,"\n")
    varList = .$getVarList(dataname)
    varTypes = varList[,.$translate("Type")]
    names(varTypes) = varList[,.$translate("Variables")]

    transExp = .$trim(svalue(.$toCompute))
    if(transExp == "") return()

    sourceVar = svalue(.$currentFactor)
    destVar   = svalue(.$newVarName)
    sourceAccess = paste(dataname,"[,'",sourceVar,"']",sep="")
    destAccess = paste(dataname,"[,'",destVar,"']",sep="")
    source.isFactor = varTypes[sourceVar] == .$translate("F")
    
    # RECODING
    if(source.isFactor) {
    
      eval(parse(text=paste(destAccess,"= as.character(",sourceAccess,")")),envir=.GlobalEnv)

      # Parse and execute commands
      commands = unlist(strsplit(transExp,";"))
      for(com in commands) {
        com = .$removeSpaces(com)
        arguments = unlist(strsplit(com,"="))
        
        # No equal sign
        if(length(arguments) != 2) {
        
          if(!length(grep(":",com,fixed=TRUE))) { gmessage(.$translate("Syntax problem in your command.")) ; return() }
          
          # Combining factors
          eval(parse(text=paste(destAccess,"<-with(",dataname,",",com,")[drop=TRUE]",sep="")),envir=.GlobalEnv)
          return()
        }
        
        # Recoding a factor
	      sourceCat = unlist(strsplit(arguments[1],","))
	      newCat = arguments[2]
	      eval(parse(text=paste(dataname,"[",destAccess," %in% c('",paste(sourceCat,collapse="','"),"'),'",destVar,"'] ='",newCat,"'",sep="")),envir=.GlobalEnv)
      }
	
	  # Recode to a factor or to a numeric (indicator) variable
	  if(any(is.na(as.numeric(eval(parse(text=destAccess),envir=.GlobalEnv)))))
	    eval(parse(text=paste(destAccess,"=factor(",destAccess,")",sep="")),envir=.GlobalEnv)
	  else
	    eval(parse(text=paste(destAccess,"=as.numeric(",destAccess,")",sep="")),envir=.GlobalEnv)
    }
    
    # NUMERIC TRANSFORM
    else {
    
      # The dot is interpreted as the source variable itself
      transExp = .$removeSpaces(transExp)
      transExp = gsub("[.]",sourceVar,transExp)

      eval(parse(text=paste(destAccess,"<-with(",dataname,",",transExp,")")),envir=.GlobalEnv)
    }
    
    # Update various variable lists
    if(svalue(.$currentData)==dataname) {
      .$updateVarList(h,...)
      .$updateWeightList(h,...)
      .$updateConstrFactor(h,...)
    }

    # Redisplay data (this also call updateTransf() via the implicit tab change
   .$showData(dataname)

    # Update the droplists
    svalue(.$currentFactor) = destVar
    svalue(.$newVarName)    = destVar
  },
  ### Update the transform fields when dataset is changed
  updateTransf = function(.,h,...) {

    displayed.grids = names(.$gridNotebook)
    if(length(displayed.grids) == 0) { 
     .$currentFactor[,]   = ""
	    svalue(.$toCompute) = ""
	   .$newVarName[,]      = ""
	    return()
    }
    
    # Get the current active data tab
    if(is.null(h$action))           dataname = displayed.grids[svalue(.$gridNotebook)]
    else if(h$action=="tabChange")  dataname = displayed.grids[h$pageno]

    if(!length(dataname)) return()

    .$currentFactor[,] <- .$newVarName[,] <- .$getVarList(dataname)[,1]
    
    # Just to trigger printCat()
    if(length(.$currentFactor[,])>1) {
      svalue(.$currentFactor,index=TRUE) = 2
      svalue(.$newVarName,index=TRUE) = 2
    }
  },
  ### Print factor categories in the recode/transform field
  printCat = function(.,h,...) {
  
    # Get the opened data tables
    displayed.grids = names(.$gridNotebook)

    # Empty transformation fields if none
    if(length(displayed.grids) == 0) { 
     .$currentFactor[,]   = ""
	    svalue(.$toCompute) = ""
	   .$newVarName[,]      = ""
	    return()
    }

    # Get the current active data tab
    dataname = displayed.grids[svalue(.$gridNotebook)]
    if(!length(dataname)) return()

    # Current source variable
    sourceVar = svalue(.$currentFactor)
    if(is.null(sourceVar)) return()
    
    # Display category names in transform field
    if(is.factor(eval(parse(text=paste(dataname,sourceVar,sep="$")),envir=.GlobalEnv))) {
      factorLevels = .$getFacLevelList(dataname,sourceVar)
      svalue(.$toCompute) = paste(factorLevels,collapse=",")
    }
    
    else svalue(.$toCompute) = ""
  },
  ### Save data grid under a new name
  saveGridAs = function(.,h,...) {

    filename = .$fileChoose("save")

    if(is.na(filename)) return()
    if(filename=="")    return()
    
    # Extract file extension
    name.ext = tolower(.$fileExtension(filename))
    if(name.ext != "csv") {
      gmessage(.$translate("Files may only be saved in CSV format."))
      return()
    }
    
    tabname =  .$baseName(filename)

    res = try(write.csv2(tabname,file=filename,row.names=FALSE))
    if(inherits(res,"try-error")) {
      gmessage("Write error while saving the file.")
      return()
    }
    
    # The working directory is implicitly set
    setwd(dirname(filename))
  },
  ### Save grid under its current name
  saveGrid = function(.,h,...) {

    displayed.grids = names(.$gridNotebook)
    dataname = displayed.grids[svalue(.$gridNotebook)]
    filename = paste(dataname,".csv",sep="")
    
    res = try(write.csv2(eval(parse(text=dataname),envir=.GlobalEnv),file=filename,row.names=FALSE))
    if(inherits(res,"try-error")) {
      gmessage(.$translate("Write error while saving the file."))
      return()
    }

    .$setStatus(.$translate("Status: CSV file successfully written to disk."))
  },
  ### Update data view in the datagrid tab
  showData = function(.,dataname) {

    if(!length(dataname) || (dataname == .$translate("No table"))) return()
    displayed.grids = names(.$gridNotebook)

    # Position of the new tab
    pos = length(displayed.grids) + 1
    
    # Data table already opened: Close it before updating
    if(dataname %in% displayed.grids) {
      pos = which(displayed.grids==dataname)
      svalue(.$gridNotebook) = pos
      dispose(.$gridNotebook)
    }
    
   .$setStatus(.$translate("Status: Data loading in progress. Please wait..."))
    
    # Add a new page to gridNotebook
    add(.$gridNotebook,tmp <- ggroup(),label = dataname,index=pos,expand=TRUE)
    addHandlerUnrealize(tmp,.$closeData)
    
    # A workaround, suggested by Tom Taverner, as DfEdit does not deal with ordered factors (to be corrected soon)
    # eval(parse(text=paste("for(var in which(sapply(",dataname,",is.ordered))) class(",dataname,"[,var])='factor'",sep="")),envir=.GlobalEnv)
    
    # Open grid
    add(tmp,  gdfedit(eval(parse(text=dataname),envir=.GlobalEnv),name=dataname), expand=TRUE)

    # Make it the active data
   .$currentDataName = dataname
    
    # Set to 2nd element of currentFactor just to trigger a change and call printCat
    svalue(.$currentFactor,index=TRUE) = 2
    svalue(.$newVarName,index=TRUE) = 2

    .$setStatus(.$translate("Status: Ready."))

  },
  ### Called when a data grid is closed
  closeData = function(.,h,...) {

    displayed.grids = names(.$gridNotebook)
    if(is.null(displayed.grids)) return()

    dataname = displayed.grids[svalue(.$gridNotebook)]
   .$currentDataName = dataname
    
    svalue(.$openGrid) = dataname
    svalue(.$currentData) = dataname
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS MODEL DEFINITION METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Dataset is changed from somewhere: Reset all dataframes droplists
  changeDataset = function(.,h,...) {
  
    # Authorized actions only
    if(is.null(h$action)) return()
    
    # From data grid selector
    if(h$action == "openGrid")  {
    
      if(.$debug) cat("Function: changeDataset, action: openGrid\n")
      
      dataName = svalue(.$openGrid)

      # Warning: a droplist may be temporarily set to NULL when changed
      if(is.null(dataName)) return()
      
      # No available table in working memory
      if(dataName == .$translate("No table")) return()
      
     .$currentDataName = dataName

      # Which tables are already displayed?
      displayed.grids = names(.$gridNotebook)
      
      # Avoid loops
      blockhandler(.$currentData,.$handler.ID['onChangeModelData'])

      # None: Then open the table
      if(is.null(displayed.grids)) {
       .$showData(dataName)
        return()
      }
      
      # Which one is visible?
      which.visible = displayed.grids[svalue(.$gridNotebook)]

      # Already selected: Useless to trigger all other handlers
      if(which.visible == dataName) return()

      # Put the data grid to the foreground if already opened
      if(dataName %in% displayed.grids) {
        svalue(.$gridNotebook) = which(displayed.grids == dataName)
        return()
      }
      
      # ... or open it
     .$showData(dataName)
           
      unblockhandler(.$currentData,.$handler.ID['onChangeModelData'])
    }
    
    # From data frame selector in the model tab
    else if(h$action == "changeModelData") {
    
      if(.$debug) cat("Function: changeDataset, action: changemodeldata\n")

      newName = svalue(.$currentData)

      # Warning: A droplist may be temporarily set to NULL when changed
      if(is.null(newName)) return()
          
     .$currentDataName = newName
     
      # This will trigger gridNotebook updating too
      blockhandler(.$gridNotebook,.$handler.ID['onTabChange'])
      svalue(.$openGrid) = newName
      unblockhandler(.$gridNotebook,.$handler.ID['onTabChange'])
    }
    
    else if(h$action == "tabChange") {
    
      if(.$debug) cat("Function: changeDataset, action: ontabchange\n")
      
      newName = names(.$gridNotebook)[h$pageno]
      if(is.null(newName)) return()
      
     .$currentDataName = newName
     
      # This will trigger the model panel reset too
      blockhandler(.$openGrid,.$handler.ID['onOpenGrid'])
      svalue(.$currentData) = newName
      unblockhandler(.$openGrid,.$handler.ID['onOpenGrid'])
    }
        
  },
  ### Get the variable names and types list in a data frame
  getVarList = function(.,dataname,type=TRUE) {

    if(.$debug) cat("Function: getVarList, File: ",dataname,"\n")
    
    if(!length(dataname) || (dataname == "") || (dataname == .$translate("No table"))) {
      vl = data.frame(Variables=.$translate("No variable"),Type="-")
      names(vl) = .$translate(names(vl))
      return(vl)
    }
    
    dataFrame = eval(parse(text=dataname),envir=.GlobalEnv)
    vList = colnames(dataFrame)
    varType = with(dataFrame,sapply(vList,function(x) is.numeric(eval(parse(text=x)))))

    if(type) {
      vl = cbind(Variables=vList,Type=ifelse(varType,.$translate("N"),.$translate("F")))
      return(vl)
    }
    else {
      vl = cbind(Variables=vList)
      names(vl) = .$translate(names(vl))
      return(vl)
    }
  },
  ### Get the vector of numeric variable names in the current dataframe
  getNumVarList = function(.,dataname) {

    if(!length(dataname)) return()
    if(dataname == "")    return()
    if(dataname == .$translate("No table")) return(NULL)

    dataFrame = eval(parse(text=dataname),envir=.GlobalEnv)
    vList = colnames(dataFrame)
    varType = with(dataFrame,sapply(vList,function(x) is.numeric(eval(parse(text=x)))))
    return(vList[varType])
  },
  ### Get the vector of categorical variable names in the current dataframe
  getCatVarList = function(.,dataname) {

    if(!length(dataname)) return()
    if(dataname == "")    return()
    if(dataname == .$translate("No table")) return(NULL)

    dataFrame = eval(parse(text=dataname),envir=.GlobalEnv)
    vList = colnames(dataFrame)
    varType = with(dataFrame,sapply(vList,function(x) is.numeric(eval(parse(text=x)))))
    return(vList[!varType])
  },
  ### Update the variable list
  updateVarList = function(.,h,...) {

    newName = svalue(.$currentData)
    if(is.null(newName)) return()
    
   .$currentDataName = newName
    if(.$currentDataName == .$translate("No table")) return()

    # Update main variable list
    variables = .$getVarList(.$currentDataName)
    if(!length(variables)) return()
   .$varList[,] = variables
    svalue(.$varList) = NULL
  },
  ### Update the weighting variable list
  updateWeightList = function(.,h,...) {

    if(.$currentDataName == .$translate("No table")) {
      svalue(.$weightList) = .$translate("No variable")
      return()
    }
  
    variables = .$getNumVarList(.$currentDataName)
    weightVar = svalue(.$weightList)
   .$weightList[] = c(.$translate("No variable"),variables)
    if(!weightVar %in% .$weightList[]) svalue(.$weightList) = .$translate("No variable")
    else svalue(.$weightList) = weightVar  
  },
  ### Update the numeric variables list on the graph tab (to be used as the x-axis)
  updateGraphNumVarList = function(.,h,...) {

    if(.$currentDataName == .$translate("No table")) {
      svalue(.$graphAxisX) = .$translate("Default")
      return()
    }
  
    variables = .$getNumVarList(.$currentDataName)
    xaxisVar = svalue(.$graphAxisX)
   .$graphAxisX[] = c(.$translate("Default"),variables)
    if(!xaxisVar %in% .$graphAxisX[]) svalue(.$graphAxisX) = .$translate("Default")
    else svalue(.$graphAxisX) = xaxisVar  
  },
  ### Update the constraint factor list
  updateConstrFactor = function(.,h,...) {

    if(.$currentDataName == .$translate("No table")) {
      svalue(.$structList) = .$translate("No factor")
      return()
    }
  
    variables = .$getCatVarList(.$currentDataName)
    structVar = svalue(.$structList)
   .$structList[] = c(.$translate("No factor"),.$translate("Constant"),variables)
    if(!structVar %in% .$structList[]) svalue(.$structList) = .$translate("No factor")
    else svalue(.$structList) = structVar  
  },
  ### Update the variable summary each time a variable name is clicked
  updateVarSummary = function(.,h,...) {

    currentVar = svalue(.$varList,drop=FALSE)
    varType    = currentVar[1,2]
    currentVar = currentVar[1,1]

    # No variable name selected
    if(!length(currentVar) || is.na(currentVar)) {
     .$varSummary[,] = data.frame(Attribute=.$translate("None"),Value=.$translate("None"),stringsAsFactors=FALSE)
      names(.$varSummary) = .$translate(names(.$varSummary))
      return()
    }
    
    varContent = eval(parse(text=paste(.$getDataName(),"[,'",currentVar,"']",sep="")),envir=.GlobalEnv)

    # A factor is selected
    if(varType == .$translate("F")) {
      vs = table(varContent)
      cn = c(names(vs),.$translate("Missing"),.$translate("Total"))
      vs = c(as.vector(vs),sum(is.na(varContent)),length(varContent))
     .$varSummary[,] = cbind(Attribute=cn,Value=vs)
      names(.$varSummary) = .$translate(names(.$varSummary))
    }
    
    # A numeric variable
    else if(varType == .$translate("N")) {
      Q = quantile(varContent,probs=c(0,.25,.5,.75,1),na.rm=TRUE)
     .$varSummary[,] = cbind(Attribute=.$translate(c("Minimum","1st quartile","Median","Mean","3rd quartile","Maximum","Std. dev.","Total","Missing")),
                             Value=round(c(Q[1],Q[2],Q[3],mean(varContent,na.rm=TRUE),Q[4],Q[5],sd(varContent,na.rm=TRUE),length(varContent),sum(is.na(varContent))),4))
      names(.$varSummary) = .$translate(names(.$varSummary))
    }
  },

  ### Update link function list upon distribution selection
  updateLink = function(.,h,...) {
   .$currentLinkList[] = linkLists[[svalue(.$distribList)]]
    svalue(.$currentLinkList) = .$currentLinkList[1]
  },
  ### Get the defined model name
  getModelName = function(.) {
    .$trim(svalue(.$modelName))
  },
  ### Get all model names
  getModelNames = function(.) {
  
    if(!length(.$models)) return("")
    sapply(.$models, function(m) m$getName())
  },
  ### Clear current model name field
  clearModelName = function(.,h,...) {
    svalue(.$modelName) = ""
  },
  ### Refill the model fields with a given model specs
  retrieveModel = function(.,h,...) {
  
    modname = .$getModelName()
    fitted.models = .$getModelNames()
    if(!(modname %in% fitted.models)) return()
    
    svalue(.$dvList)  = .$models[[modname]]$dvField
    svalue(.$fivList) = .$models[[modname]]$ivField
    svalue(distribList,index=TRUE) = .$models[[modname]]$family
    svalue(currentLinkList,index=TRUE) = .$models[[modname]]$link
    svalue(weightList) = ifelse(.$models[[modname]]$weights=="NULL",.$translate("No variable"),.$models[[modname]]$weights)
    svalue(subsetVar) = ifelse(.$models[[modname]]$subset=="NULL", "",.$models[[modname]]$subset)
    svalue(structList) = .$models[[modname]]$constrFactor
  },
  ### Update the list of model names when a new one is fitted
  updateModelNameList = function(.,h,...) {
    .$modelName[,] = .$getModelNames()
  },
  ### Add dependent variable
  updateDVField = function(.,h,...) {
    dv = .$trim(svalue(.$dvList))
    model = paste(svalue(.$varList),collapse=",")
    if(nchar(dv)) svalue(.$dvList) = paste(dv,", ",model,sep="")
    else          svalue(.$dvList) = model
  },
  ### Clear dependent variable field
  clearDVField = function(.,h,...) {
    svalue(.$dvList) = ""
  },
  ### Edit model field
  updateFIVField = function(.,h,...) {

    op = h$action
    fiv = .$trim(svalue(.$fivList))
    vl = svalue(.$varList)
    
    if(op == "-1") {
      svalue(.$fivList) = paste(fiv,"+0",sep="")
      return()
    }
    
    else if(op == "1") {
      if(nchar(fiv)) svalue(.$fivList) = paste(fiv,"+1",sep="")
      else           svalue(.$fivList) = "1"
      return()
    }
    
    else if(op == "()") {
      pat = svalue(.$fivList,drop=TRUE)
      repl = paste("(",pat,")",sep="")

      if(nchar(pat)) svalue(.$fivList) = sub(pat,repl,fiv,fixed=TRUE)
      else           svalue(.$fivList) = paste(fiv,"()",sep="")
      return()
    }
    
    else if(op == "offset") {

      # Is there selected text ?
      pat = svalue(.$fivList,drop=TRUE)
      if(nchar(pat)) {
        repl = paste("offset(",pat,")",sep="")
        svalue(.$fivList) = sub(pat,repl,fiv,fixed=TRUE)
        return()
      }

      # Any variable selected in the list ?
      pat = vl[1]
      if(length(pat)) {
        if(nchar(fiv)) svalue(.$fivList) = paste(fiv,"+offset(",pat,")",sep="")
        else           svalue(.$fivList) = paste(fiv,"offset(",pat,")",sep="")
        return()
      }

      svalue(.$fivList) = paste(fiv,"+offset()",sep="")
      return()
    }
    
    else if(op == "|") {
      svalue(.$fivList) = paste(fiv,"|",sep="")
      return()
    }
    
    else if(op == "(1|)") {
	    if(length(vl))  svalue(.$fivList) = paste(fiv,"+(1|",vl[1],")",sep="")
	    else            svalue(.$fivList) = paste(fiv,"+(1|.)",sep="")
      return()
    }
    
    else if(op == "(.|)") {
	    if(length(vl))  svalue(.$fivList) = paste(fiv,"+(",vl[1],"|.)",sep="")
	    else            svalue(.$fivList) = paste(fiv,"+(.|.)",sep="")
      return()
    }
    
    else if(op == "add") {
    
      # Only relevant when a variable is selected
      if(!nchar(vl[1])) return()
      
      # Is there selected text ?
      pat = svalue(.$fivList,drop=TRUE)
      if(nchar(pat)) {
      
        # Only replace the first occurrence of '.'
        if(.$trim(pat) == '.') svalue(.$fivList) = sub(pat,vl[1],fiv,fixed=TRUE)

        # Replace all occurrences of a variable name
        else                   svalue(.$fivList) = gsub(pat,vl[1],fiv,fixed=TRUE)
        
        return()
      }
      else {
	      svalue(.$fivList) = paste(fiv,vl[1],sep="")
        return()
      }
    }
    
    else if(op == "/") {
	    paste(fiv,"/",sep="")
      return()
    }
    
    else {

      # no variable selected: just print the operator
      if(length(vl)==0) { svalue(.$fivList) = paste(fiv,op,sep="")    ; return() }

      # one variable selected: append it with the operator
      if(length(vl)==1) { 
        if(nchar(fiv)) svalue(.$fivList) = paste(fiv,op,vl,sep="")
        else           svalue(.$fivList) = vl
        return() 
      }

      # At least two variables selected
      model = paste(vl,collapse=op)
      if(nchar(fiv)) svalue(.$fivList) = paste(fiv,"+",model,sep="")
      else           svalue(.$fivList) = model
    }
  },
  ### Extract random variable name from a (1|.) term
  extractRandomVar = function(.,h,...) {
    fterms = unlist(strsplit(.$removeSpaces(.$ivField)),"+",fixed=TRUE)
    which.rand = grep("|",fterms,fixed=TRUE)
    if(length(which.rand)>1) {
      gmessage(.$translate("Only one random factor is accepted."))
      return("error")
    }
    
    # grep("^\\(1\\|.*)$",fterms,value=T) # Look exactly for the pattern (1|.)
    randvar = .$removeParentheses(fterms[which.rand])
    randvar = unlist(strsplit(randvar,"|",fixed=T))
    if(randvar[1] != "1") {
      gmessage(.$translate("Only random intercept models are accepted."))
      return("error")
    }
    
    randvar[2]
  },
  ### Clear model field
  clearFIVField = function(.,h,...) {
    svalue(.$fivList) = ""
  },
  ### Clear the case subsetting field
  clearSelectField = function(.,h,...) {
    svalue(.$subsetVar) = ""
  },
  ### Clear all definition fields
  clearModelFields = function(.,h,...) {
    .$clearModelName(h,...)
    .$clearDVField(h,...)
    .$clearFIVField(h,...)
    .$clearSelectField(h,...)
  },
  ### Get the active R2STATS model
  getCurrentModel = function(.,h,...) {
    if(.$currentModelName == .$translate("No model")) return(.$translate("No model"))
    .$models[[.$currentModelName]]
  },
  ### Get a named R2STATS model
  getModelByName = function(.,modelname) {
    .$models[[modelname]]
  },
  ### Get variables
  getDVField = function(.) {
    .$trim(svalue(.$dvList))
  },
  getDV = function(.) {
    unlist(strsplit(.$strip.cbind(.$trim(.$getDVField())),","))
  },
  getIVField = function(.) {
    ivf = .$trim(svalue(.$fivList))
    if(ivf=="") ivf = "1"
    ivf
  },
  getIV = function(.) {
    all.vars(as.formula(paste("~",.$getIVField(),sep="")))
  },
  getModelVars = function(.) {
    all.vars(as.formula(.$getFormulaAsString()))
  },
  ### Get model formula
  getFormulaAsString = function(.) {
    dvs = .$getDV()
    dv  = .$getDVField()
    iv  = .$getIVField()
    
    if(length(dvs) == 0)    return("")
    else if(length(dvs)==1) return(paste(dv,"~",iv))
    else                    return(paste("cbind(",dv,")~",iv))
  },
  ### Get link index as specified in the interface
  getLink = function(.) {
    svalue(.$currentLinkList,index=TRUE)
  },
  ### Get distribution family index as specified in interface
  getFamily = function(.) {
    svalue(.$distribList,index=TRUE)
  },
  ### Get subset command as text
  getSubset = function(.) {
    # subset=NULL is the default in glm()
    subset = "NULL"
    if(nchar(svalue(.$subsetVar))) subset = svalue(.$subsetVar)
    subset
  },
  ### Get the name of weighting variable
  getWeights = function(.) {
    weights = svalue(.$weightList)
    if(weights == .$translate("No variable")) return("NULL")
    weights
  },
  ### Get constraint factor
  getConstrFactor = function(.) {
    svalue(.$structList)
  },
  ### Model identification and estimation
  run = function(.,h,...) {

    modname = .$getModelName()
    tableau = .$getDataName()

    ##------------------------------- Some checks about model specs
    
    # No model name provided
    if(modname %in% c("",.$translate("New"))) {
      gmessage(.$translate("Please give a name for this model."))
      return()
    }

    # Don't give the model the dataframe name
    if(modname == tableau) {
      gmessage(.$translate("You can't give the same name to both data and model."))
      return()
    }
	
    # This model name already exists
    if(modname %in% .$getModelNames()) {
      if(!gconfirm(.$translate("This model already exists.\nDo you want to replace it?"))) {
        return()
      }
    }
    
    # Check dependent variable field
    vd = .$getDVField()
    if(!nchar(vd)) {
      gmessage(.$translate("You must specify at least one dependent variable."))
      return()
    }

    # Check independent variable field
    vif = .$getIVField()
    if(!nchar(vif)) {
      vif = "1"
    }
    
    # Do variables exist?
    listvar = .$getVarList(tableau,type=FALSE)   
    if(!all(.$getModelVars() %in% listvar)) {
      gmessage(.$translate("Error: One of the variables does not exist."))
      return()
    }
    
    # Distribution misspecification
    vds   = .$getDV()
    ndv   =  length(vds)
    distr = .$getFamily()

    current.data = eval(parse(text=tableau),envir=.GlobalEnv)
    if( (length(ndv)==1) && is.factor(current.data[,vds]) && ( !(distr %in% c(2,6,7))) ) {
      gmessage(.$translate("This distribution is not suited\nfor a categorical variable."))
      return()
    }
    
   .$setStatus(.$translate("Status: Parameter estimation in progress..."))

    # Model definition
    vifs    = .$getIV()
    link    = .$getLink()
    
    subset  = .$getSubset()
    if(.$trim(subset) == "") subset = "NULL"
    
    weights = .$getWeights()
    if(weights == .$translate("No variable")) weights = "NULL"
    
    constrFactor = .$getConstrFactor()

    # Random effect model
    if(length(grep("|",vif,fixed=T))) {

      # Multinomial mixed model
      if(distr %in% 6:7) {
        randvar = .$extractRandomVar()
        if(randvar == "error") return()
        
        model = r2sCLMM$new(name=modname,class="clmm",func="clmm",dvField=vd,dv=vds,ivField=vif,iv=vifs,random=randvar,
                            data=tableau,family=distr,link=link,weights=weights,constrFactor=constrFactor,subset=subset)      
      }
      
      else {
        model = r2sGLMM$new(name=modname,class="mer",func="glmer",dvField=vd,dv=vds,ivField=vif,iv=vifs,data=tableau,
                             family=distr,link=link,weights=weights,constrFactor=constrFactor,subset=subset)
      }
	  }
     	   
    # Fixed effects model
    else {

      # Multivariate model
      if(ndv>1) {

        # Multivariate gaussian
        if(distr==1) {
          model = r2sMANOVA$new(name=modname,class="mlm",func="manova",dvField=vd,dv=vds,ivField=vif,iv=vifs,data=tableau,
                               family=distr,link=link,weights=weights,constrFactor=constrFactor,subset=subset)
	      }
	      
	      # Binomial with success/failure counts in two columns
	      else if(distr==2) {
	        if(ndv>2) {
	          gmessage(.$translate("For a binomial model, counts\nmust appear in two columns."))
	          return()
	        }
          model = r2sGLM$new(name=modname,class="glm",func="glm",dvField=vd,dv=vds,ivField=vif,iv=vifs,data=tableau,
                             family=distr,link=link,weights=weights,constrFactor=constrFactor,subset=subset)    
	      }
      }
      
      # Multinomial with multicolumn counts
      else if(distr %in% 6:7) {
        model = r2sCLM$new(name=modname,class="clm",func="clm",dvField=vd,dv=vds,ivField=vif,iv=vifs,data=tableau,
                            family=distr,link=link,weights=weights,constrFactor=constrFactor,subset=subset)
      }
        
      # Standard univariate GLM
      else {
        model = r2sGLM$new(name=modname,class="glm",func="glm",dvField=vd,dv=vds,ivField=vif,iv=vifs,data=tableau,
                           family=distr,link=link,weights=weights,constrFactor=constrFactor,subset=subset)
      }
    }
    
    # Model estimation
    res = model$estimate()
    if(inherits(res,"try-error")) {
     .$setStatus(.$translate("Status: Ready."))
      return()
    }
    
    # Save model in R2STATS model list
   .$models[[model$name]] = model
   .$updateModelNameList(h)
   
    # Make it the active model
   .$currentModelName = modname
    
    # Print results
    model$Summary()

    # Add to model lists (in graph and compare tabs)
   .$updateModelList(h)
       
    # Update list of available plots (if necessary)
   .$setStatus(.$translate("Status: Building plots..."))
    possiblePlots = model$getAvailablePlots()
    lastPlotType = svalue(.$plotType)
    if(any(.$plotType[,] != possiblePlots)) {
     .$plotType[,] = possiblePlots
      if(lastPlotType %in% possiblePlots) svalue(.$plotType) = lastPlotType
      else                                svalue(.$plotType,index=TRUE) = 1
    }

    # Highlight in plot list (this change automatically triggers the plotCurrentModel() handler)
    whichToPlot = which(.$graphModelList[,1] == modname)
    svalue(.$graphModelList,index=TRUE) = whichToPlot

    # Switch to result tab
   .$setStatus(.$translate("Status: Ready."))
    svalue(.$mainNotebook) = 4
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS RESULT PRINTING METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  clearResults = function(.,h,...) {
    svalue(.$results)=""
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS MODEL PLOTTING METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Generates group colors (borrowed from fBasics and RColorBrewer)
  getColors = function(.,n,name=c("Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3"),alpha=255) {

    Accent = rgb(c(127, 190, 253, 255, 56, 240, 191, 102),
                 c(201, 174, 192, 255, 108, 2, 91, 102), 
                 c(127, 212, 134, 153, 176, 127, 23, 102),maxColorValue = 255)
    Dark2 = rgb(c(27, 217, 117, 231, 102, 230, 166, 102), 
                c(158, 95, 112, 41, 166, 171, 118, 102),
                c(119, 2, 179, 138, 30, 2, 29, 102), maxColorValue = 255)
    Paired = rgb(c(166, 31, 178, 51, 251, 227, 253, 255, 202, 106, 255, 177),
                 c(206, 120, 223, 160, 154, 26, 191, 127, 178, 61, 255, 89),
                 c(227, 180, 138, 44, 153, 28, 111,   0, 214, 154, 153, 40), maxColorValue = 255)
    Pastel1 = rgb(c(179, 251, 204, 222, 254, 255, 229, 253, 242), 
                  c(205, 180, 235, 203, 217, 255, 216, 218, 242), 
                  c(227, 174, 197, 228, 166, 204, 189, 236, 242), maxColorValue = 255)
    Pastel2 = rgb(c(179, 253, 203, 244, 230, 255, 241, 204), 
                  c(226, 205, 213, 202, 245, 242, 226, 204),
                  c(205, 172, 232, 228, 201, 174, 204, 204), maxColorValue = 255)
    Set1 =    rgb(c(55, 228, 77, 152, 255, 255, 166, 247, 153), 
                  c(126, 26, 175, 78, 127, 255, 86, 129, 153),
                  c(184, 28, 74, 163, 0, 51, 40, 191, 153), maxColorValue = 255)
    Set2 =    rgb(c(102, 252, 141, 231, 166, 255, 229, 179),
                  c(194, 141, 160, 138, 216, 217, 196, 179),
                  c(165, 98, 203, 195, 84, 47, 148, 179), maxColorValue = 255)
    Set3 =    rgb(c(141, 255, 190, 251, 128, 253, 179, 252, 217, 188, 204, 255),
                  c(211, 255, 186, 128, 177, 180, 222, 205, 217, 128, 235, 237),
                  c(199, 179, 218, 114, 211, 98, 105, 229, 217, 189, 197, 111), maxColorValue = 255)
    name = match.arg(name)
    orig = eval(parse(text = name))
    
    if(n<10) return(orig[1:n])

    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0,1, ,length(orig))
    xg = seq(0,1,,n)
    for (k in 1:3) {
        hold = spline(x, rgb[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[,k] = round(hold)
    }
    
    rgb(temp[,1],temp[,2],temp[,3],alpha=alpha,maxColorValue = 255)
  },
  ### Set graphical parameters for a given model, depending upon the number of groups or classes
  setPlotParams = function(.,n=9) {

    fullColors   = .$getColors(n,"Set1")
    pastelColors = .$getColors(n,"Pastel1")

    trellis.par.set(    
      plot.symbol       = list(col = fullColors,fill = pastelColors,cex = rep(.8,n),font = rep(1,n),pch = rep(21,n),alpha = rep(1,n)),
      superpose.symbol  = list(col = fullColors,fill = pastelColors,cex = rep(.8,n),font = rep(1,n),pch = rep(21,n),alpha = rep(1,n)),
      plot.line         = list(alpha = rep(1,n),col = fullColors,lty = rep(1,n),lwd = rep(1,n)),
      superpose.line    = list(alpha = rep(1,n),col = fullColors,lty = rep(1,n),lwd = rep(1,n)),
      plot.polygon      = list(border = fullColors,col = fullColors),
      superpose.polygon = list(alpha = rep(1,n),col = fullColors,border="black",lty = rep(1,n),lwd = rep(1,n)),
      add.line          = list(alpha=rep(1,n),col=fullColors,lty=rep(1,n),lwd=rep(2,n)),
      add.text          = list(alpha=rep(1,n),cex=rep(1,n),col=fullColors,font=rep(1,n),lineheight=rep(1.2,n))
    )
  },
  ### Plot current model
  plotCurrentModel = function(.,h,...) {

    # No plot if no model exists
    if(!length(.$models)) return()

    # Current model is defined by the state of graphModelList
    selectedModel = svalue(.$graphModelList)
    
    # No plot if no model selected
    if(!length(selectedModel)) return()

    # Make it the active model
   .$currentModelName = svalue(.$graphModelList)
   
    # Set colors
    currentModel = .$models[[.$currentModelName]]
    if(length(currentModel$designFactors)) .$setPlotParams(nlevels(currentModel$groupLabels))
    
    # The model constructs its own plot...
    currentModel$Plot(h)
    
    # ... but R2STATS prints it on the device
    print(.$currentPlot)
    
    # Reset observation locator
    svalue(.$obsId) = ""
  },
  ### A blank plot when none is sensible
  emptyPlot = function(.) {
    xyplot(1:5~1:5,type="n",bty="n",scales=list(draw=F),xlab="",ylab="",
           panel = function(x,y,...) {
             panel.text(3,3,.$translate("No plot available."))
           })
  },
  ### Get desired plot type
  getPlotType = function(.) {
    svalue(.$plotType)
  },
  getPlotTypeAsIndex = function(.) {
    svalue(.$plotType,index=TRUE)
  },
  ### Get the selected group for plotting
  getSelectedGroup = function(.) {
    svalue(.$groupList)
  },
  ### Get the desired legend placement
  getLegendLocation = function(.) {
    c("none","right","left","top","bottom")[svalue(.$legendLoc,index=T)]  
  },
  ### Get the desired legend number of columns
  getLegendCols = function(.) {
    c(1:10)[svalue(.$legendCols,index=T)]  
  },
  ### Get limits on the x-axis
  getXLim = function(.) {
    
    xlim = NULL
    graph.xlim = .$trim(svalue(.$graphLimitsX))

    if(graph.xlim != "") {
      if(length(grep(",",graph.xlim))) {
        gmessage(.$translate("Please use the dot as the decimal separator\nand a blank space as a value separator."))
      }
      else {
        xlim = as.numeric(unlist(strsplit(graph.xlim,split=" ")))
        xlim = xlim[!is.na(xlim)]
      }
    }
    xlim
  },
  ### Get limits on the x-axis
  getYLim = function(.) {
    
    ylim = NULL
    graph.ylim = .$trim(svalue(.$graphLimitsY))

    if(graph.ylim != "") {
      if(length(grep(",",graph.ylim))) {
        gmessage(.$translate("Please use the dot as the decimal separator\nand a blank space as a value separator."))
      }
      else {
        ylim = as.numeric(unlist(strsplit(graph.ylim,split=" ")))
        ylim = ylim[!is.na(ylim)]
      }
    }
    ylim
  },
  ### Get the flag value for adding the Y=X line
  getAddLine01 = function(.) {
    svalue(.$addLine01)
  },
  ### Display observation id. on a click
  getXY = function(.,h,...) {
  
    # Only relevant for XY scatterplots
    if(is.null(.$currentPlot$panel.args[[1]]$y)) return()

    # Which obs. is closest to mouse pointer
    currentPlot.XY = data.frame(x=.$currentPlot$panel.args[[1]]$x,y=.$currentPlot$panel.args[[1]]$y)
    
    # Argh.. this returns (Inf,Inf) with lattice plots
    target = c(h$x,h$y)

    D2 = (currentPlot.XY$x - h$x)**2 + (currentPlot.XY$y - h$y)**2
    closest = which(D2 == min(D2))[1]

    # Display the observation number/name that is closest to the mouse pointer
    # svalue(.$obsId) = paste("Obs.",names(.$currentPlot$panel.args[[1]]$y)[closest])
  },
  ### Save plot in a PNG file
  savePlot = function(.,h,...) {
  
    fn = .$fileChoose("save")
    if(is.na(fn) || (fn == "")) return()
    
    fn = unlist(strsplit(fn,"\\."))[1]
    fn = paste(fn,".png",sep="")
    png(filename=fn)
   .$plotCurrentModel(h,...)
    dev.off()
  },
  ### Copy plot in clipboard
  copyPlot = function(.,h,...) {
    if(Sys.info()["sysname"] != "Windows") {
      gmessage(.$translate("This function is only available under Windows."))
    	return()
    }
    
    win.metafile()
   .$plotCurrentModel(h,...)
    dev.off()
    gmessage(.$translate("The graph has been saved to clipboard."))
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                           R2STATS MODEL COMPARISON METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Return the full list of all R2STATS models and their main attributes
  getModelList = function(.) {

    if(!length(.$models)) {
      df = data.frame(Name=.$translate("No model"),Formula=.$translate("No formula"),Constraint=.$translate("No factor"),Distribution=.$translate("Unspecified"),Link=.$translate("Unspecified"),stringsAsFactors=FALSE)
      names(df) = .$translate(names(df))
      return(df)
    }

    mnames = sapply(.$models, function(model) model$getName())
    formul = sapply(.$models, function(model) model$getFormula())
    constr = sapply(.$models, function(model) model$getConstrFactor())
    distr  = sapply(.$models, function(model) model$getFamilyAsIndex())
    distr  = names(.$linkLists)[distr]
    link   = sapply(.$models, function(model) .$linkLists[[model$getFamilyAsIndex()]][model$getLinkAsIndex()])

    df = cbind(Name=mnames,Formula=formul,Constraint=constr,Distribution=distr,Link=link)
    names(df) = .$translate(names(df))
    return(df)
  },
  ### Update model list in tab 5
  updateModelList = function(.,h,...) {

    ml = .$getModelList()
   .$modelList[,]      = ml
   .$graphModelList[,] = ml
  },
  ### Delete models in model list
  deleteModels = function(.,h,...) {
   
    # Vector of selected model names
    modelsToDelete = svalue(.$modelList)
    if(!length(modelsToDelete) || (modelsToDelete == .$translate("No model"))) {
      .$currentModelName = .$translate("No model")
      return()
    }
    
    # Remove from model list in the Compare tab
   .$models[modelsToDelete] = NULL
   
    # Remove from model list in the Model tab
   .$updateModelNameList(h,...)
    
    # Refresh model lists
   .$updateModelList(h,...)
  },
  ### Select all entries in R2STATS model list (tab 5)
  selectAllModels = function(.,h,...) {
    if(.$modelList[1,1]==.$translate("No model")) return()
    nmodels = nrow(as.data.frame(.$modelList[,]))
    svalue(.$modelList,index=TRUE) = 1:nmodels
  },
  ### Compare selected models
  compareModels = function(.,h,...) {

    # Names and indicies of selected models
    comp = svalue(.$modelList)
    idx  = svalue(.$modelList,index=TRUE)

    # Reorder by complexity
    Dfs = sapply(.$models[comp], function(m) m$df())
    k = order(Dfs)
    if(length(k)>1) {
      comp = comp[k]
      idx = idx[k]
    }

    # Model characteristics (as vectors or matrices)
    classes  = unique(lapply(.$models[comp], function(m) m$getClass()))
    dvs      = lapply(.$models[comp], function(m) m$getDV())
    distribs = sapply(.$models[comp], function(m) m$getFamily())
    sizes    = lapply(.$models[comp], function(m) length(m$getY()))
    
    # Do model have the same dependent variable?
    same.depvar = length(unique(dvs)) == 1
    if(!same.depvar) {
      gmessage(.$translate("The dependent variable is not the same across these models.\nThe comparison is meaningless."))
      return()
    }

    # Do models bear upon the same number of observations?
    same.size = length(unique(sizes)) == 1
    if(!same.size) {
      gmessage(.$translate("These models do not bear upon the same observation set."))
      return()
    }
    
    # If all models are from the same family, deviance is analyzed (AIC-BIC otherwise)
    same.distrib = length(unique(distribs)) == 1

    # Model statistics
    AICs     = sapply(.$models[comp], function(m) m$Aic())
    BICs     = sapply(.$models[comp], function(m) m$Bic())
    Expl     = sapply(.$models[comp], function(m) m$devExplained())
    logLiks  = sapply(.$models[comp], function(m) m$LogLik())
    Dfs      = sapply(.$models[comp], function(m) m$df())

    # Print titles
    add(.$results,.$translate("Model comparison and summary"),font.attr=c(style="normal",weights="bold",size="large",col="blue"))
    add(.$results,"")
    add(.$results,.$translate("Analysis of variance/deviance table"),font.attr=c(style="normal",weights="bold",col="black"))
    add(.$results,"")
    
    if( (length(classes)==1) && (classes == "glm")) {
    
      # If all models are gaussian or Gamma, then F tests are used (chi-square otherwise)
      all.gaussian = all(distribs == "gaussian")
      all.gamma    = all(distribs == "Gamma")

      # Several models selected: Compare them
     	if(length(comp)>1) {
     	
        if(same.distrib)  {

          # Use 'F' in case of gaussian or Gamma models, chi-square otherwise
          cmd = paste(".$tabdev = anova(",paste(".$models$",comp,"$Rmodel",collapse=",",sep=""),ifelse(all.gaussian || all.gamma, ",test='F')", ",test='Chisq')"))

          # Execute command
          eval(parse(text=cmd))
          
          # Adapt column headers
          if( (all.gaussian) || (all.gamma) )
            attr(.$tabdev,"names")=.$translate(c("Resid. Df", "Resid. Dev","Diff. Df","LR","F","Pr(>F)"))
          else 
            attr(.$tabdev,"names")=.$translate(c("Resid. Df", "Resid. Dev","Diff. Df","LR","Pr(>Chi2)"))

          # Model names appear as row names
          attr(.$tabdev,"row.names") = comp
          
          # Remind me of model formulae
          attr(.$tabdev,"heading") = c(capture.output(.$modelList[idx,2:5]),"")

          # Add AIC and BIC when several models are compared (anova also works for a single one)
          if(length(comp)>1)  {
            # TODO: cet ajout de colonnes modifie l'affichage des décimales (!?)
            .$tabdev[["AIC"]]      = AICs
            .$tabdev[["BIC"]]      = BICs
            .$tabdev[["Expl.(%)"]] = Expl
          }
        }

        # Just AIC and BIC in case of different distributions or non-nested models
        else  {
        
          .$tabdev = data.frame(Distribution=.$modelList[idx,4],Lien=.$modelList[idx,5],AIC=AICs,BIC=BICs)

          # Model names appear as row names
          attr(.$tabdev,"row.names") = comp
          
          # Remind me of model formulae
          attr(.$tabdev,"heading") = c(paste(.$modelList[idx,1],.$modelList[idx,2],sep=" : ",collapse="\n"),"")
        }
      }
      
      # Only one model selected: Print summary
     	if(length(k)==1) {
     	
     	  # Gaussian models: Standard anova table
     	  if(distribs == "gaussian") {
       	    cmd = paste(".$tabdev = summary(aov(.$models$",comp,"$Rmodel))[[1]]",sep="")
      	    eval(parse(text=cmd))
            attr(.$tabdev,"names")=.$translate(c("Ddl.","SC","CM","F","Pr(>F)"))
     	  }
     	  
     	  # Non-gaussian models: Compare to the null
     	  else {
       	    cmd = paste(".$tabdev = anova(.$models$",comp,"$Rmodel,test='Chisq')",sep="")
      	    eval(parse(text=cmd))
      	  
      	    # For some reason, columns are not in the same order with a single model!
           .$tabdev = .$tabdev[2:1,c(3,4,1,2,5)]
            attr(.$tabdev,"names")=.$translate(c("Resid. Df", "Resid. Dev","Diff. Df","LR","Pr(>Chi2)"))
     	  }

        # Remind me of model formulae
        attr(.$tabdev,"heading") = c(paste(.$modelList[idx,1],.$modelList[idx,2],sep=" : ",collapse="\n"),"")
     	}
	  }
	  
	  # Compare GLMMs
	  else if( (length(classes)==1) && (classes == "mer") ) {
	  
	    # anova() for a single model not implemented
     	if(length(comp)==1) {
     	  gmessage(.$translate("You must select at least two models."))
     	  return()
      }
      
      chisq = 2 * pmax(0, c(NA, diff(logLiks)))
      dfChisq = c(NA, diff(Dfs))
      pchi2 = pchisq(chisq, dfChisq, lower = FALSE)
     .$tabdev = data.frame(logL=logLiks,Df=Dfs,Chi2=chisq,"Diff. Df"=dfChisq,Prob.=pchi2,AIC=AICs,BIC=BICs)
      names(.$tabdev) = .$translate(names(.$tabdev))

      # Model names appear as row names
      attr(.$tabdev,"row.names") = comp
      
      # Remind me of model formulae
      attr(.$tabdev,"heading") = c(paste(.$modelList[idx,1],.$modelList[idx,2],sep=" : ",collapse="\n"),"")
      
      # Just to get figures nicely printed
      class(.$tabdev) = c("anova",class(.$tabdev))
    }

    # Compare different models from different classes: use AIC and BIC
	  else {
	  
     .$tabdev = data.frame(Distribution=.$modelList[idx,4],Lien=.$modelList[idx,5],AIC=AICs,BIC=BICs,stringsAsFactors=FALSE)

      # Model names appear as row names
      attr(.$tabdev,"row.names") = comp
      
      # Remind me of model formulae
      heading = c(paste(.$modelList[idx,1],.$modelList[idx,2],sep=" : "))

      # print.data.frame() does not print headings
      add(.$results,heading,font.attr=c(family="monospace",size="medium"))
      add(.$results,"")
	  }
	  
    # Output
    add(.$results,capture.output(.$tabdev),font.attr=c(family="monospace",size="medium"))
    add(.$results,"")
    
    # Automatically switch to result page
    svalue(.$mainNotebook) = 4
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS GENERIC METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  getVersion = function(.) {
    return(.$version)
  },
  setStatus = function(.,text) {
    svalue(.$status)
    svalue(.$status) = text
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS TOOLS
  #
  #------------------------------------------------------------------------------------------------------------------------
  loadModels = function(.,h,...) {
    load(.$fileChoose("open"),envir=.)
   .$updateModelNameList(h)
   .$updateModelList(h)
  },
  saveModels = function(.,h,...) {
    save(models,envir=.,file=.$fileChoose("open"))
  },
  editOptions = function(.,h,...) {
    gmessage(.$translate("Function not yet available."))
  },
  aboutR2stats = function(.,h,...) {
    aboutMessage = gbasicdialog(title="About...",do.buttons=FALSE)
    messageFrame = gframe(cont=aboutMessage,horizontal=FALSE)
    add(messageFrame,glabel("<span foreground='blue' size='x-large' weight='ultrabold'>R2STATS</span>",markup=TRUE))
    add(messageFrame,glabel(paste("version",.$version)))
    add(messageFrame,glabel("\n   A GTK GUI for fitting GLM and GLMM in R   \n"))
    add(messageFrame,glabel("<b>Yvonnick Noel</b>",markup=TRUE))
    add(messageFrame,glabel("University of Brittany at Rennes, France"))
    add(messageFrame,glabel("<i>yvonnick.noel@uhb.fr</i>\n",markup=TRUE))
    visible(aboutMessage,set=TRUE)
  },
  getAllTables = function(.) {
    l = ls(envir=.GlobalEnv)
    classes = unlist(sapply(l,function(x) class(get(x))))
    data.tabs = which(classes %in% c("data.frame","matrix"))
    names(classes)[data.tabs]
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS STRING MANAGEMENT METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Remove all spaces from a string
  removeSpaces = function(.,x) {
    gsub("[ ]","",x)
  },
  ### Trim leading and trailing spaces from a string
  trim = function(.,x) {

    # Replace all multiple spaces with a single space
    x <- gsub("[ ]+", " ", x)
    # Remove all trailing spaces
    x <- gsub("[ ]+$", "", x)
    # Remove all leading spaces
    x <- gsub("^[ ]+", "", x)
    x
  },
  ### Get file base name
  getBaseName = function(.,x) {
    sub("(?x) # allow embedded comments
         (.+) # match and remember at least one arbitrary character
         [.] # match a dot
         [^.]+ # match at least one non-dot character
         $", # end of string anchor
         "\\1",basename(x), perl=TRUE)
  },
  ### Get file extension
  getFileExtension = function(.,x) {
    sub(".*\\.", "", x)
  },
  ### Remove cbind() from a string
  strip.cbind = function(.,l) {
    l = sapply(l,function(x) sub("cbind","",x))
    .$removeParentheses(l)
  },
  ### Remove parentheses from an expression
  removeParentheses = function(.,l) {
    l = sapply(l,function(x) sub("\\(","",x))
    sapply(l,function(x) sub("\\)","",x))  
  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-R2STATS")
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 R2STATS SLOTS
  #
  #------------------------------------------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #------------------------------------------------------------------------------------------------------------------------
  #----- General slots
  version                    = "??",
  mainWindow                 = NULL,
  mainWindowTitle            = "R2STATS",
  menu                       = NULL,
  mainNotebook               = NULL,
  status                     = "",
  debug                      = FALSE,
  handler.ID                 = list(),
  #----- Data slots
  dataUrl                    = NULL,
  hasHeader                  = NULL,
  hasRowNames                = NULL,
  libList                    = NULL,
  libDataList                = NULL,
  loadDataBut                = NULL,
  #----- Grid slots
  currentFactor              = NULL,
  toCompute                  = NULL,
  newVarName                 = NULL,
  gridNotebook               = NULL,
  openGrid                   = NULL,
  #----- Model slots
  models                     = list(),
  panedGroup                 = NULL,
  currentData                = NULL,
  currentDataName            = "",
  varList                    = NULL,
  varSummary                 = NULL,
  modelName                  = NULL,                  # Name of the model currently defined in the model tab
  currentModelName           = "Aucun",               # Name of the model currently active (from the model or graph tabs)
  dvList                     = NULL,
  fivList                    = NULL,
  distribList                = NULL,
  linkLists                  = NULL,
  currentLinkList            = NULL,
  weightList                 = NULL,
  structList                 = NULL,
  subsetVar                  = "NULL",                # Must be a string
  #----- Result slots
  results                    = NULL,
  #----- Graphics slots
  currentPlot                = NULL,                  # Current trellis/lattice plot
  plotType                   = NULL,
  graphModelList             = NULL,
  legendLoc                  = NULL,
  legendCols                 = NULL,
  groupList                  = NULL,
  graphAxisX                 = NULL,
  graphLimitsX               = NULL,
  graphLimitsY               = NULL,
  addData                    = NULL,
  addModel                   = NULL,
  addCondMeans               = NULL,
  addRefLine                 = NULL,
  addGrid                    = NULL,
  addNoise                   = NULL,
  addSmooth                  = NULL,
  addRandCurves              = NULL,
  obsId                      = NULL,
  #----- Model comparison slots
  modelList                  = NULL,
  tabdev                     = NULL
  #----- Option setting
)

# Main
R2STATS = function() {

  # Create R2STATS main window
  r2stats$create()
  
  # Show R2STATS interface
  r2stats$show()
}

