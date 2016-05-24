# ---------------------------------------
# Author: Daniel Schopfhauser
#         Vienna University of Technology
# ---------------------------------------



#' GUI for Visualization and Imputation of Missing Values
#' 
#' Graphical user interface for visualization and imputation of missing values.
#' 
#' Details about handling survey objects follow soon.
#' 
#' @param startupObject Object loaded at the start of the GUI
#' @author Daniel Schopfhauser
#' @keywords multivariate hplot
#' @export VIMGUI
VIMGUI <- function(startupObject=NULL){
  #fixate underlying GUI-Toolkit for gWidgets to GTK
  #as there are some specific used features
  options("guiToolkit"="RGtk2")
  
  ####
  #HANDLER FUNCTIONS
  
  #handler for notebook
  #called when changing the main tab (data, imputation,...)
  mainNotebook.handler <- function(h,...){
    updatePanels(pageno=h$pageno)
  }
  
  
  #handler for the plot notebook 
  #called when changing the tabs for the different plots
  imputationPlotHandler <- function(h,...){
    makeImputationPlot(h$pageno)
  }
  
  #main method for drawing the plot graphics
  #is called in cases where the plot graphic needs a redraw
  #like switching to the plot tab, switching the plot type
  #or change some parameters
  makeImputationPlot <- function(index = svalue(impVis.plotBook), savePlot=FALSE){
	#select which dataset is plotted (with or without imputed values) depending on the 
	#users choice in the plot tab
    if (svalue(impVis.plotImputed) == "original"){
      plotData <- getVm("activeDataSetOriginal")
	  #use no delimiter in case of not imputed data to prevent unnecessary warnings
      delimiter <- NULL
    }
    else{
      plotData <- getVm("activeDataSetImputed")
      delimiter <- "_imp"
    }
	
	#create plot depending on which tab is active
	#index is normally the index of the plot tab
    if (index == 1){
			#use the buffered plot function to create a aggr-plot
      bufferedPlot(aggr(plotData, bars=svalue(impVis.aggr.bars), numbers = svalue(impVis.aggr.numbers),
                        prop = svalue(impVis.aggr.prop), combined = svalue(impVis.aggr.combined),
                        only.miss = svalue(impVis.aggr.only.miss), sortVars = svalue(impVis.aggr.sortVars),
                        sortCombs = svalue(impVis.aggr.sortCombs),
                        col = getVm("plotColors"), delimiter=delimiter
                        #,weighted=svalue(impVis.aggr.weighted)
                        ), savePlot=savePlot)
    }
    else if (index == 2){
			#use the buffered plot function to create a barMiss-plot
      bufferedPlot(barMiss(plotData, pos = svalue(impVis.barMiss.pos, index=TRUE), 
              selection = svalue(impVis.barMiss.selection),
              only.miss = svalue(impVis.barMiss.only.miss),
              col = getVm("plotColors"),interactive=FALSE,
                           delimiter=delimiter
                           #,weighted=svalue(impVis.barMiss.weighted)
                           ), savePlot=savePlot)
    }
    else if (index == 3){
			#use the buffered plot function to create a histMiss-plot
			#select if breaks represents the number of breaks or a algorithm
      breaks <- svalue(impVis.histMiss.breaks)
      if (isNumber(breaks) == TRUE){
        breaks <- as.numeric(breaks)
      }
      bufferedPlot(histMiss(plotData, pos = svalue(impVis.histMiss.pos, index=TRUE),
               selection = svalue(impVis.histMiss.selection),
               breaks = breaks, right = svalue(impVis.histMiss.right),
               only.miss = svalue(impVis.histMiss.only.miss),
               col = getVm("plotColors"),interactive=FALSE,
                            delimiter=delimiter
                            #,weighted=svalue(impVis.histMiss.weighted)
                            ), savePlot=savePlot)
    }
    else if (index == 4){
			#use the buffered plot function to create a marginmatrix-plot
      drawNames <- as.character(svalue(impVis.marginMatrix.plotvars))
			#remove the delimiter variables from the actual variable selection
      impNames <- intersect(sapply(drawNames, FUN=function(s) paste(s,"_imp", sep="")),
                            names(plotData))
      drawData <- plotData[,c(drawNames, impNames)]
      bufferedPlot(marginmatrix(drawData, delimiter=delimiter,
                   col = getVm("plotColors"), alpha = getVm("plotAlpha")), savePlot=savePlot)
    }
    else if (index == 5){
			#use the buffered plot function to create a scattmatrixMiss-plot
      plotvars <- NULL
			#prevent error if no variables or highlights are selected
      v <- as.character(svalue(impVis.scattmatrixMiss.plotvars))
      if (length(v) > 0){
        plotvars <- v
      }
      highlight <- NULL
      v <- as.character(svalue(impVis.scattmatrixMiss.highlight))
      if (length(v) > 0){
        highlight <- v
      }
      bufferedPlot(scattmatrixMiss(plotData,
																	col = getVm("plotColors"), alpha = getVm("plotAlpha"),interactive=FALSE,
                                  delimiter=delimiter, plotvars=plotvars, highlight=highlight,
                                  selection=svalue(impVis.scattmatrixMiss.selection),
                                  diagonal=svalue(impVis.scattmatrixMiss.diagonal)
                                  #,weighted=svalue(impVis.scattmatrixMiss.weighted)
                                   ), savePlot=savePlot)
    }
    else if (index == 6){
			#use the buffered plot function to create a mosaicMiss-plot
			#prevent error if no variables or highlights are selected
      highlight <- as.character(svalue(impVis.mosaicMiss.highlight))
      if (length(highlight) < 1){
        highlight <- NULL
      }
      plotvars <- as.character(svalue(impVis.mosaicMiss.plotvars))
      if (length(plotvars) < 1){
        plotvars <- NULL
      }
      bufferedPlot(mosaicMiss(plotData, highlight=highlight, plotvars=plotvars, 
									selection=svalue(impVis.mosaicMiss.selection),
									col = getVm("plotColors"),
									delimiter=delimiter
                  #,weighted=svalue(impVis.mosaicMiss.weighted)
                  ), savePlot=savePlot)
    }
    else if (index == 7){
			#use the buffered plot function to create a parcoordMiss-plot
			#prevent error if no variables or highlights are selected
      highlight <- as.character(svalue(impVis.parcoordMiss.highlight))
      if (length(highlight) < 1){
        highlight <- NULL
      }
      plotvars <-as.character(svalue(impVis.parcoordMiss.plotvars))
      if (length(plotvars) < 1){
        plotvars <- NULL
      }
      bufferedPlot(parcoordMiss(plotData, highlight=highlight, plotvars=plotvars, 
                   selection=svalue(impVis.parcoordMiss.selection),
                   plotNA = svalue(impVis.parcoordMiss.plotNA),
                   col = getVm("plotColors"), alpha = getVm("plotAlpha"),interactive=FALSE,
                   delimiter=delimiter
                                #,weighted=svalue(impVis.parcoordMiss.weighted)
                                ), savePlot=savePlot)
    }
    else if (index == 8){
			#use the buffered plot function to create a pbox-plot
      bufferedPlot(pbox(plotData, pos=svalue(impVis.pbox.pos, index=TRUE),
           selection = svalue(impVis.pbox.selection),
           numbers = svalue(impVis.pbox.numbers),
           col = getVm("plotColors"),interactive=FALSE,
                        delimiter=delimiter
                        #,weighted=svalue(impVis.pbox.weighted)
                        ), savePlot=savePlot)
    }
    else if (index == 9){
      #use the buffered plot function to create a matrixplot-plot
      suppressWarnings({
        bufferedPlot(matrixplot(plotData, sortby=svalue(impVis.matrixplot.sortby, index=TRUE),
                                col = getVm("plotColors"),interactive=FALSE,
                                delimiter=delimiter
                                #,weighted=svalue(impVis.matrixplot.weighted)
                                ), savePlot=savePlot)
      })
    }
  }
  
  #called after using the contextual menu items of the plot window
  #uses the regular plot function but saves the created content 
	#to a permanent file on the hard-drive
  savePlotToFile <- function(FileFormat="PNG"){
		#extract the size of the plot inside the window
		#mostly to preserve the proportion of width and height
    tgtk <- getToolkitWidget(impVis.plot)
    a <- tgtk$getAllocation()
		#choose the file format the user wants
    if (FileFormat=="PNG"){
			#opens a file save dialog for selected file format
      location <- gfile("Save as PNG", type="save", filter=list("PNG","png"))
			#user input was OK
      if (!is.na(location)){
				#substitute missing file-extension
        if (!endsWithText(location,c("PNG","png"))){
          location <- paste(location,".png",sep="")
        }
				#use Cairo-function for saving 
        dev <- CairoPNG(filename=location, width=a$allocation$width, height=a$allocation$height)
        makeImputationPlot(savePlot=TRUE)
        dev.off(dev)
      }
    }
    else if (FileFormat=="PDF"){
			#opens a file save dialog for selected file format
      location <- gfile("Save as PDF", type="save", filter=list("PDF","pdf"))
			#user input was OK
      if (!is.na(location)){
				#pdf needs only proportion of weight and height not pixel-size
				#but absolute value of size changes font size
        w <- a$allocation$width / max(a$allocation$width,a$allocation$height)
        h <- a$allocation$height / max(a$allocation$width,a$allocation$height)
				#substitute missing file-extension
        if (!endsWithText(location,c("PDF","pdf"))){
          location <- paste(location,".pdf",sep="")
        }
        dev <- CairoPDF(file=location, width=w*12, height=h*12)
        makeImputationPlot(savePlot=TRUE)
        dev.off(dev)
      }
    }
    else if (FileFormat=="PS"){
			#opens a file save dialog for selected file format
      location <- gfile("Save as PS", type="save", filter=list("PS","ps"))
			#user input was OK
      if (!is.na(location)){
				#substitute missing file-extension
        if (!endsWithText(location,c("PS","ps"))){
          location <- paste(location,".ps",sep="")
        }
        dev <- CairoPS(file=location)
        makeImputationPlot(savePlot=TRUE)
        dev.off(dev)
      }
    }
    else if (FileFormat=="JPEG"){
			#opens a file save dialog for selected file format
      location <- gfile("Save as JPEG", type="save", filter=list("JPEG","jpeg","jpg","JPG"))
			#user input was OK
      if (!is.na(location)){
				#substitute missing file-extension
        if (!endsWithText(location,c("JPEG","jpeg","jpg","JPG"))){
          location <- paste(location,".jpeg",sep="")
        }
        dev <- CairoJPEG(filename=location, width=a$allocation$width, height=a$allocation$height)
        makeImputationPlot(savePlot=TRUE)
        dev.off(dev)
      }
    }
    else if (FileFormat=="SVG"){
			#opens a file save dialog for selected file format
      location <- gfile("Save as SVG", type="save", filter=list("SVG","svg"))
			#user input was OK
      if (!is.na(location)){
				#svg needs only proportion of weight and height not pixel-size
				#but absolute value of size changes font size
        w <- a$allocation$width / max(a$allocation$width,a$allocation$height)
        h <- a$allocation$height / max(a$allocation$width,a$allocation$height)
				#substitute missing file-extension
        if (!endsWithText(location,c("SVG","svg"))){
          location <- paste(location,".svg",sep="")
        }
        dev <- CairoSVG(file=location, width=w*12, height=h*12)
        makeImputationPlot(savePlot=TRUE)
        dev.off(dev)
      }
    }
  }
  
  #sets a new dataset as active i.e. sets the activeDataSetOriginal and
  #activeDataSetImputed Variables
  #after that initializes the GUI with the new values (meaning puts values in
	#different tables and widgets, resets settings, ...)
	#firstPage 	...  	resets the main notebook
	#prepare 		...		used after setting dataset after using the prepare command
	#									doesn't remove the old dataset to enable usage of undo
	#loadScript ...		R-source used to load variable, first line in script dialog
	#parent			...		parent for the small dialog windows
	#adjustTypes...		should adjusting of types be performed
  #deleteScript     should the previous accumulated scripts be deleted
  setActiveDataset <- function(x, firstPage=TRUE, prepare=FALSE, loadScript="", parent=NULL,
                               adjustTypes=TRUE, deleteScript=TRUE){
    if (firstPage){
      svalue(mainNotebook) <- 1
    }
    if(adjustTypes) x <- adjustTypesDialog(x)
		#open loading window and save its ID for later
    putVm("loadingWindowID",loadingWindow(parent=parent))
		#reset different settings
    putVm("activeDataSetOriginal", x)
    putVm("activeDataSetImputed", NULL)
    enabled(menu.Undo) <- FALSE
    #if called after preparing data keep the old dataset for possible undo 
    if (prepare == FALSE){
      putVm("undoDataSetOriginal", NULL)
    }
    
    #delete script history
    if (deleteScript== TRUE){
      putVm("ScriptHistory", list(loadScript, "originaldataset <- activedataset"))
    }
    
    initPanels()
    updatePanels(firstTime=TRUE)
		#destroy the loading window after loading is done
    dispose(getVm("loadingWindowID"))
  }
  
  #loads a R-dataset from a file
	#used after clicking on the corresponding menu entry
  loadDataSet <- function(...) {
		#open file dialog
    xname <- gfile("Select file to load", parent=window, type="open" ,filter=list("R-Data"=list(patterns=c("*.rda", "*.RData","*.RDA","*.rdata","*.RDATA")), "All files" = list(patterns = c("*"))))
    #is valid name load to global name space and open select dataset window
		if (is.na(xname) == FALSE) {
      putVm("importFilename",xname)
      load(xname, envir=.GlobalEnv)
      setDataSet()
    }
  }
  

  
  #lets the user select a dataset from the global environment to load it into the application
	#used after clicking on the corresponding menu entry
  setDataSet <- function(...) {
		#get variable names in global environment
    vardt <- ls(envir = .GlobalEnv, all.names=TRUE)
		vards <- character(0)
		#filter to only use dataframes and surveys
		if (length(vardt) != 0){
			vards <- names(which(sapply(vardt, function(.x) is.data.frame(get(.x)))))
			vards <- c(vards,names(which(sapply(vardt, function(.x) is.survey(get(.x))))))
		}
		#if there aren't any usable objects inform the user which error
    if( length(vards)==0 ) {
      gmessage("No datasets loaded.", title="Information", icon="warning")
    } else {
			#present the user with possible objects
      gbasicdialog(title="Choose Dataset",
                   x<-gdroplist(vards), parent=NULL,
                    handler=function(x, ...) {
											#set new dataset after confirmation by user
                      setActiveDataset(get(svalue(x$obj)), parent=mainWindow,
                                       , loadScript=paste("activedataset <- ",svalue(x$obj)))
                    })
    }
  }
  
  #lets the user save the current active dataset to a file
	#called after clicking the corresponding menu item
  saveToFile <- function(...) {
		#there is no imputed data -> save original
    if( is.null(getVm("activeDataSetImputed")) == FALSE) {
      xname <- gfile("Choose a file to save the Dataset", type="save", parent=window)
      if( xname != "" ) {
        data <- getVm("activeDataSetImputed")
        save(data, file=paste(xname,".RData", sep=""))
      }
			#there is imputed data -> save imputed data
    } else if (is.null(getVm("activeDataSetOriginal")) == FALSE){
      xname <- gfile("Choose a file to save the Dataset", type="save", parent=window)
      if( xname != "" ) {
        data <- getVm("activeDataSetOriginal")
        save(data, file=paste(xname,".RData", sep=""))
      }
			#catch possible error situation by not allowing to save empty data
    } else {
      gmessage("No active Dataset found.", title="Information", icon="warning",
               parent=window)
    }
  
  }
  
	#let the user save the currently active dataset to a variable in the global environment
	#called after clicking the corresponding menu entry
  saveToVariable <- function(...){
		#find and filter already exising variable names
    #get variable names in global environment
    vardt <- ls(envir = .GlobalEnv, all.names=TRUE)
    vards <- character(0)
    #filter to only use dataframes and surveys
    if (length(vardt) != 0){
      vards <- names(which(sapply(vardt, function(.x) is.data.frame(get(.x)))))
      vards <- c(vards,names(which(sapply(vardt, function(.x) is.survey(get(.x))))))
    }
    else{
      vards <- " "
    }
		#build small window with a combobox and two buttons for accept and discard
    saveVariable.window <- gwindow("Save Dataset to Variable", width=300, height=75)
    gg <- glayout(container=saveVariable.window)
    saveVariable.list <- gcombobox(c("",vards), editable=TRUE)
    saveVariable.accept <- gbutton(" Accept")
    saveVariable.discard <- gbutton(" Discard")
    gg[1,1:2, expand=TRUE] <- saveVariable.list
    gg[2,1] <- saveVariable.discard 
    gg[2,2] <- saveVariable.accept
		#destroy window when clicking discard button
    addHandlerClicked(saveVariable.discard, handler=function(h,...){
      dispose(saveVariable.window)
    })
		#accept button handler
    addHandlerClicked(saveVariable.accept, handler=function(h,...){
			#get save-to name and original data
      varName <- svalue(saveVariable.list)
      dataObject <- getVm("activeDataSetOriginal")
			#is imputed data exists ask user hwo to proceed
      if (is.null(getVm("activeDataSetImputed")) == FALSE){
        w <- gconfirm("Do you want to use the imputed values?", title="Imputed Values", icon="question")
				#user wants to save imputed data
        if (w == TRUE){
          dataObject <- getVm("activeDataSetImputed")
					#is dataset is survey remove the delimtier variables, as they are possible interfering 
					#which the survey methods
          if (is.survey(dataObject)){
            dataobject <- dataObject$variables
            dataobject <- dataobject[,grep("_imp", colnames(dataobject), invert=TRUE)]
            dataObject$variables <- dataobject
          }
        }
      }
      #if user selected a name which already existed in the global environment
			#ask the user if he wants to override it
      if( exists(varName, envir=.GlobalEnv) ) {
                gconfirm("Variable already exists, do you want to replace it?",
                    title="Information",
                    handler=function(h, ...) {
                      en <- as.environment(1)
                      assign(varName, dataObject, envir=en)
                      dispose(saveVariable.window)
                      })
              }
      else{
        en <- as.environment(1)
        assign(varName, dataObject, envir=en)
        dispose(saveVariable.window)
      }
      
    })
  }
  
  #a small dialog which a "Please Wait"-Text.
	#always called if something needs a few seconds to build up (loading data, init widgets,...)
  WaitingDialog <- function(parent, text="<b><big>Importing Data, Please Wait!</big></b>", 
                            header="Importing!", Parent=NULL){
    window <- gwindow(header, parent=Parent, width=100, height=50)
    glabel(text, markup=TRUE,container=window)
    return(window)
  }
  
  #window for importing CSV files with a interactive preview 
	#called after clicking the corresponding menu entry
  importCSV <- function(...){
		#init different window options and widgets
    importDialog <- gwindow("Import CSV", width=400, height=800)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    putVm("changedTypes", FALSE)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    csvfilename <- gedit()
    enabled(csvfilename) <- FALSE
		#handler for the file selection button handler
		#lets the user select a file and uses a simple heuristic for choosing default parameters
    buttonHandler <- function(...){
      gfile(text = "Open CSV File", type = "open", 
            filter=list("CSV files"=list(patterns=c("*.csv", "*.CSV")), "All files" = list(patterns = c("*"))),
            handler=function(h,...){
              svalue(csvfilename) <- h$file
							#reads the beginning of the CSV file for preview and parameter estimation
              tryCatch({
                fl <- readLines(svalue(csvfilename), n=2)
                comma <- sapply(strsplit(as.character(fl[1]), ","), length)
                semicolon <- sapply(strsplit(as.character(fl[1]), ";"), length)
                dot <- sum(sapply(strsplit(as.character(fl[2]), "."), length))
                comma2 <- sum(sapply(strsplit(as.character(fl[2]), ","), length))
                if(comma > semicolon){
                  svalue(csvseperator) <- ","
                  svalue(csvdecimal) <- "."
                }else{
                  svalue(csvseperator) <- ";"
                  if(comma2 > dot){
                    svalue(csvdecimal) <- ","
                  }
                  else{
                    svalue(csvdecimal) <- "."
                  }
                }
              },
							 error=function(e){
								 gmessage(paste("There was a problem while preparing your data: '",e,"'"), "Problem",
													icon="error")       
							 })
              previewCSV()
            })
    }
    
    #creates the actual preview inside the table, also the handler for all gui elements
    #beside the OK-button
    previewCSV <- function(...){
      f <- gframe("Preview:")
      g <- ggroup(use.scrollwindow = TRUE)
      testimport <- NULL
      error <- FALSE
      if(svalue(csvfilename)==''){
        testimport <- data.frame(column="preview loading ...")
      }
      else{
        svalue(statusbar) <- "compiling preview!"
        tryCatch({testimport <- read.table(svalue(csvfilename), nrows=10,
                                           fill=svalue(csvfill),
                                           header=svalue(csvheader),
                                           strip.white=svalue(csvstrip.white),
                                           stringsAsFactors=svalue(csvstringsAsFactors),
                                           blank.lines.skip=svalue(csvblank.lines.skip), 
                                           sep=svalue(csvseperator),
                                           dec=svalue(csvdecimal),
                                           quote=svalue(csvquotes),
                                           skip=svalue(csvskip),
                                           na.strings=strsplit(svalue(csvnastrings),",")[[1]])
                  putVm("colclasses",NA)
                  putVm("changedTypes",FALSE)}, 
                 error=function(e){svalue(statusbar) <- "read.table was not successful, please check your settings";
                                   error<-TRUE})
      }
      if(is.null(testimport)==FALSE){
        svalue(statusbar) <- "preview complete!"
      }
      else{
        testimport <- data.frame(column="preview loading ...")
      }
      add(g, gtable(testimport), expand=TRUE)
      add(f, g, expand=TRUE)
      layout[6:10, 1:7, expand=TRUE] <- f
      putVm("dframe",testimport)
      
    }
    
    #creates the layout for the parameter widgets in the CSV import window
		#and setups the handler
    statusbar <- gstatusbar("")
    csvfilebutton <- gbutton("...", handler=buttonHandler)
    csvheader <- gcheckbox("header", checked=TRUE, handler=previewCSV)
    csvfill <- gcheckbox("fill", checked=TRUE, handler=previewCSV)
    csvstrip.white <- gcheckbox("strip white", , handler=previewCSV)
    csvblank.lines.skip <- gcheckbox("blank line skip", handler=previewCSV)
    csvstringsAsFactors <- gcheckbox("strings As Factors", handler=previewCSV)
    csvseperator <- gedit(",", handler=previewCSV)
    addHandlerKeystroke(csvseperator, previewCSV)
    csvdecimal <- gedit(".", handler=previewCSV)
    addHandlerKeystroke(csvdecimal, previewCSV)
    csvquotes <- gedit("\"", handler=previewCSV)
    addHandlerKeystroke(csvquotes, previewCSV)
    csvskip <- gedit("0")
    addHandlerKeystroke(csvskip, previewCSV)
    csvnastrings <- gedit("")
    addHandlerKeystroke(csvnastrings, previewCSV)
    csvaccept <- gbutton("Accept", handler=function(...){
      ###real CSV import after pressing the accept button
      tryCatch({testimport <- getVm("dframe")
                if(getVm("changedTypes")==TRUE){
                  colclasses <- getVm("colclasses")
                  colclassesSTR <- parseVarStr(colclasses)
                }  
                else{
                  colclasses <- NA
                  colclassesSTR <- "NA"
                }
                wd <- WaitingDialog(Parent=importDialog)
                focus(wd) <- TRUE
                putVm("importFilename",svalue(csvfilename))
                filename=gsub("\\\\","/",svalue(csvfilename))
								#actual import of CSV
                df <- read.table(svalue(csvfilename),
                                 fill=svalue(csvfill),
                                 header=svalue(csvheader),
                                 strip.white=svalue(csvstrip.white),
                                 stringsAsFactors=svalue(csvstringsAsFactors),
                                 blank.lines.skip=svalue(csvblank.lines.skip), 
                                 sep=svalue(csvseperator),
                                 dec=svalue(csvdecimal),
                                 quote=svalue(csvquotes),
                                 skip=svalue(csvskip),
                                 colClasses=colclasses,
                                 na.strings=strsplit(svalue(csvnastrings),",")[[1]])
                dname <- format(Sys.time(), "importedCSV_%H_%M")
								#build a command line script for the script window
                cmdimp <- paste("activedataset <- read.table(\"",filename,"\"",
                                ",fill=",svalue(csvfill), 
                                ",header=",svalue(csvheader),
                                ",strip.white=",svalue(csvstrip.white),
                                ",stringsAsFactors=",svalue(csvstringsAsFactors),
                                ",blank.lines.skip=",svalue(csvblank.lines.skip),
                                ",sep=",parseVarStr(svalue(csvseperator)),
                                ",dec=",parseVarStr(svalue(csvdecimal)),
                                ",quote=\"\\",svalue(csvquotes),"\"",
                                ",skip=",parseVarStr(svalue(csvskip)),
                                ",colClasses=",colclassesSTR,
                                ",na.strings=",parseVarStr(svalue(strsplit(svalue(csvnastrings),",")[[1]])),
                                ")", sep="")
                setActiveDataset(df, loadScript=cmdimp)
                putVm("importFileName", svalue(csvfilename))
                #save import parameters for later export
                csvimportparams <- list(fill=svalue(csvfill),
                                        header=svalue(csvheader),
                                        strip.white=svalue(csvstrip.white),
                                        stringsAsFactors=svalue(csvstringsAsFactors),
                                        blank.lines.skip=svalue(csvblank.lines.skip), 
                                        sep=svalue(csvseperator),
                                        dec=svalue(csvdecimal),
                                        quote=svalue(csvquotes),
                                        skip=svalue(csvskip),
                                        colClasses=colclasses,
                                        na.strings=strsplit(svalue(csvnastrings),",")[[1]])
                putVm("csvimportparameters", csvimportparams)
                dispose(wd)
                dispose(importDialog)
      },
               error=function(e){gmessage(paste("There was a problem while importing your data: '",e,"'"), "Problem",
                                          icon="error")})
      
    })
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
    #more CSV window layout
    ftop <- gframe("Choose CSV-File:")
    gtop <- ggroup(horizontal=TRUE, container=ftop)
    add(ftop, csvfilename, expand=TRUE)
    add(ftop, csvfilebutton)
    layout[1,1:7] <- ftop
    fparams <- gframe("CSV-Parameters:")
    glayout <- glayout(container=fparams)
    glayout[2,1] <- csvheader
    glayout[3,1] <- csvfill
    glayout[4,1] <- csvstrip.white
    glayout[5,1] <- csvstringsAsFactors
    glayout[2,2] <- csvblank.lines.skip
    glayout[2,3, anchor=c(0,0)] <- glabel("seperator:")
    glayout[2,4] <- csvseperator
    glayout[3,3, anchor=c(0,0)] <- glabel("decimal:")
    glayout[3,4] <- csvdecimal
    glayout[4,3, anchor=c(0,0)] <- glabel("quotes:")
    glayout[4,4] <- csvquotes
    glayout[5,3, anchor=c(0,0)] <- glabel("skip:")
    glayout[5,4] <- csvskip
    glayout[2,5, anchor=c(0,0)] <- glabel("NA-strings:")
    glayout[2,6, expand=FALSE] <- csvnastrings
    layout[2:5, 1:7] <- fparams
		#init window after layouting
    previewCSV()
    layout[11,6, expand=FALSE] <- csvaccept
    layout[11,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
    add(importDialogFrame, statusbar)
    buttonHandler()
  }
  
  #opens a small window for importing SPSS files
	#called after clicking the corresponding menu entry
  importSPSS <- function(...){
    importDialog <- gwindow("Import SPSS", parent=window, width=100, height=100)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    filename <- gedit()
    enabled(filename) <- FALSE
    #handler for open file button
		#opens opens file window
    buttonHandler <- function(...){
      gfile(text = "Open SPSS File", type = "open", , 
            filter=list("SPSS files"=list(patterns=c("*.sav", "*.SAV")),"All files" = list(patterns = c("*"))),
            handler=function(h,...){
              svalue(filename) <- h$file
            })
    }
    
    #setup SPSS import gui
    statusbar <- gstatusbar("")
    filebutton <- gbutton("...", handler=buttonHandler)
    check.use.value.labels <- gcheckbox("convert value labels to factors")
    check.lowernames <- gcheckbox("convert variable names to lower case")
    check.force.single <- gcheckbox("force storage mode double to single", checked=TRUE)
    check.charfactor <- gcheckbox("convert character variables to factors")
    csvaccept <- gbutton("Accept", handler=function(...){
      #try to import spss file, if not message error
      tryCatch({
        wd <- WaitingDialog(Parent=importDialog)
        focus(wd) <- TRUE
				#actual import
        df <- spss.get(svalue(filename),
                       use.value.labels = svalue(check.use.value.labels),
                       lowernames = svalue(check.lowernames),
                       force.single = svalue(check.force.single),
                       charfactor= svalue(check.charfactor),
                       to.data.frame = TRUE)
        putVm("importFilename",svalue(filename))
        filename=gsub("\\\\","/",svalue(filename))
				#build command line script for script window
        cmdimp <- paste("activedataset <- spss.get(\"",gsub("\\\\","/",filename),"\"",
                        ",use.value.labels=",svalue(check.use.value.labels), 
                        ",lowernames=",svalue(check.lowernames),
                        ",force.single=",svalue(check.force.single),
                        ",charfactor=",svalue(check.charfactor),
                        ",to.data.frame = TRUE)", sep="")
				#set dataset close window
        setActiveDataset(df, loadScript=cmdimp)
        dispose(wd)
        dispose(importDialog)
      },error=function(e){
        gmessage(paste("There was a problem while importing your SPSS file: '",e,"'"),"Import Error!",icon="error")
      })
      
    })
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
		#more layout setup
    ftop <- gframe("Choose SPSS-File:")
    gtop <- ggroup(horizontal=TRUE, container=ftop)
    add(ftop, filename, expand=TRUE)
    add(ftop, filebutton)
    layout[1,1:7] <- ftop
    fparams <- gframe("SPSS-Parameters:")
    glayout <- glayout(container=fparams)
    glayout[1,1] <- check.use.value.labels
    glayout[1,2] <- check.lowernames
    glayout[2,1] <- check.force.single
    glayout[2,2] <- check.charfactor
    layout[2:3, 1:7] <- fparams
    layout[4,6, expand=FALSE] <- csvaccept
    layout[4,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
    buttonHandler()
  }
  
  #opens a small window for importing STATA files
	#called after clicking the corresponding menu entry
  importSTATA <- function(...){
		#create window layout and init variables
    importDialog <- gwindow("Import STATA", parent=window, width=100, height=100)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    filename <- gedit()
    enabled(filename) <- FALSE
    
		#handler for the choose file button
		#opens file selection window
    buttonHandler <- function(...){
      gfile(text = "Open STATA File", 
            filter=list("STATA files"=list(patterns=c("*.dta", "*.DTA")),"All files" = list(patterns = c("*"))),
            type = "open", handler=function(h,...){
              svalue(filename) <- h$file
            })
    }
    
    #setup STATA import gui
    statusbar <- gstatusbar("")
    filebutton <- gbutton("...", handler=buttonHandler)
    check.use.value.labels <- gcheckbox("convert value labels to factors", checked=TRUE)
    csvaccept <- gbutton("Accept", handler=function(...){
      #try to import stata file, if not message error
      tryCatch({
        wd <- WaitingDialog(Parent=importDialog)
        focus(wd) <- TRUE
        #actual loading
        df <- read.dta(svalue(filename),
                       convert.factors = svalue(check.use.value.labels))
        putVm("importFilename",svalue(filename))
        filename=gsub("\\\\","/",svalue(filename))
				#build command line script for script window
        cmdimp <- paste("activedataset <- read.dta(\"",gsub("\\\\","/",filename),"\"",
                        ",convert.factors=",svalue(check.use.value.labels),")",sep="")
				#set dataset active and close window
        setActiveDataset(df, loadScript=cmdimp)
        dispose(wd)
        dispose(importDialog)  
      },error=function(e){
        gmessage(paste("There was a problem while importing your STATA file: ",e,"'"),"Import Error!",icon="error")
      })
      
    })
		#more window layout 
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
    ftop <- gframe("Choose STATA-File:")
    gtop <- ggroup(horizontal=TRUE, container=ftop)
    add(ftop, filename, expand=TRUE)
    add(ftop, filebutton)
    layout[1,1:7] <- ftop
    fparams <- gframe("STATA-Parameters:")
    glayout <- glayout(container=fparams)
    glayout[1,1] <- check.use.value.labels
    layout[2, 1:7] <- fparams
    layout[3,6, expand=FALSE] <- csvaccept
    layout[3,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
    buttonHandler()
  }
  
  #opens a small window for importing SAS files
	#called after clicking the corresponding menu entry
  importSAS <- function(...){
		#init window and variables
    importDialog <- gwindow("Import SAS", parent=window, width=100, height=100)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    filename <- gedit()
    enabled(filename) <- FALSE
    
		#handler for file choose button
		#opens file choose dialog
    buttonHandler <- function(...){
      gfile(text = "Open SAS Export File", 
            filter=list("SAS XPORT"=list(patterns=c("*.xpt", "*.XPT")),"All files" = list(patterns = c("*"))),
            type = "open", handler=function(h,...){
              svalue(filename) <- h$file
            })
    }
    
    #setup SAS import gui
    statusbar <- gstatusbar("")
    filebutton <- gbutton("...", handler=buttonHandler)
    csvaccept <- gbutton("Accept", handler=function(...){
      #try to import sas file, if not message error
      tryCatch({
        wd <- WaitingDialog(Parent=importDialog)
        focus(wd) <- TRUE
				#actual import
        df <- sasxport.get(svalue(filename))
        putVm("importFilename",svalue(filename))
        filename=gsub("\\\\","/",svalue(filename))
				#build command line script for script window
        cmdimp <- paste("activedataset <- sasxport.get(\"",gsub("\\\\","/",filename),"\")",sep="")
        setActiveDataset(df, loadScript=cmdimp)
        
      },error=function(e){
        gmessage(paste("There was a problem while importing your SAS file: '",e,"'"),"Import Error!",icon="error")
      })
      
    })
		#more window layout
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
    ftop <- gframe("Choose SAS-File:")
    gtop <- ggroup(horizontal=TRUE, container=ftop)
    add(ftop, filename, expand=TRUE)
    add(ftop, filebutton)
    layout[1,1:7] <- ftop
    layout[2,6, expand=FALSE] <- csvaccept
    layout[2,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
    buttonHandler()
  }
  
  #opens a small window for exporting CSV files
	#called after clicking the corresponding menu entry
  exportCSV <- function(...){
		#init window and variables
    importDialog <- gwindow("Export CSV", parent=window, width=200, height=200)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    csvfilename <- gedit()
    enabled(csvfilename) <- FALSE
    
		#handler for file chooser button
		#opens file choose dialog
    buttonHandler <- function(...){
      gfile(text = "Save CSV File", type = "save", filter=list("CSV-Files"=list("*.csv")),handler=function(h,...){
        if(grepl("^.*\\.(csv|CSV)$", h$file)){
          svalue(csvfilename) <- h$file
        }
        else{
          svalue(csvfilename) <- paste(h$file, ".csv", sep="")
        }
      })
    }
    
    #retrieves previously saved CSV import settings
		#to allow simple export of similar CSV
		#or use default settings
    if(existsVm("csvimportparameters")){
      ip <- getVm("csvimportparameters")
      ip$na.strings<-"NA"
    }
    else{
      ip <- list(fill=TRUE,
                 header=TRUE,
                 strip.white=TRUE,
                 stringsAsFactors=TRUE,
                 blank.lines.skip=TRUE, 
                 sep=",",
                 dec=".",
                 quote="'",
                 skip=0,
                 colClasses=NULL,
                 na.strings="NA")
    }
    
		#more layouting
    statusbar <- gstatusbar("")
    csvfilebutton <- gbutton("...", handler=buttonHandler)
    csvheader <- gcheckbox("header", checked=ip$header)
    csvseperator <- gedit(ip$sep)
    csvdecimal <- gedit(ip$dec)
    csvnastrings <- gedit(ip$na.strings)
		#accept button handler
		#saves CSV file
    csvaccept <- gbutton("Accept", handler=function(...){
      dataobject <- getVm("activeDataSetOriginal")
			#if data was imputed, ask user which version to export
      if (is.null(getVm("activeDataSetImputed")) == FALSE){
        w <- gconfirm("Do you want to use the imputed values?", title="Imputed Values", icon="question")
        if (w == TRUE){
          dataobject <- getVm("activeDataSetImputed")
        }
      }
			#is data object is survey, only save real dataset
      if (is.survey(dataobject)){
        dataobject <- dataobject$variables
      }
			#remove delimiter variables
      dataobject <- dataobject[,grep("_imp", colnames(dataobject), invert=TRUE)]
			#actual export
      tryCatch({write.table(dataobject, file=svalue(csvfilename),
                            sep=svalue(csvseperator),
                            na=svalue(csvnastrings),
                            dec=svalue(csvdecimal),
                            row.names=svalue(csvheader))
                putVm("exportFileName", svalue(csvfilename))
								#create command line script for script browser
                cmdimp <- paste("write.table(activedataset, file='",gsub("\\\\","/",svalue(csvfilename)),"',",
                                "sep='",svalue(csvseperator),"',",
                                "na='",svalue(csvnastrings),"'',",
                                "dec='",svalue(csvdecimal),"',",
                                "row.names=",svalue(csvheader),")",sep="")
                addScriptLine(cmdimp)
                },
               error=function(e){gmessage(paste("There was a problem while exporting your data: '",e,"'"), "Problem",
                                          icon="error")})
      dispose(importDialog)
    })
		#more layouting
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
    ftop <- gframe("Choose CSV-File:")
    gtop <- ggroup(horizontal=TRUE, container=ftop)
    add(ftop, csvfilename, expand=TRUE)
    add(ftop, csvfilebutton)
    layout[1,1:7] <- ftop
    fparams <- gframe("CSV-Parameters:")
    glayout <- glayout(container=fparams)
    glayout[2,1] <- csvheader
    glayout[2,3, anchor=c(0,0)] <- glabel("seperator:")
    glayout[2,4] <- csvseperator
    glayout[3,3, anchor=c(0,0)] <- glabel("decimal:")
    glayout[3,4] <- csvdecimal
    glayout[2,5, anchor=c(0,0)] <- glabel("NA-strings:")
    glayout[2,6, expand=FALSE] <- csvnastrings
    layout[2:5, 1:7] <- fparams
    layout[9,6, expand=FALSE] <- csvaccept
    layout[9,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
		#init window
    buttonHandler()
  }
  
  #opens a small window for exporting SPSS files
	#called after clicking the corresponding menu entry
  exportSPSS <- function(...){
    #init window and different variables
    importDialog <- gwindow("Export SPSS", parent=window, width=100, height=100)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    datafilename <- gedit()
    codefilename <- gedit()
    enabled(datafilename) <- FALSE
    enabled(codefilename) <- FALSE
    
		#handler for choose file button
		#opens file selection dialog
    databuttonHandler <- function(...){
      gfile(text = "Save Data File", type = "save", handler=function(h,...){
        if(grepl("^.*\\.(dat)$", h$file)){
          svalue(datafilename) <- h$file
        }
        else{
          svalue(datafilename) <- paste(h$file, ".dat", sep="")
        }
      })
    }
    
		#handler for choose codebook file button
		#opens file selection dialog
    codebuttonHandler <- function(...){
      gfile(text = "Save SPS File", type = "save", filter=list(".sps"=list("*.sps")),handler=function(h,...){
        if(grepl("^.*\\.(sps|SPS)$", h$file)){
          svalue(codefilename) <- h$file
        }
        else{
          svalue(codefilename) <- paste(h$file, ".sps", sep="")
        }
      })
    }
    
    #setup sas export gui
    datafilebutton <- gbutton("...", handler=databuttonHandler)
    codefilebutton <- gbutton("...", handler=codebuttonHandler)
    csvaccept <- gbutton("Accept", handler=function(...){
      dataobject <- getVm("activeDataSetOriginal")
			#if imputed data exists ask user if he wants to export it
      if (is.null(getVm("activeDataSetImputed")) == FALSE){
        w <- gconfirm("Do you want to use the imputed values?", title="Imputed Values", icon="question")
        if (w == TRUE){
          dataobject <- getVm("activeDataSetImputed")
        }
      }
			#for surveys only export the actual data
      if (is.survey(dataobject)){
        dataobject <- dataobject$variables
      }
      dataobject <- dataobject[,grep("_imp", colnames(dataobject), invert=TRUE)]
			#acutal export
      tryCatch({write.foreign(dataobject, datafile=svalue(datafilename),
                              codefile=svalue(codefilename),
                              package = "SPSS")
                putVm("exportFileName", svalue(datafilename))
                #create command line for script window
                cmdimp <- paste("write.foreign(activedataset, datafile='",gsub("\\\\","/",svalue(datafilename)),"',",
                                              "codefile='",gsub("\\\\","/",svalue(codefilename)),"',",
                                              "package = 'SPSS')",sep="")
                addScriptLine(cmdimp)
                },
               error=function(e){gmessage(paste("There was a problem while exporting your data: '",e,"'"), "Problem",
                                          icon="error")})
      dispose(importDialog)
    })
		#more layout and init
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
    fdata <- gframe("Choose Data-File (Contains exported data as freetext):")
    gdata <- ggroup(horizontal=TRUE, container=fdata)
    add(fdata, datafilename, expand=TRUE)
    add(fdata, datafilebutton)
    layout[1,1:7] <- fdata
    fcode <- gframe("Choose Code-File (Contains SPSS Code for import):")
    gcode <- ggroup(horizontal=TRUE, container=fcode)
    add(fcode, codefilename, expand=TRUE)
    add(fcode, codefilebutton)
    layout[2,1:7] <- fcode
    layout[4,6, expand=FALSE] <- csvaccept
    layout[4,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
  }
  
  #creates a small window for exporting STATA files
	#called after clicking the corresponding menu entry
  exportSTATA <- function(...){
		#init window and variables
    importDialog <- gwindow("Export STATA", parent=window, width=100, height=100)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    filename <- gedit()
    enabled(filename) <- FALSE
    edit.version <- gedit("7")
    check.dates <- gcheckbox("convert dates to STATA-dates", checked=TRUE)
    combo.convert.factors <- gcombobox(c("labels","string","numeric","codes"))
    
		#handler for the file choose button
		#open file selection dialog
    buttonHandler <- function(...){
      gfile(text = "Save STATA File", type = "save", filter=list(".dta"=list("*.dta")),handler=function(h,...){
        if(grepl("^.*\\.(dta|DTA)$", h$file)){
          svalue(filename) <- h$file
        }
        else{
          svalue(filename) <- paste(h$file, ".dta", sep="")
        }
      })
    }
    
    #setup stata export gui
    filebutton <- gbutton("...", handler=buttonHandler)
		#handler for accept button
    csvaccept <- gbutton("Accept", handler=function(...){
      dataobject <- getVm("activeDataSetOriginal")
			#if imputed data exists, ask user how to proceed
      if (is.null(getVm("activeDataSetImputed")) == FALSE){
        w <- gconfirm("Do you want to use the imputed values?", title="Imputed Values", icon="question")
        if (w == TRUE){
          dataobject <- getVm("activeDataSetImputed")
        }
      }
			#for surveys only export actual data
      if (is.survey(dataobject)){
        dataobject <- dataobject$variables
      }
      dataobject <- dataobject[,grep("_imp", colnames(dataobject), invert=TRUE)]
			#actual export
      tryCatch({write.dta(dataobject, file=svalue(filename),
                          version=as.numeric(svalue(edit.version)),
                          convert.dates=svalue(check.dates),
                          convert.factors=svalue(combo.convert.factors))
                putVm("exportFileName", svalue(filename))
								#create command line script for script browser
                cmdimp <- paste("write.dta(activedataset, file='",gsub("\\\\","/",svalue(filename)),"',",
                                "version=",as.numeric(svalue(edit.version)),",",
                                "convert.dates=",svalue(check.dates),",",
                                "convert.factors = '",svalue(combo.convert.factors),"')",sep="")
                addScriptLine(cmdimp)
                },
               error=function(e){gmessage(paste("There was a problem while exporting your data: '",e,"'"), "Problem",
                                          icon="error")})
      dispose(importDialog)
    })
		#more layout and init
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
    ftop <- gframe("Choose STATA-File:")
    gtop <- ggroup(horizontal=TRUE, container=ftop)
    add(ftop, filename, expand=TRUE)
    add(ftop, filebutton)
    layout[1,1:7] <- ftop
    fparams <- gframe("STATA-Parameters:")
    glayout <- glayout(container=fparams)
    glayout[1,1] <- check.dates
    glayout[1,2, anchor=c(0,0)] <- glabel("version:")
    glayout[1,6, expand=FALSE] <- edit.version
    glayout[2,2, anchor=c(0,0)] <- glabel("handle factors as:")
    glayout[2,6, expand=FALSE] <- combo.convert.factors
    layout[2:3, 1:7] <- fparams
    layout[5,6, expand=FALSE] <- csvaccept
    layout[5,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
  }
  
  #opens small window for exporting SAS files
	#called after clicking the corresponding menu item 	
  exportSAS <- function(...){
    #init window and variables
    importDialog <- gwindow("Export SAS", parent=window, width=100, height=100)
    putVm("importDialog",importDialog)
    putVm("dframe", NULL)
    importDialogFrame <- ggroup(container=importDialog, horizontal=FALSE)
    layout <- glayout()
    datafilename <- gedit()
    codefilename <- gedit()
    enabled(datafilename) <- FALSE
    enabled(codefilename) <- FALSE
    edit.dataname <- gedit("rdata")
    combo.validvarname <- gcombobox(c("<=6",">=7"), selected=2)
    
		#handler for the data file chooser
		#opens file selection dialog
    databuttonHandler <- function(...){
      gfile(text = "Save Data File", type = "save", handler=function(h,...){
        if(grepl("^.*\\.(dat)$", h$file)){
          svalue(datafilename) <- h$file
        }
        else{
          svalue(datafilename) <- paste(h$file, ".dat", sep="")
        }
      })
    }
		
    #handler for the code file chooser
		#opens file selection dialog
    codebuttonHandler <- function(...){
      gfile(text = "Save SAS File", type = "save", filter=list(".sas"=list("*.sas")),handler=function(h,...){
        if(grepl("^.*\\.(sas|SAS)$", h$file)){
          svalue(codefilename) <- h$file
        }
        else{
          svalue(codefilename) <- paste(h$file, ".sas", sep="")
        }
      })
    }
    
    #setup sas export gui
    datafilebutton <- gbutton("...", handler=databuttonHandler)
    codefilebutton <- gbutton("...", handler=codebuttonHandler)
		#accept button handler
    csvaccept <- gbutton("Accept", handler=function(...){
      dataobject <- getVm("activeDataSetOriginal")
			#in case of imputed data export ask user which to export
      if (is.null(getVm("activeDataSetImputed")) == FALSE){
        w <- gconfirm("Do you want to use the imputed values?", title="Imputed Values", icon="question")
        if (w == TRUE){
          dataobject <- getVm("activeDataSetImputed")
        }
      }
			#for surveys only save the actual data
      if (is.survey(dataobject)){
        dataobject <- dataobject$variables
      }
      dataobject <- dataobject[,grep("_imp", colnames(dataobject), invert=TRUE)]
			#actual export
      tryCatch({version <- paste("V",substr(svalue(combo.validvarname), 3,3), sep="")
                write.foreign(dataobject, datafile=svalue(datafilename),
                              codefile=svalue(codefilename),
                              package = "SAS",
                              dataname = svalue(edit.dataname),
                              validvarname = version)
                putVm("exportFileName", svalue(datafilename))
								#create command line script for script browser
                cmdimp <- paste("write.foreign(activedataset, datafile='",gsub("\\\\","/",svalue(datafilename)),"',",
                                "codefile='",gsub("\\\\","/",svalue(codefilename)),"',",
                                "package = 'SAS', validvarname='",version,"')",sep="")
                addScriptLine(cmdimp)
                },
               error=function(e){gmessage(paste("There was a problem while exporting your data: '",e,"'"), "Problem",
                                          icon="error")})
      dispose(importDialog)
    })
		#more layout and init
    csvdiscard <- gbutton("Discard ", handler=function(...){dispose(importDialog)})
    fdata <- gframe("Choose Data-File (Contains exported data as free-text):")
    gdata <- ggroup(horizontal=TRUE, container=fdata)
    add(fdata, datafilename, expand=TRUE)
    add(fdata, datafilebutton)
    layout[1,1:7] <- fdata
    
    fcode <- gframe("Choose Code-File (Contains SAS Code for import):")
    gcode <- ggroup(horizontal=TRUE, container=fcode)
    add(fcode, codefilename, expand=TRUE)
    add(fcode, codefilebutton)
    layout[2,1:7] <- fcode
    
    fparams <- gframe("SAS-Parameters:")
    glayout <- glayout(container=fparams)
    glayout[1,1, anchor=c(-1,0)] <- glabel("future SAS data set name:")
    glayout[1,2, expand=FALSE] <- edit.dataname
    glayout[2,1, anchor=c(-1,0)] <- glabel("SAS version :")
    glayout[2,2, expand=FALSE] <- combo.validvarname
    layout[3:4, 1:7] <- fparams
    layout[6,6, expand=FALSE] <- csvaccept
    layout[6,7, expand=FALSE] <- csvdiscard
    add(importDialogFrame, layout, expand=TRUE)
  }
  
	#handler for the main notebook
	#does some init and update (e.i. plots)
  #called after switching the main tabs 
  updatePanels <- function(pageno=svalue(mainNotebook), firstTime=FALSE){
		#data tab
    if(pageno==1){
			#allow user different functionality for imputed data
      if (is.null(getVm("activeDataSetImputed"))) {
        svalue(dataPanel.imputationSelection) <- "original"
        enabled(dataPanel.imputationSelection) <- FALSE
      }
      else{
        svalue(dataPanel.imputationSelection) <- "imputed"
        enabled(dataPanel.imputationSelection) <- TRUE
      }
      #adapt widgets above table to dataset
			#block change handler to improve performance
      blockHandler(dataPanel.variableSelection, dataPanel.variableSelectionHandler)
      blockHandler(dataPanel.bySelection, dataPanel.bySelectionHandler)
      blockHandler(dataPanel.statisticsSelection, dataPanel.statisticsSelectionHandler)
      dataPanel.variableSelection[]<- getVariableNames(getVm("activeDataSetOriginal"))
      dataPanel.bySelection[]<- c(" ", getVariableNames(getVm("activeDataSetOriginal")))
      svalue(dataPanel.bySelection) <- " "
      #change available statistics depending on data is survey or not
      if (is.survey(getVm("activeDataSetOriginal"))){
        dataPanel.statisticsSelection[] <- c("svymean", "svyvar", "svytotal")
        svalue(dataPanel.statisticsSelection) <- "svymean"
      }
      else{
        dataPanel.statisticsSelection[] <- c("mean", "var")
        svalue(dataPanel.statisticsSelection) <- "mean"
      }
			#if is first pass (after program start) make additional inits
      if (firstTime == TRUE){
        svalue(dataPanel.variableSelection, index=TRUE) <- 1
        svalue(dataPanel.bySelection) <- " "
        dataPanelChangeHandler()
      }
			#reinstate different handlers
      unblockHandler(dataPanel.statisticsSelection, dataPanel.statisticsSelectionHandler)
      unblockHandler(dataPanel.variableSelection, dataPanel.variableSelectionHandler)
      unblockHandler(dataPanel.bySelection, dataPanel.bySelectionHandler)
    }
		#imputation tab
    if(pageno==2){
      #for possible future additions
    }
		#visualization tab
    if(pageno==3){
			#different settings if imputed data is available
      if (is.null(getVm("activeDataSetImputed"))) {
        svalue(impVis.plotImputed) <- "original"
        enabled(impVis.plotImputed) <- FALSE
      }
      else{
        svalue(impVis.plotImputed) <- "imputed"
        enabled(impVis.plotImputed) <- TRUE
      }
			#make sure that no empty plot window exists
      visible(impVis.plot) <- TRUE
      makeImputationPlot()
    }
  }
  
  #adaptes all widgets on all panels to a new dataset
	#called while loading a new dataset (setActiveDataset())
	#fills a lot of tables with the possible variable names of the dataset
  initPanels <- function(){
    clearTable(dataPanel.table)
    
    #init the imputation tab
    variables <- getVariableNames(getVm("activeDataSetOriginal"))
    variables <- data.frame(var=variables)
    insertTable(imputation.hotdeck.variable, variables)
    insertTable(imputation.hotdeck.ord_var, variables)
    insertTable(imputation.hotdeck.domain_var, variables)
    insertTable(imputation.irmi.variable, variables)
    insertTable(imputation.irmi.count, variables)
    insertTable(imputation.kNN.variable, variables)
    insertTable(imputation.kNN.dist_var, variables)
    insertTable(imputation.regression.variables, variables)
    clearTable(imputation.kNN.mixed)
    clearTable(imputation.irmi.mixed)
    
		#init the plot tab
		#block handlers to improve performance
    sapply(handlerList, FUN=function(s){blockHandler(s[[1]],s[[2]])})
    try(impVis.barMiss.pos[] <- variables, silent = TRUE)
    try(impVis.histMiss.pos[] <- variables, silent = TRUE)
    try(impVis.pbox.pos[] <- variables, silent = TRUE)
    try(impVis.matrixplot.sortby[] <- variables, silent = TRUE)
    svalue(impVis.barMiss.pos, index=TRUE) <- 1
    svalue(impVis.histMiss.pos, index=TRUE) <- 1
    svalue(impVis.pbox.pos, index=TRUE) <- 1
    svalue(impVis.matrixplot.sortby, index=TRUE) <- 1
    insertTable(impVis.scattmatrixMiss.highlight, variables)
    insertTable(impVis.scattmatrixMiss.plotvars, variables)
    insertTable(impVis.marginMatrix.plotvars, variables)
		#default shown variables for scatterplot matrices
		#the first 5 existing variables
    svalue(impVis.marginMatrix.plotvars, index=TRUE) <- 1:min(dim(variables)[1],5)
    svalue(impVis.scattmatrixMiss.plotvars, index=TRUE) <- 1:min(dim(variables)[1],5)
    insertTable(impVis.mosaicMiss.highlight, variables)
    insertTable(impVis.mosaicMiss.plotvars, variables)
    insertTable(impVis.parcoordMiss.highlight, variables)
    insertTable(impVis.parcoordMiss.plotvars, variables)
    svalue(impVis.parcoordMiss.plotvars, index=TRUE) <- 1:min(dim(variables)[1],5)
		#choose default shown variables for mosaic plot
		#heuristic: the first 2 discrete (less then 25 different values) variables
		#otherwise just the first two variables
    dat <- getVm("activeDataSetOriginal")
    if(is.survey(dat)){
      dat <- dat$variables
    }
    vc <- sapply(dat, is.categorical)
    vc <- which(vc, TRUE)
    if (length(vc) < 2){
      vc <- 1:2
    }
    else{
      vc <- vc[1:2]
    }
    svalue(impVis.mosaicMiss.plotvars, index=TRUE) <- vc
		#reinstate handlers
    sapply(handlerList, FUN=function(s){unblockHandler(s[[1]],s[[2]])})
    
    #default focus for regression tabs
    putVm("regressionFocus", 0)
    
    #if surveyobject allow weighted plots
    if (is.survey(getVm("activeDataSetOriginal")) == TRUE){
      enabled(impVis.aggr.weighted) <- TRUE
      enabled(impVis.histMiss.weighted) <- TRUE
      enabled(impVis.barMiss.weighted) <- TRUE
      enabled(impVis.scattmatrixMiss.weighted) <- TRUE
      enabled(impVis.mosaicMiss.weighted) <- TRUE
      enabled(impVis.parcoordMiss.weighted) <- TRUE
      enabled(impVis.pbox.weighted) <- TRUE
      enabled(impVis.matrixplot.weighted) <- TRUE
    }
    else{
      enabled(impVis.aggr.weighted) <- FALSE
      enabled(impVis.histMiss.weighted) <- FALSE
      enabled(impVis.barMiss.weighted) <- FALSE
      enabled(impVis.scattmatrixMiss.weighted ) <- FALSE
      enabled(impVis.mosaicMiss.weighted) <- FALSE
      enabled(impVis.parcoordMiss.weighted) <- FALSE
      enabled(impVis.pbox.weighted) <- FALSE
      enabled(impVis.matrixplot.weighted) <- FALSE
    }
    clearTable(imputation.undo.variables)
    
    #adaped summary page to dataset
    svalue(dataPanel.summarytext) <- ""
    sumdat <- data.frame()
		#if data is survey
    if(is.survey(getVm("activeDataSetOriginal"))){
      d <- dim(getVm("activeDataSetOriginal"))
			#put survey summary + head count into summary field
      svalue(dataPanel.summarytext) <- paste(paste(capture.output(print(getVm("activeDataSetOriginal"))), collapse="\n"),"\n",d[1],' Observations  -  ',d[2],' Variables',collapse="")
      sumdat <- getVm("activeDataSetOriginal")$variables
    }
    else{
      d <- dim(getVm("activeDataSetOriginal"))
			#test if a non-empty data set is loaded
      if (d[1] == 1 & d[2] == 1){
        svalue(dataPanel.summarytext) <- "No Dataset loaded!"
      }
      else{
				#write head-count into summary field
        svalue(dataPanel.summarytext) <- paste("Dataframe",paste(d[1],' Observations  -  ',d[2],' Variables',sep=""),sep="\n")
      }
      sumdat <- getVm("activeDataSetOriginal")
    }
    sumtable <- createSummaryDataframe(sumdat)
    insertTable(dataPanel.summaryTable, sumtable)
    insertTable(imputation.summarytable, sumtable)

    #init the data structure for the undo imputation tab
    vars <- rep(FALSE, length(getVariableNames(getVm("activeDataSetOriginal"))))
    methods <- rep(" ", length(getVariableNames(getVm("activeDataSetOriginal"))))
    times <- rep(" ", length(getVariableNames(getVm("activeDataSetOriginal"))))
    vars <- data.frame(varnames=vars, methods=methods, times=times, stringsAsFactors = FALSE)
    rownames(vars) <- getVariableNames(getVm("activeDataSetOriginal"))
    putVm("ImputedVariables", vars)
  }
  
  #creates the option window which allows the user to select the colors for the plots
	#called after clicking the corresponding menu entry
  OptionsHandler <- function(h,...){
		#init window, varaibles and layout
    options.window <- gwindow(title="Options", width=50, height=50)
		#the init takes rather long, so put a "please wait"-message there
    temp <- ggroup(container=options.window)
    glabel("<b><big>LOADING!</big></b>", container=temp, markup = TRUE)
    options.layout <- glayout()
    c <- colors()
    colors <- getVm("plotColors")
    options.color1 <- gdroplist(c)
    options.color2 <- gdroplist(c)
    options.color3 <- gdroplist(c)
    options.color4 <- gdroplist(c)
    options.color5 <- gdroplist(c)
    options.color6 <- gdroplist(c)
    svalue(options.color1) <- colors[1]
    svalue(options.color2) <- colors[2]
    svalue(options.color3) <- colors[3]
    svalue(options.color4) <- colors[4]
    svalue(options.color5) <- colors[5]
    svalue(options.color6) <- colors[6]
    options.alpha <- gslider(from=0, to=1, by=0.01)
    svalue(options.alpha) <- getVm("plotAlpha")
    options.layout[1,1] <- glabel("Color 1:")
    options.layout[1,2] <- options.color1
    options.layout[2,1] <- glabel("Color 2:")
    options.layout[2,2] <- options.color2
    options.layout[3,1] <- glabel("Color 3:")
    options.layout[3,2] <- options.color3
    options.layout[4,1] <- glabel("Color 4:")
    options.layout[4,2] <- options.color4
    options.layout[5,1] <- glabel("Color 5:")
    options.layout[5,2] <- options.color5
    options.layout[6,1] <- glabel("Color 6:")
    options.layout[6,2] <- options.color6
    options.layout[7,1:2, anchor <- c(-1, -1)] <- glabel("Alpha:")
    options.layout[8,1:2] <- options.alpha
    options.smallgroup <- ggroup()
		#handler for discard button, closes window
    gbutton(" Discard", container=options.smallgroup, handler=function(h,...){
      dispose(options.window)
    })
		#handler for accept button
		#saves colors and redraws plot
    gbutton(" Accept", container=options.smallgroup, handler=function(h,...){
      color <- c(svalue(options.color1),
                 svalue(options.color2),
                 svalue(options.color3),
                 svalue(options.color4),
                 svalue(options.color5),
                 svalue(options.color6))
      putVm("plotColors", color)
      putVm("plotAlpha", svalue(options.alpha))
      makeImputationPlot()
      dispose(options.window)
    })
    options.layout[9,2] <- options.smallgroup
		#loading is done now, remove "please wait" and substitute real widgets
    delete(options.window, temp)
    add(options.window, options.layout)
  }
  
  #shows a copyable list of the R code of the called operations
	#called after somebody clicked the script menu item
  ScriptHandler <- function(h,...){
		#init window and layout
    script.window <- gwindow("Script View", width=1024, height=300)
    script.code <- gtext(container=script.window)
    hist <- getVm("ScriptHistory")
		#piece together all saved script lines and present them in textbox
    script_text <- character(0)
    for (line in hist){
      if (!isEmpty(line)){
        script_text <- paste(script_text,as.character(line), sep="\n")
      }
    }
    svalue(script.code) <- script_text
  }
  
  #small helper which simplifies the addition of a new line of code
	#puts "line" at the end of a list of all added script lines 
  addScriptLine <- function(line){
    putVm("ScriptHistory", c(getVm("ScriptHistory"), line))
  }
  
  #draws a plot, created by expression plotExpr, onto a gimage widget
  #this plot is temporary created as image file inside the working directory
  #this resolves some bugs with the as.cairodevice-methode and allows for better image quality
	#plotExpr 	...		ordinary R-code for to-be-drawn plot
	#target			...		widget to put the created image
	#savePlot		...		small workaround to double use this function for also permanent image saving
  bufferedPlot <- function(plotExpr, target=impVis.plot, savePlot=FALSE){
    #if the plot graphics has to be saved to disk with user defined criteria
    #outsource this to other function
    #done to reduce overhead while rewriting plot functionality
    if (savePlot){
      eval(plotExpr)
    }
    else{
      try({
				#get size of widget to preserve plot proportions
        tgtk <- getToolkitWidget(target)
        a <- tgtk$getAllocation()
				#save temporary as PNG in working directory
        dev <- CairoPNG(filename="current.tmp", a$allocation$width-3, a$allocation$height-3)
        eval(plotExpr)
        dev.off(dev)
        svalue(target) <- "current.tmp"
        file.remove("current.tmp")
      }, silent=TRUE)
    }

  }
  
  #sets the variable to its original value and removes the entry from the table
	#called after somebody double-clicks a row to undo the corresponding imputation
  undoImputation <- function(h,...){
		#get the to-undo variable and ask user if he is sure
    varname <- as.character(svalue(imputation.undo.variables))
    w <- gconfirm(paste("Undo imputation of",varname,"?"))
    if (w==TRUE){
      orig <- getVm("activeDataSetOriginal")
      imp <- getVm("activeDataSetImputed")
      #change value back to original, mind if survey or not
      if (is.survey(orig)){
        imp$variables[,varname] <- orig$variables[,varname]
        imp$variables[,paste(varname,"_imp",sep="")] <- NULL
        #build command line script for script browser
        addScriptLine(paste("activedataset$variables$",varname," <- originaldataset$variables$",varname, sep=""))
        addScriptLine(paste("activedataset$variables$",varname,"_imp <- NULL", sep=""))
      }
      else{
        imp[,varname] <- orig[,varname]
        imp[,paste(varname,"_imp",sep="")] <- NULL
        #build command line script for script browser
        addScriptLine(paste("activedataset$",varname," <- originaldataset$",varname, sep=""))
        addScriptLine(paste("activedataset$",varname,"_imp <- NULL", sep=""))
      }
      putVm("activeDataSetImputed",imp)
      #remove entry from table
      vars <- getVm("ImputedVariables")
      vars[varname,1] <- FALSE
      vars[varname,2] <- " "
      vars[varname,3] <- " "
      putVm("ImputedVariables", vars)
      sumtable <- createSummaryDataframe(imp)
			#update tables
      insertTable(imputation.summarytable, sumtable)
      updateUndoVariablesTable()
    }
  }
  
	#main working method for the imputation
	#handler for the "apply imputation" button
  #applies selected imputation method, updates dataset and tables
	#called after somebody presses the "Apply Imputation"-Button
  Imputation <- function(h,...){
		#if there was already a previous imputation, apply imputation to that dataset
    if (is.null(getVm("activeDataSetImputed"))){
      dataset <- getVm("activeDataSetOriginal")
    }
    else{
      dataset <- getVm("activeDataSetImputed")
    }
    
		#select imputation method depending on which tab the user is on
    #perform kNN imputation
    if (svalue(imputation.notebook)==1){
      #convert widgets to useable function parameters, i.e. convert from character to numeric,
			#substitute NULL for empty elements,...
      variable <- getVariableNames(getVm("activeDataSetOriginal"))
      v <- as.character(svalue(imputation.kNN.variable))
      if (length(v) > 0){
        variable <- v
      }
      
      k <- 5
      v <- svalue(imputation.kNN.k)
      if (v != ""){
        k <- v
      }
      
      dist_var <- getVariableNames(getVm("activeDataSetOriginal"))
      v <- as.character(svalue(imputation.kNN.dist_var))
      if (length(v) > 0){
        dist_var <- v
      }
            
      weights <- NULL
      v <- cutParam(svalue(imputation.kNN.weights))
      if (length(v) > 0) {
        weights <- v
      }
      
      mixed <- vector()
      v <- as.character(svalue(imputation.kNN.mixed))
      if (length(v) > 0){
        mixed <- v
      }
      
      mixed.constant <- NULL
      v <- cutParam(svalue(imputation.kNN.mixed.constant))
      if (length(v) > 0){
        mixed.constant <- v
      }
      
			#perform actual imputation
			#capture different errors and warnings
      sumText <- capture.output(impData <- tryCatch({kNN(dataset,
                                                        variable = variable,
                                                        k = k,
                                                        dist_var = dist_var,
                                                        weights = weights,
                                                        numFun = get(svalue(imputation.kNN.numFun)),
                                                        catFun = get(svalue(imputation.kNN.catFun)),
                                                        #makeNA = makeNA,
                                                        impNA = svalue(imputation.kNN.impNA),
                                                        addRandom = svalue(imputation.kNN.addRandom),
                                                        mixed = mixed,
                                                        mixed.constant = mixed.constant)},
                                                    error=function(e){
                                                      message(e$message)
                                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Warning", icon="error")
                                                      return(NULL)
                                                    },
                                                    warning=function(e){
                                                      message(e$message)
                                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Warning", icon="error")
                                                      return(NULL)
                                                    }))
      
      #no error while imputation
      if(!is.null(impData)) {
				#save imputed data and update different tables 
        putVm("activeDataSetImputed", impData)  
				enabled(impVis.plotImputed) <- TRUE 
				#test if there was a actual imputations, i.e. if there are fewer missing values
				vars <- getVm("ImputedVariables")
				if (is.survey(impData)){
					variable <- compareImputations(dataset$variables, impData$variables)
				}
				else{
					variable <- compareImputations(dataset, impData)
				}
				#add the newly imputed variables to the undo tab
				vars[variable,1] <- TRUE
				vars[variable,2] <- "kNN"
				vars[variable,3] <- format(Sys.time(), "%H:%M:%S")
				putVm("ImputedVariables", vars)
				updateUndoVariablesTable()
				sumtable <- createSummaryDataframe(impData)
				insertTable(imputation.summarytable, sumtable)
				#create command-line script for script browser
				#convert to different parameter from widgets to readable strings
				variable <- getVariableNames(getVm("activeDataSetOriginal"))
				v <- as.character(svalue(imputation.kNN.variable))
				if (length(v) > 0){
					variable <- v
				}
				variable <- paste("c(",paste(sapply(variable, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				
				k <- "5"
				v <- svalue(imputation.kNN.k)
				if (v != ""){
					k <- as.character(v)
				}
				
				dist_var <- getVariableNames(getVm("activeDataSetOriginal"))
				v <- as.character(svalue(imputation.kNN.dist_var))
				if (length(v) > 0){
					dist_var <- v
				}
				dist_var <- paste("c(",paste(sapply(dist_var, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				
				weights <- "NULL"
				v <- cutParam(svalue(imputation.kNN.weights))
				if (length(v) > 0) {
					weights <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				
				mixed <- "NULL"
				v <- as.character(svalue(imputation.kNN.mixed))
				if (length(v) > 0){
					mixed <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				
				mixed.constant <- "NULL"
				v <- cutParam(svalue(imputation.kNN.mixed.constant))
				if (length(v) > 0){
					mixed.constant <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				cmd <- paste('activedataset <- kNN(activedataset, variable=',variable,',',
										 'k=',as.character(k),", ",
										 'dist_var=',dist_var,", ",
										 'weights=',weights,", ",
										 'numFun=',svalue(imputation.kNN.numFun),", ",
										 'catFun=',svalue(imputation.kNN.catFun),", ",
										 'impNA=',svalue(imputation.kNN.impNA),", ",
										 'addRandom=',svalue(imputation.kNN.addRandom),", ",
										 'mixed=',mixed,", ",
										 'mixed.constant=',mixed.constant,") ", collapse="")
				addScriptLine(cmd)
				gmessage("Imputation successful!", title="Success", icon="info")
        
      }
    }
    #perform irmi imputation
    if (svalue(imputation.notebook)==2){
      #convert widgets to useable function parameters, i.e. convert from character to numeric,
			#substitute NULL for empty elements,...
      variable <- getVariableNames(getVm("activeDataSetOriginal"))
      v <- as.character(svalue(imputation.irmi.variable))
      if (length(v) > 0){
        variable <- v
      }
      
      mixed <- NULL
      v <- as.character(svalue(imputation.irmi.mixed))
      if (length(v) > 0){
        mixed <- v
      }
      
      mixed.constant <- NULL
      v <- cutParam(svalue(imputation.irmi.mixed.constant))
      if (length(v) > 0){
        mixed.constant <- v
      }
      
      count <- NULL
      v <- as.character(svalue(imputation.irmi.count))
      if (length(v) > 0){
        count <- v
      }
      
			#perform actual imputation
			#capture different errors and warnings
      sumText <- capture.output(impData <- tryCatch({irmi(dataset,
                                                #variable = variable,
                                                mixed = mixed,
                                                mixed.constant = mixed.constant,
                                                count = count,
                                                robust = svalue(imputation.irmi.robust),
                                                noise = svalue(imputation.irmi.noise),
                                                noise.factor = svalue(imputation.irmi.noise.factor))},
                                               error=function(e){
                                                 message(e$message)
                                                 gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                                 return(NULL)
                                               },
                                               warning=function(e){
                                                 message(e$message)
                                                 gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                                 return(NULL)
                                               }))
      #no error while imputation
      if(!is.null(impData) ) {
				#save imputed data
        putVm("activeDataSetImputed", impData)  
				enabled(impVis.plotImputed) <- TRUE 
				vars <- getVm("ImputedVariables")
				#FOR FUTURE REFERENCE: IRMI VARIABLE SELECTION
				#vars[getVariableNames(getVm("activeDataSetOriginal"))] <- TRUE
				#variable <- getVariableNames(getVm("activeDataSetOriginal"))
				
				#find actual imputed variables and save them for undo operation
				if (is.survey(impData)){
					variable <- compareImputations(dataset$variables, impData$variables)
				}
				else{
					variable <- compareImputations(dataset, impData)
				}
				vars[variable,1] <- TRUE
				vars[variable,2] <- "irmi"
				vars[variable,3] <- format(Sys.time(), "%H:%M:%S")
				putVm("ImputedVariables", vars)
				updateUndoVariablesTable()
				sumtable <- createSummaryDataframe(impData)
				insertTable(imputation.summarytable, sumtable)
				#build command line string
				#convert the content of widgets to a actual string representation
				mixed <- "NULL"
				v <- as.character(svalue(imputation.irmi.mixed))
				if (length(v) > 0){
					mixed <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				
				mixed.constant <- "NULL"
				v <- cutParam(svalue(imputation.irmi.mixed.constant))
				if (length(v) > 0){
					mixed.constant <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				
				count <- "NULL"
				v <- as.character(svalue(imputation.irmi.count))
				if (length(v) > 0){
					count <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				cmd <- paste('activedataset <- irmi(activedataset, mixed=',mixed,',',
										 'mixed.constant=',mixed.constant,", ",
										 'count=',count,", ",
										 'robust=',svalue(imputation.irmi.robust),", ",
										 'noise=',svalue(imputation.irmi.noise),", ",
										 'noise.factor=',svalue(imputation.irmi.noise.factor),") ", collapse="")
				addScriptLine(cmd)
				gmessage("Imputation successful!", title="Success", icon="info") 
      }
    }
    #perform hotdeck imputation
    if (svalue(imputation.notebook)==3){
      #convert widgets to useable function parameters, i.e. convert from character to numeric,
			#substitute NULL for empty elements,...
      variable <- getVariableNames(getVm("activeDataSetOriginal"))
      v <- as.character(svalue(imputation.hotdeck.variable))
      if (length(v) > 0){
        variable <- v
      }
      
      ord_var <- NULL
      v <- as.character(svalue(imputation.hotdeck.ord_var))
      if (length(v) > 0){
        ord_var <- v
      }
      
      domain_var <- NULL
      v <- as.character(svalue(imputation.hotdeck.domain_var))
      if (length(v) > 0){
        domain_var <- v
      }
      
      #perform actual imputation
			#capture different errors and warnings
      sumText <- capture.output(impData <- tryCatch({hotdeck(dataset,
                                                   variable = variable,
                                                   ord_var = ord_var,
                                                   domain_var = domain_var,
                                                   #makeNA = makeNA,
                                                   impNA = svalue(imputation.hotdeck.impNA))},
                                                    error=function(e){
                                                      message(e$message)
                                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                                      return(NULL)
                                                    },
                                                    warning=function(e){
                                                      message(e$message)
                                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                                      return(NULL)
                                                    }))
      #no error while imputation
      if(!is.null(impData)) {
				#save imputed data
        putVm("activeDataSetImputed", impData)  
				#find actual imputed variables and save them for future reference (undo)
				if (is.survey(impData)){
					variable <- compareImputations(dataset$variables, impData$variables)
				}
				else{
					variable <- compareImputations(dataset, impData)
				}
				enabled(impVis.plotImputed) <- TRUE 
				vars <- getVm("ImputedVariables")
				vars[variable,1] <- TRUE
				vars[variable,2] <- "hotdeck"
				vars[variable,3] <- format(Sys.time(), "%H:%M:%S")
				putVm("ImputedVariables", vars)
				updateUndoVariablesTable()
				sumtable <- createSummaryDataframe(impData)
				insertTable(imputation.summarytable, sumtable)
				#build command-line script for script browser
				#convert content of widgets to readable string
				variable <- getVariableNames(getVm("activeDataSetOriginal"))
				v <- as.character(svalue(imputation.hotdeck.variable))
				if (length(v) > 0){
					variable <- v
				}
				variable <- paste("c(",paste(sapply(variable, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				
				ord_var <- "NULL"
				v <- as.character(svalue(imputation.hotdeck.ord_var))
				if (length(v) > 0){
					ord_var <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				
				domain_var <- "NULL"
				v <- as.character(svalue(imputation.hotdeck.domain_var))
				if (length(v) > 0){
					domain_var <- paste("c(",paste(sapply(v, FUN=function(s)paste('"',s,'"',sep="")), collapse=","),")")
				}
				cmd <- paste('activedataset <- hotdeck(activedataset, variable=',variable,',',
										 'ord_var=',ord_var,", ",
										 'domain_var=',domain_var,", ",
										 'impNA=',svalue(imputation.hotdeck.impNA),")", collapse="")
				addScriptLine(cmd)
				gmessage("Imputation successful!", title="Success", icon="info")
          
        
      }
    }
    
    #perform regression imputation
    if (svalue(imputation.notebook)==4){
      #convert widgets to useable function parameters, i.e. convert from character to numeric,
      #substitute NULL for empty elements,...
      formula <- paste(svalue(imputation.regression.dependent),
                                  "~",
                                  svalue(imputation.regression.independent))
      
      #perform actual imputation
      #capture different errors and warnings
      sumText <- capture.output(impData <- tryCatch({regressionImp(as.formula(formula),
                                                                   family=svalue(imputation.regression.family),
                                                                   robust=svalue(imputation.regression.robust),
                                                                   data=dataset)},
                                                    error=function(e){
                                                      message(e$message)
                                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                                      return(NULL)
                                                    },
                                                    warning=function(e){
                                                      message(e$message)
                                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                                      return(NULL)
                                                    }))
      #no error while imputation
      if(!is.null(impData) ) {
        #save imputed data
        putVm("activeDataSetImputed", impData)  
        enabled(impVis.plotImputed) <- TRUE 
        vars <- getVm("ImputedVariables")
        
        #find actual imputed variables and save them for undo operation
        if (is.survey(impData)){
          variable <- compareImputations(dataset$variables, impData$variables)
        }
        else{
          variable <- compareImputations(dataset, impData)
        }
        vars[variable,1] <- TRUE
        vars[variable,2] <- "regression"
        vars[variable,3] <- format(Sys.time(), "%H:%M:%S")
        putVm("ImputedVariables", vars)
        updateUndoVariablesTable()
        sumtable <- createSummaryDataframe(impData)
        insertTable(imputation.summarytable, sumtable)
        #build command line string
        #convert the content of widgets to a actual string representation
        cmd <- paste('activedataset <- regressionImp(',as.character(formula), ', family="',svalue(imputation.regression.family),'",',
                     'robust=',svalue(imputation.regression.robust),", ",
                     'data=activedataset)"',collapse="")
        addScriptLine(cmd)
        gmessage("Imputation successful!", title="Success", icon="info") 
      }
    }
  }
  
	#recalculates the content of the data overview table in the data tab
  #handler for the different dropdown-lists of the data overview
  dataPanelChangeHandler <- function(...){
		#chooses if to use the original or the imputed dataset, depending on the users choice in gui
    curdat <- numeric()
    if (svalue(dataPanel.imputationSelection)=="original"){
      curdat <- getVm("activeDataSetOriginal")
    }
    else{
      curdat <- getVm("activeDataSetImputed")
    }
    #if a variable is selected, for which there should be a calculation
    if (is.null(svalue(dataPanel.variableSelection))==FALSE){
      #if a by-operation (i.e. calculating statistic for different domains) should be performed
      if (!is.Empty(svalue(dataPanel.bySelection))){
				#data is survey, use survey package for calculation
        if (is.survey(curdat)) {
          tab <- as.data.frame(svyby(as.formula(paste("~",svalue(dataPanel.variableSelection))),
                                     as.formula(paste("~",svalue(dataPanel.bySelection))),
                                     curdat, 
                                     get(svalue(dataPanel.statisticsSelection))), stringsAsFactors=FALSE)
          tab <- data.frame(lapply(tab, as.character), stringsAsFactors=FALSE)
          insertTable(dataPanel.table, tab)
        }
        else{
          #not a survey, use aggregate for calculation
          tab <- aggregate(as.formula(paste(svalue(dataPanel.variableSelection),"~",svalue(dataPanel.bySelection))),
                           curdat, get(svalue(dataPanel.statisticsSelection)))
					#calculate the standard deviation if the mean was calculated
          if (svalue(dataPanel.statisticsSelection) == "mean"){
						#again use aggregate for the SE of mean calculation
            se <- aggregate(as.formula(paste(svalue(dataPanel.variableSelection),"~",svalue(dataPanel.bySelection))),
                            curdat, FUN = function(s){sd(as.numeric(s),na.rm = TRUE)/sqrt(length(s))})
            se <- se[,2]
          }
          else{
            se <- rep(" ", nrow(tab))
          }

          insertTable(dataPanel.table, cbind(tab, se))
        }
      }#a single value (not a by-oepration) shall be performed
      else{
				#perform single calculation
				#for surveys use survey package to do so
        if (is.survey(curdat)) {
          tab <- do.call(svalue(dataPanel.statisticsSelection), list(as.formula(paste("~",svalue(dataPanel.variableSelection))), curdat))
          insertTable(dataPanel.table, data.frame(domain="all",value=coef(tab), se=SE(tab)))
        }
        else{
          curdat <- curdat[,svalue(dataPanel.variableSelection)]
          tab <- do.call(svalue(dataPanel.statisticsSelection), list(curdat, na.rm = TRUE))
					#calculate SA if mean was previously calculated
          if (svalue(dataPanel.statisticsSelection) == "mean"){
            se <- sd(as.numeric(curdat), na.rm = TRUE)/sqrt(length(curdat))
          }
          else{
            se <- " "
          }
          insertTable(dataPanel.table, data.frame(domain="all",value=tab, se=se))
        }
      }
    }
  }
  
  #opens a dialog to covert a normal data.frame to a survey object
  #allows to select the different weighting variables
  #called when clicking on the corresponding menu item
  createSurveyDialog <- function(){
		#find all datasets in the global environment which can be used as base for new survey
		#alternative use the currently loaded dataset
    vardt <- ls(envir = .GlobalEnv, all.names=TRUE)
    vards <- names(which(sapply(vardt, function(.x) is.data.frame(get(.x)))))
    vards <- c(vards,names(which(sapply(vardt, function(.x) is.survey(get(.x))))))
    vards <- c("Current Dataset",vards)
		#set default focus to the ID textfield, uses to determine the target of the text input actions
    putVm("CreateSurveyDialogFocus", "ids")
		#use per default the already loaded dataset
    dataobject <- getVm("activeDataSetOriginal")
		#build layout and init widgets
    variableNames <- getVariableNames(dataobject)
    survey.window <- gwindow("Create Survey Object from Dataset")
    survey.layout <- glayout(container=survey.window)
    survey.variables <- gtable(variableNames)
    size(survey.variables) <- c(120,-1)
    names(survey.variables) <- "Variables"
    survey.dataset <- gdroplist(vards)
    bg <- ggroup()
    survey.button1 <- gbutton("+", container=bg)
    survey.button2 <- gbutton(":", container=bg)
    survey.button3 <- gbutton("*", container=bg)
    survey.button4 <- gbutton(",", container=bg)
    survey.button5 <- gbutton("^", container=bg)
    survey.button6 <- gbutton("-", container=bg)
    s <- c(25,-1)
    size(survey.button1) <- s
    size(survey.button2) <- s
    size(survey.button3) <- s
    size(survey.button4) <- s
    size(survey.button5) <- s
    size(survey.button6) <- s
    survey.ids <- gedit("1")
    survey.probs <- gedit()
    survey.strata <- gedit()
    survey.vars <- gedit(paste(variableNames, sep=" + "))
    survey.fpc <- gedit()
    survey.weights <- gedit()
    gg <- ggroup()
    survey.discard <- gbutton(" Discard", container=gg)
    survey.accept <- gbutton("Accept", container=gg)
    size(survey.discard) <- c(100, -1 )
    size(survey.accept) <- c(100, -1 )
    survey.layout[1,1, anchor=c(0,0)] <- glabel("original Dataset:")
    survey.layout[1,2:3] <- survey.dataset
    survey.layout[2:11,1] <- survey.variables
    survey.layout[2,3] <- bg
    survey.layout[3,2, anchor=c(0,0)] <- glabel("ids ~")
    survey.layout[4,2, anchor=c(0,0)] <- glabel("props ~")
    survey.layout[5,2, anchor=c(0,0)] <- glabel("strata ~")
    survey.layout[6,2, anchor=c(0,0)] <- glabel("variables ~")
    survey.layout[7,2, anchor=c(0,0)] <- glabel("fpc ~")
    survey.layout[8,2, anchor=c(0,0)] <- glabel("weights ~")
    survey.layout[3,3] <- survey.ids
    survey.layout[4,3] <- survey.probs
    survey.layout[5,3] <- survey.strata
    survey.layout[6,3] <- survey.vars
    survey.layout[7,3] <- survey.fpc
    survey.layout[8,3] <- survey.weights
    survey.layout[9,3] <- gg
    
		#set the background of the in-foxus field to a light yellow
    setWidgetBgColor(survey.ids, "palegoldenrod")
    
    #resets the background color of all entry fields
		#small helper
    resetSurveyColors <- function(){
      setWidgetBgColor(survey.ids, "white")
      setWidgetBgColor(survey.probs, "white")
      setWidgetBgColor(survey.strata, "white")
      setWidgetBgColor(survey.vars, "white")
      setWidgetBgColor(survey.fpc, "white")
      setWidgetBgColor(survey.weights, "white")
    }
    
		#use focus handler of all text-fields to determine current focus
		#save it and recolor backgrounds dependently
    addHandlerFocus(survey.ids, handler=function(h,...){
      putVm("CreateSurveyDialogFocus", "ids")
      resetSurveyColors()
      setWidgetBgColor(survey.ids, "palegoldenrod")
      })
    addHandlerFocus(survey.probs, handler=function(h,...){
      putVm("CreateSurveyDialogFocus", "probs")
      resetSurveyColors()
      setWidgetBgColor(survey.probs, "palegoldenrod")
      })
    addHandlerFocus(survey.strata, handler=function(h,...){
      putVm("CreateSurveyDialogFocus", "strata")
      resetSurveyColors()
      setWidgetBgColor(survey.strata, "palegoldenrod")
      })
    addHandlerFocus(survey.vars, handler=function(h,...){
      putVm("CreateSurveyDialogFocus", "vars")
      resetSurveyColors()
      setWidgetBgColor(survey.vars, "palegoldenrod")
      })
    addHandlerFocus(survey.fpc, handler=function(h,...){
      putVm("CreateSurveyDialogFocus", "fpc")
      resetSurveyColors()
      setWidgetBgColor(survey.fpc, "palegoldenrod")
      })
    addHandlerFocus(survey.weights, handler=function(h,...){
      putVm("CreateSurveyDialogFocus", "weights")
      resetSurveyColors()
      setWidgetBgColor(survey.weights, "palegoldenrod")
      })
			
    #if double click on table add selected variable to textfield which had focus last
		#handler for the variable table
    addHandlerChanged(survey.variables, handler=function(h,...){
			#retrieve saved focus
      field <- getVm("CreateSurveyDialogFocus")
			#insert the variable name to field in a way that a valid R formula is created
			#i.e. for each variable after the first add a "+" before
      if (field=="ids") {
        old <- svalue(survey.ids)
        t <- ""
        if (!isEmpty(old) & !endsWithSymbol(old)){
          t <- "+"
        } 
        svalue(survey.ids) <- paste(svalue(survey.ids),t, svalue(survey.variables), sep="")
      }
      if (field=="probs"){
        old <- svalue(survey.probs)
        t <- ""
        if (!isEmpty(old) & !endsWithSymbol(old)){
          t <- "+"
        } 
        svalue(survey.probs) <- paste(svalue(survey.probs),t, svalue(survey.variables), sep="")
      } 
      if (field=="strata") {
        old <- svalue(survey.strata)
        t <- ""
        if (!isEmpty(old) & !endsWithSymbol(old)){
          t <- "+"
        } 
        svalue(survey.strata) <- paste(svalue(survey.strata),t, svalue(survey.variables), sep="")
      }
      if (field=="vars"){
        old <- svalue(survey.vars)
        t <- ""
        if (!isEmpty(old) & !endsWithSymbol(old)){
          t <- "+"
        } 
        svalue(survey.vars) <- paste(svalue(survey.vars),t, svalue(survey.variables), sep="")
      } 
      if (field=="fpc") {
        old <- svalue(survey.fpc)
        t <- ""
        if (!isEmpty(old) & !endsWithSymbol(old)){
          t <- "+"
        } 
        svalue(survey.fpc) <- paste(svalue(survey.fpc),t, svalue(survey.variables), sep="")
      }
      if (field=="weights") {
        old <- svalue(survey.weights)
        t <- ""
        if (!isEmpty(old) & !endsWithSymbol(old)){
          t <- "+"
        } 
        svalue(survey.weights) <- paste(svalue(survey.weights),t, svalue(survey.variables), sep="")
      }
    })
		
    #if new dataset is selected in the dropdown menu
		#update the variable table
    addHandlerChanged(survey.dataset, handler=function(h,...){
      #current VIMGUI dataset
      if(svalue(survey.dataset, index=TRUE) == 1){
        dataobject <- getVm("activeDataSetOriginal")
        variableNames <- getVariableNames(dataobject)
        insertTable(survey.variables, variableNames)
				svalue(survey.vars) <- paste(variableNames, collapse=" + ")
      }
      else{
        dataobject <- get(svalue(survey.dataset))
        variableNames <- getVariableNames(dataobject)
        insertTable(survey.variables, variableNames)
				svalue(survey.vars) <- paste(variableNames, collapse=" + ")
      }
    })
		
    #discard button, close window
    addHandlerChanged(survey.discard, handler=function(h,...){dispose(survey.window)})
		
    #accept button was clicked, creates survey object and sets it as new active object for GUI
    addHandlerChanged(survey.accept, handler=function(h,...){
			#if survey object is based on current active dataset and there are imputed values
			#use them?
      if(svalue(survey.dataset, index=TRUE) == 1){
        if (is.null(getVm("activeDataSetImputed")) == TRUE){
          dataobject <- getVm("activeDataSetOriginal")
        }
        else{
          w <- gconfirm("Do you want to use the imputed values?", title="Imputed Values", icon="question")
          if (w == TRUE){
            dataobject <- getVm("activeDataSetImputed")
          }
          else{
            dataobject <- getVm("activeDataSetOriginal")
          }
        }
      }
      else{
        dataobject <- get(svalue(survey.dataset))
      }
      if(is.survey(dataobject)){
        dataobject <- dataobject$variables
      }
      #in case of already imputed values, remove "imp" columns
      dataobject <- dataobject[,grep("_imp", colnames(dataobject), invert=TRUE)]
      ids <- NULL
      probs <- NULL
      strata <- NULL
      variables <- NULL
      fpc <- NULL
      weights <- NULL
			#build formulas out of widget strings
      if (is.Empty(svalue(survey.ids)) == FALSE){
        ids <- as.formula(paste("~", svalue(survey.ids)))
      }
      if (is.Empty(svalue(survey.probs)) == FALSE){
        probs <- as.formula(paste("~", svalue(survey.probs)))
      }
      if (is.Empty(svalue(survey.strata)) == FALSE){
        strata <- as.formula(paste("~", svalue(survey.strata)))
      }
      if (is.Empty(svalue(survey.vars)) == FALSE){
        variables <- as.formula(paste("~", svalue(survey.vars)))
      }
      if (is.Empty(svalue(survey.fpc)) == FALSE){
        fpc <- as.formula(paste("~", svalue(survey.fpc)))
      }
      if (is.Empty(svalue(survey.weights)) == FALSE){
        weights <- as.formula(paste("~", svalue(survey.weights)))
      }
			#perform actual survey creation
      surveyObject <- tryCatch({svydesign(ids=ids, probs = probs, strata = strata, variables = variables,
                                         fps = fpc, weights = weights, data = dataobject)},
                               error=function(e){
                                 message(e$message)
                                 gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                 return(NULL)
                               },
                               warning=function(e){
                                 message(e$message)
                                 gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                 return(NULL)
                               })
      if(!is.null(surveyObject)){
        #survey could be created
				#build commandline string for script dialog
        if(svalue(survey.dataset, index=TRUE) == 1){
          datasetname <- "activedataset"
        }
        else{
          datasetname <- svalue(survey.dataset)
        }
				#convert widget content to readable strings
        ids <- "NULL"
        probs <- "NULL"
        strata <- "NULL"
        variables <- "NULL"
        fpc <- "NULL"
        weights <- "NULL"
        if (is.Empty(svalue(survey.ids)) == FALSE){
          ids <- paste("~", svalue(survey.ids))
        }
        if (is.Empty(svalue(survey.probs)) == FALSE){
          probs <- paste("~", svalue(survey.probs))
        }
        if (is.Empty(svalue(survey.strata)) == FALSE){
          strata <- paste("~", svalue(survey.strata))
        }
        if (is.Empty(svalue(survey.vars)) == FALSE){
          variables <- paste("~", svalue(survey.vars))
        }
        if (is.Empty(svalue(survey.fpc)) == FALSE){
          fpc <- paste("~", svalue(survey.fpc))
        }
        if (is.Empty(svalue(survey.weights)) == FALSE){
          weights <- paste("~", svalue(survey.weights))
        }
        cmd <- paste('activedataset <- svydesign(ids=',ids,' ,',
                     'probs=',probs,' ,',
                     'strata=',strata,' ,',
                     'variables=',variables,' ,',
                     'fpc=',fpc,' ,',
                     'weights=',weights,' ,',
                     'data=',datasetname,' )')
        setActiveDataset(surveyObject, loadScript=cmd)
        dispose(survey.window)
      }

    })
    
		#handler for the small symbol buttons, to add strings 
		#just adds the corresponding string to the focused field
    surveyButtonHandler <- function(h,...){
      field <- getVm("CreateSurveyDialogFocus")
      if (field=="ids") {
        svalue(survey.ids) <- paste(svalue(survey.ids), h$action, sep="")
      }
      if (field=="probs"){
        svalue(survey.probs) <- paste(svalue(survey.probs), h$action, sep="")
      } 
      if (field=="strata") {
        svalue(survey.strata) <- paste(svalue(survey.strata), h$action, sep="")
      }
      if (field=="vars"){
        svalue(survey.vars) <- paste(svalue(survey.vars), h$action, sep="")
      } 
      if (field=="fpc") {
        svalue(survey.fpc) <- paste(svalue(survey.fpc), h$action, sep="")
      }
      if (field=="weights") {
        svalue(survey.weights) <- paste(svalue(survey.weights), h$action, sep="")
      }
    }
    
		#add handlers for small buttons
    addHandlerClicked(survey.button1, action="+", handler=surveyButtonHandler)
    addHandlerClicked(survey.button2, action=":", handler=surveyButtonHandler)
    addHandlerClicked(survey.button3, action="*", handler=surveyButtonHandler)
    addHandlerClicked(survey.button4, action=",", handler=surveyButtonHandler)
    addHandlerClicked(survey.button5, action="^", handler=surveyButtonHandler)
    addHandlerClicked(survey.button6, action="-", handler=surveyButtonHandler)
  }
  
  #opens a dialog which allows the usage of the prepare method 
  #of the VIM package onto the current dataset
  #called after pressing the corresponding menu item
  createPrepareDialog <- function(){
		#inits the window, its layout and the different widgets
    prepare.window <- gwindow("Prepare Dataset", width=100, height=100)
    prepare.scaling <- gdroplist(c("none","classical","MCD","robust","onestep"))
    prepare.transformation <- gdroplist(c("none","minus","reciprocal","logarithm",
                                          "exponential","boxcox","clr","ilr","alr"))
    prepare.alpha <- gedit()
    prepare.powers <- gedit()
    prepare.start <- gedit()
    prepare.apply <- gbutton(" Apply")
    prepare.alrVar <- gdroplist(getVariableNames(getVm("activeDataSetOriginal")))
    prepare.layout <- glayout(container=prepare.window)
    prepare.layout[1,1] <- glabel("scaling:")
    prepare.layout[1,2] <- prepare.scaling
    prepare.layout[1,3] <- glabel("transform:")
    prepare.layout[1,4] <- prepare.transformation
    prepare.layout[2,1] <- glabel("start:")
    prepare.layout[2,2] <- prepare.start
    prepare.layout[2,3] <- glabel("alpha:")
    prepare.layout[2,4] <- prepare.alpha
    prepare.layout[3,1] <- glabel("powers:")
    prepare.layout[3,2] <- prepare.powers
    prepare.layout[3,3] <- glabel("alrVar:")
    prepare.layout[3,4] <- prepare.alrVar
    prepare.layout[4,4] <- prepare.apply
    
		#disable elements not working with current setting
    enabled(prepare.alpha) <- FALSE
    enabled(prepare.powers) <- FALSE
    enabled(prepare.start) <- FALSE
    enabled(prepare.alrVar) <- FALSE
    
    #gray out not usable parameter
    addHandlerChanged(prepare.scaling, function(h,..){
      if (svalue(prepare.scaling) == "MCD"){
        enabled(prepare.alpha) <- TRUE
      }
      else{
        enabled(prepare.alpha) <- FALSE
      }
    })
    
    #gray out not usable parameter
    addHandlerChanged(prepare.transformation, function(h,..){
      enabled(prepare.powers) <- FALSE
      enabled(prepare.start) <- FALSE
      enabled(prepare.alrVar) <- FALSE
      if (svalue(prepare.transformation) == "boxcox"){
        enabled(prepare.powers) <- TRUE
        enabled(prepare.start) <- TRUE
      } else if (svalue(prepare.transformation) == "alr"){
        enabled(prepare.alrVar) <- TRUE
      }
    })
    
    #handler for apply button
		#does the actual prepare
    addHandlerClicked(prepare.apply, function(h,...){
      proceed <- TRUE
			w <- FALSE
      #dataset already imputed, proceed?
      if (is.null(getVm("activeDataSetImputed")) == FALSE){
        w <- gconfirm("Preparing the dataset will remove imputation, proceed?", title="Imputed Values", icon="question")
				if (w == TRUE){
					proceed <- TRUE
				}
      }
      if(proceed){
				#depending on chosen modes, select parameters from gui widgets
        if (w == TRUE) {
					currentDataset <- getVm("activeDataSetImputed")
				}
				else
				{
					currentDataset <- getVm("activeDataSetOriginal")
				}
        alpha <- NULL
        powers <- NULL
        start <- NULL
        alrVar <- NULL
        if (svalue(prepare.scaling) == "MCD"){
          alpha <- as.numeric(svalue(prepare.alpha))
        }
        if (svalue(prepare.transformation) == "boxcox"){
          powers <- as.numeric(svalue(prepare.powers))
          start <- as.numeric(svalue(prepare.start))
        } 
        if (svalue(prepare.transformation) == "alr"){
          alrVar <- svalue(prepare.alrVar) 
        }
        #try to prepare dataset
        preparedDataset <- tryCatch({prepare(currentDataset, scaling=svalue(prepare.scaling),
                                    transformation = svalue(prepare.transformation),
                                    alpha=alpha, powers=powers, start=start, alrVar=alrVar)},
                                    error=function(e){
                                      message(e$message)
                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                      return(NULL)
                                    },
                                    warning=function(e){
                                      message(e$message)
                                      gmessage(paste("A problem occurred (see also in console):",e$message), title="Problem", icon="error")
                                      return(NULL)
                                    })
        #there was no problem
        if(!is.null(preparedDataset)) {
					#save old dataset for possible undo operation
          putVm("undoDataSetOriginal", currentDataset)
					#build command line script for script browser
          if (is.null(alpha)) alpha <- "NULL"
          if (is.null(powers)) powers <- "NULL"
          if (is.null(start)) start <- "NULL"
          if (is.null(alrVar)) {
            alrVar <- "NULL"
          }
          else{
            alrVar <- paste('"',alrVar,'"', sep="")
          }
          cmdimp <- paste("undoDataset <- activedataset \n","activedataset <- prepare(activedataset, ",
                          'scaling="',svalue(prepare.scaling),'", ',
                          'transformation="',svalue(prepare.transformation),'", ',
                          'alpha=',alpha,', ',
                          'powers=',powers,', ',
                          'start=',start,', ',
                          'alrVar=',alrVar,')', sep="")
          #addScriptLine(cmdimp)
					#set new prepared dataset as the currently active dataset
					#dont adjustTypes as the old dataset must have been
          setActiveDataset(preparedDataset, prepare=TRUE, 
                           adjustTypes=FALSE, loadScript=cmdimp)
          enabled(menu.Undo) <- TRUE
          dispose(prepare.window)
        }
        
      }
    })
  }
  
  #undo a previous prepare action
  #called from corresponding menu entry
  undoPrepare <- function(){
		#Do you really want to?
    w <- gconfirm("Undoing the previous prepare action removes all subsequent changes, proceed?", title="Undo prepare", icon="question")
    if (w == TRUE){
      undoDataset <- getVm("undoDataSetOriginal")
      setActiveDataset(undoDataset, adjustTypes=FALSE, deleteScript=FALSE)
      addScriptLine("activedataset <- undoDataset")
    }
  }
  
  #updates the content of the undo variables table
  #called after imputation occurred
  updateUndoVariablesTable <- function(){
    vars <- getVm("ImputedVariables")
    elements <- vars[vars[,1]==TRUE,]
    insertTable(imputation.undo.variables, cbind(rownames(elements), elements[,2], elements[,3]))
  }
  
  #opens a dialog to convert specific values in specific variables to NA
	#called after clicking the corresponding menu entry
  createNaDialog <- function(){
    setNA.window <- gwindow("Set Values to NA")
    setNA.panel <- glayout(container=setNA.window, expand=TRUE)
    setNA.variables <- gtable(getVariableNames(getVm("activeDataSetOriginal")), multiple=TRUE)
    names(setNA.variables) <- "variables"
    setNA.values <- gedit()
    setNA.button <- gbutton("Set Values to NA")
    setNA.panel[1:8,1, expand=TRUE] <- setNA.variables
    setNA.panel[9,1, anchor=c(-1,-1)] <- glabel("Values: ")
    setNA.panel[10,1] <- setNA.values
    setNA.panel[11,1] <- setNA.button
    #add Handler to button
    addHandlerClicked(setNA.button, function(h,...){
      curDat <- getVm("activeDataSetOriginal")
      vals <- cutParam(svalue(setNA.values))
      if (is.survey(curDat)){
        cols <- curDat$variables[,svalue(setNA.variables)]
        for (s in vals){
          cols[cols==s] <- NA
        }
        curDat$variables[,svalue(setNA.variables)] <- cols
        #putVm("activeDataSetOriginal", curDat)
        setActiveDataset(curDat, adjustTypes=FALSE)
      }
      else{
        cols <- curDat[,svalue(setNA.variables)]
        for (s in vals){
          cols[cols==s] <- NA
        }
        curDat[,svalue(setNA.variables)] <- cols
        #putVm("activeDataSetOriginal", curDat)
        setActiveDataset(curDat, adjustTypes=FALSE)
      }
    })
  }
  
	#updates some tables with possible values combinations
  #called if certain widgets in the imputation tabs are changed like the tables
  imputationChangeHandler <- function(h,...){
    if (svalue(imputation.notebook)==1){
      #kNN tab changed
      insertTable(imputation.kNN.mixed, as.character(svalue(imputation.kNN.variable)))
    }
    if (svalue(imputation.notebook)==2){
      #irmi tab changed
      insertTable(imputation.irmi.mixed, as.character(svalue(imputation.irmi.variable)))
    }
  }
  
  #writes selected variable names into the formula widgets of the regression imputation tab
	#called  after double click on table for regression imputation
  regressionTableHandler <- function(h,...){
		#retrieve previously saved focus (dependent or independent variable)
    focus <- getVm("regressionFocus")
		#insert variable into last focused field
    if (focus == 0){
      old <- svalue(imputation.regression.dependent)
      t <- ""
			#test if a "+" is necessary for a valid formula representation
			#i.e. if this is the first variable in the line or the character before is special 
      if (!isEmpty(old) & !endsWithSymbol(old)){
        t <- "+"
      }
      text <- paste(svalue(imputation.regression.dependent),t,
                    svalue(imputation.regression.variables), sep = "")
                    svalue(imputation.regression.dependent) <- text
    }
    else{
      old <- svalue(imputation.regression.independent)
      t <- ""
			#test if a "+" is necessary for a valid formula representation
			#i.e. if this is the first variable in the line or the character before is special 
      if (!isEmpty(old) & !endsWithSymbol(old)){
        t <- "+"
      }
      text <- paste(svalue(imputation.regression.independent),t,
                    svalue(imputation.regression.variables), sep = "")
                    svalue(imputation.regression.independent) <- text
    }
  }
  
  #called when clicking on the symbol buttons in the regression tab
  #adds the specified symbol
  regressionButtonHandler <- function(h,...){
    focus <- getVm("regressionFocus")
    if (focus == 0){
      text <- paste(svalue(imputation.regression.dependent), h$action, sep = "")
      svalue(imputation.regression.dependent) <- text
    }
    else{
      text <- paste(svalue(imputation.regression.independent),h$action, sep = "")
      svalue(imputation.regression.independent) <- text
    }
  }
  
  
  ########
  #Initializations for the most part of the user interface
	#mostly layout and widget creation
	
  #data(sleep)
  #putVm("activeDataSetOriginal", sleep)
  #putVm("activeDataSetImputed", NULL)
	
	#save default colors
  putVm("plotColors", c("skyblue","red","skyblue4","red4","orange","orange4"))
  putVm("plotAlpha", 1)
  
  #create menu structure for main window
  menu.LoadDataset <- gaction(label="Load Dataset", handler=loadDataSet)
  menu.ChooseDataset <- gaction(label="Choose Dataset", handler=setDataSet)
  menu.SaveFile <- gaction(label="to File", handler=saveToFile)
  menu.SaveVariable <- gaction(label="to Variable", handler=saveToVariable)
  menu.ImportCSV <- gaction(label="CSV", handler=importCSV)
  menu.ExportCSV <- gaction(label="CSV", handler=exportCSV)
  menu.ImportSPSS <- gaction(label="SPSS", handler=importSPSS)
  menu.ExportSPSS <- gaction(label="SPSS", handler=exportSPSS)
  menu.ImportSTATA <- gaction(label="STATA", handler=importSTATA)
  menu.ExportSTATA <- gaction(label="STATA", handler=exportSTATA)
  menu.ImportSAS <- gaction(label="SAS", handler=importSAS)
  menu.ExportSAS <- gaction(label="SAS", handler=exportSAS)
  menu.CreateSurvey <- gaction(label="Create Survey", handler=function(h,...) createSurveyDialog())
  menu.SetNA <- gaction(label="Set Value NA", handler = function(h,...) createNaDialog())  
  menu.Prepare <- gaction(label="Prepare Dataset", handler = function(h,...) createPrepareDialog()) 
  menu.Undo <- gaction(label="Undo Prepare", handler = function(h,...) undoPrepare()) 
  menu.Preferences <- gaction(label="Preferences", handler = OptionsHandler) 
  menu.Script <- gaction(label="Script", handler = ScriptHandler) 
  
	#create main window
  mainWindow <- gwindow(title="VIM GUI", width=1024, height=768)
 
	#build menu bar
  ml <- list(
    Data = list(
    load=menu.LoadDataset,
    choose=menu.ChooseDataset,
    Save=list(
      file=menu.SaveFile,
      variable=menu.SaveVariable),
    Import=list(
      impcsv=menu.ImportCSV,
      impspss=menu.ImportSPSS,
      impsas=menu.ImportSAS,
      impsata=menu.ImportSTATA),
    Export=list(
      expcsv=menu.ExportCSV,
      expspss=menu.ExportSPSS,
      expsas=menu.ExportSAS,
      expstata=menu.ExportSTATA)),
    Survey=list(create=menu.CreateSurvey),
    Edit=list(na=menu.SetNA, prep=menu.Prepare, undo=menu.Undo),
    Miscellaneous=list(script=menu.Script, prefs=menu.Preferences)
    )
  menuBar <- gmenu(ml, container=mainWindow)
  
  #add notebook with tabs for different tasks to main window
  mainNotebook <- gnotebook(container=mainWindow)
  dataGroup <- ggroup(container=mainNotebook, label="Data")
  imputationGroup <- ggroup(horizontal=FALSE, container=mainNotebook, label="Imputation")
  imputationVisGroup <- ggroup(container=mainNotebook, label="Visualization")
  addHandlerChanged(mainNotebook, mainNotebook.handler)
  
  #####
  #init and layout of the Data tab
  dataPanel.summaryframe <- gframe("Summary of original dataset: ", container=dataGroup, expand=TRUE)
  dataPanel.summarypanel <- ggroup(horizontal=FALSE, container=dataPanel.summaryframe, expand=TRUE)
  dataPanel.summarytext <- gtext(container=dataPanel.summarypanel)
  size(dataPanel.summarytext) <- c(-1,100)
  df <- data.frame(a="", b="", c="", d="", e="", f="", g="", h="", stringsAsFactors=FALSE)
  dataPanel.summaryTable <- gtable(df, container=dataPanel.summarypanel, expand=TRUE)
  names(dataPanel.summaryTable) <- c("Name", "Class", "NA", "Min","Lower","Median","Upper","Max")
  enabled(dataPanel.summarytext) <- FALSE
  dataPanel.overviewframe <- gframe("Overview: ", container=dataGroup)
  dataPanel.overviewpanel <- glayout(container=dataPanel.overviewframe)
  dataPanel.variableSelection <- gdroplist("")
  dataPanel.statisticsSelection <- gdroplist(c("svymean", "svyvar", "svytotal"))
  dataPanel.bySelection <- gdroplist(" ")
  size(dataPanel.variableSelection) <- c(120,-1)
  size(dataPanel.bySelection) <- c(120,-1)
  dataPanel.imputationSelection <- gradio(c("original","imputed"), horizontal=TRUE)
  dataPanel.table <- gtable(data.frame(Domain=rep(" ",1000), Value=rep(" ",1000), SE=rep(" ",1000)))
  dataPanel.overviewpanel[1,1:3] <- dataPanel.imputationSelection 
  dataPanel.overviewpanel[2,1] <- glabel("Variable:")
  dataPanel.overviewpanel[3,1] <- dataPanel.variableSelection
  dataPanel.overviewpanel[2,2] <- glabel("Statistic:")
  dataPanel.overviewpanel[3,2] <- dataPanel.statisticsSelection
  dataPanel.overviewpanel[2,3] <- glabel("By:")
  dataPanel.overviewpanel[3,3] <- dataPanel.bySelection
  dataPanel.overviewpanel[4,1:3, expand=TRUE] <- dataPanel.table
  dataPanel.variableSelectionHandler <- addHandlerChanged(dataPanel.variableSelection, handler=dataPanelChangeHandler)
  dataPanel.statisticsSelectionHandler <- addHandlerChanged(dataPanel.statisticsSelection, handler=dataPanelChangeHandler)
  dataPanel.bySelectionHandler <- addHandlerChanged(dataPanel.bySelection, handler=dataPanelChangeHandler)
  dataPanel.imputationSelectionHandler <- addHandlerChanged(dataPanel.imputationSelection, handler=dataPanelChangeHandler)
 
  
  ####
  #init and layout Imputation tab
  imputation.notebook <- gnotebook(tab.pos = 2, container=imputationGroup, expand=TRUE)
  g <- ggroup(container=imputationGroup)
  f <- gframe("Summary of imputed dataset:", container=g, expand=TRUE)
  df <- data.frame(a="", b="", c="", d="", e="", f="", g="", h="", stringsAsFactors=FALSE)
  imputation.summarytable <- gtable(df, expand=TRUE, container=f)
  names(imputation.summarytable) <- c("Name", "Class", "NA", "Min","Lower","Median","Upper","Max")
  size(imputation.summarytable) <- c(-1,300)
  f2 <- gframe("", container=g)
  gg <- ggroup(horizontal=FALSE, container=f2)
  imputation.ok <- gbutton("Apply Imputation", container=gg)
  addHandlerClicked(imputation.ok, handler=Imputation)
  imputation.kNN.group <- glayout(container=imputation.notebook, label="kNN", expand=TRUE)
  imputation.irmi.group <- glayout(container=imputation.notebook, label="irmi", expand=TRUE)
  imputation.hotdeck.group <- glayout(container=imputation.notebook, label="hotdeck", expand=TRUE)
  imputation.regression.group <- glayout(container=imputation.notebook, label="regression", expand=TRUE)
  imputation.undo.group <- glayout(container=imputation.notebook, label="undo", expand=TRUE)
	
  ###kNN-imputation-Tab
  imputation.kNN.variable <- gtable(rep("",100), multiple=TRUE)
  size(imputation.kNN.variable) <- c(120,200)
  names(imputation.kNN.variable) <- "variables"
  imputation.kNN.k <- gdroplist(1:15)
  imputation.kNN.dist_var <- gtable(rep("",100), multiple=TRUE)
  size(imputation.kNN.dist_var) <- c(120,200)
  names(imputation.kNN.dist_var) <- "dist_var"
  imputation.kNN.weights <- gedit()
  imputation.kNN.numFun <- gdroplist(c("median","mean"))
  imputation.kNN.catFun <- gdroplist(c("maxCat","sampleCat"))
  imputation.kNN.impNA <- gcheckbox("impNA")
  svalue(imputation.kNN.impNA) <- TRUE
  imputation.kNN.addRandom <- gcheckbox("addRandom")
  imputation.kNN.mixed <- gtable(rep("",100), multiple=TRUE)
  size(imputation.kNN.mixed) <- c(120,200)
  names(imputation.kNN.mixed) <- "mixed"
  imputation.kNN.mixed.constant <- gedit()
  imputation.kNN.group[1:8,1] <- imputation.kNN.variable
  imputation.kNN.group[1:8,2] <- imputation.kNN.dist_var
  imputation.kNN.group[1:8,3] <- imputation.kNN.mixed
  imputation.kNN.group[1,4, anchor = c(-1, 0)] <- glabel("k:")
  imputation.kNN.group[1,5] <- imputation.kNN.k
  imputation.kNN.group[1,6, anchor = c(-1, 0)] <- glabel("weights:")
  imputation.kNN.group[1,7] <- imputation.kNN.weights
  imputation.kNN.group[2,4, anchor = c(-1, 0)] <- glabel("numFun:")
  imputation.kNN.group[2,5] <- imputation.kNN.numFun
  imputation.kNN.group[2,6, anchor = c(-1, 0)] <- glabel("catFun:")
  imputation.kNN.group[2,7] <- imputation.kNN.catFun
  imputation.kNN.group[3,4, anchor = c(-1, 0)] <- glabel("mixed.constant:")
  imputation.kNN.group[3,5] <- imputation.kNN.mixed.constant
  imputation.kNN.group[4,4, anchor = c(-1, 0)] <- imputation.kNN.impNA
  imputation.kNN.group[4,5, anchor = c(-1, 0)] <- imputation.kNN.addRandom
  addHandlerClicked(imputation.kNN.variable, handler=imputationChangeHandler)
	
  ###irmi-imputation-Tab
  imputation.irmi.variable <- gtable(rep("",100), multiple=TRUE)
  size(imputation.irmi.variable) <- c(120,200)
  names(imputation.irmi.variable) <- "variables"
  imputation.irmi.mixed <- gtable(rep("",100), multiple=TRUE)
  size(imputation.irmi.mixed) <- c(120,200)
  names(imputation.irmi.mixed) <- "mixed"
  imputation.irmi.mixed.constant <- gedit()
  imputation.irmi.count <- gtable(rep("",100), multiple=TRUE)
  size(imputation.irmi.count) <- c(120,200)
  names(imputation.irmi.count) <- "count"
  imputation.irmi.robust <- gcheckbox("robust")
  imputation.irmi.noise <- gcheckbox("noise")
  imputation.irmi.noise.factor <- gslider(from = 0, to = 1, by = 0.01)
  size(imputation.irmi.noise.factor) <- c(120,-1)
  imputation.irmi.group[1:8,1] <- imputation.irmi.variable
  imputation.irmi.group[1:8,2] <- imputation.irmi.count
  imputation.irmi.group[1:8,3] <- imputation.irmi.mixed
  imputation.irmi.group[1,4, anchor = c(-1, 0)] <- glabel("mixed.constant:")
  imputation.irmi.group[1,5] <- imputation.irmi.mixed.constant
  imputation.irmi.group[2,4] <- imputation.irmi.robust
  imputation.irmi.group[3,4, anchor = c(-1, 0)] <- imputation.irmi.noise
  imputation.irmi.group[4,4, anchor = c(-1, -0.5)] <- glabel("noise.factor:")
  imputation.irmi.group[4,5, anchor = c(-1, 0)] <- imputation.irmi.noise.factor
  addHandlerClicked(imputation.irmi.variable, handler=imputationChangeHandler)
	
  ###hotdeck-imputation-Tab
  imputation.hotdeck.variable <- gtable(rep("",100), multiple=TRUE)
  size(imputation.hotdeck.variable) <- c(120,200)
  names(imputation.hotdeck.variable) <- "variables"
  imputation.hotdeck.ord_var <- gtable(rep("",100), multiple=TRUE)
  size(imputation.hotdeck.ord_var) <- c(120,200)
  names(imputation.hotdeck.ord_var) <- "ord_var"
  imputation.hotdeck.domain_var<- gtable(rep("",100), multiple=TRUE)
  size(imputation.hotdeck.domain_var) <- c(120,200)
  names(imputation.hotdeck.domain_var) <- "domain_var"
  imputation.hotdeck.impNA <- gcheckbox("impNA")
  svalue(imputation.hotdeck.impNA) <- TRUE
  imputation.hotdeck.group[1:8,1] <- imputation.hotdeck.variable
  imputation.hotdeck.group[1:8,2] <- imputation.hotdeck.ord_var
  imputation.hotdeck.group[1:8,3] <- imputation.hotdeck.domain_var
  imputation.hotdeck.group[1,4] <- imputation.hotdeck.impNA
	
  ###regression-imputation-tab
  imputation.regression.variables <- gtable("")
  names(imputation.regression.variables) <- "variables"
  size(imputation.regression.variables) <- c(120,200)
  bg <- ggroup()
  imputation.regression.button1 <- gbutton("+", container=bg)
  imputation.regression.button2 <- gbutton(":", container=bg)
  imputation.regression.button3 <- gbutton("*", container=bg)
  imputation.regression.button4 <- gbutton(",", container=bg)
  imputation.regression.button5 <- gbutton("^", container=bg)
  imputation.regression.button6 <- gbutton("-", container=bg)
  s <- c(25,-1)
  size(imputation.regression.button1) <- s
  size(imputation.regression.button2) <- s
  size(imputation.regression.button3) <- s
  size(imputation.regression.button4) <- s
  size(imputation.regression.button5) <- s
  size(imputation.regression.button6) <- s
  bgg <- ggroup()
  imputation.regression.dependent <- gedit(container=bgg)
  glabel("~", container=bgg)
  imputation.regression.independent <- gedit(container=bgg)
  size(imputation.regression.independent) <- c(250, -1)
  imputation.regression.family <- gdroplist(c("AUTO", "normal", "binomial", "multinomial",
                                              "poisson", "lognormal"))
  imputation.regression.robust <- gcheckbox("robust")
  imputation.regression.group[1:8,1] <- imputation.regression.variables
  imputation.regression.group[1,2] <- bg
  imputation.regression.group[2,2:6] <- bgg
  imputation.regression.group[3,2] <- imputation.regression.family
  imputation.regression.group[4,2] <- imputation.regression.robust
  addHandlerChanged(imputation.regression.variables, handler=regressionTableHandler)
  setWidgetBgColor(imputation.regression.dependent, "palegoldenrod")
	#add handler for the dependent and independent formula fields
	#these save the last focus and change their colors
  addHandlerFocus(imputation.regression.dependent, handler=function(h,...){
    putVm("regressionFocus", 0)
    setWidgetBgColor(imputation.regression.dependent, "palegoldenrod")
    setWidgetBgColor(imputation.regression.independent, "white")
  })
  addHandlerFocus(imputation.regression.independent, handler=function(h,...){
    putVm("regressionFocus", 1)
    setWidgetBgColor(imputation.regression.independent, "palegoldenrod")
    setWidgetBgColor(imputation.regression.dependent, "white")
  })
  addHandlerClicked(imputation.regression.button1, action="+", handler=regressionButtonHandler)
  addHandlerClicked(imputation.regression.button2, action=":", handler=regressionButtonHandler)
  addHandlerClicked(imputation.regression.button3, action="*", handler=regressionButtonHandler)
  addHandlerClicked(imputation.regression.button4, action=",", handler=regressionButtonHandler)
  addHandlerClicked(imputation.regression.button5, action="^", handler=regressionButtonHandler)
  addHandlerClicked(imputation.regression.button6, action="-", handler=regressionButtonHandler)
	
  ###Undo-imputation-Tab
  imputation.undo.variables <- gtable(data.frame(variable="", method="", time=""))
  names(imputation.undo.variables) <- c("variable","method", "time")
  #size(imputation.undo.variables) <- c(120,200)
  imputation.undo.group[1,1:2] <- glabel("Double click on a variable to undo imputation:")
  imputation.undo.group[2:10,1:2, expand=TRUE] <- imputation.undo.variables
  addHandlerChanged(imputation.undo.variables, handler=undoImputation) 

  ####
  #init and layout Visualization - tab
	#save references of handlers in a single list to allow a quick deactivation
	#done for performance reasons
	handlerList <- list()
  g <- ggroup(horizontal=FALSE, container=imputationVisGroup)
  impVis.plotImputed <- gradio(c("original","imputed"), horizontal=TRUE, container=g)
  enabled(impVis.plotImputed) <- FALSE
  addHandlerChanged(impVis.plotImputed, function(h,...){makeImputationPlot()})
  impVis.plotBook <- gnotebook(container=g, tab.pos = 2, expand=TRUE, fill="x")
  size(impVis.plotBook) <- c(260,-1)
  impVis.aggr.group <- glayout(container=impVis.plotBook)
  impVis.barMiss.group <- glayout(container=impVis.plotBook)
  impVis.histMiss.group <- glayout(container=impVis.plotBook)
  impVis.marginMatrix.group <- glayout(container=impVis.plotBook)
  impVis.scattmatrixMiss.group <- glayout(container=impVis.plotBook)
  impVis.mosaicMiss.group <- glayout(container=impVis.plotBook)
  impVis.parcoordMiss.group <- glayout(container=impVis.plotBook)
  impVis.pbox.group <- glayout(container=impVis.plotBook)
  impVis.matrixplot.group <- glayout(container=impVis.plotBook)
  impVis.plot <- gimage(container=imputationVisGroup, expand=TRUE)
    
  #aggregation of imputed values plot
  impVis.aggr.bars <- gcheckbox("bars")
  impVis.aggr.numbers <- gcheckbox("numbers")
  impVis.aggr.prop <- gcheckbox("prop")
  impVis.aggr.combined <- gcheckbox("combined")
  impVis.aggr.only.miss <- gcheckbox("only.miss")
  impVis.aggr.sortVars <- gcheckbox("sortVars")
  impVis.aggr.sortCombs <- gcheckbox("sortCombs")
  impVis.aggr.weighted <- gcheckbox("weighted")
  impVis.aggr.group[1,1] <- impVis.aggr.bars
  impVis.aggr.group[1,2] <- impVis.aggr.numbers
  impVis.aggr.group[2,1] <- impVis.aggr.prop
  impVis.aggr.group[2,2] <- impVis.aggr.combined
  impVis.aggr.group[3,1] <- impVis.aggr.only.miss
  impVis.aggr.group[3,2] <- impVis.aggr.sortVars
  impVis.aggr.group[4,1] <- impVis.aggr.sortCombs
  impVis.aggr.group[4,2] <- impVis.aggr.weighted
  addHandlerChanged(impVis.aggr.bars, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.aggr.numbers, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.aggr.prop, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.aggr.combined, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.aggr.only.miss, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.aggr.sortVars, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.aggr.sortCombs, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.aggr.weighted, function(h,...){makeImputationPlot()})
  
  #barMiss
  impVis.barMiss.pos <- gdroplist("")
  impVis.barMiss.selection <- gdroplist(c("any","all"))
  impVis.barMiss.only.miss <- gcheckbox("only.miss")
  size(impVis.barMiss.selection) <- c(140,-1)
  impVis.barMiss.weighted <- gcheckbox("weighted")
  impVis.barMiss.group[1,1] <- glabel("pos: ")
  impVis.barMiss.group[1,2] <- impVis.barMiss.pos
  impVis.barMiss.group[2,1] <- glabel("selection: ")
  impVis.barMiss.group[2,2] <- impVis.barMiss.selection
  impVis.barMiss.group[3,1] <- impVis.barMiss.only.miss
  impVis.barMiss.group[3,2] <- impVis.barMiss.weighted
  addHandlerChanged(impVis.barMiss.selection, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.barMiss.only.miss, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.barMiss.weighted, function(h,...){makeImputationPlot()})
  handlerList[[length(handlerList)+1]] <- c(impVis.barMiss.pos, 
                                            addHandlerChanged(impVis.barMiss.pos, function(h,...){makeImputationPlot()}))
  #histMiss
  impVis.histMiss.pos <- gdroplist("")
  impVis.histMiss.selection <- gdroplist(c("any","all"))
  impVis.histMiss.breaks <- gcombobox(c("Sturges", "Scott", "Freedman-Diaconis"),editable = TRUE)
  impVis.histMiss.right <- gcheckbox("right")
  impVis.histMiss.only.miss <- gcheckbox("only.miss")
  impVis.histMiss.weighted <- gcheckbox("weighted")
  size(impVis.histMiss.breaks) <- c(140,-1)
  impVis.histMiss.group[1,1] <- glabel("pos: ")
  impVis.histMiss.group[1,2] <- impVis.histMiss.pos
  impVis.histMiss.group[2,1] <- glabel("selection: ")
  impVis.histMiss.group[2,2] <- impVis.histMiss.selection
  impVis.histMiss.group[3,1] <- glabel("breaks: ")
  impVis.histMiss.group[3,2] <- impVis.histMiss.breaks
  impVis.histMiss.group[4,1] <- impVis.histMiss.right
  impVis.histMiss.group[4,2] <- impVis.histMiss.only.miss
  impVis.histMiss.group[5,1] <- impVis.histMiss.weighted
  addHandlerChanged(impVis.histMiss.selection, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.histMiss.breaks, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.histMiss.right, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.histMiss.weighted, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.histMiss.only.miss, function(h,...){makeImputationPlot()})
  handlerList[[length(handlerList)+1]] <- c(impVis.histMiss.pos,
                                            addHandlerChanged(impVis.histMiss.pos, function(h,...){makeImputationPlot()}))
  
  #marginmatrix
  impVis.marginMatrix.plotvars <- gtable(rep("",100),  multiple = TRUE)
  names(impVis.marginMatrix.plotvars) <- "plotvars"
  impVis.marginMatrix.group[1,1:2, fill="y", expand="TRUE"] <- impVis.marginMatrix.plotvars
  addHandlerClicked(impVis.marginMatrix.plotvars, function(h,...){makeImputationPlot()})
  
  #scattmatrixMiss
  impVis.scattmatrixMiss.highlight <- gtable(rep("",100), multiple=TRUE)
  names(impVis.scattmatrixMiss.highlight) <- "highlight"
  impVis.scattmatrixMiss.plotvars <- gtable(rep("",100), multiple=TRUE)
  names(impVis.scattmatrixMiss.plotvars) <- "plotvars"
  impVis.scattmatrixMiss.selection <- gdroplist(c("any","all"))
  impVis.scattmatrixMiss.diagonal<- gdroplist(c("density","none"))
  impVis.scattmatrixMiss.weighted <- gcheckbox("weighted")
  impVis.scattmatrixMiss.group[1,1:2, fill="y", expand="TRUE"] <- impVis.scattmatrixMiss.plotvars
  impVis.scattmatrixMiss.group[2,1:2, fill="y", expand="TRUE"] <- impVis.scattmatrixMiss.highlight
  impVis.scattmatrixMiss.group[3,1] <- glabel("selection: ")
  impVis.scattmatrixMiss.group[3,2] <- impVis.scattmatrixMiss.selection
  impVis.scattmatrixMiss.group[4,1] <- glabel("diagonal: ")
  impVis.scattmatrixMiss.group[4,2] <- impVis.scattmatrixMiss.diagonal
  impVis.scattmatrixMiss.group[5,1] <- impVis.scattmatrixMiss.weighted
  addHandlerClicked(impVis.scattmatrixMiss.plotvars, function(h,...){makeImputationPlot()})
  addHandlerClicked(impVis.scattmatrixMiss.highlight, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.scattmatrixMiss.selection, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.scattmatrixMiss.diagonal, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.scattmatrixMiss.weighted, function(h,...){makeImputationPlot()})
  
  #mosaicplot
  impVis.mosaicMiss.highlight <- gtable(rep("",100), multiple=TRUE)
  names(impVis.mosaicMiss.highlight) <- "highlight"
  impVis.mosaicMiss.plotvars <- gtable(rep("",100), multiple=TRUE)
  names(impVis.mosaicMiss.plotvars) <- "plotvars"
  impVis.mosaicMiss.selection <- gdroplist(c("any","all"))
  impVis.mosaicMiss.weighted <- gcheckbox("weighted")
  impVis.mosaicMiss.group[1,1] <- glabel("selection: ")
  impVis.mosaicMiss.group[1,2] <- impVis.mosaicMiss.selection
  impVis.mosaicMiss.group[2,1] <- impVis.mosaicMiss.weighted
  impVis.mosaicMiss.group[3,1:2, fill="y", expand="TRUE"] <- impVis.mosaicMiss.plotvars
  impVis.mosaicMiss.group[4,1:2, fill="y", expand="TRUE"] <- impVis.mosaicMiss.highlight
  addHandlerChanged(impVis.mosaicMiss.selection, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.mosaicMiss.weighted, function(h,...){makeImputationPlot()})
  h <- addHandlerClicked(impVis.mosaicMiss.highlight, function(h,...){makeImputationPlot()})
  handlerList[[length(handlerList)+1]] <- c(impVis.mosaicMiss.highlight,h)
  handlerList[[length(handlerList)+1]] <- c(impVis.mosaicMiss.plotvars,
                                            addHandlerClicked(impVis.mosaicMiss.plotvars, function(h,...){makeImputationPlot()}))
  #parcoord
  impVis.parcoordMiss.highlight <- gtable(rep("",100), multiple=TRUE)
  names(impVis.parcoordMiss.highlight) <- "highlight"
  impVis.parcoordMiss.plotvars <- gtable(rep("",100), multiple=TRUE)
  names(impVis.parcoordMiss.plotvars) <- "plotvars"
  impVis.parcoordMiss.selection <- gdroplist(c("any","all"))
  impVis.parcoordMiss.plotNA <- gcheckbox("plotNA")
  impVis.parcoordMiss.weighted <- gcheckbox("weighted")
  impVis.parcoordMiss.group[1,1] <- glabel("selection: ")
  impVis.parcoordMiss.group[1,2] <- impVis.parcoordMiss.selection
  impVis.parcoordMiss.group[2,1] <- impVis.parcoordMiss.plotNA
  impVis.parcoordMiss.group[2,2] <- impVis.parcoordMiss.weighted
  impVis.parcoordMiss.group[3,1:2, fill="y", expand="TRUE"] <- impVis.parcoordMiss.plotvars
  impVis.parcoordMiss.group[4,1:2, fill="y", expand="TRUE"] <- impVis.parcoordMiss.highlight
  addHandlerChanged(impVis.parcoordMiss.selection, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.parcoordMiss.plotNA, function(h,...){makeImputationPlot()})
  handlerList[[length(handlerList)+1]] <- c(impVis.parcoordMiss.highlight,
                                            addHandlerClicked(impVis.parcoordMiss.highlight, function(h,...){makeImputationPlot()}))
  handlerList[[length(handlerList)+1]] <- c(impVis.parcoordMiss.plotvars,
                                            addHandlerClicked(impVis.parcoordMiss.plotvars, function(h,...){makeImputationPlot()}))
  
  
  #pbox
  impVis.pbox.pos <- gdroplist("")
  impVis.pbox.selection <- gdroplist(c("any","all","none"))
  size(impVis.pbox.selection) <- c(140,-1)
  impVis.pbox.numbers <- gcheckbox("numbers")
  impVis.pbox.weighted<- gcheckbox("weighted")
  impVis.pbox.group[1,1] <- glabel("pos: ")
  impVis.pbox.group[1,2] <- impVis.pbox.pos
  impVis.pbox.group[2,1] <- glabel("selection: ")
  impVis.pbox.group[2,2] <- impVis.pbox.selection
  impVis.pbox.group[3,1] <- impVis.pbox.numbers
  impVis.pbox.group[3,2] <- impVis.pbox.weighted
  addHandlerChanged(impVis.pbox.selection, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.pbox.weighted, function(h,...){makeImputationPlot()})
  addHandlerChanged(impVis.pbox.numbers, function(h,...){makeImputationPlot()})
  handlerList[[length(handlerList)+1]] <- c(impVis.pbox.pos,
                                            addHandlerChanged(impVis.pbox.pos, function(h,...){makeImputationPlot()}))
  
  #matrixplot
  impVis.matrixplot.sortby <- gdroplist("")
  size(impVis.matrixplot.sortby) <- c(140,-1)
  impVis.matrixplot.weighted <- gcheckbox("weighted")
  impVis.matrixplot.group[1,1] <- glabel("sortby: ")
  impVis.matrixplot.group[1,2] <- impVis.matrixplot.sortby
  impVis.matrixplot.group[2,1] <- impVis.matrixplot.weighted
  handlerList[[length(handlerList)+1]] <- c(impVis.matrixplot.sortby,
                                            addHandlerChanged(impVis.matrixplot.sortby, function(h,...){makeImputationPlot()}))
  addHandlerChanged(impVis.matrixplot.weighted, function(h,...){makeImputationPlot()})
  
  
  #popup menu for saving images
	#activates by right click onto plot 
  sml <- list()
  sml$"Save as JPEG"$handler = function(h,...) savePlotToFile("JPEG")
  sml$"Save as PDF"$handler = function(h,...) savePlotToFile("PDF")
  sml$"Save as PNG"$handler = function(h,...) savePlotToFile("PNG")
  sml$"Save as PS"$handler = function(h,...) savePlotToFile("PS")
  sml$"Save as SVG"$handler = function(h,...) savePlotToFile("SVG")
  add3rdMousePopupmenu(impVis.plot, sml)
  
  #add rotated labels to the imputation plot notebook
	#this is done by using native GTK-methods
	#obtained by using the RGtk2 package
  notebook <- getToolkitWidget(impVis.plotBook)
  label <- gtkLabel("aggr")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 0), label)
  label <- gtkLabel("barMiss")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 1), label)
  label <- gtkLabel("histMiss")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 2), label)
  label <- gtkLabel("marginMatrix")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 3), label)
  label <- gtkLabel("scattMatrix")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 4), label)
  label <- gtkLabel("mosaicMiss")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 5), label)
  label <- gtkLabel("parcoordMiss")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 6), label)
  label <- gtkLabel("pbox")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 7), label)
  label <- gtkLabel("matrixplot")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 8), label)
  svalue(impVis.plotBook) <- 1
  addHandlerChanged(impVis.plotBook, imputationPlotHandler)
  
  #add rotated labels to the imputation notebook
	#this is done by using native GTK-methods
	#obtained by using the RGtk2 package
  notebook <- getToolkitWidget(imputation.notebook)
  label <- gtkLabel("knn")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 0), label)
  label <- gtkLabel("irmi")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 1), label)
  label <- gtkLabel("hotdeck")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 2), label)
  label <- gtkLabel("regression")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 3), label)
  label <- gtkLabel("undo imputation")
  gtkLabelSetAngle(label, 90)
  gtkNotebookSetTabLabel(notebook, gtkNotebookGetNthPage(notebook, 4), label)
  svalue(imputation.notebook) <- 1
  
  #clean up after closing the main window
	#in case of closing the window while creating a new plot
	#tries to delete a possible remaining temporary image in the working directory
  addHandlerUnrealize(mainWindow, handler=function(h,...){
    try({
      if (file.exists("current.tmp")){
        file.remove("current.tmp")
      }
    }, silent=TRUE)
  })
  
  #create a artificial dataset without real values in case the application
	#was started with a dataset
	#this has the purpose of allowing the widgets to initialize
	#but do this much faster than with a real dataset
  if (is.null(startupObject)){
    startupObject <- data.frame(empty="")
    setActiveDataset(startupObject, firstPage=FALSE, adjustTypes=FALSE)
  }
  else{
		#parse the name of the dataset used to start the application
		#used for the script window
    cmdimp <- paste("activedataset <- ",deparse(substitute(startupObject)))
    setActiveDataset(startupObject, firstPage=FALSE, loadScript=cmdimp)
  }
  
	#init all the tabs in the interface
  svalue(mainNotebook) <- 3
  #svalue(mainNotebook) <- 4
  svalue(mainNotebook) <- 1
}
