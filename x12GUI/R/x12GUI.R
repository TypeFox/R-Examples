#' Graphical User Interface for the S4 implementation of the X12-Arima wrapper
#' in package x12
#' 
#' GUI
#' 
#' @param x12orig object of class x12Batch or x12Single
#' @param \dots further arguments (currently ignored).
#' @author Daniel Schopfhauser
#' @seealso \code{\link{x12}}, \code{\link{x12env}},
#' \code{\linkS4class{x12Single}}, \code{\linkS4class{x12Batch}},
#' \code{\linkS4class{x12Parameter}}, \code{\linkS4class{x12List}},
#' \code{\linkS4class{x12Output}}, \code{\linkS4class{x12BaseInfo}},
#' \code{\link{summary}}, %\code{\link{summary.x12}}, \code{\link{x12work}}
#' @examples
#' 
#' \dontrun{
#' x12path("../x12a.exe")
#' data(AirPassengersX12Batch)
#' xbn <- x12GUI(AirPassengersX12Batch)
#' 
#' ##Create new x12Single and x12Batch objects
#' s1 <- new("x12Single",ts=AirPassengers,tsName="air")
#' s1 <- setP(s1,list(slidingspans=TRUE,history=TRUE,
#'   history.estimates=c("sadj","sadjchng","trend","trendchng","aic"),
#'   history.sadjlags=c(1,12),automdl=TRUE))
#' s2 <- new("x12Single",ts=ldeaths,tsName="ldeaths")
#' s3 <- new("x12Single",ts=UKgas,tsName="UKgas")
#' b <- new("x12Batch",list(s1,s2,s3))
#' ##Use GUI to handle single object
#' s1 <- x12GUI(s1)
#' ##Use GUI to handle batch object
#' b <- x12GUI(b)   
#' }
#' 
#' @import Hmisc
#' @import cairoDevice
#' @import grid
#' @import lattice
#' @import stringr
#' @import x12
#' @export x12GUI
x12GUI <- function(x12orig,...){
  tmpfsink <- tempfile()##tmpfile for catching the output to the console of x12
  window.main <- gtkWindow("x12 GUI")
  if(existd("x12path")==FALSE){
    dialog <- gtkMessageDialog(window.main, "destroy-with-parent",
                               "question", "yes-no", "The 'x12path' variable is not set, but is necessary for this application to work. Would you like to specify it?")
    ret <- dialog$Run()
    dialog$Destroy()
    if(ret==-8){
      dialog <- gtkFileChooserDialog("Open X12 Executable", window.main, "open",
                                     "gtk-cancel", GtkResponseType["cancel"], 
                                     "gtk-open", GtkResponseType["accept"])
      
      if (dialog$run() == GtkResponseType["accept"]) {
        x12path(dialog$getFilename())
      }
      
      dialog$destroy()
    }
  }
  
	if(class(x12orig)=="x12Single"){
		xl <- new("x12List")
		xl <- list(x12orig)
		xb <- new("x12Batch",list(x12orig@ts))
		xb@x12List[[1]] <- x12orig
		x12orig <- x12(xb)
	}else if(class(x12orig)=="ts"){
		x12orig <- x12(new("x12Batch",list(x12orig)))
	}else if(class(x12orig)=="x12Batch"){
    x12orig <- x12(x12orig)
  }
  
	########################################################
	##   Variables
	########################################################
	object <- x12orig
	indices <- c(1)
	
  regression <- rep(FALSE, length(object@x12List))
  
	locate <- FALSE
	clickedX <- 100
  
  zooming<-FALSE
	zoomX <- 0
  zoomY <- 0
  handler.moving <- 0
  pixbuf.plot <- 0
  
	context.plot1 <- 0
	context.plot2 <- 0
	context.plot3 <- 0
	context.plot4 <- 0
# context.plot5 <- 0
	
	notebook.plot <- gtkNotebook()
	statusbar <- gtkStatusbar()
	button.update <- gtkButton("Update")
	
	menubar.main <- gtkMenuBar()
	menuitem.x12update <- gtkMenuItemNewWithMnemonic("x12 _Update")
	menuitem.expplotaspdf <- gtkMenuItem("Export current Plot as PDF")
	menuitem.expplotaspng <- gtkMenuItem("Export current Plot as PNG")
	menuitem.expsummarycsv <- gtkMenuItem("Export summary as CSV")
	menuitem.expsummaryclipboard <- gtkMenuItem("Copy summary to clipboard")
	menuitem.x12loadp <- gtkMenuItem("load Parameter")
	menuitem.x12savep <- gtkMenuItem("save Parameter")
	menuitem.x12load <- gtkMenuItem("load")
	menuitem.x12save <- gtkMenuItem("save")
	menuitem.export <- gtkMenuItem("Export")
	menuitem.x12 <- gtkMenuItem("x12")
  menuitem.path <- gtkMenuItem(paste("x12Path: ",capPath(x12path())))
	menu.export <- gtkMenu()
	menu.x12 <- gtkMenu()
	
	#ts table
	table.ts <- gtkTreeView()
	renderer.ts <- gtkCellRendererText()
	column.ts <- gtkTreeViewColumn()
	table.model <- gtkListStore("character")
	handler.tstable <- 0
	
	frame.history <- gtkFrame("Archive")
	panel.history <- gtkVBox()
	label.history <- gtkLabel("revert to:")
	combobox.history <- gtkComboBoxNewText()
	button.revert <- gtkButton("Revert")
	button.clearhistory <- gtkButton("Clean Archive")
	count.history <- 0
	
	#Panels
	panel.main <- gtkHPaned()
	panel.window <- gtkVBox()
	#Left half of GUI
	panel.gui <- gtkHBox()
	panel.ts <- gtkVBox()
  panel.params <- gtkVBox()
	#panel.params <- gtkTable(60, 6)
	panel.scrolledparams <- gtkScrolledWindow()
	panel.plotp <- gtkVBox()
	
	#plotting frames und areas
	frame.plot1 <- gtkHBox()
	area.plot1 <- gtkDrawingArea()
	frame.plot2 <- gtkHBox()
	area.plot2 <- gtkDrawingArea()
	frame.plot3 <- gtkHBox()
	area.plot3 <- gtkDrawingArea()
	frame.plot4 <- gtkVBox()
	area.plot4 <- gtkDrawingArea()
# frame.plot5 <- gtkHBox()
# area.plot5 <- gtkDrawingArea()
	
	#plotcontextmenu
	menu.contextplotall <- gtkMenu()
	menuitem.saveaspdf <- gtkMenuItemNewWithLabel("Save as PDF")
	
	#plotcontextmenu with manual outlier
	menu.contextplotwithoutlier <- gtkMenu()
	menuitem.saveaspdfwithoutlier <- gtkMenuItemNewWithLabel("Save as PDF")
	menuitem.addAO <- gtkMenuItemNewWithLabel("Add AO Outlier (Regression must be active)")
	menuitem.addTC <- gtkMenuItemNewWithLabel("Add TC Outlier (Regression must be active)")
	menuitem.addLS <- gtkMenuItemNewWithLabel("Add LS Outlier (Regression must be active)")
	
	#plotslider
	slider.plotmin <- gtkHScale(min = 0, max= 100, step=5)
	slider.plotmax <- gtkHScale(min = 0, max= 100, step=5)
	
	#summary tab
	frame.summary <- gtkScrolledWindow()
	buffer.summary <- gtkTextBuffer()
	textview.summary <- gtkTextView()
	table.summary <- gtkTreeView()
	columns.summary <- 0
	model.summary <- gtkListStore(rep("character", length(object)+3))
	
	#summarytotal tab
	frame.summarytotal <- gtkScrolledWindow()
	buffer.summarytotal <- gtkTextBuffer()
	textview.summarytotal <- gtkTextView()
	
	#summary parameter
	frame.summaryparameter <- gtkFrame("Summary")
	panel.summaryparameter <- gtkVBox()
	checkb.fullSummary <- gtkCheckButtonNewWithLabel("show full summary")
	checkb.spectraldetail <- gtkCheckButtonNewWithLabel("spectral detail")
	checkb.almostout <- gtkCheckButtonNewWithLabel("almostout")
	checkb.rsdautocorr <- gtkCheckButtonNewWithLabel("rsd autocorr")
	checkb.quality.stat <- gtkCheckButtonNewWithLabel("quality stat")
	checkb.likelihoodstat <- gtkCheckButtonNewWithLabel("likelihood stat")
	checkb.aape <- gtkCheckButtonNewWithLabel("aape")
	checkb.idrsdseas <- gtkCheckButtonNewWithLabel("id rsdseas")
	checkb.summaryslidingspans <- gtkCheckButtonNewWithLabel("slidingspans")
	checkb.summaryhistory <- gtkCheckButtonNewWithLabel("history")
	checkb.summaryidentify <- gtkCheckButtonNewWithLabel("identify")
	
	#plotparameter panel plot
	panel.scrolledplotparams <- gtkScrolledWindow()
	frame.plotparams <- gtkFrame("Plot")
	panel.plotparams <- gtkVBox()
  checkb.orig <- gtkCheckButtonNewWithLabel("original")
  checkb.orig$SetActive(TRUE)
  checkb.sa <- gtkCheckButtonNewWithLabel("seasonally adjusted")
	checkb.trend <- gtkCheckButtonNewWithLabel("trend")
	checkb.logtransform <- gtkCheckButtonNewWithLabel("log-transformation")
	checkb.showAllout <- gtkCheckButtonNewWithLabel("show all outliers")
	frame.showout <- gtkFrame("show specific outlier")
	checkb.showout <- gtkCheckButton();
	panel.showout <- gtkTable(rows=3, columns=2)
	label.showoutyear <- gtkLabel("year:") 
	entry.showoutyear <- gtkEntry()
	label.showoutperiod <- gtkLabel("period:")
	entry.showoutperiod <- gtkEntry()
	#label.showouttype <- gtkLabel("type:")
	#combobox.showouttype <- gtkComboBoxNewText()
# checkb.showAlloutLines <- gtkCheckButtonNewWithLabel("Show Allout Lines")
# checkb.annComp <- gtkCheckButtonNewWithLabel("AnnComp")
# checkb.annCompTrend <- gtkCheckButtonNewWithLabel("AnnComp Trend")
	
	#plotparameter panel plotFbcast
	frame.plotFbcastparams <- gtkFrame("Plot Fore- & Backcast")
	panel.plotFbcastparams <- gtkVBox()
# checkb.forecast <- gtkCheckButtonNewWithLabel("Forecast")
# checkb.backcast <- gtkCheckButtonNewWithLabel("Backcast")
	checkb.showCI <- gtkCheckButtonNewWithLabel("show CI")
	checkb.logtransform_fb <- gtkCheckButtonNewWithLabel("Log-Transform")
	#checkb.showLine <- gtkCheckButtonNewWithLabel("Show Line")
	checkb.pointsOriginal <- gtkCheckButtonNewWithLabel("original points")
	
	#plotparameter panel plotSeasFac
	frame.plotSeasFacparams <- gtkFrame("Seasonal Factors")
	panel.plotSeasFacparams <- gtkVBox()
	checkb.SIratios <- gtkCheckButtonNewWithLabel("SI-Ratios")
	checkb.SIratiosreplaced <- gtkCheckButtonNewWithLabel("SI-Ratios replaced")
	
	#spectral plot
	frame.spectral <- gtkFrame("Spectral Plot")
	panel.spectral <- gtkVBox()
	radiob.spectralsa <- gtkRadioButtonNewWithLabel(label="sa")
	radiob.spectraloriginal <- gtkRadioButtonNewWithLabel("original", group=radiob.spectralsa$GetGroup())
	radiob.spectralirregular <- gtkRadioButtonNewWithLabel("irregular", group=radiob.spectralsa$GetGroup())
	radiob.spectralresiduals <- gtkRadioButtonNewWithLabel("residuals", group=radiob.spectralsa$GetGroup())
	#sa, original, irregular und residuals
	
	#autocorrelation plot
	frame.rsdacf <- gtkFrame("Autocorrelation Plot")
	panel.rsdacf <- gtkVBox()
	radiob.rsdacfacf <- gtkRadioButtonNewWithLabel(label="acf")
	radiob.rsdacfpacf <- gtkRadioButtonNewWithLabel("pacf", group=radiob.rsdacfacf$GetGroup())
	radiob.rsdacfacf2 <- gtkRadioButtonNewWithLabel("acf2", group=radiob.rsdacfacf$GetGroup())
	
	
	#manual outliers
	frame.manualoutlier <- gtkFrame("Manual Outliers")
	panel.manualoutlier <- gtkTable(rows=6, columns=2)
	scroll.manualoutlier <- gtkScrolledWindow()
	label.manualoutliertype <- gtkLabel("Type:")
	label.manualoutlieryear <- gtkLabel("Year:")
	label.manualoutlierperiod <- gtkLabel("Period:")
	table.manualoutlier <- gtkTreeView()
	renderer.manualoutliertype <- gtkCellRendererText()
	renderer.manualoutlieryear <- gtkCellRendererText()
	renderer.manualoutlierperiod <- gtkCellRendererText()
	column.manualoutliertype <- gtkTreeViewColumn()
	column.manualoutlieryear <- gtkTreeViewColumn()
	column.manualoutlierperiod <- gtkTreeViewColumn()
	tablemodel.manualoutlier <- gtkListStore("character", "character", "character")
	button.manualoutlierremove <- gtkButton("Remove")
	button.manualoutlieraddclick <- gtkToggleButton("Add by Click")
	button.manualoutlieradd <- gtkButton("Add")
	combobox.manualoutliertype <- gtkComboBoxNewText()
	entry.manualoutlieryear <- gtkEntry()
	entry.manualoutlierperiod <- gtkEntry()
	#List of OutlierLists for each ts, if no manual outliers the list element should be NA
	outlierlist <- list()
	
	#x12 parameters########################
	#span
	checkb.spanactive <- gtkCheckButtonNewWithLabel("span")
	frame.series <- gtkFrame("Series")
  #label.series.span <- gtkLabel("span")
	panel.series <- gtkTable(12, 6)
	checkb.spanstart <- gtkCheckButtonNewWithLabel("start")
	label.spanstartyear <- gtkLabel("year:")
	entry.spanstartyear <- gtkEntry()
	label.spanstartperiod <- gtkLabel("period:")
	entry.spanstartperiod <- gtkEntry()
	panel.spanstart <- gtkHBox()
	checkb.spanend <- gtkCheckButtonNewWithLabel("end")
	label.spanendyear <- gtkLabel("year:")
	entry.spanendyear <- gtkEntry()
	label.spanendperiod <- gtkLabel("period:")
	entry.spanendperiod <- gtkEntry()
	panel.spanend <- gtkHBox()
	handler.spanstartyear <- 1
	handler.spanstartperiod <- 1
	handler.spanendyear <- 1
	handler.spanendperiod <- 1
	handler.spanactive <- 1
	
	#modelspan
	checkb.modelspanactive <- gtkCheckButtonNewWithLabel("modelspan")
	frame.modelspan <- gtkFrame("Modelspan")
	panel.modelspan <- gtkVBox()
	checkb.modelspanstart <- gtkCheckButtonNewWithLabel("start")
	label.modelspanstartyear <- gtkLabel("year:")
	entry.modelspanstartyear <- gtkEntry()
	label.modelspanstartperiod <- gtkLabel("period:")
	entry.modelspanstartperiod <- gtkEntry()
	panel.modelspanstart <- gtkHBox()
	checkb.modelspanend <- gtkCheckButtonNewWithLabel("end")
	label.modelspanendyear <- gtkLabel("year:")
	entry.modelspanendyear <- gtkEntry()
	label.modelspanendperiod <- gtkLabel("period:")
	entry.modelspanendperiod <- gtkEntry()
	panel.modelspanend <- gtkHBox()
	handler.modelspanstartyear <- 1
	handler.modelspanstartperiod <- 1
	handler.modelspanendyear <- 1
	handler.modelspanendperiod <- 1
  
  #series.type
#   checkb.series.type <- gtkCheckButtonNewWithLabel("type")
#   combobox.series.type <- gtkComboBoxNewText()
#   label.series.type <- gtkLabel("type:") 
#   handler.series.type <- 1
	
	
#	#decimals
## checkb.decimalsactive <- gtkCheckButton()
##	panel.decimals <- gtkHBox()
#	label.decimals <- gtkLabel("Decimals:")
#	entry.decimals <- gtkEntry()
#	handler.decimals <- 1
	
	#transform
  frame.transform <- gtkFrame("Transform")
  panel.transform <- gtkTable(2,3)
  #transform.function
	label.transform <- gtkLabel("function:")
	combobox.transform <- gtkComboBoxNewText()
	handler.transform <- 1
  #transform.power
  checkb.transform.power <- gtkCheckButtonNewWithLabel("power:")
  entry.transform.power <- gtkEntry()
  handler.transform.power <- 1
  #transform.adjust
	checkb.transform.adjust <- gtkCheckButtonNewWithLabel("adjust:")
	combobox.transform.adjust <- gtkComboBoxNewText()
	handler.transform.adjust <- 1
  
#  ######Seats-Frame
#  frame.seats <- gtkFrame("Seats")
#  panel.seats <- gtkTable(3,3)
#	checkb.seats <- gtkCheckButton()
#	handlercheckb.seats <- 1
#	#seatsparameter
#	checkb.seatsparameteractive <- gtkCheckButtonNewWithLabel("seatsparameter")
#	entry.seatsparameter <- gtkEntry()
#	handler.seatsparameter <- 1
  
  
	## Arima
	frame.arima <- gtkFrame("Arima")
	panel.arima <- gtkTable(5, 6)	
	#arima
	checkb.arimaactive <- gtkCheckButtonNewWithLabel("model:")
	entry.arima1 <- gtkEntry()
	entry.arima2 <- gtkEntry()
	entry.arima3 <- gtkEntry()
	handler.arima1 <- 1
	handler.arima2 <- 1
	handler.arima3 <- 1	
	#sarima
	checkb.sarimaactive <- gtkCheckButtonNewWithLabel("smodel:")
	entry.sarima1 <- gtkEntry()
	entry.sarima2 <- gtkEntry()
	entry.sarima3 <- gtkEntry()
	handler.sarima1 <- 1
	handler.sarima2 <- 1
	handler.sarima3 <- 1
  #as
  checkb.arima.ar <- gtkCheckButtonNewWithLabel("ar:")
  entry.arima.ar <- gtkEntry()
	#ma
	checkb.arima.ma <- gtkCheckButtonNewWithLabel("ma:")
	entry.arima.ma <- gtkEntry()
  
  ########Automdl-Frame
	checkb.automdl <- gtkCheckButton()
  panel.automdl <- gtkTable(3,5)
	frame.automdl <- gtkFrame("Automdl")
	checkb.acceptdefault <- gtkCheckButtonNewWithLabel("Acceptdefault")
	checkb.balanced <- gtkCheckButtonNewWithLabel("Balanced")
	handlercheckb.acceptdefault <- 1
	handlercheckb.balanced <- 1
	handlercheckb.automdl <- 1
	#maxorder
	#checkb.maxorderactive <- gtkCheckButtonNewWithLabel("maxorder:")
  label.maxorder <- gtkLabel("maxorder:")
	entry.maxorder1 <- gtkEntry()
  entry.maxorder2 <- gtkEntry()
	handler.maxorder1 <- 1
  handler.maxorder2 <- 1
  #maxdiff
	#checkb.maxdiffactive <- gtkCheckButtonNewWithLabel("maxdiff:")
  label.maxdiff <- gtkLabel("maxdiff:")
	entry.maxdiff1 <- gtkEntry()
	handler.maxdiff1 <- 1
  entry.maxdiff2 <- gtkEntry()
  handler.maxdiff2 <- 1
  
  #########Estimate-Frame
  panel.estimate <- gtkTable(2,2)
  frame.estimate <- gtkFrame("Estimate")
	checkb.estimate <- gtkCheckButton()
	checkb.estOutofsample <- gtkCheckButtonNewWithLabel("outofsample")
	handlercheckb.estimate <- 1
	handlercheckb.estOutofsample <- 1
  
	#########Check-Frame
	panel.check <- gtkTable(2,2)
	frame.check <- gtkFrame("Check")
	checkb.check <- gtkCheckButton()
	checkb.check.maxlag <- gtkCheckButtonNewWithLabel("maxlag")
  entry.check.maxlag <- gtkEntry()
	handlercheckb.check <- 1
	handlercheckb.check.maxlag <- 1
  handler.check.maxlag <- 1
  
  #########Slidingspans-Frame
  panel.slidingspans <- gtkTable(6, 8)
  frame.slidingspans <- gtkFrame("Slidingspans")
	checkb.slidingspans <- gtkCheckButton()
	handlercheckb.slidingspans <- 1
  combobox.slidingspans.fixmdl <- gtkComboBoxNewText()
	handler.slidingspans.fixmdl <- 1
  checkb.slidingspans.fixmdl <- gtkCheckButtonNewWithLabel("fixmdl:")
	handlercheckb.slidingspans.fixmdl <- 1
	checkb.slidingspans.fixreg <- gtkCheckButtonNewWithLabel("fixreg:")
	handlercheckb.slidingspans.fixreg <- 1
	checkb.slidingspans.fixreg1 <- gtkCheckButtonNewWithLabel("td")
	checkb.slidingspans.fixreg2 <- gtkCheckButtonNewWithLabel("holiday")
	checkb.slidingspans.fixreg3 <- gtkCheckButtonNewWithLabel("outlier")
	checkb.slidingspans.fixreg4 <- gtkCheckButtonNewWithLabel("user")
	handler.slidingspans.fixreg1 <- 1
	handler.slidingspans.fixreg2 <- 1
	handler.slidingspans.fixreg3 <- 1
	handler.slidingspans.fixreg4 <- 1
	checkb.slidingspans.length <- gtkCheckButtonNewWithLabel("length:")
  entry.slidingspans.length <- gtkEntry()
	handler.slidingspans.length <- 1
	handlercheckb.slidingspans.length <- 1
	checkb.slidingspans.numspans <- gtkCheckButtonNewWithLabel("numspans:")
	entry.slidingspans.numspans <- gtkEntry()
	handler.slidingspans.numspans <- 1
	handlercheckb.slidingspans.numspans <- 1
	checkb.slidingspans.outlier <- gtkCheckButtonNewWithLabel("outlier:")
	combobox.slidingspans.outlier <- gtkComboBoxNewText()
	handler.slidingspans.outlier <- 1
	handlercheckb.slidingspans.numspans <- 1
	checkb.slidingspans.additivesa <- gtkCheckButtonNewWithLabel("additivesa:")
	combobox.slidingspans.additivesa <- gtkComboBoxNewText()
	handler.slidingspans.additivesa <- 1
	handlercheckb.slidingspans.numspans <- 1
	checkb.slidingspans.start <- gtkCheckButtonNewWithLabel("start:")
	entry.slidingspans.start1 <- gtkEntry()
	entry.slidingspans.start2 <- gtkEntry()
	label.slidingspans.start1 <- gtkLabel("year:")
	label.slidingspans.start2 <- gtkLabel("period:")
	handlerhandlercheck.slidingspans.numspans <- 1
	handler.slidingspans.start1 <- 1
	handler.slidingspans.start2 <- 1
  
	#div checkboxes
	frame.divboxes <- gtkFrame("Settings")
	panel.divboxes <- gtkVBox()
	#handlercheckb.onlytd <- 1
	
  ####Regression-Frame
#   checkb.x11regress <- gtkCheckButtonNewWithLabel("x11regress")
# 	handlercheckb.x11regress <- 1
  
  panel.regression <- gtkTable(10, 3)
  frame.regression <- gtkFrame("Regression")
  
  checkb.regressionactive <- gtkCheckButton()
  handler.regressionactive <- 1
  radiob.regression <- gtkRadioButtonNewWithLabel(label="Regression")
	radiob.x11regression <- gtkRadioButtonNewWithLabel("X11 Regression", group=radiob.regression$GetGroup())
  handler.x11regression <- 1
  
	#regvariables
	checkb.regvariablesactive <- gtkCheckButtonNewWithLabel("variables:")
#	panel.regvariables <- gtkHBox()
	#label.regvariables <- gtkLabel("Regvariables:")
	entry.regvariables <- gtkEntry()
	handler.regvariables <- 1
	
	#reguser
	checkb.reguseractive <- gtkCheckButtonNewWithLabel("user:")
#	panel.reguser <- gtkHBox()
	#label.reguser <- gtkLabel("Reguser:")
	entry.reguser <- gtkEntry()
	handler.reguser <- 1
	
	#regfile
	checkb.regfileactive <- gtkCheckButtonNewWithLabel("file:")
#	panel.regfile <- gtkHBox()
	#label.regfile <- gtkLabel("Regfile:")
	filebutton.regfile <- gtkFileChooserButton("RegFile","GTK_FILE_CHOOSER_ACTION_OPEN")
	
	#usertype
	checkb.usertypeactive <- gtkCheckButtonNewWithLabel("usertype:")
#	panel.usertype <- gtkHBox()
	#label.usertype <- gtkLabel("usertype:")
	entry.usertype <- gtkEntry()
	handler.usertype <- 1
	
	#centeruser
	checkb.centeruseractive <- gtkCheckButtonNewWithLabel("centeruser:")
#	panel.centeruser <- gtkHBox()
	#label.centeruser <- gtkLabel("centeruser:")
	combobox.centeruser <- gtkComboBoxNewText()
	handler.centeruser <- 1
	
	#regfilestart
	checkb.regfilestartactive <- gtkCheckButtonNewWithLabel("start:")
	#frame.regfilestart <- gtkFrame("RegFileStart")
	panel.regfilestart <- gtkVBox()
	label.regfilestartstartyear <- gtkLabel("Year:")
	entry.regfilestartstartyear <- gtkEntry()
	label.regfilestartstartperiod <- gtkLabel("Period:")
	entry.regfilestartstartperiod <- gtkEntry()
	panel.regfilestartstart <- gtkHBox()
	handler.regfilestartstartyear <- 1
	handler.regfilestartstartperiod <- 1
  
	#aictest
	checkb.aictestactive <- gtkCheckButtonNewWithLabel("aictest:")
	#	panel.aictest <- gtkHBox()
	#label.aictest <- gtkLabel("aictest:")
	entry.aictest <- gtkEntry()
	handler.aictest <- 1
	
	#########Outlier-Frame
  frame.outlier <- gtkFrame("Outlier")
  panel.outlier <- gtkTable(5, 12)
	#outlier
	checkb.outlieractive <- gtkCheckButton()
	label.outlier <- gtkLabel("types:")
	checkb.outlierall <- gtkCheckButtonNewWithLabel("all")
	checkb.outlierTC <- gtkCheckButtonNewWithLabel("TC")
	checkb.outlierAO <- gtkCheckButtonNewWithLabel("AO")
	checkb.outlierLS <- gtkCheckButtonNewWithLabel("LS")
	handler.outlier <- 1
	handler.outlierall <- 1
	handler.outlierAO <- 1
	handler.outlierLS <- 1
	handler.outlierTC <- 1
	#critical
	checkb.criticalactive <- gtkCheckButtonNewWithLabel("critical:")
	radiob.criticalall <- gtkRadioButtonNewWithLabel(label="all")
	radiob.criticalspecific <- gtkRadioButtonNewWithLabel(group = radiob.criticalall$GetGroup(), label="specific")
	entry.criticalall <- gtkEntry()
	entry.criticalAO <- gtkEntry()
	entry.criticalLS <- gtkEntry()
	entry.criticalTC <- gtkEntry()
	label.criticalAO <- gtkLabel("AO")
	label.criticalLS <- gtkLabel("LS")
	label.criticalTC <- gtkLabel("TC")
	handler.criticalall <- 1
	handler.criticalAO <- 1
	handler.criticalTC <- 1
	handler.criticalLS <- 1
	#outlierspan
	checkb.outlierspanactive <- gtkCheckButtonNewWithLabel("span:")
  checkb.outlierspan.start <- gtkCheckButtonNewWithLabel("start:")
	checkb.outlierspan.end <- gtkCheckButtonNewWithLabel("end:")
  label.outlierspan.year <- gtkLabel("year:")
	label.outlierspan.period <- gtkLabel("period:")
	entry.outlierspan.start1 <- gtkEntry()
	entry.outlierspan.start2 <- gtkEntry()
	handler.outlierspan.start1 <- 1
	handler.outlierspan.start2 <- 1
	entry.outlierspan.end1 <- gtkEntry()
	entry.outlierspan.end2 <- gtkEntry()
	handler.outlierspan.end1 <- 1
	handler.outlierspan.end2 <- 1
	#outliermethod
	checkb.outliermethodactive <- gtkCheckButtonNewWithLabel("method:")
	combobox.outliermethod <- gtkComboBoxNewText()
	handler.outliermethod <- 1
  
# #tblnames 
# checkb.tblnamesactive <- gtkCheckButton()
# panel.tblnames <- gtkHBox()
# label.tblnames <- gtkLabel("tblnames:")
# entry.tblnames <- gtkEntry()
# 
# #Rtblnames 
# checkb.Rtblnamesactive <- gtkCheckButton()
# panel.Rtblnames <- gtkHBox()
# label.Rtblnames <- gtkLabel("Rtblnames:")
# entry.Rtblnames <- gtkEntry()
# 
# #x12path
# checkb.x12pathactive <- gtkCheckButton()
# panel.x12path <- gtkHBox()
# label.x12path <- gtkLabel("x12path:")
# filebutton.x12path <- gtkFileChooserButton("x12path","GTK_FILE_CHOOSER_ACTION_OPEN")
# 
# #x13path
# checkb.x13pathactive <- gtkCheckButton()
# panel.x13path <- gtkHBox()
# label.x13path <- gtkLabel("x13path:")
# filebutton.x13path <- gtkFileChooserButton("x13path","GTK_FILE_CHOOSER_ACTION_OPEN")
# 
# #use
# checkb.useactive <- gtkCheckButton()
# panel.use <- gtkHBox()
# label.use <- gtkLabel("use:")
# combobox.use <- gtkComboBoxNewText()
  
  #######Identify-Frame
  checkb.identify <- gtkCheckButton()
  panel.identify <- gtkTable(3,4)
  frame.identify <- gtkFrame("Identify")
  checkb.identify.diff <- gtkCheckButtonNewWithLabel("diff:")
	checkb.identify.sdiff <- gtkCheckButtonNewWithLabel("sdiff:")
	checkb.identify.maxlag <- gtkCheckButtonNewWithLabel("maxlag:")
  entry.identify.diff <- gtkEntry()
  entry.identify.sdiff <- gtkEntry()
  entry.identify.maxlag <- gtkEntry()
  handler.identify.sdiff <- 1 
	handler.identify.diff <- 1
	handler.identify.maxlag <- 1
  
# #file
# checkb.fileactive <- gtkCheckButton()
# panel.file <- gtkHBox()
# label.file <- gtkLabel("file:")
# entry.file <- gtkEntry()
	
	########History-Frame
	frame.historyparam <- gtkFrame("History")
	panel.historyparam <- gtkTable()
  checkb.historyactive <- gtkCheckButton()
  checkb.history.estimates <- gtkCheckButtonNewWithLabel("estimates:")
	checkb.history.estimatessadj <- gtkCheckButtonNewWithLabel("sadj")
	checkb.history.estimatessadjchng <- gtkCheckButtonNewWithLabel("sadjchng")
	checkb.history.estimatestrend <- gtkCheckButtonNewWithLabel("trend")
	checkb.history.estimatestrendchng <- gtkCheckButtonNewWithLabel("trendchng")
	checkb.history.estimatesseasonal <- gtkCheckButtonNewWithLabel("seasonal")
	checkb.history.estimatesfcst <- gtkCheckButtonNewWithLabel("fcst")
	checkb.history.estimatesaic <- gtkCheckButtonNewWithLabel("aic")
	handler.history.estimatessadj <- 1
	handler.history.estimatessadjchng <- 1
	handler.history.estimatestrend <- 1
	handler.history.estimatestrendchng <- 1
	handler.history.estimatesseasonal <- 1
	handler.history.estimatesfcst <- 1
	handler.history.estimatesaic <- 1
  handler.history <- 1
	checkb.history.fixmdl <- gtkCheckButtonNewWithLabel("fixmdl")
	checkb.history.fixreg <- gtkCheckButtonNewWithLabel("fixreg:")
	checkb.history.fixreg1 <- gtkCheckButtonNewWithLabel("td")
	checkb.history.fixreg2 <- gtkCheckButtonNewWithLabel("holiday")
	checkb.history.fixreg3 <- gtkCheckButtonNewWithLabel("outlier")
	checkb.history.fixreg4 <- gtkCheckButtonNewWithLabel("user")
	handler.history.fixreg1 <- 1
	handler.history.fixreg2 <- 1
	handler.history.fixreg3 <- 1
	handler.history.fixreg4 <- 1
  #outlier
  checkb.history.outlier <- gtkCheckButtonNewWithLabel("outlier:")
  combobox.history.outlier <- gtkComboBoxNewText()
  handler.history.outlier <- 1
  #target
	checkb.history.target <- gtkCheckButtonNewWithLabel("target:")
	combobox.history.target <- gtkComboBoxNewText()
	handler.history.target <- 1
  #sadjlags
  checkb.history.sadjlags <- gtkCheckButtonNewWithLabel("sadjlags:")
  entry.history.sadjlags <- gtkEntry()
  handler.history.sadjlags <- 1
	#trendlags
	checkb.history.trendlags <- gtkCheckButtonNewWithLabel("trendlags:")
	entry.history.trendlags <- gtkEntry()
	handler.history.trendlags <- 1
  #start
  checkb.history.start <- gtkCheckButtonNewWithLabel("start:")
  label.history.startyear <- gtkLabel("year:")
  label.history.startperiod <- gtkLabel("period:")
  entry.history.startyear <- gtkEntry()
  entry.history.startperiod <- gtkEntry()
	handler.history.startyear <- 1
	handler.history.startperiod <- 1
  
  ######Forecast-Frame
  frame.forecast <- gtkFrame("Forecast")
  panel.forecast <- gtkTable(2,3)
	#forecast_years
	checkb.forecast_yearsactive <- gtkCheckButtonNewWithLabel("forecast_years")
	entry.forecast_years <- gtkEntry()
	handler.forecast_years <- 1
	#backcast_years
	checkb.backcast_yearsactive <- gtkCheckButtonNewWithLabel("backcast_years")
	entry.backcast_years <- gtkEntry()
	handler.backcast_years <- 1
	#forecast_conf
  checkb.forecast_confactive <- gtkCheckButton("forecast_conf")
	entry.forecast_conf <- gtkEntry()
	handler.forecast_conf <- 1
  
  ######X11-Frame
  frame.x11 <- gtkFrame("X11")
  panel.x11 <- gtkTable(3,11)
	#sigmalim
	checkb.sigmalimactive <- gtkCheckButtonNewWithLabel("sigmalim:")
	entry.sigmalim1 <- gtkEntry()
	entry.sigmalim2 <- gtkEntry()
	handler.sigmalim1 <- 1
	handler.sigmalim2 <- 1
  #type
  checkb.x11.type <- gtkCheckButtonNewWithLabel("type:")
  combobox.x11.type <- gtkComboBoxNewText()
  handler.x11type <- 1
  #div
	checkb.sfshort <- gtkCheckButtonNewWithLabel("sfshort")
	checkb.x11appendfcst <- gtkCheckButtonNewWithLabel("x11appendfcst")
	checkb.x11appendbcst <- gtkCheckButtonNewWithLabel("x11appendbcst")
	checkb.x11excludefcst <- gtkCheckButtonNewWithLabel("x11excludefcst")
	handlercheckb.sfshort <- 1
	handlercheckb.x11appendfcst <- 1
	handlercheckb.x11appendbcst <- 1
	handlercheckb.x11excludefcst <- 1
	#samode
	checkb.samodeactive <- gtkCheckButtonNewWithLabel("samode:")
	combobox.samode <- gtkComboBoxNewText()
	handler.samode <- 1
	#seasonalma
	checkb.seasonalmaactive <- gtkCheckButtonNewWithLabel("seasonalma:")
	entry.seasonalma <- gtkEntry()
	handler.seasonalma <- 1
	#trendma
	checkb.trendmaactive <- gtkCheckButtonNewWithLabel("trendma:")
	entry.trendma <- gtkEntry()
	handler.trendma <- 1
	#x11calendarsigma
	checkb.x11calendarsigmaactive <- gtkCheckButtonNewWithLabel("x11calendarsigma:")
	combobox.x11calendarsigma <- gtkComboBoxNewText()
	handler.x11calendarsigma <- 1
	#x11final
	checkb.x11.finalactive <- gtkCheckButtonNewWithLabel("final:")
	checkb.x11.finalAO <- gtkCheckButtonNewWithLabel("AO")
  handler.x11.finalAO <- 1
	checkb.x11.finalLS <- gtkCheckButtonNewWithLabel("LS")
	handler.x11.finalLS <- 1
	checkb.x11.finalTC <- gtkCheckButtonNewWithLabel("TC")
	handler.x11.finalTC <- 1
	checkb.x11.finaluser <- gtkCheckButtonNewWithLabel("user")
	handler.x11.finaluser <- 1
	checkb.x11.finalnone <- gtkCheckButtonNewWithLabel("none")
	handler.x11.finalnone <- 1
  
	###DELETE ME SOMEDAY
  entry.x11final <- gtkEntry()
	handler.x11final <- 1
  ###
  
	
	####x12TOOLTIPS
	string.span <- "text <i>kursiv</i><b>fett</b>"
	string.modelspan <- "text <i>kursiv</i><b>fett</b>"
	string.decimals <- "text <i>kursiv</i><b>fett</b>"
	string.transform <- "text <i>kursiv</i><b>fett</b>"
	string.arima <- "text <i>kursiv</i><b>fett</b>"
	string.sarima <- "text <i>kursiv</i><b>fett</b>"
	string.automdl <- "text <i>kursiv</i><b>fett</b>"
	string.acceptdefault <- "text <i>kursiv</i><b>fett</b>"
	string.balanced <- "text <i>kursiv</i><b>fett</b>"
#	string.seats <- "text <i>kursiv</i><b>fett</b>"
	string.estimate <- "text <i>kursiv</i><b>fett</b>"
	string.estimateoutofsamples <- "text <i>kursiv</i><b>fett</b>"
	string.slidingspans <- "text <i>kursiv</i><b>fett</b>"
	string.onlytd <- "text <i>kursiv</i><b>fett</b>"
	string.sfshort <- "text <i>kursiv</i><b>fett</b>"
	string.x11appendfcst <- "text <i>kursiv</i><b>fett</b>"
	string.x11appendfbst <- "text <i>kursiv</i><b>fett</b>"
	string.x11excludefcst <- "text <i>kursiv</i><b>fett</b>"
	string.x11regress <- "text <i>kursiv</i><b>fett</b>"
	string.maxorder <- "text <i>kursiv</i><b>fett</b>"
	string.maxdiff <- "text <i>kursiv</i><b>fett</b>"
	string.regvariables <- "text <i>kursiv</i><b>fett</b>"
	string.reguser <- "text <i>kursiv</i><b>fett</b>"
	string.regfile <- "text <i>kursiv</i><b>fett</b>"
	string.usertype <- "text <i>kursiv</i><b>fett</b>"
	string.centeruser <- "text <i>kursiv</i><b>fett</b>"
	string.regfilestart <- "text <i>kursiv</i><b>fett</b>"
#	string.seatsparameter <- "text <i>kursiv</i><b>fett</b>"
	string.sigmalim <- "text <i>kursiv</i><b>fett</b>"
	string.critical <- "text <i>kursiv</i><b>fett</b>"
	string.outlier <- "text <i>kursiv</i><b>fett</b>"
	string.outlierspan <- "text <i>kursiv</i><b>fett</b>"
	string.outliermethod <- "text <i>kursiv</i><b>fett</b>"
	string.forecast_years <- "text <i>kursiv</i><b>fett</b>"
	string.backcast_years <- "text <i>kursiv</i><b>fett</b>"
	string.forecast_conf <- "text <i>kursiv</i><b>fett</b>"
	string.aictest <- "text <i>kursiv</i><b>fett</b>"
	string.samode <- "text <i>kursiv</i><b>fett</b>"
	string.seasonalma <- "text <i>kursiv</i><b>fett</b>"
	string.trendma <- "text <i>kursiv</i><b>fett</b>"
	string.x11calendarsigma <- "text <i>kursiv</i><b>fett</b>"
	string.x11final <- "text <i>kursiv</i><b>fett</b>"
	
	####END Variables###################################
	
	####################################################
	#FUNCTIONS
	####################################################
	#Handler function for the notepad with the plots
	notebookhandler <- function(notebook, page, page.num, user.data){
		update_notebook(page.num)
	}
	
	#Handler function for the NULL-checkboxes of the x12 panel 
	checkbuttonhandler <- function(togglebutton, user_data){
		update_toggle(button=togglebutton, system=FALSE)
	}
	
	#Handler function for the checkbuttons for the boolean x12 values
	checkbuttonx12handler <- function(togglebutton, user_data){
		#automdl
		if(togglebutton == checkb.automdl){
		  toggle(c(checkb.acceptdefault, checkb.balanced,
               entry.maxorder1, entry.maxdiff1,entry.maxorder2, entry.maxdiff2), togglebutton)
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(automdl=TRUE),s)})
			else{
			  lapply(indices, FUN=function(s){object <<- setP(object,list(automdl=FALSE,
                                                                    automdl.acceptdefault=FALSE,
			                                                              automdl.balanced=TRUE,
			                                                              automdl.maxorder=c(3,2),
			                                                              automdl.maxdiff=c(1,1)),s)})
        checkb.acceptdefault$SetActive(FALSE)
        checkb.balanced$SetActive(FALSE)
        entry.maxorder1$SetText("3")
			  entry.maxdiff1$SetText("1")
			  entry.maxorder2$SetText("2")
			  entry.maxdiff2$SetText("1")
			} 
			status_print("automdl changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#regression
		if(togglebutton == checkb.regressionactive){
# 		  toggle(c(checkb.acceptdefault, checkb.balanced, checkb.maxorderactive, checkb.maxdiffactive,
# 		           entry.maxorder1, entry.maxdiff1), togglebutton)
		  if(togglebutton$GetActive()==TRUE) {
        status_print("regressionactive changed!")
        regression[indices] <<- TRUE
        toggle(c(radiob.regression, radiob.x11regression,
                 checkb.regvariablesactive, 
                 checkb.reguseractive, 
                 checkb.regfileactive, 
                 checkb.usertypeactive, 
                 checkb.centeruseractive, 
                 checkb.regfilestartactive, 
                 checkb.aictestactive), checkb.regressionactive)
        #removes the possibility to add manual outlier if the regression is not active
        #would lead to some questionable situations
        toggle(c(menuitem.addAO, menuitem.addLS, menuitem.addTC, 
                 button.manualoutlieraddclick, button.manualoutlieradd), checkb.regressionactive)
		  }
		  else{
		    ####Regvariables deletion?
		    capture.output(v <- getP(object, list("regression.variables")))
        remove <- TRUE
		    v <- v[indices]
		    t <- sapply(lapply(v, FUN=function(s){s$regression.variables}), FUN=function(s){containsNULL(list(s))})
		    if(FALSE %in% t){
		      dialog <- gtkMessageDialog(window.main, "destroy-with-parent",
		                                 "warning", "yes-no", "Removing regression will also remove specified regression.variables and manually defined outliers. Proceed?")
		      if(dialog$run()==GtkResponseType["yes"]){
		        remove <- TRUE
		      }
		      else{
		        checkb.regressionactive$SetActive(TRUE)
            remove <- FALSE
		      }
		      dialog$destroy()
		    }  
        if(remove==TRUE){
          regression[indices] <<- FALSE
          lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=NULL,
                                                                      regression.user=NULL,
                                                                      regression.file=NULL,
                                                                      regression.usertype=NULL,
                                                                      regression.centeruser=NULL,
                                                                      regression.start=NULL,
                                                                      regression.aictest=NULL,
                                                                      x11regression=FALSE),s)})
          toggle(c(radiob.regression, radiob.x11regression,
                   checkb.regvariablesactive, entry.regvariables,
                   checkb.reguseractive, entry.reguser,
                   checkb.regfileactive, filebutton.regfile,
                   checkb.usertypeactive, entry.usertype,
                   checkb.centeruseractive, combobox.centeruser,
                   checkb.regfilestartactive, label.regfilestartstartyear,
                   entry.regfilestartstartyear, label.regfilestartstartperiod,
                   entry.regfilestartstartperiod,
                   checkb.aictestactive, entry.aictest), checkb.regressionactive)
          #removes the possibility to add manual outlier if the regression is not active
          #would lead to some questionable situations
          toggle(c(menuitem.addAO, menuitem.addLS, menuitem.addTC, 
                   button.manualoutlieraddclick, button.manualoutlieradd), checkb.regressionactive)
          radiob.regression$SetActive(TRUE) 
          radiob.x11regression$SetActive(FALSE)
          checkb.regvariablesactive$SetActive(FALSE)
          entry.regvariables$SetText("")
          checkb.reguseractive$SetActive(FALSE)
          entry.reguser$SetText("")
          checkb.regfileactive$SetActive(FALSE) 
          filebutton.regfile$UnselectAll()
          checkb.usertypeactive$SetActive(FALSE) 
          entry.usertype$SetText("")
          checkb.centeruseractive$SetActive(FALSE) 
          combobox.centeruser$SetActive(-1)
          checkb.regfilestartactive$SetActive(FALSE)
          entry.regfilestartstartyear$SetText("") 
          entry.regfilestartstartperiod$SetText("")
          checkb.aictestactive$SetActive(FALSE) 
          entry.aictest$SetText("")
        }
		    
		  }
      ######
		    
		  status_print("regression changed!")
		  togglebutton$SetInconsistent(FALSE)
		}
    #identify
		if(togglebutton == checkb.identify){
		  if(togglebutton$GetActive()==TRUE){
        lapply(indices, FUN=function(s){object <<- setP(object,list(identify=TRUE),s)})
        toggle(c(checkb.identify.diff, checkb.identify.sdiff, checkb.identify.maxlag), togglebutton)
		  }
		  else{
		    toggle(c(checkb.identify.diff, checkb.identify.sdiff, checkb.identify.maxlag,
		             entry.identify.diff, entry.identify.sdiff, entry.identify.maxlag), togglebutton)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(identify=FALSE,
                                                                    identify.diff=NULL,
                                                                    identify.sdiff=NULL,
                                                                    identify.maxlag=NULL),s)})
		    checkb.identify.diff$SetActive(FALSE)
		    entry.identify.diff$SetText("")
		    checkb.identify.sdiff$SetActive(FALSE)
		    entry.identify.sdiff$SetText("")
		    checkb.identify.maxlag$SetActive(FALSE)
		    entry.identify.maxlag$SetText("")
		  } 
		  status_print("identify changed!")
		  togglebutton$SetInconsistent(FALSE)
		}
		#check
		if(togglebutton == checkb.check){
		  if(togglebutton$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(check=TRUE),s)})
		    toggle(c(checkb.check.maxlag), togglebutton)
		  }
		  else{
		    toggle(c(checkb.check.maxlag, entry.check.maxlag), togglebutton)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(check=FALSE,
		                                                                check.maxlag=NULL),s)})
		    checkb.check.maxlag$SetActive(FALSE)
		    entry.check.maxlag$SetText("")
		  } 
		  status_print("check changed!")
		  togglebutton$SetInconsistent(FALSE)
		}
		#history
		if(togglebutton == checkb.historyactive){
		  if(togglebutton$GetActive()==TRUE){
        lapply(indices, FUN=function(s){object <<- setP(object,list(history=TRUE),s)})
        toggle(c(checkb.history.estimates,
                 checkb.history.fixmdl,
                 checkb.history.fixreg,
                 checkb.history.outlier, 
				 checkb.history.target, 
                 checkb.history.sadjlags, 
                 checkb.history.trendlags, 
                 checkb.history.start), togglebutton)
		  }
		  else{
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history=FALSE,
		                                                                history.estimates=NULL,
		                                                                history.fixmdl=FALSE,
		                                                                history.fixreg=NULL,
		                                                                history.outlier=NULL,
																		history.target=NULL,
		                                                                history.sadjlags=NULL,
		                                                                history.trendlags=NULL,
		                                                                history.start=NULL),s)})
        toggle(c(checkb.history.estimates, checkb.history.estimatessadj,
                 checkb.history.estimatessadjchng, checkb.history.estimatestrend,
                 checkb.history.estimatestrendchng, checkb.history.estimatesseasonal,
                 checkb.history.estimatesfcst, checkb.history.estimatesaic,
                 checkb.history.fixmdl,
                 checkb.history.fixreg, checkb.history.fixreg1,
                 checkb.history.fixreg2, checkb.history.fixreg3,
                 checkb.history.fixreg4,
                 checkb.history.outlier, combobox.history.outlier,
				 checkb.history.target, combobox.history.target,
                 checkb.history.sadjlags, entry.history.sadjlags,
                 checkb.history.trendlags, entry.history.trendlags,
                 checkb.history.start, label.history.startyear,
                 label.history.startperiod, entry.history.startyear,
                 entry.history.startperiod), togglebutton)
        checkb.history.estimates$SetActive(FALSE)
		    checkb.history.estimatessadj$SetActive(FALSE)
		    checkb.history.estimatessadjchng$SetActive(FALSE)
        checkb.history.estimatestrend$SetActive(FALSE)
		    checkb.history.estimatestrendchng$SetActive(FALSE)
        checkb.history.estimatesseasonal$SetActive(FALSE)
		    checkb.history.estimatesfcst$SetActive(FALSE) 
        checkb.history.estimatesaic$SetActive(FALSE)
		    checkb.history.fixmdl$SetActive(FALSE)
		    checkb.history.fixreg$SetActive(FALSE)
        checkb.history.fixreg1$SetActive(FALSE)
		    checkb.history.fixreg2$SetActive(FALSE)
        checkb.history.fixreg3$SetActive(FALSE)
		    checkb.history.fixreg4$SetActive(FALSE)
		    checkb.history.outlier$SetActive(FALSE)
        combobox.history.outlier$SetActive(-1)		
		checkb.history.target$SetActive(FALSE)
		combobox.history.target$SetActive(-1)
		checkb.history.sadjlags$SetActive(FALSE)
        entry.history.sadjlags$SetText("")
		    checkb.history.trendlags$SetActive(FALSE)
        entry.history.trendlags$SetText("")
		    checkb.history.start$SetActive(FALSE)
		    entry.history.startyear$SetText("")
		    entry.history.startperiod$SetText("")
		  } 
		  status_print("history changed!")
		  togglebutton$SetInconsistent(FALSE)
		}
		#acceptdefault
		if(togglebutton == checkb.acceptdefault){
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(automdl.acceptdefault=TRUE),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(automdl.acceptdefault=FALSE),s)})
			status_print("automdl.acceptdefault changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#balanced
		if(togglebutton == checkb.balanced){
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(automdl.balanced=TRUE),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(automdl.balanced=FALSE),s)})
			status_print("automdl.balanced changed!")
			togglebutton$SetInconsistent(FALSE)
		}
#		#seats
#		if(togglebutton == checkb.seats){
#			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(seats=TRUE),s)})
#			else lapply(indices, FUN=function(s){object <<- setP(object,list(seats=FALSE),s)})
#			status_print("seats changed!")
#			togglebutton$SetInconsistent(FALSE)
#		}
		#estimate
		if(togglebutton == checkb.estimate){
		  toggle(c(checkb.estOutofsample), togglebutton)
			if(togglebutton$GetActive()==TRUE){
        lapply(indices, FUN=function(s){object <<- setP(object,list(estimate=TRUE),s)})
			}
			else{
			  lapply(indices, FUN=function(s){object <<- setP(object,list(estimate=FALSE,
                                                                    estimate.outofsample=TRUE),s)})
			  gSignalHandlerBlock(checkb.estOutofsample, handlercheckb.estOutofsample)
        checkb.estOutofsample$SetActive(TRUE)
			  gSignalHandlerUnblock(checkb.estOutofsample, handlercheckb.estOutofsample)
			} 
			status_print("estimate changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#estOutofsample
		if(togglebutton == checkb.estOutofsample){
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(estimate.outofsample=TRUE),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(estimate.outofsample=FALSE),s)})
			status_print("estimate.outofsample changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#slidingspans
		if(togglebutton == checkb.slidingspans){
			if(togglebutton$GetActive()==TRUE){
        lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans=TRUE),s)})
        toggle(c(checkb.slidingspans.fixmdl, checkb.slidingspans.fixreg, 
                 checkb.slidingspans.length, checkb.slidingspans.numspans, 
                 checkb.slidingspans.outlier, checkb.slidingspans.additivesa,
                 checkb.slidingspans.start), togglebutton)
			}
			else{
			  lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans=FALSE,  		
                                                                    slidingspans.fixmdl=NULL,
			                                                              slidingspans.fixreg=NULL,
			                                                              slidingspans.length=NULL,
			                                                              slidingspans.numspans=NULL,
			                                                              slidingspans.outlier=NULL,
			                                                              slidingspans.additivesa=NULL,
			                                                              slidingspans.start=NULL),s)})
			  toggle(c(combobox.slidingspans.fixmdl,checkb.slidingspans.fixmdl, 
                 checkb.slidingspans.fixreg, checkb.slidingspans.fixreg1,
			           checkb.slidingspans.fixreg2, checkb.slidingspans.fixreg3,
			           checkb.slidingspans.fixreg4, 
                 checkb.slidingspans.length, entry.slidingspans.length,
			           checkb.slidingspans.numspans, entry.slidingspans.numspans,
			           checkb.slidingspans.outlier, combobox.slidingspans.outlier,
			           checkb.slidingspans.additivesa, combobox.slidingspans.additivesa, 
			           checkb.slidingspans.start, entry.slidingspans.start1,
			           entry.slidingspans.start2, label.slidingspans.start1,
			           label.slidingspans.start2), togglebutton)
        combobox.slidingspans.fixmdl$SetActive(-1)
			  combobox.slidingspans.outlier$SetActive(-1)
			  combobox.slidingspans.additivesa$SetActive(-1)
        checkb.slidingspans.fixreg1$SetActive(FALSE)
			  checkb.slidingspans.fixreg2$SetActive(FALSE)
			  checkb.slidingspans.fixreg3$SetActive(FALSE)
			  checkb.slidingspans.fixreg4$SetActive(FALSE)
        entry.slidingspans.length$SetText("")
			  entry.slidingspans.numspans$SetText("")
			  entry.slidingspans.start1$SetText("")
			  entry.slidingspans.start2$SetText("")
			  checkb.slidingspans.fixmdl$SetActive(FALSE)
        checkb.slidingspans.fixreg$SetActive(FALSE) 
			  checkb.slidingspans.length$SetActive(FALSE) 
        checkb.slidingspans.numspans$SetActive(FALSE) 
			  checkb.slidingspans.outlier$SetActive(FALSE) 
        checkb.slidingspans.additivesa$SetActive(FALSE)
			  checkb.slidingspans.start$SetActive(FALSE)
			} 
			status_print("slidingspans changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#slidingspans.fixreg
		if(togglebutton == checkb.slidingspans.fixreg1 || togglebutton == checkb.slidingspans.fixreg2 ||
		   togglebutton == checkb.slidingspans.fixreg3 || togglebutton == checkb.slidingspans.fixreg4){
      values <- NULL
      if(checkb.slidingspans.fixreg1$GetActive()==TRUE)values <- c(values,"td")
      if(checkb.slidingspans.fixreg2$GetActive()==TRUE)values <- c(values,"holiday")
      if(checkb.slidingspans.fixreg3$GetActive()==TRUE)values <- c(values,"outlier")
      if(checkb.slidingspans.fixreg4$GetActive()==TRUE)values <- c(values,"user")
		  lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.fixreg=values),s)})
		  status_print("slidingspans.fixreg changed!")
      checkb.slidingspans.fixreg1$SetInconsistent(FALSE)
      checkb.slidingspans.fixreg2$SetInconsistent(FALSE)
      checkb.slidingspans.fixreg3$SetInconsistent(FALSE)
      checkb.slidingspans.fixreg4$SetInconsistent(FALSE)
      checkb.slidingspans.fixreg$SetInconsistent(FALSE)
		}
    
		#history.fixreg
		if(togglebutton == checkb.history.fixreg1 || togglebutton == checkb.history.fixreg2 ||
		  togglebutton == checkb.history.fixreg3 || togglebutton == checkb.history.fixreg4){
		  values <- NULL
		  if(checkb.history.fixreg1$GetActive()==TRUE)values <- c(values,"td")
		  if(checkb.history.fixreg2$GetActive()==TRUE)values <- c(values,"holiday")
		  if(checkb.history.fixreg3$GetActive()==TRUE)values <- c(values,"outlier")
		  if(checkb.history.fixreg4$GetActive()==TRUE)values <- c(values,"user")
		  lapply(indices, FUN=function(s){object <<- setP(object,list(history.fixreg=values),s)})
		  status_print("history.fixreg changed!")
		  checkb.history.fixreg1$SetInconsistent(FALSE)
		  checkb.history.fixreg2$SetInconsistent(FALSE)
		  checkb.history.fixreg3$SetInconsistent(FALSE)
		  checkb.history.fixreg4$SetInconsistent(FALSE)
		  checkb.history.fixreg$SetInconsistent(FALSE)
		}
		#history.estimates
		if(togglebutton == checkb.history.estimatessadj || togglebutton == checkb.history.estimatessadjchng ||
		  togglebutton == checkb.history.estimatestrend || togglebutton == checkb.history.estimatestrendchng ||
		  togglebutton == checkb.history.estimatesseasonal || togglebutton == checkb.history.estimatesfcst|| 
      togglebutton == checkb.history.estimatesaic){
		  values <- NULL
		  if(checkb.history.estimatessadj$GetActive()==TRUE)values <- c(values,"sadj")
		  if(checkb.history.estimatessadjchng$GetActive()==TRUE)values <- c(values,"sadjchng")
		  if(checkb.history.estimatestrend$GetActive()==TRUE)values <- c(values,"trend")
		  if(checkb.history.estimatestrendchng$GetActive()==TRUE)values <- c(values,"trendchng")
		  if(checkb.history.estimatesseasonal$GetActive()==TRUE)values <- c(values,"seasonal")
		  if(checkb.history.estimatesfcst$GetActive()==TRUE)values <- c(values,"fcst")
		  if(checkb.history.estimatesaic$GetActive()==TRUE)values <- c(values,"aic")
		  lapply(indices, FUN=function(s){object <<- setP(object,list(history.estimates=values),s)})
		  status_print("history.estimates changed!")
		  checkb.history.estimates$SetInconsistent(FALSE)
		  checkb.history.estimatessadj$SetInconsistent(FALSE)
		  checkb.history.estimatessadjchng$SetInconsistent(FALSE)
		  checkb.history.estimatestrend$SetInconsistent(FALSE)
		  checkb.history.estimatestrendchng$SetInconsistent(FALSE)
		  checkb.history.estimatesseasonal$SetInconsistent(FALSE)
		  checkb.history.estimatesfcst$SetInconsistent(FALSE)
		  checkb.history.estimatesaic$SetInconsistent(FALSE)
		}
# 		#onlytd
# 		if(togglebutton == checkb.onlytd){
# 			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(onlytd=TRUE),s)})
# 			else lapply(indices, FUN=function(s){object <<- setP(object,list(onlytd=FALSE),s)})
# 			status_print("onlytd changed!")
# 			togglebutton$SetInconsistent(FALSE)
# 		}
		#sfshort
		if(togglebutton == checkb.sfshort){
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.sfshort=TRUE),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(x11.sfshort=FALSE),s)})
			status_print("x11.sfshort changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#x11appendfcst
		if(togglebutton == checkb.x11appendfcst){
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.appendfcst=TRUE),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(x11.appendfcst=FALSE),s)})
			status_print("x11.appendfcst changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#x11appendbcst
		if(togglebutton == checkb.x11appendbcst){
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.appendbcst=TRUE),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(x11.appendbcst=FALSE),s)})
			status_print("x11.appendbcst changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#x11excludefcst
		if(togglebutton == checkb.x11excludefcst){
			if(togglebutton$GetActive()==TRUE)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.excludefcst=TRUE),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(x11.excludefcst=FALSE),s)})
			status_print("x11.excludefcst changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#x11.final
		if(togglebutton == checkb.x11.finalAO || togglebutton == checkb.x11.finalLS ||
		  togglebutton == checkb.x11.finalTC || togglebutton == checkb.x11.finaluser ||
		  togglebutton == checkb.x11.finalnone){
		  values <- NULL
		  if(checkb.x11.finalAO$GetActive()==TRUE)values <- c(values,"AO")
		  if(checkb.x11.finalLS$GetActive()==TRUE)values <- c(values,"LS")
		  if(checkb.x11.finalTC$GetActive()==TRUE)values <- c(values,"TC")
		  if(checkb.x11.finaluser$GetActive()==TRUE)values <- c(values,"user")
		  if(checkb.x11.finalnone$GetActive()==TRUE)values <- c(values,"none")
		  lapply(indices, FUN=function(s){object <<- setP(object,list(x11.final=values),s)})
		  status_print("x11.final changed!")
		  checkb.x11.finalAO$SetInconsistent(FALSE)
		  checkb.x11.finalTC$SetInconsistent(FALSE)
		  checkb.x11.finalLS$SetInconsistent(FALSE)
		  checkb.x11.finaluser$SetInconsistent(FALSE)
		  checkb.x11.finalnone$SetInconsistent(FALSE)
		  checkb.x11.finalactive$SetInconsistent(FALSE)
		}
		#x11regress
		if(togglebutton == radiob.x11regression){
		  regression[indices] <<- TRUE
			if(togglebutton$GetActive()==TRUE){
        lapply(indices, FUN=function(s){object <<- setP(object,list(x11regression=TRUE),s)})
        entry.aictest$SetText("")
        entry.aictest$SetSensitive(FALSE)
        checkb.aictestactive$SetSensitive(FALSE)
        checkb.aictestactive$SetActive(FALSE)
			}
			else{
			  lapply(indices, FUN=function(s){object <<- setP(object,list(x11regression=FALSE),s)}) 
			  entry.aictest$SetText("")
			  entry.aictest$SetSensitive(FALSE)
			  checkb.aictestactive$SetSensitive(TRUE)
			  checkb.aictestactive$SetActive(FALSE)
			}
			status_print("x11regression changed!")
			togglebutton$SetInconsistent(FALSE)
		}
		#outlier
		if(togglebutton == checkb.outlierall || togglebutton == checkb.outlierAO ||
				togglebutton == checkb.outlierLS ||togglebutton == checkb.outlierTC){
      
			if(togglebutton == checkb.outlierall)toggle(c(checkb.outlierAO, checkb.outlierLS, checkb.outlierTC),checkb.outlierall, invert=TRUE)
			checkb.outlierall$SetInconsistent(FALSE)
			checkb.outlierAO$SetInconsistent(FALSE)
			checkb.outlierTC$SetInconsistent(FALSE)
			checkb.outlierLS$SetInconsistent(FALSE)
			checkb.outlieractive$SetInconsistent(FALSE)
			i <- numeric(0)
			if(checkb.outlierAO$GetActive()==TRUE)i <- append(i, "AO")
			if(checkb.outlierLS$GetActive()==TRUE)i <- append(i, "LS")
			if(checkb.outlierTC$GetActive()==TRUE)i <- append(i, "TC")
			if(checkb.outlierall$GetActive()==TRUE)i <- append(i, "all")
			if(length(i)>0)lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.types=i),s)})
			else lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.types=NULL),s)})
		}
	}
	
	comboboxx12handler <- function(widget, user_data){
		#transform
		if(widget == combobox.transform){
			if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(transform.function="auto"),s)})
			if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(transform.function="log"),s)})
			if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(transform.function="none"),s)})
			status_print("transform.function changed!")
		}
		
    #series.type
# 		if(widget == combobox.series.type && checkb.series.type$GetActive()==1){
# 		  if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(series.type="flow"),s)})
# 		  if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(series.type="stock"),s)})
# 		  status_print("series.type changed!")
# 		}
    
		#transform.adjust
		if(widget == combobox.transform.adjust && checkb.transform.adjust$GetActive()==1){
		  if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(transform.adjust="lom"),s)})
		  if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(transform.adjust="loq"),s)})
		  if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(transform.adjust="lpyear"),s)})
      status_print("transform.adjust changed!")
		}
    
		#centeruser
		if(widget == combobox.centeruser){
		  regression[indices] <<- TRUE
			if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(regression.centeruser="mean"),s)})
			if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(regression.centeruser="seasonal"),s)})
			status_print("regression.centeruser changed!")
		}
		
		#outlier_method
		if(widget == combobox.outliermethod){
			if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.method="addone"),s)})
			if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.method="addall"),s)})
			status_print("outlier.method changed!")
		}
    
		#slidingspans.fixmdl
		if(widget == combobox.slidingspans.fixmdl){
		  if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.fixmdl="yes"),s)})
		  if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.fixmdl="no"),s)})
		  if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.fixmdl="clear"),s)})
      checkb.slidingspans.fixmdl$SetInconsistent(FALSE)
		  status_print("slidingspans.fixmdl changed!")
		}
		#slidingspans.outlier
		if(widget == combobox.slidingspans.outlier){
		  if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.outlier="keep"),s)})
		  if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.outlier="remove"),s)})
		  if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.outlier="yes"),s)})
		  checkb.slidingspans.outlier$SetInconsistent(FALSE)
		  status_print("slidingspans.outlier changed!")
		}
		#slidingspans.additivesa
		if(widget == combobox.slidingspans.additivesa){
		  if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.additivesa="differences"),s)})
		  if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.additivesa="percent"),s)})
		  checkb.slidingspans.additivesa$SetInconsistent(FALSE)
		  status_print("slidingspans.additivesa changed!")
		}
		
		#history.outlier
		if(widget == combobox.history.outlier){
		  if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(history.outlier="keep"),s)})
		  if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(history.outlier="remove"),s)})
		  if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(history.outlier="yes"),s)})
		  checkb.history.outlier$SetInconsistent(FALSE)
		  status_print("history.outlier changed!")
		}

		#history.target
		if(widget == combobox.history.target){
			if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(history.target="final"),s)})
			if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(history.target="concurrent"),s)})
#			if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(history.target="yes"),s)})
			checkb.history.target$SetInconsistent(FALSE)
			status_print("history.target changed!")
		}
		
		#samode
		if(widget == combobox.samode){
			if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.samode="mult"),s)})
			if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.samode="add"),s)})
			if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.samode="pseudoadd"),s)})
			if(widget$GetActive()==3)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.samode="logadd"),s)})
			status_print("x11.samode changed!")
		}
		
		#x11.type
		if(widget == combobox.x11.type){
		  if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.type="summary"),s)})
		  if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.type="trend"),s)})
		  if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.type="sa"),s)})
		  checkb.x11.type$SetInconsistent(FALSE)
		  status_print("x11.type changed!")
		}
    
		#x11calendarsigma
		if(widget == combobox.x11calendarsigma){
			if(widget$GetActive()==0)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.calendarsigma="all"),s)})
			if(widget$GetActive()==1)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.calendarsigma="signif"),s)})
			if(widget$GetActive()==2)lapply(indices, FUN=function(s){object <<- setP(object,list(x11.calendarsigma="select"),s)})
			status_print("x11.calendarsigma changed!")
		}
	}
	
	#Handler for timeseries table, called after new selection occurred
	#retrieves selection indices und uses them for updating plots and x12 parametergui
	tablehandler <- function(treeselection, userdata,...){
		indices <<- sapply(treeselection$GetSelectedRows()$retval, FUN=function(s){s$GetIndices()})+1
#		print(indices)
		update_notebook()
		make_history()
    text <- capture.output(read_x12(object, indices))
		status_print(text)	
# 		read_x12(object, indices)
		update_outliertable()
	}
	
	#Handler responsible for the mouseclicks on plot surfaces
	mousehandlerdrawing <- function(widget,event,userdata){
		if(widget == area.plot1){
			if(event$button == 3){
				menu.contextplotall$Popup(button=3, activate.time=0)
			}
		}
		if(widget == area.plot2){
			if(event$button == 3){
				menu.contextplotall$Popup(button=3, activate.time=0)
			}
		}
		if(widget == area.plot3){
			if(event$button == 3){
				menu.contextplotall$Popup(button=3, activate.time=0)
			}
		}
		if(widget == area.plot4){
			if(event$button == 1){
				if(locate == TRUE){
					locate <<-  FALSE
					button.manualoutlieraddclick$SetActive(FALSE)
					oc <- convertCoordinates(event$x)
					entry.manualoutlieryear$SetText(oc$year)
					entry.manualoutlierperiod$SetText(oc$period)
				}
        else{
          #lets assume that when somebody left clicked the plot and has not choosen to locate a outlier
          #then he wants to zoom in
          if(event$state==0){
            #mousedown:
            zooming <<- TRUE
            zoomX <<- event$X
            zoomY <<- event$Y
            size <- gdkDrawableGetSize(widget[["window"]])
            pixbuf.plot <<- gdkPixbufGetFromDrawable(NULL, widget[["window"]], NULL, 
                                                     0 ,0 ,0, 0, size$width, size$height)
            handler.moving <<- gSignalConnect(area.plot4, "motion_notify_event", f=mousehandlermoving)
          }
          else{
            #mouseup
            zooming <<- FALSE
            #minimal zooming rectangle size to reduce change for accidental zooming
            #exact size has to be adjusted!!
            gSignalHandlerDisconnect(area.plot4, handler.moving)
            gdkDrawPixbuf(widget[["window"]], NULL, pixbuf.plot, 0 ,0, 0, 0)
            if(abs((zoomX-event$x))*abs((zoomY-event$y)) > 10){
              start <- convertCoordinates(min(zoomX, event$x), zoom=TRUE)
              end <- convertCoordinates(max(zoomX, event$x), zoom=TRUE)
              t <- times(object@x12List[[min(indices)]])
              if(!is.null(t$backcast)){
                slider.plotmin$SetValue(calcPeriods(c(t$backcast[1], t$backcast[2], 
                                                      start$year, start$period), 
                                                    object@x12List[[indices[1]]]@ts))
                slider.plotmax$SetValue(calcPeriods(c(t$backcast[1], t$backcast[2], 
                                                      end$year, end$period), 
                                                    object@x12List[[indices[1]]]@ts)) 
              }
              else{
                slider.plotmin$SetValue(calcPeriods(c(as.numeric(t$original[1]), t$original[2], 
                                                      start$year, start$period), 
                                                    object@x12List[[indices[1]]]@ts))
                slider.plotmax$SetValue(calcPeriods(c(as.numeric(t$original[1]), t$original[2], 
                                                      end$year, end$period), 
                                                    object@x12List[[indices[1]]]@ts))                
              }
              
            }
          }
        }
			}
			if(event$button == 3 && event$state!=0){
				clickedX <<- event$x
				menu.contextplotwithoutlier$Popup(button=3, activate.time=0)
			}
		}
	}
	
  mousehandlermoving <- function(widget,event,userdata){
    gdkDrawPixbuf(widget[["window"]], NULL, pixbuf.plot, 0 ,0, 0, 0)
    gc <- gdkGCNew(widget[["window"]])
    gc$SetLineAttributes(line.width=1, line.style=GdkLineStyle["on-off-dash"],"round","round")
    gdkDrawRectangle(widget[["window"]], gc,
                     FALSE, min(zoomX, event$x), min(zoomY, event$y), abs(event$x-zoomX), abs(event$y-zoomY))
  }
  
	#Trys to Interpolate pixelcoordinates/usercoordinates to timeseries-values
	#could use some tweaking
	#######################################
	convertCoordinates <- function(x, zoom=FALSE){
		click <- x/area.plot4$GetAllocation()$allocation$width
		usract <- par("usr")[1:2]
		pltact <- par("plt")[1:2]
		xnew <- (((click -pltact[1])/(pltact[2] -pltact[1]))*(usract[2] -usract[1]))+usract[1]
		start_year <- floor(usract[1])
		end_year <- ceiling(usract[2])
		mm <- (0:(frequency(object@x12List[[min(indices)]]@ts)-1))/frequency(object@x12List[[min(indices)]]@ts)
		start_year:end_year
		tt <- expand.grid(start_year:end_year,mm)
		tt <- tt[,1]+tt[,2]
		tt <- tt[which.min(abs(xnew-tt))]
		year <- floor(tt)
		month <- round(1+(tt-year)*frequency(object@x12List[[min(indices)]]@ts))
		timess <- time(object@x12List[[min(indices)]]@ts)
		timeend <- timess[length(timess)]
		timestart <- timess[1]
		time <- timeend[length(timeend)]
		if(!zoom){
		  if(tt>timeend){
		    ee <- end(object@x12List[[min(indices)]]@ts)
		    year <- ee[1]
		    month <- ee[2]
		  }else if(tt<timestart){
		    ss <- start(object@x12List[[min(indices)]]@ts)
		    year <- ss[1]
		    month <- ss[2]
		  } 
		}
		
		ret <- list(year=year, period=month)
	}
	#######################################
	
	#reads values of the manualoutlierlist into the corresponding databuffer 
	update_outliertable <- function(){
		tablemodel.manualoutlier$Clear()
		if(length(outlierlist)>0){
			sapply(outlierlist,
					function(string) {
						## Add a new row to the model
						iter <- tablemodel.manualoutlier$Append()$iter
						tablemodel.manualoutlier$Set(iter, 0, string[1], 1, string[2], 2, string[3])
					})
		}
	}
	
	#encapsulates all NULL-checkbox updates in one function for simplified usage
	update_toggle <- function(button=NULL, system=TRUE){
		if(button==checkb.spanactive){
			toggle(c(entry.spanstartyear, entry.spanstartperiod, entry.spanendyear, 
							entry.spanendperiod, checkb.spanstart, checkb.spanend, 
               label.spanstartyear, label.spanstartperiod,
			         label.spanendyear, label.spanendperiod), button)
			if(button$GetActive()==FALSE){
				lapply(indices, FUN=function(s){object <<- setP(object,list(series.span=NULL),s)})
        entry.spanstartyear$SetText("")
				entry.spanstartperiod$SetText("")
				entry.spanendyear$SetText("")
				entry.spanendperiod$SetText("")
        checkb.spanstart$SetActive(FALSE)
				checkb.spanend$SetActive(FALSE)
				button$SetInconsistent(FALSE)
				if(system==FALSE)status_print("series.span changed!")
			}
			else{
				gSignalEmit(entry.spanstartyear, "changed")
        gSignalEmit(checkb.spanstart, "toggled")
				gSignalEmit(checkb.spanend, "toggled")
			}
		}
#     if(button==checkb.series.type){
#       toggle(c(combobox.series.type), button)
#       if(button$GetActive()==FALSE){
#         lapply(indices, FUN=function(s){object <<- setP(object,list(series.type=NULL),s)})
#         combobox.series.type$SetActive(-1)
#         button$SetInconsistent(FALSE)
#         if(system==FALSE)status_print("series.type changed!")
#       }
#       gSignalEmit(combobox.series.type, "changed")
#     }
    
    #transform.power
		if(button==checkb.transform.power){
		  toggle(c(entry.transform.power), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(transform.power=NULL),s)})
        entry.transform.power$SetText("")
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("transform.power changed!")
		  }
		  gSignalEmit(entry.transform.power, "changed")
		}
    
    #transform.adjust
		if(button==checkb.transform.adjust){
		  toggle(c(combobox.transform.adjust), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(transform.adjust=NULL),s)})
        combobox.transform.adjust$SetActive(-1)
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("transform.adjust changed!")
		  }
		  gSignalEmit(combobox.transform.adjust, "changed")
		}
    
		if(button==checkb.modelspanactive){
			toggle(c(entry.modelspanstartyear, entry.modelspanstartperiod, entry.modelspanendyear, 
							entry.modelspanendperiod, checkb.modelspanstart, checkb.modelspanend,
			         label.modelspanstartyear, label.modelspanstartperiod,
			         label.modelspanendyear, label.modelspanendperiod), button)
			if(button$GetActive()==FALSE){
				lapply(indices, FUN=function(s){object <<- setP(object,list(series.modelspan=NULL),s)})
				button$SetInconsistent(FALSE)
				entry.modelspanstartyear$SetText("")
				entry.modelspanstartperiod$SetText("")
				entry.modelspanendyear$SetText("")
				entry.modelspanendperiod$SetText("")
				checkb.modelspanstart$SetActive(FALSE)
				checkb.modelspanend$SetActive(FALSE)
				button$SetInconsistent(FALSE)
				if(system==FALSE)status_print("series.modelspan changed!")
			}
			else{
				gSignalEmit(entry.modelspanstartyear, "changed")
				gSignalEmit(checkb.modelspanstart, "toggled")
				gSignalEmit(checkb.modelspanend, "toggled")
			}
		}
		if(button==checkb.spanstart){
			toggle(c(entry.spanstartyear, entry.spanstartperiod, 
               label.spanstartyear, label.spanstartperiod), button)
      if(checkb.spanstart$GetActive()==FALSE){
        entry.spanstartyear$SetText("")
        entry.spanstartperiod$SetText("")
      }
			gSignalEmit(entry.spanstartyear, "changed")
		}
		if(button==checkb.spanend){
			toggle(c(entry.spanendyear, entry.spanendperiod, 
			         label.spanendyear, label.spanendperiod), button)
			if(checkb.spanend$GetActive()==FALSE){
			  entry.spanendyear$SetText("")
			  entry.spanendperiod$SetText("")
			}
			gSignalEmit(entry.spanstartyear, "changed")
		}
		if(button==checkb.modelspanstart){
			toggle(c(entry.modelspanstartyear, entry.modelspanstartperiod,
               label.modelspanstartyear, label.modelspanstartperiod), button)
			if(checkb.modelspanstart$GetActive()==FALSE){
			  entry.modelspanstartyear$SetText("")
			  entry.modelspanstartperiod$SetText("")
			}
			gSignalEmit(entry.modelspanstartyear, "changed")
		}
		if(button==checkb.modelspanend){
			toggle(c(entry.modelspanendyear, entry.modelspanendperiod,
               label.modelspanendyear, label.modelspanendperiod), button)
			if(checkb.modelspanend$GetActive()==FALSE){
			  entry.modelspanendyear$SetText("")
			  entry.modelspanendperiod$SetText("")
			}
			gSignalEmit(entry.modelspanstartyear, "changed")
		}

		#check.maxlag
		if(button==checkb.check.maxlag){
		  toggle(c(entry.check.maxlag), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(check.maxlag=NULL),s)})
        entry.check.maxlag$SetText("")
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("check.maxlag changed!")
		  }
		  gSignalEmit(entry.check.maxlag, "changed")
		}
    
		if(button==checkb.arimaactive){
			toggle(c(entry.arima1, entry.arima2, entry.arima3), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
        entry.arima1$SetText("")
				entry.arima2$SetText("")
				entry.arima3$SetText("")
				lapply(indices, FUN=function(s){object <<- setP(object,list(arima.model=NULL),s)})
				if(system==FALSE)status_print("arima.model changed!")
			}
			else{
				gSignalEmit(entry.arima1, "changed")
			}
		}
		if(button==checkb.sarimaactive){
			toggle(c(entry.sarima1, entry.sarima2, entry.sarima3), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				entry.sarima1$SetText("")
				entry.sarima2$SetText("")
				entry.sarima3$SetText("")
				lapply(indices, FUN=function(s){object <<- setP(object,list(arima.smodel=NULL),s)})
				if(system==FALSE)status_print("arima.smodel changed!")
			}
			else{
				gSignalEmit(entry.sarima1, "changed")
			}
		}
		#arima.ar
		if(button==checkb.arima.ar){
		  toggle(c(entry.arima.ar), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(arima.ar=NULL),s)})
        entry.arima.ar$SetText("")
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("arima.ar changed!")
		  }
		  gSignalEmit(entry.arima.ar, "changed")
		}
		#arima.ma
		if(button==checkb.arima.ma){
		  toggle(c(entry.arima.ma), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(arima.ma=NULL),s)})
		    entry.arima.ma$SetText("")
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("arima.ma changed!")
		  }
		  gSignalEmit(entry.arima.ma, "changed")
		}
    ######REGRESSSION
		if(button==checkb.regvariablesactive){
		  toggle(c(entry.regvariables), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    entry.regvariables$SetText("")
		    lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=NULL),s)})
        outlierlist <<- list()
		    update_outliertable()
		    status_print("regression.variables changed!")
		  }
		  else{
		    regression[indices] <<- TRUE
		    gSignalEmit(entry.regvariables, "changed")
		  }
		}
		if(button==checkb.reguseractive){
			toggle(c(entry.reguser), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.user=NULL),s)})
        entry.reguser$SetText("")
				if(system==FALSE)status_print("regression.user changed!")
			}
			else{
			  regression[indices] <<- TRUE
				gSignalEmit(entry.reguser, "changed")
			}
		}
		if(button==checkb.regfileactive){
			toggle(c(filebutton.regfile), button)
      if(button$GetActive()==FALSE){
        filebutton.regfile$UnselectAll()
        button$SetInconsistent(FALSE)
        lapply(indices, FUN=function(s){object <<- setP(object,list(regression.file=NULL),s)})
        if(system==FALSE)status_print("regression.file changed!")
      }
      else{
        regression[indices] <<- TRUE
      }
		}
    #Regression END
    
		#identify.diff
		if(button==checkb.identify.diff){
		  toggle(c(entry.identify.diff), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(identify.diff=NULL),s)})
        entry.identify.diff$SetText("")
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("identify.diff changed!")
		  }
		  gSignalEmit(entry.identify.diff, "changed")
		}
		#identify.sdiff
		if(button==checkb.identify.sdiff){
		  toggle(c(entry.identify.sdiff), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(identify.sdiff=NULL),s)})
		    entry.identify.sdiff$SetText("")
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("identify.sdiff changed!")
		  }
		  gSignalEmit(entry.identify.sdiff, "changed")
		}
		#identify.maxlag
		if(button==checkb.identify.maxlag){
		  toggle(c(entry.identify.maxlag), button)
		  if(button$GetActive()==FALSE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(identify.maxlag=NULL),s)})
		    entry.identify.maxlag$SetText("")
		    button$SetInconsistent(FALSE)
		    if(system==FALSE)status_print("identify.maxlag changed!")
		  }
		  gSignalEmit(entry.identify.maxlag, "changed")
		}
    
		if(button==checkb.usertypeactive){
			toggle(c(entry.usertype), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.usertype=NULL),s)})
				entry.usertype$SetText("")
				if(system==FALSE)status_print("regression.usertype changed!")
			}
			else{
				gSignalEmit(entry.usertype, "changed")
			}
		}
		if(button==checkb.centeruseractive){
			toggle(c(combobox.centeruser), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.centeruser=NULL),s)})
        combobox.centeruser$SetActive(-1)
				if(system==FALSE)status_print("regression.centeruser changed!")
			}
			else{
				gSignalEmit(combobox.centeruser, "changed")
			}
		}
    
    #slidingspans.fixmdl
		if(button==checkb.slidingspans.fixmdl){
		  toggle(c(combobox.slidingspans.fixmdl), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.fixmdl=NULL),s)})
        combobox.slidingspans.fixmdl$SetActive(-1)
		    if(system==FALSE)status_print("slidingspans.fixmdl changed!")
		  }
		  else{
		    gSignalEmit(combobox.slidingspans.fixmdl, "changed")
		  }
		}
		#slidingspans.outlier
		if(button==checkb.slidingspans.outlier){
		  toggle(c(combobox.slidingspans.outlier), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.outlier=NULL),s)})
		    combobox.slidingspans.outlier$SetActive(-1)
		    if(system==FALSE)status_print("slidingspans.outlier changed!")
		  }
		  else{
		    gSignalEmit(combobox.slidingspans.outlier, "changed")
		  }
		}
		#slidingspans.additivesa
		if(button==checkb.slidingspans.additivesa){
		  toggle(c(combobox.slidingspans.additivesa), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.additivesa=NULL),s)})
		    combobox.slidingspans.additivesa$SetActive(-1)
		    if(system==FALSE)status_print("slidingspans.additivesa changed!")
		  }
		  else{
		    gSignalEmit(combobox.slidingspans.additivesa, "changed")
		  }
		}
		#slidingspans.fixreg
		if(button==checkb.slidingspans.fixreg){
		  toggle(c(checkb.slidingspans.fixreg1,checkb.slidingspans.fixreg2,
               checkb.slidingspans.fixreg3,checkb.slidingspans.fixreg4), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    checkb.slidingspans.fixreg1$SetActive(FALSE)
		    checkb.slidingspans.fixreg2$SetActive(FALSE)
		    checkb.slidingspans.fixreg3$SetActive(FALSE)
		    checkb.slidingspans.fixreg4$SetActive(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.fixreg=NULL),s)})
		    if(system==FALSE)status_print("slidingspans.fixreg changed!")
		  }
		}
		#slidingspans.length
		if(button==checkb.slidingspans.length){
		  toggle(c(entry.slidingspans.length), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.length=NULL),s)})
		    entry.slidingspans.length$SetText("")
		    if(system==FALSE)status_print("slidingspans.length changed!")
		  }
		  else{
		    gSignalEmit(entry.slidingspans.length, "changed")
		  }
		}
		#slidingspans.numspans
		if(button==checkb.slidingspans.numspans){
		  toggle(c(entry.slidingspans.numspans), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.numspans=NULL),s)})
		    entry.slidingspans.numspans$SetText("")
		    if(system==FALSE)status_print("slidingspans.numspans changed!")
		  }
		  else{
		    gSignalEmit(entry.slidingspans.numspans, "changed")
		  }
		}
		#slidingspans.start
		if(button==checkb.slidingspans.start){
		  toggle(c(entry.slidingspans.start1, entry.slidingspans.start2, 
               label.slidingspans.start1, label.slidingspans.start2), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.start=NULL),s)})
		    entry.slidingspans.start1$SetText("")
		    entry.slidingspans.start2$SetText("")
		    if(system==FALSE)status_print("slidingspans.start changed!")
		  }
		  else{
		    gSignalEmit(entry.slidingspans.start1, "changed")
		    gSignalEmit(entry.slidingspans.start2, "changed")
		  }
		}
		
		#history.estimates
		if(button==checkb.history.estimates){
		  toggle(c(checkb.history.estimatessadj, checkb.history.estimatessadjchng,
		           checkb.history.estimatestrend ,checkb.history.estimatestrendchng,
		           checkb.history.estimatesseasonal, checkb.history.estimatesfcst,
		           checkb.history.estimatesaic), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.estimates=NULL),s)})
		    checkb.history.estimatessadj$SetActive(FALSE)
		    checkb.history.estimatessadjchng$SetActive(FALSE)
		    checkb.history.estimatestrend$SetActive(FALSE)
		    checkb.history.estimatestrendchng$SetActive(FALSE)
		    checkb.history.estimatesseasonal$SetActive(FALSE)
		    checkb.history.estimatesfcst$SetActive(FALSE)
		    checkb.history.estimatesaic$SetActive(FALSE)
		    if(system==FALSE)status_print("history.estimates changed!")
		  }
		}
		#history.fixreg
		if(button==checkb.history.fixreg){
		  toggle(c(checkb.history.fixreg1, checkb.history.fixreg2,
		           checkb.history.fixreg3, checkb.history.fixreg4), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.fixreg=NULL),s)})
		    checkb.history.fixreg$SetActive(FALSE)
		    checkb.history.fixreg1$SetActive(FALSE)
		    checkb.history.fixreg2$SetActive(FALSE)
		    checkb.history.fixreg3$SetActive(FALSE)
		    checkb.history.fixreg4$SetActive(FALSE)
		    if(system==FALSE)status_print("history.fixreg changed!")
		  }
		}
		#history.outlier
		if(button==checkb.history.outlier){
		  toggle(c(combobox.history.outlier), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.outlier=NULL),s)})
		    combobox.history.outlier$SetActive(-1)
		    if(system==FALSE)status_print("history.outlier changed!")
		  }
		  else{
		    gSignalEmit(combobox.history.outlier, "changed")
		  }
		}
		#history.target
		if(button==checkb.history.target){
			toggle(c(combobox.history.target), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(history.target=NULL),s)})
				combobox.history.target$SetActive(-1)
				if(system==FALSE)status_print("history.target changed!")
			}
			else{
				gSignalEmit(combobox.history.target, "changed")
			}
		}
		
		#history.sadjlags
		if(button==checkb.history.sadjlags){
		  toggle(c(entry.history.sadjlags), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.sadjlags=NULL),s)})
		    entry.history.sadjlags$SetText("")
		    if(system==FALSE)status_print("history.sadjlags changed!")
		  }
		  else{
		    gSignalEmit(entry.history.sadjlags, "changed")
		  }
		}
		#history.trendlags
		if(button==checkb.history.trendlags){
		  toggle(c(entry.history.trendlags), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.trendlags=NULL),s)})
		    entry.history.trendlags$SetText("")
		    if(system==FALSE)status_print("history.trend changed!")
		  }
		  else{
		    gSignalEmit(entry.history.trendlags, "changed")
		  }
		}
		#history.start
		if(button==checkb.history.start){
		  toggle(c(entry.history.startyear, entry.history.startperiod, 
               label.history.startyear, label.history.startperiod), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.start=NULL),s)})
		    entry.history.startyear$SetText("")
		    entry.history.startperiod$SetText("")
		    if(system==FALSE)status_print("history.start changed!")
		  }
		  else{
		    gSignalEmit(entry.history.startyear, "changed")
		    gSignalEmit(entry.history.startperiod, "changed")
		  }
		}
# 		#maxorder
# 		if(button==checkb.maxorderactive){
# 		  toggle(c(entry.maxorder1, entry.maxorder2), button)
# 		  if(button$GetActive()==FALSE){
# 		    lapply(indices, FUN=function(s){object <<- setP(object,list(automdl.maxorder=NULL),s)})
# 		    button$SetInconsistent(FALSE)
# 		    status_print("automdl.maxorder changed!")
# 		  }
# 		  gSignalEmit(entry.maxorder1, "changed")
# 		  gSignalEmit(entry.maxorder2, "changed")
# 		}
    
		if(button==checkb.regfilestartactive){
			toggle(c(entry.regfilestartstartyear, entry.regfilestartstartperiod,
			         label.regfilestartstartyear, label.regfilestartstartperiod), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.start=NULL),s)})
				entry.regfilestartstartyear$SetText("")
				entry.regfilestartstartperiod$SetText("")
				if(system==FALSE)status_print("regression.start changed!")
			}
			else{
				gSignalEmit(entry.regfilestartstartyear, "changed")
			}
		}
#		if(button==checkb.seatsparameteractive){
#			toggle(c(entry.seatsparameter), button)
#			if(button$GetActive()==FALSE){
#				button$SetInconsistent(FALSE)
#				lapply(indices, FUN=function(s){object <<- setP(object,list(seatsparameter=NULL),s)})
#				entry.seatsparameter$SetText("")
#				if(system==FALSE)status_print("seatsparameter changed!")
#			}
#			else{
#				gSignalEmit(entry.seatsparameter, "changed")
#			}
#		}
		if(button==checkb.sigmalimactive){
			toggle(c(entry.sigmalim1, entry.sigmalim2), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.sigmalim=NULL),s)})
        entry.sigmalim1$SetText("")
				entry.sigmalim2$SetText("")
				if(system==FALSE)status_print("x11.sigmalim changed!")
			}
			else{
				gSignalEmit(entry.sigmalim1, "changed")
			}
		}
		if(button==checkb.criticalactive){
			toggle(c(radiob.criticalall, radiob.criticalspecific), button)
			gSignalEmit(radiob.criticalall, "toggled")
			gSignalEmit(radiob.criticalspecific, "toggled")
			if(button$GetActive()==FALSE){
				toggle(c(entry.criticalAO, entry.criticalTC, entry.criticalLS, 
								radiob.criticalspecific,entry.criticalall, radiob.criticalall), checkb.criticalactive)
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.critical=NULL),s)})
        entry.criticalall$SetText("")
				entry.criticalAO$SetText("")
				entry.criticalTC$SetText("")
				entry.criticalLS$SetText("")
        label.criticalAO$SetSensitive(FALSE)
				label.criticalLS$SetSensitive(FALSE)
				label.criticalTC$SetSensitive(FALSE)
				if(system==FALSE)status_print("outlier.critical changed!")
			}
     # print(object@x12List[[1]]@x12Parameter@outlier.critical)
		}
		if(button==radiob.criticalall){
			if(button$GetActive()==TRUE){
				toggle(c(entry.criticalall), radiob.criticalall)
				toggle(c(entry.criticalAO, entry.criticalTC, entry.criticalLS,
				         label.criticalAO, label.criticalTC, label.criticalLS), radiob.criticalspecific)
				entry.criticalAO$SetText("")
				entry.criticalTC$SetText("")
				entry.criticalLS$SetText("")
			}
		}
		if(button==radiob.criticalspecific){
			if(button$GetActive()==TRUE){
				toggle(c(entry.criticalall), radiob.criticalall)
				toggle(c(entry.criticalAO, entry.criticalTC, entry.criticalLS,
				         label.criticalAO, label.criticalTC, label.criticalLS), radiob.criticalspecific)
				entry.criticalall$SetText("")
			}
		}
		if(button==checkb.outlieractive){
			if(button$GetActive()==FALSE){
			  toggle(c(checkb.outlierall, checkb.outlierAO, checkb.outlierTC, checkb.outlierLS,
			           checkb.criticalactive, checkb.outlierspanactive, checkb.outliermethodactive,
			           entry.outlierspan.start1,entry.outlierspan.start2,
			           entry.outlierspan.end1,entry.outlierspan.end2,
			           combobox.outliermethod, label.outlier,
			           entry.criticalAO, entry.criticalTC, entry.criticalLS, 
			           radiob.criticalspecific, radiob.criticalspecific, 
			           entry.criticalall, radiob.criticalall,
			           label.criticalAO, label.criticalTC,label.criticalLS), button)
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.types=NULL, 
                                                                    outlier.span=NULL,
                                                                    outlier.method=NULL,
                                                                    outlier.critical=NULL),s)})
				gSignalHandlerBlock(checkb.outlierall, handler.outlierall)
				gSignalHandlerBlock(checkb.outlierTC, handler.outlierTC)
				gSignalHandlerBlock(checkb.outlierAO, handler.outlierAO)
				gSignalHandlerBlock(checkb.outlierLS, handler.outlierLS)
        checkb.criticalactive$SetActive(FALSE)
        checkb.outlierspanactive$SetActive(FALSE)
        checkb.outliermethodactive$SetActive(FALSE)
        checkb.outlierspan.start$SetActive(FALSE)
			  checkb.outlierspan.end$SetActive(FALSE)
				checkb.outlierall$SetActive(FALSE)
				checkb.outlierAO$SetActive(FALSE)
				checkb.outlierTC$SetActive(FALSE)
				checkb.outlierLS$SetActive(FALSE)
        entry.criticalall$SetText("")
				entry.criticalAO$SetText("")
				entry.criticalLS$SetText("")
				entry.criticalTC$SetText("")
        combobox.outliermethod$SetActive(-1)
				gSignalHandlerUnblock(checkb.outlierall, handler.outlierall)
				gSignalHandlerUnblock(checkb.outlierTC, handler.outlierTC)
				gSignalHandlerUnblock(checkb.outlierAO, handler.outlierAO)
				gSignalHandlerUnblock(checkb.outlierLS, handler.outlierLS)
				if(system==FALSE)status_print("outlier.types changed!")
			}
			else{
			  toggle(c( checkb.outlierall, checkb.outlierAO, checkb.outlierTC, checkb.outlierLS,
			            checkb.criticalactive, checkb.outlierspanactive, 
                  checkb.outliermethodactive, label.outlier), button)
  ############
  if(is.null(unlist(getP(object,list("outlier.types"))))){

  button$SetInconsistent(FALSE)
  lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.types="all", 
							  outlier.span=NULL,
							  outlier.method=NULL,
							  outlier.critical=NULL),s)})
  gSignalHandlerBlock(checkb.outlierall, handler.outlierall)
  gSignalHandlerBlock(checkb.outlierTC, handler.outlierTC)
  gSignalHandlerBlock(checkb.outlierAO, handler.outlierAO)
  gSignalHandlerBlock(checkb.outlierLS, handler.outlierLS)
  checkb.criticalactive$SetActive(FALSE)
  checkb.outlierspanactive$SetActive(FALSE)
  checkb.outliermethodactive$SetActive(FALSE)
  checkb.outlierspan.start$SetActive(FALSE)
  checkb.outlierspan.end$SetActive(FALSE)
  checkb.outlierall$SetActive(TRUE)
  checkb.outlierAO$SetActive(FALSE)
  checkb.outlierTC$SetActive(FALSE)
  checkb.outlierLS$SetActive(FALSE)
  entry.criticalall$SetText("")
  entry.criticalAO$SetText("")
  entry.criticalLS$SetText("")
  entry.criticalTC$SetText("")
  combobox.outliermethod$SetActive(-1)
  gSignalHandlerUnblock(checkb.outlierall, handler.outlierall)
  gSignalHandlerUnblock(checkb.outlierTC, handler.outlierTC)
  gSignalHandlerUnblock(checkb.outlierAO, handler.outlierAO)
  gSignalHandlerUnblock(checkb.outlierLS, handler.outlierLS)
  if(system==FALSE)status_print("outlier.types changed!")
}  
  ############
			}
		}
    if(button==checkb.outlierspan.start){
      toggle(c(entry.outlierspan.start1,entry.outlierspan.start2), checkb.outlierspan.start)
      if(button$GetActive()==FALSE){
        entry.outlierspan.start1$SetText("")
        entry.outlierspan.start2$SetText("")
      }
      gSignalEmit(entry.outlierspan.start1, "changed")
    }
    
		if(button==checkb.outlierspan.end){
		  toggle(c(entry.outlierspan.end1,entry.outlierspan.end2), checkb.outlierspan.end)
		  if(button$GetActive()==FALSE){
		    entry.outlierspan.end1$SetText("")
		    entry.outlierspan.end2$SetText("")
		  }
		  gSignalEmit(entry.outlierspan.end1, "changed")
		}
    
		if(button==checkb.outlierspanactive){
			#toggle(c(entry.outlierspan.start1, entry.outlierspan.start2), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
        checkb.outlierspan.start$SetSensitive(FALSE)
        checkb.outlierspan.end$SetSensitive(FALSE)
        entry.outlierspan.start1$SetSensitive(FALSE)
				entry.outlierspan.start2$SetSensitive(FALSE)
				entry.outlierspan.end1$SetSensitive(FALSE)
				entry.outlierspan.end2$SetSensitive(FALSE)
        label.outlierspan.year$SetSensitive(FALSE)
				label.outlierspan.period$SetSensitive(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.span=NULL),s)})
				entry.outlierspan.start1$SetText("")
				entry.outlierspan.start2$SetText("")
				entry.outlierspan.end1$SetText("")
				entry.outlierspan.end2$SetText("")
				if(system==FALSE)status_print("outlier.span changed!")
			}
			else{
			  checkb.outlierspan.start$SetSensitive(TRUE)
			  checkb.outlierspan.end$SetSensitive(TRUE)
			  label.outlierspan.year$SetSensitive(TRUE)
			  label.outlierspan.period$SetSensitive(TRUE)
				gSignalEmit(entry.outlierspan.start1, "changed")
        gSignalEmit(checkb.outlierspan.start, "toggled")
        gSignalEmit(checkb.outlierspan.end, "toggled")
			}
		}
		if(button==checkb.outliermethodactive){
			toggle(c(combobox.outliermethod), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.method=NULL),s)})
        combobox.outliermethod$SetActive(-1)
				if(system==FALSE)status_print("outlier.method changed!")
			}
			else{
				
				gSignalEmit(combobox.outliermethod, "changed")
			}
		}
#   if(button==checkb.fileactive){
#     toggle(c(entry.file), button)
#   }
		if(button==checkb.forecast_yearsactive){
			toggle(c(entry.forecast_years), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(forecast_years=NULL),s)})
        entry.forecast_years$SetText("")
				if(system==FALSE)status_print("forecast_years changed!")
			}
			else{
				gSignalEmit(entry.forecast_years, "changed")
			}
		}
		if(button==checkb.backcast_yearsactive){
			toggle(c(entry.backcast_years), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(backcast_years=NULL),s)})
				entry.backcast_years$SetText("")
				if(system==FALSE)status_print("backcast_years changed!")
			}
			else{
				gSignalEmit(entry.backcast_years, "changed")
			}
		}
#   if(button==checkb.forecast_confactive){
#     toggle(c(entry.forecast_conf), button)
#     if(button$GetActive()==FALSE){
#       button$SetInconsistent(FALSE)
#       lapply(indices, FUN=function(s){x12 <<- setP(x12,list(forecast_conf=NULL),s)})
#       status_print("forecast_conf changed!")
#     }
#     else{
#       gSignalEmit(entry.forecast_conf, "changed")
#     }
#   }
		if(button==checkb.aictestactive){
			toggle(c(entry.aictest), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.aictest=NULL),s)})
        entry.aictest$SetText("")
        entry.aictest$SetSensitive(FALSE)
				if(system==FALSE)status_print("regression.aictest changed!")
			}
			else{
				gSignalEmit(entry.aictest, "changed")
			}
		}
		if(button==checkb.seasonalmaactive){
			toggle(c(entry.seasonalma), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.seasonalma=NULL),s)})
        entry.seasonalma$SetText("")
				if(system==FALSE)status_print("x11.seasonalma changed!")
			}
			else{
				gSignalEmit(entry.seasonalma, "changed")
			}
		}
		if(button==checkb.trendmaactive){
			toggle(c(entry.trendma), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.trendma=NULL),s)})
        entry.trendma$SetText("")
				if(system==FALSE)status_print("x11.trendma changed!")
			}
			else{
				gSignalEmit(entry.trendma, "changed")
			}
		}
		#x11.type
		if(button==checkb.x11.type){
		  toggle(c(combobox.x11.type), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(x11.type=NULL),s)})
		    combobox.x11.type$SetActive(-1)
		    if(system==FALSE)status_print("x11.type changed!")
		  }
		  else{
		    gSignalEmit(combobox.x11.type, "changed")
		  }
		}
		#x11.final
		if(button==checkb.x11.finalactive){
		  toggle(c(checkb.x11.finalAO, checkb.x11.finalLS,
		           checkb.x11.finalTC, checkb.x11.finalnone, checkb.x11.finaluser), button)
		  if(button$GetActive()==FALSE){
		    button$SetInconsistent(FALSE)
		    lapply(indices, FUN=function(s){object <<- setP(object,list(x11.final=NULL),s)})
#TODO: Kommentar eingefuegt, weil ein nicht vorhandes Feld aufgeruefen wird        
#		    checkb.x11.final$SetActive(FALSE)
		    checkb.x11.finalAO$SetActive(FALSE)
		    checkb.x11.finalLS$SetActive(FALSE)
		    checkb.x11.finalTC$SetActive(FALSE)
		    checkb.x11.finalnone$SetActive(FALSE)
		    checkb.x11.finaluser$SetActive(FALSE)
		    if(system==FALSE)status_print("x11.final changed!")
		  }
		}
		if(button==checkb.x11calendarsigmaactive){
			toggle(c(combobox.x11calendarsigma), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.calendarsigma=NULL),s)})
				if(system==FALSE)status_print("x11.calendarsigma changed!")
			}
			else{
				gSignalEmit(combobox.x11calendarsigma, "changed")
			}
		}
		if(button==checkb.samodeactive){
			toggle(c(combobox.samode), button)
			if(button$GetActive()==FALSE){
				button$SetInconsistent(FALSE)
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.samode=NULL),s)})
        combobox.samode$SetActive(-1)
				if(system==FALSE)status_print("x11.samode changed!")
			}
			else{
				gSignalEmit(combobox.samode, "changed")
			}
		}
	}
	
	#handler for the update button
	#calls x12 and updates summaries,plots and tables
	update_handler <- function(button, userdata){
		object <<- x12(object) 
		update_notebook(reset=FALSE)
		update_outliertable()
		make_history()
	}
	
	#handler for the buttons in the manual outlier panel
	manualoutlierhandler <- function(button, userdata){
		if(button==button.manualoutlieradd){
			if(!isEmpty(entry.manualoutlieryear$GetText())&&!isEmpty(entry.manualoutlierperiod$GetText())&&
					isNumber(entry.manualoutlieryear$GetText(), integer=TRUE) &&
					isPeriod(entry.manualoutlierperiod$GetText()) &&
					combobox.manualoutliertype$GetActive()!=-1){
				if(combobox.manualoutliertype$GetActive()==0) typ <- "TC"
				if(combobox.manualoutliertype$GetActive()==1) typ <- "LS"
				if(combobox.manualoutliertype$GetActive()==2) typ <- "AO"
				outlierlist[[length(outlierlist)+1]] <<- c(typ,as.numeric(entry.manualoutlieryear$GetText()), entry.manualoutlierperiod$GetText())
				checkb.regvariablesactive$SetActive(TRUE);
				
				update_regvariables()
			}
		}
		
		if(button==button.manualoutlierremove){
			selection <- table.manualoutlier$GetSelection()
			if(selection$CountSelectedRows()>0){
				rows <- sapply(selection$GetSelectedRows()$retval, FUN=function(s){s$GetIndices()})+1
				outlierlist[[rows[1]]] <<- NULL
				update_regvariables()
			}
			
		}
		
		if(button == button.manualoutlieraddclick){
			if(button$GetActive()==TRUE){
				locate <<- TRUE
			}
		}
	}
	
	#updates the regvariables parameter with values from the regvariables entry and the
	#manual outlier panel
	update_regvariables <- function(){
		vars <- sapply(outlierlist, FUN=function(s){
					paste(s[1],s[2],".",s[3], sep="")
				})
		if(trim(entry.regvariables$GetText())=="" && length(vars)!=0){
			lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=as.vector(vars)),s)})
		}
		else if((trim(entry.regvariables$GetText())!="" && trim(entry.regvariables$GetText())!="*") && length(vars)==0){
			lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=cutParam(entry.regvariables$GetText())),s)})
		}
		else if((trim(entry.regvariables$GetText())=="") && length(vars)==0){
			lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=NULL),s)})
		}
		else{
			if(trim(entry.regvariables$GetText())!="*")lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=append(cutParam(entry.regvariables$GetText()), as.vector(vars))),s)})
		}
		update_outliertable()
	}
	
	#Handler for the change event of all x12parameter entrys
	#validates input, if input is valid setP is called
	x12input_handler <- function(editable, user_data){
		text <- editable$GetChars(0, 30)
		
		if(editable==entry.spanstartyear||editable==entry.spanstartperiod||editable==entry.spanendyear||editable==entry.spanendperiod){
			startyear <- entry.spanstartyear$GetChars(0, -1)
			startperiod <- entry.spanstartperiod$GetChars(0, -1)
			endyear <- entry.spanendyear$GetChars(0, -1)
			endperiod <- entry.spanendperiod$GetChars(0, -1)
			values <- c(NA,NA,NA,NA)
			problem = FALSE
			if(checkb.spanstart$GetActive()){
				if(isNumber(startyear)&&isPeriod(startperiod)){
					values[1] <- as.numeric(startyear)
					values[2] <- startperiod
				}
				else{
					problem = TRUE
					status_print("possible wrong input for series.span!")
				}
			}
			if(checkb.spanend$GetActive()){
				if(isNumber(endyear)&&isPeriod(endperiod)){
					values[3] <- as.numeric(endyear)
					values[4] <- endperiod
				}
				else{
					problem=TRUE
					status_print("possible wrong input for series.span!")
				}
			}
#     print(values)
			if(problem==FALSE){
				status_print("series.span changed!")
				if(length(unique(values))==1 && is.na(values[1]))lapply(indices, FUN=function(s){object <<- setP(object,list(series.span=NULL),s)})
				else lapply(indices, FUN=function(s){object <<- setP(object,list(series.span=values),s)})
			}
		}
		
		if(editable==entry.modelspanstartyear||editable==entry.modelspanstartperiod||editable==entry.modelspanendyear||editable==entry.modelspanendperiod){
			startyear <- entry.modelspanstartyear$GetChars(0, -1)
			startperiod <- entry.modelspanstartperiod$GetChars(0, -1)
			endyear <- entry.modelspanendyear$GetChars(0, -1)
			endperiod <- entry.modelspanendperiod$GetChars(0, -1)
			values <- c(NA,NA,NA,NA)
			problem = FALSE
			if(checkb.modelspanstart$GetActive()){
				if(isNumber(startyear)&&isPeriod(startperiod)){
					values[1] <- as.numeric(startyear)
					values[2] <- (startperiod)
				}
				else{
					problem = TRUE
					status_print("possible wrong input for series.modelspan!")
				}
			}
			if(checkb.modelspanend$GetActive()){
				if(isNumber(endyear)&&isPeriod(endperiod)){
					values[3] <- as.numeric(endyear)
					values[4] <- (endperiod)
				}
				else{
					problem=TRUE
					status_print("possible wrong input for series.modelspan!")
				}
			}
#     print(values)
			if(problem==FALSE){
				status_print("series.modelspan changed!")
				if(length(unique(values))==1 && is.na(values[1]))lapply(indices, FUN=function(s){object <<- setP(object,list(series.modelspan=NULL),s)})
				else lapply(indices, FUN=function(s){object <<- setP(object,list(series.modelspan=values),s)})
			}
		}
		
#		if(editable==entry.decimals){
#			if(isNumber(text)){
#				lapply(indices, FUN=function(s){object <<- setP(object,list(decimals=as.numeric(text)),s)})
#				status_print("decimals changed!")
#			}
#			else{
#				status_print("possible wrong input for decimals!")
#			}
#		}
		
		if(editable==entry.arima1||editable==entry.arima2||editable==entry.arima3){
			text <- entry.arima1$GetChars(0, -1)
			text2 <- entry.arima2$GetChars(0, -1)
			text3 <- entry.arima3$GetChars(0, -1)
			if(isNumber(text)&&isNumber(text2)&&isNumber(text3)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(arima.model=c(as.numeric(text), as.numeric(text2), as.numeric(text3))),s)})
				status_print("mode.arima changed!")
			}
			else{
				status_print("possible wrong input for arima!")
			}
		}
		
    #transform.power
		if(editable==entry.transform.power){
		  text <- entry.transform.power$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE) && checkb.transform.power$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(transform.power=as.numeric(text)),s)})
		    status_print("transform.power changed!")
		  }
		  else{
		    status_print("possible wrong input for transform.power!")
		  }
		}
    
		#check.maxlag
		if(editable==entry.check.maxlag){
		  text <- entry.check.maxlag$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE) && checkb.check.maxlag$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(check.maxlag=as.numeric(text)),s)})
        checkb.check.maxlag$SetInconsistent(FALSE)
		    status_print("check.maxlag changed!")
		  }
		  else{
		    status_print("possible wrong input for check.maxlag!")
		  }
		}
    
		#identify.diff
		if(editable==entry.identify.diff){
		  text <- entry.identify.diff$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE) && checkb.identify.diff$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(identify.diff=as.numeric(text)),s)})
		    status_print("identify.diff changed!")
		  }
		  else{
		    status_print("possible wrong input for identify.diff!")
		  }
		}
		#identify.sdiff
		if(editable==entry.identify.sdiff){
		  text <- entry.identify.sdiff$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE) && checkb.identify.sdiff$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(identify.sdiff=as.numeric(text)),s)})
		    status_print("identify.sdiff changed!")
		  }
		  else{
		    status_print("possible wrong input for identify.sdiff!")
		  }
		}
		#identify.maxlag
		if(editable==entry.identify.maxlag){
		  text <- entry.identify.maxlag$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE) && checkb.identify.maxlag$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(identify.maxlag=as.numeric(text)),s)})
		    status_print("identify.maxlag changed!")
		  }
		  else{
		    status_print("possible wrong input for identify.maxlag!")
		  }
		}
    
		if(editable==entry.sarima1||editable==entry.sarima2||editable==entry.sarima3){
			text <- entry.sarima1$GetChars(0, -1)
			text2 <- entry.sarima2$GetChars(0, -1)
			text3 <- entry.sarima3$GetChars(0, -1)
			if(isNumber(text)&&isNumber(text2)&&isNumber(text3)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(arima.smodel=c(as.numeric(text), as.numeric(text2), as.numeric(text3))),s)})
				status_print("arima.smodel changed!")
#				print(c(as.numeric(text), as.numeric(text2), as.numeric(text3)))
			}
			else{
				status_print("possible wrong input for sarima!")
			}
		}
		
		#arima.ar
		if(editable==entry.arima.ar){
		  text <- entry.arima.ar$GetChars(0, -1)
      params <- cutParam(text)
		  if(unique(isNumber(params,integer=FALSE))[1]==TRUE && length(unique(isNumber(params,integer=FALSE)))==1){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(arima.ar=params),s)})
		    status_print("arima.ar changed!")
		  }
		  else{
		    status_print("possible wrong input for arima.ar!")
		  }
		}
		#arima.ma
		if(editable==entry.arima.ma){
		  text <- entry.arima.ma$GetChars(0, -1)
		  params <- cutParam(text)
		  if(unique(isNumber(params,integer=FALSE))[1]==TRUE && length(unique(isNumber(params,integer=FALSE)))==1){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(arima.ma=params),s)})
		    status_print("arima.ma changed!")
		  }
		  else{
		    status_print("possible wrong input for arima.ma!")
		  }
		}
    
		if(editable==entry.maxorder1||editable==entry.maxorder2){
			text <- entry.maxorder1$GetChars(0, -1)
			text2 <- entry.maxorder2$GetChars(0, -1)
			if(isNumber(text) && isNumber(text2)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(automdl.maxorder=c(as.numeric(text),as.numeric(text2))),s)})
				status_print("automdl.maxorder changed!")
			}
			else{
				status_print("possible wrong input for automdl.maxorder!")
			}
		}
		
		if(editable==entry.maxdiff1 || editable==entry.maxdiff2){
			text <- entry.maxdiff1$GetChars(0, -1)
			text2 <- entry.maxdiff2$GetChars(0, -1)
			if(isNumber(text) && isNumber(text2)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(automdl.maxdiff=c(as.numeric(text),as.numeric(text2))),s)})
				status_print("automdl.maxdiff changed!")
			}
			else{
				status_print("possible wrong input for automdl.maxdiff!")
			}
		}
		
		if(editable==entry.regvariables){
			vars <- sapply(outlierlist, FUN=function(s){
						paste(s[1],s[2],".",s[3], sep="")
					})
			if(trim(entry.regvariables$GetText())=="" && length(vars)!=0){
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=as.vector(vars)),s)})
				status_print("regression.variables changed!")
			}
			else if(trim(entry.regvariables$GetText())!="" && cutParam(entry.regvariables$GetText())[1]!="*" && length(vars)==0){
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=cutParam(entry.regvariables$GetText())),s)})
				status_print("regression.variables changed!")
			}
			else if(trim(entry.regvariables$GetText())=="" && length(vars)==0){
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=NULL),s)})
				status_print("regression.variables changed!")
			}
			else if(cutParam(entry.regvariables$GetText())[1]!="*"){
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.variables=append(cutParam(entry.regvariables$GetText()), as.vector(vars))),s)})
				status_print("regression.variables changed!")
			}
		}
		
		if(editable==entry.reguser){
			if(!isEmpty(text) && trim(text) != "*"){
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.user=(text)),s)})
				status_print("regression.user changed!")
			}
			else{
				status_print("possible wrong input for regression.user!")
			}
		}
		
		if(editable==entry.aictest){
			if(!isEmpty(text) && trim(text) != "*"){
			  regression[indices] <<- TRUE
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.aictest=(text)),s)})
				status_print("regression.aictest changed!")
			}
			else{
				status_print("possible wrong input for regression.aictest!")
			}
		}
    
		#slidingspans.length
		if(editable==entry.slidingspans.length){
		  text <- entry.slidingspans.length$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE) && checkb.slidingspans.length$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.length=as.numeric(text)),s)})
        checkb.slidingspans.length$SetInconsistent(FALSE)
		    status_print("slidingspans.length changed!")
		  }
		  else{
		    status_print("possible wrong input for slidingspans.length!")
		  }
		}
		#slidingspans.numspan
		if(editable==entry.slidingspans.numspans){
		  text <- entry.slidingspans.numspans$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE) && checkb.slidingspans.numspans$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.numspans=as.numeric(text)),s)})
		    checkb.slidingspans.numspans$SetInconsistent(FALSE)
		    status_print("slidingspans.numspans changed!")
		  }
		  else{
		    status_print("possible wrong input for slidingspans.numspans!")
		  }
		}
		#slidingspans.start
		if(editable==entry.slidingspans.start1 || editable==entry.slidingspans.start2){
		  text <- entry.slidingspans.start1$GetChars(0, -1)
		  text2 <- entry.slidingspans.start2$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE)==TRUE && isPeriod(text2)==TRUE && checkb.slidingspans.start$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(slidingspans.start=c(as.numeric(text),text2)),s)})
		    checkb.slidingspans.start$SetInconsistent(FALSE)
		    status_print("slidingspans.start changed!")
		  }
		  else{
		    status_print("possible wrong input for slidingspans.start!")
		  }
		}
    
		#history.sadjlags
		if(editable==entry.history.sadjlags){
		  text <- entry.history.sadjlags$GetChars(0, -1)
		  params <- cutParam(text)
		  if(unique(isNumber(params,integer=FALSE))[1]==TRUE && length(unique(isNumber(params,integer=FALSE)))==1
         && checkb.history.sadjlags$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.sadjlags=as.numeric(params)),s)})
		    checkb.history.sadjlags$SetInconsistent(FALSE)
		    status_print("history.sadjlags changed!")
		  }
		  else{
		    status_print("possible wrong input for history.sadjlags!")
		  }
		}
		#history.trendlags
		if(editable==entry.history.trendlags){
		  text <- entry.history.trendlags$GetChars(0, -1)
		  params <- cutParam(text)
		  if(unique(isNumber(params,integer=FALSE))[1]==TRUE && length(unique(isNumber(params,integer=FALSE)))==1
		     && checkb.history.trendlags$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.trendlags=as.numeric(params)),s)})
		    checkb.history.trendlags$SetInconsistent(FALSE)
		    status_print("history.trendlags changed!")
		  }
		  else{
		    status_print("possible wrong input for history.trendlags!")
		  }
		}
		#history.start
		if(editable==entry.history.startyear || editable==entry.history.startperiod){
		  text <- entry.history.startyear$GetChars(0, -1)
		  text2 <- entry.history.startperiod$GetChars(0, -1)
		  if(isNumber(text, integer=FALSE)==TRUE && isPeriod(text2)==TRUE && checkb.history.start$GetActive()==TRUE){
		    lapply(indices, FUN=function(s){object <<- setP(object,list(history.start=c(as.numeric(text),text2)),s)})
		    checkb.history.start$SetInconsistent(FALSE)
		    status_print("history.start changed!")
		  }
		  else{
		    status_print("possible wrong input for history.start!")
		  }
		}
		
		if(editable==entry.usertype){
      if(trim(text)!=""){
        lapply(indices, FUN=function(s){object <<- setP(object,list(regression.usertype=cutParam(text)),s)})
        status_print("regression.usertype changed!")
      }
			else{
			  lapply(indices, FUN=function(s){object <<- setP(object,list(regression.usertype=NULL),s)})
			  status_print("possible wrong input for regression.usertype!")
			}
		}
		
		if(editable==entry.regfilestartstartyear||editable==entry.regfilestartstartperiod){
			text <- entry.regfilestartstartyear$GetChars(0, -1)
			text2 <- entry.regfilestartstartperiod$GetChars(0, -1)
			if(isNumber(text)&&isPeriod(text2)){
			  regression[indices] <<- TRUE
				lapply(indices, FUN=function(s){object <<- setP(object,list(regression.start=c(text, text2)),s)})
				status_print("regression.start changed!")
			}
			else{
				status_print("possible wrong input for regression.start!")
			}
		}
		
#		if(editable==entry.seatsparameter){
#			if(!isEmpty(text) && cutParam(text)!="*"){
#				lapply(indices, FUN=function(s){object <<- setP(object,list(seatsparameter=cutParam(text)),s)})
#				status_print("seatsparameter changed!")
#			}
#			else{
#				status_print("possible wrong input for seatsparameter!")
#			}
#		}
		
		if(editable==entry.sigmalim1||editable==entry.sigmalim2){
			text <- entry.sigmalim1$GetChars(0, -1)
			text2 <- entry.sigmalim2$GetChars(0, -1)
			if(isNumber(text)&&isNumber(text2)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.sigmalim=c(as.numeric(text), as.numeric(text2))),s)})
				status_print("x11.sigmalim changed!")
			}
			else{
				status_print("possible wrong input for x11.sigmalim!")
			}
		}
		
		if(editable==entry.criticalall){
			if(isNumber(text, integer=FALSE)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.critical=as.numeric(text)),s)})
				status_print("outlier.critical changed!")
			}
			else{
				status_print("possible wrong input for outlier.critical!")
			}
		}
		
		if(editable==entry.criticalAO || editable==entry.criticalLS || editable==entry.criticalTC){
			text <- entry.criticalAO$GetChars(0, -1)
			text2 <- entry.criticalLS$GetChars(0, -1)
			text3 <- entry.criticalTC$GetChars(0, -1)
			criticallist <- list()
			if(isNumber(text)){
				criticallist[["AO"]] <- as.numeric(text)
			}
			if(isNumber(text2)){
				criticallist[["LS"]] <- as.numeric(text2)
			}
			if(isNumber(text3)){
				criticallist[["TC"]] <- as.numeric(text3)
			}
      if(length(criticallist)>0)lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.critical=criticallist),s)})
      else lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.critical=NULL),s)})
		}
				
		if(editable==entry.outlierspan.start1||editable==entry.outlierspan.start2 || editable==entry.outlierspan.end1
       ||editable==entry.outlierspan.end2){
		  startyear <- entry.outlierspan.start1$GetChars(0, -1)
		  startperiod <- entry.outlierspan.start2$GetChars(0, -1)
		  endyear <- entry.outlierspan.end1$GetChars(0, -1)
		  endperiod <- entry.outlierspan.end2$GetChars(0, -1)
		  values <- c(NA,NA,NA,NA)
		  problem = FALSE
		  if(checkb.outlierspan.start$GetActive()){
		    if(isNumber(startyear)&&isPeriod(startperiod)){
		      values[1] <- as.numeric(startyear)
		      values[2] <- startperiod
		    }
		    else{
		      problem = TRUE
		      status_print("possible wrong input for outlier.span!")
		    }
		  }
		  if(checkb.outlierspan.end$GetActive()){
		    if(isNumber(endyear)&&isPeriod(endperiod)){
		      values[3] <- as.numeric(endyear)
		      values[4] <- endperiod
		    }
		    else{
		      problem=TRUE
		      status_print("possible wrong input for outlier.span!")
		    }
		  }
		  #     print(values)
		  if(problem==FALSE){
		    status_print("outlier.span changed!")
		    if(length(unique(values))==1 && is.na(values[1]))lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.span=NULL),s)})
		    else lapply(indices, FUN=function(s){object <<- setP(object,list(outlier.span=values),s)})
		  }
		}
		
		if(editable==entry.forecast_years){
			if(isNumber(text)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(forecast_years=as.numeric(text)),s)})
				status_print("forecast_years changed!")
			}
			else{
				status_print("possible wrong input for forecast_years!")
			}
		}
		
		if(editable==entry.backcast_years){
			if(isNumber(text)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(backcast_years=as.numeric(text)),s)})
				status_print("backcast_years changed!")
			}
			else{
				status_print("possible wrong input for backcast_years!")
			}
		}
		
		if(editable==entry.forecast_conf){
			if(isNumber(text, integer=FALSE)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(forecast_conf=as.numeric(text)),s)})
				status_print("forecast_conf changed!")
			}
			else{
				status_print("possible wrong input for forecast_conf!")
			}
		}
		
		if(editable==entry.seasonalma){
			if(isEmpty(text)==FALSE && cutParam(text)[1] != "*"){
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.seasonalma=as.character(cutParam(text))),s)})
				status_print("x11.seasonalma changed!")
			}
			else{
				status_print("possible wrong input for x11.seasonalma!")
			}
		}
		
		if(editable==entry.trendma){
			if(isNumber(text)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.trendma=as.numeric(text)),s)})
				status_print("x11.trendma changed!")
			}
			else{
				status_print("possible wrong input for x11.trendma!")
			}
		}
		
		if(editable==entry.x11final){
			if(grepl("^(\\s)*(AO|LS|TC|user|all|none)(\\s)*((\\s)*(AO|LS|user|TC|all|none)(\\s)*)*$",text)){
				lapply(indices, FUN=function(s){object <<- setP(object,list(x11.final=cutParam(text)),s)})
				status_print("x11.final changed!")
			}
			else{
				status_print("possible wrong input for x11.final!")
			}
		}
	}
	
	#activates the right plottingcontext(cairoDevice) and plots depending on selected tab
	#at the first run the context for each tab will be created
	update_notebook <- function(page.num=notebook.plot$GetCurrentPage(), file=FALSE, onlyplot=FALSE, reset=TRUE){
		if(page.num==2){
			if(context.plot1==0){
				asCairoDevice(area.plot1)
				context.plot1 <<- dev.cur()
			}
			if(!file)dev.set(context.plot1)
			spect <- "acf"
			if(radiob.rsdacfacf2$GetActive() == TRUE)spect <- "acf2"
			if(radiob.rsdacfpacf$GetActive() == TRUE)spect <- "pacf"
      plotRsdAcfGUI(object@x12List[[indices[1]]], which=spect,main=object@x12List[[indices[1]]]@tsName)
		}
		if(page.num==1){
			if(context.plot2==0){
				asCairoDevice(area.plot2)
				context.plot2 <<- dev.cur()
			}
			if(!file)dev.set(context.plot2)
			spect <- "sa"
			if(radiob.spectraloriginal$GetActive() == TRUE)spect <- "original"
			if(radiob.spectralresiduals$GetActive() == TRUE)spect <- "residuals"
			if(radiob.spectralirregular$GetActive() == TRUE)spect <- "irregular"
      if(spect=="sa")
        main <- "Spectrum of the Seasonally Adjusted Series"
      else if(spect=="original")
        main <- "Spectrum of the Original Series"      
      else if(spect=="irregular")
        main <- "Spectrum of the Irregular"
      else if(spect=="residuals")
        main <- "Spectrum of the RegARIMA Residuals"
      main <- paste(object@x12List[[indices[1]]]@tsName,"-",main)
			plotSpec(object@x12List[[indices[1]]], which=spect,main=main,legend_bty="n",legend_horiz=FALSE)
		}
		if(page.num==3){
			if(context.plot3==0){
				asCairoDevice(area.plot3)
				context.plot3 <<- dev.cur()
			}
			if(!file)dev.set(context.plot3)
        main <- paste(object@x12List[[indices[1]]]@tsName,"-","Seasonal Factors by period and SI Ratios")
			  plotSeasFac(object@x12List[[indices[1]]],main=main,legend_bty="n",legend_horiz=FALSE)
		}
		if(page.num==0){
			if(context.plot4==0){
				asCairoDevice(area.plot4)
				context.plot4 <<- dev.cur()
			}
			if(!file)dev.set(context.plot4)
			make_plot(object)
			t <- times(object@x12List[[indices[1]]])
			if(!onlyplot){
				#set values(amount of seasonal periods since begin of backcast or begin of timeseries if no backcast) for sliders
				slider.plotmin$SetRange(0, calcPeriods(t$original, object@x12List[[indices[1]]]@ts))
				slider.plotmax$SetRange(0, calcPeriods(t$original, object@x12List[[indices[1]]]@ts))
				min <- 0
				max <- calcPeriods(t$original, object@x12List[[indices[1]]]@ts)
				if(!is.null(t$forecast) && !is.null(t$backcast)){
					max <-  calcPeriods(c(t$backcast[1], t$backcast[2], t$forecast[3], t$forecast[4]), object@x12List[[indices[1]]]@ts)
				}
				if(!is.null(t$backcast) && is.null(t$forecast)){
					max <-  calcPeriods(c(t$backcast[1], t$backcast[2], t$original[3], t$original[4]), object@x12List[[indices[1]]]@ts)
				}
				if(is.null(t$backcast) && !is.null(t$forecast)){
					max <-  calcPeriods(c(t$original[1], t$original[2], t$forecast[3], t$forecast[4]), object@x12List[[indices[1]]]@ts)
				}
				slider.plotmin$SetRange(min, max)
				slider.plotmax$SetRange(min, max)
        if(reset == TRUE){
          if(!is.null(t$backcast)){
            slider.plotmin$SetValue(calcPeriods(t$backcast, object@x12List[[indices[1]]]@ts))
          }
          else{
            slider.plotmin$SetValue(0)
          }
          slider.plotmax$SetValue(calcPeriods(t$original, object@x12List[[indices[1]]]@ts))
        }
			} 
		}
		if(page.num==4){
			if(onlyplot==FALSE)make_summary(object, table=FALSE)
		}
		if(page.num==5){
			if(onlyplot==FALSE)make_summary(object, text=FALSE)
		}
#   if(page.num==4){
#     if(context.plot5==0){
#       asCairoDevice(area.plot5)
#       context.plot5 <<- dev.cur()
#     }
#     if(!file)dev.set(context.plot5)
#     make_plot(object)
#   }
	}
	
	#Handler responsible for the contextmenus of the plots
	menuhandler <- function(menuitem, userdata){
		if(menuitem == menuitem.x12update){
      dialog <- gtkMessageDialog(window.main, "destroy-with-parent","warning","none", "Please wait, X12-ARIMA is runnning!")
      sink(tmpfsink)
      tt <- try(object <<- x12(object))
      sink()
      
      runOutput <- readLines(tmpfsink)
      runOutput <- runOutput[runOutput!=" "]
      if(length(grep("ERROR",runOutput))>0){
        dialog$destroy()
        runOutput <- paste(runOutput,collapse="\n")
        dialog <- gtkMessageDialog(window.main, "destroy-with-parent","error","cancel", runOutput)
        if(dialog$run()){
          dialog$destroy()
        }
      }else{
			  update_notebook(reset=FALSE)
			  update_outliertable()
			  make_history()
        dialog$destroy()
      }
    
		}
		
		if(menuitem == menuitem.saveaspdf || menuitem == menuitem.saveaspdfwithoutlier || menuitem == menuitem.expplotaspdf){
			dialog <- gtkFileChooserDialog("Save Plot as PDF", window.main, action="save","gtk-cancel", GtkResponseType["cancel"], 
					"gtk-save", GtkResponseType["ok"])
			dialog$SetDoOverwriteConfirmation(TRUE)
			if(dialog$run()==GtkResponseType["ok"]){
				cur <- dev.cur()
				if(!grepl("(.)*(\\.pdf)",dialog$getFilename())){
					filename <- paste(dialog$getFilename(),".pdf", sep = "")
				}
				else{
					filename <- dialog$getFilename()
				}
				pdf(filename)
				update_notebook(file=TRUE)
				dev.off()
				dev.set(cur)
			}
			
			dialog$Destroy()
		}
		
		if(menuitem == menuitem.expplotaspng){
			dialog <- gtkFileChooserDialog("Save Plot as PNG", window.main, action="save","gtk-cancel", GtkResponseType["cancel"], 
					"gtk-save", GtkResponseType["ok"])
			dialog$SetDoOverwriteConfirmation(TRUE)
			if(dialog$run()==GtkResponseType["ok"]){
				cur <- dev.cur()
				if(!grepl("(.)*(\\.png)",dialog$getFilename())){
					filename <- paste(dialog$getFilename(),".png", sep = "")
				}
				else{
					filename <- dialog$getFilename()
				}
				png(filename, width=800, height=800)
				update_notebook(file=TRUE)
				dev.off()
				dev.set(cur)
			}
			
			dialog$Destroy()
		}
		
		if(menuitem == menuitem.addAO || menuitem == menuitem.addLS || menuitem == menuitem.addTC){
			oc <- convertCoordinates(clickedX)
			typ <- 0
			if(menuitem == menuitem.addAO) typ <- 'AO'
			if(menuitem == menuitem.addLS) typ <- 'LS'
			if(menuitem == menuitem.addTC) typ <- 'TC'
			
			ee <- end(object@x12List[[min(indices)]]@ts)
			ss <- start(object@x12List[[min(indices)]]@ts)
			if((oc$year==ee[1]&&oc$period==ee[2])||(oc$year==ss[1]&&oc$period==ss[2]))
				typ <- 'AO'
			TFexistOut <- FALSE
			if(length(outlierlist)>0){
				for(io in 1:length(outlierlist)){
					if(all(outlierlist[[io]]==c(typ,oc$year,oc$period)))
						TFexistOut <- TRUE
				}
			}
			if(!TFexistOut)
				outlierlist[[length(outlierlist)+1]] <<- c(typ,oc$year,oc$period)
			checkb.regvariablesactive$SetActive(TRUE);
      update_regvariables()
      
		}
		
		if(menuitem == menuitem.expsummarycsv || menuitem == menuitem.expsummaryclipboard){
			if(menuitem == menuitem.expsummarycsv){
				dialog <- gtkFileChooserDialog("Save CSV", window.main, action="save","gtk-cancel", GtkResponseType["cancel"], 
						"gtk-save", GtkResponseType["ok"])
				dialog$SetDoOverwriteConfirmation(TRUE)
				if(dialog$run()==GtkResponseType["ok"]){
					
					if(!grepl("(.)*(\\.csv)",dialog$getFilename())){
						filename <- paste(dialog$getFilename(),".csv", sep = "")
					}
					else{
						filename <- dialog$getFilename()
					}
				}
				nl <- max(sapply(object@x12List, function(x) length(x@x12OldOutput)))
				if(checkb.rsdautocorr$GetActive()==TRUE)
					rsdSummary<-c("acf","pacf","acf2")
				else
					rsdSummary<-FALSE	
				
				write.csv2(getMethod("summary","x12Batch")(object, print=FALSE, oldOutput=nl,
								fullSummary=checkb.fullSummary$GetActive(), spectra.detail=checkb.spectraldetail$GetActive(),
								almostout=checkb.almostout$GetActive(), rsd.autocorr=rsdSummary,
								quality.stat=checkb.quality.stat$GetActive(), likelihood.stat=checkb.likelihoodstat$GetActive(),
								aape=checkb.aape$GetActive(), id.rsdseas=checkb.idrsdseas$GetActive(), slidingspans=checkb.summaryslidingspans$GetActive(), history=checkb.summaryhistory$GetActive(), identify=checkb.summaryidentify$GetActive()),file=filename, quote=FALSE)
				dialog$Destroy()
			}
			
			if(menuitem == menuitem.expsummaryclipboard){
				clipboard <- table.summary$GetClipboard(GDK_SELECTION_CLIPBOARD)
				nl <- max(sapply(object@x12List, function(x) length(x@x12OldOutput)))
				if(checkb.rsdautocorr$GetActive()==TRUE)
					rsdSummary<-c("acf","pacf","acf2")
				else
					rsdSummary<-FALSE	
				
				clipboard$SetText(paste(capture.output(write.csv2(getMethod("summary","x12Batch")(object, print=FALSE, oldOutput=nl,
														fullSummary=checkb.fullSummary$GetActive(), spectra.detail=checkb.spectraldetail$GetActive(),
														almostout=checkb.almostout$GetActive(), rsd.autocorr=rsdSummary,
														quality.stat=checkb.quality.stat$GetActive(), likelihood.stat=checkb.likelihoodstat$GetActive(),
														aape=checkb.aape$GetActive(), id.rsdseas=checkb.idrsdseas$GetActive(), slidingspans=checkb.summaryslidingspans$GetActive(), history=checkb.summaryhistory$GetActive(), identify=checkb.summaryidentify$GetActive()), file="", quote=FALSE)),collapse="\n"))
			}
			
		}
		
		if(menuitem == menuitem.x12savep){
			dialog <- gtkFileChooserDialog("Save x12 Parameter", window.main, action="save","gtk-cancel", GtkResponseType["cancel"], 
					"gtk-save", GtkResponseType["ok"])
			dialog$SetDoOverwriteConfirmation(TRUE)
			if(dialog$run()==GtkResponseType["ok"]){
				
				if(!grepl("(.)*(\\.RData)",dialog$getFilename())){
					filename <- paste(dialog$getFilename(),".RData", sep = "")
				}
				else{
					filename <- dialog$getFilename()
				}
				saveP(object, filename)
			}
			
			dialog$Destroy()
		}
    if(menuitem == menuitem.path){
      dialog <- gtkFileChooserDialog("Open X12 Executable", window.main, "open",
                                     "gtk-cancel", GtkResponseType["cancel"], 
                                     "gtk-open", GtkResponseType["accept"])
      
      if (dialog$run() == GtkResponseType["accept"]) {
        x12path(dialog$getFilename())
      }
      
      dialog$destroy()
      menuitem$SetLabel(paste("x12Path: ",capPath(x12path())))
    }
		if(menuitem == menuitem.x12save){
			dialog <- gtkFileChooserDialog("Save x12 Object", window.main, action="save","gtk-cancel", GtkResponseType["cancel"], 
					"gtk-save", GtkResponseType["ok"])
			dialog$SetDoOverwriteConfirmation(TRUE)
			if(dialog$run()==GtkResponseType["ok"]){
				
				if(!grepl("(.)*(\\.RData)",dialog$getFilename())){
					filename <- paste(dialog$getFilename(),".RData", sep = "")
				}
				else{
					filename <- dialog$getFilename()
				}
				save(object, file=filename)
			}
			dialog$Destroy()
		}
		
		if(menuitem == menuitem.x12load){
			dialog <- gtkFileChooserDialog("Load x12 Object", window.main, action="open","gtk-cancel", GtkResponseType["cancel"], 
					"gtk-save", GtkResponseType["ok"])
			if(dialog$run()==GtkResponseType["ok"]){
				
				if(!grepl("(.)*(\\.RData)",dialog$getFilename())){
					filename <- paste(dialog$getFilename(),".RData", sep = "")
				}
				else{
					filename <- dialog$getFilename()
				}
				gSignalHandlerBlock(table.ts$GetSelection(), handler.tstable)
				x12o <- get(load(filename))
				if(class(x12o)=="x12Single" || class(x12o)=="x12Batch"){
					if(class(x12o)=="x12Single"){
						xl <- new("x12List")
						xl <- list(x12o)
						xb <- new("x12Batch",list(x12o@ts))
						xb@x12List[[1]] <- x12o
						x12o <- x12(xb)
					}else if(class(x12o)=="ts"){
						x12o <- x12(new("x12Batch",list(x12o)))
					}	
					object <<- x12o
					elements <- sapply(object@x12List, function(x) x@tsName)
					table.ts$GetSelection()$UnselectAll()
					indices <<- c(1)
					table.model$Clear()
					sapply(elements,
							function(string) {
								iter <- table.model$Append()$iter
								table.model$Set(iter, 0, string)
							})
					setup_summarytable(object, remove=TRUE)
					read_x12(object, indices)
					update_notebook()
					update_outliertable()
					make_history()
					gSignalHandlerUnblock(table.ts$GetSelection(), handler.tstable)
				}
			}
			dialog$Destroy()
		}
		
		if(menuitem == menuitem.x12loadp){
			dialog <- gtkFileChooserDialog("Load x12 Parameter", window.main, action="open","gtk-cancel", GtkResponseType["cancel"], 
					"gtk-save", GtkResponseType["ok"])
			if(dialog$run()==GtkResponseType["ok"]){
				
				if(!grepl("(.)*(\\.RData)",dialog$getFilename())){
					filename <- paste(dialog$getFilename(),".RData", sep = "")
				}
				else{
					filename <- dialog$getFilename()
				}
				paramlist <- get(load(filename))
				dialog$Destroy()
				if(length(paramlist)>0 && class(paramlist[[1]])=="x12Parameter"){
					panel.loadplist <- gtkTable(rows=length(object@x12List)+1, columns=2)
					panel.loadplist$AttachDefaults(gtkLabel("timeseries"), 0, 1, 0, 1)
					panel.loadplist$AttachDefaults(gtkLabel("x12Parameter"), 1, 2, 0, 1)
					i <- 1
					combolist <- 1
					numparams <- length(paramlist)
					sapply(object@x12List, FUN=function(s){
								combo <- gtkComboBoxNewText()
								combo$SetSizeRequest(-1, 15)
								combo$AppendText(" ")
								sapply(seq(1,numparams,1), FUN<-function(s){combo$AppendText(s)})
								panel.loadplist$AttachDefaults(combo, 1, 2, i, i+1)
								panel.loadplist$AttachDefaults(gtkLabel(s@tsName), 0, 1, i, i+1)
								if(class(combolist)!="list")combolist <<- list(combo)
								else combolist <<- append(combolist, combo)
								i <<- i + 1
							})
					panel.scrolledload <- gtkScrolledWindow()
					panel.scrolledload$SetPolicy("GTK_POLICY_NEVER","GTK_POLICY_ALWAYS")
					panel.scrolledload$AddWithViewport(panel.loadplist)
					panel.window <- gtkVBox()
					button.accept <- gtkButton("Accept")
					button.discard <- gtkButton("Discard")
					panel.buttons <- gtkHBox()
					panel.buttons$PackStart(button.accept)
					panel.buttons$PackStart(button.discard)
					panel.window$PackStart(panel.scrolledload, expand=TRUE, fill=TRUE)
					panel.window$PackStart(panel.buttons, expand=FALSE)
					window <- gtkWindow()
					window$SetModal(TRUE)
					window$SetTitle("load Parameter")
					window$Add(panel.window)
					window$Show()
					gSignalConnect(button.accept, "released", f=function(...){
								i <- sapply(combolist, FUN<-function(s){
											if(is.null(s$GetActiveText())==FALSE)s$GetActiveText()
											else " "
										})
								for(k in 1:length(i)){
									if(isNumber(i[k]))object@x12List[[k]]@x12Parameter <<- paramlist[[as.numeric(i[k])]]
								}
								read_x12(object, indices)
								update_outliertable()
								window$Destroy()
								gtkMainQuit()
							})
					gSignalConnect(button.discard, "released", f=function(...){
								window$Destroy()
								gtkMainQuit()
							})
					gSignalConnect(window, "destroy", f=function(...){gtkMainQuit()})
					gtkMain()
				}	
			}
		}
	}
	
	filebuttonhandler <- function(widget, user_data){
#		print(filebutton.regfile$GetFile()$GetPath())
		lapply(indices, FUN=function(s){object <<- setP(object,list(regression.file=filebutton.regfile$GetFile()$GetPath()),s)})
		status_print("regfile changed!")
	}
	
	#changehandler for the sliders under plot responsible for span
	sliderhandler <- function(range, user_data){
		if(range==slider.plotmin){
			if(slider.plotmin$GetValue()>slider.plotmax$GetValue()){
				slider.plotmin$SetValue(slider.plotmax$GetValue())
			}
		}
		if(range==slider.plotmax){
			if(slider.plotmax$GetValue()<slider.plotmin$GetValue()){
				slider.plotmax$SetValue(slider.plotmin$GetValue())
			}
		}
		make_plot(object)
	}
	
	#responsible for formating the seasonalperiods count of sliders to year.period formation
	sliderformat <- function(scale, value){
		t <- calcSpan(times(object@x12List[[min(indices)]]),object@x12List[[min(indices)]]@ts)
		if(scale == slider.plotmin){
			return(paste(t[1],'.',t[2]))
		}
		else if(scale == slider.plotmax){
			return(paste(t[3],'.',t[4]))
		}
	}
	
	setup_summarytable <- function(x12o, remove=FALSE){
		if(remove==TRUE){
			model.summary$Clear()
			lapply(columns.summary, FUN=function(s){table.summary$RemoveColumn(s)})
			model.summary <<- gtkListStore(rep("character", length(x12o)+3))
			table.summary$SetModel(model.summary)
		}
		i <- 0
		sumnames <- names(getMethod("summary","x12Batch")(x12o, print=FALSE))
		sumnames[1] <- "Value"
		for(s in sumnames){
			renderer <- gtkCellRendererText()
			column <- gtkTreeViewColumn()
			renderer$SetAlignment(0.5, 0.5)
			column$SetTitle(s)
			column$PackStart(renderer)
			column$SetAlignment(0.5)
			column$SetExpand(TRUE)
			column$AddAttribute(renderer, "text", i)
			if(i==0)column$AddAttribute(renderer, "background", length(sumnames))
			else column$AddAttribute(renderer, "background", length(sumnames)+1)
			i <- i + 1
			table.summary$AppendColumn(column)
			if(class(columns.summary)!="list")columns.summary <<- list(column)
			else columns.summary <<- append(columns.summary, column)
		}
	}
	
	make_summary <- function(objectS, text = TRUE, table = TRUE){
		#textform
		nl <- max(sapply(objectS@x12List, function(x) length(x@x12OldOutput)))
#		buffer.summary$SetText(paste(capture.output(getMethod("summary","x12Single")(objectS@x12List[[indices[1]]], oldOutput=nl)), collapse="\n"))
		if(text == TRUE){
			#buffer.summarytotal$SetText(paste(capture.output(getMethod("summary","x12Batch")(object, oldOutput=nl,
			#								fullSummary=checkb.fullSummary$GetActive(), spectra.detail=checkb.spectraldetail$GetActive(),
			#								almostout=checkb.almostout$GetActive(), rsd.autocorr=checkb.rsdautocorr$GetActive(),
			#								quality.stat=checkb.quality.stat$GetActive(), likelihood.stat=checkb.likelihoodstat$GetActive(),
			#								aape=checkb.aape$GetActive(), id.rsdseas=checkb.idrsdseas$GetActive())), collapse="\n"))
		if(checkb.rsdautocorr$GetActive()==TRUE)
			rsdSummary<-c("acf","pacf","acf2")
		else
			rsdSummary<-FALSE	
		
		buffer.summarytotal$SetText(paste(capture.output(getMethod("summary","x12Single")(objectS@x12List[[min(indices)]], oldOutput=nl,
											fullSummary=checkb.fullSummary$GetActive(), spectra.detail=checkb.spectraldetail$GetActive(),
											almostout=checkb.almostout$GetActive(), 
											rsd.autocorr=rsdSummary,
											quality.stat=checkb.quality.stat$GetActive(), likelihood.stat=checkb.likelihoodstat$GetActive(),
											aape=checkb.aape$GetActive(), id.rsdseas=checkb.idrsdseas$GetActive(), slidingspans=checkb.summaryslidingspans$GetActive(), history=checkb.summaryhistory$GetActive(), identify=checkb.summaryidentify$GetActive())), collapse="\n"))
		}
		#tableform
		if(table == TRUE){
			model.summary$Clear()
			if(checkb.rsdautocorr$GetActive()==TRUE)
				rsdSummary<-c("acf","pacf","acf2")
			else
				rsdSummary<-FALSE	
			
			sum <- getMethod("summary","x12Batch")(objectS, print=FALSE, oldOutput=nl,
					fullSummary=checkb.fullSummary$GetActive(), spectra.detail=checkb.spectraldetail$GetActive(),
					almostout=checkb.almostout$GetActive(), rsd.autocorr=rsdSummary,
					quality.stat=checkb.quality.stat$GetActive(), likelihood.stat=checkb.likelihoodstat$GetActive(),
					aape=checkb.aape$GetActive(), id.rsdseas=checkb.idrsdseas$GetActive(), slidingspans=checkb.summaryslidingspans$GetActive(), history=checkb.summaryhistory$GetActive(), identify=checkb.summaryidentify$GetActive())
			iter <- model.summary$Append()$iter
			n <- 0
			nam <- names(sum)
			nam[1] <- ""
			sapply(nam, FUN=function(k){
						model.summary$Set(iter, n, k)
						n <<- n + 1
					})
			model.summary$Set(iter, n, "grey90")
			model.summary$Set(iter, n+1, "grey90")
			for(i in 1:dim(sum)[1]){
				if(grepl(".*OLD OUTPUT.*",sum[i,1])){
					iter <- model.summary$Append()$iter
					n <- 0
					sapply(sum[i,], FUN=function(k){
								model.summary$Set(iter, n, "")
								n <<- n + 1
							})
					model.summary$Set(iter, n, "white")
				}
				iter <- model.summary$Append()$iter
				n <- 0
				sapply(sum[i,], FUN=function(k){
							model.summary$Set(iter, n, k)
							n <<- n + 1
						})
				model.summary$Set(iter, n, "grey90")
				model.summary$Set(iter, n+1, "grey96")
#			model.summary$Set(iter, 0, string[1], 1, string[2], 2, string[3])
			}
		}
	}
#	setupSummarytable <- function(table, x12sum){
#		i <- 0
#		sapply(names(x12sum), FUN=function(s){
#					renderer <- gtkCellRendererText()
#					column <- gtkTreeViewColumn()
#					renderer$SetAlignment(0.5, 0.5)
#					column$SetTitle(s)
#					column$PackStart(renderer)
#					column$SetAlignment(0.5)
#					column$AddAttribute(renderer.manualoutliertype, "text", i)
#					i <<- i + 1
#					table$AppendColumn(column)
#				})
#		table
#	}
	
	#reads the x12-parameters of a x12batch from the selected indices into the gui
	#is also responsible for checking if specific aprameter is equal/different in multiple x12singles 
	#selected is a vector of the ts indices which should be used c(1,3,4) if only ts 1,3,4 should be used
	read_x12 <- function(x12batch, selected){
    #old variable names, kept for the moment for possible debugging
# 		v <- getP(x12batch, list("span", "modelspan", "decimals", "transform", "arima", "sarima", "automdl",
# 						"maxorder", "maxdiff", "regvariables", "reguser", "regfile", "automdl", "balanced", "acceptdefault", 
# 						"usertype", "centeruser", "regfilestart", "seats", "seatsparameter", "sigmalim", "outlier", "critical", "outlier_span", 
# 						"outlier_method", "forecast_years","forecast_conf", "backcast_years", "estimate", "estOutofsample",
# 						"slidingspans", "aictest", "onlytd", "sfshort", "samode", "seasonalma", "trendma", "x11appendfcst", "x11appendbcst",
# 						"x11calendarsigma", "x11excludefcst", "x11final", "x11regress"))
	  v <- getP(x12batch, list("series.span", "series.modelspan", "transform.function", "transform.power", "transform.adjust",
                             "arima.model", "arima.smodel","arima.ar", "arima.ma", "automdl",
                             "check", "check.maxlag",
	                           "automdl.maxorder", "automdl.maxdiff", "regression.variables", "regression.user", 
                             "regression.file", "automdl.balanced", "automdl.acceptdefault", 
                             "identify", "identify.diff", "identify.sdiff", "identify.maxlag",
	                           "regression.usertype", "regression.centeruser", "regression.start", "regression.aictest",
                             #"seats", "seatsparameter", 
							 "x11.sigmalim", "outlier.types", "outlier.critical", "outlier.span", 
	                           "outlier.method", "forecast_years","forecast_conf", "backcast_years", "estimate", "estimate.outofsample",
	                           "slidingspans", "slidingspans.fixmdl","slidingspans.fixreg","slidingspans.length","slidingspans.numspans",
	                           "slidingspans.outlier","slidingspans.additivesa","slidingspans.start",
                             "history", "history.estimates", "history.fixmdl", "history.fixreg", "history.outlier","history.target",
                             "history.sadjlags", "history.trendlags", "history.start",
                             "regression.aictest", "x11.type", "x11.sfshort", "x11.samode", "x11.seasonalma", 
                             "x11.trendma", "x11.appendfcst", "x11.appendbcst",
	                           "x11.calendarsigma", "x11.excludefcst", "x11.final", "x11regression"))
		v <- v[selected]
    #regression[selected] <<- sapply(v[selected],function(x)any(!sapply(x[c("regression.variables", "regression.user", "regression.file")],is.null)))
	  
    vreg <- regression[selected]
	  status_print(cat(vreg))
		####span
		gSignalHandlerBlock(entry.spanstartyear, handler.spanstartyear)
		gSignalHandlerBlock(entry.spanstartperiod, handler.spanstartperiod)
		gSignalHandlerBlock(entry.spanendyear, handler.spanendyear)
		gSignalHandlerBlock(entry.spanendperiod, handler.spanendperiod)
		if(length(unique(lapply(v, FUN=function(s){s$series.span})))==1){
			#span equal in all ts objects
			if(is.null(v[[1]]$series.span)){
				checkb.spanactive$SetInconsistent(FALSE)
				checkb.spanactive$SetActive(FALSE)
				update_toggle(checkb.spanactive)
			}
			else{
				checkb.spanactive$SetInconsistent(FALSE)
				checkb.spanactive$SetActive(TRUE)
				update_toggle(checkb.spanactive)
				#start
				if((is.na(v[[1]]$series.span)[1])||is.na((v[[1]]$series.span)[2])){
					checkb.spanstart$SetActive(FALSE)
					update_toggle(checkb.spanactive)
				}
				else{
					checkb.spanstart$SetActive(TRUE)
					update_toggle(checkb.spanactive)
					entry.spanstartyear$SetText(v[[1]]$series.span[1])
					entry.spanstartperiod$SetText(v[[1]]$series.span[2])
				}
				#end
				if((is.na(v[[1]]$series.span)[3])||is.na((v[[1]]$series.span)[4])){
					checkb.spanend$SetActive(FALSE)
					update_toggle(checkb.spanactive)
				}
				else{
					checkb.spanend$SetActive(TRUE)
					update_toggle(checkb.spanactive)
					entry.spanendyear$SetText(v[[1]]$series.span[3])
					entry.spanendperiod$SetText(v[[1]]$series.span[4])
				}
			}
		}
		else{
			checkb.spanactive$SetInconsistent(TRUE)
			checkb.spanactive$SetActive(TRUE)
			update_toggle(checkb.spanactive)
		}
		gSignalHandlerUnblock(entry.spanstartyear, handler.spanstartyear)
		gSignalHandlerUnblock(entry.spanstartperiod, handler.spanstartperiod)
		gSignalHandlerUnblock(entry.spanendyear, handler.spanendyear)
		gSignalHandlerUnblock(entry.spanendperiod, handler.spanendperiod)
	  
	  ###regression
	  #not really a x12-parameter but a variable inside the gui that works like a
	  #virtual x12-parameter i.e. its at startup FALSE for all time series and can be changed
	  #with the regression checkbox
	  gSignalHandlerBlock(checkb.regressionactive, handler.regressionactive)
	  if(length(unique(vreg))==1){
	    checkb.regressionactive$SetInconsistent(FALSE)
	    if(vreg[[1]]==TRUE){
	      checkb.regressionactive$SetActive(TRUE)
	    }else {
	      checkb.regressionactive$SetActive(FALSE)
	      radiob.regression$SetActive(TRUE) 
	      radiob.x11regression$SetActive(FALSE)
	    }
	    #removes the possibility to add manual outlier if the regression is not active
	    #would lead to some questionable situations
	    toggle(c(menuitem.addAO, menuitem.addLS, menuitem.addTC, 
	             button.manualoutlieraddclick, button.manualoutlieradd), checkb.regressionactive)
	    toggle(c(radiob.regression, radiob.x11regression,
	             checkb.regvariablesactive, 
	             checkb.reguseractive, 
	             checkb.regfileactive, 
	             checkb.usertypeactive, 
	             checkb.centeruseractive, 
	             checkb.regfilestartactive, 
	             checkb.aictestactive), checkb.regressionactive)
	  }
	  else{
	    checkb.regressionactive$SetInconsistent(TRUE)
	    checkb.regressionactive$SetActive(TRUE)
	    toggle(c(menuitem.addAO, menuitem.addLS, menuitem.addTC, 
	             button.manualoutlieraddclick, button.manualoutlieradd), checkb.regressionactive)
	    toggle(c(radiob.regression, radiob.x11regression,
	             checkb.regvariablesactive, 
	             checkb.reguseractive, 
	             checkb.regfileactive, 
	             checkb.usertypeactive, 
	             checkb.centeruseractive, 
	             checkb.regfilestartactive, 
	             checkb.aictestactive), checkb.regressionactive)
	    update_toggle(checkb.regressionactive)
	  }
	  gSignalHandlerUnblock(checkb.regressionactive, handler.regressionactive)
    
	  ###x11regress
	  gSignalHandlerBlock(radiob.x11regression, handler.x11regression)
	  allowAIC <- FALSE
	  if(length(unique(lapply(v, FUN=function(s){s$x11regression})))==1){
	    radiob.x11regression$SetInconsistent(FALSE)
	    radiob.regression$SetInconsistent(FALSE)
	    if(vreg[1]==TRUE){
	      if(v[[1]]$x11regression==TRUE){
	        radiob.x11regression$SetActive(TRUE)
	        radiob.regression$SetActive(FALSE)
	        entry.aictest$SetSensitive(FALSE)
	        checkb.aictestactive$SetSensitive(FALSE)
	      }
	      else{
	        radiob.regression$SetActive(TRUE) 
	        radiob.x11regression$SetActive(FALSE)
	        allowAIC <- TRUE
	        #entry.aictest$SetSensitive(FALSE)
	        checkb.aictestactive$SetSensitive(TRUE)
	      }
	    }
	  }
	  else{
	    radiob.x11regression$SetInconsistent(TRUE)
	    radiob.regression$SetInconsistent(TRUE)
	  }
	  gSignalHandlerUnblock(radiob.x11regression, handler.x11regression)
	  
    
	  ###aictest
	  gSignalHandlerBlock(entry.aictest, handler.aictest)
	  if(length(unique(lapply(v, FUN=function(s){s$regression.aictest})))==1){
	    if(is.null(v[[1]]$regression.aictest)|| vreg[[1]]==FALSE){
	      checkb.aictestactive$SetInconsistent(FALSE)
	      checkb.aictestactive$SetActive(FALSE)
	      entry.aictest$SetText("")
	      entry.aictest$SetSensitive(FALSE)
	      #update_toggle(checkb.aictestactive)
	    }
	    else{
	      if (allowAIC){
	        checkb.aictestactive$SetInconsistent(FALSE)
	        checkb.aictestactive$SetActive(TRUE)
	        entry.aictest$SetSensitive(TRUE)
	        #update_toggle(checkb.aictestactive)
	        entry.aictest$SetText(concParam((v[[1]]$regression.aictest)))
	      }
	    }
	  }
	  else{
	    entry.aictest$SetText("*")
	    checkb.aictestactive$SetActive(TRUE)
	    checkb.aictestactive$SetInconsistent(TRUE)
	    update_toggle(checkb.aictestactive)
	  }
	  gSignalHandlerUnblock(entry.aictest, handler.aictest)
		
		###modelspan
		gSignalHandlerBlock(entry.modelspanstartyear, handler.modelspanstartyear)
		gSignalHandlerBlock(entry.modelspanstartperiod, handler.modelspanstartperiod)
		gSignalHandlerBlock(entry.modelspanendyear, handler.modelspanendyear)
		gSignalHandlerBlock(entry.modelspanendperiod, handler.modelspanendperiod)
		if(length(unique(lapply(v, FUN=function(s){s$series.modelspan})))==1){
			#span equal in all ts objects
			if(is.null(v[[1]]$series.modelspan)){
				checkb.modelspanactive$SetInconsistent(FALSE)
				checkb.modelspanactive$SetActive(FALSE)
				update_toggle(checkb.modelspanactive)
			}
			else{
				checkb.modelspanactive$SetInconsistent(FALSE)
				checkb.modelspanactive$SetActive(TRUE)
				update_toggle(checkb.modelspanactive)
				#start
				if((is.na(v[[1]]$series.modelspan)[1])||is.na((v[[1]]$series.modelspan)[2])){
					checkb.modelspanstart$SetActive(FALSE)
					update_toggle(checkb.modelspanactive)
				}
				else{
					checkb.modelspanstart$SetActive(TRUE)
					update_toggle(checkb.modelspanactive)
					entry.modelspanstartyear$SetText(v[[1]]$series.modelspan[1])
					entry.modelspanstartperiod$SetText(v[[1]]$series.modelspan[2])
				}
				#end
				if((is.na(v[[1]]$series.modelspan)[3])||is.na((v[[1]]$series.modelspan)[4])){
					checkb.modelspanend$SetActive(FALSE)
					update_toggle(checkb.modelspanactive)
				}
				else{
					checkb.modelspanend$SetActive(TRUE)
					update_toggle(checkb.modelspanactive)
					entry.modelspanendyear$SetText(v[[1]]$series.modelspan[3])
					entry.modelspanendperiod$SetText(v[[1]]$series.modelspan[4])
				}
			}
		}
		else{
			checkb.modelspanactive$SetInconsistent(TRUE)
			checkb.modelspanactive$SetActive(TRUE)
			update_toggle(checkb.modelspanactive)
		}
		gSignalHandlerUnblock(entry.modelspanstartyear, handler.modelspanstartyear)
		gSignalHandlerUnblock(entry.modelspanstartperiod, handler.modelspanstartperiod)
		gSignalHandlerUnblock(entry.modelspanendyear, handler.modelspanendyear)
		gSignalHandlerUnblock(entry.modelspanendperiod, handler.modelspanendperiod)
    
    ###series.type
# 	  gSignalHandlerBlock(combobox.series.type, handler.series.type)
# 	  if(length(unique(lapply(v, FUN=function(s){s$series.type})))==1){
# 	    if(is.null(v[[1]]$series.type)){
# 	      checkb.series.type$SetInconsistent(FALSE)
# 	      checkb.series.type$SetActive(FALSE)
# 	      update_toggle(checkb.series.type)
# 	    }
# 	    else{
# 	      checkb.series.type$SetInconsistent(FALSE)
# 	      checkb.series.type$SetActive(TRUE)
# 	      update_toggle(checkb.series.type)
# 	      if((v[[1]]$series.type)=="flow")combobox.series.type$SetActive(0)
# 	      if((v[[1]]$series.type)=="stock")combobox.series.type$SetActive(1)
# 	    }
# 	  }
# 	  else{
# 	    checkb.series.type$SetInconsistent(TRUE)
# 	    combobox.series.type$SetActive(-1)
# 	    checkb.series.type$SetActive(TRUE)
# 	    update_toggle(checkb.series.type)
# 	  }
# 	  gSignalHandlerUnblock(combobox.series.type, handler.series.type)
    
    #transform.power
	  gSignalHandlerBlock(entry.transform.power, handler.transform.power)
	  if(length(unique(lapply(v, FUN=function(s){s$transform.power})))==1){
	    if(is.null(v[[1]]$transform.power)){
	      checkb.transform.power$SetInconsistent(FALSE)
	      checkb.transform.power$SetActive(FALSE)
	      update_toggle(checkb.transform.power)
	    }
	    else{
	      checkb.transform.power$SetActive(TRUE)
	      checkb.transform.power$SetInconsistent(FALSE)
	      update_toggle(checkb.transform.power)
	      entry.transform.power$SetText(v[[1]]$transform.power)
	    }
	  }
	  else{
	    entry.transform.power$SetText("*")
	    checkb.transform.power$SetActive(TRUE)
	    checkb.transform.power$SetInconsistent(TRUE)
	    update_toggle(checkb.transform.power)
	  }
	  gSignalHandlerUnblock(entry.transform.power, handler.transform.power)
    
	  #transform.adjust
	  gSignalHandlerBlock(combobox.transform.adjust, handler.transform.adjust)
	  if(length(unique(lapply(v, FUN=function(s){s$transform.adjust})))==1){
	    if(is.null(v[[1]]$transform.adjust)){
	      checkb.transform.adjust$SetActive(FALSE)
	      update_toggle(checkb.transform.adjust)
	    }
	    else{
	      checkb.transform.adjust$SetActive(TRUE)
	      update_toggle(checkb.transform.adjust)
	      if((v[[1]]$transform.adjust)=="lom")combobox.transform.adjust$SetActive(0)
	      if((v[[1]]$transform.adjust)=="loq")combobox.transform.adjust$SetActive(1)
	      if((v[[1]]$transform.adjust)=="lpyear")combobox.transform.adjust$SetActive(2)
	    }
	  }
	  else{
	    combobox.transform.adjust$SetActive(-1)
	    checkb.transform.adjust$SetActive(TRUE)
	    update_toggle(checkb.transform.adjust)
	  }
	  gSignalHandlerUnblock(combobox.transform.adjust, handler.transform.adjust)
    
#		###decimals
#		#supressing change signal while setting new values
#		gSignalHandlerBlock(entry.decimals, handler.decimals)
#		if(length(unique(lapply(v, FUN=function(s){s$decimals})))==1){
#			if(is.null(v[[1]]$decimals)){
##       update_toggle(checkb.decimalsactive)
#			}
#			else{
##       update_toggle(checkb.decimalsactive)
#				entry.decimals$SetText((v[[1]]$decimals))
#			}
#		}
#		else{
#			entry.decimals$SetText("*")
##     checkb.decimalsactive$SetInconsistent(TRUE)
##     checkb.decimalsactive$SetActive(TRUE)
##     update_toggle(checkb.decimalsactive)
#		}
#		gSignalHandlerUnblock(entry.decimals, handler.decimals)
		
		###transform
		gSignalHandlerBlock(combobox.transform, handler.transform)
		if(length(unique(lapply(v, FUN=function(s){s$transform.function})))==1){
			if(is.null(v[[1]]$transform.function)){
#       checkb.transformactive$SetActive(FALSE)
#       update_toggle(checkb.transformactive)
			}
			else{
#       checkb.transformactive$SetActive(TRUE)
#       update_toggle(checkb.transformactive)
				if((v[[1]]$transform.function)=="auto")combobox.transform$SetActive(0)
				if((v[[1]]$transform.function)=="log")combobox.transform$SetActive(1)
				if((v[[1]]$transform.function)=="none")combobox.transform$SetActive(2)
			}
		}
		else{
			combobox.transform$SetActive(-1)
#     checkb.transformactive$SetActive(TRUE)
#     update_toggle(checkb.transformactive)
		}
		gSignalHandlerUnblock(combobox.transform, handler.transform)
    
	  ###identify
	  if(length(unique(lapply(v, FUN=function(s){s$identify})))==1){
	    checkb.identify$SetInconsistent(FALSE)
	    if(v[[1]]$identify==TRUE)checkb.identify$SetActive(TRUE)
	    else checkb.identify$SetActive(FALSE)
      toggle(c(checkb.identify.diff, checkb.identify.sdiff, checkb.identify.maxlag), checkb.identify)
	  }
	  else{
	    checkb.identify$SetInconsistent(TRUE)
	  }
	  #identify.diff
	  gSignalHandlerBlock(entry.identify.diff, handler.identify.diff)
	  if(length(unique(lapply(v, FUN=function(s){s$identify.diff})))==1){
	    if(is.null(v[[1]]$identify.diff) || v[[1]]$identify==FALSE){
	      checkb.identify.diff$SetInconsistent(FALSE)
	      checkb.identify.diff$SetActive(FALSE)
	      update_toggle(checkb.identify.diff)
	    }
	    else{
	      checkb.identify.diff$SetActive(TRUE)
	      checkb.identify.diff$SetInconsistent(FALSE)
	      update_toggle(checkb.identify.diff)
	      entry.identify.diff$SetText(v[[1]]$identify.diff)
	    }
	  }
	  else{
	    entry.identify.diff$SetText("*")
	    checkb.identify.diff$SetActive(TRUE)
	    checkb.identify.diff$SetInconsistent(TRUE)
	    update_toggle(checkb.identify.diff)
	  }
	  gSignalHandlerUnblock(entry.identify.diff, handler.identify.diff)
	  #identify.sdiff
	  gSignalHandlerBlock(entry.identify.sdiff, handler.identify.sdiff)
	  if(length(unique(lapply(v, FUN=function(s){s$identify.sdiff})))==1){
	    if(is.null(v[[1]]$identify.sdiff)|| v[[1]]$identify==FALSE){
	      checkb.identify.sdiff$SetInconsistent(FALSE)
	      checkb.identify.sdiff$SetActive(FALSE)
	      update_toggle(checkb.identify.sdiff)
	    }
	    else{
	      checkb.identify.sdiff$SetActive(TRUE)
	      checkb.identify.sdiff$SetInconsistent(FALSE)
	      update_toggle(checkb.identify.sdiff)
	      entry.identify.sdiff$SetText(v[[1]]$identify.sdiff)
	    }
	  }
	  else{
	    entry.identify.sdiff$SetText("*")
	    checkb.identify.sdiff$SetActive(TRUE)
	    checkb.identify.sdiff$SetInconsistent(TRUE)
	    update_toggle(checkb.identify.sdiff)
	  }
	  gSignalHandlerUnblock(entry.identify.sdiff, handler.identify.sdiff)
	  #identify.maxlag
	  gSignalHandlerBlock(entry.identify.maxlag, handler.identify.maxlag)
	  if(length(unique(lapply(v, FUN=function(s){s$identify.maxlag})))==1){
	    if(is.null(v[[1]]$identify.maxlag)|| v[[1]]$identify==FALSE){
	      checkb.identify.maxlag$SetInconsistent(FALSE)
	      checkb.identify.maxlag$SetActive(FALSE)
	      update_toggle(checkb.identify.maxlag)
	    }
	    else{
	      checkb.identify.maxlag$SetActive(TRUE)
	      checkb.identify.maxlag$SetInconsistent(FALSE)
	      update_toggle(checkb.identify.maxlag)
	      entry.identify.maxlag$SetText(v[[1]]$identify.maxlag)
	    }
	  }
	  else{
	    entry.identify.maxlag$SetText("*")
	    checkb.identify.maxlag$SetActive(TRUE)
	    checkb.identify.maxlag$SetInconsistent(TRUE)
	    update_toggle(checkb.identify.maxlag)
	  }
	  gSignalHandlerUnblock(entry.identify.maxlag, handler.identify.maxlag)
		
	  ###slidingspans.fixreg
	  if(length(unique(lapply(v, FUN=function(s){s$slidingspans.fixreg})))==1){
	    if(is.null(v[[1]]$slidingspans.fixreg) || v[[1]]$slidingspans==FALSE){
	      checkb.slidingspans.fixreg$SetInconsistent(FALSE)
	      checkb.slidingspans.fixreg$SetActive(FALSE)
	      checkb.slidingspans.fixreg1$SetActive(FALSE)
	      checkb.slidingspans.fixreg2$SetActive(FALSE)
	      checkb.slidingspans.fixreg3$SetActive(FALSE)
	      checkb.slidingspans.fixreg4$SetActive(FALSE)
	      update_toggle(checkb.slidingspans.fixreg)
	    }
	    else{
	      checkb.slidingspans.fixreg$SetInconsistent(FALSE)
	      checkb.slidingspans.fixreg$SetActive(TRUE)
	      update_toggle(checkb.slidingspans.fixreg)
	      if('td' %in% v[[1]]$slidingspans.fixreg)checkb.slidingspans.fixreg1$SetActive(TRUE)
        else checkb.slidingspans.fixreg1$SetActive(FALSE)
	      if('holiday' %in% v[[1]]$slidingspans.fixreg)checkb.slidingspans.fixreg2$SetActive(TRUE)
	      else checkb.slidingspans.fixreg2$SetActive(FALSE)
	      if('outlier' %in% v[[1]]$slidingspans.fixreg)checkb.slidingspans.fixreg3$SetActive(TRUE)
	      else checkb.slidingspans.fixreg3$SetActive(FALSE)
	      if('user' %in% v[[1]]$slidingspans.fixreg)checkb.slidingspans.fixreg4$SetActive(TRUE)
	      else checkb.slidingspans.fixreg4$SetActive(FALSE)
	    }
	  }
	  else{
	    checkb.slidingspans.fixreg$SetInconsistent(TRUE)
	    checkb.slidingspans.fixreg$SetActive(TRUE)
	    update_toggle(checkb.slidingspans.fixreg)
	  }
	  ###slidingspans.fixmdl
	  gSignalHandlerBlock(combobox.slidingspans.fixmdl, handler.slidingspans.fixmdl)
	  if(length(unique(lapply(v, FUN=function(s){s$slidingspans.fixmdl})))==1){
	    if(is.null(v[[1]]$slidingspans.fixmdl) || v[[1]]$slidingspans==FALSE){
	      checkb.slidingspans.fixmdl$SetInconsistent(FALSE)
	      checkb.slidingspans.fixmdl$SetActive(FALSE)
	      combobox.slidingspans.fixmdl$SetActive(-1)
	      update_toggle(checkb.slidingspans.fixmdl)
	    }
	    else{
	      checkb.slidingspans.fixmdl$SetInconsistent(FALSE)
	      checkb.slidingspans.fixmdl$SetActive(TRUE)
	      update_toggle(checkb.slidingspans.fixmdl)
	      if((v[[1]]$slidingspans.fixmdl)=="yes")combobox.slidingspans.fixmdl$SetActive(0)
	      if((v[[1]]$slidingspans.fixmdl)=="no")combobox.slidingspans.fixmdl$SetActive(1)
	      if((v[[1]]$slidingspans.fixmdl)=="clear")combobox.slidingspans.fixmdl$SetActive(2)
	    }
	  }
	  else{
	    checkb.slidingspans.fixmdl$SetInconsistent(TRUE)
	    combobox.slidingspans.fixmdl$SetActive(-1)
	    checkb.slidingspans.fixmdl$SetActive(TRUE)
	    update_toggle(checkb.slidingspans.fixmdl)
	  }
	  gSignalHandlerUnblock(combobox.slidingspans.fixmdl, handler.slidingspans.fixmdl)
	  ###slidingspans.outlier
	  gSignalHandlerBlock(combobox.slidingspans.outlier, handler.slidingspans.outlier)
	  if(length(unique(lapply(v, FUN=function(s){s$slidingspans.outlier})))==1){
	    if(is.null(v[[1]]$slidingspans.outlier) || v[[1]]$slidingspans==FALSE){
	      checkb.slidingspans.outlier$SetInconsistent(FALSE)
	      checkb.slidingspans.outlier$SetActive(FALSE)
        combobox.slidingspans.outlier$SetActive(-1)
	      update_toggle(checkb.slidingspans.outlier)
	    }
	    else{
	      checkb.slidingspans.outlier$SetInconsistent(FALSE)
	      checkb.slidingspans.outlier$SetActive(TRUE)
	      update_toggle(checkb.slidingspans.outlier)
	      if((v[[1]]$slidingspans.outlier)=="keep")combobox.slidingspans.outlier$SetActive(0)
	      if((v[[1]]$slidingspans.outlier)=="remove")combobox.slidingspans.outlier$SetActive(1)
	      if((v[[1]]$slidingspans.outlier)=="yes")combobox.slidingspans.outlier$SetActive(2)
	    }
	  }
	  else{
	    checkb.slidingspans.outlier$SetInconsistent(TRUE)
	    combobox.slidingspans.outlier$SetActive(-1)
	    checkb.slidingspans.outlier$SetActive(TRUE)
	    #update_toggle(checkb.slidingspans.outlier)
	  }
	  gSignalHandlerUnblock(combobox.slidingspans.outlier, handler.slidingspans.outlier)
	  ###slidingspans.additativesa
	  gSignalHandlerBlock(combobox.slidingspans.additivesa, handler.slidingspans.additivesa)
	  if(length(unique(lapply(v, FUN=function(s){s$slidingspans.additivesa})))==1){
	    if(is.null(v[[1]]$slidingspans.additivesa) || v[[1]]$slidingspans==FALSE){
	      checkb.slidingspans.additivesa$SetInconsistent(FALSE)
	      combobox.slidingspans.additivesa$SetActive(-1)
	      checkb.slidingspans.additivesa$SetActive(FALSE)
	      update_toggle(checkb.slidingspans.additivesa)
	    }
	    else{
	      checkb.slidingspans.additivesa$SetInconsistent(FALSE)
	      checkb.slidingspans.additivesa$SetActive(TRUE)
	      update_toggle(checkb.slidingspans.additivesa)
	      if((v[[1]]$slidingspans.additivesa)=="differences")combobox.slidingspans.additivesa$SetActive(0)
	      if((v[[1]]$slidingspans.additivesa)=="percent")combobox.slidingspans.additivesa$SetActive(1)
	    }
	  }
	  else{
	    checkb.slidingspans.additivesa$SetInconsistent(TRUE)
	    combobox.slidingspans.additivesa$SetActive(-1)
	    checkb.slidingspans.additivesa$SetActive(TRUE)
	    update_toggle(checkb.slidingspans.additivesa)
	  }
	  gSignalHandlerUnblock(combobox.slidingspans.additivesa, handler.slidingspans.additivesa)
	  ###slidingspans.length
	  gSignalHandlerBlock(entry.slidingspans.length, handler.slidingspans.length)
	  if(length(unique(lapply(v, FUN=function(s){s$slidingspans.length})))==1){
	    if(is.null(v[[1]]$slidingspans.length) || v[[1]]$slidingspans==FALSE){
	      checkb.slidingspans.length$SetInconsistent(FALSE)
	      entry.slidingspans.length$SetText("")
	      checkb.slidingspans.length$SetActive(FALSE)
	      update_toggle(checkb.slidingspans.length)
	    }
	    else{
	      checkb.slidingspans.length$SetInconsistent(FALSE)
	      checkb.slidingspans.length$SetActive(TRUE)
	      entry.slidingspans.length$SetText(v[[1]]$slidingspans.length)
	      update_toggle(checkb.slidingspans.length)
	    }
	  }
	  else{
	    checkb.slidingspans.length$SetInconsistent(TRUE)
	    entry.slidingspans.length$SetText("*")
	    checkb.slidingspans.length$SetActive(TRUE)
	    update_toggle(checkb.slidingspans.length)
	  }
	  gSignalHandlerUnblock(entry.slidingspans.length, handler.slidingspans.length)
	  ###slidingspans.numspans
	  gSignalHandlerBlock(entry.slidingspans.numspans, handler.slidingspans.numspans)
	  if(length(unique(lapply(v, FUN=function(s){s$slidingspans.numspans})))==1){
	    if(is.null(v[[1]]$slidingspans.numspans) || v[[1]]$slidingspans==FALSE){
	      checkb.slidingspans.numspans$SetInconsistent(FALSE)
	      checkb.slidingspans.numspans$SetActive(FALSE)
	      entry.slidingspans.numspans$SetText("")
	      update_toggle(checkb.slidingspans.numspans)
	    }
	    else{
	      checkb.slidingspans.numspans$SetInconsistent(FALSE)
	      checkb.slidingspans.numspans$SetActive(TRUE)
	      entry.slidingspans.numspans$SetText(v[[1]]$slidingspans.numspans)
	      update_toggle(checkb.slidingspans.numspans)
	    }
	  }
	  else{
	    checkb.slidingspans.numspans$SetInconsistent(TRUE)
	    entry.slidingspans.numspans$SetText("*")
	    checkb.slidingspans.numspans$SetActive(TRUE)
	    update_toggle(checkb.slidingspans.numspans)
	  }
	  gSignalHandlerUnblock(entry.slidingspans.numspans, handler.slidingspans.numspans)
	  ###slidingspans.start
	  gSignalHandlerBlock(entry.slidingspans.start1, handler.slidingspans.start1)
	  gSignalHandlerBlock(entry.slidingspans.start2, handler.slidingspans.start2)
	  if(length(unique(lapply(v, FUN=function(s){s$slidingspans.start})))==1){
	    if(is.null(v[[1]]$slidingspans.start) || v[[1]]$slidingspans==FALSE){
	      checkb.slidingspans.start$SetInconsistent(FALSE)
	      entry.slidingspans.start1$SetText("")
	      entry.slidingspans.start2$SetText("")
	      checkb.slidingspans.start$SetActive(FALSE)
	      update_toggle(checkb.slidingspans.start)
	    }
	    else{
	      checkb.slidingspans.start$SetInconsistent(FALSE)
	      checkb.slidingspans.start$SetActive(TRUE)
	      entry.slidingspans.start1$SetText(v[[1]]$slidingspans.start[1])
	      entry.slidingspans.start2$SetText(v[[1]]$slidingspans.start[2])
	      update_toggle(checkb.slidingspans.start)
	    }
	  }
	  else{
	    checkb.slidingspans.start$SetInconsistent(TRUE)
	    entry.slidingspans.start1$SetText("*")
	    entry.slidingspans.start2$SetText("*")
	    checkb.slidingspans.start$SetActive(TRUE)
	    update_toggle(checkb.slidingspans.start)
	  }
	  gSignalHandlerUnblock(entry.slidingspans.start1, handler.slidingspans.start1)
	  gSignalHandlerUnblock(entry.slidingspans.start2, handler.slidingspans.start2)
    
		###arima
		gSignalHandlerBlock(entry.arima1, handler.arima1)
		gSignalHandlerBlock(entry.arima2, handler.arima2)
		gSignalHandlerBlock(entry.arima3, handler.arima3)
		if(length(unique(lapply(v, FUN=function(s){s$arima.model})))==1){
			if(is.null(v[[1]]$arima.model)){
				checkb.arimaactive$SetInconsistent(FALSE)
				checkb.arimaactive$SetActive(FALSE)
				entry.arima1$SetText("")
				entry.arima2$SetText("")
				entry.arima3$SetText("")
				update_toggle(checkb.arimaactive)
			}
			else{
				checkb.arimaactive$SetInconsistent(FALSE)
				checkb.arimaactive$SetActive(TRUE)
				update_toggle(checkb.arimaactive)
				entry.arima1$SetText(v[[1]]$arima.model[1])
				entry.arima2$SetText(v[[1]]$arima.model[2])
				entry.arima3$SetText(v[[1]]$arima.model[3])
			}
		}
		else{
			entry.arima1$SetText("*")
			entry.arima2$SetText("*")
			entry.arima3$SetText("*")
			checkb.arimaactive$SetActive(TRUE)
			checkb.arimaactive$SetInconsistent(TRUE)
			update_toggle(checkb.arimaactive)
		}
		gSignalHandlerUnblock(entry.arima1, handler.arima1)
		gSignalHandlerUnblock(entry.arima2, handler.arima2)
		gSignalHandlerUnblock(entry.arima3, handler.arima3)
		
	  ###check
	  if(length(unique(lapply(v, FUN=function(s){s$check})))==1){
	    checkb.check$SetInconsistent(FALSE)
	    if(v[[1]]$check==TRUE)checkb.check$SetActive(TRUE)
	    else checkb.check$SetActive(FALSE)
      toggle(c(checkb.check.maxlag), checkb.check)
	  }
	  else{
	    checkb.check$SetInconsistent(TRUE)
	  }
	  #check.maxlag
	  gSignalHandlerBlock(entry.check.maxlag, handler.check.maxlag)
	  if(length(unique(lapply(v, FUN=function(s){s$check.maxlag})))==1){
	    if(is.null(v[[1]]$check.maxlag) || v[[1]]$check==FALSE){
	      checkb.check.maxlag$SetInconsistent(FALSE)
	      checkb.check.maxlag$SetActive(FALSE)
	      update_toggle(checkb.check.maxlag)
	    }
	    else{
	      checkb.check.maxlag$SetActive(TRUE)
	      checkb.check.maxlag$SetInconsistent(FALSE)
	      update_toggle(checkb.check.maxlag)
	      entry.check.maxlag$SetText(v[[1]]$check.maxlag)
	    }
	  }
	  else{
	    entry.check.maxlag$SetText("*")
	    checkb.check.maxlag$SetActive(TRUE)
	    checkb.check.maxlag$SetInconsistent(TRUE)
	    update_toggle(checkb.check.maxlag)
	  }
	  gSignalHandlerUnblock(entry.check.maxlag, handler.check.maxlag)
    
		###sarima
		gSignalHandlerBlock(entry.sarima1, handler.sarima1)
		gSignalHandlerBlock(entry.sarima2, handler.sarima2)
		gSignalHandlerBlock(entry.sarima3, handler.sarima3)
		if(length(unique(lapply(v, FUN=function(s){s$arima.smodel})))==1){
			if(is.null(v[[1]]$arima.smodel)){
				checkb.sarimaactive$SetInconsistent(FALSE)
				checkb.sarimaactive$SetActive(FALSE)
				entry.sarima1$SetText("")
				entry.sarima2$SetText("")
				entry.sarima3$SetText("")
				update_toggle(checkb.sarimaactive)
			}
			else{
				checkb.sarimaactive$SetInconsistent(FALSE)
				checkb.sarimaactive$SetActive(TRUE)
				update_toggle(checkb.sarimaactive)
				entry.sarima1$SetText(v[[1]]$arima.smodel[1])
				entry.sarima2$SetText(v[[1]]$arima.smodel[2])
				entry.sarima3$SetText(v[[1]]$arima.smodel[3])
			}
		}
		else{
			entry.sarima1$SetText("*")
			entry.sarima2$SetText("*")
			entry.sarima3$SetText("*")
			checkb.sarimaactive$SetInconsistent(TRUE)
			checkb.sarimaactive$SetActive(TRUE)
			update_toggle(checkb.sarimaactive)
		}
		gSignalHandlerUnblock(entry.sarima1, handler.sarima1)
		gSignalHandlerUnblock(entry.sarima2, handler.sarima2)
		gSignalHandlerUnblock(entry.sarima3, handler.sarima3)
    
	  #arima.ar
	  gSignalHandlerBlock(entry.arima.ar, handler.arima.ar)
	  if(length(unique(lapply(v, FUN=function(s){s$arima.ar})))==1){
	    if(is.null(v[[1]]$arima.ar)){
	      checkb.arima.ar$SetInconsistent(FALSE)
	      checkb.arima.ar$SetActive(FALSE)
	      update_toggle(checkb.arima.ar)
	    }
	    else{
	      checkb.arima.ar$SetInconsistent(FALSE)
	      checkb.arima.ar$SetActive(TRUE)
	      update_toggle(checkb.arima.ar)
#arima.ar in entry.arima.ar geaendert #alex250713
        entry.arima.ar$SetText(concParam(v[[1]]$arima.ar))
	    }
	  }
	  else{
      checkb.arima.ar$SetInconsistent(TRUE)
	    checkb.arima.ar$SetActive(TRUE)
	    update_toggle(checkb.arima.ar)
	  }
	  gSignalHandlerUnblock(entry.arima.ar, handler.arima.ar)
	  #arima.ma
	  gSignalHandlerBlock(entry.arima.ma, handler.arima.ma)
	  if(length(unique(lapply(v, FUN=function(s){s$arima.ma})))==1){
	    if(is.null(v[[1]]$arima.ma)){
	      checkb.arima.ma$SetInconsistent(FALSE)
	      checkb.arima.ma$SetActive(FALSE)
	      update_toggle(checkb.arima.ma)
	    }
	    else{
	      checkb.arima.ma$SetInconsistent(FALSE)
	      checkb.arima.ma$SetActive(TRUE)
	      update_toggle(checkb.arima.ma)
#arima.ar auf entry arima.ar geaendert #alex250713
	      entry.arima.ar$SetText(concParam(v[[1]]$arima.ma))
	    }
	  }
	  else{
	    checkb.arima.ma$SetInconsistent(TRUE)
	    checkb.arima.ma$SetActive(TRUE)
	    update_toggle(checkb.arima.ma)
	  }
	  gSignalHandlerUnblock(entry.arima.ma, handler.arima.ma)
		
		###maxorder
		gSignalHandlerBlock(entry.maxorder1, handler.maxorder1)
	  gSignalHandlerBlock(entry.maxorder2, handler.maxorder2)
		if(length(unique(lapply(v, FUN=function(s){s$automdl.maxorder})))==1){
			if(is.null(v[[1]]$automdl.maxorder)){
      #checkb.maxorderactive$SetActive(FALSE)
      #update_toggle(checkb.maxorderactive)
			}
			else{
      #checkb.maxorderactive$SetActive(TRUE)
      #update_toggle(checkb.maxorderactive)
				entry.maxorder1$SetText(v[[1]]$automdl.maxorder[1])
				entry.maxorder2$SetText(v[[1]]$automdl.maxorder[2])
			}
		}
		else{
			entry.maxorder1$SetText("*")
			entry.maxorder2$SetText("*")
    #checkb.maxorderactive$SetActive(TRUE)
    #update_toggle(checkb.maxorderactive)
		}
		gSignalHandlerUnblock(entry.maxorder1, handler.maxorder1)
	  gSignalHandlerUnblock(entry.maxorder2, handler.maxorder2)
		
		###maxdiff
		gSignalHandlerBlock(entry.maxdiff1, handler.maxdiff1)
	  gSignalHandlerBlock(entry.maxdiff2, handler.maxdiff2)
		if(length(unique(lapply(v, FUN=function(s){s$automdl.maxdiff})))==1){
			if(is.null(v[[1]]$automdl.maxdiff)){
      #checkb.maxdiffactive$SetActive(FALSE)
      #update_toggle(checkb.maxdiffactive)
			}
			else{
			  #checkb.maxdiffactive$SetActive(TRUE)
			  #update_toggle(checkb.maxdiffactive)
			  entry.maxdiff1$SetText(v[[1]]$automdl.maxdiff[1])
			  entry.maxdiff2$SetText(v[[1]]$automdl.maxdiff[2])
			}
		}
		else{
			entry.maxdiff1$SetText("*")
			entry.maxdiff2$SetText("*")
#     checkb.maxdiffactive$SetActive(TRUE)
#     update_toggle(checkb.maxdiffactive)
		}
		gSignalHandlerUnblock(entry.maxdiff1, handler.maxdiff1)
	  gSignalHandlerUnblock(entry.maxdiff2, handler.maxdiff2)
    
		###regvariables
		gSignalHandlerBlock(entry.regvariables, handler.regvariables)
		if(length(unique(lapply(v, FUN=function(s){s$regression.variables})))==1){
			if(is.null(v[[1]]$regression.variables) || vreg[[1]]==FALSE){
				checkb.regvariablesactive$SetInconsistent(FALSE)
				checkb.regvariablesactive$SetActive(FALSE)
				entry.regvariables$SetText("")
				update_toggle(checkb.regvariablesactive)
				outlierlist <<- list()
			}
			else{
				checkb.regvariablesactive$SetInconsistent(FALSE)
				checkb.regvariablesactive$SetActive(TRUE)
				update_toggle(checkb.regvariablesactive)
				retval <- splitRegvariables(v[[1]]$regression.variables)
        #				entry.regvariables$SetText(concParam(v[[1]]$regvariables))
				entry.regvariables$SetText(concParam(retval$regvariables))
				outlierlist <<- splitOulierstring(retval$outliers)
				update_outliertable()	
			}
		}
		else{
			entry.regvariables$SetText("*")
			outlierlist <<- list()
			checkb.regvariablesactive$SetInconsistent(TRUE)
			checkb.regvariablesactive$SetActive(TRUE)
			update_toggle(checkb.regvariablesactive)
		}
		gSignalHandlerUnblock(entry.regvariables, handler.regvariables)
		
		###reguser
		gSignalHandlerBlock(entry.reguser, handler.reguser)
		if(length(unique(lapply(v, FUN=function(s){s$regression.user})))==1){
			if(is.null(v[[1]]$regression.user)|| vreg[[1]]==FALSE){
				checkb.reguseractive$SetInconsistent(FALSE)
				checkb.reguseractive$SetActive(FALSE)
				entry.reguser$SetText("")
				update_toggle(checkb.reguseractive)
			}
			else{
				checkb.reguseractive$SetInconsistent(FALSE)
				checkb.reguseractive$SetActive(TRUE)
				update_toggle(checkb.reguseractive)
				entry.reguser$SetText((v[[1]]$regression.user))
			}
		}
		else{
			entry.reguser$SetText("*")
			checkb.reguseractive$SetInconsistent(TRUE)
			checkb.reguseractive$SetActive(TRUE)
			update_toggle(checkb.reguseractive)
		}
		gSignalHandlerUnblock(entry.reguser, handler.reguser)
		
		###regfile
		if(length(unique(lapply(v, FUN=function(s){s$regression.file})))==1){
			if(is.null(v[[1]]$regression.file)|| vreg[[1]]==FALSE){
				checkb.regfileactive$SetActive(FALSE)
				update_toggle(checkb.regfileactive)
			}
			else{
				checkb.regfileactive$SetActive(TRUE)
				update_toggle(checkb.regfileactive)
				filebutton.regfile$SetFilename((v[[1]]$regression.file))
			}
		}
		else{
			filebutton.regfile$UnselectAll()
			checkb.regfileactive$SetActive(TRUE)
      checkb.regfileactive$SetInconsistent(TRUE)
			update_toggle(checkb.regfileactive)
		}
		
		###usertype
		gSignalHandlerBlock(entry.usertype, handler.usertype)
		if(length(unique(lapply(v, FUN=function(s){s$regression.usertype})))==1){
			if(is.null(v[[1]]$regression.usertype)|| vreg[[1]]==FALSE){
				checkb.usertypeactive$SetInconsistent(FALSE)
				checkb.usertypeactive$SetActive(FALSE)
				entry.usertype$SetText("")
				update_toggle(checkb.usertypeactive)
			}
			else{
				checkb.usertypeactive$SetInconsistent(FALSE)
				checkb.usertypeactive$SetActive(TRUE)
				update_toggle(checkb.usertypeactive)
				entry.usertype$SetText(concParam((v[[1]]$regression.usertype)))
			}
		}
		else{
			entry.usertype$SetText("*")
			checkb.usertypeactive$SetActive(TRUE)
			checkb.usertypeactive$SetInconsistent(TRUE)
			update_toggle(checkb.usertypeactive)
		}
		gSignalHandlerUnblock(entry.usertype, handler.usertype)
		
		###centeruser
		gSignalHandlerBlock(combobox.centeruser, handler.centeruser)
		if(length(unique(lapply(v, FUN=function(s){s$regression.centeruser})))==1){
			if(is.null(v[[1]]$regression.centeruser)|| vreg[[1]]==FALSE){
				checkb.centeruseractive$SetInconsistent(FALSE)
				checkb.centeruseractive$SetActive(FALSE)
				update_toggle(checkb.centeruseractive)
			}
			else{
				checkb.centeruseractive$SetInconsistent(FALSE)
				checkb.centeruseractive$SetActive(TRUE)
				update_toggle(checkb.centeruseractive)
				if((v[[1]]$regression.centeruser)=="mean")combobox.centeruser$SetActive(0)
				if((v[[1]]$regression.centeruser)=="seasonal")combobox.centeruser$SetActive(1)
			}
		}
		else{
			combobox.centeruser$SetActive(-1)
			checkb.centeruseractive$SetActive(TRUE)
			checkb.centeruseractive$SetInconsistent(TRUE)
			update_toggle(checkb.centeruseractive)
		}
		gSignalHandlerUnblock(combobox.centeruser, handler.centeruser)
		
		###regfilestart
		gSignalHandlerBlock(entry.regfilestartstartyear, handler.regfilestartstartyear)
		gSignalHandlerBlock(entry.regfilestartstartperiod, handler.regfilestartstartperiod)
		if(length(unique(lapply(v, FUN=function(s){s$regression.start})))==1){
			if(is.null(v[[1]]$regression.start)|| vreg[[1]]==FALSE){
				checkb.regfilestartactive$SetInconsistent(FALSE)
				checkb.regfilestartactive$SetActive(FALSE)
				entry.regfilestartstartyear$SetText("")
				entry.regfilestartstartperiod$SetText("")
				update_toggle(checkb.regfilestartactive)
			}
			else{
				checkb.regfilestartactive$SetInconsistent(FALSE)
				checkb.regfilestartactive$SetActive(TRUE)
				label.regfilestartstartyear$SetSensitive(TRUE)
				label.regfilestartstartperiod$SetSensitive(TRUE)
				update_toggle(checkb.regfilestartactive)
				entry.regfilestartstartyear$SetText(v[[1]]$regression.start[1])
				entry.regfilestartstartperiod$SetText(v[[1]]$regression.start[2])
			}
		}
		else{
			entry.regfilestartstartyear$SetText("*")
			entry.regfilestartstartperiod$SetText("*")
			checkb.regfilestartactive$SetActive(TRUE)
			checkb.regfilestartactive$SetInconsistent(TRUE)
			update_toggle(checkb.regfilestartactive)
		}
		gSignalHandlerUnblock(entry.regfilestartstartyear, handler.regfilestartstartyear)
		gSignalHandlerUnblock(entry.regfilestartstartperiod, handler.regfilestartstartperiod)
		
#		###seatsparameter
#		gSignalHandlerBlock(entry.seatsparameter, handler.seatsparameter)
#		if(length(unique(lapply(v, FUN=function(s){s$seatsparameter})))==1){
#			if(is.null(v[[1]]$seatsparameter)){
#				checkb.seatsparameteractive$SetInconsistent(FALSE)
#				checkb.seatsparameteractive$SetActive(FALSE)
#				entry.seatsparameter$SetText("")
#				update_toggle(checkb.seatsparameteractive)
#			}
#			else{
#				checkb.seatsparameteractive$SetInconsistent(FALSE)
#				checkb.seatsparameteractive$SetActive(TRUE)
#				update_toggle(checkb.seatsparameteractive)
#				entry.seatsparameter$SetText((v[[1]]$seatsparameter))
#			}
#		}
#		else{
#			entry.seatsparameter$SetText("*")
#			checkb.seatsparameteractive$SetActive(TRUE)
#			checkb.seatsparameteractive$SetInconsistent(TRUE)
#			update_toggle(checkb.seatsparameteractive)
#		}
#		gSignalHandlerUnblock(entry.seatsparameter, handler.seatsparameter)
		
		###sigmalim
		gSignalHandlerBlock(entry.sigmalim1, handler.sigmalim1)
		gSignalHandlerBlock(entry.sigmalim2, handler.sigmalim2)
		if(length(unique(lapply(v, FUN=function(s){s$x11.sigmalim})))==1){
			if(is.null(v[[1]]$x11.sigmalim)){
				checkb.sigmalimactive$SetInconsistent(FALSE)
				checkb.sigmalimactive$SetActive(FALSE)
				entry.sigmalim1$SetText("")
				entry.sigmalim2$SetText("")
				update_toggle(checkb.sigmalimactive)
			}
			else{
				checkb.sigmalimactive$SetInconsistent(FALSE)
				checkb.sigmalimactive$SetActive(TRUE)
				update_toggle(checkb.sigmalimactive)
				entry.sigmalim1$SetText(v[[1]]$x11.sigmalim[1])
				entry.sigmalim2$SetText(v[[1]]$x11.sigmalim[2])
			}
		}
		else{
			entry.sigmalim1$SetText("*")
			entry.sigmalim2$SetText("*")
			checkb.sigmalimactive$SetActive(TRUE)
			checkb.sigmalimactive$SetInconsistent(TRUE)
			update_toggle(checkb.sigmalimactive)
		}
		gSignalHandlerUnblock(entry.sigmalim1, handler.sigmalim1)
		gSignalHandlerUnblock(entry.sigmalim2, handler.sigmalim2)
		
		###critical
		gSignalHandlerBlock(entry.criticalall, handler.criticalall)
		gSignalHandlerBlock(entry.criticalTC, handler.criticalTC)
		gSignalHandlerBlock(entry.criticalLS, handler.criticalLS)
		gSignalHandlerBlock(entry.criticalAO, handler.criticalAO)
		if(length(unique(lapply(v, FUN=function(s){s$outlier.critical})))==1){
			if(is.null(v[[1]]$outlier.critical)){
				checkb.criticalactive$SetInconsistent(FALSE)
				checkb.criticalactive$SetActive(FALSE)
				entry.criticalAO$SetText("")
				entry.criticalLS$SetText("")
				entry.criticalTC$SetText("")
				entry.criticalall$SetText("")
				update_toggle(checkb.criticalactive)
			}
			else{
				checkb.criticalactive$SetInconsistent(FALSE)
				checkb.criticalactive$SetActive(TRUE)
				update_toggle(checkb.criticalactive)
				if(length(v[[1]]$outlier.critical)==1 && is.null(names(v[[1]]$outlier.critical))){
					entry.criticalall$SetText(v[[1]]$outlier.critical[1])
					radiob.criticalall$SetActive(TRUE)
					radiob.criticalspecific$SetActive(FALSE)
				}
				else{
					radiob.criticalall$SetActive(FALSE)
					radiob.criticalspecific$SetActive(TRUE)
					if(is.null(v[[1]]$outlier.critical[["TC"]])==FALSE)entry.criticalTC$SetText(v[[1]]$outlier.critical["TC"])
					else entry.criticalTC$SetText(" ")
					if(is.null(v[[1]]$outlier.critical[["LS"]])==FALSE)entry.criticalLS$SetText(v[[1]]$outlier.critical["LS"])
					else entry.criticalLS$SetText(" ")
					if(is.null(v[[1]]$outlier.critical[["AO"]])==FALSE)entry.criticalAO$SetText(v[[1]]$outlier.critical["AO"])
					else entry.criticalAO$SetText(" ")
				}
			}
		}
		else{
			entry.criticalAO$SetText("*")
			entry.criticalLS$SetText("*")
			entry.criticalTC$SetText("*")
			entry.criticalall$SetText("*")
			checkb.criticalactive$SetActive(TRUE)
			checkb.criticalactive$SetInconsistent(TRUE)
			update_toggle(checkb.criticalactive)
		}
		gSignalHandlerUnblock(entry.criticalall, handler.criticalall)
		gSignalHandlerUnblock(entry.criticalTC, handler.criticalTC)
		gSignalHandlerUnblock(entry.criticalLS, handler.criticalLS)
		gSignalHandlerUnblock(entry.criticalAO, handler.criticalAO)
		
		###outlier
#		gSignalHandlerBlock(entry.outlier, handler.outlier)
		gSignalHandlerBlock(checkb.outlierall, handler.outlierall)
		gSignalHandlerBlock(checkb.outlierTC, handler.outlierTC)
		gSignalHandlerBlock(checkb.outlierAO, handler.outlierAO)
		gSignalHandlerBlock(checkb.outlierLS, handler.outlierLS)
		if(length(unique(lapply(v, FUN=function(s){s$outlier.types})))==1){
			if(is.null(v[[1]]$outlier.types)){
				checkb.outlieractive$SetInconsistent(FALSE)
				checkb.outlieractive$SetActive(FALSE)
#				entry.outlier$SetText("")
				update_toggle(checkb.outlieractive)
			}
			else{
				checkb.outlieractive$SetInconsistent(FALSE)
				checkb.outlierall$SetInconsistent(FALSE)
				checkb.outlierTC$SetInconsistent(FALSE)
				checkb.outlierAO$SetInconsistent(FALSE)
				checkb.outlierLS$SetInconsistent(FALSE)
				checkb.outlieractive$SetActive(TRUE)
				update_toggle(checkb.outlieractive)
#				entry.outlier$SetText(concParam(v[[1]]$outlier))
				checkb.outlierall$SetActive(FALSE)
				checkb.outlierTC$SetActive(FALSE)
				checkb.outlierAO$SetActive(FALSE)
				checkb.outlierLS$SetActive(FALSE)
				for(i in v[[1]]$outlier.types){
					if(tolower(i)=="all")checkb.outlierall$SetActive(TRUE)
					if(toupper(i)=="TC")checkb.outlierTC$SetActive(TRUE)
					if(toupper(i)=="AO")checkb.outlierAO$SetActive(TRUE)
					if(toupper(i)=="LS")checkb.outlierLS$SetActive(TRUE)
				}
				toggle(c(checkb.outlierTC, checkb.outlierAO, checkb.outlierLS), checkb.outlierall, invert=TRUE)
			}
		}
		else{
			#entry.outlier$SetText("*")
			checkb.outlierall$SetInconsistent(TRUE)
			checkb.outlierTC$SetInconsistent(TRUE)
			checkb.outlierAO$SetInconsistent(TRUE)
			checkb.outlierLS$SetInconsistent(TRUE)
			checkb.outlieractive$SetActive(TRUE)
			checkb.outlieractive$SetInconsistent(TRUE)
			update_toggle(checkb.outlieractive)
		}
		gSignalHandlerUnblock(checkb.outlierall, handler.outlierall)
		gSignalHandlerUnblock(checkb.outlierTC, handler.outlierTC)
		gSignalHandlerUnblock(checkb.outlierAO, handler.outlierAO)
		gSignalHandlerUnblock(checkb.outlierLS, handler.outlierLS)
#		gSignalHandlerUnblock(entry.outlier, handler.outlier)
		
		###outlierspan
		gSignalHandlerBlock(entry.outlierspan.start1, handler.outlierspan.start1)
		gSignalHandlerBlock(entry.outlierspan.start2, handler.outlierspan.start2)
	  gSignalHandlerBlock(entry.outlierspan.end1, handler.outlierspan.end1)
	  gSignalHandlerBlock(entry.outlierspan.end2, handler.outlierspan.end2)
		if(length(unique(lapply(v, FUN=function(s){s$outlier.span})))==1){
			if(is.null(v[[1]]$outlier.span) || is.null(v[[1]]$outlier.types) == TRUE){
				checkb.outlierspanactive$SetInconsistent(FALSE)
				checkb.outlierspanactive$SetActive(FALSE)
				entry.outlierspan.start1$SetText("")
				entry.outlierspan.start2$SetText("")
				entry.outlierspan.end1$SetText("")
				entry.outlierspan.end2$SetText("")
				update_toggle(checkb.outlierspanactive)
			}
			else{
				checkb.outlierspanactive$SetInconsistent(FALSE)
				checkb.outlierspanactive$SetActive(TRUE)
				update_toggle(checkb.outlierspanactive)
				if(!is.na(v[[1]]$outlier.span[1]))entry.outlierspan.start1$SetText(v[[1]]$outlier.span[1])
				if(!is.na(v[[1]]$outlier.span[2]))entry.outlierspan.start2$SetText(v[[1]]$outlier.span[2])
				if(!is.na(v[[1]]$outlier.span[3]))entry.outlierspan.end1$SetText(v[[1]]$outlier.span[3])
				if(!is.na(v[[1]]$outlier.span[4]))entry.outlierspan.end2$SetText(v[[1]]$outlier.span[4])
			}
		}
		else{
		  entry.outlierspan.start1$SetText("*")
		  entry.outlierspan.start1$SetText("*")
			checkb.outlierspanactive$SetActive(TRUE)
			checkb.outlierspanactive$SetInconsistent(TRUE)
			update_toggle(checkb.outlierspanactive)
		}
		gSignalHandlerUnblock(entry.outlierspan.start1, handler.outlierspan.start1)
		gSignalHandlerUnblock(entry.outlierspan.start2, handler.outlierspan.start2)
	  gSignalHandlerUnblock(entry.outlierspan.end1, handler.outlierspan.end1)
	  gSignalHandlerUnblock(entry.outlierspan.end2, handler.outlierspan.end2)
    
		###outliermethod
		gSignalHandlerBlock(combobox.outliermethod, handler.outliermethod)
		if(length(unique(lapply(v, FUN=function(s){s$outlier.method})))==1){
			if(is.null(v[[1]]$outlier.method) || is.null(v[[1]]$outlier.types) == TRUE){
				checkb.outliermethodactive$SetInconsistent(FALSE)
				checkb.outliermethodactive$SetActive(FALSE)
				update_toggle(checkb.outliermethodactive)
			}
			else{
				checkb.outliermethodactive$SetInconsistent(FALSE)
				checkb.outliermethodactive$SetActive(TRUE)
				update_toggle(checkb.outliermethodactive)
				if((v[[1]]$outlier.method)=="addone")combobox.outliermethod$SetActive(0)
				if((v[[1]]$outlier.method)=="addall")combobox.outliermethod$SetActive(1)
			}
		}
		else{
			combobox.outliermethod$SetActive(-1)
			checkb.outliermethodactive$SetActive(TRUE)
			checkb.outliermethodactive$SetInconsistent(TRUE)
			update_toggle(checkb.outliermethodactive)
		}
		gSignalHandlerUnblock(combobox.outliermethod, handler.outliermethod)
		
		###samode
		gSignalHandlerBlock(combobox.samode, handler.samode)
		if(length(unique(lapply(v, FUN=function(s){s$x11.samode})))==1){
			if(is.null(v[[1]]$x11.samode)){
				checkb.samodeactive$SetInconsistent(FALSE)
				checkb.samodeactive$SetActive(FALSE)
				update_toggle(checkb.samodeactive)
			}
			else{
				checkb.samodeactive$SetInconsistent(FALSE)
				checkb.samodeactive$SetActive(TRUE)
				update_toggle(checkb.samodeactive)
				if((v[[1]]$x11.samode)=="mult")combobox.samode$SetActive(0)
				if((v[[1]]$x11.samode)=="add")combobox.samode$SetActive(1)
				if((v[[1]]$x11.samode)=="pseudoadd")combobox.samode$SetActive(2)
				if((v[[1]]$x11.samode)=="logadd")combobox.samode$SetActive(3)
			}
		}
		else{
			combobox.samode$SetActive(-1)
			checkb.samodeactive$SetActive(TRUE)
			checkb.samodeactive$SetInconsistent(TRUE)
			update_toggle(checkb.samodeactive)
		}
		gSignalHandlerUnblock(combobox.samode, handler.samode)
		
		###forecast_years
		gSignalHandlerBlock(entry.forecast_years, handler.forecast_years)
		if(length(unique(lapply(v, FUN=function(s){s$forecast_years})))==1){
			if(is.null(v[[1]]$forecast_years)){
				checkb.forecast_yearsactive$SetInconsistent(FALSE)
				checkb.forecast_yearsactive$SetActive(FALSE)
				entry.forecast_years$SetText("")
				update_toggle(checkb.forecast_yearsactive)
			}
			else{
				checkb.forecast_yearsactive$SetInconsistent(FALSE)
				checkb.forecast_yearsactive$SetActive(TRUE)
				update_toggle(checkb.forecast_yearsactive)
				entry.forecast_years$SetText((v[[1]]$forecast_years))
			}
		}
		else{
			entry.forecast_years$SetText("*")
			checkb.forecast_yearsactive$SetActive(TRUE)
			checkb.forecast_yearsactive$SetInconsistent(TRUE)
			update_toggle(checkb.forecast_yearsactive)
		}
		gSignalHandlerUnblock(entry.forecast_years, handler.forecast_years)
		
		###forecast_conf
		gSignalHandlerBlock(entry.forecast_conf, handler.forecast_conf)
		if(length(unique(lapply(v, FUN=function(s){s$forecast_conf})))==1){
			if(is.null(v[[1]]$forecast_conf)){
#       checkb.forecast_confactive$SetActive(FALSE)
#       update_toggle(checkb.forecast_confactive)
			}
			else{
#       checkb.forecast_confactive$SetActive(TRUE)
#       update_toggle(checkb.forecast_confactive)
				entry.forecast_conf$SetText((v[[1]]$forecast_conf))
			}
		}
		else{
			entry.forecast_conf$SetText("*")
#     checkb.forecast_confactive$SetActive(TRUE)
#     update_toggle(checkb.forecast_confactive)
		}
		gSignalHandlerUnblock(entry.forecast_conf, handler.forecast_conf)
		
		###backcast_years
		gSignalHandlerBlock(entry.backcast_years, handler.backcast_years)
		if(length(unique(lapply(v, FUN=function(s){s$backcast_years})))==1){
			if(is.null(v[[1]]$backcast_years)){
				checkb.backcast_yearsactive$SetInconsistent(FALSE)
				checkb.backcast_yearsactive$SetActive(FALSE)
				entry.backcast_years$SetText("")
				update_toggle(checkb.backcast_yearsactive)
			}
			else{
				checkb.backcast_yearsactive$SetInconsistent(FALSE)
				checkb.backcast_yearsactive$SetActive(TRUE)
				update_toggle(checkb.backcast_yearsactive)
				entry.backcast_years$SetText((v[[1]]$backcast_years))
			}
		}
		else{
			entry.backcast_years$SetText("*")
			checkb.backcast_yearsactive$SetActive(TRUE)
			checkb.backcast_yearsactive$SetInconsistent(TRUE)
			update_toggle(checkb.backcast_yearsactive)
		}
		gSignalHandlerUnblock(entry.backcast_years, handler.backcast_years)
		
	
		###seasonalma
		gSignalHandlerBlock(entry.seasonalma, handler.seasonalma)
		if(length(unique(lapply(v, FUN=function(s){s$x11.seasonalma})))==1){
			if(is.null(v[[1]]$x11.seasonalma)){
				checkb.seasonalmaactive$SetInconsistent(FALSE)
				checkb.seasonalmaactive$SetActive(FALSE)
				entry.seasonalma$SetText("")
				update_toggle(checkb.seasonalmaactive)
			}
			else{
				checkb.seasonalmaactive$SetInconsistent(FALSE)
				checkb.seasonalmaactive$SetActive(TRUE)
				update_toggle(checkb.seasonalmaactive)
				entry.seasonalma$SetText(concParam((v[[1]]$x11.seasonalma)))
			}
		}
		else{
			entry.seasonalma$SetText("*")
			checkb.seasonalmaactive$SetActive(TRUE)
			checkb.seasonalmaactive$SetInconsistent(TRUE)
			update_toggle(checkb.seasonalmaactive)
		}
		gSignalHandlerUnblock(entry.seasonalma, handler.seasonalma)
		
		###trendma
		gSignalHandlerBlock(entry.trendma, handler.trendma)
		if(length(unique(lapply(v, FUN=function(s){s$x11.trendma})))==1){
			if(is.null(v[[1]]$x11.trendma)){
				checkb.trendmaactive$SetInconsistent(FALSE)
				checkb.trendmaactive$SetActive(FALSE)
				entry.trendma$SetText("")
				update_toggle(checkb.trendmaactive)
			}
			else{
				checkb.trendmaactive$SetInconsistent(FALSE)
				checkb.trendmaactive$SetActive(TRUE)
				update_toggle(checkb.trendmaactive)
				entry.trendma$SetText((v[[1]]$x11.trendma))
			}
		}
		else{
			entry.trendma$SetText("*")
			checkb.trendmaactive$SetActive(TRUE)
			checkb.trendmaactive$SetInconsistent(TRUE)
			update_toggle(checkb.trendmaactive)
		}
		gSignalHandlerUnblock(entry.trendma, handler.trendma)
		
		
		###x11calendarsigma
		gSignalHandlerBlock(combobox.x11calendarsigma, handler.x11calendarsigma)
		if(length(unique(lapply(v, FUN=function(s){s$x11.calendarsigma})))==1){
			if(is.null(v[[1]]$x11.calendarsigma)){
				checkb.x11calendarsigmaactive$SetActive(FALSE)
				update_toggle(checkb.x11calendarsigmaactive)
			}
			else{
				checkb.x11calendarsigmaactive$SetActive(TRUE)
				update_toggle(checkb.x11calendarsigmaactive)
				if((v[[1]]$x11.calendarsigma)=="all")combobox.x11calendarsigma$SetActive(0)
				if((v[[1]]$x11.calendarsigma)=="signif")combobox.x11calendarsigma$SetActive(1)
				if((v[[1]]$x11.calendarsigma)=="select")combobox.x11calendarsigma$SetActive(2)
			}
		}
		else{
			combobox.x11calendarsigma$SetActive(-1)
			checkb.x11calendarsigmaactive$SetActive(TRUE)
			update_toggle(checkb.x11calendarsigmaactive)
		}
		gSignalHandlerUnblock(combobox.x11calendarsigma, handler.x11calendarsigma)
		
		###x11final
		gSignalHandlerBlock(entry.x11final, handler.x11final)
		if(length(unique(lapply(v, FUN=function(s){s$x11.final})))==1){
			if(is.null(v[[1]]$x11.final)){
#       checkb.x11finalactive$SetInconsistent(FALSE)
#       checkb.x11finalactive$SetActive(FALSE)
#       update_toggle(checkb.x11finalactive)
			}
			else{
#       checkb.x11finalactive$SetInconsistent(FALSE)
#       checkb.x11finalactive$SetActive(TRUE)
#       update_toggle(checkb.x11finalactive)
				entry.x11final$SetText(concParam((v[[1]]$x11.final)))
			}
		}
		else{
			entry.x11final$SetText("*")
#     checkb.x11finalactive$SetActive(TRUE)
#     checkb.x11finalactive$SetInconsistent(TRUE)
#     update_toggle(checkb.x11finalactive)
		}
		gSignalHandlerUnblock(entry.x11final, handler.x11final)
		
		###automdl
		gSignalHandlerBlock(checkb.automdl, handlercheckb.automdl)
		if(length(unique(lapply(v, FUN=function(s){s$automdl})))==1){
			checkb.automdl$SetInconsistent(FALSE)
			if(v[[1]]$automdl==TRUE){
        checkb.automdl$SetActive(TRUE)
			}
			else {
        checkb.automdl$SetActive(FALSE)
			}
      #checkb.maxorderactive, checkb.maxdiffactive,
			toggle(c(checkb.acceptdefault, checkb.balanced, 
			         entry.maxorder1, entry.maxdiff1,entry.maxorder2, entry.maxdiff2), checkb.automdl)
		}
		else{
			checkb.automdl$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.automdl, handlercheckb.automdl)
    
	  ###history
	  gSignalHandlerBlock(checkb.historyactive, handler.history)
	  if(length(unique(lapply(v, FUN=function(s){s$history})))==1){
	    checkb.historyactive$SetInconsistent(FALSE)
	    if(v[[1]]$history==TRUE)checkb.historyactive$SetActive(TRUE)
	    else checkb.historyactive$SetActive(FALSE)
	    toggle(c(checkb.history.estimates,
	             checkb.history.fixmdl,
	             checkb.history.fixreg,
	             checkb.history.outlier, 
				 checkb.history.target,
	             checkb.history.sadjlags, 
	             checkb.history.trendlags, 
	             checkb.history.start), checkb.historyactive)
	  }
	  else{
	    checkb.historyactive$SetInconsistent(TRUE)
	  }
	  gSignalHandlerUnblock(checkb.historyactive, handler.history)
	  ###history.fixreg
	  if(length(unique(lapply(v, FUN=function(s){s$history.fixreg})))==1){
	    if(is.null(v[[1]]$history.fixreg) || v[[1]]$history==FALSE){
	      checkb.history.fixreg$SetInconsistent(FALSE)
	      checkb.history.fixreg$SetActive(FALSE)
	      checkb.history.fixreg1$SetActive(FALSE)
	      checkb.history.fixreg2$SetActive(FALSE)
	      checkb.history.fixreg3$SetActive(FALSE)
	      checkb.history.fixreg4$SetActive(FALSE)
	      update_toggle(checkb.history.fixreg)
	    }
	    else{
	      checkb.history.fixreg$SetInconsistent(FALSE)
	      checkb.history.fixreg$SetActive(TRUE)
	      update_toggle(checkb.history.fixreg)
	      if('td' %in% v[[1]]$history.fixreg)checkb.history.fixreg1$SetActive(TRUE)
	      else checkb.history.fixreg1$SetActive(FALSE)
	      if('holiday' %in% v[[1]]$history.fixreg)checkb.history.fixreg2$SetActive(TRUE)
	      else checkb.history.fixreg2$SetActive(FALSE)
	      if('outlier' %in% v[[1]]$history.fixreg)checkb.history.fixreg3$SetActive(TRUE)
	      else checkb.history.fixreg3$SetActive(FALSE)
	      if('user' %in% v[[1]]$history.fixreg)checkb.history.fixreg4$SetActive(TRUE)
	      else checkb.history.fixreg4$SetActive(FALSE)
	    }
	  }
	  else{
	    checkb.history.fixreg$SetInconsistent(TRUE)
	    checkb.history.fixreg$SetActive(TRUE)
	    update_toggle(checkb.history.fixreg)
	  }
	  ###history.estimates
	  if(length(unique(lapply(v, FUN=function(s){s$history.estimates})))==1){
	    if(is.null(v[[1]]$history.estimates) || v[[1]]$history==FALSE){
	      checkb.history.estimates$SetInconsistent(FALSE)
	      checkb.history.estimates$SetActive(FALSE)
	      checkb.history.estimatessadj$SetActive(FALSE)
	      checkb.history.estimatessadjchng$SetActive(FALSE)
	      checkb.history.estimatestrend$SetActive(FALSE)
	      checkb.history.estimatestrendchng$SetActive(FALSE)
	      checkb.history.estimatesseasonal$SetActive(FALSE)
	      checkb.history.estimatesfcst$SetActive(FALSE)
	      checkb.history.estimatesaic$SetActive(FALSE)
	      update_toggle(checkb.history.estimates)
	    }
	    else{
	      checkb.history.estimates$SetInconsistent(FALSE)
	      checkb.history.estimates$SetActive(TRUE)
	      update_toggle(checkb.history.estimates)
	      if('sadj' %in% v[[1]]$history.estimates)checkb.history.estimatessadj$SetActive(TRUE)
	      else checkb.history.estimatessadj$SetActive(FALSE)
	      if('sadjchng' %in% v[[1]]$history.estimates)checkb.history.estimatessadjchng$SetActive(TRUE)
	      else checkb.history.estimatessadjchng$SetActive(FALSE)
	      if('trend' %in% v[[1]]$history.estimates)checkb.history.estimatestrend$SetActive(TRUE)
	      else checkb.history.estimatestrend$SetActive(FALSE)
	      if('trendchng' %in% v[[1]]$history.estimates)checkb.history.estimatestrendchng$SetActive(TRUE)
	      else checkb.history.estimatestrendchng$SetActive(FALSE)
	      if('seasonal' %in% v[[1]]$history.estimates)checkb.history.estimatesseasonal$SetActive(TRUE)
	      else checkb.history.estimatesseasonal$SetActive(FALSE)
	      if('fcst' %in% v[[1]]$history.estimates)checkb.history.estimatesfcst$SetActive(TRUE)
	      else checkb.history.estimatesfcst$SetActive(FALSE)
	      if('aic' %in% v[[1]]$history.estimates)checkb.history.estimatesaic$SetActive(TRUE)
	      else checkb.history.estimatesaic$SetActive(FALSE)
	    }
	  }
	  else{
	    checkb.history.estimates$SetInconsistent(TRUE)
	    checkb.history.estimates$SetActive(TRUE)
	    update_toggle(checkb.history.estimates)
	  }
	  ###history.outlier
	  gSignalHandlerBlock(combobox.history.outlier, handler.history.outlier)
	  if(length(unique(lapply(v, FUN=function(s){s$history.outlier})))==1){
	    if(is.null(v[[1]]$history.outlier) || v[[1]]$history==FALSE){
	      checkb.history.outlier$SetInconsistent(FALSE)
	      checkb.history.outlier$SetActive(FALSE)
	      combobox.history.outlier$SetActive(-1)
	      update_toggle(checkb.history.outlier)
	    }
	    else{
	      checkb.history.outlier$SetInconsistent(FALSE)
	      checkb.history.outlier$SetActive(TRUE)
	      update_toggle(checkb.history.outlier)
	      if((v[[1]]$history.outlier)=="keep")combobox.history.outlier$SetActive(0)
	      if((v[[1]]$history.outlier)=="remove")combobox.history.outlier$SetActive(1)
	      if((v[[1]]$history.outlier)=="yes")combobox.history.outlier$SetActive(2)
	    }
	  }
	  else{
	    checkb.history.outlier$SetInconsistent(TRUE)
	    combobox.history.outlier$SetActive(-1)
	    checkb.history.outlier$SetActive(TRUE)
	    #update_toggle(checkb.slidingspans.outlier)
	  }
	  gSignalHandlerUnblock(combobox.history.outlier, handler.history.outlier)
	  ###history.target
	  gSignalHandlerBlock(combobox.history.target, handler.history.target)
	  if(length(unique(lapply(v, FUN=function(s){s$history.target})))==1){
		  if(is.null(v[[1]]$history.target) || v[[1]]$history==FALSE){
			  checkb.history.target$SetInconsistent(FALSE)
			  checkb.history.target$SetActive(FALSE)
			  combobox.history.target$SetActive(-1)
			  update_toggle(checkb.history.target)
		  }
		  else{
			  checkb.history.target$SetInconsistent(FALSE)
			  checkb.history.target$SetActive(TRUE)
			  update_toggle(checkb.history.target)
			  if((v[[1]]$history.target)=="final")combobox.history.target$SetActive(0)
			  if((v[[1]]$history.target)=="concurrent")combobox.history.target$SetActive(1)
#			  if((v[[1]]$history.target)=="yes")combobox.history.target$SetActive(2)
		  }
	  }
	  else{
		  checkb.history.target$SetInconsistent(TRUE)
		  combobox.history.target$SetActive(-1)
		  checkb.history.target$SetActive(TRUE)
		  #update_toggle(checkb.slidingspans.target)
	  }
	  gSignalHandlerUnblock(combobox.history.target, handler.history.target)
	  ###history.sadjlags
	  gSignalHandlerBlock(entry.history.sadjlags, handler.history.sadjlags)
	  if(length(unique(lapply(v, FUN=function(s){s$history.sadjlags})))==1){
	    if(is.null(v[[1]]$history.sadjlags) || v[[1]]$history==FALSE){
	      checkb.history.sadjlags$SetInconsistent(FALSE)
	      entry.history.sadjlags$SetText("")
	      checkb.history.sadjlags$SetActive(FALSE)
	      update_toggle(checkb.history.sadjlags)
	    }
	    else{
	      checkb.history.sadjlags$SetInconsistent(FALSE)
	      checkb.history.sadjlags$SetActive(TRUE)
	      entry.history.sadjlags$SetText(concParam(v[[1]]$history.sadjlags))
	      update_toggle(checkb.history.sadjlags)
	    }
	  }
	  else{
	    checkb.history.sadjlags$SetInconsistent(TRUE)
	    entry.history.sadjlags$SetText("*")
	    checkb.history.sadjlags$SetActive(TRUE)
	    update_toggle(checkb.history.sadjlags)
	  }
	  gSignalHandlerUnblock(entry.history.sadjlags, handler.history.sadjlags)
	  ###history.trendlags
	  gSignalHandlerBlock(entry.history.trendlags, handler.history.trendlags)
	  if(length(unique(lapply(v, FUN=function(s){s$history.trendlags})))==1){
	    if(is.null(v[[1]]$history.trendlags) || v[[1]]$history==FALSE){
	      checkb.history.trendlags$SetInconsistent(FALSE)
	      entry.history.trendlags$SetText("")
	      checkb.history.trendlags$SetActive(FALSE)
	      update_toggle(checkb.history.trendlags)
	    }
	    else{
	      checkb.history.trendlags$SetInconsistent(FALSE)
	      checkb.history.trendlags$SetActive(TRUE)
	      entry.history.trendlags$SetText(concParam(v[[1]]$history.trendlags))
	      update_toggle(checkb.history.trendlags)
	    }
	  }
	  else{
	    checkb.history.trendlags$SetInconsistent(TRUE)
	    entry.history.trendlags$SetText("*")
	    checkb.history.trendlags$SetActive(TRUE)
	    update_toggle(checkb.history.trendlags)
	  }
	  gSignalHandlerUnblock(entry.history.trendlags, handler.history.trendlags)
	  ###history.start
	  gSignalHandlerBlock(entry.history.startyear, handler.history.start1)
	  gSignalHandlerBlock(entry.history.startperiod, handler.history.start2)
	  if(length(unique(lapply(v, FUN=function(s){s$history.start})))==1){
	    if(is.null(v[[1]]$history.start) || v[[1]]$history==FALSE){
	      checkb.history.start$SetInconsistent(FALSE)
	      entry.history.startyear$SetText("")
	      entry.history.startperiod$SetText("")
	      checkb.history.start$SetActive(FALSE)
	      update_toggle(checkb.history.start)
	    }
	    else{
	      checkb.history.start$SetInconsistent(FALSE)
	      checkb.history.start$SetActive(TRUE)
	      entry.history.startyear$SetText(v[[1]]$history.start[1])
	      entry.history.startperiod$SetText(v[[1]]$history.start[2])
	      update_toggle(checkb.history.start)
	    }
	  }
	  else{
	    checkb.history.start$SetInconsistent(TRUE)
	    entry.history.startyear$SetText("*")
	    entry.history.startperiod$SetText("*")
	    checkb.history.start$SetActive(TRUE)
	    update_toggle(checkb.history.start)
	  }
	  gSignalHandlerUnblock(entry.history.startyear, handler.history.start1)
	  gSignalHandlerUnblock(entry.history.startperiod, handler.history.start2)
    
		###balanced
		gSignalHandlerBlock(checkb.balanced, handlercheckb.balanced)
		if(length(unique(lapply(v, FUN=function(s){s$automdl.balanced})))==1){
			checkb.balanced$SetInconsistent(FALSE)
			if(v[[1]]$automdl.balanced==TRUE)checkb.balanced$SetActive(TRUE)
			else checkb.balanced$SetActive(FALSE)
		}
		else{
			checkb.balanced$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.balanced, handlercheckb.balanced)
		
		###acceptdefault
		gSignalHandlerBlock(checkb.acceptdefault, handlercheckb.acceptdefault)
		if(length(unique(lapply(v, FUN=function(s){s$automdl.acceptdefault})))==1){
			checkb.acceptdefault$SetInconsistent(FALSE)
			if(v[[1]]$automdl.acceptdefault==TRUE)checkb.acceptdefault$SetActive(TRUE)
			else checkb.acceptdefault$SetActive(FALSE)
		}
		else{
			checkb.acceptdefault$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.acceptdefault, handlercheckb.acceptdefault)
		
#		###seats
#		gSignalHandlerBlock(checkb.seats, handlercheckb.seats)
#		if(length(unique(lapply(v, FUN=function(s){s$seats})))==1){
#			checkb.seats$SetInconsistent(FALSE)
#			if(v[[1]]$seats==TRUE)checkb.seats$SetActive(TRUE)
#			else checkb.seats$SetActive(FALSE)
#		}
#		else{
#			checkb.seats$SetInconsistent(TRUE)
#		}
#		gSignalHandlerUnblock(checkb.seats, handlercheckb.seats)
		
		###estimate
		gSignalHandlerBlock(checkb.estimate, handlercheckb.estimate)
		if(length(unique(lapply(v, FUN=function(s){s$estimate})))==1){
			checkb.estimate$SetInconsistent(FALSE)
			if(v[[1]]$estimate==TRUE)checkb.estimate$SetActive(TRUE)
			else checkb.estimate$SetActive(FALSE)
      toggle(c(checkb.estOutofsample), checkb.estimate)
		}
		else{
			checkb.estimate$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.estimate, handlercheckb.estimate)
		
		###estOutofsample
		gSignalHandlerBlock(checkb.estOutofsample, handlercheckb.estOutofsample)
		if(length(unique(lapply(v, FUN=function(s){s$estimate.outofsample})))==1){
			checkb.estOutofsample$SetInconsistent(FALSE)
			if(v[[1]]$estimate.outofsample==FALSE){
        checkb.estOutofsample$SetActive(FALSE)
			}
			else checkb.estOutofsample$SetActive(TRUE)
		}
		else{
			checkb.estOutofsample$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.estOutofsample, handlercheckb.estOutofsample)
		
		###slidingspans
		gSignalHandlerBlock(checkb.slidingspans, handlercheckb.slidingspans)
		if(length(unique(lapply(v, FUN=function(s){s$slidingspans})))==1){
			checkb.slidingspans$SetInconsistent(FALSE)
			if(v[[1]]$slidingspans==TRUE){
        checkb.slidingspans$SetActive(TRUE)
			}
			else checkb.slidingspans$SetActive(FALSE)
			toggle(c(checkb.slidingspans.fixmdl, checkb.slidingspans.fixreg, 
			         checkb.slidingspans.length, checkb.slidingspans.numspans, 
			         checkb.slidingspans.outlier, checkb.slidingspans.additivesa,
			         checkb.slidingspans.start), checkb.slidingspans)
		}
		else{
			checkb.slidingspans$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.slidingspans, handlercheckb.slidingspans)
		
		###onlytd
# 		gSignalHandlerBlock(checkb.onlytd, handlercheckb.onlytd)
# 		if(length(unique(lapply(v, FUN=function(s){s$onlytd})))==1){
# 			checkb.onlytd$SetInconsistent(FALSE)
# 			if(v[[1]]$onlytd==TRUE)checkb.onlytd$SetActive(TRUE)
# 			else checkb.onlytd$SetActive(FALSE)
# 		}
# 		else{
# 			checkb.onlytd$SetInconsistent(TRUE)
# 		}
# 		gSignalHandlerUnblock(checkb.onlytd, handlercheckb.onlytd)
		
		###sfshort
		gSignalHandlerBlock(checkb.sfshort, handlercheckb.sfshort)
		if(length(unique(lapply(v, FUN=function(s){s$x11.sfshort})))==1){
			checkb.sfshort$SetInconsistent(FALSE)
			if(v[[1]]$x11.sfshort==TRUE)checkb.sfshort$SetActive(TRUE)
			else checkb.sfshort$SetActive(FALSE)
		}
		else{
			checkb.sfshort$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.sfshort, handlercheckb.sfshort)
    
	  ###x11.type
	  gSignalHandlerBlock(combobox.x11.type, handler.x11.type)
	  if(length(unique(lapply(v, FUN=function(s){s$x11.type})))==1){
	    if(is.null(v[[1]]$x11.type)){
	      checkb.x11.type$SetInconsistent(FALSE)
	      checkb.x11.type$SetActive(FALSE)
	      combobox.x11.type$SetActive(-1)
	      update_toggle(checkb.x11.type)
	    }
	    else{
	      checkb.x11.type$SetInconsistent(FALSE)
	      checkb.x11.type$SetActive(TRUE)
	      update_toggle(checkb.x11.type)
	      if((v[[1]]$x11.type)=="summary")combobox.x11.type$SetActive(0)
	      if((v[[1]]$x11.type)=="trend")combobox.x11.type$SetActive(1)
	      if((v[[1]]$x11.type)=="sa")combobox.x11.type$SetActive(2)
	    }
	  }
	  else{
	    checkb.x11.type$SetInconsistent(TRUE)
	    combobox.x11.type$SetActive(-1)
	    checkb.x11.type$SetActive(TRUE)
	    #update_toggle(checkb.slidingspans.outlier)
	  }
	  gSignalHandlerUnblock(combobox.x11.type, handler.x11.type)
	  ###x11.final
	  if(length(unique(lapply(v, FUN=function(s){s$x11.final})))==1){
	    if(is.null(v[[1]]$x11.final)){
	      checkb.x11.finalactive$SetInconsistent(FALSE)
	      checkb.x11.finalactive$SetActive(FALSE)
	      checkb.x11.finalAO$SetActive(FALSE)
	      checkb.x11.finalLS$SetActive(FALSE)
	      checkb.x11.finalTC$SetActive(FALSE)
	      checkb.x11.finalnone$SetActive(FALSE)
	      checkb.x11.finaluser$SetActive(FALSE)
	      update_toggle(checkb.x11.finalactive)
	    }
	    else{
	      checkb.x11.finalactive$SetInconsistent(FALSE)
	      checkb.x11.finalactive$SetActive(TRUE)
	      update_toggle(checkb.x11.finalactive)
	      if('AO' %in% v[[1]]$x11.final)checkb.x11.finalAO$SetActive(TRUE)
	      else checkb.x11.finalAO$SetActive(FALSE)
	      if('LS' %in% v[[1]]$x11.finalLS)checkb.x11.finalLS$SetActive(TRUE)
	      else checkb.x11.finalLS$SetActive(FALSE)
	      if('TC' %in% v[[1]]$x11.final)checkb.x11.finalTC$SetActive(TRUE)
	      else checkb.x11.finalTC$SetActive(FALSE)
	      if('user' %in% v[[1]]$x11.final)checkb.x11.finaluser$SetActive(TRUE)
	      else checkb.x11.finaluser$SetActive(FALSE)
	      if('none' %in% v[[1]]$x11.final)checkb.x11.finalnone$SetActive(TRUE)
	      else checkb.x11.finalnone$SetActive(FALSE)
	    }
	  }
	  else{
	    checkb.x11.finalactive$SetInconsistent(TRUE)
	    checkb.x11.finalactive$SetActive(TRUE)
	    update_toggle(checkb.x11.finalactive)
	  }
		###x11appendfcst
		gSignalHandlerBlock(checkb.x11appendfcst, handlercheckb.x11appendfcst)
		if(length(unique(lapply(v, FUN=function(s){s$x11.appendfcst})))==1){
			checkb.x11appendfcst$SetInconsistent(FALSE)
			if(v[[1]]$x11.appendfcst==TRUE)checkb.x11appendfcst$SetActive(TRUE)
			else checkb.x11appendfcst$SetActive(FALSE)
		}
		else{
			checkb.x11appendfcst$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.x11appendfcst, handlercheckb.x11appendfcst)
		
		###x11appendbcst
		gSignalHandlerBlock(checkb.x11appendbcst, handlercheckb.x11appendbcst)
		if(length(unique(lapply(v, FUN=function(s){s$x11.appendbcst})))==1){
			checkb.x11appendbcst$SetInconsistent(FALSE)
			if(v[[1]]$x11.appendbcst==TRUE)checkb.x11appendbcst$SetActive(TRUE)
			else checkb.x11appendbcst$SetActive(FALSE)
		}
		else{
			checkb.x11appendbcst$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.x11appendbcst, handlercheckb.x11appendbcst)
		
		###x11excludefcst
		gSignalHandlerBlock(checkb.x11excludefcst, handlercheckb.x11excludefcst)
		if(length(unique(lapply(v, FUN=function(s){s$x11.excludefcst})))==1){
			checkb.x11excludefcst$SetInconsistent(FALSE)
			if(v[[1]]$x11.excludefcst==TRUE)checkb.x11excludefcst$SetActive(TRUE)
			else checkb.x11excludefcst$SetActive(FALSE)
		}
		else{
			checkb.x11excludefcst$SetInconsistent(TRUE)
		}
		gSignalHandlerUnblock(checkb.x11excludefcst, handlercheckb.x11excludefcst)
		
	  
		
		make_history()
	}
##########end read_x12
	#updates the history combobox to values correspondig the old outputs
	make_history <- function(){
		###History Combobox
		for(i in 1:(count.history+1))combobox.history$RemoveText(0)
		k <- length(object@x12List[[min(indices)]]@x12OldParameter)
		if(k>0){
			combobox.history$AppendText("previous")
			count.history <<- 0
			for(i in k:1){
				combobox.history$AppendText(i)
				count.history <<- count.history + 1
			}
			combobox.history$SetActive(0)
		}
	}
	
# make_plotFbcast <- function(...){
#   s <- capture.output(plotFbcast(x12@x12List[[indices[1]]], backcast=checkb.backcast$GetActive(), forecast=checkb.forecast$GetActive(),
#       showCI=checkb.showCI$GetActive(), log_transform=checkb.logtransform_fb$GetActive(), 
#       showLine=checkb.showLine$GetActive(), 
#       points_original=checkb.pointsOriginal$GetActive()))
#   if(is.character(s) & length(s)>0){
#     status_print(s)
#   }
# }
	
	#draws plot depending on parameters set in gui
	make_plot <- function(objectP,...){
#		print(times(x12@x12List[[indices[1]]]))
#		print(calcSpan(times(x12@x12List[[indices[1]]]),x12@x12List[[indices[1]]]@ts))
		showallout = FALSE;
		capture.output(v <- getP(objectP, list("regression.variables")))
		v <- v[indices]
		showout <- NULL
		if(checkb.showout$GetActive() == TRUE){
			if(isNumber(entry.showoutyear$GetText()) == TRUE && isPeriod(entry.showoutperiod$GetText()) == TRUE &&
					as.Period(entry.showoutperiod$GetText()) <= frequency(objectP@x12List[[min(indices)]]@ts) &&
					as.numeric(entry.showoutyear$GetText()) >= start(objectP@x12List[[min(indices)]]@ts) &&
					as.numeric(entry.showoutyear$GetText()) <= end(objectP@x12List[[min(indices)]]@ts)){
				showout <- paste(entry.showoutyear$GetText(),".",as.Period(entry.showoutperiod$GetText()),sep="")
			}
		}
#		if(length(v[[1]])>0 && checkb.showAllout$GetActive()==TRUE)showallout=TRUE;
		s <- capture.output(plotgui(objectP@x12List[[min(indices)]]@x12Output, original=checkb.orig$GetActive(),sa=checkb.sa$GetActive(), trend=checkb.trend$GetActive(), 
						log_transform=checkb.logtransform$GetActive(),
						showCI=checkb.showCI$GetActive(), 
						points_original=checkb.pointsOriginal$GetActive(),showAllout=checkb.showAllout$GetActive(),
						span=calcSpan(times(objectP@x12List[[min(indices)]]),objectP@x12List[[min(indices)]]@ts),showOut=showout,tsName=objectP@x12List[[min(indices)]]@tsName                
                ))
		if(is.character(s) & length(s)>0){
			status_print(s)
		}
	}
	
	#calculates a vector of form c(startyear, startperiod, endyear, endperiod) from
	#the timeselement of the x12Object, a timeseries(for the frequency) and the values of the sliders
	#in the gui
	calcSpan <- function(t, tss){
#		startdate <- c(floor(slider.plotmin$getValue()/(frequency(tss)+1)), slider.plotmin$getValue()%%(frequency(tss)+1))
#		enddate <- c(floor(slider.plotmax$getValue()/(frequency(tss)+1)), slider.plotmax$getValue()%%(frequency(tss)+1))
#		print(enddate)
#		if(!is.null(t$backcast)){
#			return(c((t$backcast[1]+startdate[1])+floor(((t$backcast[2]+startdate[2])/(frequency(tss)+1))),
#							((t$backcast[2]+startdate[2])%%(frequency(tss)+1)),
#							(t$backcast[1]+enddate[1])+floor(((t$backcast[2]+enddate[2])/(frequency(tss)+1))),
#							((t$backcast[2]+enddate[2])%%(frequency(tss)+1))+1))
#		}
#		else{
#			return(c((t$original[1]+startdate[1])+floor(((t$original[2]+startdate[2])/(frequency(tss)+1))),
#							((t$original[2]+startdate[2])%%(frequency(tss)+1)+1),
#							(t$original[1]+enddate[1])+floor(((t$original[2]+enddate[2])/(frequency(tss)+1))),
#							((t$original[2]+enddate[2])%%(frequency(tss)+1))+1))
#		}
		min <- slider.plotmin$getValue()
		max <- slider.plotmax$getValue()
		if(!is.null(t$backcast)){
			return(c((t$backcast[1]+floor(min/frequency(tss))),
							((t$backcast[2]-1+min)%%frequency(tss)+1),
							(t$backcast[1]+floor(max/frequency(tss))),
							((t$backcast[2]-1+max)%%frequency(tss)+1)))
		}
		else{
			return(c((t$original[1]+floor(min/frequency(tss))),
							((t$original[2]-1+min)%%frequency(tss)+1),
							(t$original[1]+floor(max/frequency(tss))),
							((t$original[2]-1+max)%%frequency(tss)+1)))
		}
	}
	
	#helper methode to make printing to statusbar more readable
	status_print <- function(s,p="",...){
		statusbar$Push(statusbar$GetContextId(p), s)
	}
	
	####################################################
	#END FUNCTIONS
	####################################################
	
	window.main$Resize(1100,700)
	
	#Table for timeseries
	renderer.ts$SetAlignment(0.5,0.5)
	column.ts$SetTitle("timeseries")
	column.ts$PackStart(renderer.ts)
	column.ts$SetAlignment(0.5)
	column.ts$AddAttribute(renderer.ts, "text", 0)
	table.ts$AppendColumn(column.ts)
	elements <- sapply(object@x12List, function(x) x@tsName)
	sapply(elements,
			function(string) {
				## Add a new row to the model
				iter <- table.model$Append()$iter
				table.model$Set(iter, 0, string)
			})
	table.ts$SetModel(table.model)
	table.ts$GetSelection()$SetMode("GTK_SELECTION_MULTIPLE")
	handler.tstable <- gSignalConnect(table.ts$GetSelection(), "changed", f=tablehandler)
	
	#history frame
	panel.history$PackStart(label.history, padding=3)
	panel.history$PackStart(combobox.history, padding=3)
	panel.history$PackStart(button.revert, padding=3)
	panel.history$PackStart(button.clearhistory, padding=3)
	frame.history$Add(panel.history)
	gSignalConnect(button.clearhistory, "released", f=function(...){
				dialog <- gtkMessageDialog(window.main, "destroy-with-parent",
						"warning", "yes-no", "History will be lost, do you really want to do this?")
				res <- dialog$Run()
				if(res==-8){
					object@x12List[[indices[1]]] <<- cleanArchive(object@x12List[[min(indices)]])
					make_history()
					update_notebook()
				}
				dialog$destroy()
			})
	gSignalConnect(button.revert, "released", f=function(...){
				dialog <- gtkMessageDialog(window.main, "destroy-with-parent",
						"warning", "yes-no", "Do you really want to revert the parameter?")
				res <- dialog$Run()
				if(res==-8){
					if(isNumber(combobox.history$GetActiveText())==TRUE) k <- as.numeric(combobox.history$GetActiveText())
					else if(combobox.history$GetActiveText()=="previous") k <- length(object@x12List[[min(indices)]]@x12OldParameter)
					if(k>0){
						object@x12List[[indices[1]]] <<- prev(object@x12List[[min(indices)]], n=k)
						make_history()
						update_notebook()
						status_print(capture.output(read_x12(object, indices)))	
						update_outliertable()
					}
					else{
						status_print("Nothing to revert!")
					}
				}
				dialog$destroy()
			})
	
	#######################################################
	#Column timeseries
	#######################################################
	gSignalConnect(button.update, "released", f=update_handler)
#	panel.ts$PackStart(button.update,expand=FALSE)
  panel.ts.scrolled <- gtkScrolledWindow()
  
  panel.ts.scrolled$SetPolicy("GTK_POLICY_NEVER","GTK_POLICY_ALWAYS")
  panel.ts.scrolled$AddWithViewport(table.ts)
  panel.ts$PackStart(panel.ts.scrolled)
  
	#panel.ts$PackStart(table.ts,expand=TRUE)
	panel.ts$PackStart(frame.history, expand=FALSE)
	
	#######################################################
	#Column x12-parameters
	#######################################################
  
	##span
	#start
	panel.series$AttachDefaults(checkb.spanactive, 0, 3, 0, 1)
	handler.spanactive <- gSignalConnect(checkb.spanactive, "toggled", f=checkbuttonhandler)
	panel.series$AttachDefaults(checkb.spanstart, 1, 3, 1, 2)
	gSignalConnect(checkb.spanstart, "toggled", f=checkbuttonhandler)
	panel.spanstart$PackStart(label.spanstartyear)
	entry.spanstartyear$SetSizeRequest(30,-1)
	panel.spanstart$PackStart(entry.spanstartyear, expand=TRUE, padding=2)
	panel.spanstart$PackStart(label.spanstartperiod)
	entry.spanstartperiod$SetSizeRequest(30,-1)
	panel.spanstart$PackStart(entry.spanstartperiod, expand=TRUE, padding=2)
	panel.series$AttachDefaults(panel.spanstart, 1, 3, 2, 3)
	handler.spanstartyear <- gSignalConnect(entry.spanstartyear, "changed", f=x12input_handler)
	handler.spanstartperiod <- gSignalConnect(entry.spanstartperiod, "changed", f=x12input_handler)
	
	#end
	panel.series$AttachDefaults(checkb.spanend, 1, 3, 3, 4)
	gSignalConnect(checkb.spanend, "toggled", f=checkbuttonhandler)
	panel.spanend$PackStart(label.spanendyear)
	entry.spanendyear$SetSizeRequest(30,-1)
	panel.spanend$PackStart(entry.spanendyear, expand=TRUE, padding=2)
	panel.spanend$PackStart(label.spanendperiod)
	entry.spanendperiod$SetSizeRequest(30,-1)
	panel.spanend$PackStart(entry.spanendperiod, expand=TRUE, padding=2)
	panel.series$AttachDefaults(panel.spanend, 1, 3, 4, 5)
	handler.spanendyear <- gSignalConnect(entry.spanendyear, "changed", f=x12input_handler)
	handler.spanendperiod <- gSignalConnect(entry.spanendperiod, "changed", f=x12input_handler)
	
	##modelspan
	#start
	panel.series$AttachDefaults(checkb.modelspanactive, 0, 3, 5, 6)
	gSignalConnect(checkb.modelspanactive, "toggled", f=checkbuttonhandler)
	panel.series$AttachDefaults(checkb.modelspanstart, 1, 3, 7, 8)
	gSignalConnect(checkb.modelspanstart, "toggled", f=checkbuttonhandler)
	panel.modelspanstart$PackStart(label.modelspanstartyear)
	entry.modelspanstartyear$SetSizeRequest(30,-1)
	panel.modelspanstart$PackStart(entry.modelspanstartyear, expand=TRUE, padding=2)
	panel.modelspanstart$PackStart(label.modelspanstartperiod)
	entry.modelspanstartperiod$SetSizeRequest(30,-1)
	panel.modelspanstart$PackStart(entry.modelspanstartperiod, expand=TRUE, padding=2)
	panel.series$AttachDefaults(panel.modelspanstart, 1, 3, 8 ,9)
	handler.modelspanstartyear <- gSignalConnect(entry.modelspanstartyear, "changed", f=x12input_handler)
	handler.modelspanstartperiod <- gSignalConnect(entry.modelspanstartperiod, "changed", f=x12input_handler)
	
	#end
	panel.series$AttachDefaults(checkb.modelspanend, 1, 3, 9, 10)
	gSignalConnect(checkb.modelspanend, "toggled", f=checkbuttonhandler)
	panel.modelspanend$PackStart(label.modelspanendyear)
	entry.modelspanendyear$SetSizeRequest(30,-1)
	panel.modelspanend$PackStart(entry.modelspanendyear, expand=TRUE, padding=2)
	panel.modelspanend$PackStart(label.modelspanendperiod)
	entry.modelspanendperiod$SetSizeRequest(30,-1)
	panel.modelspanend$PackStart(entry.modelspanendperiod, expand=TRUE, padding=2)
	panel.series$AttachDefaults(panel.modelspanend, 1, 3, 10, 11)
	handler.modelspanendyear <- gSignalConnect(entry.modelspanendyear, "changed", f=x12input_handler)
	handler.modelspanendperiod <- gSignalConnect(entry.modelspanendperiod, "changed", f=x12input_handler)
	
  #type
  #panel.series$AttachDefaults(checkb.series.type, 0, 1, 11, 12)
#   gSignalConnect(checkb.series.type, "toggled", f=checkbuttonhandler)
#   combobox.series.type$AppendText("flow")
# 	combobox.series.type$AppendText("stock")
#   panel.series$SetRowSpacing(10, 5)
#   panel.series$AttachDefaults(combobox.series.type, 2, 3, 11, 12)
# 	combobox.series.type$SetSizeRequest(30,-1)
# 	handler.series.type <- gSignalConnect(combobox.series.type, "changed", f=comboboxx12handler)
  
	frame.series$Add(panel.series)
	panel.params$PackStart(frame.series, expand=FALSE)
	
#	#####decimals
## panel.decimals$PackStart(checkb.decimalsactive)
## gSignalConnect(checkb.decimalsactive, "toggled", f=checkbuttonhandler)
#	panel.params$AttachDefaults(label.decimals, 1, 3, 10, 11)
#	label.decimals$SetAlignment(0,0.5)
#	entry.decimals$SetSizeRequest(50,-1)
#	panel.params$AttachDefaults(entry.decimals, 4, 6, 10, 11)
##	panel.params$AttachDefaults(panel.decimals)
#	handler.decimals <- gSignalConnect(entry.decimals, "changed", f=x12input_handler)
	
	#####transform
  #function
# panel.transform$PackStart(checkb.transformactive) 
# gSignalConnect(checkb.transformactive, "toggled", f=checkbuttonhandler)
	panel.transform$Attach(label.transform, 0 , 1, 0, 1)
	#label.transform$SetAlignment(0,0.5)
	combobox.transform$SetSizeRequest(20,-1)
	panel.transform$Attach(combobox.transform, 1, 2, 0, 1, xpadding=5)
	combobox.transform$AppendText("auto")
	combobox.transform$AppendText("log")
	combobox.transform$AppendText("none")
#	panel.params$AttachDefaults(panel.transform)
	handler.transform <- gSignalConnect(combobox.transform, "changed", f=comboboxx12handler)
  
  #power
  panel.transform$Attach(checkb.transform.power, 0, 1, 1, 2)
  panel.transform$Attach(entry.transform.power, 1, 2, 1, 2, xpadding=5)
  handler.transform.power <- gSignalConnect(entry.transform.power, "changed", f=x12input_handler)
  gSignalConnect(checkb.transform.power, "toggled", f=checkbuttonhandler)
  
	#adjust
  combobox.transform.adjust$AppendText("lom")
	combobox.transform.adjust$AppendText("loq")
	combobox.transform.adjust$AppendText("lpyear")
	panel.transform$Attach(checkb.transform.adjust, 0, 1, 2, 3)
	panel.transform$Attach(combobox.transform.adjust, 1, 2, 2, 3, xpadding=5)
	handler.transform.adjust <- gSignalConnect(combobox.transform.adjust, "changed", f=comboboxx12handler)
	gSignalConnect(checkb.transform.adjust, "toggled", f=checkbuttonhandler)
  
	frame.transform$Add(panel.transform)
  panel.params$PackStart(frame.transform)

  
	########Regression-Frame
	panel.regression$Attach(checkb.regressionactive, 0, 1, 0, 1)
	panel.regression$Attach(radiob.regression, 1, 2, 1, 2)
	panel.regression$Attach(radiob.x11regression, 2, 3, 1, 2)
	handler.x11regression <- gSignalConnect(radiob.x11regression, "toggled", f=checkbuttonx12handler)
  handler.regressionactive <- gSignalConnect(checkb.regressionactive, "toggled", f=checkbuttonx12handler)
	
	####regvariables
	panel.regression$Attach(checkb.regvariablesactive, 1, 2, 2, 3)
	gSignalConnect(checkb.regvariablesactive, "toggled", f=checkbuttonhandler)
	#panel.params$PackStart(label.regvariables)
	#label.regvariables$SetAlignment(0, 0.5)
	entry.regvariables$SetSizeRequest(75,-1)
	panel.regression$Attach(entry.regvariables, 2, 3, 2, 3, xpadding=5)
	#	panel.params$AttachDefaults(panel.regvariables, expand=FALSE, padding=5)
	handler.regvariables <- gSignalConnect(entry.regvariables, "changed", f=x12input_handler)
	
	####reguser
	panel.regression$Attach(checkb.reguseractive, 1, 2, 3, 4)
	gSignalConnect(checkb.reguseractive, "toggled", f=checkbuttonhandler)
	# 	panel.params$PackStart(label.reguser)
	# 	label.reguser$SetAlignment(0, 0.5)
	entry.reguser$SetSizeRequest(53,-1)
	panel.regression$Attach(entry.reguser, 2, 3, 3, 4, xpadding=5)
	#	panel.params$AttachDefaults(panel.reguser, expand=FALSE, padding=5)
	handler.reguser <- gSignalConnect(entry.reguser, "changed", f=x12input_handler)
	
	#####regfile
	panel.regression$Attach(checkb.regfileactive, 1, 2, 4, 5)
	gSignalConnect(checkb.regfileactive, "toggled", f=checkbuttonhandler)
	# 	panel.params$PackStart(label.regfile)
	# 	label.regfile$SetAlignment(0, 0.5)
	panel.regression$Attach(filebutton.regfile, 2, 3, 4, 5, xpadding=5)
	gSignalConnect(filebutton.regfile, "file-set", f=filebuttonhandler)
	
	####usertype
	panel.regression$Attach(checkb.usertypeactive, 1, 2, 5, 6)
	gSignalConnect(checkb.usertypeactive, "toggled", f=checkbuttonhandler)
	# 	panel.params$PackStart(label.usertype)
	# 	label.usertype$SetAlignment(0, 0.5)
	entry.usertype$SetSizeRequest(53,-1)
	panel.regression$Attach(entry.usertype, 2, 3, 5, 6, xpadding=5)
	#	panel.params$AttachDefaults(panel.usertype, expand=FALSE, padding=5)
	handler.usertype <- gSignalConnect(entry.usertype, "changed", f=x12input_handler)
	
	#####centeruser
	panel.regression$Attach(checkb.centeruseractive, 1, 2, 6, 7)
	gSignalConnect(checkb.centeruseractive, "toggled", f=checkbuttonhandler)
	# 	panel.params$PackStart(label.centeruser)
	# 	label.centeruser$SetAlignment(0,0.5)
	combobox.centeruser$SetSizeRequest(57,-1)
	panel.regression$Attach(combobox.centeruser, 2, 3, 6, 7, xpadding=5)
	combobox.centeruser$AppendText("mean")
	combobox.centeruser$AppendText("seasonal")
	#	panel.params$AttachDefaults(panel.centeruser, expand=FALSE, padding=5)
	handler.centeruser <- gSignalConnect(combobox.centeruser, "changed", f=comboboxx12handler)
	
	#####regfilestart
	panel.regression$Attach(checkb.regfilestartactive, 1, 2, 7, 8)
	gSignalConnect(checkb.regfilestartactive, "toggled", f=checkbuttonhandler)
	panel.regfilestartstart$PackStart(label.regfilestartstartyear)
	entry.regfilestartstartyear$SetSizeRequest(30,-1)
	panel.regfilestartstart$PackStart(entry.regfilestartstartyear)
	panel.regfilestartstart$PackStart(label.regfilestartstartperiod)
	entry.regfilestartstartperiod$SetSizeRequest(30,-1)
	panel.regfilestartstart$PackStart(entry.regfilestartstartperiod)
	panel.regfilestart$PackStart(panel.regfilestartstart, padding=5)
	#frame.regfilestart$Add(panel.regfilestart)
	panel.regression$Attach(panel.regfilestart, 1, 3, 8, 9, xpadding=5)
	handler.regfilestartstartyear <- gSignalConnect(entry.regfilestartstartyear, "changed", f=x12input_handler)
	handler.regfilestartstartperiod <- gSignalConnect(entry.regfilestartstartperiod, "changed", f=x12input_handler)
	
	####aictest
	panel.regression$Attach(checkb.aictestactive, 1, 2, 9, 10)
	gSignalConnect(checkb.aictestactive, "toggled", f=checkbuttonhandler)
	#panel.params$PackStart(label.aictest)
# 	label.aictest$SetAlignment(0, 0.5)
	entry.aictest$SetSizeRequest(53,-1)
	panel.regression$Attach(entry.aictest, 2, 3, 9, 10, xpadding=5)
	#	panel.params$AttachDefaults(panel.aictest, expand=FALSE, padding=5)
	handler.aictest <- gSignalConnect(entry.aictest, "changed", f=x12input_handler)
	
	frame.regression$Add(panel.regression)
	panel.params$PackStart(frame.regression)
	
  ########Outlier-Frame
	#outlier
	panel.outlier$AttachDefaults(checkb.outlieractive, 0, 1, 1, 2)
  panel.outlier$AttachDefaults(label.outlier, 1, 2, 2, 3)
	panel.outlier$AttachDefaults(checkb.outlierall, 2, 3, 3, 4)
	panel.outlier$AttachDefaults(checkb.outlierAO, 3, 4, 3, 4)
	panel.outlier$AttachDefaults(checkb.outlierTC, 2, 3, 4, 5)
	panel.outlier$AttachDefaults(checkb.outlierLS, 3, 4, 4, 5)
	gSignalConnect(checkb.outlieractive, "toggled", f=checkbuttonhandler)
	handler.outlierall <- gSignalConnect(checkb.outlierall, "toggled", checkbuttonx12handler)
	handler.outlierTC <- gSignalConnect(checkb.outlierTC, "toggled", checkbuttonx12handler)
	handler.outlierAO <- gSignalConnect(checkb.outlierAO, "toggled", checkbuttonx12handler)
	handler.outlierLS <- gSignalConnect(checkb.outlierLS, "toggled", checkbuttonx12handler)
	#critical
	entry.criticalall$SetSizeRequest(30, -1)
	entry.criticalTC$SetSizeRequest(10, -1)
	entry.criticalAO$SetSizeRequest(10, -1)
	entry.criticalLS$SetSizeRequest(10, -1)
	gSignalConnect(radiob.criticalspecific, "toggled", f=checkbuttonhandler)
	gSignalConnect(radiob.criticalall, "toggled", f=checkbuttonhandler)
	gSignalConnect(checkb.criticalactive, "toggled", f=checkbuttonhandler)
	panel.outlier$Attach(checkb.criticalactive, 1, 2, 5, 6)
	panel.outlier$AttachDefaults(radiob.criticalall, 2, 3, 5, 6)
	panel.outlier$AttachDefaults(entry.criticalall, 3, 5, 5, 6)
	panel.outlier$AttachDefaults(radiob.criticalspecific, 2, 3, 6, 7)
	panel.outlier$AttachDefaults(label.criticalAO, 2, 3, 7, 8)
	panel.outlier$AttachDefaults(label.criticalLS, 2, 3, 8, 9)
	panel.outlier$AttachDefaults(label.criticalTC, 2, 3, 9, 10)
	panel.outlier$AttachDefaults(entry.criticalAO, 3, 5, 7, 8)
	panel.outlier$AttachDefaults(entry.criticalLS, 3, 5, 8, 9)
	panel.outlier$AttachDefaults(entry.criticalTC, 3, 5, 9, 10)
	handler.criticalall <- gSignalConnect(entry.criticalall, "changed", f=x12input_handler)
	handler.criticalAO <- gSignalConnect(entry.criticalAO, "changed", f=x12input_handler)
	handler.criticalLS <- gSignalConnect(entry.criticalLS, "changed", f=x12input_handler)
	handler.criticalTC <- gSignalConnect(entry.criticalTC, "changed", f=x12input_handler)
	#outlierspan
	panel.outlier$Attach(checkb.outlierspanactive, 1, 2, 10, 11)
  panel.outlier$Attach(label.outlierspan.year, 2, 3, 11, 12)
	panel.outlier$Attach(label.outlierspan.period, 3, 4, 11, 12)
	gSignalConnect(checkb.outlierspanactive, "toggled", f=checkbuttonhandler)
  panel.outlier$Attach(checkb.outlierspan.start, 2, 3, 10, 11)
	entry.outlierspan.start1$SetSizeRequest(10,-1)
	panel.outlier$Attach(entry.outlierspan.start1, 2, 3, 12, 13)
	entry.outlierspan.start2$SetSizeRequest(10,-1)
	panel.outlier$Attach(entry.outlierspan.start2, 3, 4, 12, 13)
	panel.outlier$Attach(checkb.outlierspan.end, 2, 3, 13, 14)
	entry.outlierspan.end1$SetSizeRequest(10,-1)
	panel.outlier$Attach(entry.outlierspan.end1, 2, 3, 14, 15)
	entry.outlierspan.end2$SetSizeRequest(10,-1)
	panel.outlier$Attach(entry.outlierspan.end2, 3, 4, 14, 15)
  gSignalConnect(checkb.outlierspan.start, "toggled", f=checkbuttonhandler)
  gSignalConnect(checkb.outlierspan.end, "toggled", f=checkbuttonhandler)
	handler.outlierspan.start1 <- gSignalConnect(entry.outlierspan.start1, "changed", f=x12input_handler)
	handler.outlierspan.start2<- gSignalConnect(entry.outlierspan.start2, "changed", f=x12input_handler)
	handler.outlierspan.end1 <- gSignalConnect(entry.outlierspan.end1, "changed", f=x12input_handler)
	handler.outlierspan.end2<- gSignalConnect(entry.outlierspan.end2, "changed", f=x12input_handler)
	#outliermethod
	panel.outlier$Attach(checkb.outliermethodactive, 1, 2, 15, 16) 
	gSignalConnect(checkb.outliermethodactive, "toggled", f=checkbuttonhandler)
	combobox.outliermethod$SetSizeRequest(57,-1)
	panel.outlier$Attach(combobox.outliermethod, 2, 4, 15, 16)
	combobox.outliermethod$AppendText("addone")
	combobox.outliermethod$AppendText("addall")
	handler.outliermethod <- gSignalConnect(combobox.outliermethod, "changed", f=comboboxx12handler)
  
  frame.outlier$Add(panel.outlier)
  panel.params$PackStart(frame.outlier)
  
	######## Arima-Frame
	frame.arima$Add(panel.arima)
	panel.params$PackStart(frame.arima)
	#####arima	
	panel.arima$AttachDefaults(checkb.arimaactive, 1, 2, 0, 1)
	gSignalConnect(checkb.arimaactive, "toggled", f=checkbuttonhandler)
	entry.arima1$SetSizeRequest(8,-1)
	panel.arima$AttachDefaults(entry.arima1, 2, 3, 0, 1)
	entry.arima2$SetSizeRequest(8,-1)
	panel.arima$AttachDefaults(entry.arima2, 3, 4, 0, 1)
	entry.arima3$SetSizeRequest(8,-1)
	panel.arima$AttachDefaults(entry.arima3, 4, 5, 0, 1)
#	panel.params$AttachDefaults(panel.arima, expand=FALSE, padding=5)
	handler.arima1 <- gSignalConnect(entry.arima1, "changed", f=x12input_handler)
	handler.arima2 <- gSignalConnect(entry.arima2, "changed", f=x12input_handler)
	handler.arima3 <- gSignalConnect(entry.arima3, "changed", f=x12input_handler)	
	#####sarima
	panel.arima$AttachDefaults(checkb.sarimaactive, 1, 2, 1, 2)
	gSignalConnect(checkb.sarimaactive, "toggled", f=checkbuttonhandler)
	entry.sarima1$SetSizeRequest(8,-1)
	panel.arima$AttachDefaults(entry.sarima1, 2, 3, 1, 2)
	entry.sarima2$SetSizeRequest(8,-1)
	panel.arima$AttachDefaults(entry.sarima2, 3, 4, 1, 2)
	entry.sarima3$SetSizeRequest(8,-1)
	panel.arima$AttachDefaults(entry.sarima3, 4, 5, 1, 2)
#	panel.params$AttachDefaults(panel.sarima, expand=FALSE, padding=5)
	handler.sarima1 <- gSignalConnect(entry.sarima1, "changed", f=x12input_handler)
	handler.sarima2 <- gSignalConnect(entry.sarima2, "changed", f=x12input_handler)
	handler.sarima3 <- gSignalConnect(entry.sarima3, "changed", f=x12input_handler)
  #####as
  panel.arima$Attach(checkb.arima.ar, 1, 2, 2, 3)
	panel.arima$Attach(entry.arima.ar, 2, 5, 2, 3)
  handler.arima.ar <- gSignalConnect(entry.arima.ar, "changed", f=x12input_handler)
	gSignalConnect(checkb.arima.ar, "toggled", f=checkbuttonhandler)
	#####as
	panel.arima$Attach(checkb.arima.ma, 1, 2, 3, 4)
	panel.arima$Attach(entry.arima.ma, 2, 5, 3, 4)
	handler.arima.ma <- gSignalConnect(entry.arima.ma, "changed", f=x12input_handler)
	gSignalConnect(checkb.arima.ma, "toggled", f=checkbuttonhandler)
  
  #######Automdl-Frame
  panel.automdl$Attach(checkb.automdl, 0, 1, 0, 1)
	panel.automdl$Attach(checkb.acceptdefault, 1, 2, 1, 2)
	panel.automdl$Attach(checkb.balanced, 1, 2, 2, 3)
	#####maxorder
	#panel.automdl$Attach(checkb.maxorderactive, 1, 2, 3, 4)
  gtkMiscSetAlignment(label.maxorder, 0, 0.5)
  panel.automdl$Attach(label.maxorder, 1, 2, 3, 4)
  #gSignalConnect(checkb.maxorderactive, "toggled", f=checkbuttonhandler)
	entry.maxorder1$SetSizeRequest(10,-1)
	panel.automdl$Attach(entry.maxorder1, 2, 3, 3, 4)
  entry.maxorder2$SetSizeRequest(10,-1)
  panel.automdl$Attach(entry.maxorder2, 3, 4, 3, 4)
	handler.maxorder1 <- gSignalConnect(entry.maxorder1, "changed", f=x12input_handler)
  handler.maxorder2 <- gSignalConnect(entry.maxorder2, "changed", f=x12input_handler)
	#####maxdiff
	#panel.automdl$Attach(checkb.maxdiffactive, 1, 2, 4, 5)
  gtkMiscSetAlignment(label.maxdiff, 0, 0.5)
  panel.automdl$Attach(label.maxdiff, 1, 2, 4, 5)
	#gSignalConnect(checkb.maxdiffactive, "toggled", f=checkbuttonhandler)
	entry.maxdiff1$SetSizeRequest(10,-1)
	panel.automdl$Attach(entry.maxdiff1, 2, 3, 4, 5)
  entry.maxdiff2$SetSizeRequest(10,-1)
  panel.automdl$Attach(entry.maxdiff2, 3, 4, 4, 5)
	handler.maxdiff1 <- gSignalConnect(entry.maxdiff1, "changed", f=x12input_handler)
  handler.maxdiff2 <- gSignalConnect(entry.maxdiff2, "changed", f=x12input_handler)
	handlercheckb.automdl <- gSignalConnect(checkb.automdl, "toggled", checkbuttonx12handler)
	handlercheckb.acceptdefault <- gSignalConnect(checkb.acceptdefault, "toggled", checkbuttonx12handler)
	handlercheckb.balanced <- gSignalConnect(checkb.balanced, "toggled", checkbuttonx12handler)
  
  frame.automdl$Add(panel.automdl)
  panel.params$PackStart(frame.automdl)
	

	
# ####tblnames
# panel.tblnames$PackStart(checkb.tblnamesactive)
# gSignalConnect(checkb.tblnamesactive, "toggled", f=checkbuttonhandler)
# panel.tblnames$PackStart(label.tblnames)
# label.tblnames$SetAlignment(0, 0.5)
# entry.tblnames$SetSizeRequest(53,-1)
# panel.tblnames$PackStart(entry.tblnames)
# panel.params$AttachDefaults(panel.tblnames, expand=FALSE, padding=5)
# 
# ####Rtblnames
# panel.Rtblnames$PackStart(checkb.Rtblnamesactive)
# gSignalConnect(checkb.Rtblnamesactive, "toggled", f=checkbuttonhandler)
# panel.Rtblnames$PackStart(label.Rtblnames)
# label.Rtblnames$SetAlignment(0, 0.5)
# entry.Rtblnames$SetSizeRequest(53,-1)
# panel.Rtblnames$PackStart(entry.Rtblnames)
# panel.params$AttachDefaults(panel.Rtblnames, expand=FALSE, padding=5)
# 
# #####x12path
# panel.x12path$PackStart(checkb.x12pathactive)
# gSignalConnect(checkb.x12pathactive, "toggled", f=checkbuttonhandler)
# panel.x12path$PackStart(label.x12path)
# label.x12path$SetAlignment(0, 0.5)
# panel.x12path$PackStart(filebutton.x12path)
# panel.params$AttachDefaults(panel.x12path, expand=FALSE, padding=5)
# 
# #####x13path
# panel.x13path$PackStart(checkb.x13pathactive)
# gSignalConnect(checkb.x13pathactive, "toggled", f=checkbuttonhandler)
# panel.x13path$PackStart(label.x13path)
# label.x13path$SetAlignment(0, 0.5)
# panel.x13path$PackStart(filebutton.x13path)
# panel.params$AttachDefaults(panel.x13path, expand=FALSE, padding=5)
# 
# #####use
# panel.use$PackStart(checkb.useactive)
# gSignalConnect(checkb.useactive, "toggled", f=checkbuttonhandler)
# panel.use$PackStart(label.use)
# label.use$SetAlignment(0, 0.5)
# combobox.use$SetSizeRequest(57,-1)
# panel.use$PackStart(combobox.use, expand=TRUE)
# combobox.use$AppendText("x12")
# combobox.use$AppendText("x13")
# panel.params$AttachDefaults(panel.use, expand=FALSE, padding=5)

  
  #####Identify-Frame
  panel.identify$Attach(checkb.identify, 0, 1, 0, 1)
  panel.identify$Attach(checkb.identify.diff, 1, 2, 1, 2)
	panel.identify$Attach(entry.identify.diff, 2, 3, 1, 2, xpadding=5)
  entry.identify.diff$SetSizeRequest(35,-1)
	panel.identify$Attach(checkb.identify.sdiff, 1, 2, 2, 3)
	panel.identify$Attach(entry.identify.sdiff, 2, 3, 2, 3, xpadding=5)
	entry.identify.sdiff$SetSizeRequest(35,-1)
	panel.identify$Attach(checkb.identify.maxlag, 1, 2, 3, 4)
	panel.identify$Attach(entry.identify.maxlag, 2, 3, 3, 4, xpadding=5)
	entry.identify.maxlag$SetSizeRequest(35,-1)
	handler.identify.sdiff <- gSignalConnect(entry.identify.sdiff, "changed", f=x12input_handler)
	handler.identify.diff <- gSignalConnect(entry.identify.diff, "changed", f=x12input_handler)
	handler.identify.maxlag <- gSignalConnect(entry.identify.maxlag, "changed", f=x12input_handler)
  gSignalConnect(checkb.identify.diff, "toggled", checkbuttonhandler)
	gSignalConnect(checkb.identify.sdiff, "toggled", checkbuttonhandler)
	gSignalConnect(checkb.identify.maxlag, "toggled", checkbuttonhandler)
  gSignalConnect(checkb.identify, "toggled", checkbuttonx12handler)
  
  frame.identify$Add(panel.identify)
  panel.params$PackStart(frame.identify)
	
# ####file
# panel.file$PackStart(checkb.fileactive)
# gSignalConnect(checkb.fileactive, "toggled", f=checkbuttonhandler)
# panel.file$PackStart(label.file)
# label.file$SetAlignment(0, 0.5)
# entry.file$SetSizeRequest(53,-1)
# panel.file$PackStart(entry.file)
# panel.params$AttachDefaults(panel.file, expand=FALSE, padding=5)
	
  ######Forecast-Frame
	####forecast_years
	panel.forecast$Attach(checkb.forecast_yearsactive, 0, 1, 0, 1)
	gSignalConnect(checkb.forecast_yearsactive, "toggled", f=checkbuttonhandler)
	entry.forecast_years$SetSizeRequest(53,-1)
	panel.forecast$Attach(entry.forecast_years, 1, 2, 0, 1)
	handler.forecast_years <- gSignalConnect(entry.forecast_years, "changed", f=x12input_handler)
	####backcast_years
	panel.forecast$Attach(checkb.backcast_yearsactive, 0 ,1, 1, 2)
	gSignalConnect(checkb.backcast_yearsactive, "toggled", f=checkbuttonhandler)
	entry.backcast_years$SetSizeRequest(53,-1)
	panel.forecast$Attach(entry.backcast_years, 1, 2, 1, 2)
	handler.backcast_years <- gSignalConnect(entry.backcast_years, "changed", f=x12input_handler)
	####forecast_conf
	panel.forecast$Attach(checkb.forecast_confactive, 0, 1, 2, 3)
  gSignalConnect(checkb.forecast_confactive, "toggled", f=checkbuttonhandler)
	entry.forecast_conf$SetSizeRequest(53,-1)
	panel.forecast$Attach(entry.forecast_conf, 1, 2, 2, 3)
	handler.forecast_conf <- gSignalConnect(entry.forecast_conf, "changed", f=x12input_handler)
  ##
  frame.forecast$Add(panel.forecast)
  panel.params$PackStart(frame.forecast)
  
  ######Estimate-Frame
	panel.estimate$Attach(checkb.estimate, 0, 1, 0, 1)
  panel.estimate$Attach(checkb.estOutofsample, 1, 2, 1, 2)
  frame.estimate$Add(panel.estimate)
  panel.params$PackStart(frame.estimate)
	handlercheckb.estimate <- gSignalConnect(checkb.estimate, "toggled", checkbuttonx12handler)
	handlercheckb.estOutofsample <- gSignalConnect(checkb.estOutofsample, "toggled", checkbuttonx12handler)
  
	######Check-Frame
	panel.check$Attach(checkb.check, 0, 1, 0, 1)
	panel.check$Attach(checkb.check.maxlag, 1, 2, 1, 2)
  entry.check.maxlag$SetSizeRequest(35,-1)
	panel.check$Attach(entry.check.maxlag, 2, 3, 1, 2)
	frame.check$Add(panel.check)
	panel.params$PackStart(frame.check)
	handlercheckb.check <- gSignalConnect(checkb.check, "toggled", checkbuttonx12handler)
	handlercheckb.check.maxlag <- gSignalConnect(checkb.check.maxlag, "toggled", checkbuttonhandler)
  handler.check.maxlag <- gSignalConnect(entry.check.maxlag, "changed", x12input_handler)
  
  ######Slidingspans-Frame
  panel.slidingspans$Attach(checkb.slidingspans, 0, 1, 0, 1)
  frame.slidingspans$Add(panel.slidingspans)
  panel.params$PackStart(frame.slidingspans)
	handlercheckb.slidingspans <- gSignalConnect(checkb.slidingspans, "toggled", checkbuttonx12handler)
  panel.slidingspans$Attach(checkb.slidingspans.fixmdl, 1, 2, 1, 2)
	handlercheckb.slidingspans.fixmdl <- gSignalConnect(checkb.slidingspans.fixmdl, "toggled", checkbuttonhandler)
	panel.slidingspans$Attach(combobox.slidingspans.fixmdl, 2, 6, 1, 2)
  combobox.slidingspans.fixmdl$AppendText("yes")
	combobox.slidingspans.fixmdl$AppendText("no")
	combobox.slidingspans.fixmdl$AppendText("clear")
  handler.slidingspans.fixmdl <- gSignalConnect(combobox.slidingspans.fixmdl, "changed", f=comboboxx12handler)
	#fixreg
	panel.slidingspans$Attach(checkb.slidingspans.fixreg, 1, 2, 2, 3)
	panel.slidingspans$Attach(checkb.slidingspans.fixreg1, 2, 3, 2, 3)
	panel.slidingspans$Attach(checkb.slidingspans.fixreg2, 3, 4, 2, 3)
	panel.slidingspans$Attach(checkb.slidingspans.fixreg3, 2, 3, 3, 4)
	panel.slidingspans$Attach(checkb.slidingspans.fixreg4, 3, 4, 3, 4)
	handlercheckb.slidingspans.fixreg <- gSignalConnect(checkb.slidingspans.fixreg, "toggled", checkbuttonhandler)
  handler.slidingspans.fixreg1 <- gSignalConnect(checkb.slidingspans.fixreg1, "toggled", checkbuttonx12handler)
	handler.slidingspans.fixreg2 <- gSignalConnect(checkb.slidingspans.fixreg2, "toggled", checkbuttonx12handler)
	handler.slidingspans.fixreg3 <- gSignalConnect(checkb.slidingspans.fixreg3, "toggled", checkbuttonx12handler)
	handler.slidingspans.fixreg4 <- gSignalConnect(checkb.slidingspans.fixreg4, "toggled", checkbuttonx12handler)
	#length
  panel.slidingspans$Attach(checkb.slidingspans.length, 1, 2, 4, 5)
	handlercheckb.slidingspans.length <- gSignalConnect(checkb.slidingspans.length, "toggled", checkbuttonhandler)
	panel.slidingspans$Attach(entry.slidingspans.length, 2, 6, 4, 5)
	handler.slidingspans.length <- gSignalConnect(entry.slidingspans.length, "changed", x12input_handler)
	#numspan
  panel.slidingspans$Attach(checkb.slidingspans.numspans, 1, 2, 5, 6)
	handlercheckb.slidingspans.numspans <- gSignalConnect(checkb.slidingspans.numspans, "toggled", checkbuttonhandler)
	panel.slidingspans$Attach(entry.slidingspans.numspans, 2, 6, 5, 6)
	handler.slidingspans.numspans <- gSignalConnect(entry.slidingspans.numspans, "changed", x12input_handler)
	#outlier
  panel.slidingspans$Attach(checkb.slidingspans.outlier, 1, 2, 6, 7)
	handlercheckb.slidingspans.outlier <- gSignalConnect(checkb.slidingspans.outlier, "toggled", checkbuttonhandler)
	panel.slidingspans$Attach(combobox.slidingspans.outlier, 2, 6, 6, 7)
	combobox.slidingspans.outlier$AppendText("keep")
	combobox.slidingspans.outlier$AppendText("remove")
	combobox.slidingspans.outlier$AppendText("yes")
	handler.slidingspans.outlier <- gSignalConnect(combobox.slidingspans.outlier, "changed", comboboxx12handler)
	#additivesa
	panel.slidingspans$Attach(checkb.slidingspans.additivesa, 1, 2, 7, 8)
	handlercheckb.slidingspans.additivesa <- gSignalConnect(checkb.slidingspans.additivesa, "toggled", checkbuttonhandler)
	panel.slidingspans$Attach(combobox.slidingspans.additivesa, 2, 6, 7, 8)
	combobox.slidingspans.additivesa$AppendText("difference")
	combobox.slidingspans.additivesa$AppendText("percent")
	handler.slidingspans.outlier <- gSignalConnect(combobox.slidingspans.additivesa, "changed", comboboxx12handler)
	#start
	panel.slidingspans$Attach(checkb.slidingspans.start, 1, 2, 9, 10)
	handlercheckb.slidingspans.start <- gSignalConnect(checkb.slidingspans.start, "toggled", checkbuttonhandler)
	entry.slidingspans.start1$SetSizeRequest(10,-1)
	entry.slidingspans.start2$SetSizeRequest(10,-1)
	panel.slidingspans$Attach(label.slidingspans.start1, 2, 3, 9, 10)
	panel.slidingspans$Attach(label.slidingspans.start2, 3, 6, 9, 10)
	panel.slidingspans$Attach(entry.slidingspans.start1, 2, 3, 10, 11)
	panel.slidingspans$Attach(entry.slidingspans.start2, 3, 6, 10, 11)
  handler.slidingspans.start1 <- gSignalConnect(entry.slidingspans.start1, "changed", x12input_handler)
	handler.slidingspans.start2 <- gSignalConnect(entry.slidingspans.start2, "changed", x12input_handler)
	handler.slidingspans.outlier <- gSignalConnect(combobox.slidingspans.additivesa, "changed", comboboxx12handler)
  
  #########History-Frame
  panel.historyparam$Attach(checkb.historyactive, 0, 1, 0, 1)
  panel.historyparam$Attach(checkb.history.estimates, 1, 2, 1, 2)
  panel.historyparam$Attach(checkb.history.estimatessadj, 2, 3, 1, 2)
	panel.historyparam$Attach(checkb.history.estimatessadjchng, 3, 4, 1, 2)
	panel.historyparam$Attach(checkb.history.estimatestrend, 2, 3, 2, 3)
	panel.historyparam$Attach(checkb.history.estimatestrendchng, 3, 4, 2, 3)
	panel.historyparam$Attach(checkb.history.estimatesseasonal, 2, 3, 3, 4)
	panel.historyparam$Attach(checkb.history.estimatesfcst, 3, 4, 3, 4)
	panel.historyparam$Attach(checkb.history.estimatesaic, 2, 3, 4, 5)
  handler.history <- gSignalConnect(checkb.historyactive, "toggled", checkbuttonx12handler)
	handler.history.estimates <- gSignalConnect(checkb.history.estimates, "toggled", checkbuttonhandler)
	handler.history.estimatessadj <- gSignalConnect(checkb.history.estimatessadj, "toggled", checkbuttonx12handler)
	handler.history.estimatessadjchng <- gSignalConnect(checkb.history.estimatessadjchng, "toggled", checkbuttonx12handler)
	handler.history.estimatestrend <- gSignalConnect(checkb.history.estimatestrend, "toggled", checkbuttonx12handler)
	handler.history.estimatestrendchng <- gSignalConnect(checkb.history.estimatestrendchng, "toggled", checkbuttonx12handler)
	handler.history.estimatesseasonal <- gSignalConnect(checkb.history.estimatesseasonal, "toggled", checkbuttonx12handler)
	handler.history.estimatesfcst <- gSignalConnect(checkb.history.estimatesfcst, "toggled", checkbuttonx12handler)
	handler.history.estimatesaic <- gSignalConnect(checkb.history.estimatesaic, "toggled", checkbuttonx12handler)
  #fixmdl
  panel.historyparam$Attach(checkb.history.fixmdl, 1, 2, 5, 6)
  #fixreg
	panel.historyparam$Attach(checkb.history.fixreg, 1, 2, 6, 7)
	panel.historyparam$Attach(checkb.history.fixreg1, 2, 3, 6, 7)
	panel.historyparam$Attach(checkb.history.fixreg2, 3, 4, 6, 7)
	panel.historyparam$Attach(checkb.history.fixreg3, 2, 3, 7, 8)
	panel.historyparam$Attach(checkb.history.fixreg4, 3, 4, 7, 8)
  handler.history.fixreg <- gSignalConnect(checkb.history.fixreg, "toggled", checkbuttonhandler)
	handler.history.fixreg1 <- gSignalConnect(checkb.history.fixreg1, "toggled", checkbuttonx12handler)
	handler.history.fixreg2 <- gSignalConnect(checkb.history.fixreg2, "toggled", checkbuttonx12handler)
	handler.history.fixreg3 <- gSignalConnect(checkb.history.fixreg3, "toggled", checkbuttonx12handler)
	handler.history.fixreg4 <- gSignalConnect(checkb.history.fixreg4, "toggled", checkbuttonx12handler)
	#outlier
	combobox.history.outlier$AppendText("keep")
	combobox.history.outlier$AppendText("remove")
	combobox.history.outlier$AppendText("auto")
	panel.historyparam$Attach(checkb.history.outlier, 1, 2, 8, 9)
	panel.historyparam$Attach(combobox.history.outlier, 2, 4, 8, 9)
  gSignalConnect(checkb.history.outlier, "toggled", checkbuttonhandler)
	handler.history.outlier <- gSignalConnect(combobox.history.outlier, "changed", comboboxx12handler)
	#target
	combobox.history.target$AppendText("final")
	combobox.history.target$AppendText("concurrent")
#	combobox.history.target$AppendText("auto")
	panel.historyparam$Attach(checkb.history.target, 1, 2, 9, 10)
	panel.historyparam$Attach(combobox.history.target, 2, 4, 9, 10)
	gSignalConnect(checkb.history.target, "toggled", checkbuttonhandler)
	handler.history.target <- gSignalConnect(combobox.history.target, "changed", comboboxx12handler)	
	#sadjlags
	panel.historyparam$Attach(checkb.history.sadjlags, 1, 2, 10, 11)
	panel.historyparam$Attach(entry.history.sadjlags, 2, 4, 10, 11)
	gSignalConnect(checkb.history.sadjlags, "toggled", checkbuttonhandler)
	handler.history.sadjlags <- gSignalConnect(entry.history.sadjlags, "changed", f=x12input_handler)
	#trendlags
	panel.historyparam$Attach(checkb.history.trendlags, 1 ,2 ,11, 12)
	panel.historyparam$Attach(entry.history.trendlags, 2 ,4, 11, 12)
	gSignalConnect(checkb.history.trendlags, "toggled", checkbuttonhandler)
	handler.history.trendlags <- gSignalConnect(entry.history.trendlags, "changed", f=x12input_handler)
  frame.historyparam$Add(panel.historyparam)
  panel.params$PackStart(frame.historyparam)
	#start
	panel.historyparam$Attach(checkb.history.start, 1, 2, 12, 13)
	panel.historyparam$Attach(label.history.startyear, 2, 3, 12, 13)
	panel.historyparam$Attach(label.history.startperiod, 3, 4, 12, 13)
  entry.history.startperiod$SetSizeRequest(30,-1)
	entry.history.startyear$SetSizeRequest(30,-1)
  panel.historyparam$Attach(entry.history.startyear, 2, 3, 13, 14)
	panel.historyparam$Attach(entry.history.startperiod, 3, 4, 13, 14)
  gSignalConnect(checkb.history.start, "toggled", checkbuttonhandler)
	handler.history.start1 <- gSignalConnect(entry.history.startyear, "changed", x12input_handler)
	handler.history.start2 <- gSignalConnect(entry.history.startperiod, "changed", x12input_handler)
	handler.history.startyear <- 1
	handler.history.startperiod <- 1
  
  #########X11-Frame
	#handlercheckb.x11regress <- gSignalConnect(checkb.x11regress, "toggled", checkbuttonx12handler)
	#sigmalim
	panel.x11$Attach(checkb.sigmalimactive, 0, 1, 0, 1)
	gSignalConnect(checkb.sigmalimactive, "toggled", f=checkbuttonhandler)
	entry.sigmalim1$SetSizeRequest(10,-1)
	panel.x11$Attach(entry.sigmalim1, 1, 2, 0, 1)
	entry.sigmalim2$SetSizeRequest(10,-1)
	panel.x11$Attach(entry.sigmalim2, 2, 3, 0, 1)
	handler.sigmalim1 <- gSignalConnect(entry.sigmalim1, "changed", f=x12input_handler)
	handler.sigmalim2<- gSignalConnect(entry.sigmalim2, "changed", f=x12input_handler)
  #type
	panel.x11$Attach(checkb.x11.type, 0, 1, 1, 2)
  panel.x11$Attach(combobox.x11.type, 1, 3, 1, 2)
  combobox.x11.type$AppendText("summary")
	combobox.x11.type$AppendText("trend")
	combobox.x11.type$AppendText("sa")
  handler.x11.type <- gSignalConnect(combobox.x11.type, "changed", f=comboboxx12handler)
  gSignalConnect(checkb.x11.type, "toggled", checkbuttonhandler)
  #sfshort
	panel.x11$Attach(checkb.sfshort, 0, 1, 2, 3)
	handlercheckb.sfshort <- gSignalConnect(checkb.sfshort, "toggled", checkbuttonx12handler)
  #samode
	panel.x11$Attach(checkb.samodeactive, 0, 1, 3, 4) 
	gSignalConnect(checkb.samodeactive, "toggled", f=checkbuttonhandler)
	combobox.samode$SetSizeRequest(57,-1)
	panel.x11$Attach(combobox.samode, 1, 3, 3, 4)
	combobox.samode$AppendText("mult")
	combobox.samode$AppendText("add")
	combobox.samode$AppendText("pseudoadd")
	combobox.samode$AppendText("logadd")
	handler.samode <- gSignalConnect(combobox.samode, "changed", f=comboboxx12handler)
	#seasonalma
	panel.x11$Attach(checkb.seasonalmaactive, 0, 1, 4, 5)
	gSignalConnect(checkb.seasonalmaactive, "toggled", f=checkbuttonhandler)
	entry.seasonalma$SetSizeRequest(53,-1)
	panel.x11$Attach(entry.seasonalma, 1, 3, 4, 5)
	handler.seasonalma <- gSignalConnect(entry.seasonalma, "changed", f=x12input_handler)
	#trendma
	panel.x11$Attach(checkb.trendmaactive, 0, 1, 5, 6)
	gSignalConnect(checkb.trendmaactive, "toggled", f=checkbuttonhandler)
	entry.trendma$SetSizeRequest(53,-1)
	panel.x11$Attach(entry.trendma, 1, 3, 5, 6)
	handler.trendma <- gSignalConnect(entry.trendma, "changed", f=x12input_handler)
  #div
	panel.x11$Attach(checkb.x11appendfcst, 0, 1, 6, 7)
	panel.x11$Attach(checkb.x11appendbcst, 0, 1, 7, 8)
	panel.x11$Attach(checkb.x11excludefcst, 0, 1, 8, 9)
	handlercheckb.x11appendfcst <- gSignalConnect(checkb.x11appendfcst, "toggled", checkbuttonx12handler)
	handlercheckb.x11appendbcst <- gSignalConnect(checkb.x11appendbcst, "toggled", checkbuttonx12handler)
	handlercheckb.x11excludefcst <- gSignalConnect(checkb.x11excludefcst, "toggled", checkbuttonx12handler)
	#x11calendarsigma
	panel.x11$Attach(checkb.x11calendarsigmaactive, 0, 1, 9, 10) 
	gSignalConnect(checkb.x11calendarsigmaactive, "toggled", f=checkbuttonhandler)
	combobox.x11calendarsigma$SetSizeRequest(57,-1)
	panel.x11$Attach(combobox.x11calendarsigma, 1, 3, 9, 10)
	combobox.x11calendarsigma$AppendText("all")
	combobox.x11calendarsigma$AppendText("signif")
	combobox.x11calendarsigma$AppendText("select")
	handler.x11calendarsigma <- gSignalConnect(combobox.x11calendarsigma, "changed", f=comboboxx12handler)
  #final
	panel.x11$Attach(checkb.x11.finalactive, 0, 1, 10, 11)
  panel.x11$Attach(checkb.x11.finalAO, 1, 2, 10, 11)
	panel.x11$Attach(checkb.x11.finalLS, 2, 3, 10, 11)
	panel.x11$Attach(checkb.x11.finalTC, 1, 2, 11, 12)
	panel.x11$Attach(checkb.x11.finaluser, 2, 3, 11, 12)
	panel.x11$Attach(checkb.x11.finalnone, 1, 2, 12, 13)
  gSignalConnect(checkb.x11.finalactive, "toggled", checkbuttonhandler)
  handler.x11.finalAO <- gSignalConnect(checkb.x11.finalAO, "toggled", checkbuttonx12handler)
	handler.x11.finalLS <- gSignalConnect(checkb.x11.finalLS, "toggled", checkbuttonx12handler)
	handler.x11.finalTC <- gSignalConnect(checkb.x11.finalTC, "toggled", checkbuttonx12handler)
	handler.x11.finaluser <- gSignalConnect(checkb.x11.finaluser, "toggled", checkbuttonx12handler)
	handler.x11.finalnone <- gSignalConnect(checkb.x11.finalnone, "toggled", checkbuttonx12handler)
  
  frame.x11$Add(panel.x11)
  panel.params$PackStart(frame.x11)
  
#   ######Seats-Frame
# 	panel.seats$Attach(checkb.seats, 0, 1, 0, 1)
# 	handlercheckb.seats <- gSignalConnect(checkb.seats, "toggled", checkbuttonx12handler)
# 	####seatsparameter
# 	panel.seats$Attach(checkb.seatsparameteractive, 1, 2, 1, 2)
# 	gSignalConnect(checkb.seatsparameteractive, "toggled", f=checkbuttonhandler)
# 	entry.seatsparameter$SetSizeRequest(53,-1)
# 	panel.seats$Attach(entry.seatsparameter, 2, 3, 1, 2)
# 	handler.seatsparameter <- gSignalConnect(entry.seatsparameter, "changed", f=x12input_handler)
#   frame.seats$Add(panel.seats)
#   panel.params$PackStart(frame.seats)
	
# 	#####TOOLTIPS
# 	#frame.span$SetTooltipMarkup(string.span)
# 	frame.modelspan$SetTooltipMarkup(string.modelspan)
# #	label.decimals$SetTooltipMarkup(string.decimals)
# 	label.transform$SetTooltipMarkup(string.transform)
# 	label.arima$SetTooltipMarkup(string.arima)
# 	label.sarima$SetTooltipMarkup(string.sarima)
# 	checkb.automdl$SetTooltipMarkup(string.automdl)
# 	checkb.balanced$SetTooltipMarkup(string.balanced)
# 	checkb.seats$SetTooltipMarkup(string.seats)
# 	checkb.estimate$SetTooltipMarkup(string.estimate)
# 	checkb.estOutofsample$SetTooltipMarkup(string.estimateoutofsamples)
# 	checkb.slidingspans$SetTooltipMarkup(string.slidingspans)
# 	checkb.onlytd$SetTooltipMarkup(string.onlytd)
# 	checkb.sfshort$SetTooltipMarkup(string.sfshort)
# 	checkb.x11appendfcst$SetTooltipMarkup(string.x11appendfcst)
# 	checkb.x11appendbcst$SetTooltipMarkup(string.x11appendfbst)
# 	checkb.x11excludefcst$SetTooltipMarkup(string.x11excludefcst)
# 	checkb.x11regress$SetTooltipMarkup(string.x11regress)
# 	label.maxorder$SetTooltipMarkup(string.maxorder)
# 	label.maxdiff$SetTooltipMarkup(string.maxdiff)
# 	label.regvariables$SetTooltipMarkup(string.regvariables)
# 	label.reguser$SetTooltipMarkup(string.reguser)
# 	label.regfile$SetTooltipMarkup(string.regfile)
# 	label.usertype$SetTooltipMarkup(string.usertype)
# 	label.centeruser$SetTooltipMarkup(string.centeruser)
# 	frame.regfilestart$SetTooltipMarkup(string.regfilestart)
# 	label.seatsparameter$SetTooltipMarkup(string.seatsparameter)
# 	label.sigmalim$SetTooltipMarkup(string.sigmalim)
# 	frame.critical$SetTooltipMarkup(string.critical)
# 	label.outlier$SetTooltipMarkup(string.outlier)
# 	label.outlierspan$SetTooltipMarkup(string.outlierspan)
# 	label.outliermethod$SetTooltipMarkup(string.outliermethod)
# 	label.forecast_years$SetTooltipMarkup(string.forecast_years)
# 	label.backcast_years$SetTooltipMarkup(string.backcast_years)
# 	label.forecast_conf$SetTooltipMarkup(string.forecast_conf)
# 	label.aictest$SetTooltipMarkup(string.aictest)
# 	label.samode$SetTooltipMarkup(string.samode)
# 	label.seasonalma$SetTooltipMarkup(string.seasonalma)
# 	label.trendma$SetTooltipMarkup(string.trendma)
# 	label.x11calendarsigma$SetTooltipMarkup(string.x11calendarsigma)
# 	label.x11final$SetTooltipMarkup(string.x11final)
	
	
	
	#######################################################
	#Column plot-parameters
	#######################################################
	#manualoutlier table
	renderer.manualoutliertype$SetAlignment(0.5, 0.5)
	column.manualoutliertype$SetTitle("Type")
	column.manualoutliertype$PackStart(renderer.manualoutliertype)
	column.manualoutliertype$SetAlignment(0.5)
	column.manualoutliertype$AddAttribute(renderer.manualoutliertype, "text", 0)
	renderer.manualoutlieryear$SetAlignment(0.5, 0.5)
	column.manualoutlieryear$SetTitle("Year")
	column.manualoutlieryear$PackStart(renderer.manualoutlieryear)
	column.manualoutlieryear$SetAlignment(0.5)
	column.manualoutlieryear$AddAttribute(renderer.manualoutlieryear, "text", 1)
	renderer.manualoutlierperiod$SetAlignment(0.5, 0.5)
	column.manualoutlierperiod$SetTitle("period")
	column.manualoutlierperiod$PackStart(renderer.manualoutlierperiod)
	column.manualoutlierperiod$SetAlignment(0.5)
	column.manualoutlierperiod$AddAttribute(renderer.manualoutlierperiod, "text", 2)
	table.manualoutlier$AppendColumn(column.manualoutliertype)
	table.manualoutlier$AppendColumn(column.manualoutlieryear)
	table.manualoutlier$AppendColumn(column.manualoutlierperiod)
#	sapply(outlierlist,
#			function(string) {
#				## Add a new row to the model
#				iter <- tablemodel.manualoutlier$Append()$iter
#				tablemodel.manualoutlier$Set(iter, 0, string[1], 1, string[2], 2, string[3])
#			})
	update_outliertable()
	table.manualoutlier$SetModel(tablemodel.manualoutlier)
	#manualoutlier panel
	table.manualoutlier$SetSizeRequest(-1, 180)
	scroll.manualoutlier$Add(table.manualoutlier)
	panel.manualoutlier$AttachDefaults(scroll.manualoutlier, 1, 3, 1, 2)
	panel.manualoutlier$AttachDefaults(button.manualoutlierremove, 1, 3, 2, 3)
	panel.manualoutlier$AttachDefaults(label.manualoutliertype, 1, 2, 3, 4)
	combobox.manualoutliertype$AppendText("TC")
	combobox.manualoutliertype$AppendText("LS")
	combobox.manualoutliertype$AppendText("AO")
	panel.manualoutlier$AttachDefaults(combobox.manualoutliertype, 2, 3, 3, 4)
	panel.manualoutlier$AttachDefaults(label.manualoutlieryear, 1, 2, 4, 5)
	entry.manualoutlieryear$SetSizeRequest(20,-1)
	entry.manualoutlierperiod$SetSizeRequest(20,-1)
	panel.manualoutlier$AttachDefaults(entry.manualoutlieryear, 2, 3, 4, 5)
	panel.manualoutlier$AttachDefaults(label.manualoutlierperiod, 1, 2, 5, 6)
	panel.manualoutlier$AttachDefaults(entry.manualoutlierperiod, 2, 3, 5, 6)
	panel.manualoutlier$AttachDefaults(button.manualoutlieradd, 1, 3, 6, 7)
	panel.manualoutlier$AttachDefaults(button.manualoutlieraddclick, 1, 3, 7, 8)
	frame.manualoutlier$Add(panel.manualoutlier)
	panel.plotp$PackStart(frame.manualoutlier, expand=FALSE)
	gSignalConnect(button.manualoutlierremove, "released", f=manualoutlierhandler)
	gSignalConnect(button.manualoutlieraddclick, "released", f=manualoutlierhandler)
	gSignalConnect(button.manualoutlieradd, "released", f=manualoutlierhandler)
	
  #plot(...)
  panel.plotparams$PackStart(checkb.orig)  
  panel.plotparams$PackStart(checkb.sa)
  panel.plotparams$PackStart(checkb.trend)
  panel.plotparams$PackStart(checkb.logtransform)
  panel.plotparams$PackStart(checkb.showCI)
  #panel.plotparams$PackStart(checkb.showLine)
  panel.plotparams$PackStart(checkb.pointsOriginal)
  panel.plotparams$PackStart(checkb.showAllout)
  # panel.plotparams$PackStart(checkb.showAlloutLines)
  # panel.plotparams$PackStart(checkb.annComp)
  # panel.plotparams$PackStart(checkb.annCompTrend)
  # gSignalConnect(checkb.original, "toggled", f=function(...) update_notebook())
  gSignalConnect(checkb.orig, "toggled", f=function(...) update_notebook(onlyplot=TRUE))  
  gSignalConnect(checkb.sa, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
  gSignalConnect(checkb.trend, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
  gSignalConnect(checkb.logtransform, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
  gSignalConnect(checkb.showCI, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
  #gSignalConnect(checkb.showLine, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
  gSignalConnect(checkb.pointsOriginal, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
  gSignalConnect(checkb.showAllout, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
  # gSignalConnect(checkb.showAlloutLines, "toggled", f=function(...) update_notebook())
  # gSignalConnect(checkb.annComp, "toggled", f=function(...) update_notebook())
  # gSignalConnect(checkb.annCompTrend, "toggled", f=function(...) update_notebook())
  ##showout panel
  #	panel.showout$SetColSpacings(5)
  #	panel.showout$SetRowSpacings(5)
  entry.showoutyear$SetSizeRequest(15,-1)
  entry.showoutperiod$SetSizeRequest(15,-1)
  #combobox.showouttype$SetSizeRequest(15,-1)
  #combobox.showouttype$AppendText("AO")
  #combobox.showouttype$AppendText("LS")
  #combobox.showouttype$AppendText("TC")
  panel.showout$AttachDefaults(checkb.showout, 1, 2, 1, 2)
  panel.showout$AttachDefaults(label.showoutyear, 1, 2, 3, 4)
  panel.showout$AttachDefaults(entry.showoutyear, 2, 3, 3, 4)
  panel.showout$AttachDefaults(label.showoutperiod, 1, 2, 4, 5)
  panel.showout$AttachDefaults(entry.showoutperiod, 2, 3, 4, 5)
  #panel.showout$AttachDefaults(label.showouttype, 1, 2, 2, 3)
  #panel.showout$AttachDefaults(combobox.showouttype, 2, 3, 2, 3)
  alignment.showout <- gtkAlignment()
  alignment.showout$SetPadding(3,3,3,3)
  alignment.showout$Add(panel.showout)
  frame.showout$Add(alignment.showout)
  panel.plotparams$PackStart(frame.showout, padding=3)
  gSignalConnect(checkb.showout, "toggled", f=function(...){
    update_notebook(onlyplot=TRUE)
    toggle(c(entry.showoutyear), checkb.showout)
    #toggle(c(combobox.showouttype), checkb.showout)
    toggle(c(entry.showoutperiod), checkb.showout)})
  gSignalConnect(entry.showoutyear, "changed", f=function(...) update_notebook(onlyplot=TRUE))
  gSignalConnect(entry.showoutperiod, "changed", f=function(...) update_notebook(onlyplot=TRUE))
  #gSignalConnect(combobox.showouttype, "changed", f=function(...) update_notebook(onlyplot=TRUE))
  entry.showoutyear$SetSensitive(FALSE)
  entry.showoutperiod$SetSensitive(FALSE)
  #combobox.showouttype$SetSensitive(FALSE)
  checkb.showout$SetActive(FALSE)
  ##
  frame.plotparams$Add(panel.plotparams)
  panel.plotp$PackStart(frame.plotparams,expand=FALSE)
  
  #spectral frame
  radiob.spectralsa$SetActive(TRUE)
  panel.spectral$PackStart(radiob.spectralsa)
  panel.spectral$PackStart(radiob.spectraloriginal)
  panel.spectral$PackStart(radiob.spectralirregular)
  panel.spectral$PackStart(radiob.spectralresiduals)
  gSignalConnect(radiob.spectralsa, "toggled", f=function(...) if(radiob.spectralsa$GetActive()==TRUE)update_notebook(onlyplot=TRUE))
  gSignalConnect(radiob.spectraloriginal, "toggled", f=function(...) if(radiob.spectraloriginal$GetActive()==TRUE)update_notebook(onlyplot=TRUE))
  gSignalConnect(radiob.spectralirregular, "toggled", f=function(...) if(radiob.spectralirregular$GetActive()==TRUE)update_notebook(onlyplot=TRUE))
  gSignalConnect(radiob.spectralresiduals, "toggled", f=function(...) if(radiob.spectralresiduals$GetActive()==TRUE)update_notebook(onlyplot=TRUE))
  frame.spectral$Add(panel.spectral)
  panel.plotp$PackStart(frame.spectral, expand=FALSE)
  
  
  #seasonal factors plot
  radiob.rsdacfacf$SetActive(TRUE)
  panel.rsdacf$PackStart(radiob.rsdacfacf)
  panel.rsdacf$PackStart(radiob.rsdacfpacf)
  panel.rsdacf$PackStart(radiob.rsdacfacf2)
  gSignalConnect(radiob.rsdacfacf, "toggled", f=function(...) if(radiob.rsdacfacf$GetActive()==TRUE)update_notebook(onlyplot=TRUE))
  gSignalConnect(radiob.rsdacfpacf, "toggled", f=function(...) if(radiob.rsdacfpacf$GetActive()==TRUE)update_notebook(onlyplot=TRUE))
  gSignalConnect(radiob.rsdacfacf2, "toggled", f=function(...) if(radiob.rsdacfacf2$GetActive()==TRUE)update_notebook(onlyplot=TRUE))
  frame.rsdacf$Add(panel.rsdacf)
  panel.plotp$PackStart(frame.rsdacf, expand=FALSE)
  
  
	#summary parameter
	panel.summaryparameter$PackStart(checkb.fullSummary)
	panel.summaryparameter$PackStart(checkb.spectraldetail)
	panel.summaryparameter$PackStart(checkb.almostout)
	panel.summaryparameter$PackStart(checkb.rsdautocorr)
	panel.summaryparameter$PackStart(checkb.quality.stat)
	panel.summaryparameter$PackStart(checkb.likelihoodstat)
	panel.summaryparameter$PackStart(checkb.aape)
	panel.summaryparameter$PackStart(checkb.idrsdseas)
	panel.summaryparameter$PackStart(checkb.summaryslidingspans)
	panel.summaryparameter$PackStart(checkb.summaryhistory)
	panel.summaryparameter$PackStart(checkb.summaryidentify)	
	gSignalConnect(checkb.fullSummary, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.spectraldetail, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.almostout, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.rsdautocorr, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.quality.stat, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.likelihoodstat, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.aape, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.idrsdseas, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.summaryslidingspans, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.summaryhistory, "toggled", f=function(...) update_notebook())
	gSignalConnect(checkb.summaryidentify, "toggled", f=function(...) update_notebook())
	
	frame.summaryparameter$Add(panel.summaryparameter)
	panel.plotp$PackStart(frame.summaryparameter)
	
	#plotFbcast(...)
# panel.plotFbcastparams$PackStart(checkb.forecast)
# panel.plotFbcastparams$PackStart(checkb.backcast)
#  panel.plotparams$PackStart(checkb.showCI)
	##  panel.plotFbcastparams$PackStart(checkb.logtransform_fb)
#  panel.plotparams$PackStart(checkb.showLine)
#  panel.plotparams$PackStart(checkb.pointsOriginal)
	## gSignalConnect(checkb.forecast,"toggled", f=function(...) update_notebook(onlyplot=TRUE))
	## gSignalConnect(checkb.backcast, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
#  gSignalConnect(checkb.showCI, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
	##  gSignalConnect(checkb.logtransform_fb, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
#  gSignalConnect(checkb.showLine, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
#  gSignalConnect(checkb.pointsOriginal, "toggled", f=function(...) update_notebook(onlyplot=TRUE))
#  frame.plotFbcastparams$Add(panel.plotFbcastparams)
#  panel.plotp$PackStart(frame.plotFbcastparams, expand=FALSE)
	
	panel.gui$PackStart(panel.ts)
	panel.scrolledparams$SetPolicy("GTK_POLICY_NEVER","GTK_POLICY_ALWAYS")
	panel.scrolledparams$AddWithViewport(panel.params)
	panel.gui$PackStart(panel.scrolledparams)
	panel.scrolledplotparams$SetPolicy("GTK_POLICY_NEVER","GTK_POLICY_ALWAYS")
	panel.scrolledplotparams$AddWithViewport(panel.plotp)
	panel.gui$PackStart(panel.scrolledplotparams)
	
	#for each plot tab: frame + drawingarea
	frame.plot1$Add(area.plot1)
	frame.plot2$Add(area.plot2)
	frame.plot3$Add(area.plot3)
	frame.plot4$PackStart(area.plot4, expand=TRUE)
# frame.plot5$Add(area.plot5)
	
	#sliders for frame plot
	gSignalConnect(slider.plotmax, "value-changed", f=sliderhandler)
	gSignalConnect(slider.plotmin, "value-changed", f=sliderhandler)
	gSignalConnect(slider.plotmin, "format-value", f=sliderformat)
	gSignalConnect(slider.plotmax, "format-value", f=sliderformat)
	frame.plot4$PackStart(slider.plotmin, expand=FALSE)
	frame.plot4$PackStart(slider.plotmax, expand=FALSE)
	
	#contextmenu for plots
	menu.contextplotall$Attach(menuitem.saveaspdf, 0, 1, 0, 1)
	gSignalConnect(menuitem.saveaspdf, "activate", f=menuhandler)
	gSignalConnect(area.plot1, "button_release_event", f=mousehandlerdrawing)
	gSignalConnect(area.plot2, "button_release_event", f=mousehandlerdrawing)
	gSignalConnect(area.plot3, "button_release_event", f=mousehandlerdrawing)
	gSignalConnect(area.plot4, "button_release_event", f=mousehandlerdrawing)
  gSignalConnect(area.plot4, "button_press_event", f=mousehandlerdrawing)
  
	menu.contextplotwithoutlier$Attach(menuitem.saveaspdfwithoutlier, 0, 1, 0, 1)
	menu.contextplotwithoutlier$Attach(menuitem.addTC, 0, 1, 1, 2)
	menu.contextplotwithoutlier$Attach(menuitem.addLS, 0, 1, 2, 3)
	menu.contextplotwithoutlier$Attach(menuitem.addAO, 0, 1, 3, 4)
	gSignalConnect(menuitem.addTC, "activate", f=menuhandler)
	gSignalConnect(menuitem.addAO, "activate", f=menuhandler)
	gSignalConnect(menuitem.addLS, "activate", f=menuhandler)
	gSignalConnect(menuitem.saveaspdfwithoutlier, "activate", f=menuhandler)
	
	
	#summary tab setup
# buffer.summary$SetText(paste(capture.output(summary(x12)), collapse="\n"))
#	make_summary(object)
#	textview.summary$SetBuffer(buffer.summary)
#	textview.summary$setEditable(FALSE)
#	frame.summary$add(textview.summary)
#	table.summary <- setupSummarytable(table.summary, s)
#	i <- 0
#	sumnames <- names(getMethod("summary","x12Batch")(object, print=FALSE))
#	sumnames[1] <- "Value"
#	for(s in sumnames){
#		renderer <- gtkCellRendererText()
#		column <- gtkTreeViewColumn()
#		renderer$SetAlignment(0.5, 0.5)
#		column$SetTitle(s)
#		column$PackStart(renderer)
#		column$SetAlignment(0.5)
#		column$SetExpand(TRUE)
#		column$AddAttribute(renderer, "text", i)
#		if(i==0)column$AddAttribute(renderer, "background", length(sumnames))
#		else column$AddAttribute(renderer, "background", length(sumnames)+1)
#		i <- i + 1
#		table.summary$AppendColumn(column)
#		if(class(columns.summary)!="list")columns.summary <- list(column)
#		else columns.summary <- append(columns.summary, column)
#	}
	setup_summarytable(object)
	table.summary$SetModel(model.summary)
	table.summary$SetHeadersVisible(FALSE)
	frame.summary$add(table.summary)
	
	#summarytotal tab setup
	make_summary(object)
#	buffer.summarytotal$SetText(paste(capture.output(getMethod("summary","x12Batch")(object)), collapse="\n"))
  
  pF <- pangoFontDescriptionNew()
  pangoFontDescriptionSetFamily(pF,"Monospace")
  textview.summarytotal$modifyFont(pF)
	textview.summarytotal$SetBuffer(buffer.summarytotal)
	textview.summarytotal$setEditable(FALSE)
	frame.summarytotal$add(textview.summarytotal)
	
  notebook.plot$AppendPage(frame.plot4, tab.label=gtkLabel("plot"))
  notebook.plot$AppendPage(frame.plot2, tab.label=gtkLabel("spectral"))
	notebook.plot$AppendPage(frame.plot1, tab.label=gtkLabel("autocorrelations"))
	notebook.plot$AppendPage(frame.plot3, tab.label=gtkLabel("seasonal factors"))
# notebook.plot$AppendPage(frame.plot5, tab.label=gtkLabel("Summary Total"))
	notebook.plot$AppendPage(frame.summarytotal, tab.label=gtkLabel("summary text"))
	notebook.plot$AppendPage(frame.summary, tab.label=gtkLabel("summary table"))
	
	gSignalConnect(menuitem.x12update, "activate", f=menuhandler)
	gSignalConnect(menuitem.x12loadp, "activate", f=menuhandler)
	gSignalConnect(menuitem.x12save, "activate", f=menuhandler)
	gSignalConnect(menuitem.x12savep, "activate", f=menuhandler)
	gSignalConnect(menuitem.x12load, "activate", f=menuhandler)
	gSignalConnect(menuitem.expplotaspdf, "activate", f=menuhandler)
	gSignalConnect(menuitem.expplotaspng, "activate", f=menuhandler)
	gSignalConnect(menuitem.expsummarycsv, "activate", f=menuhandler)
	gSignalConnect(menuitem.expsummaryclipboard, "activate", f=menuhandler)
  gSignalConnect(menuitem.path, "activate", f=menuhandler)
	accgroup <- gtkAccelGroupNew()
	window.main$AddAccelGroup(accgroup)
	menuitem.x12update$AddAccelerator("activate", accgroup, GDK_U, "GDK_CONTROL_MASK", "GTK_ACCEL_VISIBLE")
	menu.export$Append(menuitem.expplotaspdf)
	menu.export$Append(menuitem.expplotaspng)
	menu.export$Append(menuitem.expsummarycsv)
	menu.export$Append(menuitem.expsummaryclipboard)
	menuitem.export$SetSubmenu(menu.export)
	menu.x12$Append(menuitem.x12update)
	menu.x12$Append(menuitem.x12loadp)
	menu.x12$Append(menuitem.x12savep)
	menu.x12$Append(menuitem.x12load)
	menu.x12$Append(menuitem.x12save)
  menu.x12$Append(menuitem.path)
	menuitem.x12$SetSubmenu(menu.x12)
	menubar.main$Append(menuitem.x12)
	menubar.main$Append(menuitem.export)
	
	panel.main$Add(panel.gui)
	panel.main$Add(notebook.plot)
	panel.window$PackStart(menubar.main, expand=FALSE)
	panel.window$PackStart(panel.main, expand=TRUE)
	panel.window$PackStart(statusbar, expand=FALSE)
	window.main$Add(panel.window)
	
	asCairoDevice(area.plot1)
	Sys.sleep(.1)
	
	#preloading the plots of the notebook tab
  notebook.plot$SetCurrentPage(0)
	notebook.plot$SetCurrentPage(1)
	notebook.plot$SetCurrentPage(2)
	notebook.plot$SetCurrentPage(3)
# notebook.plot$SetCurrentPage(4)
	notebook.plot$SetCurrentPage(0)
	notebook.plot$SetCurrentPage(5)
	
	gSignalConnect(notebook.plot, "switch-page", f=notebookhandler)
	window.main$Resize(1100,700)
	window.main$Show()
	capture.output(read_x12(object, c(1)))
	window.main$SetFocus(table.ts)
  gSignalHandlerBlock(table.ts$GetSelection(), handler.tstable)
	table.ts$GetSelection()$SelectPath(gtkTreePathNewFirst())
  gSignalHandlerUnblock(table.ts$GetSelection(), handler.tstable)
	#table.ts$SetCursor(gtkTreePathNewFirst())
	gSignalConnect(window.main, "destroy", f=function(...){gtkMainQuit()})
	status_print("Programm started!")
	gtkMain()
	invisible(object)
}







