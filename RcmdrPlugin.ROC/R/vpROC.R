
#=========================================================================================================================================

fncvpROC <- function(){
  varPosListn <- function(vars, var){
    if (is.null(var)) return(NULL)
    if (any(!var %in% vars)) NULL
    else apply(outer(var, vars, "=="), 1, which) - 1
  }
  #require(pROC)
  #Daniel
  defaults <- list(# data tab
                   initial.prediction = NULL, initial.label = NULL, 
                   initial.narm = 1, initial.percent = 0, initial.direction = "auto", 
                   # smoothing tab
                   initial.smooth = 0, initial.smoothingmethod = "binormal", 
                   initial.smoothinbandwidth = "nrd0", initial.bandwidthnumeric = "", 
                   initial.bandwidthadjustment = "1", initial.bandwidthwindow = "gaussian",
                   initial.distributioncontrols = "normal", initial.distributioncases = "normal",
                   # ci tab
                   initial.cicompute = 1, initial.cilevel = "0.95", initial.cimethod = "bootstrap", 
                   initial.cibootn = "2000", initial.cibootstratified = 0, 
                   initial.citype = "se", initial.cithresholds = "local maximas", 
                   initial.civalues = "seq(0, 1, 0.05)", initial.ciplottype = "shape", 
                   # auc tab
                   initial.auc = 1, initial.partialauc = 0, 
                   initial.partialfrom = 0, initial.partialto = 1,
                   initial.partialfocus = "specificity", initial.partialcorrect = 0,
                   # plot tab
                   initial.plot = 1, initial.add = 0, 
                   initial.printauc = 0, initial.aucpolygon = 0, initial.maxaucpolygon = 0, 
                   initial.grid = 0, initial.identity = 1, initial.ciplot = 0, initial.values = 0, 
                   initial.printthresrb = "no", initial.customthres = "c(0.5, 1, 10, 100)",
                   initial.colorroc=palette()[1], #initial.colorroc="#1c61b6", 
                   initial.ltyroc="solid", 
                   initial.xlab=gettextRcmdr("<auto>"), initial.ylab=gettextRcmdr("<auto>"), 
                   initial.main=gettextRcmdr("<auto>"),
                   initial.tab=0) # tab
  dialog.values <- getDialog("pROC", defaults)

  initializeDialog(title=gettext("Plot ROC curve", domain="R-RcmdrPlugin.ROC"), use.tabs=TRUE, tabs=c("dataTab", "smoothingTab", "aucTab", "ciTab", "optionsTab")) # tab
  
  #Daniel
  #General/data tab:
  generalFrame <- tkframe(dataTab)# tab
  generaldataFrame <- ttklabelframe(generalFrame, text = gettext("Data", domain="R-RcmdrPlugin.ROC"))
  predictionBox <- variableListBox(generaldataFrame, Numeric(), title=gettext("Predictions variable (pick one)", domain="R-RcmdrPlugin.ROC"),
                                   initialSelection=varPosn(dialog.values$initial.prediction, "numeric"))
  labelBox <- variableListBox(generaldataFrame, Factors(), title=gettext("Outcome variable (pick one)", domain="R-RcmdrPlugin.ROC"),
                              initialSelection=varPosn(dialog.values$initial.label, "factor"))
 
  checkBoxes(window = generalFrame, 
             frame = "dataoptionsFrame", # tab
             boxes = c("narm", "percent"), 
             initialValues = c(dialog.values$initial.narm, dialog.values$initial.percent),
             labels = c(gettext("Remove NAs", domain="R-RcmdrPlugin.ROC"), gettext("Show/input % instead of 0-1", domain="R-RcmdrPlugin.ROC")), 
             title = gettext("Options", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
 
  radioButtons(dataoptionsFrame, 
               name="directionrb", 
               buttons=c("auto", "gt", "lt"), 
               values=c("auto", ">", "<"),
               labels=gettextRcmdr(c("auto", "Control > cases", "Control <= cases")), 
               title=gettext("Direction", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.direction)  
  
  # Smoothing tab:
  smoothingFrame <- tkframe(smoothingTab)# tab
  smoothingleftpaneFrame <- tkframe(smoothingFrame)#ttklabelframe(smoothingFrame, text = "")
  smoothinggeneralFrame <- ttklabelframe(smoothingleftpaneFrame, text = gettext("General", domain="R-RcmdrPlugin.ROC"))
  smoothingdensityFrame <- ttklabelframe(smoothingleftpaneFrame, text = gettext("Density options", domain="R-RcmdrPlugin.ROC"))
  smoothingdistributionFrame <- ttklabelframe(smoothingFrame, text = gettext("Distributions options", domain="R-RcmdrPlugin.ROC"))
  
  radioButtons(smoothinggeneralFrame, 
               name="smoothingmethodrb", 
               buttons=c("binormal", "density", "fitdistr", "logcondens", "logcondens.smooth"), 
               values=c("binormal", "density", "fitdistr", "logcondens", "logcondens.smooth"),
               labels=gettextRcmdr(c("binormal", "density", "fit distribution", "logcondens", "logcondens.smooth")), 
               title=gettext("Smoothing method", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.smoothingmethod)  

  radioButtons(smoothingdensityFrame, 
               name="smoothinbandwidthrb", 
               buttons=c("nrd0", "nrd", "ucv", "bcv", "SJ", "numeric"), 
               values=c("nrd0", "nrd", "ucv", "bcv", "SJ", "numeric"),
               labels=gettextRcmdr(c("nrd0", "nrd", "ucv", "bcv", "SJ", "<numeric>")), 
               title=gettext("Bandwidth", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.smoothinbandwidth)  
  
  bandwidthnumericVar <- tclVar(dialog.values$initial.bandwidthnumeric) # tab
  bandwidthnumericEntry <- ttkentry(smoothingdensityFrame, width = "25", textvariable = bandwidthnumericVar)# tab
  bandwidthnumericScroll <- ttkscrollbar(smoothingdensityFrame, orient = "horizontal",
                                    command = function(...) tkxview(bandwidthnumericEntry, ...))
  tkconfigure(bandwidthnumericEntry, xscrollcommand = function(...) tkset(bandwidthnumericScroll,
                                                                     ...))
  tkbind(bandwidthnumericEntry, "<FocusIn>", function() tkselection.clear(bandwidthnumericEntry))
  
  bandwidthadjustmentVar <- tclVar(dialog.values$initial.bandwidthadjustment) # tab
  bandwidthadjustmentEntry <- ttkentry(smoothingdensityFrame, width = "25", textvariable = bandwidthadjustmentVar)# tab
 
  radioButtons(smoothingdensityFrame, # kernel for density 
               name="bandwidthwindowrb", 
               buttons=c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"), 
               values=c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"),
               labels=gettextRcmdr(c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")), 
               title=gettext("Kernel", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.bandwidthwindow) 
  
  radioButtons(smoothingdistributionFrame, #"beta", "chi-squared", "f", "geometric", "negative binomial",, "poisson", "t"
               name="distributioncontrolsrb", 
               buttons=c("normal", "lognormal", "logistic", "exponential", "weibull", "gamma", "cauchy"), 
               values=c("normal", "lognormal", "logistic", "exponential", "weibull", "gamma", "cauchy"),
               labels=gettextRcmdr(c("normal", "lognormal", "logistic", "exponential", "weibull", "gamma", "cauchy")), 
               title=gettext("Distribution of controls", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.distributioncontrols) 
  radioButtons(smoothingdistributionFrame, #"beta", "chi-squared", "f", "geometric", "negative binomial",, "poisson", "t"
               name="distributioncasesrb", 
               buttons=c("normal", "lognormal", "logistic", "exponential", "weibull", "gamma", "cauchy"), 
               values=c("normal", "lognormal", "logistic", "exponential", "weibull", "gamma", "cauchy"),
               labels=gettextRcmdr(c("normal", "lognormal", "logistic", "exponential", "weibull", "gamma", "cauchy")), 
               title=gettext("Distribution of cases", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.distributioncases) 
  
  # CI tab:
  ciFrame <- tkframe(ciTab)# tab
  #cigeneralFrame <- ttklabelframe(ciFrame, text = gettext("General", domain="R-RcmdrPlugin.ROC"))
  #cibootstrapFrame <- ttklabelframe(ciFrame, text = gettext("Bootstrap options", domain="R-RcmdrPlugin.ROC"))
  checkBoxes(window = ciFrame, frame = "cibootstrapFrame", # tab
             boxes = c("cibootstratified"), initialValues = c(
               dialog.values$initial.cibootstratified
             ),labels = gettextRcmdr(c(
               "Stratified")), title = gettext("Bootstrap options", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  
  checkBoxes(window = ciFrame, frame = "cigeneralFrame", # tab
             boxes = c("cicompute"), initialValues = c(
               dialog.values$initial.cicompute
             ),labels = gettextRcmdr(c(
               "Compute Confidence Interval (CI)")), title = gettext("General", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  
  cilevelVar <- tclVar(dialog.values$initial.cilevel) # tab
  cilevelEntry <- ttkentry(cigeneralFrame, width = "25", textvariable = cilevelVar)# tab

  radioButtons(cigeneralFrame, name="cimethodrb", buttons=c("delong", "bootstrap", "auto"), values=c("delong", "bootstrap", "auto"),
               labels=gettextRcmdr(c("delong", "bootstrap", "auto")), title=gettext("Method", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.cimethod)  

  radioButtons(cigeneralFrame, name="cityperb", buttons=c("auc", "se", "sp", "thresholds"), values=c("auc", "se", "sp", "thresholds"),
               labels=gettextRcmdr(c("auc", "se", "sp", "thresholds")), title=gettext("Type of CI", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.citype)
  
  radioButtons(cigeneralFrame, name="cithresholdsrb", buttons=c("all", "localmaximas", "custom"), values=c("all", "local maximas", "custom"),
               labels=gettextRcmdr(c("all", "local maximas", "<custom>")), title=gettext("Thresholds", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.cithresholds)
  
  civaluesVar <- tclVar(dialog.values$initial.civalues) # tab
  civaluesEntry <- ttkentry(cigeneralFrame, width = "25", textvariable = civaluesVar)# tab
  civaluesScroll <- ttkscrollbar(cigeneralFrame, orient = "horizontal",
                                    command = function(...) tkxview(civaluesEntry, ...))
  tkconfigure(civaluesEntry, xscrollcommand = function(...) tkset(civaluesScroll,
                                                                     ...))
  tkbind(civaluesEntry, "<FocusIn>", function() tkselection.clear(civaluesEntry))


  cibootnVar <- tclVar(dialog.values$initial.cibootn) # tab
  cibootnEntry <- ttkentry(cibootstrapFrame, width = "5", textvariable = cibootnVar)# tab
  tkgrid(labelRcmdr(cibootstrapFrame, text = gettext("Confidence level number of replicates", domain="R-RcmdrPlugin.ROC")), cibootnEntry, sticky = "ew", padx=6)
  
  # AUC tab:
  aucFrame <- tkframe(aucTab)# tab
  #generalaucFrame <- ttklabelframe(aucFrame, text = gettext("General", domain="R-RcmdrPlugin.ROC"))
  #partialaucFrame <- ttklabelframe(aucFrame, text = gettext("Partial AUC", domain="R-RcmdrPlugin.ROC"))

  checkBoxes(window = aucFrame, frame = "generalaucFrame", # tab
             boxes = c("auc", "partialauc"), initialValues = c(
               dialog.values$initial.auc, dialog.values$initial.partialauc
               ),labels = gettextRcmdr(c(
                 "Compute Area Under Curve (AUC)", "Compute partial AUC")), title = gettext("General", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  
    checkBoxes(window = aucFrame, frame = "partialaucFrame", # tab             !!!!!!!!!!  NU merge inca!!!
               boxes = c("partialcorrect"), initialValues = c(
                 dialog.values$initial.partialcorrect),labels = gettextRcmdr(c(
                   "Correct partial AUC")), title = gettext("Partial AUC", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  
  partialfromVar <- tclVar(dialog.values$initial.partialfrom) # tab
  partialfromEntry <- ttkentry(partialaucFrame, width = "25", textvariable = partialfromVar)# tab
  tkgrid(labelRcmdr(partialaucFrame, text = gettext("From:", domain="R-RcmdrPlugin.ROC")), partialfromEntry, sticky = "ew", padx=6)
  partialtoVar <- tclVar(dialog.values$initial.partialto) # tab
  partialtoEntry <- ttkentry(partialaucFrame, width = "25", textvariable = partialtoVar)# tab
  tkgrid(labelRcmdr(partialaucFrame, text = gettext("To:", domain="R-RcmdrPlugin.ROC")), partialtoEntry, sticky = "ew", padx=6)
  radioButtons(partialaucFrame, name="partialfocus", buttons=c("specificity", "sensitivity"), values=c("specificity", "sensitivity"),
               labels=c(gettext("specificity", domain="R-RcmdrPlugin.ROC"), gettext("sensitivity", domain="R-RcmdrPlugin.ROC")), title=gettext("Focus", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.partialfocus)    

  
  # Plot tab:
  optionsParFrame <- tkframe(optionsTab)# tab
  optFrame <- ttklabelframe(optionsParFrame, text = gettext("Plot Options", domain="R-RcmdrPlugin.ROC"))
  parFrame <- ttklabelframe(optionsParFrame, text = gettext("Plot Labels", domain="R-RcmdrPlugin.ROC"))
  #paletteFrame <- tkframe(optionsTab)# tab

  checkBoxes(window = optFrame, frame = "optionsFrame", # tab
             boxes = c("plot", "add", "smooth", "grid","identity","ciplot","values"), initialValues = c(
               dialog.values$initial.plot, dialog.values$initial.add, dialog.values$initial.smooth, 
               dialog.values$initial.grid, dialog.values$initial.identity, dialog.values$initial.ciplot, dialog.values$initial.values),labels = gettextRcmdr(c(
                 "Plot", "Add curve to existing plot", "Smooth","Display grid","Display identity line",
                 "Display confidence interval","Display values (Se, Sp, Thresholds)")), title = gettext("General", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  checkBoxes(window = optFrame, frame = "aucpolygonFrame", # tab
             boxes = c("aucpolygon", "maxaucpolygon"), initialValues = c(
               dialog.values$initial.aucpolygon, dialog.values$initial.maxaucpolygon),labels = gettextRcmdr(c(
                 "Polygon of AUC", "Polygon of maximal AUC")), title = gettext("Display area as polygon", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  checkBoxes(window = optFrame, frame = "informationFrame", # tab
             boxes = c("printauc"), initialValues = c(
               dialog.values$initial.printauc),labels = gettextRcmdr(c(
                 "AUC")), title = gettext("Display information on plot", domain="R-RcmdrPlugin.ROC"), ttk=TRUE)
  radioButtons(informationFrame, name="printthresrb", buttons=c("no", "best", "all", "localmaximas", "customthres"), values=c("no", "best", "all", "local maximas", "customthres"),
               labels=gettextRcmdr(c("no", "best: max(sum(Se + Sp))", "all", "local maximas", "<custom>")), title=gettext("Display threshold(s)", domain="R-RcmdrPlugin.ROC"),
               initialValue = dialog.values$initial.printthresrb)  

  customthresVar <- tclVar(dialog.values$initial.customthres) # tab
  customthresEntry <- ttkentry(informationFrame, width = "25", textvariable = customthresVar)# tab
  customthresScroll <- ttkscrollbar(informationFrame, orient = "horizontal",
                                command = function(...) tkxview(customthresEntry, ...))
  tkconfigure(customthresEntry, xscrollcommand = function(...) tkset(customthresScroll,
                                                                 ...))
  tkbind(customthresEntry, "<FocusIn>", function() tkselection.clear(customthresEntry))

  
  xlabVar <- tclVar(dialog.values$initial.xlab) # tab
  ylabVar <- tclVar(dialog.values$initial.ylab)
  mainVar <- tclVar(dialog.values$initial.main)
  xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
  xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                             command = function(...) tkxview(xlabEntry, ...))
  tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
                                                              ...))
  tkbind(xlabEntry, "<FocusIn>", function() tkselection.clear(xlabEntry))
  tkgrid(labelRcmdr(parFrame, text = gettext("x-axis label", domain="R-RcmdrPlugin.ROC")), xlabEntry, sticky = "ew", padx=6)
  tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
  ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
  ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                             command = function(...) tkxview(ylabEntry, ...))
  tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
                                                              ...))
  tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
  tkgrid(labelRcmdr(parFrame, text=""), ylabScroll, sticky = "ew", padx=6)
  mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
  mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                             command = function(...) tkxview(mainEntry, ...))
  tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
                                                              ...))
  tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
  tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)

  radioButtons(parFrame, name="ciplottyperb", buttons=c("shape", "bars"), values=c("shape", "bars"),
             labels=c(gettext("shape", domain="R-RcmdrPlugin.ROC"), gettext("bars", domain="R-RcmdrPlugin.ROC")), title=gettext("CI plot type", domain="R-RcmdrPlugin.ROC"),
             initialValue = dialog.values$initial.ciplottype) 

  colorrocBox <- variableListBox(parFrame, palette(), title=gettext("Color of ROC (from Palette)", domain="R-RcmdrPlugin.ROC"),
                                 initialSelection=varPosListn(palette(), dialog.values$initial.colorroc))
  
  ltys <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "blank")
  ltyrocBox <- variableListBox(parFrame, ltys, title=gettext("Line type of ROC (from Palette)", domain="R-RcmdrPlugin.ROC"),
                               initialSelection=varPosListn(ltys, dialog.values$initial.ltyroc))

  
  onOK <- function(){
    tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1 # tab            de modificat!!!!!!!!!!!
    #Daniel
    #general tab
    prediction <- getSelection(predictionBox)
    label <- getSelection(labelBox)    
    narm <- as.character("1" == tclvalue(narmVariable)) 
    percent <- as.character("1" == tclvalue(percentVariable)) 
    direction <- as.character(tclvalue(directionrbVariable)) 
    
    #smoothing tab
    smoothingmethod <- as.character(tclvalue(smoothingmethodrbVariable)) 
    smoothinbandwidth <- as.character(tclvalue(smoothinbandwidthrbVariable)) 
    bandwidthnumeric <- as.character(tclvalue(bandwidthnumericVar))
    bandwidthadjustment <- as.character(tclvalue(bandwidthadjustmentVar))    
    bandwidthwindow <- as.character(tclvalue(bandwidthwindowrbVariable)) 
    distributioncases <- as.character(tclvalue(distributioncasesrbVariable)) 
    distributioncontrols <- as.character(tclvalue(distributioncontrolsrbVariable)) 
    
    #ci tab
    cicompute <- as.character("1" == tclvalue(cicomputeVariable)) 
    cilevel <- as.numeric(as.character(tclvalue(cilevelVar)))
    cimethod <- as.character(tclvalue(cimethodrbVariable)) 
    citype <- as.character(tclvalue(cityperbVariable)) 
    cithresholds <- as.character(tclvalue(cithresholdsrbVariable)) 
    civalues <- as.character(tclvalue(civaluesVar)) 
    cibootn <- as.integer(as.character(tclvalue(cibootnVar)))
    cibootstratified <- as.character("1" == tclvalue(cibootstratifiedVariable))          
    
    #auc tab
    auc <- as.character("1" == tclvalue(aucVariable))
    partialauc <- as.character("1" == tclvalue(partialaucVariable))
    partialfrom <- as.character(tclvalue(partialfromVar))
    partialto <- as.character(tclvalue(partialtoVar))
    partialfocus <- as.character(tclvalue(partialfocusVariable))
    partialcorrect <- as.character("1" == tclvalue(partialcorrectVariable))
    
    #plot tab
    add <- as.character("1" == tclvalue(addVariable))
    plot <- as.character("1" == tclvalue(plotVariable))
    smooth <- as.character("1" == tclvalue(smoothVariable))
    printauc <- as.character("1" == tclvalue(printaucVariable))
    aucpolygon <- as.character("1" == tclvalue(aucpolygonVariable))
    maxaucpolygon <- as.character("1" == tclvalue(maxaucpolygonVariable))
    grid <- as.character("1" == tclvalue(gridVariable))
    identity <- as.character("1" == tclvalue(identityVariable))
    ciplot <- as.character("1" == tclvalue(ciplotVariable))
    values <- as.character("1" == tclvalue(valuesVariable))
    
    printthresrb <- as.character(tclvalue(printthresrbVariable))
    customthres <- as.character(tclvalue(customthresVar))
    
    xlab <- trim.blanks(tclvalue(xlabVar))
    xlab <- if (xlab == gettextRcmdr("<auto>"))
      ""
    else paste(", xlab=\"", xlab, "\"", sep = "")
    ylab <- trim.blanks(tclvalue(ylabVar))
    ylab <- if (ylab == gettextRcmdr("<auto>"))
      ""
    else paste(", ylab=\"", ylab, "\"", sep = "")
    main <- trim.blanks(tclvalue(mainVar))
    main <- if (main == gettextRcmdr("<auto>"))
      ""
    else paste(", main=\"", main, "\"", sep = "")
 
    ciplottype <- as.character(tclvalue(ciplottyperbVariable))
    colorroc <- getSelection(colorrocBox)
    convert <- function (color){
      f=col2rgb(color)
      rgb(f[1],f[2],f[3],maxColorValue=255)
    }
    if(substr(colorroc,1,1) != "#") colorroc <- convert(colorroc)
    ltyroc <- as.character(getSelection(ltyrocBox))
    
    putDialog ("pROC", list(# data tab
                            initial.prediction = prediction, initial.label = label,
                            initial.narm = tclvalue(narmVariable), initial.percent = tclvalue(percentVariable), 
                            initial.direction = as.character(tclvalue(directionrbVariable)),
                            # smoothing tab
                            initial.smooth = tclvalue(smoothVariable), initial.smoothingmethod = tclvalue(smoothingmethodrbVariable),
                            initial.smoothinbandwidth = tclvalue(smoothinbandwidthrbVariable), initial.bandwidthnumeric = tclvalue(bandwidthnumericVar),
                            initial.bandwidthadjustment = "1", initial.bandwidthwindow = tclvalue(bandwidthwindowrbVariable),
                            initial.distributioncontrols = tclvalue(distributioncontrolsrbVariable), initial.distributioncases = tclvalue(distributioncasesrbVariable),
                            # ci tab
                            initial.cicompute = tclvalue(cicomputeVariable), initial.cilevel = tclvalue(cilevelVar), initial.cimethod = tclvalue(cimethodrbVariable), 
                            initial.cibootn = tclvalue(cibootnVar), initial.cibootstratified = tclvalue(cibootstratifiedVariable), 
                            initial.citype = tclvalue(cityperbVariable), initial.cithresholds = tclvalue(cithresholdsrbVariable), 
                            initial.civalues = tclvalue(civaluesVar), initial.ciplottype = tclvalue(ciplottyperbVariable),   
                            # auc tab
                            initial.auc = tclvalue(aucVariable), initial.partialauc = tclvalue(partialaucVariable), 
                            initial.partialfrom = tclvalue(partialfromVar), initial.partialto = tclvalue(partialtoVar),
                            initial.partialfocus = tclvalue(partialfocusVariable), initial.partialcorrect = tclvalue(partialcorrectVariable),
                            # plot tab
                            initial.plot = tclvalue(plotVariable), initial.add = tclvalue(addVariable), 
                            initial.printauc = tclvalue(printaucVariable), initial.aucpolygon = tclvalue(aucpolygonVariable),  initial.maxaucpolygon = tclvalue(maxaucpolygonVariable),
                            initial.grid = tclvalue(gridVariable), initial.identity = tclvalue(identityVariable), 
                            initial.ciplot = tclvalue(ciplotVariable), initial.values = tclvalue(valuesVariable), initial.printthresrb = tclvalue(printthresrbVariable), initial.customthres = as.character(tclvalue(customthresVar)),
                            initial.colorroc = getSelection(colorrocBox), 
                            initial.ltyroc = getSelection(ltyrocBox), 
                            initial.xlab=tclvalue(xlabVar), initial.ylab=tclvalue(ylabVar), 
                            initial.main=tclvalue(mainVar),
                            initial.tab=tab)) # tab
closeDialog()

# Checking input ============================================
# data tab
if (0 == length(prediction)) {
  errorCondition(recall=fncvpROC, message=gettext("You must select a prediction variable.", domain="R-RcmdrPlugin.ROC"))
  return()
}
if (0 == length(label)) {
  errorCondition(recall=fncvpROC, message=gettext("No outcome variable selected.", domain="R-RcmdrPlugin.ROC"))
  return()
}
if (percent == "TRUE") {
  percentupper = 100
} else {
  percentupper = 1  
}
# ci tab
if (cicompute == "TRUE") {
  if (0 == length(cilevel)) {
    errorCondition(recall=fncvpROC, message=gettext("You must set a confidence interval level.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
  cilevel = as.numeric(cilevel)
  if ((cilevel < 0) || (cilevel > 1)) {
    errorCondition(recall=fncvpROC, message=gettext("Confidence interval level outside of range.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
  if (0 == length(cibootn)) {
    errorCondition(recall=fncvpROC, message=gettext("You must set a confidence interval number of replicates.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
  if (cibootn < 0) {
    errorCondition(recall=fncvpROC, message=gettext("Confidence interval number of replicates should be a pozitive number.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
}
# auc tab
if (partialauc == "TRUE") {
  if (0 == length(partialto)) {
    errorCondition(recall=fncvpROC, message=gettext("You must set a partial AUC 'to' limit.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
  partialto = as.numeric(partialto)
  partialfrom = as.numeric(partialfrom)
  if ((partialto < 0) | (partialto > percentupper)) {
    errorCondition(recall=fncvpROC, message=gettext("Partial AUC 'to' limit outside of range.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
  if (0 == length(partialfrom)) {
    errorCondition(recall=fncvpROC, message=gettext("You must set a partial AUC 'from' limit.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
  if ((partialfrom < 0) | (partialfrom > percentupper)) {
    errorCondition(recall=fncvpROC, message=gettext("Partial AUC 'from' limit outside of range.", domain="R-RcmdrPlugin.ROC"))
    return()
  }
  if ((max(c(partialfrom, partialto)) <= 1) & (percent=="TRUE")) {
    Message(message="Maybe you didn't specified well the values, you probably wanted to set the values between 0-100 instead of between 0-1, since percent is checked", type="warning")    
  }
  if ((max(c(partialfrom, partialto)) > 1) & (percent=="FALSE")) {
    Message(message="Maybe you didn't specified well the values, you probably wanted to set the values between 0-1 instead of between 0-100, since percent is not checked", type="warning")    
  }
  
}
# plot tab
if ((printthresrb == "custom") & (0 == length(customthres))) {
  errorCondition(recall=fncvpROC, message=gettext("Custom threshold should not be empty.", domain="R-RcmdrPlugin.ROC"))
  return()
}

# transformations
.activeDataSet <- ActiveDataSet()
if (printthresrb == "customthres") {
  threshold = customthres
} else {
  threshold = paste("'", printthresrb, "'", sep="") 
}
if (partialauc == "TRUE") {
  partialauc = paste("c(", partialfrom, ", ", partialto, ")", sep="") 
}

#Daniel 
# command <- paste("roc.obj <- pROC::roc(", label, " ~ ", prediction, ", data=", .activeDataSet, ", na.rm=", narm, ", percent=", percent, ", direction='", direction, "'",  
#                  ", auc=", auc, ", partial.auc=", partialauc, ", partial.auc.focus='", partialfocus, "'", ", partial.auc.correct=", partialcorrect, 
#                  ", plot=", plot, ", add=", add,
#                  ", print.auc=", printauc, ", auc.polygon=", aucpolygon, ", max.auc.polygon=", maxaucpolygon, 
#                  ", grid=", grid, ", identity=", identity, 
#                  ", print.thres=", threshold, xlab, ylab, main, ")", sep = "")
# doItAndPrint(command)
command <- paste("roc.obj <- pROC::roc(", label, " ~ ", prediction, ", data=", .activeDataSet, ", na.rm=", narm, ", percent=", percent, ", direction='", direction, "'",  
                 ", partial.auc=", partialauc, ", partial.auc.focus='", partialfocus, "'", ", partial.auc.correct=", partialcorrect, 
                 ", auc=", auc, ", plot=FALSE, ci=TRUE, of='auc', conf.level=", cilevel, ", ci.method='", cimethod,"', boot.n=", cibootn, ", boot.stratified=", cibootstratified,")", sep = "")
doItAndPrint(command)
if (plot == "TRUE") {
  command <- paste("plot(roc.obj, add=", add,
                   ", print.auc=", printauc, ", auc.polygon=", aucpolygon, ", max.auc.polygon=", maxaucpolygon, 
                   ", print.auc.x=ifelse(roc.obj$percent, 50, .5), print.auc.y=ifelse(roc.obj$percent, 45, .45), print.auc.pattern='AUC: %.2f (%.2f, %.2f)'",
                   #", auc.polygon.col='", colorroc, "AA'", ", max.auc.polygon.col='", colorroc, "22'",
                   ", grid=", grid, ", identity=", identity, ", col='", colorroc, "', lty='", ltyroc, "'", #", col='", colorroc, "', lty='", ltyroc, "'",
                   ", print.thres=", threshold, ", print.thres.adj=c(0,0.5), print.thres.cex=0.7, print.thres.pattern='%.2f (%.2f, %.2f)'", 
                   xlab, ylab, main, ")", sep = "")
  doItAndPrint(command)
}

command <- paste("roc.obj$levels[1] # The controls are:", sep = "")
doItAndPrint(command)
command <- paste("roc.obj$levels[2] # The cases are:", sep = "")
doItAndPrint(command)

if (cicompute == "TRUE") {
  cilevel = paste(", conf.level=", cilevel, sep="") 
  cimethod = paste(", method='", cimethod, "'", sep="") 
}
if (ciplot == "TRUE") {
  if (citype == "thresholds") {
    if (cithresholds == "custom") {
      threshold = civalues
    } else {
      threshold = paste("'", cithresholds, "'", sep="") 
    }
    command <- paste("roc.ci.obj <- ci(roc.obj, of='thresholds', thresholds=", threshold, cilevel, cimethod,", boot.n=", cibootn, ", boot.stratified=", cibootstratified,")", sep = "")
    doItAndPrint(command) 
    command <- paste("plot(roc.ci.obj, type='", ciplottype, "', col='#1c61b6AA')", sep = "")
    doItAndPrint(command) 
  } else {
    #check if civalues are probably correct (ex. if the max of them is <=1 and percent was selected then the specification is incorrect it should have been seq(0,100,5)))
    if ((citype == "se") & (citype == "sp")) {
      if ((max(eval(parse(text=as.character(civalues)))) <= 1) & (percent=="TRUE")) {
        Message(message="Maybe you didn't specified well the values, since percent is selected you probably wanted to set seq(0,100,5) (or values between 0-100%) instead of seq(0,1,0.05), since percent is checked", type="warning")    
      }
      if ((max(eval(parse(text=as.character(civalues)))) > 1) & (percent=="FALSE")) {
        Message(message="Maybe you didn't specified well the values, since percent is selected you probably wanted to set seq(0,1,0.05) (or values between 0-1) instead of seq(0,100,5), since percent is not checked", type="warning")    
      }
    }
    if (citype == "se") {
      command <- paste("roc.ci.obj <- ci(roc.obj, of='se', specificities=", civalues, cilevel, cimethod,", boot.n=", cibootn, ", boot.stratified=", cibootstratified,")", sep = "")
      doItAndPrint(command) 
    }
    if (citype == "sp") {
      command <- paste("roc.ci.obj <- ci(roc.obj, of='sp', sensitivities=", civalues, cilevel, cimethod,", boot.n=", cibootn, ", boot.stratified=", cibootstratified,")", sep = "")
      doItAndPrint(command) 
    }
    if (citype == "auc") {
      command <- paste("roc.ci.obj <- ci(roc.obj, of='auc'", cilevel, cimethod,", boot.n=", cibootn, ", boot.stratified=", cibootstratified,")", sep = "")
      doItAndPrint(command) 
      doItAndPrint("roc.ci.obj") 
    }
    command <- paste("plot(roc.ci.obj, type='", ciplottype, "', col='#1c61b6AA')", sep = "")
    doItAndPrint(command) 
  }
}

if (smooth == "TRUE") {
  bandwidth = ""
  density = ""
  if (smoothingmethod == "density") {
    if (smoothinbandwidth == "numeric") {
      bandwidth = paste(", bw=", bandwidthnumeric, "", sep="")    
    } else {
      bandwidth = paste(", bw='", smoothinbandwidth, "'", sep="")      
    }
    #density = paste(", density='", bandwidthwindow, "'", sep="") #!!!!!!!!!!!!!!!!!!       !!!!!!!!!!            !    de modificat! nu stiu cum...
  }
  if (smoothingmethod == "fitdistr") {
      density = paste(", density.cases='", distributioncases, "', density.controls='", distributioncontrols, "'", sep="")      
  }
  command <- paste("lines(smooth(roc.obj, method = '", smoothingmethod, "'", bandwidth, density, "), col='#1c61b6')", sep = "")
  doItAndPrint(command)
}

if (values == "TRUE") {
  doItAndPrint("roc.obj$sensitivities")
  doItAndPrint("roc.obj$specificities")
  doItAndPrint("roc.obj$thresholds")
}

# removing variables
command <- paste("remove(roc.obj)", sep = "")
doItAndPrint(command)
if (ciplot == "TRUE") {
  command <- paste("remove(roc.ci.obj)", sep = "")
  doItAndPrint(command)
}

activateMenus()
tkfocus(CommanderWindow())
  }



OKCancelHelp(helpSubject="plot.roc", reset = "fncvpROC", apply="fncvpROC")

# general tab
tkgrid(getFrame(predictionBox), getFrame(labelBox), sticky = "nw", padx=6, pady=c(6, 6))
tkgrid(directionrbFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(generaldataFrame , dataoptionsFrame, sticky = "nswe", padx=6, pady=6)
tkgrid(generalFrame, sticky = "we")

# smoothing tab 
tkgrid(smoothingmethodrbFrame, sticky = "w", padx=6, pady=c(6, 6))
tkgrid(smoothinbandwidthrbFrame, sticky = "w", padx=6, pady=c(6, 0))
tkgrid(labelRcmdr(smoothingdensityFrame, text = gettext("Numeric bandwidth", domain="R-RcmdrPlugin.ROC")), bandwidthnumericEntry, sticky = "ew", padx=6, pady=c(6, 0))
tkgrid(labelRcmdr(smoothingdensityFrame, text =""), bandwidthnumericScroll, sticky = "ew", padx=6)
 tkgrid(labelRcmdr(smoothingdensityFrame, text = gettext("Adjustment", domain="R-RcmdrPlugin.ROC")), bandwidthadjustmentEntry, sticky = "ew", padx=6, pady=c(6, 0)) #adaugat!!!
tkgrid(bandwidthwindowrbFrame, sticky = "w", padx=6, pady=c(6, 6))
tkgrid(distributioncontrolsrbFrame, sticky = "w", padx=6, pady=c(6, 6))
tkgrid(distributioncasesrbFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(smoothinggeneralFrame, sticky = "w")
tkgrid(smoothingdensityFrame, sticky = "w")
tkgrid(smoothingleftpaneFrame , smoothingdistributionFrame, sticky = "nswe", padx=6, pady=6)
tkgrid(smoothingFrame, sticky = "we")

# ci tab
tkgrid(labelRcmdr(cigeneralFrame, text = gettext("Confidence level", domain="R-RcmdrPlugin.ROC")), cilevelEntry, sticky = "ew", padx=6)
tkgrid(cimethodrbFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(cityperbFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(cithresholdsrbFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(labelRcmdr(cigeneralFrame, text = gettext("Values (Se/Sp/Custom thres.)", domain="R-RcmdrPlugin.ROC")), civaluesEntry, sticky = "ew", padx=6, pady=c(0, 6))
tkgrid(labelRcmdr(cigeneralFrame, text =""), civaluesScroll, sticky = "ew", padx=6, pady=c(0, 6))
tkgrid(cigeneralFrame , cibootstrapFrame, sticky = "nswe", padx=6, pady=6)
tkgrid(ciFrame, sticky = "we")

# auc tab
tkgrid(partialfocusFrame, sticky = "w", padx=6, pady=c(6, 6))
#tkgrid(partialcorrectFrame, sticky = "w")
tkgrid(generalaucFrame , partialaucFrame, sticky = "nswe", padx=6, pady=6)
tkgrid(aucFrame, sticky = "we")

# plot tab
tkgrid(optionsFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(aucpolygonFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(printthresrbFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(labelRcmdr(informationFrame, text = gettext("Custom threshold", domain="R-RcmdrPlugin.ROC")), customthresEntry, sticky = "ew", padx=6)
tkgrid(labelRcmdr(informationFrame, text =""), customthresScroll, sticky = "ew", padx=6)
tkgrid(informationFrame, sticky = "w", padx=6, pady=c(0, 6))
tkgrid(ciplottyperbFrame, sticky = "w", padx=6, pady=c(6, 6))
tkgrid(getFrame(colorrocBox), sticky = "w", padx=6, pady=c(6, 0))
tkgrid(getFrame(ltyrocBox), sticky = "w", padx=6, pady=c(6, 18))

tkgrid(optFrame , parFrame, sticky = "nswe", padx=6, pady=6)
tkgrid(optionsParFrame, sticky = "we")
tkgrid(ttklabel(dataTab, text=""))
tkgrid(ttklabel(dataTab, text=""))
tkgrid(labelRcmdr(top, text = " "), padx=6)
dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tabs=c("dataTab", "smoothingTab", "aucTab", "ciTab", "optionsTab"), 
             tab.names=c("General", "Smoothing", "AUC", "CI", "Plot")) #
}    
# library(pROC)
# data(aSAH)
# roc.obj <- pROC::roc(outcome ~ age, data=aSAH, conf.level=0.95,  auc=TRUE, 
#                      partial.auc=FALSE, partial.auc.focus='specificity', plot=TRUE, add=FALSE, smooth=FALSE, 
#                      print.auc=FALSE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE, identity=TRUE, 
#                      print.thres='no')
# ci.obj <- ci(roc.obj, of="se", specificities=seq(0, 1, 0.05), method="bootstrap", boot.n=500, boot.stratified=F)
# plot(ci.obj, type='shape', col='#5c6126ff')
#=========================================================================================================================================
