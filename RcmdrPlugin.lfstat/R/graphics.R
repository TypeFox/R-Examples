#########################
# Seasonal Barchart     #
#########################

sbplotcalc <- function(){
  command <- paste('sbplot(',ActiveDataSet(),')',sep = "")
  doItAndPrint(command)
}

#########################
# Baseflow Plot         #
#########################

bfplotcalc <- function(){
  initializeDialog(title = gettextRcmdr("Baseflow Plot"))
years <- justDoIt(paste('min(',ActiveDataSet(),'$hyear):max(',ActiveDataSet(),'$hyear)',sep = ""))
  choice <-  c(gettextRcmdr("Whole period"),gettextRcmdr(years))
yearsBox <- variableListBox(top, choice, title=gettextRcmdr("Hydrological year"),initialSelection=which(choice ==getOption("RcmdrPlugin.lfstat")$hydroyearlist)-1)

  onOK <- function(){
        year <- getSelection(yearsBox)
        options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(hydroyearlist = year)))
        if(year == "Whole period"){year <- '"any"'}

        closeDialog()

        command <- paste('bfplot(lfobj = ',ActiveDataSet(),', year = ',year,')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "bfplot")
tkgrid(getFrame(yearsBox), sticky="nw")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=2, columns=1)
}

#########################
# Flow Duration Curve   #
#########################
fdccalc <- function(){
  initializeDialog(title = gettextRcmdr("Flow Duration Curve"))
  optionsFrame <- tkframe(top)
  breakFrame <- tkframe(top)
years <- justDoIt(paste('min(',ActiveDataSet(),'$hyear):max(',ActiveDataSet(),'$hyear)',sep = ""))
yearsBox <- variableListBox(top, c(gettextRcmdr("Whole period"),years), title=gettextRcmdr("Hydrological year (start)"),initialSelection=0)
yearsBox2 <- variableListBox(top, years, title=gettextRcmdr("Hydrological year (end)"),initialSelection=0)

breakseas <- tclVar(gettextRcmdr(getlfopt("fdcseason")))
entrybreakseas <- ttkentry(breakFrame, width="25", textvariable=breakseas)
  
log <- tclVar(getlfopt("fdclog"))
logCheckBox <- tkcheckbutton(optionsFrame, variable = log)
xnorm  <- tclVar(getlfopt("fdcxnorm"))
xnormCheckBox <- tkcheckbutton(optionsFrame, variable = xnorm)
leg <- tclVar(getlfopt("fdclegend"))
legCheckBox <- tkcheckbutton(optionsFrame, variable = leg)
color <- tclVar(getlfopt("fdccolors"))
colCheckBox <- tkcheckbutton(optionsFrame, variable = color)
sep <- tclVar(getlfopt("fdcseparate"))
sepCheckBox <- tkcheckbutton(optionsFrame, variable = sep)

  
  onOK <- function(){
        year1 <- getSelection(yearsBox)
        if(year1 == "Whole period"){year1 <- '"any"'}
        year2 <- getSelection(yearsBox2)
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(fdcyear1 = year1)))
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(fdcyear2 = year2)))
        if(year1 == '"any"'){
          years = '"any"'} else {
            years = paste("c(",year1,",",year2,")",sep = "")}
        closeDialog()
        log <- tclvalue(log) == "1"
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(fdclog = log)))
        xnorm <- tclvalue(xnorm) == "1"
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(fdcxnorm = xnorm))) 
        leg <- tclvalue(leg)=="1"
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(fdclegend = leg))) 
        col <- tclvalue(color) == "1"
        options("RcmdrPlugin.lfstat" =
                      modifyList(getOption("RcmdrPlugin.lfstat"),list(fdccolors = col))) 
        sep <- tclvalue(sep) == "1"
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(fdcseparate = sep)))
       seas <- tclvalue(breakseas)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(fdcseason = seas)))
        command <- paste('fdc(lfobj =',ActiveDataSet(),', year = ',years,',breakdays = c(',seas,'),ylog =',log,',xnorm =',xnorm,',colors=',col,', legend =',leg,',separate =',sep,')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
      }
OKCancelHelp(helpSubject = "fdc")
tkgrid(getFrame(yearsBox),getFrame(yearsBox2), sticky="nw")

tkgrid(labelRcmdr(breakFrame, text=gettextRcmdr("Breakdays (optional): ")),entrybreakseas , sticky="w")
tkgrid(breakFrame, sticky = "w")
#tkgrid(getFrame(yearsBox2),sticky="ne")
tkgrid(optionsFrame,sticky="w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Log transform y-axes:")),
       logCheckBox, sticky = "w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Use normal scale on x-axes:")),
       xnormCheckBox, sticky = "w")
  tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Plot legend:")),
       legCheckBox, sticky = "w")
    tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Use colors:")),
       colCheckBox, sticky = "w")
    tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Get a separate plot for every period:")),
       sepCheckBox, sticky = "w")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=7, columns=2)
}

#############################
# Streamflow deficite plot  #
#############################

streamdefplotcalc <- function(){
  initializeDialog(title = gettextRcmdr("Streamflow deficit plot"))
  optionsFrame <- tkframe(top)
  levelFrame <- tkframe(top)
  breakFrame <- tkframe(top)
  paraFrame <- tkframe(top)

years <- justDoIt(paste('min(',ActiveDataSet(),'$hyear):max(',ActiveDataSet(),'$hyear)',sep = ""))

yearsBox <- variableListBox(top,
                            years,
                            title=gettextRcmdr("Hydrological year"),
                            initialSelection=if(is.null(getOption("RcmdrPlugin.lfstat")$strplotyear)){
                              0}else{
                                which(years ==getOption("RcmdrPlugin.lfstat")$strplotyear)-1})
tlevel <- tclVar(gettextRcmdr(getlfopt("strlevel")))  
entrytlevel <- ttkentry(levelFrame, width="2", textvariable=tlevel)

 #Maybe do this like QQ-Plot!!! (Breakdays) 
 radioButtons(levelFrame,
               "thresbreaks",
               buttons = c("fixed", "monthly", "daily", "seasonal"),
               labels = c("fixed", "monthly", "daily", "seasonal"),
               title = gettextRcmdr("Threshold selection:"),
               initialValue=getlfopt("strthresbreaks"))

 breakseas <- tclVar(gettextRcmdr(getlfopt("season")))
 entrybreakseas <- ttkentry(breakFrame, width="25", textvariable=breakseas)
 
  
onOK <- function(){
         year <- getSelection(yearsBox)
         options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strplotyear = as.numeric(year))))
		closeDialog()

threshold <- tclvalue(tlevel)
         options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strlevel = as.numeric(threshold))))        
breaks <- tclvalue(thresbreaksVariable)
         options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strthresbreaks = breaks)))        
seas <- tclvalue(breakseas)
         options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(season = seas)))
                
command <- paste("streamdefplot(lfobj = ", ActiveDataSet(),",year = ",year, ", threslevel =", threshold, ", thresbreaks = \"",breaks,"\", breakdays =c(",seas,"))",sep = "")

  doItAndPrint(command)
  #Class?! Deftable needed?!              
  tkfocus(CommanderWindow())
              } #END ONok

OKCancelHelp(helpSubject="streamdefplot")
tkgrid(getFrame(yearsBox), sticky="nw")  
tkgrid(labelRcmdr(levelFrame, text=gettextRcmdr("Threshold: Q")), entrytlevel, sticky="w")
tkgrid(levelFrame, sticky="w")
tkgrid(thresbreaksFrame, sticky="w")

tkgrid(labelRcmdr(breakFrame, text=gettextRcmdr("Breakdays (seasonal only): ")),entrybreakseas , sticky="w")
tkgrid(breakFrame, sticky = "w")

tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=5, columns=1)
}

#########################
#Double Mass Curve      #
#########################
dmcurvecalc <- function(){
	initializeDialog(title=gettextRcmdr("Double mass curve"))
        optionsFrame <-  tkframe(top)

	lfobjBox1 <- variableListBox(top, listlfobj(),
		selectmode="single", title=gettextRcmdr("Select low flow object x-Axis: "))
	lfobjBox2 <- variableListBox(top, listlfobj(),
		selectmode="single", title=gettextRcmdr("Select low flow object y-Axis: "))

onOK <-  function(){
  lfobj1 <- getSelection(lfobjBox1)
  lfobj2 <- getSelection(lfobjBox2)
  if (length(lfobj1)+length(lfobj2) < 2 ) {
     errorCondition(recall=Recode, message=gettextRcmdr("You must select an object."))
			return()
		}
  command <- paste('dmcurve(x =',lfobj1,', y =', lfobj2,')')
  doItAndPrint(command)  
} #End on OK


OKCancelHelp(helpSubject="dmcurve")
	tkgrid(getFrame(lfobjBox1),getFrame(lfobjBox2), sticky="nw")
        tkgrid(buttonsFrame, sticky = "w")
        dialogSuffix(rows=2, columns=2)
      }

#########################
# Hydrograph            #
#########################

hydrocalc<-function(){
initializeDialog(title=gettextRcmdr("Hydrograph"))
levelFrame <- tkframe(top)
optionsFrame <- tkframe(top)

radioButtons(top,
               "period",
               buttons = c("A", "B"),
               labels = c("Use whole period", "Enter dates"),
               title = gettextRcmdr("Period:"),
               initialValue=getlfopt("hchoice"))
  
sdate <- tclVar(gettextRcmdr(getlfopt("hsdate")))  
entrysdate <- ttkentry(levelFrame, width="12", textvariable=sdate)
edate <- tclVar(gettextRcmdr(getlfopt("hedate")))  
entryedate <- ttkentry(levelFrame, width="12", textvariable=edate)
min <- tclVar(getlfopt("hminpoints"))
mincheckbox <- tkcheckbutton(optionsFrame, variable = min)
onOK <- function(){
  closeDialog()
  per <- tclvalue(periodVariable)
            options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(hchoice = per)))
  start <- "NULL"
  end <- "NULL"
  if(per == "B"){
    start <- tclvalue(sdate)
    end <- tclvalue(edate)
            options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(hsdate = start)))
            options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(hedate = end)))
    start = paste('"',start,'"',sep = "")
    end = paste('"',end,'"',sep = "")
  }
 minima <- tclvalue(min) == "1"
              options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(hminpoints =minima)))

  command <- paste('hydrograph(lfobj =',ActiveDataSet(),',startdate = ',start,',enddate = ',end,',amin =',minima,')',sep = "")
  doItAndPrint(command)
  tkfocus(CommanderWindow())



}

OKCancelHelp(helpSubject="hydrograph")
tkgrid(periodFrame, sticky = "w")
tkgrid(labelRcmdr(levelFrame, text=gettextRcmdr("Startdate: ")),entrysdate , sticky="w")
tkgrid(labelRcmdr(levelFrame, text=gettextRcmdr("Enddate: ")),entryedate , sticky="w")
tkgrid(levelFrame,sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Mark annual minima:")),
       mincheckbox, sticky = "w")
tkgrid(optionsFrame,sticky = "w")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=5, columns=2)
}

