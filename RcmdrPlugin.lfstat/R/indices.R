###########################
# MEAN FLOW               #
###########################
meanflowcalc <- function(){
initializeDialog(title = gettextRcmdr("Mean Flow"))
optionsFrame <- tkframe(top)
years <- justDoIt(paste('min(',ActiveDataSet(),'$hyear):max(',ActiveDataSet(),'$hyear)',sep = ""))
choice <-  c(gettextRcmdr("Whole period"),gettextRcmdr(years))
init <-NULL
for(ii in seq_along(getlfopt("hydroyearlist"))){
  init[ii] <- which(getlfopt("hydroyearlist")[ii] == choice)-1
}


yearsBox <- variableListBox(top,
                            choice,
                            title=gettextRcmdr("Hydrological year"),
                            initialSelection = init,
                            selectmode = "multiple")
#monthly <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$monthsep+0)))
#monthlyCheckBox <- tkcheckbutton(optionsFrame, variable = monthly)
yearly <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$yearsep+0)))
yearlyCheckBox <- tkcheckbutton(optionsFrame, variable = yearly)
season <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$seasonsep+0)))
seasonCheckBox <- tkcheckbutton(optionsFrame, variable = season)
breakseas <- tclVar(gettextRcmdr(getlfopt("season")))
entrybreakseas <- ttkentry(optionsFrame, width="25", textvariable=breakseas)


onOK <- function(){
        year <- getSelection(yearsBox)
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(hydroyearlist = year)))
        if("Whole period" %in% year){yearvec <- '"any"'}else{
          yearvec <- NULL
          for(ii in seq_along(year)){
            yearvec <- paste(yearvec,year[ii],sep =if(ii == 1){""}else{","})
          }
        yearvec <- paste("c(",yearvec,")",sep="")
        }
        closeDialog()
        #month <- tclvalue(monthly) == "1"
        #  options("RcmdrPlugin.lfstat" =
        #                modifyList(getOption("RcmdrPlugin.lfstat"),list(monthsep = month)))
        yearl <- tclvalue(yearly) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(yearsep = yearl)))
        seas <- tclvalue(season) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(seasonsep = seas)))
        bday <- tclvalue(breakseas)
          options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(season = bday)))
        command <- paste('meanflow(lfobj =',ActiveDataSet(),
                         ', year = ',yearvec,
                         #',monthly =',month,
                         ',yearly =', yearl,
                         if(seas){paste(',breakdays = c(', bday,')',sep="")}else{" "},
                         ')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "meanflow")

tkgrid(getFrame(yearsBox), sticky="nw")
tkgrid(optionsFrame, sticky="w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each year separately:")),
       yearlyCheckBox, sticky = "w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each season separately:")),
       seasonCheckBox, sticky = "w")

tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter Breakdays:")),entrybreakseas , sticky="w")


tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)                   
}
###########################
# Q95                     #
###########################
Q95calc <- function(){
initializeDialog(title = gettextRcmdr("Q95"))

optionsFrame <- tkframe(top)
years <- justDoIt(paste('min(',ActiveDataSet(),'$hyear):max(',ActiveDataSet(),'$hyear)',sep = ""))
choice <-  c(gettextRcmdr("Whole period"),gettextRcmdr(years))
init <-NULL
for(ii in seq_along(getlfopt("hydroyearlist"))){
  init[ii] <- which(getlfopt("hydroyearlist")[ii] == choice)-1
}


yearsBox <- variableListBox(top,
                            choice,
                            title=gettextRcmdr("Hydrological year"),
                            initialSelection = init,
                            selectmode = "multiple")
#monthly <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$monthsep+0)))
#monthlyCheckBox <- tkcheckbutton(optionsFrame, variable = monthly)
yearly <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$yearsep+0)))
yearlyCheckBox <- tkcheckbutton(optionsFrame, variable = yearly)
season <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$seasonsep+0)))
seasonCheckBox <- tkcheckbutton(optionsFrame, variable = season)
breakseas <- tclVar(gettextRcmdr(getlfopt("season")))
entrybreakseas <- ttkentry(optionsFrame, width="25", textvariable=breakseas)

onOK <- function(){
        year <- getSelection(yearsBox)
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(hydroyearlist = year)))
        if("Whole period" %in% year){yearvec <- '"any"'}else{
          yearvec <- NULL
          for(ii in seq_along(year)){
            yearvec <- paste(yearvec,year[ii],sep =if(ii == 1){""}else{","})
          }
        yearvec <- paste("c(",yearvec,")",sep="")
        }
        closeDialog()
        #month <- tclvalue(monthly) == "1"
        #  options("RcmdrPlugin.lfstat" =
        #                modifyList(getOption("RcmdrPlugin.lfstat"),list(monthsep = month)))
        yearl <- tclvalue(yearly) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(yearsep = yearl)))
        seas <- tclvalue(season) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(seasonsep = seas)))
        bday <- tclvalue(breakseas)
          options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(season = bday)))
        command <- paste('Q95(lfobj =',ActiveDataSet(),
                         ', year = ',yearvec,
                         #',monthly =',month,
                         ',yearly =', yearl,
                         if(seas){paste(',breakdays = c(', bday,')',sep="")}else{" "},
                         ')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "Q95")

tkgrid(getFrame(yearsBox), sticky="nw")
tkgrid(optionsFrame, sticky="w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each year separately:")),
       yearlyCheckBox, sticky = "w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each season separately:")),
       seasonCheckBox, sticky = "w")

tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter Breakdays:")),entrybreakseas , sticky="w")


tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)                   
}



###########################
# MAM                     #
###########################

MAMcalc <- function(){
initializeDialog(title = gettextRcmdr("Mean Annual Minima (n-day)"))

#   sliderValue <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$MAM)) #Som#ething better then a slider????
#    slider <- tkscale(top, from=1, to=90, showvalue=TRUE, variable=sliderValue,
#        resolution=1, orient="horizontal")

optionsFrame <- tkframe(top)

ndays <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$MAM))
entryn <- ttkentry(optionsFrame, width="2", textvariable=ndays)


years <- justDoIt(paste('min(',ActiveDataSet(),'$hyear):max(',ActiveDataSet(),'$hyear)',sep = ""))
choice <-  c(gettextRcmdr("Whole period"),gettextRcmdr(years))
init <-NULL
for(ii in seq_along(getlfopt("hydroyearlist"))){
  init[ii] <- which(getlfopt("hydroyearlist")[ii] == choice)-1
}


yearsBox <- variableListBox(top,
                            choice,
                            title=gettextRcmdr("Hydrological year"),
                            initialSelection = init,
                            selectmode = "multiple")
#monthly <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$monthsep+0)))
#monthlyCheckBox <- tkcheckbutton(optionsFrame, variable = monthly)
yearly <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$yearsep+0)))
yearlyCheckBox <- tkcheckbutton(optionsFrame, variable = yearly)
season <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$seasonsep+0)))
seasonCheckBox <- tkcheckbutton(optionsFrame, variable = season)
breakseas <- tclVar(gettextRcmdr(getlfopt("season")))
entrybreakseas <- ttkentry(optionsFrame, width="25", textvariable=breakseas)


onOK <- function(){
        year <- getSelection(yearsBox)
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(hydroyearlist = year)))
        if("Whole period" %in% year){yearvec <- '"any"'}else{
          yearvec <- NULL
          for(ii in seq_along(year)){
            yearvec <- paste(yearvec,year[ii],sep =if(ii == 1){""}else{","})
          }
        yearvec <- paste("c(",yearvec,")",sep="")
        }
        closeDialog()
 nday <- tclvalue(ndays)
          options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(MAM = nday)))
        
        yearl <- tclvalue(yearly) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(yearsep = yearl)))
        seas <- tclvalue(season) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(seasonsep = seas)))
        bday <- tclvalue(breakseas)
          options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(season = bday)))
        command <- paste('MAM(lfobj =',ActiveDataSet(),
                         ',n = ',nday,
                         ', year = ',yearvec,
                         #',monthly =',month,
                         ',yearly =', yearl,
                         if(seas){paste(',breakdays = c(', bday,')',sep="")}else{" "},
                         ')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "MAM")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("n-days:")),entryn , sticky="w")
tkgrid(getFrame(yearsBox), sticky="nw")
tkgrid(optionsFrame, sticky="w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each year separately:")),
       yearlyCheckBox, sticky = "w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each season separately:")),
       seasonCheckBox, sticky = "w")

tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter Breakdays:")),entrybreakseas , sticky="w")


tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)                   
}


###########################
# BFI                     #
###########################
BFIcalc <- function(){
initializeDialog(title = gettextRcmdr("Baseflow Index"))
optionsFrame <- tkframe(top)
years <- justDoIt(paste('min(',ActiveDataSet(),'$hyear):max(',ActiveDataSet(),'$hyear)',sep = ""))
choice <-  c(gettextRcmdr("Whole period"),gettextRcmdr(years))

init <-NULL
for(ii in seq_along(getlfopt("hydroyearlist"))){
  init[ii] <- which(getlfopt("hydroyearlist")[ii] == choice)-1
}
yearsBox <- variableListBox(top,
                            choice,
                            title=gettextRcmdr("Hydrological year"),
                            initialSelection = init,
                            selectmode = "multiple")
yearly <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$yearsep+0)))
yearlyCheckBox <- tkcheckbutton(optionsFrame, variable = yearly)
season <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$seasonsep+0)))
seasonCheckBox <- tkcheckbutton(optionsFrame, variable = season)
breakseas <- tclVar(gettextRcmdr(getlfopt("season")))
entrybreakseas <- ttkentry(optionsFrame, width="25", textvariable=breakseas)

onOK <- function(){
        year <- getSelection(yearsBox)
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(hydroyearlist = year)))
        if("Whole period" %in% year){yearvec <- '"any"'}else{
          yearvec <- NULL
          for(ii in seq_along(year)){
            yearvec <- paste(yearvec,year[ii],sep =if(ii == 1){""}else{","})
          }
        yearvec <- paste("c(",yearvec,")",sep="")
        }
        closeDialog()
        #month <- tclvalue(monthly) == "1"
        #  options("RcmdrPlugin.lfstat" =
        #                modifyList(getOption("RcmdrPlugin.lfstat"),list(monthsep = month)))
        yearl <- tclvalue(yearly) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(yearsep = yearl)))
        seas <- tclvalue(season) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(seasonsep = seas)))
        bday <- tclvalue(breakseas)
          options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(season = bday)))
        command <- paste('BFI(lfobj =',ActiveDataSet(),
                         ', year = ',yearvec,
                         #',monthly =',month,
                         ',yearly =', yearl,
                         if(seas){paste(',breakdays = c(', bday,')',sep="")}else{" "},
                         ')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "BFI")

tkgrid(getFrame(yearsBox), sticky="nw")
tkgrid(optionsFrame, sticky="w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each year separately:")),
       yearlyCheckBox, sticky = "w")
tkgrid(labelRcmdr(optionsFrame,
                  text = gettextRcmdr("Calculate each season separately:")),
       seasonCheckBox, sticky = "w")

tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter Breakdays:")),entrybreakseas , sticky="w")


tkgrid(buttonsFrame, sticky="w")
	dialogSuffix(rows=4, columns=1)                   
}

###########################
# RECESSION               #
###########################
recessionanalysis <- function(){
  initializeDialog(title=gettextRcmdr("Recession Analysis"))
  optionsFrame <- tkframe(top)
  seasonsFrame <- tkframe(top)
  options2Frame <- tkframe(top)
  
 radioButtons(optionsFrame,
             "method",
            buttons=c("MRC", "IRS"), 
             labels=gettextRcmdr(c("MRC", "IRS")),
             title=gettextRcmdr("Method"),
              initialValue = getOption("RcmdrPlugin.lfstat")$recessionmethod)
  seglen <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$segLength))
  entryseglength <- ttkentry(optionsFrame, width="4", textvariable=seglen)

  plotyear <- tclVar(gettextRcmdr("0"))
  entryplotyear <- ttkentry(optionsFrame, width="4", textvariable = plotyear)

  thrlevel <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$threslevelrec))
  entrythrlevel <- ttkentry(optionsFrame, width = "4", textvariable = thrlevel)
  rainlevel <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$rainpeaklevel))
  entryrain <- ttkentry(optionsFrame, width = "4", textvariable = rainlevel)
  flowplot <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$flowplot + 0))
  plotCheckBox <- tkcheckbutton(options2Frame, variable = flowplot)
  tmean <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$tmean))
  entrytmean <- ttkentry(options2Frame, width= "4", textvariable = tmean)

 radioButtons(seasonsFrame,
               "thresbreaks",
               buttons = c("fixed", "monthly", "seasonal"),
               labels = c("fixed", "monthly", "seasonal"),
               title = gettextRcmdr("Threshold selection:"),
              initialValue = getlfopt("strthresbreaks"))
  
season <- tclVar(gettextRcmdr(paste(getOption("RcmdrPlugin.lfstat")$seasonsep+0)))
seasonCheckBox <- tkcheckbutton(seasonsFrame, variable = season)
breakseas <- tclVar(gettextRcmdr(getlfopt("season")))
entrybreakseas <- ttkentry(seasonsFrame, width="25", textvariable=breakseas)

  
Pressedseglenplot <- function()
{   method <- tclvalue(methodVariable)
    threshold <- as.numeric(tclvalue(thrlevel))
    breaks <- tclvalue(thresbreaksVariable)
    bday <- tclvalue(breakseas)
    rain <- tclvalue(rainlevel)
        command <- paste('seglenplot(',activeDataSet(),
                         ', threslevel = ',threshold,
                         ', thresbreaks = "', breaks,'"',
                         ', thresbreakdays = c(',bday,')',
                         ', rainpeaklevel = ',rain,')',sep = "")
   doItAndPrint(command)
}

Pressedrainplot <- function()
{rain <- as.numeric(tclvalue(rainlevel))
 pyear <- as.numeric(tclvalue(plotyear))
    threshold <- as.numeric(tclvalue(thrlevel))
    breaks <- tclvalue(thresbreaksVariable)
    seg <- tclvalue(seglen)
    bday <- tclvalue(breakseas)
    if(pyear == 0){
    pyear <- sample(eval(parse(text = paste0(activeDataSet(),'$hyear'))),1)}
    command <- paste0('recessionplot(',activeDataSet(),
                      ', peaklevel = ',rain,
                      ', startdate = ',pyear,
                      ', enddate =', pyear,
                      ', seglength =', seg,
                      ', threshold =', threshold,
                      ', thresbreaks = "', breaks,'"',
                      ', thresbreakdays = c(', bday,'))')
     doItAndPrint(command)
  }

onOK <- function(){
 breaks <- tclvalue(thresbreaksVariable)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strthresbreaks = breaks)))
 
 seas <- tclvalue(season) == "1"
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(seasonsep = seas)))
 
 bday <- tclvalue(breakseas)
          options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(season = bday)))
  
closeDialog()
method <- tclvalue(methodVariable)
rain <- tclvalue(rainlevel)
seg <- tclvalue(seglen)
thr <- tclvalue(thrlevel)
plot <- tclvalue(flowplot)=="1"
trim <- tclvalue(tmean)

options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(threslevelrec = thr)))
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(rainpeaklevel = rain)))

options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(recessionmethod = method)))
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(segLength = seg)))
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(flowplot = plot)))
options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(tmean = trim)))
sdays <- NULL
if(seas){sdays <- bday}
 
command <- paste('recession(lfobj = ', activeDataSet(),
                 ',method = "', method,
                 '", seglength = ',seg,
                 ',threshold = ',thr,
                 ',peaklevel = ',rain,
                 ',seasonbreakdays = c(',sdays,')',
                 ',thresbreaks = "',breaks,'"',
                 ',thresbreakdays = c(', bday,')',
                 ',plotMRC = ',plot,
                 ',trimIRS = ',trim,')',sep = "")
doItAndPrint(command)
tkfocus(CommanderWindow())
} # END onOK

seglenplot <- buttonRcmdr(top,text=gettextRcmdr("Barchart of recession length"),command=Pressedseglenplot)
rainplot <- buttonRcmdr(top,text=gettextRcmdr('Recession plot ("0" for random year)'),command=Pressedrainplot)
OKCancelHelp(helpSubject="recession")

tkgrid(optionsFrame,sticky="w")
tkgrid(methodFrame,sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Peaklevel:")), entryrain,rainplot, entryplotyear, sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Threshold level:")), entrythrlevel, sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Segment length:")), entryseglength,seglenplot, sticky = "w")


tkgrid(seasonsFrame,sticky = "w")
tkgrid(thresbreaksFrame, sticky = "w")
tkgrid(labelRcmdr(seasonsFrame, text=gettextRcmdr("Enter Breakdays:")),entrybreakseas , sticky="w")
tkgrid(labelRcmdr(seasonsFrame,
                  text = gettextRcmdr("Calculate each season separately:")),
       seasonCheckBox, sticky = "w")
tkgrid(options2Frame, sticky = "w")
tkgrid(labelRcmdr(options2Frame, text = gettextRcmdr("Plot pairs (MRC only):")),plotCheckBox, sticky = "w")
tkgrid(labelRcmdr(options2Frame, text = gettextRcmdr("Trim level (IRS only):")), entrytmean, sticky = "w")
tkgrid(buttonsFrame,sticky = "w")



dialogSuffix(rows=6, columns=2)                   
}
#########################
# Seasonality Index     #
#########################


seasindexcalc <- function(){
initializeDialog(title = gettextRcmdr("Seasonality Index"))
optionsFrame <- tkframe(top)
qlevel <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$seasq))
entryqlevel <- ttkentry(optionsFrame, width = "3", textvariable = qlevel)


onOK <- function(){
        closeDialog()
        
        q <- as.numeric(tclvalue(qlevel))
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(seasq = q)))
        
        command <- paste('seasindex(lfobj =',ActiveDataSet(),',Q = ',q,')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "seasindex")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Q:")), entryqlevel, sticky = "w")
tkgrid(optionsFrame, sticky="nw")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=2, columns=1)                                     
}
#########################
# Seasonal Ratio        #
#########################


seasratiocalc <- function(){

initializeDialog(title = gettextRcmdr("Seasonality Ratio"))
optionsFrame <- tkframe(top)
qlevel <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$seasq))
entryqlevel <- ttkentry(optionsFrame, width = "3", textvariable = qlevel)
breakdays <- tclVar(gettextRcmdr(getOption("RcmdrPlugin.lfstat")$seasbreak))
entrybd <- ttkentry(optionsFrame, width = "25", textvariable = breakdays)

onOK <- function(){
        closeDialog()
        
        q <- as.numeric(tclvalue(qlevel))
        bd <- tclvalue(breakdays)
        options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(seasq = q)))
        
        command <- paste('seasratio(lfobj =',ActiveDataSet(),',breakdays = c(',bd,'),Q = ',q,')',sep = "")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "seasratio")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Q:")), entryqlevel, sticky = "w")
tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Breakday(s):")), entrybd, sticky = "w")
tkgrid(optionsFrame, sticky="nw")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=2, columns=1)
}
