#########################
# Streamflow deficite   #
#########################

streamdefcalc <- function(){
  initializeDialog(title = gettextRcmdr("Streamflow deficite"))
  optionsFrame <- tkframe(top)
  levelFrame <- tkframe(top)
  breakFrame <- tkframe(top)
  paraFrame <- tkframe(top)
  tablename <- tclVar(gettextRcmdr(paste("StrdefTableof",ActiveDataSet(),sep = "")))
  entrytablename <- ttkentry(optionsFrame, width="25", textvariable=tablename)

 radioButtons(optionsFrame,
             "pooling",
             buttons=c("none", "MA","IT","IC"), 
             labels=gettextRcmdr(c("none", "MA", "IT", "IC" )),
             title=gettextRcmdr("Pooling method:"),
             initialValue=getlfopt("strmethod"))
  tlevel <- tclVar(gettextRcmdr(getlfopt("strlevel")))  
entrytlevel <- ttkentry(levelFrame, width="2", textvariable=tlevel)

 #Maybe do this like QQ-Plot!!! (Breakdays) 
 radioButtons(levelFrame,
               "thresbreaks",
               buttons = c("fixed", "monthly", "daily", "seasonal"),
               labels = c("fixed", "monthly", "daily", "seasonal"),
               title = gettextRcmdr("Threshold selection:"),
              initialValue = getlfopt("strthresbreaks"))

 breakseas <- tclVar(gettextRcmdr(getlfopt("season")))
 entrybreakseas <- ttkentry(breakFrame, width="25", textvariable=breakseas)
 
 Madays <- tclVar(gettextRcmdr(getlfopt("strMadays")))
 entrymadays <- ttkentry(paraFrame, width = "3", textvariable = Madays)

 tmin <- tclVar(gettextRcmdr(getlfopt("strtmin")))
 entrytmin <- ttkentry(paraFrame, width = "3", textvariable = tmin)

 IClevel <- tclVar(gettextRcmdr(getlfopt("strIClevel")))
 entryIClevel <- ttkentry(paraFrame, width = "3", textvariable = IClevel)

 mindur <- tclVar(gettextRcmdr(getlfopt("strmindur")))
  entrymindur <- ttkentry(paraFrame, width = "3", textvariable = mindur)
 minvol <- tclVar(gettextRcmdr(getlfopt("strminvol")))
  entryminvol <- ttkentry(paraFrame, width = "3", textvariable = minvol)

radioButtons(top,
               "table",
               buttons = c("all","volmax", "durmax"),
               labels = gettextRcmdr(c("All deficites", "Annual extremes (vol)", "Annual extremes (dur)")),
               title = gettextRcmdr("Streamflow table contains:"),
               initialValue = getlfopt("strtable"))

  
onOK <- function(){
          pool <- tclvalue(poolingVariable)
           options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(strmethod = pool)))

          threshold <- tclvalue(tlevel)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strlevel = as.numeric(threshold))))

          breaks <- tclvalue(thresbreaksVariable)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strthresbreaks = breaks)))

          seas <- tclvalue(breakseas)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(season = seas)))

          MAdays <- tclvalue(Madays)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strMadays = as.numeric(MAdays))))

          Tmin <- tclvalue(tmin)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strtmin = as.numeric(Tmin))))

          icLevel <- tclvalue(IClevel)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strIClevel = as.numeric(icLevel))))

          minDur <- tclvalue(mindur)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strmindur = as.numeric(minDur))))

          minVol <- tclvalue(minvol)
           options("RcmdrPlugin.lfstat" =
              modifyList(getOption("RcmdrPlugin.lfstat"),list(strminvol = as.numeric(minVol))))

          table <- tclvalue(tableVariable)
          options("RcmdrPlugin.lfstat" =
                  modifyList(getOption("RcmdrPlugin.lfstat"),list(strtable = table)))
          
		closeDialog()
                nameValue <- trim.blanks(tclvalue(tablename))
                if (nameValue == ""){
			errorCondition(recall=readDataSet,
				message=gettextRcmdr("You must enter a name for the data set."))
			return()
		}
		if (!is.valid.name(nameValue)){
			errorCondition(recall=readDataSet,
				message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
			return()
		}
                if (is.element(nameValue, listDataSets())) {
			if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
				readDataSet()
				return()
			}
		}
        

                
command <- paste("streamdef(lfobj = ", ActiveDataSet(), ", pooling = \"",pool,"\", threslevel =",
                 threshold, ", thresbreaks = \"",breaks,"\", breakdays = c(",seas,"), MAdays = ",
                 MAdays, ", tmin = ",Tmin,", IClevel = ", icLevel,", mindur = ", minDur, ",minvol = ",minVol,",table =\"",table,"\")",sep = "")
  doItAndPrint(paste(nameValue, " <- ", command, sep=""))
  #result <- data.frame(justDoIt(command))
  #Class?! Deftable needed?!              
  #gassign(nameValue, result)
  tkfocus(CommanderWindow())
              } #END ONok

OKCancelHelp(helpSubject="streamdef")
tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter name for streamflow table:")), entrytablename, sticky="w")
tkgrid(optionsFrame, sticky="w")
tkgrid(poolingFrame, sticky="w")
  
tkgrid(labelRcmdr(levelFrame, text=gettextRcmdr("Threshold: Q")), entrytlevel, sticky="w")
tkgrid(levelFrame, sticky="w")
tkgrid(thresbreaksFrame, sticky="w")

tkgrid(labelRcmdr(breakFrame, text=gettextRcmdr("Breakdays (seasonal only): ")),entrybreakseas , sticky="w")
tkgrid(breakFrame, sticky = "w")

tkgrid(labelRcmdr(paraFrame, text= gettextRcmdr("MAdays (MA only): ")), entrymadays,sticky = "w")
tkgrid(labelRcmdr(paraFrame, text = gettextRcmdr("tmin (IT, IC only): ")), entrytmin, sticky = "w")
tkgrid(labelRcmdr(paraFrame, text = gettextRcmdr("IC-Level (IC only): ")), entryIClevel, sticky = "w")
tkgrid(labelRcmdr(paraFrame, text= gettextRcmdr("Minimal duration (IT, IC only) ")), entrymindur,sticky = "w")
tkgrid(labelRcmdr(paraFrame, text= gettextRcmdr("Minimal volume (IT, IC only): ")), entryminvol,sticky = "w")
  

tkgrid(paraFrame, sticky = "w")
tkgrid(tableFrame, sticky = "w")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=11, columns=2)
}
