getlfopt <- function(x){
  justDoIt(paste('getOption("RcmdrPlugin.lfstat")$',x,sep = ""))}

resetlfoptions <- function(){
para <- list(
             datasheet = "GRDC",
             hyear = 1,
             baseflow = TRUE,
             hydroyearlist = "Whole period",
             monthsep = FALSE,
             yearsep = FALSE,
             seasonsep = FALSE,
             MAM = 7,
             recessionmethod = "MRC",
             rainpeaklevel = 0.95,
             segLength = 7,
             threslevelrec = 70,
             flowplot = TRUE,
             tmean = .1,
             strmethod = "none",
             strlevel = 70,
             strthresbreaks = "fixed",
             strtable = "all",
             season = '"01/06","01/10"',
             strMadays = 7,
             strtmin = 5,
             strIClevel = .1,
             strmindur = 0,
             strminvol = 0,
             strplotyear = NULL,
             extyears = 100,
             extn = 5,
             extdist = "Weibull",
             fdclog = TRUE,
             fdcxnorm = FALSE,
             fdclegend = TRUE,
             fdcseparate = FALSE,
             fdccolors = TRUE,
             fdcseason = NULL,
             fdcyear1 = "Whole period",
             fdcyear2 = NULL,
             seasq = 95,
             seasbreak = '"01/06","01/10"',
             dmyear = "Whole period",
             hsdate = "01/01/1960",
             hedate = "31/12/2010",
             hchoice = "A",
             hminpoints = TRUE
             )
 
  options("RcmdrPlugin.lfstat" = para)
}

savelfopt <- function(){
 	file <- tclvalue(tkgetSaveFile(filetypes=
				gettextRcmdr('{"R Data Files" {".RData"}}'),
			defaultextension=".RData", initialfile=paste("lfoptions.RData", sep="")))
	if (file == ""){ return()}
	command <- paste('lfoptions <- list(RC = getOption("RcmdrPlugin.lfstat"),R = getOption("lfstat"));save(lfoptions, file="', file, '");rm(lfoptions)', sep="")
	doItAndPrint(command)
}

loadlfopt <- function(){
     	file <- tclvalue(tkgetOpenFile(filetypes=
				gettextRcmdr('{"R Data Files" {".RData"}}')))
	if (file == "") return()
	command <- paste('load("', file,'")', sep="")
	lf <- justDoIt(command)
        if("lfoptions" == lf)
        {command2 <- paste('options("RcmdrPlugin.lfstat" =', lf,'$RC)',sep = "")
        justDoIt(command2)
         command3 <- paste('options("lfstat" =', lf,'$R)',sep = "")
        justDoIt(command3)
        justDoIt("rm(lfoptions)")
        Message(message=gettextRcmdr("Options loaded"),type = "note")}
        else{stop("This file does not contain LF-options")}
	tkfocus(CommanderWindow())
}

setunitcalc <- function(){
 initializeDialog(title = gettextRcmdr("Set Unit in Plots"))
 optionsFrame <- tkframe(top)
 radioButtons(top,
               "unit",
               buttons = c(LETTERS[1:12],"Other"),
               labels = gettextRcmdr(c(LETTERS[1:12],"Custom")),
               title = gettextRcmdr("Select Unit:"))

self <- tclVar(gettextRcmdr(getOption("lfstat")$unit))
entryunit <- ttkentry(optionsFrame, width="20", textvariable=self)
 
 expre <- function(unit){
        switch(which(unit == LETTERS),
               "m^3/s",
               "m^{3}*s^{-1}",
               "scriptscriptstyle(frac(m^3,s))",
               "l/s",
               "l*s^{-1}",
               "scriptscriptstyle(frac(l,s))",
               "m^3/s/km^2",
               "m^3*s^{-1}*km^{-2}",
               "scriptscriptstyle(frac(m^3,s%.%km^2))",
               "l/s/km^2",
               "l*s^{-1}*km^{-2}",
               "scriptscriptstyle(frac(l,s%.%km^2))"
       )
      }
 onOK <- function(){
 closeDialog()
 unit <- tclvalue(unitVariable)
 if(unit == "Other"){
   expout <- tclvalue(self)
 } else{
 expout <- expre(unit=unit)
 }
 command <- paste('setlfunit("',expout,'")',sep = "")
 doItAndPrint(command)
 tkfocus(CommanderWindow())
       }

 Pressedshowlabels <- function(){
   plot(1:3, c(1,2,4), type = "n",xaxt = "n",yaxt = "n",xlab = "",ylab = "", xlim = c(0,4), ylim = c(0,5))
   lab <- NULL
     for(ii in 1:12){
       lab[ii] <- parse(text =expre(LETTERS[ii]))}
     oth <- trim.blanks(tclvalue(self))
     lab[13] <- if(oth == ""){""}else{parse(text = oth)}
    
   text(c(rep(1:3,4),2), c(rep(4,3),rep(3,3),rep(2,3),rep(1,3),.3),labels =lab,pos = 4)
   text(c(rep(1:3,4),2), c(rep(4,3),rep(3,3),rep(2,3),rep(1,3),.3),labels =c(LETTERS[1:12],"Custom"),pos = 2,cex = 1,col = "blue")
   
 }
 unitplot <- buttonRcmdr(top,text=gettextRcmdr("Show labels"),command=Pressedshowlabels)
 OKCancelHelp(helpSubject = "setlfunit")
 tkgrid(unitFrame,sticky = "w")
 tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Enter unit:")), entryunit, sticky="w")
 tkgrid(optionsFrame,sticky = "w")
 tkgrid(unitplot,sticky = "w")
 tkgrid(buttonsFrame)
 dialogSuffix(rows=4, columns=2)
}
  
