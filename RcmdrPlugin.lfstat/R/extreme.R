#########################
# T-YEARS               #
#########################

tyearscalc <- function(){
initializeDialog(title=gettextRcmdr("T year event"))
optionsFrame <-  tkframe(top)
textFrame <- tkframe(top)

MAdays <- tclVar(gettextRcmdr(getlfopt("extn")))
MAentryframe <- ttkentry(textFrame, width="2", textvariable=MAdays)
tyears <- tclVar(gettextRcmdr(getlfopt("extyears")))
tyearsentryframe <- ttkentry(textFrame, width="4", textvariable=tyears)
choice <- c("Weibull", "GEV", "Lognormal","Gumbel","Pearson Type 3")
init <-NULL

for(ii in seq_along(getlfopt("extdist"))){
  init[ii] <- which(getlfopt("extdist")[ii] == choice)-1
}

distBox <- variableListBox(top, choice,
		selectmode="multiple", title=gettextRcmdr("Select distributions: "),
                           initialSelection = init,
                           listHeight = 5)

onOK <- function(){
  dist <- getSelection(distBox)
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(extdist = dist)))
  closeDialog()
  distname <- NULL
  names <- c("wei","gev","ln3","gum","pe3")
  for(ii in seq_along(dist)){
    distname[ii] <-names[which(dist[ii] == choice)]
    }
  distname2 <- NULL
  for (ii in seq_along(distname)){
    distname2 <- paste(distname2,distname[ii],sep =if(ii == 1){"\""}else{"\",\""})
  }

  event <- tclvalue(tyears)
   options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(extyears = event)))
  n <- tclvalue(MAdays)
   options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(extn = n)))
  command <- paste('tyears(lfobj = ',ActiveDataSet(),', event = ', event,', n = ',n,',dist = c(',distname2,'"))',sep="")
  doItAndPrint(command)
  tkfocus(CommanderWindow())
}#end onOK

OKCancelHelp(helpSubject="tyears")
	tkgrid(textFrame,sticky = "w")
        tkgrid(labelRcmdr(textFrame, text = gettextRcmdr("Annual n-day minima, n:")),MAentryframe,
               sticky = "w")
        tkgrid(labelRcmdr(textFrame, text = gettextRcmdr("Return period T (years):")),tyearsentryframe,
               sticky = "w")
        tkgrid(getFrame(distBox), sticky="nw")
       
        tkgrid(buttonsFrame, sticky = "w")
dialogSuffix(rows=4, columns=2)
}


listlfobj <- function(envir=.GlobalEnv, ...) {
	objects <- ls(envir=envir, ...)
	if (length(objects) == 0) NULL
	else objects[sapply(objects,
		function(.x) "lfobj" == (class(get(.x, envir=envir))[1]))]
}

#########################
#Regional frequency     #
#########################


rfap <- function(){
	initializeDialog(title=gettextRcmdr("Regional frequency analysis - Plot"))
        optionsFrame <-  tkframe(top)
	lfobjBox <- variableListBox(top, listlfobj(),
		selectmode="multiple", title=gettextRcmdr("Select low flow objects: "))
        Madays <- tclVar(gettextRcmdr("1"))
        MAentryframe <- ttkentry(optionsFrame, width = "2", textvariable =Madays)
                  
onOK <-  function(){
  lfobjs <- getSelection(lfobjBox)
  if (length(lfobjs) == 0) {
     errorCondition(recall=Recode, message=gettextRcmdr("You must select an object."))
			return()
		}
  a<- NULL
  for(ii in lfobjs){
    if(is.null(a)){a <- ii} else{
    a <- paste(a,ii,sep = ",")}
  }
  closeDialog()
  n <- tclvalue(Madays)
  command <- paste("rfaplot(lflist = list(",a,"),n =",n,")",sep = "")
  doItAndPrint(command)
  tkfocus(CommanderWindow())
} #End on OK


OKCancelHelp(helpSubject="rfaplot")
	tkgrid(getFrame(lfobjBox), sticky="nw")
        tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Annual n-day minima, n:")),MAentryframe,
               sticky = "w")
        tkgrid(optionsFrame, sticky = "w")
        tkgrid(buttonsFrame, sticky = "w")
        dialogSuffix(rows=3, columns=2)
      }
  
#Verteilung - Eingeben fehlt
rfacalc <- function(){
	initializeDialog(title=gettextRcmdr("Regional frequency analysis - Fitting procedure"))
        .activeModel <- ActiveModel()
        UpdateModelNumber()
        modelName <- tclVar(paste("RfaModel.", getRcmdr("modelNumber"), sep=""))
        modelFrame <- tkframe(top)
        model <- ttkentry(modelFrame, width="20", textvariable=modelName)
        optionsFrame <-  tkframe(top)
	lfobjBox <- variableListBox(top, listlfobj(),
		selectmode="multiple", title=gettextRcmdr("Select low flow objects: "))
        Madays <- tclVar(gettextRcmdr("1"))
        MAentryframe <- ttkentry(optionsFrame, width = "2", textvariable =Madays)
        choice <- c("Weibull", "GEV", "Lognormal","Gumbel","Pearson Type 3")
        init <-NULL

        for(ii in seq_along(getlfopt("extdist"))){
          init[ii] <- which(getlfopt("extdist")[ii] == choice)-1
        }

        distBox <- variableListBox(top, choice,
		selectmode="single", title=gettextRcmdr("Select distribution: "),
                           initialSelection = init,
                           listHeight = 5)
                  
onOK <-  function(){
  lfobjs <- getSelection(lfobjBox)
  if (length(lfobjs) == 0) {
    UpdateModelNumber(-1)
     errorCondition(recall=Recode, message=gettextRcmdr("You must select an object."))
			return()
		}
  a<- NULL
  for(ii in lfobjs){
    if(is.null(a)){a <- paste0(ii,"=",ii)} else{
    a <- paste0(a,",",ii,"=",ii)}
  }
  dist <- getSelection(distBox)
          options("RcmdrPlugin.lfstat" =
                        modifyList(getOption("RcmdrPlugin.lfstat"),list(extdist = dist)))
  modelValue <- trim.blanks(tclvalue(modelName))
        if (!is.valid.name(modelValue)){
            UpdateModelNumber(-1)
            errorCondition(recall=rfacalc, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue))
            return()
            }
        if (is.element(modelValue, listLinearModels())) {
            if ("no" == tclvalue(checkReplace(modelValue, type=gettextRcmdr("Model")))){
                UpdateModelNumber(-1)
                linearRegressionModel()
                return()
                }}
  closeDialog()
    names <- c("wei","gev","ln3","gum","pe3")
  distname <-paste0('"',names[which(dist == choice)],'"')
    
  n <- tclvalue(Madays)
  command <- paste0("rfa(lflist = list(",a,"),n =",n,", dist = ",distname,")")
  doItAndPrint(paste0(modelValue, " <- ", command))

  #gassign(modelValue, justDoIt(command))

  activateMenus()
  tkfocus(CommanderWindow())
 }#End on OK


OKCancelHelp(helpSubject="rfa")
        tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
        tkgrid(modelFrame, sticky="w")
	tkgrid(getFrame(lfobjBox), sticky="nw")
        tkgrid(labelRcmdr(optionsFrame, text = gettextRcmdr("Annual n-day minima, n:")),MAentryframe,
               sticky = "w")
        tkgrid(optionsFrame, sticky = "w")
        tkgrid(getFrame(distBox), sticky="nw")
        tkgrid(buttonsFrame, sticky = "w")
        dialogSuffix(rows=3, columns=2)
      }



      
rfaindex <- function(){
	initializeDialog(title=gettextRcmdr("Regional frequency analysis - Index values"))

rfaBox <- variableListBox(top, listrfd(),
		selectmode="single", title=gettextRcmdr("Select rfd object: "))
                   
onOK <-  function(){
  rfdobj <- getSelection(rfaBox)
  doItAndPrint(paste0(rfdobj,"$index"))
  closeDialog()
  tkfocus(CommanderWindow())
 }#End on OK
OKCancelHelp(helpSubject="rfa")
tkgrid(getFrame(rfaBox), sticky="nw")
tkgrid(buttonsFrame, sticky = "w")        
dialogSuffix(rows=3, columns=2)        
}
                       
listrfd <- function(envir=.GlobalEnv, ...) {
	objects <- ls(envir=envir, ...)
	if (length(objects) == 0) NULL
	else objects[sapply(objects,
		function(.x) "rfd" == (class(get(.x, envir=envir))[1]))]
}

rcgquantiles <- function(){
initializeDialog(title=gettextRcmdr("Regional Growth Curve - Quantiles"))

rfaBox <- variableListBox(top, listrfd(),
		selectmode="single", title=gettextRcmdr("Select rfd object: "))
textFrame <- tkframe(top)  
tyears <- tclVar(gettextRcmdr(getlfopt("extyears")))
tyearsentryframe <- ttkentry(textFrame, width="4", textvariable=tyears)

onOK <-  function(){
  rfdobj <- getSelection(rfaBox)
   event <- tclvalue(tyears)
   options("RcmdrPlugin.lfstat" =
                    modifyList(getOption("RcmdrPlugin.lfstat"),list(extyears = event)))
   doItAndPrint(paste0("regquant(",1/as.numeric(event),",",rfdobj,")"))
  
  closeDialog()
  tkfocus(CommanderWindow())
 }#End on OK
OKCancelHelp(helpSubject="regquant")
tkgrid(textFrame,sticky = "w")
tkgrid(labelRcmdr(textFrame, text = gettextRcmdr("Return period T (years):")),tyearsentryframe,sticky = "w")
tkgrid(getFrame(rfaBox), sticky="nw")
tkgrid(buttonsFrame, sticky = "w")
dialogSuffix(rows=3, columns=2)        
}

rcgsitequantiles <- function(){
initializeDialog(title=gettextRcmdr("Regional Growth Curve - Site quantiles"))

rfaBox <- variableListBox(top, listrfd(),
		selectmode="single", title=gettextRcmdr("Select rfd object: "))
textFrame <- tkframe(top)  
tyears <- tclVar(gettextRcmdr(getlfopt("extyears")))
tyearsentryframe <- ttkentry(textFrame, width="4", textvariable=tyears)

onOK <-  function(){
  rfdobj <- getSelection(rfaBox)
   event <- tclvalue(tyears)
   options("RcmdrPlugin.lfstat" =
                    modifyList(getOption("RcmdrPlugin.lfstat"),list(extyears = event)))
   doItAndPrint(paste0("sitequant(",1/as.numeric(event),",",rfdobj,")"))
  
  closeDialog()
  tkfocus(CommanderWindow())
 }#End on OK
OKCancelHelp(helpSubject="regquant")
tkgrid(textFrame,sticky = "w")
tkgrid(labelRcmdr(textFrame, text = gettextRcmdr("Return period T (years):")),tyearsentryframe,sticky = "w")
tkgrid(getFrame(rfaBox), sticky="nw")
tkgrid(buttonsFrame, sticky = "w")
dialogSuffix(rows=3, columns=2)     
}
   
