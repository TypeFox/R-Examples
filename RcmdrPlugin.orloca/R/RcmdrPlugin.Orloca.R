# Some Rcmdr dialogs for the orloca package (non graphical functions)

.onAttach <- function(libname, pkgname){
   if (!interactive()) return()
   Rcmdr <- options()$Rcmdr
   plugins <- Rcmdr$plugins
   options(list(".RcmdrPlugin.orloca.l2" = T))
   options(list(".RcmdrPlugin.orloca.lp"  = NA))
   if (!pkgname %in% plugins) {
       Rcmdr$plugins <- c(plugins, pkgname)
       options(Rcmdr=Rcmdr)
       if("package:Rcmdr" %in% search()) {
           if(!getRcmdr("autoRestart")) {
               closeCommander(ask=FALSE, ask.save=TRUE)
               Commander()
           }
       }
       else {
           Commander()
       }
   }
} 

.RcmdrPlugin.orloca.get.norma <- function(sep=",")
  {
    l2 <- options(".RcmdrPlugin.orloca.l2")
    command <- ""
    if (l2 != TRUE)
       {
       lp <- options(".RcmdrPlugin.orloca.lp")
       command <- paste(sep, " lp = ", lp, sep="")
       }
    command
  }

###
# New loca.p object as data.frame
###
Rcmdr.new.loca.p <- function(){
initializeDialog(title=gettext("New loca.p", domain="R-RcmdrPlugin.orloca"))
nameVar <- tclVar(gettextRcmdr("Data"))
nameEntry <- tkentry(top, width="8", textvariable=nameVar)
onOK <- function(){
	closeDialog()
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)) {
            errorCondition(recall=Rcmdr.new.loca.p, message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (is.element(name, listDataSets()))
          {
          if ("no" == tclvalue(checkReplace(name, gettextRcmdr("Data set"))))
            {
              errorCondition(recall=Rcmdr.new.loca.p, message=gettextRcmdr("Introduce the name (another) for the new data.frame."))
              return()
             }
          }
        command <- paste(name, "<- edit(data.frame(x=numeric(0), y=numeric(0), w=numeric(0)))")
        doItAndPrint(command)
        if (nrow(get(name)) == 0){
            errorCondition(recall=Rcmdr.new.loca.p, message=gettextRcmdr("empty data set."))
            return()
            }
        activeDataSet(name)
        closeDialog()
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject="loca.p")
tkgrid(tklabel(top, text=gettext("Name of new loca.p object", domain="R-RcmdrPlugin.orloca")), nameEntry, sticky="e")
tkgrid(buttonsFrame, sticky="w", columnspan=2)
tkgrid.configure(nameEntry, sticky="w")
dialogSuffix(rows=2, columns=2, focus=nameEntry)
}
###
# End of Rcmdr.new.loca.p
###

Rcmdr.rloca.p <- function(){
initializeDialog(title=gettext("New loca.p Random Instance", domain="R-RcmdrPlugin.orloca"))
nameVar <- tclVar(gettextRcmdr("Data"))
nameEntry <- tkentry(top, width="8", textvariable=nameVar)
nVar <- tclVar("100")
nEntry <- tkentry(top, width="8", textvariable=nVar)
xminVar <- tclVar("0")
xminEntry <- tkentry(top, width="8", textvariable=xminVar)
xmaxVar <- tclVar("1")
xmaxEntry <- tkentry(top, width="8", textvariable=xmaxVar)
yminVar <- tclVar("0")
yminEntry <- tkentry(top, width="8", textvariable=yminVar)
ymaxVar <- tclVar("1")
ymaxEntry <- tkentry(top, width="8", textvariable=ymaxVar)
groupsVar <- tclVar("0")
groupsEntry <- tkentry(top, width="8", textvariable=groupsVar)
onOK <- function(){
	closeDialog()
        name <- tclvalue(nameVar)
        if (!is.valid.name(name)) {
            errorCondition(recall=Rcmdr.rloca.p, message=paste('"', name, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
          }
        if (is.element(name, listDataSets()))
          {
          if ("no" == tclvalue(checkReplace(name, gettextRcmdr("Data set"))))
            {
              errorCondition(recall=Rcmdr.rloca.p, message=gettextRcmdr("Introduce the name (another) for the new data.frame."))
              return()
             }
          }
	n <- round(as.numeric(tclvalue(nVar)))
        if (is.na(n) || n <= 0){
            errorCondition(recall=Rcmdr.rloca.p, message=gettext("The number of demand points must be a positive integer.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	xmin <- as.numeric(tclvalue(xminVar))
        if (is.na(xmin)){
            errorCondition(recall=Rcmdr.rloca.p, message=gettext("xmin must be a real number.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	xmax <- as.numeric(tclvalue(xmaxVar))
        if (is.na(xmax) || xmax < xmin){
            errorCondition(recall=Rcmdr.rloca.p, message=gettext("xmax must be a real number bigger that xmin.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	ymin <- as.numeric(tclvalue(yminVar))
        if (is.na(ymin)){
            errorCondition(recall=Rcmdr.rloca.p, message=gettext("ymin must be a real number.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	ymax <- as.numeric(tclvalue(ymaxVar))
        if (is.na(ymax) || ymax < ymin){
            errorCondition(recall=Rcmdr.rloca.p, message=gettext("ymax must be a real number bigger that ymin.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	groups <- round(as.numeric(tclvalue(groupsVar)))
        if (is.na(groups) || groups < 0){
            errorCondition(recall=Rcmdr.rloca.p, message=gettext("groups must be a non negative integer.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	command <- paste(name, " <- rloca.p(n = ", n,", xmin = ", xmin,", xmax = ", xmax,", ymin = ", ymin,", ymax = ", ymax,", groups = ", groups,")", sep="")
	doItAndPrint(command)
#        command <- paste("summary(", name, ")")
#        doItAndPrint(command)
	command <- paste(name, " <- as(", name, ", \"data.frame\")", sep="")
        doItAndPrint(command)
        activeDataSet(name)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject="rloca.p")
tkgrid(tklabel(top, text=gettext("Name of new loca.p object", domain="R-RcmdrPlugin.orloca")), nameEntry, sticky="e")
tkgrid(tklabel(top, text=gettext("Number of demand points", domain="R-RcmdrPlugin.orloca")), nEntry, sticky="e")
tkgrid(tklabel(top, text=gettext("x Minimum", domain="R-RcmdrPlugin.orloca")), xminEntry, sticky="e")
tkgrid(tklabel(top, text=gettext("x Maximum", domain="R-RcmdrPlugin.orloca")), xmaxEntry, sticky="e")
tkgrid(tklabel(top, text=gettext("y Minimum", domain="R-RcmdrPlugin.orloca")), yminEntry, sticky="e")
tkgrid(tklabel(top, text=gettext("y Maximum", domain="R-RcmdrPlugin.orloca")), ymaxEntry, sticky="e")
tkgrid(tklabel(top, text=gettext("Number of groups", domain="R-RcmdrPlugin.orloca")), groupsEntry, sticky="e")
tkgrid(buttonsFrame, sticky="w", columnspan=2)
tkgrid.configure(nameEntry, sticky="w")
tkgrid.configure(nEntry, sticky="w")
tkgrid.configure(xminEntry, sticky="w")
tkgrid.configure(xmaxEntry, sticky="w")
tkgrid.configure(yminEntry, sticky="w")
tkgrid.configure(ymaxEntry, sticky="w")
tkgrid.configure(groupsEntry, sticky="w")
dialogSuffix(rows=7, columns=2, focus=nEntry)
}


Rcmdr.zsum <- function(){
initializeDialog(title=gettext("Evaluation of Objective Function for weighted sum Location Problem", domain="R-RcmdrPlugin.orloca"))
xVar <- tclVar("0")
xEntry <- tkentry(top, width="6", textvariable=xVar)
yVar <- tclVar("0")
yEntry <- tkentry(top, width="6", textvariable=yVar)
onOK <- function(){
	closeDialog()
	x <- as.numeric(tclvalue(xVar))
        if (is.na(x)){
            errorCondition(recall=Rcmdr.zsum, message=gettext("x-axis must be a number.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	y <- as.numeric(tclvalue(yVar))
        if (is.na(y)){
            errorCondition(recall=Rcmdr.zsum, message=gettext("y-axis must be a number.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
        command <- paste("zsum(as(", ActiveDataSet(), ", \"loca.p\") , x = ", x,", y = ", y, sep="")
        command <- paste(command, .RcmdrPlugin.orloca.get.norma(), sep="")
        command <- paste(command, ") # ", sep="")
        command <- paste(command, gettext("Weighted sum of distances", domain="R-RcmdrPlugin.orloca"), sep="")
	doItAndPrint(command)
        command <- paste("zsumgra(as(", ActiveDataSet(), ", \"loca.p\") , x = ", x,", y = ", y, sep="")
        command <- paste(command, .RcmdrPlugin.orloca.get.norma(), sep="")
        command <- paste(command, ") # ", sep="")
        command <- paste(command, gettext("Gradient of the weighted sum of distances function", domain="R-RcmdrPlugin.orloca"), sep="")
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject="zsum")
tkgrid(tklabel(top, text="x-axis"), xEntry, sticky="e")
tkgrid(tklabel(top, text="y-axis"), yEntry, sticky="e")
tkgrid(buttonsFrame, sticky="w", columnspan=2)
tkgrid.configure(xEntry, sticky="w")
tkgrid.configure(yEntry, sticky="w")
dialogSuffix(rows=2, columns=2, focus=xEntry)
}

Rcmdr.zsummin <- function(){
initializeDialog(title=gettext("Solve weighted sum Location Problem", domain="R-RcmdrPlugin.orloca"))
xVar <- tclVar("0")
xEntry <- tkentry(top, width="6", textvariable=xVar)
yVar <- tclVar("0")
yEntry <- tkentry(top, width="6", textvariable=yVar)
nVar <- tclVar("100")
nEntry <- tkentry(top, width="6", textvariable=nVar)
epsVar <- tclVar("0.001")
epsEntry <- tkentry(top, width="6", textvariable=epsVar)
radioButtons(name="algorithm", buttons=c("w", "g", "s", "u"), values=c("w", "g", "s", "u"), initialValue="w", labels=gettext(c("Weiszfeld", "Gradient", "Search method", "ucminf"), domain="R-RcmdrPlugin.orloca"), title=gettext("Select algorithm", domain="R-RcmdrPlugin.orloca"))

#tkgrid(labelRcmdr(statisticFrame), sticky="w")
onOK <- function(){
	closeDialog()
	x <- as.numeric(tclvalue(xVar))
        if (is.na(x)){
            errorCondition(recall=Rcmdr.zsummin, message=gettext("x-axis must be a number.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	y <- as.numeric(tclvalue(yVar))
        if (is.na(y)){
            errorCondition(recall=Rcmdr.zsummin, message=gettext("y-axis must be a number.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	n <- as.numeric(tclvalue(nVar))
        if (is.na(n) || n <= 0){
            errorCondition(recall=Rcmdr.zsummin, message=gettext("The maximum number of iterations must be a positive integer", domain="R-RcmdrPlugin.orloca"))
            return()
            }
	eps <- as.numeric(tclvalue(epsVar))
        if (is.na(eps) || eps <= 0){
            errorCondition(recall=Rcmdr.zsummin, message=gettext("The norm of the gradient must be positive.", domain="R-RcmdrPlugin.orloca"))
            return()
            }
        algorithm <- tclvalue(algorithmVariable)
        command <- paste(".sol <- zsummin(as(", ActiveDataSet(), ", \"loca.p\") , x = ", x,", y = ", y,", eps =", eps, ", algorithm =\"", algorithm, "\"", sep="")
        command <- paste(command, .RcmdrPlugin.orloca.get.norma(), " ) # ", gettext("Solve the minsum location problem", domain="R-RcmdrPlugin.orloca"), sep="")
	doItAndPrint(command)
        doItAndPrint(paste(".sol #", gettext("Show the solution", domain="R-RcmdrPlugin.orloca")))
        command <- paste("zsum(as(", ActiveDataSet(), ", \"loca.p\") , x =", .sol[1], ", y = ", .sol[2], wep="")
        command <- paste(command, .RcmdrPlugin.orloca.get.norma(), sep="")
        command <- paste(command, ") # ", gettext("Weighted sum of distances", domain="R-RcmdrPlugin.orloca"), sep="")
        doItAndPrint(command)
        doItAndPrint("remove(.sol)")
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject="zsummin")
tkgrid(tklabel(top, text="Iter. Max."), nEntry, sticky="w")
tkgrid.configure(nEntry, sticky="e")
tkgrid(tklabel(top, text="x-axis"), xEntry, sticky="w")
tkgrid.configure(xEntry, sticky="e")
tkgrid(tklabel(top, text="y-axis"), yEntry, sticky="w")
tkgrid.configure(yEntry, sticky="e")
tkgrid(tklabel(top, text="Grad. norm."), epsEntry, sticky="w")
tkgrid.configure(epsEntry, sticky="e")
tkgrid.configure(algorithmFrame, sticky="w")
tkgrid(algorithmFrame, sticky="w", columnspan=2)
tkgrid.configure(buttonsFrame, sticky="w")
tkgrid(buttonsFrame, sticky="w", columnspan=2)
dialogSuffix(rows=20, columns=2, focus=xEntry)
}



Rcmdr.help.orloca <- function(){
   # To ensure that menu name is included in pot file
   gettext("Help about orloca", domain="R-RcmdrPlugin.orloca")
   command <- paste("help(\"", gettext("orloca", domain="R-orloca"), sep="")
   command <- paste(command, "\")", sep="")
   doItAndPrint(command)
   invisible(NULL)
}

Rcmdr.help.RcmdrPlugin.orloca <- function(){
   # To ensure that menu name is included in pot file
   gettext("Help about RcmdrPlugin.orloca", domain="R-RcmdrPlugin.orloca")
   command <- paste("help(\"", gettext("RcmdrPlugin.orloca", domain="R-RcmdrPlugin.orloca"), sep="")
   command <- paste(command, "\")", sep="")
   doItAndPrint(command)
   invisible(NULL)
}

Rcmdr.summary.loca.p <- function(){
   # To ensure that menu name is included in pot file
   gettext("Summary", domain="R-RcmdrPlugin.orloca")
   command <- paste("summary(as(", ActiveDataSet(), ", \"loca.p\"))", sep="")
   doItAndPrint(command)
   invisible(NULL)
}

activeDataSetLocaP <- function() activeDataSetP() && validObject(new("loca.p",x=get(ActiveDataSet())$x, y=get(ActiveDataSet())$y, w=get(ActiveDataSet())$w))

activeDataSetLocaP <- function()
  {
    if (activeDataSetP())
      {
      .activeDataSet <- get(ActiveDataSet())
      (nrow(.activeDataSet)==length(.activeDataSet$x)) && (nrow(.activeDataSet)==length(.activeDataSet$y)) && (nrow(.activeDataSet)==length(.activeDataSet$w)) && (sum(is.na(.activeDataSet$x))+sum(is.na(.activeDataSet$y)+sum(is.na(.activeDataSet$w)))==0)
      }
    else FALSE
  }

Rcmdr.orloca.norm <- function(){
  # To ensure that menu name is include in pot file
  gettext("Show/Set norm", domain="R-RcmdrPlugin.orloca")
  # Leemos los valores actuales
  l2 <- options(".RcmdrPlugin.orloca.l2")
  if (l2 == TRUE)
    {
      lp <- ""
      iv <- "l2"
    }
  else
    {
      lp <- as.character(options(".RcmdrPlugin.orloca.lp"))
      iv <- 'lp'
    }
initializeDialog(title=gettext("Selection of the norm", domain="R-RcmdrPlugin.orloca"))
radioButtons(name="norma", title= gettext("Select the norm", domain="R-RcmdrPlugin.orloca"), buttons=c("l2", "lp"), labels=gettext(c("l_2 ", "l_p "), domain="R-RcmdrPlugin.orloca"), values=c("l2", "lp"), initialValue=iv)

nameVar <- tclVar(lp)
nameEntry <- tkentry(top, width="8", textvariable=nameVar)
onOK <- function(){
	closeDialog()
        name <- as.numeric(tclvalue(nameVar))
        on <- tclvalue(normaVariable)
        if (identical(on, 'l2')) options(list(".RcmdrPlugin.orloca.l2" = T))
        else if (name >= 1)
            {
            options(list(".RcmdrPlugin.orloca.l2" = F))
            options(list(".RcmdrPlugin.orloca.lp" = name))
            tkfocus(CommanderWindow())
            }
          else
            {
            errorCondition(recall=Rcmdr.orloca.norm, message=paste('"', name, '" ', gettext("is not a valid l_p norm.", domain="R-RcmdrPlugin.orloca"), sep=""))
            }
        return()
	}
OKCancelHelp(helpSubject="zsum")
tkgrid(normaFrame, sticky="w")
tkgrid(tklabel(top, text=gettext("p = ", domain="R-RcmdrPlugin.orloca")), nameEntry, sticky="e")
tkgrid(buttonsFrame, sticky="w", columnspan=2)
tkgrid.configure(nameEntry, sticky="w")
dialogSuffix(rows=3, columns=2, focus=normaFrame)
}
