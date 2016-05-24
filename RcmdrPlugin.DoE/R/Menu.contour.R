Menu.contour <- function(){
#### Aktualisierung der gewählten Auswahl-Form klappt noch nicht

putRcmdr("resphilf", all.vars(formula(ActiveModel()))[1])

## all x variables in model; unresolved issue: minus in formula
putRcmdr("varlistt", setdiff(all.vars(formula(ActiveModel())), resphilf))

## all numeric x variables in model
putRcmdr("varlistshortt", intersect(listNumeric(), varlistt))
                              ### numeric variables for plotting
                              ### inclusion of response possible (not good, but 
                              ###   perhaps better than alternatives)

if (length(varlistshortt)<=1) 
    tkmessageBox(title=gettextRcmdr("Too few numeric factors"), 
        text=gettextRcmdr("Model not suitable for contour plots"), 
        icon="error", type="ok")


### keep existing slice choices, if names from current model
if (exists("contours.at", where="RcmdrEnv"))
   if (!all(names(getRcmdr("contours.at") %in% varlistt))) 
           rm("contours.at", pos="RcmdrEnv")

### default slice choices
    putRcmdr("contours.default", lapply(get(getRcmdr(".activeDataSet"))[,varlistt], function(var) {
        if (is.factor(var))
            factor(levels(var)[1], levels = levels(var))
        else round(mean(var),2)
    }))

## set default choices
if (!exists("contours.at", where="RcmdrEnv")) putRcmdr("contours.at", contours.default)
if (!exists("comprrbVariable", where="RcmdrEnv")) putRcmdr("comprrbVariable", tclVar("compr"))
if (!exists("plottyperbVar", where="RcmdrEnv")) putRcmdr("plottyperbVar", tclVar("contour"))
if (!exists("imgcbVariable", where="RcmdrEnv")) putRcmdr("imgcbVariable", tclVar("0"))
if (!exists("colcbVariable", where="RcmdrEnv")) putRcmdr("colcbVariable", tclVar("0"))

onOK <- function(){
     closeDialog(window=top)
                ## compromise plans = structured plots
                if (tclvalue(comprrbVariable) == "compr") {
                    ## G1 is notest2fislist, G2 est2fislist
                    command <- as.character(parse(text=paste("fpairsc(c(",paste(shQuote(varlistshortt),collapse=","),") , c(",
                       paste(which(varlistshortt %in% notest2fislist),collapse=","),
                       "), ",substr(tclvalue(comprclassVar),1,1),")")))
                    estimable <- justDoItDoE(command)
                    if (class(estimable)[1]=="try-error") {
                          Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
                          errorCondition(window=top,recall=Menu.contour, message=gettextRcmdr(hilf))
                          return()
                    }
                }
                else estimable <- est2fislist
        ## determine slice positions
        hilf.cont <- vector("list", length(varlistt))
        names(hilf.cont) <- names(contours.default)
        for (i in 1:length(varlistt)){
             hilf.cont[[i]] <- tclvalue(getRcmdr(paste("cVar",i,sep="")))
             hilf <- contours.at[[i]]
             if (is.numeric(hilf)) hilf.cont[[i]] <- as.numeric(hilf.cont[[i]])
             if (is.factor(hilf)){
                 interim <- hilf.cont[[i]]
                 if (is.numeric(levels(hilf))) interim <- as.numeric(interim)
                 hilf.cont[[i]] <- hilf
                 hilf.cont[[i]][1] <- interim
                 }
             if (is.na(hilf.cont[[i]])){ 
                Message(gettextRcmdr("invalid slice positions, please correct"), type="error")
                errorCondition(window=top,recall=Menu.contour, message=gettextRcmdr("perhaps wrong decimal separator?"))
                 return()
                }
             }
        putRcmdr("contours.at", hilf.cont)
        nplots <- as.numeric(tclvalue(nrowVar))*as.numeric(tclvalue(ncolVar))

        ## where to plot slice info
        atpos <- 1                   ## sub title for each plot
        if (nplots>1) atpos <- 9     ## common sub title 

        ## preparation for multiple plots on one page
        if (nplots>1){
            command <- paste("par(mfrow=c(", tclvalue(nrowVar), ",", tclvalue(ncolVar),"), oma=c(3,0,2,0))")
            ## prepare graphics window for requested number row / column layout
            hilf <- justDoItDoE(command)
            if (class(hilf)[1]=="try-error") {
                Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
                errorCondition(window=top,recall=Menu.contour, message=gettextRcmdr(hilf))
                 return()
                }
            logger(command)
            }
        
        ## command for specifying slice positions
        if (!identical(getRcmdr("contours.default"), getRcmdr("contours.at")))
        textforcontoursat <- paste(", at = list(", paste(paste(names(contours.default),
                   paste(as.character(sapply(contours.at, function(obj) 
                             if (is.factor(obj)){if(is.character(levels(obj))) 
                                    shQuote(obj) else obj} else obj))), 
                                    sep="="), collapse=","), ")")
        else textforcontoursat <- ""

        ## define plot-type-specific commands
        if (tclvalue(plottyperbVar)=="contour") 
        command <- paste("contour(", activeModel(), ", as.list(c(", paste( paste("~", estimable, sep=""), collapse=","),
             ")), image=", as.logical(as.numeric(tclvalue(imgcbVariable))), textforcontoursat, ", atpos=", atpos,")")
        else
        if (tclvalue(plottyperbVar)=="persp"){
            if (as.logical(as.numeric(tclvalue(colcbVariable))))
              command <- paste("persp(", activeModel(), ", as.list(c(", paste( paste("~", estimable, sep=""), 
                                    collapse=","), ")), contours=\"col\", col=rainbow(50, end=5/6), atpos=", atpos, 
                                    textforcontoursat,")")
            else 
              command <- paste("persp(", activeModel(), ", as.list(c(", paste( paste("~", estimable, sep=""), 
                             collapse=","),")), atpos=", atpos, textforcontoursat,")")
             }
        else
        command <- paste("image(", activeModel(), ", as.list(c(", paste( paste("~", estimable, sep=""), 
                            collapse=","),")), atpos=", atpos, textforcontoursat,")")
        ## plot
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=top,recall=Menu.contour, message=gettextRcmdr(hilf))
             return()
            }
        logger(command)
        ## annotation for multiple plots
        if (nplots>1){
        command <- paste("mtext(",shQuote(paste("Response surface plots for ", resphilf, sep="")), ", outer=TRUE)")
        logger(command)
        justDoItDoE(command)
        command <- paste("mtext(",shQuote(paste("Slice at: ", paste(paste(names(contours.at), contours.at, sep="="), 
                                collapse=", "))), ", outer=TRUE, side=1)")
        logger(command)
        justDoItDoE(command)
        }
        ## remove variables for which there is a risk of confusion with other dialogs
        rm("resphilf", pos="RcmdrEnv")
        rm("varlistt", pos="RcmdrEnv")
        rm("varlistshortt", pos="RcmdrEnv")
 ### omitted following two rows in order to have graph on top
        ## tkwm.deiconify(CommanderWindow())
        ## tkfocus(CommanderWindow())
  }


dquote <- function(obj){
    ## quote vector elements for use as character vector in a command
    aus <- rep("",length(obj))
    wopt <- options("warn")[[1]]
    options(warn=-1)
    for (i in 1:length(obj)) if (is.na(as.numeric(obj[i]))) {
            if (length(grep('"',obj[i])>0))
            aus[i] <- paste("'",obj[i],"'",sep="") 
            else
            aus[i] <- paste('"',obj[i],'"',sep="") 
            }
          else aus[i] <- obj[i]
    options(warn=wopt)
    aus
}



oncomprestrb <- function(){
        oncomprestrb.worefresh
    #    tkconfigure(notest2fis$listbox, listvariable=tclVar(paste(notest2fislist,collapse=" ")))
    #    notest2fis$varlist <- notest2fislist
    #    tkconfigure(est2fis$listbox, listvariable=tclVar(paste(est2fislist,collapse=" ")))
    #    est2fis$varlist <- est2fislist
          tkfocus(CommanderWindow())
          tkdestroy(top)
          Menu.contour()

}
oncomprestrb.worefresh <- function(){
        if (tclvalue(comprrbVariable)=="compr"){
             tkconfigure(comprclassEntry, state="normal")
             }
        else 
             tkconfigure(comprclassEntry, state="disabled")
        
        if (tclvalue(comprrbVariable)=="manual"){
          hilf <- combn(length(varlistshortt),2)
          putRcmdr("intaclistt", paste(varlistshortt[hilf[1,]],varlistshortt[hilf[2,]],sep="*"))
          }
        else{
          putRcmdr("intaclistt", varlistshortt)
        }        
          putRcmdr("est2fislist", "")
          putRcmdr("notest2fislist", intaclistt)
}


 onSelect <- function(){
     if (tclvalue(tcl(notest2fis$listbox, "curselection"))=="") return()
     ## curselection is a character string with blank separated selection positions
     add <- notest2fislist[as.numeric(unlist(strsplit(tclvalue(tcl(notest2fis$listbox, "curselection")), " ")))+1]
     putRcmdr("notest2fislist", setdiff(notest2fislist,add))
     putRcmdr("est2fislist", intaclistt[intaclistt %in% c(est2fislist,add)])
     add <- NULL
     tcl(notest2fis$listbox, "selection", "clear", "0", "999")
     tkconfigure(notest2fis$listbox, listvariable=tclVar(paste(notest2fislist,collapse=" ")))
     notest2fis$varlist <- notest2fislist
     tkconfigure(est2fis$listbox, listvariable=tclVar(paste(est2fislist,collapse=" ")))
     est2fis$varlist <- est2fislist
 }
 onDeselect <- function(){
     ## curselection is a character string with blank separated selection positions
     if (tclvalue(tcl(est2fis$listbox, "curselection"))=="") return()
     add <- est2fislist[as.numeric(unlist(strsplit(tclvalue(tcl(est2fis$listbox, "curselection")), " ")))+1]
     putRcmdr("est2fislist", setdiff(est2fislist,add))
     putRcmdr("notest2fislist", intaclistt[intaclistt %in% c(notest2fislist,add)])
     add <- NULL
     tcl(est2fis$listbox, "selection", "clear", "0", "999")
     tkconfigure(notest2fis$listbox, listvariable=tclVar(paste(notest2fislist,collapse=" ")))
     notest2fis$varlist <- notest2fislist
     tkconfigure(est2fis$listbox, listvariable=tclVar(paste(est2fislist,collapse=" ")))
     est2fis$varlist <- est2fislist
 }
######## end define functions                          


##### define userform
initializeDialog(title=gettextRcmdr("Contour plots ..."))   

putRcmdr("tn",ttknotebook(top))
#tn <- ttknotebook(topdes2)

putRcmdr("tab1",ttkframe(tn))
putRcmdr("tab2",ttkframe(tn))

tkadd(tn,tab1,text="Base Settings")   ### tabid=0
tkadd(tn,tab2,text="Modify slice positions")  ### tabid=1

tkconfigure(tn, takefocus=0)

baseFrame <- ttklabelframe(tab1,text=gettextRcmdr("Properties of plots"))

## inform about response
tkgrid(ttklabel(top, text=gettextRcmdr("The response analysed:")), ttklabel(top, text=resphilf), sticky="w")

## for mfrow
layoutFrame <- ttklabelframe(baseFrame, text=gettextRcmdr("Select plot layout (rows and columns)"))
if (!exists("nrowVar", where="RcmdrEnv")) putRcmdr("nrowVar", tclVar("1"))
if (!exists("ncolVar", where="RcmdrEnv")) putRcmdr("ncolVar", tclVar("1"))
nrowLabel <- tklabel(layoutFrame, text=gettextRcmdr("Number of rows"))
nrowEntry <- tkentry(layoutFrame, width="8", textvariable=nrowVar)
ncolLabel <- tklabel(layoutFrame, text=gettextRcmdr("Number of columns"))
ncolEntry <- tkentry(layoutFrame, width="8", textvariable=ncolVar)
tkgrid(nrowLabel, nrowEntry, sticky="w")
tkgrid(ncolLabel, ncolEntry, sticky="w")


## choice of plot function
typeFrame <- ttklabelframe(baseFrame, text=gettextRcmdr("Select plot type"))
contourrb <- tkradiobutton(typeFrame,text=gettextRcmdr("Contour plot"),variable=plottyperbVar,value="contour")
persprb <- tkradiobutton(typeFrame,text=gettextRcmdr("3D perspective plot"),variable=plottyperbVar,value="persp")
imagerb  <- tkradiobutton(typeFrame,text=gettextRcmdr("Image plot"),variable=plottyperbVar,value="image")
imgcb <- tkcheckbutton(typeFrame, text=gettextRcmdr("Overlay image on contours"), variable=imgcbVariable)
colcb <- tkcheckbutton(typeFrame, text=gettextRcmdr("Colour the surface"), variable=colcbVariable)
tkgrid(contourrb, imgcb, sticky="w")
tkgrid(persprb, colcb, sticky="w")
tkgrid(imagerb, sticky="w")

## choice of plotted pairs

selFrame <- ttklabelframe(baseFrame, text=gettextRcmdr("Select pairs for plotting"))
comprradioFrame <- ttklabelframe(selFrame, text=gettextRcmdr("Type of specification"))

## design from older version of RcmdrPlugin.DoE
comprestrb <- tkradiobutton(comprradioFrame,text=gettextRcmdr("Pairs formed from two groups of variables"),
      variable=comprrbVariable,value="compr",wraplength="500",justify="left",command=oncomprestrb)
manualestrb <- tkradiobutton(comprradioFrame,text=gettextRcmdr("Select manually"),
      variable=comprrbVariable,value="manual",command=oncomprestrb)

## design from older version of RcmdrPlugin.DoE
if (!exists("comprclassVar", where="RcmdrEnv")) 
    putRcmdr("comprclassVar", tclVar("1: all pairs within group 1"))

#resEntry <- tkentry(despropframe, textvariable=resVar)
    putRcmdr("comprclassEntry", ttkcombobox(comprradioFrame, textvariable=comprclassVar, width=50, 
         values=c("1: pairs within group1",
                  "2: pairs within groups 1 and groups 2",
                  "3: all pairs of group 1", 
                  "4: pairs between groups 1 and 2"), state="readonly"))
    #tkbind(comprclassEntry, "<<ComboboxSelected>>", onRefresh)

tkgrid(comprestrb, comprclassEntry, sticky="w")
tkgrid(manualestrb, sticky="w")
tkgrid(comprradioFrame, sticky="w", columnspan=6)

## initialize variables intaclistt, est2fislist and notest2fislist
## configure 
oncomprestrb.worefresh()

estbuttonFrame <- ttkframe(selFrame)
selectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr(">"), 
        foreground = "darkgreen", command = onSelect, 
        default = "normal", borderwidth = 3)
tkgrid(selectButton)
deselectButton <- buttonRcmdr(estbuttonFrame, text = gettextRcmdr("<"), 
        foreground = "darkgreen", command = onDeselect, 
        default = "normal", borderwidth = 3)
tkgrid(deselectButton)


## define variable list boxes for selection
## intaclistt is the master list

## make sure that both lists are of equal length and long enough
## ## make selection offering depend on choice of compromise settings
     if (tclvalue(comprrbVariable)=="compr"){
        TITEL.LINKS <- "Group 1"
        TITEL.RECHTS <- "Group 2"
     }
     else {
        TITEL.LINKS <- "Available pairs"
        TITEL.RECHTS <- "Selected pairs"
     }
putRcmdr("est2fis", variableListBox(selFrame, variableList=intaclistt, listHeight=15, 
    title=TITEL.RECHTS, selectmode="multiple"))
putRcmdr("notest2fis", variableListBox(selFrame, variableList=intaclistt, listHeight=15, 
    title=TITEL.LINKS, selectmode="multiple"))
     tkconfigure(notest2fis$listbox, listvariable=tclVar(paste(notest2fislist,collapse=" ")))
     notest2fis$varlist <- notest2fislist
     tkconfigure(est2fis$listbox, listvariable=tclVar(paste(est2fislist,collapse=" ")))
     est2fis$varlist <- est2fislist


tkgrid(layoutFrame, typeFrame, sticky="w")
tkgrid(baseFrame, sticky="w", columnspan=2)
tkgrid(notest2fis$frame, estbuttonFrame, est2fis$frame, sticky="w")
tkgrid(selFrame, sticky="ew", columnspan=6, columnspan=2)

## tab 2

contFrame <- ttklabelframe(tab2, text=gettextRcmdr("Modify values at which to fix variables that are not in the current plot"))
tkgrid(ttklabel(tab2, text=gettextRcmdr("Variable"), foreground="blue"),
       ttklabel(tab2, text=gettextRcmdr("Value (one in each field)"), foreground="blue"),
       sticky="w",pady=c(10,5))
for (i in 1:length(getRcmdr("varlistt"))){
    putRcmdr(paste("lab", i, sep=""), tklabel(tab2, text=varlistt[i]))
    putRcmdr(paste("cVar", i, sep=""), tclVar(contours.at[[i]]))
    if (!is.factor(contours.at[[i]])) 
        putRcmdr(paste("contEntry", i, sep=""), tkentry(tab2, textvariable= getRcmdr(paste("cVar",i,sep=""))))
    else putRcmdr(paste("contEntry", i, sep=""), 
          ttkcombobox(tab2, textvariable=getRcmdr(paste("cVar",i,sep="")), values=levels(contours.at[[i]]), state="readonly"))
    tkgrid(getRcmdr(paste("lab",i,sep="")), getRcmdr(paste("contEntry",i,sep="")), sticky="w")
    tkgrid.configure(getRcmdr(paste("lab",i,sep="")), sticky="e")
        }
tkgrid(tn)


OKCancelHelp(window=top, helpSubject="Menu.contour")

tkgrid(buttonsFrame, sticky="ew", columnspan=2)


dialogSuffix(window=top, rows=2, columns=2, bindReturn=FALSE)

}
