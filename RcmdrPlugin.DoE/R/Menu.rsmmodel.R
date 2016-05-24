Menu.rsmmodel <- function(){
     ## menu for setting the overall scene (degree, response, coding)
    .activeDataSet <- ActiveDataSet()
     di <- design.info(eval(parse(text=.activeDataSet)))
     putRcmdr("resphilf", di$response.names[1])
     putRcmdr("resp.list", di$response.names)
     putRcmdr("factor.list", names(di$factor.names))
     putRcmdr("coded", as.numeric(!is.null(di$coding)))
     ### needed (among other things) for correction of automatic formula, 
                          ### if individual response of wide design is to be analyzed
    codedcbVariable <- tclVar(getRcmdr("coded"))
    degreeVar <- tclVar(2)
    initializeDialog(window=top, title=gettextRcmdr("Select Design, Degree, and Coding"))
    degreeEntry <- ttkentry(top, width=7, textvariable=degreeVar)
    degreelab <- ttklabel(top, text="degree (1, 1.5 or 2)")
    
    codingcb <- ttkcheckbutton(top,text=gettextRcmdr("Use coded values for calculations"),
                variable=codedcbVariable)
    
    putRcmdr("sel.resps", variableListBox(top, variableList=resp.list, listHeight=10, 
        title=gettextRcmdr("Response to be analysed (select one)"), 
        selectmode="single", initialSelection=0))
    putRcmdr("sel.factors", variableListBox(top, variableList=factor.list, listHeight=10, 
        title=gettextRcmdr("Factors for response surface model (select at least two)"), 
        selectmode="multiple", initialSelection=0:(di$nfactors-1)))
        
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
    onOK <- function(){
        response <- getSelection(sel.resps)
        fn <- getSelection(sel.factors)
        #fn is something like x1<A
        #the following gets the A from there
        #if (length(grep("<",fn)>0)) 
        #fn[grep("<",fn)] <- sapply(strsplit(fn[grep("<",fn)],"<", fixed=TRUE), function(obj) obj[2])
        putRcmdr("resphilf", response)
        putRcmdr("facthilf", paste("c(", paste(dquote(fn),collapse=","), ")", sep=""))

        putRcmdr("coded", as.numeric(tclvalue(codedcbVariable)))

        command <- paste("rsmformula(",.activeDataSet,", response=", dquote(response), 
              ", factor.names=c(", paste(dquote(fn),collapse=","), "), degree=", tclvalue(degreeVar), 
              ", coded=", as.logical(getRcmdr("coded")), ")")

        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error") {
            errorCondition(window=top,recall=Menu.rsmmodel, message=gettextRcmdr(hilf))
             return()
            }

        if (!is.null(di$responselist)){
          if (di$response.names[1]==di$responselist[1,1]){
          tk_messageBox(type="ok","This is a wide format repeated measurement or parameter design.\n
           You may want to aggregate it before applying a default analysis.\n
           Cancel the dialogue and go to \n
           Combine or Modify Designs --> Change from long to wide format \n
           for this purpose.")
          }
          } 
        if (di$repeat.only | (length(grep("param",di$type))>0 & is.null(di$responselist))){
          tk_messageBox(type="ok","This is a long format repeated measurement or parameter design.\n
           It has to be brought into wide format and subsequently should be aggregated before default analysis.\n
           Go to Combine or Modify Designs --> Change from long to wide format for this purpose.")
          closeDialog()
          tkfocus(CommanderWindow())
          } 
        else{
        if (length(grep("splitplot", di$type))>0) 
          warning("The default analysis can be misleading for splitplot designs!
          Whole-plot factors are often subject to much more variability than split-plot factors.")
          
          
          putRcmdr("degree", tclvalue(degreeVar))
          closeDialog()
          Menu.rsm(response=resphilf, factor.names=eval(parse(text=facthilf)))
        }
        
        }
    OKCancelHelp(helpSubject="Menu.rsmmodel")
    tkgrid(getFrame(sel.resps), getFrame(sel.factors), sticky="n")
    tkgrid(codingcb, columnspan=2)
    tkgrid(degreelab, degreeEntry, sticky="w")
    #tkgrid.configure(degreeEntry, sticky="nse")
    tkgrid.configure(degreelab, sticky="e")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(window=top, rows=2, columns=2)
    }
