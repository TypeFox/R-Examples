Menu.summarize <- function(){
    ### is currently not used! 
    designs <- listDesigns()
    .activeDataSet <- ActiveDataSet()
  #  if ((length(designs) == 1) && !is.null(.activeDataSet) && designs==.activeDataSet) {
  #      Message(message=gettextRcmdr("There is only one design in memory."),
  #              type="warning")
  #      tkfocus(CommanderWindow())
  #      return()
  #      }
    if (length(designs) == 0){
        Message(message=gettextRcmdr("There are no design from which to choose."),
                type="error")
        tkfocus(CommanderWindow())
        return()
        }
    initializeDialog(title=gettextRcmdr("Summarize Design"))
    designsBox <- variableListBox(top, designs, title=gettextRcmdr("Designs (pick one)"),
        initialSelection=if (is.null(.activeDataSet) || !.activeDataSet %in% designs) NULL 
             else which(.activeDataSet == designs) - 1)
    plotcheck <- tclVar(FALSE)
    tablecheck <- tclVar(FALSE)
    
    cbFrame <- ttklabelframe(top, text="Additional output?")
    cbPlot <- ttkcheckbutton(cbFrame, text="Also create plot ?", variable=plotcheck)
    cbTable <- ttkcheckbutton(cbFrame, text="Also create table ?", variable=tablecheck)
    tkgrid(cbPlot, sticky="w")
    tkgrid(cbTable, sticky="w")
    
    
    onOK <- function(){   
        if (as.logical(as.numeric(tclvalue(plotcheck)))){ 
                    command <- paste("plot(",getSelection(designsBox),")")
        hilf <- justDoItDoE(command)
        if (class(hilf)[1]=="try-error"){
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=tab6, recall=Menu.summarize, message=gettextRcmdr(hilf))
             return()
            }
            }
        ## this is always done
        command <- paste("summary(",getSelection(designsBox),")")
        hilf <- doItAndPrint(command)
        if (class(hilf)[1]=="try-error"){
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=tab6, recall=Menu.summarize, message=gettextRcmdr(hilf))
             return()
            }
        
        if (as.logical(as.numeric(tclvalue(tablecheck)))){ 
                    name <- getSelection(designsBox) 
                    command <- paste("table(",name,"[,names(design.info(",name,")$factor.names)])")
        hilf <- doItAndPrint(command)
        if (class(hilf)[1]=="try-error"){
            Message(paste(gettextRcmdr("Offending command:"), "\n", command), type="error")
            errorCondition(window=tab6, recall=Menu.summarize, message=gettextRcmdr(hilf))
             return()
            }
        #logger(command)
        #print(hilf)
        }
        closeDialog()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Menu.summarize")
    tkgrid(getFrame(designsBox), cbFrame, sticky="n")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
    }
