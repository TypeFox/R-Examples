Menu.display <- function(){
    initializeDialog(title=gettextRcmdr("Display Active Design"))
    
    rbFrame <- ttklabelframe(top, text="Actual run order or standard order ?")
    rbVariable <- tclVar("actual")
    rbactual <- ttkradiobutton(rbFrame, text=gettextRcmdr("actual run order"), variable=rbVariable, value="actual")
    rbstd <- ttkradiobutton(rbFrame, text=gettextRcmdr("standard order"), variable=rbVariable, value="standard")
    tkgrid(rbactual, sticky="w")
    tkgrid(rbstd, sticky="w")
    
    onOK <- function(){   
        if (tclvalue(rbVariable)=="actual") display.design()
           else display.design.std.order()
        closeDialog()
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="Menu.display")
    tkgrid(rbFrame, sticky="n")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=2, columns=1)
    }
