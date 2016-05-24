fix.design <- function(x,..., prompt=FALSE){
    tkmessageBox(title="Info for editing design objects",
       message="It is not permitted to edit any factor variables.\nIf you change anything for any factor variable, ALL changes will be lost.", 
       type="ok", icon="info")
    ## attributes are preserved
    ##    add response names element to design.info
    ##    bugfix prompt
    pickresps <- function (varliste){
        result <- NULL
        ## functions internal in the prompt function
        getSelection <- function(object) varliste[as.numeric(tkcurselection(yBox)) + 1]
        onOK <- function(){
          hilf <- getSelection(yBox)
          command <- paste("result <- hilf")          
          ## assign("result", hilf, envir = .GlobalEnv)
          justDoIt(command)
          tkdestroy(selresp)
        }
        onCancel <- function() tkdestroy(selresp)
        onHelp <- function() {
                if (.Platform$OS.type != "windows") tkgrab.release(window)
                print(help("fix.design"))
            }

        ## open the prompt window
        selresp <- tktoplevel(borderwidth=10)
#        tkwm.withdraw(window)
        tkwm.title(selresp, gettext("Select response(s)"))
        yBoxlab <- tklabel(selresp, text=gettext("Which of the variables are response (pick one or more)?"))
        yBox <- tklistbox(selresp, listvariable=tclVar(varliste), height=min(10, length(varliste)), selectmode = "extended",
              width=max(20, nchar(varliste)))
        scrollbar <- ttkscrollbar(selresp, command=function(...) tkyview(yBox, ...))
        tkconfigure(yBox, yscrollcommand=function(...) tkset(scrollbar, ...))
        for (var in varliste) tkinsert(yBox, "end", var)
        
        tkgrid(yBoxlab)
        tkgrid(yBox)

        buttonsFrame <- tkframe(selresp, borderwidth=5)
        OKbutton <- tkbutton(buttonsFrame, text="OK", width="12", 
            command=onOK, default="active", borderwidth=3)
        cancelButton <- tkbutton(buttonsFrame, text=gettext("Cancel"), 
            width="12", command=onCancel, borderwidth=3)
        helpButton <- tkbutton(buttonsFrame, text=gettext("Help"), width="12", command=onHelp, borderwidth=3)
        tkgrid(OKbutton, tklabel(buttonsFrame, text="  "), cancelButton, tklabel(buttonsFrame, text="            "),
            helpButton, sticky="w")

        tkgrid(buttonsFrame)
        tkwm.resizable(selresp, 1, 1)
        for (row in 0:(5 - 1)) tkgrid.rowconfigure(selresp, row, 
                  weight = 0)
        for (col in 0:(3 - 1)) tkgrid.columnconfigure(selresp, 
                  col, weight = 0)
              .Tcl("update idletasks")
        tkwm.deiconify(selresp)
        tkgrab.set(selresp)
        tkfocus(yBox)
        tkwait.window(selresp)
        tclServiceMode(on=TRUE)
        result
    }
    
    y <- x
    subx <- substitute(x)
    if (is.name(subx))
        subx <- deparse(subx)
    if (!is.character(subx) || length(subx) != 1L)
        stop("'fix' requires a name")
    parent <- parent.frame()
    if (exists(subx, envir = parent, inherits = TRUE))
        x <- edit(get(subx, envir = parent), title = subx, ...)
    else {
        x <- edit(function() {
        }, title = subx, ...)
        environment(x) <- .GlobalEnv
    }
    if (identical(rownames(x),rownames(y))){
            newnam <- setdiff(colnames(x),colnames(y))
            if (length(newnam)>0){
              numnam <- newnam[which(sapply(as.list(x)[newnam], "is.numeric"))]
              if (prompt & length(numnam)>0) respnam <- pickresps(numnam)
                else if (length(numnam) > 0) respnam <- numnam
                  else respnam <- NULL
              xnumnew <- data.frame(x[,newnam])
              colnames(xnumnew) <- newnam
              if (length(newnam) > length(numnam)) 
                  xnumnew[,setdiff(newnam,numnam)] <- 
                  lapply(xnumnew[setdiff(newnam,numnam)], function(obj) as.numeric(factor(obj)))
              attr(x,"desnum") <- cbind(desnum(y), as.matrix(xnumnew))
              hilf <- design.info(y)
              if (length(respnam)>0) hilf$response.names <- c(hilf$response.names, respnam)
            }
            else{
                attr(x,"desnum") <- desnum(y)
                hilf <- design.info(y)
            }
            attr(x,"run.order") <- run.order(y)
            attr(x,"design.info") <- hilf
            class(x) <- class(y)
            ## I did not manage to make justDoIt work here (instead of assign)
            ## perhaps because fix is used inside of a justDoIt
            ## the two ## assign(subx, x, envir = .GlobalEnv) ##
            ## were changed to hidden assigments to the global environment for now
            ## G <- .GlobalEnv; assign(subx, x, envir = G) works
            ## equivalent to gassign in Rcmdr
            .G <- .GlobalEnv
            if (all(names(factor.names(y)) %in% colnames(x) & names(factor.names(y)) %in% colnames(y))){
              if (all(x[,names(factor.names(y))]==y[,names(factor.names(y))]))
                    assign(subx, x, envir = .G) 
              else warning("changes have not been stored, as it is not permitted to edit factor variables")
              }
            else warning("changes have not been stored, as it is not permitted to rename or delete factor variables")
         }
    else{
        antwort <- tkmessageBox(title=gettext("Really save changes ?"),
            type="yesno", message="The row names have changed. The design will loose all its properties. Do you really want to save these changes?")
        if (tclvalue(antwort)=="yes") assign(subx, x, envir = .G)
             else (message("The changes were not saved due to a user decision."))
    }
}