arMODEL= function (input="",...)   {

listDataSets = function (envir = .GlobalEnv, ...){
    Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
    out=names(which(sapply(Vars, function(.x) is.data.frame(get(.x, envir = envir))||is.matrix(get(.x, envir = envir)))))
    out
    }



if (input=="<No active dataset>") input<-""

    output=""
    if(input!="") output<-paste(input, "res", sep=".")   
   filenamevar <- tclVar(input)
    tabnamevar <- tclVar(output)
    fictrt <- tclVar()
    arValue <- tclVar(arMAX)

    tfARparent <- tktoplevel()
 #   size=c(240,120,0,132)
 #  geo <- paste(size[1], "x", size[2], "+", size[3],"+", size[4], sep = "")
 #   tkwm.geometry(tfAR, geo)
    tkwm.geometry(tfARparent, paste("+0+",.heigth, sep=""))
    
    tkwm.title(tfARparent, "AR model")   
    tkwm.resizable(tfARparent, 0, 0)  
    tkwm.deiconify(tfARparent)
    tkgrab.set(tfARparent)
    tkfocus(tfARparent)
    
    choose.data = function() {
        input <- tk_select.list(sort(listDataSets ()), title="Select one")
        output <- paste(input, "res", sep=".")
        tkgrab.set(tfARparent)
        if(input!=""){
        tkdelete(file.entry, 0, "end")
        tkinsert(file.entry, "end", input)
        tkdelete(tab.entry, 0, "end")
        tkinsert(tab.entry, "end",output)
        }
        }

 
    done <- tclVar(0)
    tfAR <- tkframe(tfARparent, relief = "groove")
    tkpack(tfAR, side="top")
    frame1.a <- tkframe(tfAR, relief = "groove")
    frame1 <- tkframe(tfAR, relief = "groove")
   
    tkgrid(tklabel(frame1.a, text = "Options:", foreground = "blue"))
    tkpack(frame1.a, fill = "x")
      
    tab.entry <- tkentry(frame1, textvariable = tabnamevar)
    file.entry <- tkentry(frame1, textvariable = filenamevar)
    choosefile.but <- tkbutton(frame1, text = "...", command = function() choose.data())
    tkgrid(tklabel(frame1, text = "Input name: "), file.entry, tklabel(frame1, text = " "), choosefile.but,  sticky = "w")
    tkgrid(tklabel(frame1, text = "Output name:"), tab.entry,  sticky = "w")
    tkpack(frame1, fill = "x")
       
    frame3 <- tkframe(tfAR, relief = "groove")
    frame3.1 <- tkframe(frame3)
   
    makeAr.value <- tclVar(makeAr)
    arMAX.value <- tclVar(arMAX)
    arMAX.lab<- tklabel(frame3.1, text = "AR model of max order:", foreground="blue")
     
    slider <- tkscale(frame3.1, from=1, to=10, showvalue=T, variable=arValue,
                   resolution=1, orient="horizontal")
tkgrid(arMAX.lab, slider)
     tkpack(frame3.1, fill = "x")
    tkpack(frame3, fill = "x")
       
     frame.exit<-  tkframe(tfAR, relief = "groove")

 OnOk = function(){
       try(exists(tclvalue(filenamevar)),silent=T)->flag
if ( flag==TRUE){
eval(parse(text=paste("ArFunction(", tclvalue( filenamevar),", order.max=",as.numeric(tclvalue(arValue )), ")",sep="")))
eval(parse(text=paste(tclvalue(tabnamevar),"<<-apply(",tclvalue( filenamevar), ", 2,ar.func, order.max=",as.numeric(tclvalue(arValue )),")",sep="")))
tclvalue(done) <- 2 
}
.assign("arMAX",as.numeric(tclvalue(arValue )))

}


    fr.exit.space         <-   tklabel(frame.exit, text = " ")
    ok.but <- tkbutton(frame.exit, text =     "      Ok      ", command = function() OnOk())
    cancel.but <- tkbutton(frame.exit, text = "    Cancel    ", command = function() tkdestroy(tfARparent))
    tkgrid( cancel.but,fr.exit.space, ok.but)
    tkpack( frame.exit, side="right")
    
   # tkpack(frame2, fill="x")
    tkbind(tfAR, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(tfAR, "<KeyPress-Return>", function() OnOk())
    tkbind(tfAR, "<KeyPress-Escape>", function() tkdestroy(tfAR))
        tkfocus(tfAR)
    tkwait.variable(done)
    if (tclvalue(done) == "2")
    tkgrab.release(tfAR)
    tkdestroy(tfARparent)
    }

 # arMODEL()
