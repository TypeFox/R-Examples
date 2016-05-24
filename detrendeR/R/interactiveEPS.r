interactiveEPS = function (input="")  {


listDataSets = function (envir = .GlobalEnv, ...){
    Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
    out=names(which(sapply(Vars, function(.x) is.data.frame(get(.x, envir = envir))||is.matrix(get(.x, envir = envir)))))
    out
    }
    
#bringToTop(which = -1)
if (input=="<No active dataset>") input<-""
    filenamevar <- tclVar(input)
    output=""
    if (input!="")  output<-paste(input, "EPS", sep=".")
    done<-tclVar(3)
    tabnamevar <- tclVar(output)
    fictrt <- tclVar()
    
    choose.data = function() {
        # input <- select.list(sort(listDataSets ()))
           input <- tk_select.list(sort(listDataSets ()), title="Select one")
        output <- paste(input, "EPS", sep=".")
        tkgrab.set(tf)
        if(input!=""){
        tkdelete(file.entry, 0, "end")
        tkinsert(file.entry, "end", input)
        tkdelete(tab.entry, 0, "end")
        tkinsert(tab.entry, "end",output)
          tkfocus(tf)
        }
        }

   
    tf <- tktoplevel()  
    tkwm.geometry(tf, paste("+0+",.heigth, sep=""))
    tkwm.title(tf, "EPS analysis")
    tkwm.resizable(tf,0,0)
      tkwm.deiconify(tf)
    tkgrab.set(tf)
    
    #       size=c(295,170,0,132)
    #geo <- paste(size[1], "x", size[2], "+", size[3],"+", size[4], sep = "")
    #tkwm.geometry(tf, geo)
    
    
    
    tkfocus(tf)
    done <- tclVar(0)
 
   
    frame1.a <- tkframe(tf, relief = "groove")
    frame1 <- tkframe(tf, relief = "groove")
    
    tkgrid(tklabel(frame1.a, text = "Options:", foreground = "blue"))
    tkpack(frame1.a, fill = "x")
      
    tab.entry <- tkentry(frame1, textvariable = tabnamevar)
    file.entry <- tkentry(frame1, textvariable = filenamevar)
    choosefile.but <- tkbutton(frame1, text = "...", command = function() choose.data())
    tkgrid(tklabel(frame1, text = "Input name: "), file.entry, tklabel(frame1, text = " "), choosefile.but,  sticky = "w")
    tkgrid(tklabel(frame1, text = "Output name:"), tab.entry,  sticky = "w")
    tkpack(frame1, fill = "x")
       
  
frame4 <- tkframe(tf, relief = "groove", borderwidth = 2)
frame4.1 <- tkframe(frame4, relief = "groove", borderwidth = 2)
frame4.1.1 <- tkframe( frame4.1, relief = "groove", borderwidth = 0)
frame4.1.2 <- tkframe( frame4.1, relief = "groove", borderwidth = 0)
frame4.1.3 <- tkframe( frame4.1, relief = "groove", borderwidth = 0)

frame4.2 <- tkframe(frame4, relief = "groove", borderwidth = 2)
frame4.2.1 <- tkframe( frame4.2, relief = "groove", borderwidth = 0)
frame4.2.2 <- tkframe( frame4.2, relief = "groove", borderwidth = 0)
frame4.2.3 <- tkframe( frame4.2, relief = "groove", borderwidth = 0)

run.win.analysis.value=tclVar(run.win.analysis)
winLength.value=tclVar(winLength)
stepWin.value=tclVar(stepWin)
make.common.EPS.value=tclVar(make.common.EPS)
make.select.period.EPS.value=tclVar(make.select.period.EPS)
first.year.common.value=tclVar(first.year.common)
last.year.common.value=tclVar(last.year.common)
 
frame4.label<-tklabel(frame4, text = " - EPS analysis - ", foreground = "blue")

run.win.analysis.cbut <- tkcheckbutton(frame4.1.1, text = "EPS window analysis", variable = run.win.analysis.value)
tkpack(frame4.label, fill = "x")
tkpack(run.win.analysis.cbut, fill = "x", side ="left")
tkpack(frame4.1.1, fill = "x")

winLength.entry <- tkentry(frame4.1.2, textvariable = winLength.value, width = 3)
winLength.lab<-tklabel(frame4.1.2, text = "Window length: ")
tkpack(winLength.lab,winLength.entry, side ="left") 
tkpack(frame4.1.2, fill = "x")

stepWin.entry <- tkentry(frame4.1.3, textvariable = stepWin.value, width = 3)
stepWin.lab<-tklabel(frame4.1.3, text = "                  Lag: ")
tkpack(stepWin.lab,stepWin.entry, side ="left") 
tkpack(frame4.1.3, fill = "x")

make.common.EPS.cbut <- tkcheckbutton(frame4.2.1, text = "Common interval", variable = make.common.EPS.value)

tkpack(make.common.EPS.cbut, fill = "x", side ="left")


   make.select.period.EPS.cbut <- tkcheckbutton(frame4.2.2, text = "Period        First year:", variable = make.select.period.EPS.value)
   first.year.common.entry <- tkentry(frame4.2.2, textvariable =  first.year.common.value, width = 4)
   tkpack(make.select.period.EPS.cbut,first.year.common.entry, fill = "x", side ="left")

last.year.common.entry <- tkentry(frame4.2.3, textvariable =  last.year.common.value, width = 4)
last.year.common.label<-tklabel(frame4.2.3, text = "                         Last year:")

tkpack(last.year.common.label, last.year.common.entry, side="left")

tkpack(frame4.2.1,frame4.2.2,frame4.2.3, fill = "x", side="top")

tkpack(frame4.1, frame4.2,fill = "x", side="left")
tkpack(frame4, fill = "x")
 
       frame.exit<-  tkframe(tf, relief = "groove")

 OnOk = function(){
               input<-tclvalue( filenamevar)
               output<-tclvalue(tabnamevar)
               winLength<<-tclvalue(winLength.value)
               run.win.analysis <<- as.logic(tclvalue(run.win.analysis.value))
               stepWin<<-tclvalue(stepWin.value)
                make.common.EPS<<-tclvalue(make.common.EPS.value)
                make.select.period.EPS<<- tclvalue(make.select.period.EPS.value)
               first.year.common<<-toNumber(tclvalue(first.year.common.value))
               last.year.common<<-toNumber(tclvalue(last.year.common.value))
       try(exists(tclvalue(filenamevar)),silent=T)->flag
      
if ( flag==TRUE){
    cat(paste("\nEPS analysis [",input,"]\n", sep=""))            
 if ( run.win.analysis ){   cat(rep("=",94),sep="")
 cat(paste( "\nRunning analysis:\n", sep=""), sep="")
       eval(parse(text=paste(output,"<<-Run.Win(", input,", winLength=",winLength, ", stc=stc,  step =", stepWin, ")",sep="")))
  eval(parse(text=paste("for (i in c(0.75, 0.80, 0.85, 0.90)) {EPS.resume(",output , ", EPS=i)}", sep="")))
          }
          
if (as.logic(tclvalue(make.common.EPS.value))) { eval(parse(text=paste("EPS.common.interval(",input ,", stc=stc, out=FALSE)"))) }

if (as.logic(tclvalue(make.select.period.EPS.value))) { if (!(is.null(first.year.common)||is.null(last.year.common.value))) {

        eval(parse(text=paste("EPS.common.interval(", input, ", first.year.common=", first.year.common, ",last.year.common=" ,last.year.common,  ", stc=stc, out=FALSE)"))) 
            } 
             }
tclvalue(done) <- 2 

}
}


    fr.exit.space         <-   tklabel(frame.exit, text = " ")
    ok.but <- tkbutton(frame.exit, text =     "      Ok      ", command = function() OnOk())
    cancel.but <- tkbutton(frame.exit, text = "    Cancel    ", command = function() tkdestroy(tf))
    tkgrid( cancel.but,fr.exit.space, ok.but)
    tkpack( frame.exit, side="right")
    tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(tf, "<KeyPress-Return>", function() OnOk())
    tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
    
    tkwait.variable(done)
    tkgrab.release(tf)
    if (tclvalue(done) == "2")
   
    tkdestroy(tf)
    #tkfocus(tt)
   # return(invisible(TRUE))
    }
    
    #interactiveEPS ()
