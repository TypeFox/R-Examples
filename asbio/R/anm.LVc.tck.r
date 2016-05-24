anm.LVc.tck<-function(){         

local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkcheckbutton <- ttkcheckbutton
        tkentry <- ttkentry
        tkframe <- ttkframe
        tklabel <- ttklabel
        tkradiobutton <- ttkradiobutton
    }
    tclServiceMode(FALSE)
    dialog.sd <- function(){
        tt <- tktoplevel()
        tkwm.title(tt,"Lotka Volterra competition")
        n1.entry <- tkentry(tt, textvariable=N1, width = 10)              
        n2.entry <- tkentry(tt, textvariable=N2, width = 10) 
        r1.entry <- tkentry(tt, textvariable=R1, width = 10) 
        r2.entry <- tkentry(tt, textvariable=R2, width = 10) 
        K1.entry <- tkentry(tt, textvariable=K1, width = 10)
        K2.entry<-tkentry(tt, textvariable=K2, width = 10)
        a2.1entry <- tkentry(tt, textvariable=A2.1, width = 10)
        a1.2entry <- tkentry(tt, textvariable=A1.2, width = 10)
        int.entry<-tkentry(tt, textvariable=Int, width = 10)
        Ts.entry<-tkentry(tt, textvariable=Ts, width = 10)
                
	done <- tclVar(0)
  
        reset <- function()
        {
            tclvalue(N1)<-"150"                                       
            tclvalue(N2)<-"50"
            tclvalue(R1)<-"0.7"
            tclvalue(R2)<-"0.8"
            tclvalue(K1)<-"750"
            tclvalue(K2)<-"1000"
            tclvalue(A2.1)<-"0.5"
            tclvalue(A1.2)<-"0.7"
            tclvalue(Int)<-"0.1"
            tclvalue(Ts)<-"seq(0,200)"
          }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)

        build <- function()
        {
            n1 <-tclvalue(N1)
            n2 <-tclvalue(N2)
            r1 <-tclvalue(R1)
            r2 <-tclvalue(R2)
            K1 <-tclvalue(K1)
            K2 <-tclvalue(K2)
            a2.1<-tclvalue(A2.1)
            a1.2<-tclvalue(A1.2)
            interval <-tclvalue(Int) 
            time<-parse(text=tclvalue(Ts))[[1]]    
           substitute(anm.LVcomp(n1=as.numeric(n1),n2=as.numeric(n2),r1=as.numeric(r1),r2=as.numeric(r2),K1=as.numeric(K1),K2=as.numeric(K2),a1.2=as.numeric(a1.2),a2.1=as.numeric(a2.1),interval=as.numeric(interval),time=time))
        }
        f0<-c("Calibri","10")
        f1<-c("Calibri","11")
        tkgrid(tklabel(tt,text="Lotka Volterra competition",font=f0),columnspan=4)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="n\u2081",font=f1), n1.entry,tklabel(tt,text="n\u2082",font=f1), n2.entry)
        tkgrid(tklabel(tt,text="r\u2081",font=f1), r1.entry,tklabel(tt,text="r\u2082",font=f1), r2.entry)
        tkgrid(tklabel(tt,text="K\u2081",font=f1), K1.entry,tklabel(tt,text="K\u2082",font=f1), K2.entry)
        tkgrid(tklabel(tt,text="\u03b1\u2081\u2082",font=f1), a1.2entry,tklabel(tt,text="\u03b1\u2082\u2081",font=f1), a2.1entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Time seq.",font=f0), Ts.entry)
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Anim. int. ",font=f0),int.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but,tklabel(tt,text=""),tklabel(tt,text=""), reset.but, sticky = "w")
 
        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)
        tkwait.variable(done)
        if(tclvalue(done)=="2") stop("aborted")
        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
    invisible(tclServiceMode(TRUE))
    }                            
      N1<-tclVar("150")   
      N2<-tclVar("50")
      R1<-tclVar("0.7")
      R2<-tclVar("0.8")
      K1<-tclVar("750")
      K2<-tclVar("1000")
      A2.1<-tclVar("0.5")
      A1.2<-tclVar("0.7")
      Int<-tclVar("0.1")
      Ts<-tclVar("seq(0,200)")
      dialog.sd()
})
}
