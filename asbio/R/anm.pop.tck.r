#-------------Geometric--------------#

anm.geo.growth.tck<-function(){         

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
        tkwm.title(tt,"Geometric population growth")       
        n0.entry <- tkentry(tt, textvariable=N0, width = 10)              
        lambda.entry <- tkentry(tt, textvariable=Lambda, width = 10) 
        time.entry<-tkentry(tt, textvariable=Time, width = 10)
        int.entry<-tkentry(tt, textvariable=Int, width = 10)        

	done <- tclVar(0)
  
        reset <- function()
        {
            tclvalue(N0)<-""                                       
            tclvalue(Lambda)<-""
            tclvalue(Time)<-""
            tclvalue(Int)<-""
        }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)
     
        build <- function()
        {
            n0 <-tclvalue(N0)
            lambda <-tclvalue(Lambda)
            time <-parse(text=tclvalue(Time))[[1]]
            interval <-tclvalue(Int) 
            
           substitute(anm.geo.growth(n0=as.numeric(n0),lambda=as.numeric(lambda),time=time, interval=as.numeric(interval)))
        }
        
        tkgrid(tklabel(tt,text="Geometric growth"),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="N\u2080"), n0.entry)
        tkgrid(tklabel(tt,text="\u03bb"), lambda.entry)
        tkgrid(tklabel(tt,text="Time seq."), time.entry)
        tkgrid(tklabel(tt,text="Anim. int. "),int.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but, reset.but, sticky="e")

        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
        invisible(tclServiceMode(TRUE))
    }                            
      N0<-tclVar("10")                          
      Lambda<-tclVar("1.3")
      Int<-tclVar("0.1")
      Time<-tclVar("seq(0,20)")
      dialog.sd()

})
}

#-------------Exponential--------------#

anm.exp.growth.tck<-function(){         

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
        tkwm.title(tt,"Exponential population growth")       
        n.entry <- tkentry(tt, textvariable=N, width = 10)              
        rmax.entry <- tkentry(tt, textvariable=Rmax, width = 10) 
        time.entry<-tkentry(tt, textvariable=Time, width = 10)
        int.entry<-tkentry(tt, textvariable=Int, width = 10)        

	done <- tclVar(0)
  
        reset <- function()
        {
            tclvalue(N)<-""                                       
            tclvalue(Rmax)<-""
            tclvalue(Time)<-""
            tclvalue(Int)<-""
        }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)
     
        build <- function()
        {
            n <-tclvalue(N)
            rmax <-tclvalue(Rmax)
            time <-parse(text=tclvalue(Time))[[1]]
            interval <-tclvalue(Int) 
            
           substitute(anm.exp.growth(n=as.numeric(n),rmax=as.numeric(rmax),time=time, interval=as.numeric(interval)))
        }
        
        tkgrid(tklabel(tt,text="Exponential growth"),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="N"), n.entry)
        tkgrid(tklabel(tt,text="r.max"), rmax.entry)
        tkgrid(tklabel(tt,text="Time seq."), time.entry)
        tkgrid(tklabel(tt,text="Anim. int. "),int.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but, reset.but, sticky="w")
  
        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
        invisible(tclServiceMode(TRUE))
    }                            
      N<-tclVar("10")                          
      Rmax<-tclVar("0.6")
      Int<-tclVar("0.1")
      Time<-tclVar("seq(0,20)")
      dialog.sd()
      
})
}


#-------------Logistic--------------#

anm.log.growth.tck<-function(){         

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
        tkwm.title(tt,"Logistic population growth")       
        n.entry <- tkentry(tt, textvariable=N, width = 10)              
        rmax.entry <- tkentry(tt, textvariable=Rmax, width = 10) 
        K.entry<-tkentry(tt, textvariable=Ke, width = 10)
        time.entry<-tkentry(tt, textvariable=Time, width = 10)
        int.entry<-tkentry(tt, textvariable=Int, width = 10)        

	done <- tclVar(0)
  
        reset <- function()
        {
            tclvalue(N)<-""                                       
            tclvalue(Rmax)<-""
            tclvalue(Ke)<-""
            tclvalue(Time)<-""
            tclvalue(Int)<-""
        }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)
     
        build <- function()
        {
            n <-tclvalue(N)
            rmax <-tclvalue(Rmax)
            time <-parse(text=tclvalue(Time))[[1]]
            K <-tclvalue(Ke)
            interval <-tclvalue(Int) 
            
           substitute(anm.log.growth(n=as.numeric(n),rmax=as.numeric(rmax), K=as.numeric(K), time=time, interval=as.numeric(interval)))
        }
        
        tkgrid(tklabel(tt,text="Logisitic growth"),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="N"), n.entry)
        tkgrid(tklabel(tt,text="r.max"), rmax.entry)
        tkgrid(tklabel(tt,text="K"), K.entry)
        tkgrid(tklabel(tt,text="Time seq."), time.entry)
        tkgrid(tklabel(tt,text="Anim. int."),int.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but, reset.but, sticky="w")
             

        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
        invisible(tclServiceMode(TRUE))
    }                            
      N<-tclVar("10")                          
      Rmax<-tclVar("0.6")
      Ke<-tclVar("1000")
      Int<-tclVar("0.1")
      Time<-tclVar("seq(0,50)")
      dialog.sd()
      
})
}
