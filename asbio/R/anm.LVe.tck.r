anm.LVe.tck<-function(){         

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
        tkwm.title(tt,"Lotka Volterra exploitation")       
        nh.entry <- tkentry(tt, textvariable=Nh, width = 10)              
        np.entry <- tkentry(tt, textvariable=Np, width = 10) 
        rh.entry <- tkentry(tt, textvariable=Rh, width= 10) 
        con.entry <- tkentry(tt, textvariable=Con, width = 10 ) 
        p.entry <- tkentry(tt, textvariable=P, width = 10)
        dp.entry<-tkentry(tt, textvariable=Dp, width = 10)
        int.entry<-tkentry(tt, textvariable=Int, width = 10)
        Ts.entry<-tkentry(tt, textvariable=Ts, width = 10)
                
	done <- tclVar(0)
  show.circle<-tclVar(0)
        reset <- function()
        {
            tclvalue(Nh)<-"300"                                       
            tclvalue(Np)<-"50"
            tclvalue(Rh)<-"0.7"
            tclvalue(Con)<-"0.4"
            tclvalue(P)<-"0.006"
            tclvalue(Dp)<-"0.2"
            tclvalue(Int)<-"0.1"
            tclvalue(Ts)<-"seq(0,200)"
            tclvalue(show.circle)<-"0"
          }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)
        c.cbut <- tkcheckbutton(tt, text="Circle", variable=show.circle)
        build <- function()
        {
            nh <-tclvalue(Nh)
            np <-tclvalue(Np)
            rh <-tclvalue(Rh)
            con <-tclvalue(Con)
            p <-tclvalue(P)
            d.p <-tclvalue(Dp)
            interval <-tclvalue(Int) 
            time<-parse(text=tclvalue(Ts))[[1]]
            circle<-as.logical(tclObj(show.circle))    
           substitute(anm.LVexp(nh=as.numeric(nh),np=as.numeric(np),rh=as.numeric(rh),con=as.numeric(con),p=as.numeric(p),d.p=as.numeric(d.p),interval=as.numeric(interval),time=time,circle=circle))
        }
        
        tkgrid(tklabel(tt,text="Lotka Volterra exploitation"),columnspan=4)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""), tklabel(tt,text="Prey"),tklabel(tt,text=""),tklabel(tt,text="Predator"))
        tkgrid(tklabel(tt,text="    n.h", width = 5), nh.entry,tklabel(tt,text="n.p"), np.entry)
        tkgrid(tklabel(tt,text="    r.h", width = 5), rh.entry,tklabel(tt,text="c"), con.entry)
        tkgrid(tklabel(tt,text="    p", width = 5), p.entry,tklabel(tt,text="d.p"), dp.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Time seq."), Ts.entry)
        tkgrid(tklabel(tt,text=""),tklabel(tt,text="Anim. int."),int.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(c.cbut)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but, tklabel(tt,text=""),tklabel(tt,text=""),reset.but, sticky="e")
      
        

        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
    invisible(tclServiceMode(TRUE))
    }                            
      Nh<-tclVar("300")                          
      Np<-tclVar("50")
      Rh<-tclVar("0.7")
      Con<-tclVar("0.4")
      P<-tclVar("0.006")
      Dp<-tclVar("0.2")
      Int<-tclVar("0.1")
      Ts<-tclVar("seq(0,200)")
      dialog.sd()
})
}
