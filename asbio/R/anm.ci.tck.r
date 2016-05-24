anm.ci.tck<-function(){

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
    dialog.ci <- function(){
        tt <- tktoplevel()
        tkwm.title(tt,"Confidence intervals")
        parent.entry <- tkentry(tt, textvariable=Parent, width =10)
        para.entry<-tkentry(tt, textvariable=True.val, width = 10)
        conf.entry <- tkentry(tt, textvariable=Conf, width = 10)
        sigma.entry <- tkentry(tt, textvariable=Sigma, width = 10)
        n.est.entry<-tkentry(tt, textvariable=N.est, width = 10)
        int.entry<-tkentry(tt, textvariable=Int, width = 10)
	done <- tclVar(0)
  Par<- tclVar("mu")

        reset <- function()
        {
            tclvalue(Parent)<-""
            tclvalue(True.val)<-"0"
            tclvalue(Par)<-"mu"
            tclvalue(Conf)<-".95"
            tclvalue(Sigma)<-"1"
            tclvalue(N.est)<-"100"
            tclvalue(Int)<-"0.1"
        }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",
                               command=function()tclvalue(done)<-0)

        build <- function()
        {
            parent  <- parse(text=tclvalue(Parent))[[1]]
            par.type  <- tclvalue(Par)
            par.val  <- tclvalue(True.val)
            sigma <- tclvalue(Sigma)
            conf<-tclvalue(Conf)
            n.est<-tclvalue(N.est)
            interval<-tclvalue(Int)
            
           substitute(anm.ci(parent,par.type=par.type,par.val=as.numeric(par.val),conf=as.numeric(conf),n.est=as.numeric(n.est),sigma=as.numeric(sigma), interval = as.numeric(interval)))
        }
        alt.rbuts <- tkframe(tt)

        tkpack(tklabel(alt.rbuts, text="Parameter"))
        for ( i in c("mu", "median","sigma.sq", "p")){
            tmp <- tkradiobutton(alt.rbuts, text=i, variable=Par, value=i)
            tkpack(tmp,anchor="w")
        }
        
        tkgrid(tklabel(tt,text="Confidence intervals"),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="Parent"), parent.entry)
        tkgrid(tklabel(tt,text="True value"), para.entry)
        tkgrid(tklabel(tt,text="Conf"), conf.entry)
        tkgrid(tklabel(tt,text='\u03c3',font=c("Helvetica","9","italic")), sigma.entry)
        tkgrid(tklabel(tt,text="Iterations"),n.est.entry)
        tkgrid(tklabel(tt,text="Animation interval"),int.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(alt.rbuts)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but,reset.but, sticky="w")

        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")
      
        tkdestroy(tt) 
        
        cmd <- build()
        eval.parent(cmd)
        invisible(tclServiceMode(TRUE))
    }                            
      Parent<-tclVar("rnorm(1000)")
      True.val<-tclVar("0")
      Conf<-tclVar(".95")
      N.est<-tclVar("100")
      Sigma<-tclVar("1")
      Int<-tclVar("0.1")
      dialog.ci()
   
})
}
