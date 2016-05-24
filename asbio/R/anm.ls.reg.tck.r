anm.ls.reg.tck<-function(){

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
        tkwm.title(tt,'OLS estimation of regression parameters')
        X.entry <- tkentry(tt, textvariable=X, width = 10)
        Y.entry <- tkentry(tt, textvariable=Y, width = 10)
        Par.entry <- tkentry(tt, textvariable=Par, width = 10)
        int.entry <- tkentry(tt, textvariable=Int, width = 10)
              
                
	done <- tclVar(0)
  
        reset <- function()
        {
            tclvalue(X)<-"rnorm(30)"
            tclvalue(Y)<-"rnorm(30)"
            tclvalue(Y)<-"Par"
            tclvalue(Int)<-"0.1"
           
          }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)

        build <- function()
        {
            X <-n.seq<-parse(text=tclvalue(X))[[1]]
            Y <-n.seq<-parse(text=tclvalue(Y))[[1]]
            interval <-tclvalue(Int) 
            par<-tclvalue(Par)
     
           substitute(anm.ls.reg(X=as.numeric(X),Y=as.numeric(Y),parameter=par,interval=as.numeric(interval)))
        }
        
        tkgrid(tklabel(tt,text='OLS estimation of regression'),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="X"), X.entry)
        tkgrid(tklabel(tt,text="Y"), Y.entry)
        tkgrid(tklabel(tt,text="Parameter"), Par.entry)
        tkgrid(tklabel(tt,text="Anim. int."),int.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but, reset.but, sticky = "w")
         

        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
    invisible(tclServiceMode(TRUE))
    }                            
      X<-tclVar("rnorm(30)")
      Y<-tclVar("rnorm(30)")
      Int<-tclVar("0.1")
      Par<-tclVar("slope")
      dialog.sd()
})
}
