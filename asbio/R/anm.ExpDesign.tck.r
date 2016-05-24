anm.ExpDesign.tck<-function(){

tclRequire("BWidget")
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
        tkwm.title(tt,"Experimental designs")
        Proc<-tclVar("all")
        int.entry <- tkentry(tt, textvariable=Int, width = 10)
        iter.entry<-tkentry(tt, textvariable=Iter, width = 10)
       
                
	done <- tclVar(0)
  
        reset <- function()
        {
          
            tclvalue(Int)<-"0.5"
            tclvalue(Iter)<-"30"
          }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)

        build <- function()
        {
       
            interval <-tclvalue(Int) 
            iter <-tclvalue(Iter)
            Proc <-tclvalue(Proc)
           substitute(anm.ExpDesign(method = Proc, interval=as.numeric(interval),iter=as.numeric(iter)))
        }
        
        tkgrid(tklabel(tt,text="Experimental designs"),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt, text="Design"),columnspan=2)
        proc <- c("all","CRD","factorial2by2","factorial2by2by2","nested","RCBD","RIBD","split","split.split", "SPRB","strip","split.block","strip.split","latin","pairs")
        comboBox <- tkwidget(tt,"ComboBox", editable=FALSE, values=proc, textvariable = Proc, width = 15)
        tkgrid(comboBox,columnspan=2)
        
        tkgrid(tklabel(tt,text="Anim. int."),int.entry)
        tkgrid(tklabel(tt,text="Iterations "), iter.entry)
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

      Iter<-tclVar("30")
      Int<-tclVar("0.5")
      dialog.sd()
})
}
