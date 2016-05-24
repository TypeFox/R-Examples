anm.TM.tck<-function(){   



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
    dialog <- function(){
        tt <- tktoplevel()
        tkwm.title(tt,"Matrix population models")
        a.entry <- tkentry(tt, textvariable=Av, width = 20)
        n.entry <- tkentry(tt, textvariable=nv, width = 20)
        int.entry <- tkentry(tt, textvariable=Int, width = 10) 
	      aint.entry <- tkentry(tt, textvariable=Aint, width = 10) 
    done <- tclVar(0)
	
        reset <- function()
          {
            tclvalue(Av)<-""
            tclvalue(nv)<-""
            tclvalue(Int)<-""
            tclvalue(Aint)<-"0.1"
          }
        
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",
                               command=function()tclvalue(done)<-1)

        build <- function()
        {
            a  <- parse(text=tclvalue(Av))[[1]]
            n  <- parse(text=tclvalue(nv))[[1]]
            inter <- tclvalue(Int)
            anim.interval <- tclvalue(Aint)
            substitute(anm.transM(a,n,inter=as.numeric(inter),anim.interval=as.numeric(anim.interval),xlab="Time intervals from present"))
        }
               
        tkgrid(tklabel(tt,text="Matrix population models"),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="A",font=c("Helvetica","9","bold"),width=3), a.entry)
        tkgrid(tklabel(tt,text="n",font=c("Helvetica","9","bold"),width=3), n.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="Time ints. "), int.entry)
        tkgrid(tklabel(tt,text="Anim. int. "),aint.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but,reset.but, sticky ="e")

        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
    invisible(tclServiceMode(TRUE))
    }

    Av <- tclVar("matrix(nrow=3,ncol=3,data=c(.672,0,.561,0.018,0.849,0,0,0.138,0.969),byrow=TRUE)")
    nv <- tclVar("c(10,2,1)")
    Int <-tclVar("100")
    Aint <-tclVar("0.1") 
    dialog()
    })
    }
