shade.F.tck<-function(){

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
        tkwm.title(tt,"Depiction of F probability")
        x.entry <- tkentry(tt, textvariable=X, width = 10)
        nu1.entry<-tkentry(tt, textvariable=Nu1, width = 10)
        nu2.entry<-tkentry(tt, textvariable=Nu2, width = 10)
        from.entry<-tkentry(tt, textvariable=From, width = 10)
        to.entry<-tkentry(tt, textvariable=To, width = 10)
       
  Tail.par<-tclVar("lower")
  done <- tclVar(0)
 	show.p<-tclVar(1)
  show.d<-tclVar(0)
  show.dist<-tclVar(1)
        reset <- function()
        {
            tclvalue(X)<-"-.2"
            tclvalue(Nu1)<-"5"
            tclvalue(Nu2)<-"5"
            tclvalue(From)<-""
            tclvalue(To)<-""
            tclvalue(show.p)<-"1"
            tclvalue(show.d)<-"0"
            tclvalue(show.dist)<-"1"
         }
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)

        build <- function()
        {
            x <- tclvalue(X)
            nu1 <-tclvalue(Nu1)
            nu2 <-tclvalue(Nu2)
            from <-tclvalue(From)
            to<-tclvalue(To)
            tail<-tclvalue(Tail.par)
            show.p <- as.logical(tclObj(show.p))
            show.d <- as.logical(tclObj(show.d))
            show.dist <- as.logical(tclObj(show.dist))
                      
         substitute(shade.F(x=as.numeric(x),nu1=as.numeric(nu1), nu2=as.numeric(nu2),from = as.numeric(from), to=as.numeric(to),tail=tail,show.p=show.p,show.d=show.d,show.dist=show.dist))
        }
        
        p.cbut <- tkcheckbutton(tt, text="Show probability", variable=show.p)
        d.cbut <- tkcheckbutton(tt, text="Show density", variable=show.d)
        dist.cbut <- tkcheckbutton(tt, text="Show distribution", variable=show.dist)
        
               
        tkgrid(tklabel(tt,text="F probability"),columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="x",font=c("Helvetica","9","italic")), x.entry)
        tkgrid(tklabel(tt,text='\u03bd\u2081'), nu1.entry)
        tkgrid(tklabel(tt,text='\u03bd\u2082'), nu2.entry)
        tkgrid(tklabel(tt,text=""))
        alt.rbuts <- tkframe(tt)

        tkpack(tklabel(alt.rbuts, text="Tail"))
        for ( i in c("lower","upper","two","middle")){
            tmp <- tkradiobutton(alt.rbuts, text=i, variable=Tail.par, value=i)
            tkpack(tmp,anchor="w")
            }
        tkgrid(alt.rbuts)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="Middle 'tail' span"))
        tkgrid(tklabel(tt,text="From"),from.entry)
        tkgrid(tklabel(tt,text="To"),to.entry)
        tkgrid(tklabel(tt,text=""))
        tkgrid(p.cbut,sticky="w", columnspan = 2)
        tkgrid(d.cbut,sticky="w", columnspan = 2)
        tkgrid(dist.cbut,sticky="w", columnspan = 2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but,reset.but, sticky = "w")

       
        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
        tclServiceMode(TRUE)
    }                            
      X<-tclVar("2")
      Nu1<-tclVar("5")
      Nu2<-tclVar("5")
      Tail<-tclVar("lower")
      From<-tclVar("")
      To<-tclVar("")
      dialog.sd()
})
}
