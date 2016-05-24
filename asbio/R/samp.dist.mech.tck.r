samp.dist.mech.tck <-function(){

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
        tkwm.title(tt,"Sampling distribution basics")
        ns.entry <- tkentry(tt, textvariable=ns, width =6) 
        int.entry <- tkentry(tt, textvariable=int, width =6) 
        done <- tclVar(0)
       

   build <- function()
        {
            ns  <- tclvalue(ns)
            int  <- tclvalue(int)
                        
           substitute(samp.dist.mech(as.numeric(ns),as.numeric(int)))
        }

#desc <- tklabel(tt,text="A population of mountain goats are \r
#randomly sampled with sample size 10.\rThe goats are weighed, 
#and a mean is computed from the 10 weights. \r 
#This is repeated a user-specified number of times \r 
#to create a distribution of mean weights.")  
desc1 <- tklabel(tt, text = "Goat weight ~ N(90.5, 225)", font=c("Helvetica","12","normal"))
desc2 <- tklabel(tt, text = "", font=c("Helvetica","1","normal"))
desc3 <- tklabel(tt, text = "n = 10", font=c("Helvetica","12","normal"))

g1 <- read.pnm(system.file("pictures/goat3.pgm", package="asbio"))
f <- function(){
old.par <- par(no.readonly = TRUE)
par(mar = c(0,2,0,2))
plot(g1)
old.par <- par(no.readonly = TRUE)
}

submit.but <- tkbutton(tt, text="Submit", command=function()tclvalue(done)<-0)


#img <-tkrplot(tt, f, hscale=.8, vscale=.8)
tkpack(desc1, desc2, desc3)
tkpack(tklabel(tt,text="Iterations"), ns.entry)
tkpack(tklabel(tt,text="Animation interval"), int.entry)
tkpack(tklabel(tt,text="",font=c("Helvetica","2","normal")))
tkpack(submit.but)
tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")
      
        tkdestroy(tt) 
        
cmd <- build()
eval.parent(cmd)
invisible(tclServiceMode(TRUE))
}
ns<-tclVar("100")
int<-tclVar(".03")
dialog.ci()
})
}