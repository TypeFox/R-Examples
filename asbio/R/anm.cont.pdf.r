anm.cont.pdf<-function(part ="norm",interval=.3){
if(part == "norm") p = rnorm(1000000)
if(part == "t") p = rt(1000000,df=10)
if(part == "exp") p = rexp(1000000)
if(part == "unif") p = runif(1000000, 0, 1)

if(part == "norm") {xlim=c(-6,6);ylim=c(0,0.5)}
if(part == "t") {xlim=c(-15,15);ylim=c(0,0.4)}
if(part =="exp") {xlim=c(0,20);ylim=c(0,0.9)}
if(part =="unif") {xlim=c(0,1);ylim=c(0,1.2)}

for(i in 1:20){
hist(p,breaks=seq(xlim[1],xlim[2],by=1/i),freq=F, ylim = ylim, xlim = xlim, main = "", xlab=expression(italic(x)), ylab="Density")
dev.flush()
Sys.sleep(interval)
}
if(part=="norm")
x <- NULL; rm(x); # Dummy to trick R CMD check 
curve(if(part=="norm")dnorm(x)
else if(part=="t") dt(x, df = 10)
else if(part=="exp") dexp(x)
else if(part=="unif") dunif(x)
,xlim[1],xlim[2],ylim=ylim, xlab="",yaxt="n",bty="n",col=1,lwd=ifelse(part == "t", 3, 7), add=TRUE) 
}



see.pdf.conc.tck<-function(){

see.pdf.conc<-function(pdf,show.cdf=TRUE){
if(pdf=="t")anm.cont.pdf("t")
if(pdf=="norm")anm.cont.pdf("norm")
if(pdf=="exp")anm.cont.pdf("exp")
if(pdf=="unif")anm.cont.pdf("unif")
}


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
        tkwm.title(tt,"Conceptualization of continuous pdfs")
        Pdf<-tclVar("normal")  
        done <- tclVar(0)
        
    
    
    reset <- function()
        {
            tclvalue(Pdf)<-"normal"
        }
        
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)
        tkgrid(tklabel(tt,text="Conceptualize pdfs"),columnspan=2)  
        tkgrid(tklabel(tt,text=""),columnspan=2)
        tkgrid(tklabel(tt, text="pdf"))
        pdfs<- c("exponential","normal","t","uniform")
            comboBox <- tkwidget(tt,"ComboBox", editable=FALSE, values=pdfs, textvariable = Pdf, width = 17)
        
        
        
        build<-function(){
        pdf<-tclvalue(Pdf)
        if(pdf =="exponential")dist = "exp"
        if(pdf == "normal") dist = "norm"
        if(pdf == "t") dist = "t"
        if(pdf == "uniform") dist = "unif"
        substitute(see.pdf.conc(dist))
        }
          
        
        tkgrid(comboBox,columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text=""))
        tkgrid(submit.but, reset.but, sticky="w")
        
        
        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)

        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")

        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
    }                            
      Pdf<- tclVar("normal")
      dialog.sd()
})
}
 
