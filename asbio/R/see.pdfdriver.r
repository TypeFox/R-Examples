see.pdfdriver<-function(pdf,show.cdf=TRUE){
if(pdf=="beta"&show.cdf==TRUE)see.betacdf.tck()
if(pdf=="beta"&show.cdf==FALSE)see.beta.tck()
if(pdf=="binomial"&show.cdf==TRUE)see.bincdf.tck()
if(pdf=="binomial"&show.cdf==FALSE)see.bin.tck()
if(pdf=="chisq"&show.cdf==TRUE)see.chicdf.tck()
if(pdf=="chisq"&show.cdf==FALSE)see.chi.tck()
if(pdf=="F"&show.cdf==TRUE)see.Fcdf.tck()
if(pdf=="F"&show.cdf==FALSE)see.F.tck()
if(pdf=="gamma"&show.cdf==TRUE)see.gamcdf.tck()
if(pdf=="gamma"&show.cdf==FALSE)see.gam.tck()
if(pdf=="geometric"&show.cdf==TRUE)see.geocdf.tck()
if(pdf=="geometric"&show.cdf==FALSE)see.geo.tck()
if(pdf=="hypergeometric"&show.cdf==TRUE)see.hypercdf.tck()
if(pdf=="hypergeometric"&show.cdf==FALSE)see.hyper.tck()
if(pdf=="exponential"&show.cdf==TRUE)see.expcdf.tck()
if(pdf=="exponential"&show.cdf==FALSE)see.exp.tck()
if(pdf=="logistic"&show.cdf==TRUE)see.logiscdf.tck()
if(pdf=="logistic"&show.cdf==FALSE)see.logis.tck()
if(pdf=="lognormal"&show.cdf==TRUE)see.lnormcdf.tck()
if(pdf=="lognormal"&show.cdf==FALSE)see.lnorm.tck()
if(pdf=="negative binomial"&show.cdf==TRUE)see.nbincdf.tck()
if(pdf=="negative binomial"&show.cdf==FALSE)see.nbin.tck()
if(pdf=="normal"&show.cdf==TRUE)see.normcdf.tck()
if(pdf=="normal"&show.cdf==FALSE)see.norm.tck()
if(pdf=="poisson"&show.cdf==TRUE)see.poiscdf.tck()
if(pdf=="poisson"&show.cdf==FALSE)see.pois.tck()
if(pdf=="t"&show.cdf==TRUE)see.tcdf.tck()
if(pdf=="t"&show.cdf==FALSE)see.t.tck()
if(pdf=="uniform"&show.cdf==TRUE)see.unifcdf.tck()
if(pdf=="uniform"&show.cdf==FALSE)see.unif.tck()
if(pdf=="weibull"&show.cdf==TRUE)see.weibcdf.tck()
if(pdf=="weibull"&show.cdf==FALSE)see.weib.tck()
}



see.pdfdriver.tck<-function(){

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
        tkwm.title(tt,"Depiction of pdfs")
        Pdf<-tclVar("normal")  
        done <- tclVar(0)
        show.cdf<-tclVar(1)
    
    
    reset <- function()
        {
            tclvalue(Pdf)<-"normal"
            tclvalue(show.cdf)<-"1"
        }
        
        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)
        tkgrid(tklabel(tt,text="Choose pdf"),columnspan=2)  
        tkgrid(tklabel(tt,text=""),columnspan=2)
        cdf.cbut <- tkcheckbutton(tt, text="Show cdf", variable=show.cdf)
        tkgrid(tklabel(tt, text="pdf"))
        pdfs<- c("beta","binomial","chisq","exponential","F","gamma","geometric", "hypergeometric","logistic","lognormal","negative binomial","normal","poisson","t","uniform","weibull")
            comboBox <- tkwidget(tt,"ComboBox", editable=FALSE, values=pdfs, textvariable = Pdf, width = 17)
        
        
        
        build<-function(){
        f<-as.logical(tclObj(show.cdf))
        pdf<-tclvalue(Pdf)
        substitute(see.pdfdriver(pdf,f))
        }
          
        
        tkgrid(comboBox,columnspan=2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(cdf.cbut)
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
 
