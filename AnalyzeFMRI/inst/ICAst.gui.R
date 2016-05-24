
require(tcltk) || stop("tcltk support is absent")

local({


    gui.file<-function(){
        tclvalue(file.name) <- tcl("tk_getOpenFile")   }
    
    
    gui.mask<-function(){
        tclvalue(mask) <- tcl("tk_getOpenFile")  }    
    
    
    
    gui.ica<-function(...){
        
        if(tclvalue(mask) == "" && tclvalue(create.mask) == 0){
            err1 <- tktoplevel()
            err1.f1 <- tkframe(err1, relief="groove", borderwidth = 2)
            err1.label <- tklabel(err1.f1, text = "Either select a mask file or select create mask", bg = "#aaaaaa", pady = 20, padx = 20)
            tkgrid(err1.label)
            tkgrid(err1.f1)
            return()
        }
        
        if(tclvalue(mask) != "" && tclvalue(create.mask) == 0) msk <- tclvalue(mask)
        if(tclvalue(mask) == "" && tclvalue(create.mask) == 1) msk <- NULL
        v.norm <- 1
        if(tclvalue(var.norm) == 0) v.norm <- 0
        sl <- NULL
        if(tclvalue(slices) == 0) sl <- "all"

        if (tclvalue(spatial.or.temporal) == "spatial") is.spatial <- 1 else is.spatial <- 0
        
        cat("Computations have begun. Please, be patient ...\n")

        tmp.ica.obj <<- f.icast.fmri(foncfile = tclvalue(file.name),
                               maskfile = msk,
                               is.spatial = is.spatial,
                               n.comp.compute = as.numeric(tclvalue(n.comp.compute)),
                               n.comp = as.numeric(tclvalue(n.comp)),
                               hp.filter=TRUE)


        
        cat("done\n")
        if (exists("tmp.ica.obj")) rm(tmp.ica.obj, envir = .GlobalEnv)
        gc(F)
        
    }
    
    
    gui.end<-function(...){
        if (exists("tmp.ica.obj")) rm(tmp.ica.obj, envir = .GlobalEnv)
        tkdestroy(base.ica)
    }




    
    
    ##set tcl variables
    file.name <- tclVar()
    mask <- tclVar()
    n.comp <- tclVar("")
    var.norm <- tclVar("0")
    n.comp.compute <- tclVar("1")
    slices <- tclVar("0")
    create.mask <- tclVar("0")
    save.r <- tclVar()
    jpeg <- tclVar()
    comp <- tclVar()
    spatial.or.temporal <- tclVar("spatial")
   
    #set up base GUI window
    if(.Platform$OS.type == "windows") flush.console()
    
    base.ica <- tktoplevel(bg="#555555")
    tkwm.title(base.ica, "Spatial or Temporal ICA for fMRI datasets")
    
    
    #frame to contain file selection
    ica.f1 <- tkframe(base.ica, relief = "groove", borderwidth = 2, bg = "#555555")
    
    ica.file.entry <- tkentry(ica.f1, textvariable = file.name, width = 50, bg = "#ffffff")
    ica.file.find.but <- tkbutton(ica.f1, text = "Select File", width = 15, command = gui.file, bg = "#aaaaaa", anchor = "c")
    tkgrid(ica.file.find.but, ica.file.entry, pady = 10, padx = 10)
    
    ica.mask.entry <- tkentry(ica.f1, textvariable = mask, width = 50, bg = "#ffffff")
    ica.mask.find.but <- tkbutton(ica.f1, text = "Select Mask File", width = 15, command = gui.mask, bg = "#aaaaaa")
    tkgrid(ica.mask.find.but, ica.mask.entry, padx = 10, pady = 10)

    tkgrid(ica.f1)
    
    #frame for number of components       
    ica.f2 <- tkframe(base.ica, relief = "groove", borderwidth = 2, bg = "#555555")
    
    ica.n.comp.radio <- tkradiobutton(ica.f2)
    ica.n.comp.label <- tklabel(ica.f2, text = "Enter number of components to extract", bg = "#aaaaaa")      
    ica.n.comp.entry <- tkentry(ica.f2, textvariable = n.comp, width = 5, bg = "#ffffff")
    ica.n.comp.compute <- tklabel(ica.f2, text = "Or automatic choice", bg = "#aaaaaa")
    ica.n.comp.compute.radio <- tkradiobutton(ica.f2)

    tkconfigure(ica.n.comp.radio,variable=n.comp.compute,value="0")
    tkconfigure(ica.n.comp.compute.radio,variable=n.comp.compute,value="1")
    
    
    tkgrid(ica.n.comp.radio,ica.n.comp.label, ica.n.comp.entry, ica.n.comp.compute,ica.n.comp.compute.radio, padx = 10, pady = 10)
    
    tkgrid(ica.f2, sticky = "ew")
     
       
    #frame for options
    ica.f3 <- tkframe(base.ica, relief = "groove", borderwidth = 2, bg = "#555555")
    
    ica.normalise.but <- tkradiobutton(ica.f3) 
    ica.slices.but <- tkradiobutton(ica.f3) 
    tkconfigure(ica.normalise.but,variable=spatial.or.temporal,value="spatial")
    tkconfigure(ica.slices.but,variable=spatial.or.temporal,value="temporal")

    
    tkgrid(tklabel(ica.f3,text="Spatial ICA ", bg = "#aaaaaa"),ica.normalise.but,tklabel(ica.f3,text="Temporal ICA ", bg = "#aaaaaa"), ica.slices.but, padx = 30, pady = 10)
    tkgrid(ica.f3, sticky = "ew")



  
    #frame for start and end buttons     
    fr3 <- tkframe(base.ica, borderwidth = 2, bg = "#555555")
    go.but<- tkbutton(fr3, text = "Start", bg = "#aaaaaa", command = gui.ica)
    q.but <- tkbutton(fr3, text = "Quit",
                      command = gui.end, bg = "#aaaaaa")
    tkgrid(go.but, q.but, padx = 30, pady = 20)
    tkgrid(fr3)
    
})
