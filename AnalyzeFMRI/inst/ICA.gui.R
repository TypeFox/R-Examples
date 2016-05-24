
require(tcltk) || stop("tcltk support is absent")

local({


    gui.file<-function(){
        tclvalue(file.name) <- tcl("tk_getOpenFile")   }
    
    
    gui.mask<-function(){
        tclvalue(mask) <- tcl("tk_getOpenFile")  }
    
    gui.save<-function(...){
        assign(tclvalue(save.r), tmp.ica.obj, envir = .GlobalEnv)
    }
    
    
    
    
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
        
        tmp.ica.obj <<- f.ica.fmri(tclvalue(file.name),
                                   n.comp = as.numeric(tclvalue(n.comp)),
                                   norm.col = v.norm,
                                   fun = "logcosh",
                                   maxit = 100,
                                   alg.type = "parallel",
                                   alpha = 1,
                                   tol = 0.0001,
                                   mask.file.name = msk,
                                   sl)
        
        print("done")
        
    }
    
    gui.jpeg<-function(...){  
        f.plot.ica.fmri.jpg(tmp.ica.obj, tclvalue(jpeg), width = 700, height = 700)}
    
    gui.end<-function(...){
        rm(tmp.ica.obj, envir = .GlobalEnv)
        tkdestroy(base.ica)
    }
    gui.plot.ica<-function(...){  
        f.plot.ica.fmri(tmp.ica.obj, as.numeric(tclvalue(comp)))}
    
    ##set tcl variables
    file.name <- tclVar()
    mask <- tclVar()
    n.comp <- tclVar("30")
    var.norm <- tclVar("0")
    slices <- tclVar("0")
    create.mask <- tclVar("0")
    save.r <- tclVar()
    jpeg <- tclVar()
    comp <- tclVar()
    
    #set up base GUI window
    if(.Platform$OS.type == "windows") flush.console()
    
    base.ica <- tktoplevel(bg="#555555")
    tkwm.title(base.ica, "Spatial ICA for fMRI datasets")
    
    
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
    
    ica.n.comp.label <- tklabel(ica.f2, text = "Number of components to extract", bg = "#aaaaaa")      
    ica.n.comp.entry <- tkentry(ica.f2, textvariable = n.comp, width = 5, bg = "#ffffff")
    tkgrid(ica.n.comp.label, ica.n.comp.entry, padx = 10, pady = 10)
    
    tkgrid(ica.f2, sticky = "ew")
    
    
    
    #frame for options
    ica.f3 <- tkframe(base.ica, relief = "groove", borderwidth = 2, bg = "#555555")
    
    ica.normalise.but <- tkcheckbutton(ica.f3, text = "Variance Normalize", bg = "#aaaaaa", variable = var.norm) 
    ica.slices.but <- tkcheckbutton(ica.f3, text = "Exclude top/bottom slices", bg = "#aaaaaa", variable = slices) 
    ica.create.mask.but <- tkcheckbutton(ica.f3, text = "Create Mask", bg = "#aaaaaa", variable = create.mask)
    tkgrid(ica.normalise.but, ica.slices.but, ica.create.mask.but, padx = 30, pady = 10)
    
    tkgrid(ica.f3, sticky = "ew")
    
    #frame for saving object to R session
    ica.f4 <- tkframe(base.ica, relief = "groove", borderwidth = 2, bg = "#555555")
    
    ica.save.entry <- tkentry(ica.f4, textvariable = save.r, width = 40, bg = "#ffffff")
    ica.save.but <- tkbutton(ica.f4, text = "Save to R object", width = 15, command = gui.save, bg = "#aaaaaa")
    tkgrid(ica.save.but, ica.save.entry, padx = 10, pady = 10)
    
    tkgrid(ica.f4,sticky = "ew")
    
    #frame for plotting components to jpeg files
    ica.f5 <- tkframe(base.ica, relief = "groove", borderwidth = 2, bg = "#555555")
    
    ica.jpeg.entry <- tkentry(ica.f5, textvariable = jpeg, width = 40, bg = "#ffffff")
    ica.jpeg.but <- tkbutton(ica.f5, text = "Save to jpeg files", width = 15, command = gui.jpeg, bg = "#aaaaaa")
    tkgrid(ica.jpeg.but, ica.jpeg.entry, padx = 10, pady = 10)
    
    tkgrid(ica.f5, sticky = "ew")
    
    
    #frame for plotting components
    ica.f6 <- tkframe(base.ica, relief = "groove", borderwidth = 2, bg = "#555555")
    
    ica.plot.entry <- tkentry(ica.f6, textvariable = comp, width = 5, bg = "#ffffff")
    ica.plot.but <- tkbutton(ica.f6, text = "Plot component", width = 15, command = gui.plot.ica, bg = "#aaaaaa")
    tkgrid(ica.plot.but, ica.plot.entry, padx = 10, pady = 10)
    
    tkgrid(ica.f6, sticky = "ew")
    
    #frame for start and end buttons     
    fr3 <- tkframe(base.ica, borderwidth = 2, bg = "#555555")
    go.but<- tkbutton(fr3, text = "Start", bg = "#aaaaaa", command = gui.ica)
    q.but <- tkbutton(fr3, text = "Quit",
                      command = gui.end, bg = "#aaaaaa")
    tkgrid(go.but, q.but, padx = 30, pady = 20)
    tkgrid(fr3)
    
})
