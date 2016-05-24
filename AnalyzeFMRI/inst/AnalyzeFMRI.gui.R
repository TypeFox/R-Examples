require(tcltk) || stop("tcltk support is absent")

local({
 
    wrap.file <- function() tclvalue(file.name) <- tcl("tk_getOpenFile")
    
    wrap.mask <- function() tclvalue(mask) <- tcl("tk_getOpenFile") 
    
    do <- function(){
        if(tclvalue(alt) == "File Summary") fs()
        if(tclvalue(alt) == "Plot Time Series") pts()
        if(tclvalue(alt) == "Plot Periodogram") period()
        if(tclvalue(alt) == "Image Slice") im.sl()
        if(tclvalue(alt) == "Image Volume") im.vol() 
        if(tclvalue(alt) == "Movie") im.mov() 
        if(tclvalue(alt) == "Spectral Summary") im.spec()
    }
    
    fs <- function(...) f.analyze.file.summary(tclvalue(file.name))
    
    pts <- function(...) {
        plot(f.read.analyze.ts(tclvalue(file.name),
                               as.numeric(tclvalue(x)),
                               as.numeric(tclvalue(y)),
                               as.numeric(tclvalue(z))),
             typ = "l", ylab = "fMRI response",
             xlab = "Scans")
    }

    period <- function(...){
        par(mfrow = c(1, 1), mar = c(4, 4, 5, 5))
        a <- f.read.analyze.ts(tclvalue(file.name),
                               as.numeric(tclvalue(x)),
                               as.numeric(tclvalue(y)),
                               as.numeric(tclvalue(z)))
        b <- fft(a) / sqrt(2 * pi * length(a))
        b <- b[10:floor(length(b) / 2) + 1]
        b <- Mod(b)^2
        plot(b, ylab = "Periodogram", xlab = "Fourier Frequency")
    }
    
    im.sl <- function(...){
        par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
        a <- f.read.analyze.slice(tclvalue(file.name),
                                  as.numeric(tclvalue(z)),
                                  as.numeric(tclvalue(t)))
        
        image(a)
        par(mfrow = c(1, 1), mar = c(4, 4, 5, 5))
    }
    
    im.vol <- function(...){
        a <- f.read.analyze.header(tclvalue(file.name))$dim
        d <- ceiling(sqrt(a[4]))
        par(mfrow =c(d, d), mar = c(0, 0, 0, 0))
        b <- array(0, dim = a[2:4])
        for(i in 1:a[4]){
            b[, , i] <- f.read.analyze.slice(tclvalue(file.name),
                                             i,
                                             as.numeric(tclvalue(t)))
        }
        for(i in 1:a[4]){
            image(b[, , i], axes = FALSE)
            box()
        }
        
        par(mfrow = c(1, 1), mar = c(4, 4, 5, 5))
    }
    
    im.mov <- function(...){
        par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
        a <- f.read.analyze.header(tclvalue(file.name))$dim
        b <- array(0, dim = c(a[2], a[3], a[5]))
        for(i in 1:a[5]){
            b[, , i] <- f.read.analyze.slice(tclvalue(file.name),
                                             as.numeric(tclvalue(z)),
                                             i)
        }
        image(b[, , 1], axes = FALSE)
        for(i in 2:a[5]){
            image(b[, , i], axes = FALSE, add = TRUE)
        }
        par(mfrow = c(1, 1), mar = c(4, 4, 5, 5))
    }
    
    im.spec <- function(...){
        par(mfrow = c(1, 1), mar = c(4, 4, 5, 5))
        if(tclvalue(mask) == "") tclvalue(mask) <- FALSE
        a <- f.spectral.summary(tclvalue(file.name), tclvalue(mask))
        par(mfrow = c(1, 1), mar = c(4, 4, 5, 5))
    }
    
    ## set up tclVar variables

    mask <- tclVar("sdasd")
    alt <- tclVar()
    x <- tclVar()
    y <- tclVar()
    z <- tclVar()
    t <- tclVar()
    
    #set up base GUI window
    if(.Platform$OS.type == "windows") flush.console()
    
    base <- tktoplevel()
    tkwm.title(base, "ANALYZE file explore")
    f1 <- tkframe(base, relief = "groove", borderwidth = 2)

    file.name <- tclVar()
    tkpack(tkentry(f1, textvariable = file.name, width = 40))
    
    file.find.but <- tkbutton(f1, text = "Select File", command = wrap.file)
    tkpack(file.find.but)
    
    mask <- tclVar()
    tkpack(tkentry(f1, textvariable = mask, width = 40))
    mask.find.but <- tkbutton(f1, text = "Select Mask File", command = wrap.mask)
    tkpack(mask.find.but)
    
    opt.rbuts <- tkframe(base, relief = "groove", borderwidth = 2)
    
    tkpack(tklabel(opt.rbuts, text = "Options"))

    alt <- tclVar()
    for (i in c("File Summary",
                "Plot Time Series",
                "Plot Periodogram",
                "Image Slice",
                "Image Volume",
                "Movie",
                "Spectral Summary")) {
        tmp <- tkradiobutton(opt.rbuts, text = i, variable = alt, value = i)
        tkpack(tmp, anchor = "w")
    }
    fr2 <- tkframe(base, relief = "groove", borderwidth = 2)
    x <- tclVar()
    x.entry <- tkentry(fr2, textvariable = x)
    y <- tclVar()
    y.entry <- tkentry(fr2, textvariable = y)
    z <- tclVar()
    z.entry <- tkentry(fr2, textvariable = z)
    t <- tclVar()
    t.entry <- tkentry(fr2, textvariable = t)
    tkgrid(f1)
    
    tkgrid(tklabel(fr2, text = "Variables"), columnspan = 2)
    tkgrid(tklabel(fr2, text = "x variable"), x.entry)
    tkgrid(tklabel(fr2, text = "y variable"), y.entry)
    tkgrid(tklabel(fr2, text = "z variable"), z.entry)
    tkgrid(tklabel(fr2, text = "t variable"), t.entry)
    tkgrid(opt.rbuts)
    tkgrid(fr2)
    
    fr3 <- tkframe(base, borderwidth = 2)
    q.but <- tkbutton(fr3, text = "Quit",
                      command = function() tkdestroy(base))
    ok.but <- tkbutton(fr3, text = "OK", command = do)
    tkgrid(ok.but, q.but)
    tkgrid(fr3)
    
})


