plotBoth = function(plotfn, filename, control = plotBoth.control(), ...){
    if(control$genPlots){
        fname = paste('Figures/',filename,".eps",sep="")

        postscript(fname, horizontal = FALSE, width = 10, height = 8, paper = "special")
        plotfn(...)
        graphics.off()
        if(control$embedF)
            embedFonts(fname)

        fname = gsub('eps','pdf',fname)

        pdf(fname, width = 10, height = 8)
        plotfn(...)
        graphics.off()
        if(control$embedF)
            embedFonts(fname, options = control$embedFoptions)
    }
}
