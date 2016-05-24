.keystroke.separate <-
function (key)
{
    if (identical(key, "q")) 
        return(invisible(1))
    expectationwindow <- get("expectationwindow", envir = parent.frame())
    selectorwindow <- get("selectorwindow", envir = parent.frame())
    if (any(vapply(c("Up", "Down", "Left", "Right"), identical, logical(1), key)
        )){
        plotxsobject <- get("plotxsobject", envir = parent.frame())
        Xc <- get("Xc", envir = parent.frame())
        Xc.cond <- plotxsobject$xc.cond
        if (plotxsobject$view3d){
            dev.hold()
                try(plotxsobject$xc.cond <- get("Xc.cond", envir = 
                    parent.frame()), silent = TRUE)
                if (identical(key, "Down"))
                    plotxsobject$phi3d <- plotxsobject$phi3d + 2
                if (identical(key, "Up"))
                    plotxsobject$phi3d <- plotxsobject$phi3d - 2
                if (identical(key, "Left"))
                    plotxsobject$theta3d <- plotxsobject$theta3d + 2
                if (identical(key, "Right"))
                    plotxsobject$theta3d <- plotxsobject$theta3d - 2
                dev.set(expectationwindow)
                par(bg = "white")
                close.screen(all.screens = TRUE)
                do.call(plotxs, plotxsobject)
                assign(x = "plotxsobject", value = plotxsobject, envir = 
                    parent.frame())
            dev.flush()    
        }    
    }
    if (identical(key, "s")){
        expectationwindow <- get("expectationwindow", envir = parent.frame())
        selectorwindow <- get("selectorwindow", envir = parent.frame())
        plotxsobject <- get("plotxsobject", envir = parent.frame())
        Xc <- get("Xc", envir = parent.frame())
        Xc.cond <- plotxsobject$xc.cond
        timenow <- Sys.time()
        dev.set(selectorwindow)
        devsizesel <- dev.size()
        filename1 <- paste("snapshot_", gsub(":", ".", gsub(" ", "_", timenow)), 
            c("-expectation.pdf", "-condition.pdf"), sep = "") 
        pdf(filename1[1], width = 8.5, height = 8)
        close.screen(all.screens = TRUE)
        do.call(plotxs, plotxsobject)
        dev.off()
        pdf(filename1[2], width = devsizesel[1], height = devsizesel[2])
        close.screen(all.screens = TRUE)
        conditionselectors(Xc = Xc, type = get("selectortype", envir = 
            parent.frame()), method = get("Corder", envir = parent.frame()), 
            Xc.cond = Xc.cond)
        dev.off()
        cat(paste("\nSnapshot saved: '", filename1,"'", sep = ""))
        cat("\n")
    }
    dev.set(selectorwindow)
    points(NULL)
}
