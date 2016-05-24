see.accPrec.tck <- function(){



    if (!exists("slider.env")) 
        slider.env <- NULL
    suppressWarnings(rm(slider.env))
    slider.env <<- new.env()
    accuracy <- 0.5
    assign("accuracy", tclVar(accuracy), envir = slider.env)
    precision <- 0.5
    assign("precision", tclVar(precision), envir = slider.env)
    angle <- 180
    assign("angle", tclVar(angle), envir = slider.env)

dev.new(height = 4, width = 4)
par(mar=c(0,0,0,0))    
norm.refresh <- function(...){
        
        dev.hold()
        accuracy <- as.numeric(evalq(tclvalue(accuracy), envir = slider.env))
        precision <- as.numeric(evalq(tclvalue(precision), envir = slider.env))
        angle <- as.numeric(evalq(tclvalue(angle), envir = slider.env))

        
        plot(-1:1, -1:1, type = "n", xlab= "", ylab = "", yaxt = "n", xaxt = "n") 
        draw.circle(rep(0, 5), rep(0, 5), c(1, .8, .6, .4, .2), col = c("grey", "white", "grey", "white", "grey"))
        
        rad <- angle * pi/180
        
        xc <- sin(rad) * accuracy
        yc <- cos(rad) * accuracy
        
        x <- rnorm(10, sd = .25 * precision) + xc
        y <- rnorm(10, sd = .25 * precision) + yc
        
        points(x, y, pch = 19) 
                dev.flush()
}


    m <- tktoplevel()
    tkwm.title(m, "")
    tkpack(tklabel(m, text = "      Accuracy and Precision      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Accuracy (low to high)", width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 0, orient = "horiz", resolution = 0.1, showvalue = FALSE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = accuracy), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Precision (low to high)", width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 0, orient = "horiz", resolution = 0.1, showvalue = FALSE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = precision), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Angle", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir = slider.env)
    evalq(tkconfigure(e, textvariable = angle), envir = slider.env)
    
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")

}

