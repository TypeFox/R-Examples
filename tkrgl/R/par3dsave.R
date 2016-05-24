par3dsave <- function(params = c("userMatrix", "scale", "zoom", "FOV"), times = FALSE, dev = rgl.cur()) {

    results <- list()
    for (n in params) results[[n]] <- list()
    if (times) {
    	start <- proc.time()[3]
    	results$times <- numeric(0)
    }
    
    RecordParms <- function() {
    	values <- par3d(params)
    	if (length(params) == 1) {
    	    values <- list(values)
    	    names(values) <- params
    	}
    	for (n in params) results[[n]] <<- c(results[[n]], list(values[[n]]))
    	if (times) results$times <<- c(results$times, proc.time()[3] - start)
    }	
    base <- tktoplevel()
    tkwm.title(base, "par3d")
    
    text <- tklabel(base, text="Click on Record to save par3d parameters.",
                          justify="left",
                          wraplength="2i")
    frame <- tkframe(base)	
    save <- tkbutton(frame, text="Record", command=RecordParms)		

    quit <- tkbutton(frame,text="Quit", command=function()tkdestroy(base))

    tkpack(save, quit, side="left")
    tkpack(text, frame)
    
    tkwait.window(base)
    
    results
}