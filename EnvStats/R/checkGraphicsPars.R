checkGraphicsPars <-
function (...) 
{
    gp.arg.list <- list(...)
    gp.names <- names(gp.arg.list)
    n.gp <- length(gp.arg.list)
    gen.gp.list <- list()
    if (n.gp > 0) {
        index <- is.na(match(gp.names, c(.Gr.pars.high.and.title, 
            .Gr.pars.general)))
        if (any(index)) 
            stop(paste("Unknown or ambiguous high-level or general graphics", 
                "parameter(s):\n\t\t", paste(gp.names[index], 
                  collapse = ", "), "\n\tNote:  You cannot abbreviate graphics parameters"))
        index <- match(gp.names, .Gr.pars.general)
        if (any(!is.na(index))) 
            gen.gp.list <- gp.arg.list[!is.na(index)]
    }
    list(gp.arg.list = gp.arg.list, gp.names = gp.names, n.gp = n.gp, 
        gen.gp.list = gen.gp.list)
}
