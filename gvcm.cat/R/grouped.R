grouped <-
function(u, ...){

    old.contr <- getOption("contrasts")

    # u == factor
    if (is.factor(u)){

    options(contrasts = c("contr.treatment", "contr.treatment")) # dummy
    design <- as.matrix(model.matrix(~ u)[,-1])

    # design
    colnames(design) <- paste(".",levels(u)[-1], sep="")
    } else {

       stop("argument 'u' must be categorical. \n")
    
    # u und Rest metrische Kovariablen
#    call <- match.call()
    
#    if (length(call)<=2)
#       stop("argument 'u' must be either categorical or further to be grouped arguments are missing. \n")

#    f <- paste( " ~ ", deparse(call[[2]]), sep="") # ~ u
#    for (i in 3:length(names(call))){
#    f <- paste(f, " + ", deparse(call[[i]]), sep="")
#    }
    
#    design <- model.matrix(eval(parse(text=f)))[,-1]
    }
    
    # old options + return
    options(contrasts = old.contr)
    return(design)
}

