p <-
function(u, n="L1"){
    if (!(n %in% c("L1", "L0", "L2")))
       stop("norm 'n' must be out of 'L0', 'L1', 'L2'. \n")

    old.contr <- getOption("contrasts")

    # u
    if (is.factor(u)){

    options(contrasts = c("contr.treatment", "contr.treatment")) # dummy
    design <- as.matrix(model.matrix(~ u)[,-1])

    # design
    colnames(design) <- paste(".",levels(u)[-1], sep="")
    } else {
    design <- matrix(model.matrix(~ u)[,-1],ncol=1)    
    }
    options(contrasts = old.contr)
    return(design)
}

elastic <-
function(u){
    old.contr <- getOption("contrasts")

    # u
    if (is.factor(u)){

    options(contrasts = c("contr.treatment", "contr.treatment")) # dummy
    design <- as.matrix(model.matrix(~ u)[,-1])

    # design
    colnames(design) <- paste(".",levels(u)[-1], sep="")
    } else {
    design <- matrix(model.matrix(~ u)[,-1],ncol=1)
    }
    options(contrasts = old.contr)
    return(design)
}

SCAD <- elastic