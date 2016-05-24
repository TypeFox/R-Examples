v <-
function(x, u, n="L1", bj=TRUE){
    if (!(n %in% c("L1", "L0", "L2", "elastic")))
       stop("norm 'n' must be out of 'L0', 'L1', 'L2'. \n")

    old.contr <- getOption("contrasts")

    # u
    if (!is.factor(u))
       stop("effect modifier must be nominal or ordinal. \n")
    int <- rep(1, times = length(u))
    options(contrasts = c("contr.treatment", "contr.treatment"))
    dummies <- as.matrix(model.matrix(~ u)[,-1])

    # x
    if (is.factor(x)&&nlevels(x)>2)
       stop("varying coefficient not well defined. \n")
    options(na.action=na.pass)
    if (is.factor(x) && nlevels(x)==2)
       x <- model.matrix(~ x, contrasts = list(x="contr.effect"))[,-1]

    # design
    design <- cbind( int-rowSums(dummies) , dummies) * as.vector(x)
    colnames(design) <- paste(".",levels(u), sep="")
    options(contrasts = old.contr)
    return(design)
}

