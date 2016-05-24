"convnum" <- function(kasc)
{
    ## Verifications
    if (!inherits(kasc, "kasc"))
        stop("should be of class kasc")

    ## converts as df
    litab<-kasc2df(kasc)
    ## performs a dudi.mix
    dud<-dudi.mix(litab$tab, scannf=FALSE)

    ## converts the resulting table as kasc
    toto <- dud$tab
    names(toto)
    cw<-dud$cw
    scores <- df2kasc(toto, litab$index, kasc)

    ## output
    return(list(kasc=scores, weight=cw))
}

