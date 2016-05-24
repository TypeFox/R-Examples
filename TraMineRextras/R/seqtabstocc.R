# Frequencies of state co-occurence patterns

seqtabstocc <- function(seqdata, ...){
    ssort <- t(apply(seqdata, MARGIN=1, sort))
    ssort.seq <- seqdef(ssort, states=alphabet(seqdata), labels=stlab(seqdata))
    sdss  <- seqdef(seqdss(ssort.seq), missing="%")

    tf <- seqtab(sdss, format="STS", ...)
    t <- attr(tf,"freq")
    rownames(t) <- gsub("-\\*","",rownames(t))
    return(t)
}
