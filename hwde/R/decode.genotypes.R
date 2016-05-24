"decode.genotypes" <-
function (genotype) 
{
    gname <- deparse(substitute(genotype))
    genotype <- factor(genotype)
    gtypes <- levels(genotype)
    charswide <- nchar(gtypes)
    nabc <- length(gtypes)
    if ((nabc != 3) | any(charswide != 2)) 
        stop(paste("Illegal codes", paste(gtypes, collapse = ", "), 
            "at locus", gname))
    ch1 <- substring(gtypes, 1, 1)
    ch2 <- substring(gtypes, 2, 2)
    AA <- gtypes[ch2 == ch1][1]
    Aa <- gtypes[ch1 != ch2]
    aa <- gtypes[ch2 == ch1][2]
    ord <- match(c(AA, Aa, aa), gtypes)
    id <- c(2, 1, 0)
    names(id) <- c(AA, Aa, aa)
    n <- length(genotype)
    idloc <- c(1, 2, 1)
    names(idloc) <- c(AA, Aa, aa)
    oset <- idloc[as.character(genotype)]
    ma <- id[as.character(genotype)]
    idaa <- c(1, 0, 0)
    names(idaa) <- c(AA, Aa, aa)
    maa <- idaa[as.character(genotype)]
    list(oset = oset, ma = ma, maa = maa, types = c(AA, Aa, aa))
}
