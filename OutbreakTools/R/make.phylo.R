

################
## make.phylo ##
################

setMethod("make.phylo", "obkData", function(x, locus=NULL, result=c("obkData","multiPhylo"),
                                            model = "N", pairwise.deletion = FALSE, method=nj,
                                            plot=FALSE, ask=TRUE, ...){
    if(get.nlocus(x)==0){
        warning("No DNA sequences in the data.")
        return(NULL)
    }
    ## if(is.null(palette)){
    ##     ##palette <- colorRampPalette(brewer.pal(11, "RdYlGn"))
    ## }
    result <- match.arg(result)

    ## GET DNA SEQUENCES ##
    if(is.null(locus)) locus <- 1:get.nlocus(x)
    N <- length(locus)

    ## GET DISTANCES ##
    D <- lapply(get.dna(x, locus=locus), function(e) dist.dna(e, model=model, pairwise.deletion = pairwise.deletion))

    ## GET TREES ##
    tre <- lapply(D, function(e) ladderize(method(e)))

    ## PLOT (OPTIONAL) ##
    if(plot){
        par(ask=ask, xpd=TRUE)
        for(i in 1:N){
            plot(tre[[i]], ...)
            title(paste("locus:", locus[i]))
        }
    }


    ## RETURN OBJECT ##
    class(tre) <- "multiPhylo"
    if(result=="multiPhylo") return(tre) # returned object is multiPhylo

    x@trees <- tre # returned object is obkData
    return(x)

}) # end make.phylo

