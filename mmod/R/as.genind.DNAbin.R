#' as.genind.DNAbin
#'
#'Convert a DNAbin object into a genind object
#'
#'@import pegas
#'@import adegenet
#'@param x object of class DNAbin
#'@param pops vector of population assignemnts for each sequence
#'@return genind
#'@examples
#'library(pegas)
#'data(woodmouse)
#'wm <- as.genind.DNAbin(woodmouse, rep(c("A", "B", "C"), each=5))
#'diff_stats(wm)
#'@export

as.genind.DNAbin <- function(x, pops){
    h <- haplotype(x)
    tab <- sapply(attr(h, 'index'), function(i) 
                 sapply(1:dim(x)[1], function(j) sum(i==j)))
    colnames(tab) <- paste("L1", 1:dim(tab)[2], sep=".")
    res <- genind(tab=tab, pop=pops)
    return(res)
}

