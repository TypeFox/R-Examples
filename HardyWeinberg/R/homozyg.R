homozyg <- function (X) {
    # indicates for the elements of X whether they are homozygotes or not
    lab <- names(X)
    if (is.null(lab)) {
        warning("Genotypes are not labelled, default labels (AA, AB, BB) assumed.")
        lab <- c("AA", "AB", "BB")
    }
    if (length(unique(lab))!=3) stop("Genotype counts are not unequivocally labelled")
    coding <- "nocoding"
    if(all(lab %in% c("0","1","2"))) coding <- "coding1"
    if(all(nchar(lab)==2)) coding <- "coding2"
    status <- switch(coding, coding1 = lab == "0" | lab == "2",
                             coding2 = allele.name(lab, 1) == allele.name(lab, 2),
                     stop("Invalid labelling of the genotype counts"))
    return(status)
}
