`plot.dlcross` <- 
function(x, chr, pheno.col, ...)
{
    dots <- list(...)

    if (missing(x)) 
        stop("x is a required argument")

    if (!inherits(x, "dlcross")) 
        stop("input$genCross is not of class \"dlcross\"")

    has.map <- sum(is.na(unlist(x$map)))==0
    nplots <- has.map+min(x$nphe, 3)

    if (x$nphe>3) cat("Will only plot 3 phenotypic variables, please use option pheno.col to restrict set of phenotypes")

    par(mfrow=c(2,2))

    if (has.map)  plot(x$map)

    # histogram of phenotypes (or barplot if categorical)
    if (missing(pheno.col)) pheno.col <- 1:(min(x$nphe, 3))

    pheno.col <- pheno.col[-match(x$idname, colnames(x$dfMerged)[pheno.col])]

    for (i in pheno.col)
    if (is.numeric(x$dfMerged[,i]))
	hist(x$dfMerged[,i], col="violetred3", xlab=paste("phe", i), main=names(x$dfMerged[,i])) else
	plot(x$dfMerged[,i], col="royalblue", xlab=paste("phe", i), main=names(x$dfMerged[,i]))
}
 
