
corrplot <- function(data, varnames=NULL, cutpts=NULL, abs=TRUE, details=TRUE,
                     n.col.legend=5, cex.col=0.7, cex.var=0.9, digits=1,
                     color=FALSE)
{

    # some check!
    if (is.matrix(data)|is.data.frame(data)){
    }
    else {
        stop ("Data must be a matrix or a data frame!")
    }
    if (sum(sapply(data, FUN=is.character))>0)
        stop ("Data contains non-numeric variables!")
    if (n.col.legend > 8)
        stop ("Suggestion: More than 8 levels of colors is difficult to read!")



    # prepare correlation matrix
    if (abs){
        z.plot <- abs(cor(data, data, use="pairwise.complete.obs"))
    }
    else{
        z.plot <- cor(data, data, use="pairwise.complete.obs")
    }

    if (is.null(varnames)){
        z.names <- dimnames(data)[[2]]
    }
    else{
        z.names <- varnames
    }

    triangleplot(x=z.plot, y=z.names, cutpts=cutpts, details=details,
                n.col.legend=n.col.legend,
                cex.col=cex.col, cex.var=cex.var,
                digits=digits, color=color)
}
