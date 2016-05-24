matrify <- function(data)
{
    if (ncol(data) != 3) stop('data frame must have three column format')
    plt <- data[,1]
    spc <- data[,2]
    abu <- data[,3]
    plt.codes <- levels(factor(plt))
    spc.codes <- levels(factor(spc))
    taxa <- matrix(0,nrow=length(plt.codes),ncol=length(spc.codes))
    row <- match(plt,plt.codes)
    col <- match(spc,spc.codes)
    for (i in 1:length(abu)) {
        taxa[row[i],col[i]] <- abu[i]
    }
    taxa <- data.frame(taxa)
    names(taxa) <- spc.codes
    row.names(taxa) <- plt.codes
    taxa
}
