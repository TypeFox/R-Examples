dropplt <- function (taxa,site,which=NULL) 
{
    if (!identical(row.names(taxa),row.names(site))) stop('data frames do not match')

    orig_taxa <- deparse(substitute(taxa))
    orig_site <- deparse(substitute(site))
    if (is.null(which)) {
        keep <- apply(site,1,function(x){!any(is.na(x))})
    } else {
        keep <- 1:nrow(taxa)
        keep <- keep[-which]
    }
    taxa <- taxa[keep,]
    site <- site[keep,]
    res <- list(taxa=taxa,site=site)
    attr(res,'call') <- match.call()
    attr(res,'orig_taxa') <- orig_taxa
    attr(res,'orig_site') <- orig_site
    res
}

dropspc <- function (taxa,minocc=0,minabu=0) 
{
    taxa <- taxa[,apply(taxa>minabu,2,sum)>minocc]
    attr(taxa,'call') <- match.call()
    attr(taxa,'minocc') <- minocc
    attr(taxa,'minabu') <- minabu
    taxa
}

