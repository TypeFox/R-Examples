print.pyRAD.loci <- function (x, ...) 
{
    print(paste("pyRAD.loci object read from", x$file.read, "on", x$timestamp))
    print(paste("Contains:"))
    print(paste("--", length(unique(x$tips)[!unique(x$tips) %in% c("", "//")]), "individuals"))
    print(paste("--", length(unique(x$locus.index)[unique(x$locus.index) != ""]), "loci"))
}
