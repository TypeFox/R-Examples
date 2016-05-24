metabolomicsChanges <- function(n=5)
{
    writeLines(readLines(system.file("doc", "metabolomicsChanges.txt",
        package="metabolomics"), n=n)
    )
}
