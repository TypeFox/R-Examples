
.read.bugs<-function (codafiles, ...) 
{
    mcmc.list(lapply(codafiles, read.coda, index.file = file.path(dirname(codafiles[1]), 
        "codaIndex.txt"), ...))
}
