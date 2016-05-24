library("tmle.npvi")
log <- Arguments$getVerbose(-8, timestamp=TRUE)

dataSet <- "tcga2012brca"
path <- file.path("geneData", dataSet)
path <- Arguments$getReadablePath(path)
files <- list.files(path)

nas <- sapply(files, function(ff) {
    obs <- loadObject(file.path(path, ff))
    sum(is.na(obs))
})

files <- files[which(nas==0)]

descr <- list(thresh=2e-2,
              f=identity,
              flavor="learning",
              nodes=1, ##3,
              iter=10,
              cvControl=2,
              stoppingCriteria=list(mic = 0.001, div = 0.001, psi = 0.01),
              nMax=30)

fileout <- sprintf("%s,tmle.npvi,%s,2014-11-21.rda", dataSet, descr$flavor)

mc.cores <- 3

nms <- unlist(strsplit(files, split=".xdr"))

TMLE <- parallel::mclapply(seq(along=files), mc.cores=mc.cores, FUN=function(ii) {
    ff <- files[ii]
    ## loading the data
    print(ii)
    print(ff)
    pathname <- file.path(path, ff)
    
    obs <- loadObject(pathname)
    nbcov <- ncol(extractW(obs))
    if (nbcov==1) {
        colnames(obs) <- c("Y", "X", "W")
    }

    ## thresholding copy number data
    whichSmall <- which(abs(obs[, "X"]) <= descr$thresh)
    obs[whichSmall, "X"] <- 0

    ##
    tmle <- try(tmle.npvi(obs=obs, f=descr$f, flavor=descr$flavor,
                          stoppingCriteria=descr$stoppingCriteria,
                          cvControl=descr$cvControl, nMax=descr$nMax))
    if (inherits(tmle, "try-error")) {
        return(attr(tmle, "condition"))
    } else {
        return(list(nbcov=nbcov, hist=getHistory(tmle)))
    }
})

names(TMLE) <- nms
save(descr, TMLE, file=fileout)
