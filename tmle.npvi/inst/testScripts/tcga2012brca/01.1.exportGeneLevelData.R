library("tmle.npvi")

rootPath <- "geneData"
rootPath <- Arguments$getReadablePath(rootPath)

dataSet <- "tcga2012brca"
path <- file.path(rootPath, dataSet)
path <- Arguments$getReadablePath(path)

files <- list.files(path, pattern=".*chr17,.*.xdr")
idxs <- 1:150

if (FALSE) {
  files <- list.files(path, pattern=".*chr21,.*.xdr")
  idxs <- seq(along=files)
}

filenames <- files[idxs]
pathnames <- file.path(path, filenames)
obsList <- lapply(pathnames, loadObject)
snames <- gsub("\\.xdr", "", filenames)
names(obsList) <- snames
str(obsList)


tcga2012brca <- obsList
save(tcga2012brca, file="tcga2012brca.rda")

if (FALSE) {
opath <- "data"
opath <- Arguments$getWritablePath(opath)

  for (ff in seq(along=idxs)) {
    filename <- files[ff]
    pathname <- file.path(path, filename)
    obs <- loadObject(pathname)
    sname <- gsub("\\.xdr", "", filename)
    ofilename <- sprintf("%s,%s.txt", dataSet, sname)
    opathname <- file.path(opath, ofilename)
    write.table(obs, opathname, quote=FALSE, row.names=FALSE)
  }
}
