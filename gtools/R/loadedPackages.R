loadedPackages <- function(silent=FALSE)
{
  packageNames    <- loadedNamespaces()
  packageVersions <- sapply(packageNames, function(package) paste(packageVersion(package), sep=".") )
  packagePaths    <- find.package(packageNames)
  inSearchPath    <- match(packageNames, gsub('^package:', '', grep('^package:', search(), value=TRUE)))
  retval <- data.frame(Name=packageNames, Version=packageVersions, Path=packagePaths, SearchPath=inSearchPath)
  retval$SearchPath <- na.replace(retval$SearchPath, '-')
  retval <- retval[order(inSearchPath),]
  if(!silent) print(retval)
  retval
}
