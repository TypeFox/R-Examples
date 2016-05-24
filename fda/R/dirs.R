dirs <- function(path='.', pattern=NULL, exclude=NULL, all.files=FALSE,
                 full.names=FALSE, recursive=FALSE, ignore.case=FALSE){
##
## 1.  mainDir <- dir(...)
##
  mD.all <- dir(path=path, pattern=pattern, all.files=all.files,
                 full.names=full.names, ignore.case=ignore.case)
  {
    if(is.null(exclude)) mainDir <- mD.all
    else {
      mD.excl <- dir(path=path, pattern=exclude, all.files=all.files,
                     full.names=full.names, ignore.case=ignore.case)
      mainDir <- mD.all[!(mD.all %in% mD.excl)]
    }
  }
  if(length(mainDir)<1)return(character(0))
##
## 2.  file.info(...)$isdir
##
  mD <- {
    if(full.names) mainDir
    else file.path(path, mainDir)
  }
  fi <- file.info(mD)
# 
  mainDirs <- mainDir[fi$isdir]
  nDirs <- length(mainDirs)
  if(nDirs<1) return(character(0))
##
## 3.  if(recursive)...
##
  if(!recursive)return(mainDirs)
# 
  mD.full <- mD[fi$isdir]
  Dirs <- vector('list', nDirs)
  names(Dirs) <- mainDirs
#
  for(i in 1:nDirs){
    id <- dirs(path=mD.full[i], pattern=pattern, exclude=exclude, 
               all.files=all.files, full.names=full.names,
               ignore.case=ignore.case)
    if(!full.names)
      id <- file.path(mainDirs[[i]], id)
    Dirs[[i]] <- c(mainDirs[i], id) 
  }
#
  unlist(Dirs, use.names=FALSE)
}
