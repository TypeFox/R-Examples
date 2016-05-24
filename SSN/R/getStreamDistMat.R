getStreamDistMat <- function(x, Name = "obs")
{
  if(Name == "Obs") Name = "obs"
  if(class(x) != "SpatialStreamNetwork")
    return("Object not of class SpatialStreamNetwork")
  path = paste0(x@path,"/distance/",Name)
  flist = list.files(path)
  distMats = vector("list",length(flist))
  for(i in 1:length(flist)) {
    path1 = paste0(path,"/",flist[i])
    file_handle <- file(path1, open="rb")
    distmat <- unserialize(file_handle)
    close(file_handle)
    ordrow <- order(as.numeric(rownames(distmat)))
    ordcol <- order(as.numeric(colnames(distmat)))
    distmat = distmat[ordrow, ordcol, drop = F] 
    distMats[[i]] = distmat
    nameSplit = unlist(strsplit(flist[i],"[.]"))
    tname = nameSplit[1]
    for(j in 2:(length(nameSplit) - 1)) 
      tname = paste(tname, nameSplit[j], sep = ".")
    flist[i] = tname
  }
  names(distMats) = flist
  return(distMats)
}
