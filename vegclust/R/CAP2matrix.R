CAP2matrix<-function(CAP, type="profile", classWidths=NULL) {
  n=length(CAP)
  if(!inherits(CAP,"CAP")) stop("Wrong input data. Use function CAP() first.")
  if(!is.null(classWidths)) {
    if(ncol(CAP[[1]])!=length(classWidths)) stop("Number of size classes in stratified species data does not match the number of elements in the vector of class widths")
  }
  else classWidths = rep(1,ncol(CAP[[1]])) #All classes have equal width
  if(n<2) stop("Wrong number of plot records (should be larger than one)")
  nstrata = length(classWidths)
  nspecies = nrow(CAP[[1]])
  spnames = row.names(CAP[[1]])
  if(type=="profile") {
    m = matrix(0, nrow=n,ncol=nstrata*nspecies)
    for(i in 1:n) {
      m[i,] = as.vector(t(CAP[[i]])*classWidths)
    }    
    colnames(m) = paste(as.character(gl(nspecies,nstrata, labels=spnames)),1:nstrata, sep="_")
  } else if(type=="volume") {
    m = matrix(0, nrow=n,ncol=nspecies)
    for(i in 1:n) {
      wm = t(CAP[[i]])*classWidths
      m[i,] = colSums(wm)
    }
    colnames(m) = spnames
  } else if(type=="abundance") {
    m = matrix(0, nrow=n,ncol=nspecies)
    for(i in 1:n) {
      for(j in 1:nspecies){
        m[i,j] = CAP[[i]][j,1]
      }
    }    
    colnames(m) = spnames
  }
  rownames(m) = names(CAP)
  return(m)
}