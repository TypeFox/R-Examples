
add_dimension<-function(x,dim=3.1,add="new", nm="dummy"){
  x<-clean_magpie(x)
  dim<-as.numeric(dim)
  olddim<-old_dim_convention(dim)
  if (olddim<3) stop("Dimensions below 3 are currently not supported by add_dimensions.")
  if (is.null(getNames(x))){getNames(x)<-"NA"}
  firstnm<-nm[1]
  
  separate<-strsplit(dimnames(x)[[3]],split = "\\.")
  newnm<-NULL
  for (i in 1:length(separate)) {
    newnm<-c(newnm,paste(append(separate[[i]],firstnm,after=olddim-3),collapse="."))
  }
  dimnames(x)[[3]]<-newnm
  names(dimnames(x))[3] <- paste(append(
    strsplit(names(dimnames(x))[3],split="\\.")[[1]]
    ,add,after=olddim-3),collapse=".")

  if (length(nm)>1) {
    x<-add_columns(x,addnm = nm[2:length(nm)],dim=dim)
    tmp<-list(neu=nm[1])
    names(tmp)<-add
    x[,,]<-mselect(x,tmp)
  }
  return(x)
}
