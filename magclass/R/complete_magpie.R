

complete_magpie<-function(x,fill=NA) {
  full<-fulldim(x)[[2]]
  permute<-full[[3]]
  repeatit<-length(permute)
  if (length(full)>3) {
    for (i in 4:length(full)) {
      permute<-paste(rep(permute,each=length(full[[i]])),full[[i]],sep=".")
      repeatit<-length(permute)
    }
  }
  missing<-permute[!(permute%in%dimnames(x)[[3]])]
  if(length(missing)>0){
    add<-new.magpie(cells_and_regions = full[[1]],years = full[[2]],names = missing,fill=fill)
    out<-mbind(x,add)
  } else {out<-x}
  out<-out[,,order(getNames(out))]
  return(out)
}
