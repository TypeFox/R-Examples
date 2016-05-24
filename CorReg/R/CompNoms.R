# ' compare les nomes pr?sents dans deux vecteurs
CompNoms<-function(Tnames=Tnames,Fnames=Fnames){
  nbbon=sum(duplicated(c(Tnames,Fnames)))
  max=length(Tnames)
  nbtrop=length(Fnames)-nbbon
  mank=max-nbbon
  return(list(nbbon=nbbon,max=max,nbtrop=nbtrop,mank=mank))
}