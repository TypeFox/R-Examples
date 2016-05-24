IDD.IRR<-function(meta_de,ind_de)
{
  if(class(ind_de)!="list")
  {stop("ind_de should be a list")}
  deindst=Reduce(union,ind_de)
  DE=length(meta_de)
  gains=meta_de[which(!(meta_de %in% deindst))]
  IDD=length(gains)
  IDR=IDD/DE*100
  perte=which(!(deindst %in% meta_de))
  Loss=length(perte)
  IRR=Loss/length(deindst)*100
  res=c(DE,IDD,Loss,round(IDR,2),round(IRR,2))
  names(res)=c("DE","IDD","Loss","IDR","IRR")
  res
}