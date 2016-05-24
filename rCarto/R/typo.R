typo <-
function(pt,var, palCol,NACol){
  pt[,"varQuali"]<-as.factor(pt[,var])
  colours<-brewer.pal(n=length(levels(pt[,"varQuali"])),name=palCol)
  legQuali<-NULL
  legQuali$val<-levels(pt[,"varQuali"])
  legQuali$col<-colours
  levels(pt[,"varQuali"])<-colours
  pt[,"varQuali"]<-as.character(pt[,"varQuali"])
  pt[is.na(pt$varQuali)==TRUE,"varQuali"]<-NACol
  if ((length(pt[is.na(pt[,var]),var])>0)==TRUE){
    pdd<-1
  }else{pdd<-0}
  return(list(pt,legQuali,pdd))
}
