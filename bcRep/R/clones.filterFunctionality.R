## Julia Bischof
## 15-10-2015

clones.filterFunctionality<-function(clones.tab=NULL, filter=c("productive","unproductive")){
  if(length(clones.tab)==0){
    stop("--> Clones.tab is missing")
  }
  if(length(filter)!=1 || !(filter %in% c("productive","unproductive"))){
    stop("--> Filter criterium is missing")
  }
  col.temp<-grep("func|Func",colnames(clones.tab))
  column<-vector()
  if(length(col.temp)>1){
    for(i in 1:length(col.temp)){
      if(length(grep("prod",clones.tab[,col.temp[i]]))>0){
        column<-col.temp[i]
      }
    }
  }else{
    column<-col.temp
  }
  if(filter=="productive"){
    filtergrep<-grep(filter, clones.tab[,column])[which(!(grep(filter, clones.tab[,column]) %in% grep("unprodutive", clones.tab[,column])))]
  }else{
    filtergrep<-grep(filter, clones.tab[,column])[which(!(grep(filter, clones.tab[,column]) %in% grep("^productive| productive", clones.tab[,column])))]
  }
  clones.filtered<-clones.tab[filtergrep,]
  
  return(clones.filtered)
}
