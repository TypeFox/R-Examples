## Julia Bischof
## 15-10-2015


clones.filterJunctionFrame<-function(clones.tab=NULL, filter=c("in-frame","out-of-frame")){
  if(length(clones.tab)==0){
    stop("--> Clones.tab is missing")
  }
  if(length(filter)!=1 || !(filter %in% c("in-frame","out-of-frame"))){
    stop("--> Filter criterium is missing")
  }
  col.temp<-grep("junc|JUNC|Junc",colnames(clones.tab))
  column<-vector()
  if(length(col.temp)>1){
    for(i in 1:length(col.temp)){
      if(length(grep("-frame",clones.tab[,col.temp[i]]))>0){
        column<-col.temp[i]
      }
    }
  }else{
    column<-col.temp
  }
  if(filter=="in-frame"){
    filtergrep<-grep(filter, clones.tab[,column])[which(!(grep(filter, clones.tab[,column]) %in% grep("out-of-frame", clones.tab[,column])))]
  }else{
    filtergrep<-grep(filter, clones.tab[,column])[which(!(grep(filter, clones.tab[,column]) %in% grep("in-frame", clones.tab[,column])))]
  }
  clones.filtered<-clones.tab[filtergrep,]
  
  return(clones.filtered)
}
