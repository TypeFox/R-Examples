candidate_genes=function(res,pct=0.01,prioritizing){
  len=floor(length(res$zi.order[,1])*pct) ###number of candidate genes#####
  idx1=as.character(names(res$genesets.clear))   ###indexes of modules####
  idx2=as.character(res$zi.ordered$gene[1:(len)])###indexes of seed genes for the candidate modules####
  idx=match(idx2,idx1)
  genome.list=res$genesets.clear[idx]             
  a=unlist(genome.list)
  if (prioritizing ==F){
  index<-duplicated(a)
  list=data.frame(a[!index])
  names(list)="gene"
  return(list)}
  if (prioritizing ==T){
    t=data.frame(table(a))
    t=t[order(t$Freq,decreasing=T),]
    list=t[t$Freq>1,]
    names(list)[1]="gene"
    return(list)
  }
}

