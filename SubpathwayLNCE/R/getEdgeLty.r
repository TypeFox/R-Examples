getEdgeLty<-function(graph){
edge.name<-E(graph)$subtype_name
edge.lty=rep("solid",len=length(edge.name))
for(i in seq(edge.name)){
  if(edge.name[i]=="indirect effect"){
     edge.lty[i]<-"longdash"
  }else if(edge.name[i]=="state change"){
     edge.lty[i]<-"longdash"
  }
}
#new!
if(length(edge.lty)==0) edge.lty="solid"
return(edge.lty)
}