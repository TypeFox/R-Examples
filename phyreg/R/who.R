who <-
function(phy,nodes,find="mm"){
get_daughters<-function(phy,node){myl<-list(); if(node==max(phy) || phy[[node]]!=0L) for (kk in 1:(nnodes-1)) if (phy[[kk]]==node) myl<-union(myl,kk); return(myl)}
get_species<-function(phy,node){
  myl<-list()
  if(node==max(phy) || phy[[node]]!=0L) 
  for (jj in 1:nspec){ 
    kk<-jj; if(phy[[kk]]!=0) while (kk<node){kk<-phy[[kk]]}
    if (kk==node) myl<-union(myl,jj)
    }
   return(myl)
  }
  get_mm<-function(phy,node){
    myll<-get_species(phy,node)
    return(range(myll))
   }
  
 temp<-phy+max(phy)*(phy==0L)
 nspec<-min(temp)-1
 nnodes<-max(phy)
 tellsp<-function(node,myll){ cat(paste0("Node ",format(node)," contains species ", paste0(myll,collapse=", "),"\n"))}
 telldg<-function(node,myll){ cat(paste0("Node ",format(node)," has daughters ", paste0(myll,collapse=", "),"\n"))}

if (missing(nodes)) nodes<-(nspec+1):nnodes
if (any(nodes!=round(nodes)) || any(nodes<=0) || any(nodes>max(phy)))
 stop("The list of nodes must be integers, positive, and less than the highest node-number in the supplied phylogeny")

if (find!="mm") {
if (length(nodes)==1)  {

   {if (find=="sp") if (nodes<=nspec) cat(paste0("Node ",format(nodes), " is a species\n")) else
     if (nodes==max(phy)) cat(paste0("Node ",format(nodes)," is the root and contains all included species")) else tellsp(nodes,get_species(phy,nodes)) else  
        telldg(nodes,get_daughters(phy,nodes))}
     } else
 for (ii in 1:length(nodes)) {
   if (nodes[[ii]]<=nspec) cat(paste0("Node ",format(nodes[[ii]])," is a species\n")) else
   if (find=="sp") {
    if (nodes[[ii]]==max(phy)) cat(paste0("Node ",format(nodes[[ii]])," is the root and contains all (included) species")) else tellsp(nodes[[ii]],get_species(phy,nodes[[ii]]))} else  telldg(nodes[[ii]],get_daughters(phy,nodes[[ii]]))}
    
} else {  #  mm case

mm<-array(0L,c(length(nodes),2))
for (ii in 1:length(nodes)) mm[ii,]<-get_mm(phy,nodes[[ii]])
 myout<-data.frame(ID = nodes, minsp=mm[,1],maxsp=mm[,2])
  return(myout)
}# finish "else"
}
