sim.equadist <- function(degree)
{

deg.sort<-sort(degree,decreasing=TRUE)

n.nodes <- length(degree)
adj.mat<-matrix(0,n.nodes,n.nodes)

erreur<-0
fin<-0

while(erreur==0&&fin==0){

for(i in 1:n.nodes)
{

if(deg.sort[i]!=0){

edge.max<-deg.sort[i]

choose<-c(1:n.nodes)
choose<-choose[deg.sort!=0]
choose<-choose[choose!=i]

link.nodes<-deg.sort[choose]

for(j in 1:edge.max){

deg.sort[i]<-deg.sort[i]-1

if(length(choose)==0){
erreur<-1
}else{
if(length(choose)==1) new.link<-choose else new.link<-sample(choose,1,prob=link.nodes)

deg.sort[new.link]<-deg.sort[new.link]-1

tmp<-c(1:length(choose))
remove<-tmp[choose==new.link]
new.nodes<-length(choose)
link.nodes[remove]<-link.nodes[remove]-1

choose<-choose[choose!=new.link]
link.nodes<-deg.sort[choose]


adj.mat[i,new.link]<-1
adj.mat[new.link,i]<-1


}

}

}

}

fin<-1

}
write.table(erreur)
return(adj.mat)

}






