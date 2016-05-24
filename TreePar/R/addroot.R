addroot<-function(tree,rootlength){
height<-max(getx(tree,sersampling=1)[,1])
root<-0
summary<-rbind(getx(tree,sersampling=1),c((height+rootlength),1))
timesorig<-summary[,1]
ttypeorig<-summary[,2]
times<-summary[,1]
ttype<-summary[,2]
phylo<-tree
cut <- phylo$edge[1, 1]
        for (i in 1:length(phylo$edge[, 1])) {
            if (phylo$edge[i, 1] >= cut) {
                phylo$edge[i, 1] <- phylo$edge[i, 1] + 1
            }
            if (phylo$edge[i, 2] >= cut) {
                phylo$edge[i, 2] <- phylo$edge[i, 2] + 1
            }
        }
        phylo$edge <- rbind(c(cut, phylo$edge[1, 1]), phylo$edge)
        phylo$edge.length <- c(rootlength, phylo$edge.length)
test<-phylo
test$Nnode<-length(tree$tip.label)
test
}