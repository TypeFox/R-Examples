LTT.plot.gen<-function(trees,bound=10^(-12)){
if (length(trees)==1) {tree<-trees[[1]]
	origin<-1:length(tree)*0
	} else {origin <- trees[[2]]
		tree<-trees[[1]]
	}
	
brall<-vector()
bravg<-vector()
for (i in 1:length(tree)){
	if (class(tree[[i]]) == "phylo"){
		if (length(tree[[i]]$tip.label)==2) {
			br3<-c(-max(tree[[i]]$edge.length),2)
			br2<-c(-max(tree[[i]]$edge.length),2)
			if ((max(tree[[i]]$edge.length)-min(tree[[i]]$edge.length))>10^(-14)){
				br3<-rbind(br3,c(-(max(tree[[i]]$edge.length)-min(tree[[i]]$edge.length)),1))
				br2<-rbind(br2,c(-(max(tree[[i]]$edge.length)-min(tree[[i]]$edge.length)),-1))
				br3<-rbind(br3,c(0,1))	
			} else {br3<-rbind(br3,c(0,2))}
			if (origin[i]>0) {br3<-rbind(c(-origin[i],1),br3)
				br2<-rbind(c(-origin[i],1),br2)
				br2[2,2]<-1
			}
	} else {		
	br<-getx(tree[[i]],sersampling=1)
	del<-which(br[,1]<bound)
	br<-br[-del,]
	br2<-br[order(br[,1],decreasing=TRUE),]
	br2[which(br2[,2]==0),2]<- -1
	if (origin[i]==0){br2[1,2]<-2} else {br2<-rbind(c(origin[i],1),br2)}
	br2[,1]<- -br2[,1]
	br3<-cbind(br2[,1],cumsum(br2[,2]))
	br3<-rbind(br3,c(0,br3[length(br3[,1]),2]))}
	brall<-c(brall,list(br3))
	bravg<-rbind(bravg,br2)
	} else if (tree[[i]] == 1) {
		if (origin[i]>0) {
			br3<-rbind(c(-origin[i],1),c(-origin[i],1))
			brall<-c(brall,list(br3))
			bravg<-rbind(bravg,c(-origin[i],1))
		}
	}}
    bravg2 <- bravg[order(bravg[, 1]), ]
    bravg3 <- cbind(bravg2[, 1], cumsum(bravg2[, 2]))
    #print(bravg3)
    bravg3 <- rbind(bravg3, c(0, bravg3[length(bravg3[, 1]), 
        2]))
    #print(bravg3)
    bravg3[,2]<-bravg3[,2]/length(tree)  
    out <- c(list(bravg3), brall)
    out
# # bravg2<-bravg[order(bravg[,1]),]
# bravg3<-cbind(bravg2[,1],cumsum(bravg2[,2]))
# bravg3<-rbind(bravg3,c(0,bravg3[length(bravg3[,1]),2]))
# del<-vector()
# mini<-min(c(2*length(tree),length(bravg3[,1])-1))
# for (j in mini:1){
		 # if (abs(bravg3[j,1]-bravg3[(j+1),1])<10^(-8)) 
			 # {del<-c(del,j)
			  # }}
# bravg4<-bravg3[-del,]
# bravg4[,2]<-bravg4[,2]/length(tree)
# out<-c(list(bravg4),brall)
# out
}