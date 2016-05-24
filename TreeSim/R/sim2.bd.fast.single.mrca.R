sim2.bd.fast.single.mrca <-
function(n,lambda,mu,rho,mrca){
		n1<-sample.int(n-1,1)
		t1<-1
		if (n1>1){
			t1<-sim2.bd.fast.single.origin(n1,lambda,mu,rho,mrca)
			#t1<-add.rootedge(t1,mrca)
			}
		t2<-1
		n2<-n-n1
		if (n2>1){
			t2<-sim2.bd.fast.single.origin(n2,lambda,mu,rho,mrca)
			#t2<-add.rootedge(t2,mrca)
			}	
	    if (class(t1) == "phylo" && class(t2) == "phylo"){
	    	t1$root.edge<-t1$edge.length[1]
			t1<-collapse.singles(t1)
	    	t2$root.edge<-t2$edge.length[1]
			t2<-collapse.singles(t2)
			t<- t1+t2
			t$tip.label <- paste("t", sample(t$Nnode+1), sep = "")
		} else if (class(t1) == "phylo" || class(t2) == "phylo")  {
			if (class(t2) == "phylo") {
				t1<-t2
				}
			t<-t1
			t$edge <- t1$edge + 1
			root<- t$edge[1,1]
			t$edge <- rbind(c(root,1),t$edge)
			t$edge.length <- c(mrca,t$edge.length)
			t$tip.label <- paste("t", sample(t$Nnode+1), sep = "")
		} else {
			t<- sim2.bd.fast(2,1,1,0,1)[[1]]
			t$edge.length <- t$edge.length*0+mrca
		} 
	phy2<-collapse.singles(t)
	phy2<-reorder(phy2)
	phy2
}

