"as.treeshape.phylo" <-
function (x, model=NULL, p, ...) {

        phy <- x
        
	if (identical(model, NULL)) {convertion=FALSE}
	else {convertion=TRUE}
	
	if (missing(p)) {p=0.3}
	
	bin=is.binary.phylo(phy)

 	randomize=FALSE
	if (!(identical(bin,TRUE))){
		if (!convertion) {return(NULL)}
		else {
			randomize=TRUE
			for (i in nrow(bin):1) {
				new.nodes=bin[i,2]-2
				
                                for (j in 1:nrow(phy$edge)) {
                                        if (phy$edge[j,1] > bin[i,1])
						{phy$edge[j,1] = phy$edge[j,1] + new.nodes}
                                        if (phy$edge[j,2] > bin[i,1])
						{phy$edge[j,2] = phy$edge[j,2] + new.nodes}
				}


				if (model=="pda") {tmp=rpda(bin[i,2])}
				if (model=="yule") {tmp=ryule(bin[i,2])}
				if (model=="biased") {tmp=rbiased(bin[i,2],p)}
				if (model=="aldous") {tmp=raldous(bin[i,2])}

				tmp <- as.phylo(tmp)
                                
                                print(tmp$edge)

                                for (j in 1:length(tmp$edge)) {
					if (tmp$edge[j] > bin[i,2]) {
						tmp$edge[j] <- tmp$edge[j] + bin[i,1] - bin[i,2] -1
					}
				}
				
                                print(tmp$edge)

                                new.values =  phy$edge[phy$edge[,1]==bin[i,1],2]
                            
                                for(j in 1:bin[i,2]){
                                        tmp$edge[tmp$edge == j] = new.values[j]        
                                }                          

                                print(tmp$edge)

				current.line=new.nodes+1
				m=matrix(nrow=nrow(phy$edge) + new.nodes, ncol=2)
				
				m[1:new.nodes,]=tmp$edge[1:new.nodes,]
				for(line in 1:nrow(phy$edge)) {
					if (phy$edge[line,1]==bin[i,1]) {
						m[line+new.nodes,]=tmp$edge[current.line,]
						current.line=current.line+1
					} else {
						m[line+new.nodes,]=phy$edge[line,]
					}
				}
				phy$edge=m
				
                                
			}
		}
	}
        
	height <- (nrow(phy$edge)/2)
	merge <- matrix(0, ncol=2, nrow=height)
        tip.number = height + 1
        total <- 2 * tip.number
	
 
	for (i in nrow(phy$edge):1) {
		if (as.numeric(phy$edge[i,2])>tip.number) {
			tmp=total - phy$edge[i,2]
		} else {
			tmp= - phy$edge[i,2]
		}
		if (merge[total - phy$edge[i,1], 2]==0) {
			merge[total - phy$edge[i,1], 2]=tmp
		} else {
			merge[total - phy$edge[i,1], 1]=tmp
		}
		
	}
	
	res=treeshape(merge, phy$tip.label)
	if (randomize){class(res)=c('treeshape', 'randomized.treeshape')}
	
	res
}

