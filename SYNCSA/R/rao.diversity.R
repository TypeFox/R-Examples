rao.diversity<-function(comm,traits=NULL,phylodist=NULL,checkdata=TRUE,...){
	comm <- as.matrix(comm)
    N <- dim(comm)[1]
    S <- dim(comm)[2]
    tij2 <- 1 - diag(x = rep(1, S))
    if (!is.null(traits)) {
		if (checkdata) {
        if (is.null(rownames(traits))) {
            stop("\n Erro in row names of traits\n")
        }
        if (is.null(colnames(comm))) {
            stop("\n Erro in row names of comm\n")
        }
		match.names <- match(colnames(comm), rownames(traits))
        if (sum(is.na(match.names)) > 0) {
            stop("\n There are species from community data that are not on traits matrix\n")
        }
    		traits<-as.data.frame(traits[match.names,])
    	}
    	D1<-as.matrix(gowdis(x=traits,...))
    	S1<-1-D1
    	tij<-sqrt(1-S1)	
    }
    if (!is.null(phylodist)) {
		if (checkdata) {
        if (is.null(rownames(phylodist))) {
            stop("\n Erro in row names of phylodist\n")
        }
		if (is.null(colnames(phylodist))) {
            stop("\n Erro in column names of phylodist\n")
        }
        if (is.null(colnames(comm))) {
            stop("\n Erro in row names of comm\n")
        }
        match.names <- match(colnames(comm), colnames(phylodist))
        if (sum(is.na(match.names)) > 0) {
            stop("\n There are species from community data that are not on phylogenetic distance matrix\n")
        }
        phylodist <- as.matrix(phylodist[match.names, match.names])
    	}	
    	D1<-as.matrix(phylodist)
    	tij3<-D1/max(D1)
    }
	comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")
	inter<-comm%*%tij2
	SD<-rowSums(sweep(comm,1,inter,"*",check.margin=F))	
	if (!is.null(traits)){
		inter<-comm%*%tij
		RD<-rowSums(sweep(comm,1,inter,"*",check.margin=F))	
	}
	if (!is.null(phylodist)){
		inter<-comm%*%tij3
		FRD<-rowSums(sweep(comm,1,inter,"*",check.margin=F))	
	}
    Res<-list(Simpson=SD)
    if (!is.null(traits)){    
		Res<-list(Simpson=SD,FunRao=RD,FunRedundancy=SD-RD)
    }
    if (!is.null(phylodist)){
		Res<-list(Simpson=SD,PhyRao=FRD,PhyRedundancy=SD-FRD)
    }
    if (!is.null(phylodist) & !is.null(traits)){
		Res<-list(Simpson=SD,FunRao=RD,FunRedundancy=SD-RD,PhyRao=FRD,PhyRedundancy=SD-FRD)
    }
return(Res)    
}