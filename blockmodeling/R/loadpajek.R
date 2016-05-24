loadpajek<-function(filename){
	if(is.character(filename)) {file<-file(description=filename,open="r")
	}else file<-filename
	res<-list(Networks=list(),Partitions=list(),Vectors=list(),Permutation=list())
	nblanklines=0
	while(TRUE){
		line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
		if(length(line)==0) {
			nblanklines=nblanklines+1
			if (nblanklines>10) break
			next
		}
		nblanklines=0
		if(sum(grep(patt="^ *$",x=as.character(line))==1)) next
		if(line[1]=="*Matrix" || line[1]=="*Network"){
			objName<-paste(line[-1],collapse=" ")
			if(line[1]=="*Matrix"){
				readObj<-loadmatrix(file)
			}else readObj<-loadnetwork2(file)

			if(objName %in% names(res[["Networks"]])){
				i<-1
				while(TRUE){
					if(paste(objName,"Ver",i) %in% names(res[["Networks"]])) break
					i<-i+1
				}
				objName<-paste(objName,"Ver",i)
			}
			res[["Networks"]]<-c(res[["Networks"]],list(readObj))
			names(res[["Networks"]])[length(res[["Networks"]])]<-objName
		} else if(line[1]=="*Vector" || line[1]=="*Permutation" || line[1]=="*Partition"){
			objName<-paste(line[-1],collapse=" ")
			readObj<-loadvector2(file)
			if(line[1]=="*Vector"){
				if(objName %in% names(res[["Vectors"]])){
					i<-1
					while(TRUE){
						if(paste(objName,"Ver",i) %in% names(res[["Vectors"]])) break
						i<-i+1
					}
					objName<-paste(objName,"Ver",i)
				}
				res[["Vectors"]]<-c(res[["Vectors"]],list(readObj))
				names(res[["Vectors"]])[length(res[["Vectors"]])]<-objName
			} else if(line[1]=="*Permutation"){
				if(objName %in% names(res[["Permutations"]])){
					i<-1
					while(TRUE){
						if(paste(objName,"Ver",i) %in% names(res[["Permutations"]])) break
						i<-i+1
					}
					objName<-paste(objName,"Ver",i)
				}
				res[["Permutations"]]<-c(res[["Permutations"]],list(readObj))
				names(res[["Permutations"]])[length(res[["Permutations"]])]<-objName
			} else if(line[1]=="*Partition"){
				if(objName %in% names(res[["Partitions"]])){
					i<-1
					while(TRUE){
						if(paste(objName,"Ver",i) %in% names(res[["Partitions"]])) break
						i<-i+1
					}
					objName<-paste(objName,"Ver",i)
				}
				res[["Partitions"]]<-c(res[["Partitions"]],list(readObj))
				names(res[["Partitions"]])[length(res[["Partitions"]])]<-objName
			}
		}
	
	}
	return(res)
	close(file)
}
