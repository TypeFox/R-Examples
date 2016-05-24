gcomb <-
function(data,index,marker=NULL,...){
        dat = data[,index]
	ok = complete.cases(dat)
	dat = dat[ok,]
	r = nrow(dat)
	c = ncol(dat)
	if (is.null(marker)){
		for(i in c) marker <- union(dat[,i],marker)
		marker <- sort(marker)
	}	
	MarkerPerm = perm(length(marker),c,v=marker)
	nr = nrow(MarkerPerm) 
	A = matrix(,r,nr)
	h=1
	for(k in 1: r){
      		for(i in 1:nr){
			for(j in 1:c)if(dat[k,j]!=MarkerPerm[i,j])h=h*0	  
			A[k,i]=h
			h = 1
		}
	}
	ColName=NULL
	MP=MarkerPerm
	for(i in 1:nr){
		for(j in 1:(c-1)){MP[i,j+1]<- paste(MP[i,j],MP[i,j+1],sep=":")}
		ColName = cbind(ColName,MP[i,j+1])      
	}	
	colnames(A)=ColName
	ComLoci <- cbind(MarkerPerm,colSums(A))
	A <- NULL
	MarkerPerm <- NULL
	MP <- NULL
	return(ComLoci)
 }
