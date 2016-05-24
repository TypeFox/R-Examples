"loadnetwork2" <-
function(filename,useSparseMatrix=NULL,minN=50){
  if(is.character(filename)){
  	file<-file(description=filename,open="r")
  } else file<-filename
  n<-read.table(file=file,nrows=1)
  if(length(n)==2){
    n<-as.numeric(n[2])
    vnames<-read.table(file=file,nrows=n,as.is =TRUE)[,2]
    if(all(is.na(vnames))){
        vnames<-NULL
    } else vnames[is.na(vnames)]<-""

    if(is.null(useSparseMatrix)){
        useSparseMatrix<- n>=50
    }
    if(useSparseMatrix){
    	if(require(Matrix)){
    		M<-Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    	}else{
	        M<-matrix(0,nrow=n,ncol=n)
	        warning("Matrix package is not installed. Ordanary (dense) matrices will be used instead of sparse onse")    	
    	}
    }else{
        M<-matrix(0,nrow=n,ncol=n)
    }
    
    type=""
    
    while(TRUE){
    	line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
    	if(length(line)==0||sum(grep(patt="^ *$",x=as.character(line))==1)) break
    	if(substr(line[1],start=1,stop=1)=="*"){
    		type=line[1]
    		next
    	}else line<-as.double(line)
    	
    	if(type=="*Arcs"){
    		M[line[1],line[2]]<-line[3]
    	}else if(type=="*Edges") {
    		M[line[1],line[2]]<-line[3]
    		M[line[2],line[1]]<-line[3]
    	}	
    }
    dimnames(M)<-list(vnames,vnames)
  } else{
    n12<-as.numeric(n[2])
    n1<-as.numeric(n[3])
    n2<-n12-n1
    vnames<-read.table(file=file,skip=1,nrows=n12,as.is =TRUE)[,2]
    if(all(is.na(vnames))){
        vnames<-NULL
    } else vnames[is.na(vnames)]<-""

    if(is.null(useSparseMatrix)){
        useSparseMatrix<- n12>50
    }
    
    if(useSparseMatrix){
    	if(require(Matrix)){
    		M<-Matrix(0,nrow=n12,ncol=n12,sparse=TRUE)
    	}else{
    		warning("Matrix package is not installed. Ordanary (dense) matrices will be used instead of sparse onse")
    		M<-matrix(0,nrow=n12,ncol=n12)
    	}
    } else {
	M<-matrix(0,nrow=n12,ncol=n12)       
    }
    
    while(TRUE){
    	line<-scan(file = file, nlines =1,what="char",quiet =TRUE)
    	if(length(line)==0||sum(grep(patt="^ *$",x=as.character(line))==1)) break
    	if(substr(line[1],start=1,stop=1)=="*"){
    		type=line[1]
    		next
    	}else line<-as.double(line)
    	
	M[line[1],line[2]]<-line[3]
	M[line[2],line[1]]<-line[3]
    }
    dimnames(M)<-list(vnames,vnames)
    M<-M[1:n1,(n1+1):n12]    
  }
  return(M)
}
