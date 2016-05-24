REGE.for<-function(
    M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
	E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("rege",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}


REGD.for<-function(
   M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 0 #initial dissimiliarity between vertices (default 0 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  E<-matrix(as.double(E),ncol=N, nrow=N)
	
	res<-.Fortran("regd",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}



REGE.ow.for<-function(
   M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
	E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regeow",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}

REGD.ow.for<-function(
   M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 0 #initial dissimiliarity between vertices (default 0 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  E<-matrix(as.double(E),ncol=N, nrow=N)
	
  res<-.Fortran("regdow",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dnM[[1]],dnM[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}





REGE.ownm.for<-function(
   M, #netowrk in form of a matrix or array (in case of two relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
  if(NR==1){
    M2<-array(NA,dim=c(N,N,2))
    M2[,,1]<-diag(1/apply(M,1,sum))%*%M
    M2[,,2]<-M%*%diag(1/apply(M,2,sum))
    M2[is.nan(M2)]<-0
    NR<-2
    if(length(dimnames(M))==2) dimN<-dimnames(M) else dimN<-c(list(NULL),list(NULL))
    dimnames(M2)<-c(dimN,list(c("out","in")))
    M<-M2
  } else{
    if(NR==2){
      cat("The first matrix will be used to evalueate outgoing arcs and the second to evaluate in ingoing arcs.\n")
    } else stop("This function is only suitable for evaluating two relations obtained as a row and column normalization of a single relation network. You have supplied more than two relations.\n")
  }

  E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regeownm",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}



REGE.ownm.diag.for<-function(
   M, #netowrk in form of a matrix or array (in case of two relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
  if(NR==1){
    M2<-array(NA,dim=c(N,N,2))
    M2[,,1]<-diag(1/apply(M,1,sum))%*%M
    M2[,,2]<-M%*%diag(1/apply(M,2,sum))
    M2[is.nan(M2)]<-0
    NR<-2
    if(length(dimnames(M))==2) dimN<-dimnames(M) else dimN<-c(list(NULL),list(NULL))
    dimnames(M2)<-c(dimN,list(c("out","in")))
    M<-M2
  } else{
    if(NR==2){
      cat("The first matrix will be used to evalueate outgoing arcs and the second to evaluate in ingoing arcs.\n")
    } else stop("This function is only suitable for evaluating two relations obtained as a row and column normalization of a single relation network. You have supplied more than two relations.\n")
  }

  E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regeownmdiag",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}





REGE.nm.for<-function(
   M, #netowrk in form of a matrix or array (in case of two relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
  if(NR==1){
    M2<-array(NA,dim=c(N,N,2))
    M2[,,1]<-diag(1/apply(M,1,sum))%*%M
    M2[,,2]<-M%*%diag(1/apply(M,2,sum))
    M2[is.nan(M2)]<-0
    NR<-2
    if(length(dimnames(M))==2) dimN<-dimnames(M) else dimN<-c(list(NULL),list(NULL))
    dimnames(M2)<-c(dimN,list(c("out","in")))
    M<-M2
  } else{
    if(NR==2){
      cat("The first matrix will be used to evalueate outgoing arcs and the second to evaluate in ingoing arcs.\n")
    } else stop("This function is only suitable for evaluating two relations obtained as a row and column normalization of a single relation network. You have supplied more than two relations.\n")
  }

  E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regenm",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}


REGE.nm.diag.for<-function(
   M, #netowrk in form of a matrix or array (in case of two relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
  if(NR==1){
    M2<-array(NA,dim=c(N,N,2))
    M2[,,1]<-diag(1/apply(M,1,sum))%*%M
    M2[,,2]<-M%*%diag(1/apply(M,2,sum))
    M2[is.nan(M2)]<-0
    NR<-2
    if(length(dimnames(M))==2) dimN<-dimnames(M) else dimN<-c(list(NULL),list(NULL))
    dimnames(M2)<-c(dimN,list(c("out","in")))
    M<-M2
  } else{
    if(NR==2){
      cat("The first matrix will be used to evalueate outgoing arcs and the second to evaluate in ingoing arcs.\n")
    } else stop("This function is only suitable for evaluating two relations obtained as a row and column normalization of a single relation network. You have supplied more than two relations.\n")
  }

  E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regenmdiag",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}






REGE.ne.for<-function(
    M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
	E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regene",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}


REGE.ow.ne.for<-function(
   M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
	E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regeowne",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}


REGE.ownm.ne.for<-function(
   M, #netowrk in form of a matrix or array (in case of two relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
  if(NR==1){
    M2<-array(NA,dim=c(N,N,2))
    M2[,,1]<-diag(1/apply(M,1,sum))%*%M
    M2[,,2]<-M%*%diag(1/apply(M,2,sum))
    M2[is.nan(M2)]<-0
    NR<-2
    if(length(dimnames(M))==2) dimN<-dimnames(M) else dimN<-c(list(NULL),list(NULL))
    dimnames(M2)<-c(dimN,list(c("out","in")))
    M<-M2
  } else{
    if(NR==2){
      cat("The first matrix will be used to evalueate outgoing arcs and the second to evaluate in ingoing arcs.\n")
    } else stop("This function is only suitable for evaluating two relations obtained as a row and column normalization of a single relation network. You have supplied more than two relations.\n")
  }

  E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regeownmne",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}





REGE.nm.ne.for<-function(
   M, #netowrk in form of a matrix or array (in case of two relations)
   iter = 3,
   E = 1 #initial similiarity between vertices (default 1 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  
  if(NR==1){
    M2<-array(NA,dim=c(N,N,2))
    M2[,,1]<-diag(1/apply(M,1,sum))%*%M
    M2[,,2]<-M%*%diag(1/apply(M,2,sum))
    M2[is.nan(M2)]<-0
    NR<-2
    if(length(dimnames(M))==2) dimN<-dimnames(M) else dimN<-c(list(NULL),list(NULL))
    dimnames(M2)<-c(dimN,list(c("out","in")))
    M<-M2
  } else{
    if(NR==2){
      cat("The first matrix will be used to evalueate outgoing arcs and the second to evaluate in ingoing arcs.\n")
    } else stop("This function is only suitable for evaluating two relations obtained as a row and column normalization of a single relation network. You have supplied more than two relations.\n")
  }

  E<-matrix(E,ncol=N, nrow=N)
  res<-.Fortran("regenmne",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}


REGD.ne.for<-function(
   M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 0 #initial dissimiliarity between vertices (default 0 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  E<-matrix(as.double(E),ncol=N, nrow=N)
	
	res<-.Fortran("regdne",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dimnames(M)[[1]],dimnames(M)[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}



REGD.ow.ne.for<-function(
   M, #netowrk in form of a matrix or array (in case of several relations)
   iter = 3,
   E = 0 #initial dissimiliarity between vertices (default 0 among all vertices).
){
	if(is.array(M)){
		dM<-dim(M)
		dnM<-dimnames(M)
		N<-dM[1]
		if (length(dM)==3) {
			NR<-dM[3]
		} else {
			if(length(dM)==2) {
				NR<-1
			} else stop("An array has wrong dimensions")
		}
	} else stop("M must be an array")
  M<-structure(as.double(M),dim=dM)
  dimnames(M)<-dnM
  E<-matrix(as.double(E),ncol=N, nrow=N)
	
  res<-.Fortran("regdowne",M = M, E = E, N = as.integer(N), NR = as.integer(NR), iter = as.integer(iter))
  Eall<-array(NA,dim=c(dim(E),2))
  Eall[,,1]<-E
  Eall[,,2]<-res$E
  dimnames(Eall)<-list(dnM[[1]],dnM[[2]],c("initial","final"))
  return(list(E=Eall[,,"final"],Eall=Eall,M=M,iter=iter))
}


