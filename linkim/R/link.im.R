link.im <- function(data,r,char=NULL,oneside=NULL,twoside=NULL,trace=NULL,...){

  if(is.null(twoside))twoside=TRUE
  if(is.null(oneside))oneside=TRUE
  if(is.null(trace))trace = FALSE

  RES <- linkim(data,r,char,oneside,twoside,trace)
  return(RES)
}



linkim <-
function(data,r,char=NULL,oneside=NULL,twoside=NULL,trace=NULL,...){ 
  if(is.null(twoside))twoside=TRUE
  if(is.null(oneside))oneside=TRUE
  if(is.null(trace))trace = FALSE

  if(is.null(char)){
	cn <- ncol(data)
  	rn <- nrow(data)
  	char <- NULL
  	for(i in 1:cn) char <- union(data[,i],char)
  	id.na = which(is.na(char))
  	if(length(id.na)==1)char <- sort(char[-id.na])
  }
  if(length(char)==2){
	cn <- ncol(data)
  	rn <- nrow(data)
	homo <- TRUE
	if(!(setequal(c(0,1),char))) data <- matrix(as.numeric(greplace(data,char,c(0,1))),rn,cn)
  }
  if(length(char)==3){
	cn <- ncol(data)
  	rn <- nrow(data)
	homo <- FALSE
	if(!(setequal(c(0,1,2),char))) data <- matrix(as.numeric(greplace(data,char,c(0,1,2))),rn,cn)
  }
  if(length(char)!=2 & length(char)!=3)stop("The marker data should be Homozygous(2 levels) or Heterozygous(3 levels)")
  
  if(max(r)>0.5)distance <- TRUE
	else distance <- FALSE
  
  dataNew = data

 if(twoside){
   for(i in 1:nrow(data)){
     if(trace)cat("individual",i,"\n")
     id1 = which(is.na(data[i,]))
     id2 = which(!is.na(data[i,]))
     if(length(id1)>0) {
	a = length(id1)
	NoStop = T
	while(a>0 & NoStop){
		if(id1[1] < id2[1]){
			id1 = id1[-1]
			a = length(id1)
		}
		if(a>0) if(id1[1] > id2[1]) NoStop = F		
	}
     }
     if(length(id1)>0) { 
	a = length(id1)
	b = length(id2)
	NoStop = T
	while(a>0 & NoStop){
		if(id1[a] > id2[b]){
			id1 = id1[-a]
			a = length(id1)
		}
		if(a>0) if(id1[a] < id2[b]) NoStop = F		
	}
     }
     while(length(id1)>0){	
	Aid = id1[1]-1
	id2del = which(id2 < id1[1])
	id2 = id2[-id2del]
	Bid = id2[1]
	id1del = which(id1 < id2[1])
	id1 = id1[-id1del]
	for(j in (Aid+1):(Bid-1)){
		vec = c(Aid,j,Bid)
		if(distance){		
			r0 = r[vec]
			r1 = r0[-1]-r0[-3]
			r2 = 0.5*(1-exp(-2*0.01*r1))
		}
		if(!distance)r2=r[vec[-1]]
	  if(homo){
		ComLoci = gcomb(data,index=vec,marker=c(0,1))
		comb = numeric()
		for(ii in 1:4)comb[ii] = ComLoci[ii,4]+ComLoci[9-ii,4]
		id = which(comb==max(comb))
		LociOrder1 = t(as.matrix(ComLoci[id,1:3])) 
		LociOrder2 = t(as.matrix(ComLoci[9-id,1:3]))
		Freq = freqmat(r2[1],r2[2],twoside=TRUE,cross=TRUE)
		ConFreq = Uu(Freq[2])				
		if(data[i,vec[1]]==LociOrder1[1]){
				if(data[i,vec[3]]==LociOrder1[3])p=ConFreq[1,1]
				if(data[i,vec[3]]!=LociOrder1[3])p=ConFreq[2,1]
		}
		if(data[i,vec[1]]!=LociOrder1[1]){
				if(data[i,vec[3]]==LociOrder1[3])p=ConFreq[3,1]
				if(data[i,vec[3]]!=LociOrder1[3])p=ConFreq[4,1]
		}		
		u = runif(1)
		if(u<=p)dataNew[i,vec[2]]=LociOrder1[2]
		if(u>p)dataNew[i,vec[2]]=LociOrder2[2]
		
	  }
	  if(!homo){
		ComLoci = gcomb(data,index=vec,marker=c(0,1,2))
		comb = numeric()
		g3 <- c(1,3,7,9)
		for(ii in 1:4)comb[ii] = ComLoci[g3[ii],4]+ComLoci[28-g3[ii],4]
		id = which(comb==max(comb))
		LociOrder1 = t(as.matrix(ComLoci[g3[id],1:3]))
		LociOrder3 = t(as.matrix(ComLoci[28-g3[id],1:3]))
		LociOrder2 <- t(as.matrix(ComLoci[14,1:3]))
		Freq = freqmat(r2[1],r2[2],twoside=TRUE,cross=TRUE,homo=homo)
		ConFreq = Uu(Freq[2])		
		if(data[i,vec[1]]==LociOrder1[1]){
				if(data[i,vec[3]]==LociOrder1[3]){
					p1=ConFreq[1,1]
					p2=ConFreq[1,2]
				}
				if(data[i,vec[3]]==LociOrder2[3]){
					p1=ConFreq[2,1]
					p2=ConFreq[2,2]
				}
				if(data[i,vec[3]]==LociOrder3[3]){
					p1=ConFreq[3,1]
					p2=ConFreq[3,2]
				}
		}
		if(data[i,vec[1]]==LociOrder2[1]){
				if(data[i,vec[3]]==LociOrder1[3]){
					p1=ConFreq[4,1]
					p2=ConFreq[4,2]
				}
				if(data[i,vec[3]]==LociOrder2[3]){
					p1=ConFreq[5,1]
					p2=ConFreq[5,2]
				}
				if(data[i,vec[3]]==LociOrder3[3]){
					p1=ConFreq[6,1]
					p2=ConFreq[6,2]
				}
		}		
		if(data[i,vec[1]]==LociOrder3[1]){
				if(data[i,vec[3]]==LociOrder1[3]){
					p1=ConFreq[7,1]
					p2=ConFreq[7,2]
				}
				if(data[i,vec[3]]==LociOrder2[3]){
					p1=ConFreq[8,1]
					p2=ConFreq[8,2]
				}
				if(data[i,vec[3]]==LociOrder3[3]){
					p1=ConFreq[9,1]
					p2=ConFreq[9,2]
				}
		}		
		u = runif(1)
		if(u<=p1)dataNew[i,vec[2]]=LociOrder1[2]
		if(u>p1 & u<=(p1+p2))dataNew[i,vec[2]]=LociOrder2[2]
		if(u>(p1+p2))dataNew[i,vec[2]]=LociOrder3[2]
		
	  }
	 }
     }
   }
   data = NULL
 }
 if(oneside){
	data = dataNew
	for(i in 1:nrow(data)){
		if(trace)cat("individual",i,"\n")
		id1 = which(is.na(data[i,]))
    		id2 = which(!is.na(data[i,]))
		if(length(id1)>0){
			if(id2[1]>id1[1]){
				for(j in 1:(id2[1]-1)){
					vec=c(j,id2[1])
					if(distance){		
						r0 = r[vec]
						r2 = r0[2]-r0[1]
						r2 = 0.5*(1-exp(-2*0.01*r2))
					}
					if(!distance)r2=r[vec[2]]
	  		  	    if(homo){
					ComLoci = gcomb(data,index=vec,marker=c(0,1))
					comb = numeric()
					for(ii in 1:2)comb[ii] = ComLoci[ii,3]+ComLoci[5-ii,3]
					id = which(comb==max(comb))
					LociOrder1 = t(as.matrix(ComLoci[id,1:2]))
					LociOrder2 = t(as.matrix(ComLoci[5-id,1:2]))
					if(data[i,vec[2]]==LociOrder1[2]) p = 1-r2
					if(data[i,vec[2]]!=LociOrder1[2]) p = r2
					u = runif(1)
					if(u<=p)dataNew[i,vec[1]]=LociOrder1[1]
					if(u>p)dataNew[i,vec[1]]=LociOrder2[1]
					
  			  	    }
	  		  	    if(!homo){
					ComLoci = gcomb(data,index=vec,marker=c(0,1,2))
					comb = numeric()
					g3 <- c(1,3)
					for(ii in 1:2)comb[ii] = ComLoci[g3[ii],3]+ComLoci[10-g3[ii],3]
					id = which(comb==max(comb))
					LociOrder1 = t(as.matrix(ComLoci[g3[id],1:2]))
					LociOrder3 = t(as.matrix(ComLoci[10-g3[id],1:2]))
					LociOrder2 = t(as.matrix(ComLoci[5,1:2]))
					Freq = freqmat(r2,r2,oneside=TRUE,homo=homo)
					ConFreq = Uu(Freq[2])		
					if(data[i,vec[2]]==LociOrder1[2]){
						p1 = ConFreq[4,1]
						p2 = ConFreq[4,2]
					}
					if(data[i,vec[2]]==LociOrder2[2]){
						p1 = ConFreq[5,1]
						p2 = ConFreq[5,2]
					} 
					if(data[i,vec[2]]==LociOrder3[2]){
						p1 = ConFreq[6,1]
						p2 = ConFreq[6,2]
					} 
					u = runif(1)
					if(u<=p1)dataNew[i,vec[1]]=LociOrder1[1]
					if(u>p1 & u<=(p1+p2))dataNew[i,vec[1]]=LociOrder2[1]
					if(u>(p1+p2))dataNew[i,vec[1]]=LociOrder3[1]
					
	  		  	    }
				}
				id1del = which(id1 < id2[1])
				id1 = id1[-id1del]
			}
			if(length(id1)>0){
				a = length(id1)
				for(j in id1[1]:id1[a]){
					vec = c((id1[1]-1),j)
					if(distance){		
						r0 = r[vec]
						r1 = r0[2]-r0[1]
						r1 = 0.5*(1-exp(-2*0.01*r1))
					}
					if(!distance)r1=r[vec[2]]
	  		  	   if(homo){
					ComLoci = gcomb(data,index=vec,marker=c(0,1))
					comb = numeric()
					for(ii in 1:2)comb[ii] = ComLoci[ii,3]+ComLoci[5-ii,3]
					id = which(comb==max(comb))
					LociOrder1 = t(as.matrix(ComLoci[id,1:2]))
					LociOrder2 = t(as.matrix(ComLoci[5-id,1:2]))		
					if(data[i,vec[1]]==LociOrder1[1]) p = 1-r1 	
					if(data[i,vec[1]]!=LociOrder1[1]) p = r1 
					u = runif(1)
					if(u<=p)dataNew[i,vec[2]]=LociOrder1[2]
					if(u>p)dataNew[i,vec[2]]=LociOrder2[2]
					
  			  	   }
	  		  	   if(!homo){
					ComLoci = gcomb(data,index=vec,marker=c(0,1,2))
					comb = numeric()
					g3 <- c(1,3)
					for(ii in 1:2)comb[ii] = ComLoci[g3[ii],3]+ComLoci[10-g3[ii],3]
					id = which(comb==max(comb))
					LociOrder1 = t(as.matrix(ComLoci[g3[id],1:2]))
					LociOrder3 = t(as.matrix(ComLoci[10-g3[id],1:2]))
					LociOrder2 = t(as.matrix(ComLoci[5,1:2]))
					Freq = freqmat(r1,r1,oneside=TRUE,homo=homo)
					ConFreq = Uu(Freq[2])		
					if(data[i,vec[1]]==LociOrder1[1]){
						p1 = ConFreq[1,1]
						p2 = ConFreq[1,2] 
					}
					if(data[i,vec[1]]==LociOrder2[1]){
						p1 = ConFreq[2,1]
						p2 = ConFreq[2,2]
					} 
					if(data[i,vec[1]]==LociOrder3[1]){
						p1 = ConFreq[3,1]
						p2 = ConFreq[3,2]
					} 
					u = runif(1)
					if(u<=p1)dataNew[i,vec[2]]=LociOrder1[2]
					if(u>p1 & u<=(p1+p2))dataNew[i,vec[2]]=LociOrder2[2]
					if(u>(p1+p2))dataNew[i,vec[2]]=LociOrder3[2]
					
	  		  	   }
				}	
			}
		}
	}
	data = NULL
 }
   
	if(length(char)==2){	
		if(!(setequal(c(0,1),char))){
			cn <- ncol(dataNew)
  			rn <- nrow(dataNew)
			dataNew <- matrix(as.numeric(greplace(dataNew,c(0,1)),char),rn,cn)
		} 
  	}
  	if(length(char)==3){
		if(!(setequal(c(0,1,2),char))){
			cn <- ncol(dataNew)
  			rn <- nrow(dataNew)
			dataNew <- matrix(as.numeric(greplace(dataNew,c(0,1,2),char)),rn,cn)
		}
  	}
	return(dataNew)

}
