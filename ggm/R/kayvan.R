#### Functions by Kayvan Sadeghi 2011-2012
## May 2012 Changed the return values of some functions to TRUE FALSE
require(graph)
require(igraph)

######################################################################
######################################################################
rem<-function(a,r){ # this is setdiff (a, r)
	k<-0
	b<-a
	for (i in a){
		k<-k+1
		for(j in r){
			if(i==j){
				b<-b[-k]
				k<-k-1
				break}
			
			}
		}
	return(b)
}
#############################################################################
#############################################################################
#SPl<-function(a,alpha){
#	a<-sort(a)
#	alpha<-sort(alpha)
#	r<-c()
#	if (length(alpha)>0){
#		for(i in 1:length(a)){
#			for(j in 1:length(alpha)){
#				if(a[i]==alpha[j]){
#					r<-c(r,i)
#					break}}}}
#	return(r)
#}
############################################################################### 
# Finds indices of b in sorted a
                        
'SPl' = function(a, b) (1:length(a))[is.element(sort(a), b)]
##############################################################################
RR<-function(a){   ## This is unique(a)
	a<-sort(a)
	r<-a[1]
	i<-1
	while(i<length(a)){
		if(a[i]==a[i+1] ){
			i<-i+1}
		else{
			r<-c(r,a[i+1])
			i<-i+1	
			}}
	return(r)
}
#################################################################################
#################################################################################
`grMAT` <- function(agr)
{
	 if (class(agr) == "graphNEL") {
        	agr<-igraph.from.graphNEL(agr)
       }
	 if (class(agr)== "igraph"){
		return(get.adjacency(agr, sparse = FALSE))
	 }
	 if(class(agr) == "character"){
		if (length(agr)%%3!=0){
			stop("'The character object' is not in a valid form")}
		seqt<-seq(1,length(agr),3)
		b<-agr[seqt]
		agrn<- agr[-seqt]
		bn<-c()
		for(i in 1:length(b)){
			if(b[i]!="a" && b[i]!="l" & b[i]!="b"){
				stop("'The numeric object' is not in a valid form")}
			if(b[i]=="l"){
				bn[i]<-10}
			if(b[i]=="a"){
				bn[i]<-1}
			if(b[i]=="b"){
				bn[i]<-100}
			}
		Ragr<-RR(agrn)
		ma<-length(Ragr)
		mat<-matrix(rep(0,(ma)^2),ma,ma)
		for(i in seq(1,length(agrn),2)){
			if((bn[(i+1)/2]==1 & mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]%%10!=1)|(bn[(i+1)/2]==10 & mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]%%100<10)|(bn[(i+1)/2]==100 & mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]<100)){
				mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]<-mat[SPl(Ragr,agrn[i]),SPl(Ragr,agrn[i+1])]+bn[(i+1)/2]
			if(bn[(i+1)/2]==10 | bn[(i+1)/2]==100){
				mat[SPl(Ragr,agrn[i+1]),SPl(Ragr,agrn[i])]<-mat[SPl(Ragr,agrn[i+1]),SPl(Ragr,agrn[i])]+bn[(i+1)/2]}}
		}
		rownames(mat)<-Ragr
		colnames(mat)<-Ragr
	}			
	return(mat)
}		

##############################################################################
##############################################################################
`RG` <- function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{	
	if(class(amat) == "igraph" || class(amat) == "graphNEL" || class(amat) == "character") {
		amat<-grMAT(amat)}
	if(class(amat) == "matrix"){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}	
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(class(amat) != "matrix") {	
	stop("'object' is not in a valid form")}


	S<-C
	St<-c()
	while(identical(S,St)==FALSE)
	{
		St<-S
		for(j in S){
			for(i in rownames(amat)){
				if(amat[i,j]%%10 == 1){
					S<-c(i,S[S!=i])}
				}
			}	
	}	
			

	amatr<-amat
	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr

############################################################1
			amat21 <- amatr
 	   for (kk in M) {
		  for (kkk in 1:ncol(amatr)) {
			amat31<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat31[kk,kkk]<-amatr[kkk,kk]
				amat31[kkk,kk]<-amatr[kk,kkk]}
        		
       		 idx <- which(amat31[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat31[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amat21[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}
    		
		amatr<-amat21	

################################################################2	

		amat22 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat22[idx[ii], idy[jj]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat22+amatr

#################################################################3

		amat33<- t(amatr)
		amat23 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat33[, kk]%%10 == 1)
			idy <- which(amat33[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]%%10==0 & idx[ii]!=idy[jj]){
                  amat23[idy[jj], idx[ii]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat23+amatr

##################################################################4
			amat24 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%100>9)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat24[idx[ii], idy[jj]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat24+amatr

####################################################################5		



		 amat35<- t(amatr)
   		 amat25 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amat35[, kk]%%10 == 1)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat25[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat25+t(amat25)+amatr

######################################################################6
		amat36<- t(amatr)
		amat26 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat36[, kk]%%10 == 1)
			idy <- which(amat36[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]<100 & idx[ii]!=idy[jj]){
                  amat26[idy[jj], idx[ii]] <- 100
      				          }}
	            		}	
        			}
 		   }
		 amatr<-amat26+t(amat26)+amatr		

#################################################################7

   		 amat27 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in S) {
     			   idx <- which(amatr[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat27[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat27+t(amat27)+amatr

################################################################8

             amat28 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1)) {
            	for (ii in 1:(lenidx - 1)) {
                	for (jj in (ii + 1):lenidx) {
			if(amatr[idx[ii], idx[jj]]%%100<10){
                  amat28[idx[ii], idx[jj]] <- 10
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat28+t(amat28)+amatr
	
#################################################################9

			amat29 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%100<10 & idx[ii]!=idy[jj]){
                  amat29[idx[ii], idy[jj]] <- 10
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat29+t(amat29)+amatr

##################################################################10

   		 amat20 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amatr[, kk]%%100>9)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]%%100<10){
				   amat20[idx[ii], idx[jj]] <- 10
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat20+t(amat20)+amatr
	
	}
	Mn<-c()
	Cn<-c()
	for(i in M){
		Mn<-c(Mn,which(rownames(amat)==i))}
	for(i in C){
		Cn<-c(Cn,which(rownames(amat)==i))}
	if(length(Mn)>0&length(Cn)>0){
		fr<-amatr[-c(Mn,Cn),-c(Mn,Cn)]}
	if(length(Mn)>0&length(Cn)==0){
		fr<-amatr[-c(Mn),-c(Mn)]}
	if(length(Mn)==0&length(Cn)>0){
		fr<-amatr[-c(Cn),-c(Cn)]}
	if(length(Mn)==0&length(Cn)==0){
		fr<-amatr}	
	if(plot==TRUE){
		plotfun(fr,...)}
	if(showmat==FALSE){
		invisible(fr)}
	else{return(fr)}
}
##############################################################################
##############################################################################
`SG`<-function (amat,M=c(),C=c(),showmat=TRUE, plot=FALSE, plotfun = plotGraph, ...)
{	
	if(class(amat) == "igraph" || class(amat) == "graphNEL" || class(amat) == "character") {
		amat<-grMAT(amat)}
	if(class(amat) == "matrix"){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}	
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(class(amat) != "matrix") {	
	stop("'object' is not in a valid form")}

	S<-C
	St<-c()
	while(identical(S,St)==FALSE)
	{
		St<-S
		for(j in S){
			for(i in rownames(amat)){
				if(amat[i,j]%%10 == 1){
					S<-c(i,S[S!=i])}
				}
			}	
	}		

	amatr<-amat
	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr

############################################################1
			amat21 <- amatr
 	   for (kk in M) {
		  for (kkk in 1:ncol(amatr)) {
			amat31<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat31[kk,kkk]<-amatr[kkk,kk]
				amat31[kkk,kk]<-amatr[kk,kkk]}
        		
       		 idx <- which(amat31[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat31[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amat21[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}
    		
		amatr<-amat21	

################################################################2	

		amat22 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat22[idx[ii], idy[jj]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat22+amatr

#################################################################3

		amat33<- t(amatr)
		amat23 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat33[, kk]%%10 == 1)
			idy <- which(amat33[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]%%10==0 & idx[ii]!=idy[jj]){
                  amat23[idy[jj], idx[ii]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat23+amatr

##################################################################4
			amat24 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%100>9)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat24[idx[ii], idy[jj]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat24+amatr

####################################################################5		



		 amat35<- t(amatr)
   		 amat25 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amat35[, kk]%%10 == 1)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat25[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat25+t(amat25)+amatr

######################################################################6
		amat36<- t(amatr)
		amat26 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat36[, kk]%%10 == 1)
			idy <- which(amat36[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]<100 & idx[ii]!=idy[jj]){
                  amat26[idy[jj], idx[ii]] <- 100
      				          }}
	            		}	
        			}
 		   }
		 amatr<-amat26+t(amat26)+amatr		

#################################################################7

   		 amat27 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in S) {
     			   idx <- which(amatr[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat27[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat27+t(amat27)+amatr

################################################################8

             amat28 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1)) {
            	for (ii in 1:(lenidx - 1)) {
                	for (jj in (ii + 1):lenidx) {
			if(amatr[idx[ii], idx[jj]]%%100<10){
                  amat28[idx[ii], idx[jj]] <- 10
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat28+t(amat28)+amatr
	
#################################################################9

			amat29 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%100<10 & idx[ii]!=idy[jj]){
                  amat29[idx[ii], idy[jj]] <- 10
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat29+t(amat29)+amatr

##################################################################10

   		 amat20 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amatr[, kk]%%100>9)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]%%100<10){
				   amat20[idx[ii], idx[jj]] <- 10
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat20+t(amat20)+amatr
	
	}

		for(i in 1:ncol(amatr)) {
		for(j in 1:ncol(amatr)) {
			if(amatr[i,j]%%100>9){
				amatr[i,j]<-10
				for(k in 1:ncol(amatr)){
					if(amatr[k,j]==100){
						amatr[j,k]<-1
						amatr[k,j]<-0
						}}
				}			
			}}
	
	SS<-S
	SSt<-c()	
	while(identical(SS,SSt)==FALSE)
	{
		SSt<-SS
		for(j in SS){
			for(i in rownames(amat)) {
				if(amatr[i,j]%%10 == 1){
					SS<-c(i,SS[SS!=i])}
				}
			}	
	}	
		
		for(i in SS){
		for(j in SS) {	
			if(amatr[i,j]%%10==1){
				amatr[i,j]<-10
				amatr[j,i]<-10}}}



	Mn<-c()
	Cn<-c()
	for(i in M){
		Mn<-c(Mn,which(rownames(amat)==i))}
	for(i in C){
		Cn<-c(Cn,which(rownames(amat)==i))}
	if(length(Mn)>0&length(Cn)>0){
		fr<-amatr[-c(Mn,Cn),-c(Mn,Cn)]}
	if(length(Mn)>0&length(Cn)==0){
		fr<-amatr[-c(Mn),-c(Mn)]}
	if(length(Mn)==0&length(Cn)>0){
		fr<-amatr[-c(Cn),-c(Cn)]}
	if(length(Mn)==0&length(Cn)==0){
		fr<-amatr}	
	if(plot==TRUE){
		plotfun(fr,...)}
	if(showmat==FALSE){
		invisible(fr)}
	else{return(fr)}

}
##############################################################################
##############################################################################
`AG`<-function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{	
	if(class(amat) == "igraph" || class(amat) == "graphNEL" || class(amat) == "character") {
		amat<-grMAT(amat)}
	if(class(amat) == "matrix"){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}	
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(class(amat) != "matrix") {	
	stop("'object' is not in a valid form")}

	S<-C
	St<-c()
	while(identical(S,St)==FALSE)
	{
		St<-S
		for(j in S){
			for(i in rownames(amat)){
				if(amat[i,j]%%10 == 1){
					S<-c(i,S[S!=i])}
				}
			}	
	}		

	amatr<-amat
	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr


############################################################1
			amat21 <- amatr
 	   for (kk in M) {
		  for (kkk in 1:ncol(amatr)) {
			amat31<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat31[kk,kkk]<-amatr[kkk,kk]
				amat31[kkk,kk]<-amatr[kk,kkk]}
        		
       		 idx <- which(amat31[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat31[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amat21[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}
    		
		amatr<-amat21	

################################################################2	

		amat22 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat22[idx[ii], idy[jj]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat22+amatr

#################################################################3

		amat33<- t(amatr)
		amat23 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat33[, kk]%%10 == 1)
			idy <- which(amat33[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]%%10==0 & idx[ii]!=idy[jj]){
                  amat23[idy[jj], idx[ii]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat23+amatr

##################################################################4
			amat24 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%100>9)
			idy <- which(amatr[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]){
                  amat24[idx[ii], idy[jj]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat24+amatr

####################################################################5		



		 amat35<- t(amatr)
   		 amat25 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amat35[, kk]%%10 == 1)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat25[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat25+t(amat25)+amatr

######################################################################6
		amat36<- t(amatr)
		amat26 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in M) {
        		idx <- which(amat36[, kk]%%10 == 1)
			idy <- which(amat36[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idy[jj], idx[ii]]<100 & idx[ii]!=idy[jj]){
                  amat26[idy[jj], idx[ii]] <- 100
      				          }}
	            		}	
        			}
 		   }
		 amatr<-amat26+t(amat26)+amatr		

#################################################################7

   		 amat27 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in S) {
     			   idx <- which(amatr[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100){
				   amat27[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat27+t(amat27)+amatr

################################################################8

             amat28 <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in S) {
        		idx <- which(amatr[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1)) {
            	for (ii in 1:(lenidx - 1)) {
                	for (jj in (ii + 1):lenidx) {
			if(amatr[idx[ii], idx[jj]]%%100<10){
                  amat28[idx[ii], idx[jj]] <- 10
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat28+t(amat28)+amatr
	
#################################################################9

			amat29 <- matrix(rep(0,length(amat)),dim(amat))
    		 	for (kk in M) {
        		idx <- which(amatr[, kk]%%10 == 1)
			idy <- which(amatr[, kk]%%100> 9)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%100<10 & idx[ii]!=idy[jj]){
                  amat29[idx[ii], idy[jj]] <- 10
      				          }}
	            		}	
        			}
 		   }
		amatr<-amat29+t(amat29)+amatr

##################################################################10

   		 amat20 <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in M) {
     			   idx <- which(amatr[, kk]%%100>9)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]%%100<10){
				   amat20[idx[ii], idx[jj]] <- 10
               					 }}
           					 }
      			  }
   		 }
		 amatr<-amat20+t(amat20)+amatr
	
	}

		for(i in 1:ncol(amatr)) {
		for(j in 1:ncol(amatr)) {
			if(amatr[i,j]%%100>9){
				amatr[i,j]<-10
				for(k in 1:ncol(amatr)){
					if(amatr[k,j]==100){
						amatr[j,k]<-1
						amatr[k,j]<-0
						}}
				}			
			}}

	
	SS<-S
	SSt<-c()	
	while(identical(SS,SSt)==FALSE)
	{
		SSt<-SS
		for(j in SS){
			for(i in rownames(amat)) {
				if(amatr[i,j]%%10 == 1){
					SS<-c(i,SS[SS!=i])}
				}
			}	
	}	
		
		for(i in SS){
		for(j in SS) {	
			if(amatr[i,j]%%10==1){
				amatr[i,j]<-10
				amatr[j,i]<-10}}}



		amatn <- amatr
 	   for (kk in 1:ncol(amatr)) {
		  for (kkk in 1:ncol(amatr)) {
			amat3n<-amatr
		 	if (amatr[kkk,kk]%%10==1|amatr[kk,kkk]%%10==1) {
				amat3n[kk,kkk]<-amatr[kkk,kk]
				amat3n[kkk,kk]<-amatr[kk,kkk]}
        		
       		 idx <- which(amat3n[, kk]%%10 == 1)
        		lenidx <- length(idx)
        		if ((lenidx > 1&amat3n[kkk,kk]%%10==1)) {
           	 	for (ii in 1:(lenidx )) {
               		 #for (jj in (ii + 1):lenidx) {
                 	 	if(amatr[idx[ii], kkk]%%10==0&idx[ii]!=kkk){
			 	amatn[idx[ii], kkk] <- TRUE
                		}}
            		}
       		 }
			}
    		
	
	amatt<-2*amat
	while(identical(amatr,amatt)==FALSE)
	{
		amatt<-amatr

		
	 amat27n <- matrix(rep(0,length(amat)),dim(amat))
  		  for (kk in 1:ncol(amatr)) {
     			   idx <- which(amatn[, kk]>99)
     			   lenidx <- length(idx)
       			 if ((lenidx > 1)) {
         			   for (ii in 1:(lenidx - 1)) {
          			   for (jj in (ii + 1):lenidx) {
				   if(amatr[idx[ii], idx[jj]]<100 && (amatn[kk,idx[ii]]%%10==1 || amatn[kk,idx[jj]]%%10==1)){
				   amat27n[idx[ii], idx[jj]] <- 100
               					 }}
           					 }
      			  }
   		 }
		amatn<-amat27n+t(amat27n)+amatn
		amatr<-amat27n+t(amat27n)+amatr

		amat22n <- matrix(rep(0,length(amat)),dim(amat))
    		 for (kk in 1:ncol(amatr)) {
        		idx <- which(amatn[, kk]%%10 == 1)
			idy <- which(amatn[, kk]> 99)
        		lenidx <- length(idx)
			lenidy <- length(idy)
        		if ((lenidx > 0 & lenidy >0)) {
            	for (ii in 1:(lenidx)) {
                	for (jj in 1:lenidy) {
			if(amatr[idx[ii], idy[jj]]%%10==0 & idx[ii]!=idy[jj]& amatn[kk,idy[jj]]%%10==1){

                  amat22n[idx[ii], idy[jj]] <- 1
      				          }}
	            		}	
        			}
 		   }
		amatn<-amat22n+amatn
		amatr<-amat22n+amatr
	}
	
	for(i in 1:ncol(amatr)) {
		for(j in 1:ncol(amatr)) {
			if(amatr[i,j]==101){
				amatr[i,j]<-1
				amatr[j,i]<-0}}}

	Mn<-c()
	Cn<-c()
	for(i in M){
		Mn<-c(Mn,which(rownames(amat)==i))}
	for(i in C){
		Cn<-c(Cn,which(rownames(amat)==i))}
	if(length(Mn)>0&length(Cn)>0){
		fr<-amatr[-c(Mn,Cn),-c(Mn,Cn)]}
	if(length(Mn)>0&length(Cn)==0){
		fr<-amatr[-c(Mn),-c(Mn)]}
	if(length(Mn)==0&length(Cn)>0){
		fr<-amatr[-c(Cn),-c(Cn)]}
	if(length(Mn)==0&length(Cn)==0){
		fr<-amatr}	
	if(plot==TRUE){
		plotfun(fr, ...)}
	if(showmat==FALSE){
		invisible(fr)}
	else{return(fr)}
}

##############################################################################
##############################################################################
Max<-function(amat)
{
	if(class(amat) == "igraph" || class(amat) == "graphNEL" || class(amat) == "character") {
		amat<-grMAT(amat)}
	if(class(amat) == "matrix"){
		if(nrow(amat)==ncol(amat)){
			if(length(rownames(amat))!=ncol(amat)){
	 			rownames(amat)<-1:ncol(amat)
	 			colnames(amat)<-1:ncol(amat)}	
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(class(amat) != "matrix") {	
	stop("'object' is not in a valid form")}

	na<-ncol(amat)
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(dim(at)[1]>0){
	for(i in 1:dim(at)[1]){	
		S<-at[i,]
		St<-c()
		while(identical(S,St)==FALSE){
			St<-S
			for(j in S){
				for(k in 1:na){
					if(amat[k,j]%%10 == 1){
						S<-c(k,S[S!=k])}
					}
				}}
	one<-at[i,1]
	two<-at[i,2]
	onearrow<-c()
	onearc<-c()	
	twoarrow<-c()
	twoarc<-c()		
	Sr<-S
	Sr<-Sr[Sr!=at[i,1]]
	Sr<-Sr[Sr!=at[i,2]]
	if(length(Sr)>0){
		for(j in Sr){
			if(amat[one,j]%%10==1){
				onearrow<-c(onearrow,j)}
 			if(amat[one,j]>99){
				onearc<-c(onearc,j)}}
		for(j in Sr){
			if(amat[two,j]%%10==1){
				for(k in onearc){
					if(j==k || amat[j,k]>99){
						amat[at[i,2],at[i,1]]<-1}}
				twoarrow<-c(twoarrow,j)}
			if(amat[two,j]>99){
				for(k in onearrow){
					if(j==k || amat[j,k]>99){
						amat[at[i,1],at[i,2]]<-1}}
				for(k in onearc){
					if(j==k || amat[j,k]>99){
						amat[at[i,1],at[i,2]]<-100
						amat[at[i,2],at[i,1]]<-100}}
				twoarc<-c(twoarc,j)}}}				
		if(length(c(onearc,onearrow,twoarc,twoarrow))>0){
			
			for(j in c(onearc,onearrow,twoarc,twoarrow)){
				Sr<-Sr[Sr!=j]}
		onearcn<-c()
		twoarcn<-c()
		onearrown<-c()
		twoarrown<-c()
		while(length(Sr)>0){
			for(l in onearc){
				for(j in Sr){
 					if(amat[l,j]>99){
						onearcn<-c(onearcn,j)}}}
			for(l in onearrow){
				for(j in Sr){
 					if(amat[l,j]>99){
						onearrown<-c(onearrown,j)}}}

			for(l in twoarc){
				for(j in Sr){
					if(amat[l,j]>99){
						for(k in onearrow){
							if(j==k || amat[j,k]>99){
								amat[at[i,1],at[i,2]]<-1}}
						for(k in onearc){
							if(j==k || amat[j,k]>99){
								amat[at[i,1],at[i,2]]<-100
								amat[at[i,2],at[i,1]]<-100}}
						twoarcn<-c(twoarcn,j)}}}
			for(l in twoarrow){
				for(j in Sr){
					if(amat[l,j]>99){
						for(k in onearc){
							if(j==k || amat[j,k]>99){
								amat[at[i,1],at[i,2]]<-100
								amat[at[i,2],at[i,1]]<-100}}
						twoarrown<-c(twoarrown,j)}}}
			if(length(c(onearcn,onearrown,twoarcn,twoarrown))==0){
				break}			
			for(j in c(onearcn,onearrown,twoarcn,twoarrown)){
				Sr<-Sr[Sr!=j]}	
			onearc<-onearcn
			twoarc<-twoarcn
			onearrow<-onearrown
			twoarrow<-twoarrown					
			onearcn<-c()
			twoarcn<-c()
			onearrown<-c()
			twoarrown<-c()}}}}
	return(amat)
}
#####################################################################################################
######################################################################################################
msep<-function(a,alpha,beta,C=c()){
	if(class(a) == "igraph" || class(a) == "graphNEL" || class(a) == "character") {
		a<-grMAT(a)}
	if(class(a) == "matrix"){
		if(nrow(a)==ncol(a)){
			if(length(rownames(a))!=ncol(a)){
	 			rownames(a)<-1:ncol(a)
	 			colnames(a)<-1:ncol(a)}	
					}
		else {
      	  stop("'object' is not in a valid adjacency matrix form")}}
	if(class(a) != "matrix") {	
	stop("'object' is not in a valid form")}

		M<-rem(rownames(a),c(alpha,beta,C))
		ar<-Max(RG(a,M,C))	

	#aralpha<-as.matrix(ar[SPl(c(alpha,beta),alpha),SPl(c(alpha,beta),beta)])
	#arbeta<-as.matrix(ar[SPl(c(alpha,beta),beta),SPl(c(alpha,beta),alpha)])
	#for(i in 1:length(alpha)){
	#	for(j in 1:length(beta)){
	#		if(aralpha[j,i]!=0 || arbeta[j,i]!=0){
	#			return("NOT separated")
	#			break
	#			break}}}
	if(max(ar[as.character(beta),as.character(alpha)]+ar[as.character(alpha),as.character(beta)]!=0)){
		return(FALSE)}
			return(TRUE)
}
############################################################################
############################################################################
`MRG` <- function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	return(Max(RG(amat,M,C,showmat,plot, plotfun = plotGraph, ...)))
}
##########################################################################
##########################################################################
`MSG` <- function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	return(Max(SG(amat,M,C,showmat,plot, plotfun = plotGraph, ...)))
}
############################################################################
###########################################################################
`MAG`<-function (amat,M=c(),C=c(),showmat=TRUE,plot=FALSE, plotfun = plotGraph, ...)
{
	return(Max(AG(amat,M,C,showmat,plot, plotfun = plotGraph, ...)))
}
############################################################################
############################################################################
#Plot<-function(a)
#{
#	if(class(a) == "igraph" || class(a) == "graphNEL" || class(a) == "character") {
#		a<-grMAT(a)}
#	if(class(a) == "matrix"){
#		if(nrow(a)==ncol(a)){
#			if(length(rownames(a))!=ncol(a)){
#	 			rownames(a)<-1:ncol(a)
#	 			colnames(a)<-1:ncol(a)}	
#			l1<-c()
#			l2<-c()
#			for (i in 1:nrow(a)){
#				for (j in i:nrow(a)){
#					if (a[i,j]==1){
#						l1<-c(l1,i,j)
#						l2<-c(l2,2)}		
#					if (a[j,i]%%10==1){
#						l1<-c(l1,j,i)
#						l2<-c(l2,2)}
#					if (a[i,j]==10){
#						l1<-c(l1,i,j)
#						l2<-c(l2,0)}
#					if (a[i,j]==11){
#						l1<-c(l1,i,j,i,j)
#						l2<-c(l2,2,0)}
#					if (a[i,j]==100){
#						l1<-c(l1,i,j)
#						l2<-c(l2,3)}
#					if (a[i,j]==101){
#						l1<-c(l1,i,j,i,j)
#						l2<-c(l2,2,3)}
#					if (a[i,j]==110){
#						l1<-c(l1,i,j,i,j)
#						l2<-c(l2,0,3)}
#					if (a[i,j]==111){
#						l1<-c(l1,i,j,i,j,i,j)
#						l2<-c(l2,2,0,3)}
#					}
#				}
#			}
#		else {
#      	  stop("'object' is not in a valid adjacency matrix form")
#    	}
#    	if(length(l1)>0){
#		l1<-l1-1
#    		agr<-graph(l1,n=nrow(a),directed=TRUE)}
#    	if(length(l1)==0){
#    		agr<-graph.empty(n=nrow(a), directed=TRUE)
#		return(plot(agr,vertex.label=rownames(a)))}
#    	return( tkplot(agr, layout=layout.kamada.kawai, edge.curved=FALSE:TRUE, vertex.label=rownames(a),edge.arrow.mode=l2))	
#	}
#	else {
#        stop("'object' is not in a valid format")}
#}
############################################################################
############################################################################
MarkEqRcg<-function(amat,bmat)
	{
	if(class(amat) == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat) == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat) == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	if(class(bmat) == "igraph"){
		bmat<-grMAT(bmat)}
	if(class(bmat) == "graphNEL"){
		bmat<-grMAT(bmat)}
	if(class(bmat) == "character"){
		bmat<-grMAT(bmat)}
	if( length(rownames(bmat))!=ncol(bmat) | length(colnames(bmat))!=ncol(bmat)){
		rownames(bmat)<-1:ncol(bmat)
		colnames(bmat)<-1:ncol(bmat)}
	bmat<-bmat[rownames(amat),colnames(amat),drop=FALSE]
	na<-ncol(amat)
	nb<-ncol(bmat)
	if(na != nb){
		return(FALSE)}
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	bt<-which(bmat+t(bmat)+diag(na)==0,arr.ind=TRUE)
	if(identical(at,bt)==FALSE){
		return(FALSE)}
	ai<-c()
	bi<-c()
	if(dim(at)[1]!=0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					ai<-c(ai,j)}
			if((bmat[bt[i,1],j]%%10==1 || bmat[bt[i,1],j]>99) &&(bmat[bt[i,2],j]%%10==1 || bmat[bt[i,2],j]>99) ){
				bi<-c(bi,j)}}
			if(identical(ai,bi)==FALSE){
				return(FALSE)}
		ai<-c()
		bi<-c()}}
	return(TRUE)
}
########################################################################################
########################################################################################
MarkEqMag<-function(amat,bmat)
	{
	if(class(amat) == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat) == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat) == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	if(class(bmat) == "igraph"){
		bmat<-grMAT(bmat)}
	if(class(bmat) == "graphNEL"){
		bmat<-grMAT(bmat)}
	if(class(bmat) == "character"){
		bmat<-grMAT(bmat)}
	if( length(rownames(bmat))!=ncol(bmat) | length(colnames(bmat))!=ncol(bmat)){
		rownames(bmat)<-1:ncol(bmat)
		colnames(bmat)<-1:ncol(bmat)}
	bmat<-bmat[rownames(amat),colnames(amat),drop=FALSE]
	na<-ncol(amat)
	nb<-ncol(bmat)
	if(na != nb){
		return(FALSE)}
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	bt<-which(bmat+t(bmat)+diag(na)==0,arr.ind=TRUE)
	if(identical(at,bt)==FALSE){
		return(FALSE)}
	a<-c()
	b<-c()
	if(length(at)>0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					a<-rbind(a,c(at[i,1],j,at[i,2]))}
			if((bmat[bt[i,1],j]%%10==1 || bmat[bt[i,1],j]>99) &&(bmat[bt[i,2],j]%%10==1 || bmat[bt[i,2],j]>99) ){
					b<-rbind(b,c(bt[i,1],j,bt[i,2]))}}}
		if(identical(a,b)==FALSE){
			return(FALSE)}}
	ar<-which(amat%%10==1,arr.ind=TRUE)
	br<-which(bmat%%10==1,arr.ind=TRUE)
	ap<-c()
	bp<-c()		
	if(length(ar)>0){
		for(i in 1:dim(ar)[1]){
			for(j in 1:na){
				if((amat[ar[i,1],j]>99) &&(amat[ar[i,2],j]>99)){
					ap<-rbind(ap,c(ar[i,1],j,ar[i,2]))}}}}								
	if(length(br)>0){
		for(i in 1:dim(br)[1]){
			for(j in 1:nb){
				if((bmat[br[i,1],j]>99) &&(bmat[br[i,2],j]>99)){
					bp<-rbind(bp,c(br[i,1],j,br[i,2]))}}}}
	if(length(ap)>0){
		aptt<-ap
		apt<-c(ap,1)
		Qonen<-c()
		Qtwon<-c()
		while(length(apt)-length(aptt)>0){
			apt<-aptt
			for(i in (1:dim(ap)[1])){
				Qone<-ap[i,1]
				Qtwo<-ap[i,2]			
				while(length(Qone)>0){
					for(l in (1:length(Qone))){
						J<-which(((amat+t(amat)+diag(na))[ap[i,3],]==0) && ((amat[Qone[l],]>99) || (amat[,Qone[l]]%%100==1)))
						for(j in J){
							for(k in 1:dim(a)[1]){
								if(min(a[k,]==c(j,Qone[l],Qtwo[l]))==1){
									a<-rbind(a,ap[i,])
									aptt<-aptt[(1:dim(aptt)[1])[-i],]
									break
									break
									break
									break}}}
						Q<-which((amat[,ap[i,3]]%%10==1) && (amat[Qone[l],]>99))
						for(q in Q){
							for(k in 1:dim(a)[1]){
								if(min(a[k,]==c(q,Qone[l],Qtwo[l]))==1){
									Qtwon<-c(Qtwon,Qone[l])
									Qonen<-c(Qonen,q)}}}}
					Qtwo<-Qtwon
					Qone<-Qonen
					Qonen<-c()
					Qtwon<-c()}}}}
	if(length(bp)>0){
		bptt<-bp
		bpt<-c(bp,1)
		Qbonen<-c()
		Qbtwon<-c()
		while(length(bpt)-length(bptt)>0){
			bpt<-bptt
			for(i in (1:dim(bp)[1])){
				Qbone<-bp[i,1]
				Qbtwo<-bp[i,2]			
				while(length(Qbone)>0){
					for(l in (1:length(Qbone))){
						J<-which(((bmat+t(bmat)+diag(nb))[bp[i,3],]==0) && ((bmat[Qbone[l],]>99) || (bmat[,Qbone[l]]%%100==1)))
						for(j in J){
							for(k in 1:dim(b)[1]){
								if(min(b[k,]==c(j,Qbone[l],Qbtwo[l]))==1){
									b<-rbind(b,bp[i,])
									bptt<-bptt[(1:dim(bptt)[1])[-i],]
									break
									break
									break
									break}}}
						Qb<-which((bmat[,bp[i,3]]%%10==1) && (bmat[Qbone[l],]>99))
						for(q in Qb){
							for(k in 1:dim(b)[1]){
								if(min(b[k,]==c(q,Qbone[l],Qbtwo[l]))==1){
									Qbtwon<-c(Qbtwon,Qbone[l])
									Qbonen<-c(Qbonen,q)}}}}
					Qbtwo<-Qbtwon
					Qbone<-Qbonen
					Qbonen<-c()
					Qbtwon<-c()}}}}
	if(length(a)!=length(b)){
		return(FALSE)}
	f<-c()
	if((length(a)>0) && (length(b)>0)){
		for(i in 1:dim(a)[1]){
			for(j in 1:dim(b)[1]){
				f<-c(f,min(a[i,]==b[j,]))}
			if(max(f)==0){
				return(FALSE)}}}	
	return(TRUE)
}
##########################################################################################
###########################################################################################
RepMarUG<-function(amat)
{
	if(class(amat) == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat) == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat) == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	na<-ncol(amat)
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(dim(at)[1]!=0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					return(list(verify = FALSE, amat = NA))}}}}
	for(i in 1:na){
		for(j in 1:na){
			if(amat[i,j]==100){
				amat[i,j]<-10}
			if(amat[i,j]==1){
				amat[i,j]<-10
				amat[j,i]<-10}}}
	return(list(verify = TRUE, amat = amat))
}
########################################################################################
#########################################################################################
RepMarBG<-function(amat)
{
	if(class(amat) == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat) == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat) == "character"){
		amat<-grMAT(amat)}
	if( length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	na<-ncol(amat)
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(dim(at)[1]!=0){
		for(i in 1:dim(at)[1]){
			for(j in 1:na){
				if(amat[at[i,1],j]%%100>9 || amat[j,at[i,1]]%%10==1 || amat[at[i,2],j]%%100>9 || amat[j,at[i,2]]%%10==1){
					return(list(verify= FALSE, amat = NA))}}}}
	for(i in 1:na){
		for(j in 1:na){
			if(amat[i,j]==10){
				amat[i,j]<-100}
			if(amat[i,j]==1){
				amat[i,j]<-100
				amat[j,i]<-100}}}
	return(list(verify = TRUE, amat = amat))
}
########################################################################################
########################################################################################
RepMarDAG<-function(amat)
{
	if(class(amat) == "igraph"){
		amat<-grMAT(amat)}
	if(class(amat) == "graphNEL"){
		amat<-grMAT(amat)}
	if(class(amat) == "character"){
		amat<-grMAT(amat)}
	if(length(rownames(amat))!=ncol(amat) | length(colnames(amat))!=ncol(amat)){
		rownames(amat)<-1:ncol(amat)
		colnames(amat)<-1:ncol(amat)}
	na<-ncol(amat)
	full<-sort(unique(which(amat%%100>9,arr.ind=TRUE)[,1]))
	arc<-sort(unique(which(amat>99,arr.ind=TRUE)[,1]))
	arrow<-sort(unique(as.vector(which(amat%%10==1,arr.ind=TRUE))))
	S<-full[full!=full[1]]
	Ma<-full[1]
	while(length(S)>0){
		dim<-c()
		for(i in S){
			dim<-c(dim,length(which(amat[i,Ma]%%100>9,arr.ind=TRUE)))}
		s<-S[which(dim==max(dim))[1]]
		ns<-which(amat[s,Ma]%%100>9,arr.ind=TRUE)
		if(min(amat[ns,ns]+diag(length(ns)))==0){
			return(FALSE)}
		Ma<-c(Ma,s)
		S<-S[S!=s]}
	at<-which(amat+t(amat)==0,arr.ind=TRUE)
	ai<-c()
	if(length(at[,1])!=0){
		for(i in (1:length(at[,1]))){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10!=1) && (amat[at[i,2],j]<100) ){
					ai<-c(ai,j)}}
			if(length(ai)==0){
				break}
			for(j in ((1:na)[-ai])){
				if((max(amat[ai,j]>99)==1) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) && (amat[at[i,1],j]%%10!=1) && (amat[at[i,1],j]<100)){
					return(list(verify = FALSE, amat = NA))}}
		ai<-c()}}	
	for(i in Ma){
		v<-which(amat[i,]%%100>9)
		for(j in v){
			amat[i,j]<-1
			amat[j,i]<-0}}
	at<-which(amat+t(amat)+diag(na)==0,arr.ind=TRUE)
	if(length(at[,1])!=0){
		for(i in (1:length(at[,1]))){
			for(j in 1:na){
				if((amat[at[i,1],j]%%10==1 || amat[at[i,1],j]>99) &&(amat[at[i,2],j]%%10==1 || amat[at[i,2],j]>99) ){
					amat[at[i,1],j]<-1
					amat[at[i,2],j]<-1}}}}
	O<-c()
	Oc<-arrow
	while(length(Oc)>0){
		for(i in Oc){		
	 		if(max(amat[i,Oc]%%10)==0){
		O<-c(O,i)
		break}}
		Oc<-Oc[Oc!=i]}			
	for(i in arc){
		if(length(which(arrow==i))==0){
			O<-c(O,i)}}
	aarc<-which(amat>99,arr.ind=TRUE)
	if(length(aarc)[1]>0){
		for(i in 1:(length(aarc[,1]))){
		#	if(length(O)==0){
		#		amat[aarc[i,1],aarc[i,2]]<-0
			if(which(O==aarc[i,1])>which(O==aarc[i,2])){
				amat[aarc[i,1],aarc[i,2]]<-1}
			else{
				amat[aarc[i,1],aarc[i,2]]<-0}}}
	return(list(verify = TRUE, amat = amat))
}	
##################################################################################
