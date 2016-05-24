# Copyright (c) 2015 by Sean Downey
# Authors: Sean Downey (sean@codexdata.com) and Guowei Sun (gwsun@umd.edu)
# This software is distributed under GPL-3.
# Includes modified C++ code originally developed by Greg Kondrak (c) 2000
aline <-
function(w1,w2,sim=FALSE,m1=NULL,m2=NULL,mark=FALSE,alignment=FALSE,...){
  n1<-which(is.na(w1)==TRUE)
  n2<-which(is.na(w2)==TRUE)
  
  if(!length(n1)==0 || !length(n2)==0){
    x1<-as.character(w1[-c(n1,n2)])
    y1<-as.character(w2[-c(n1,n2)])
  }else{
    x1<-as.character(w1)
    y1<-as.character(w2)
  }
  ind<-as.integer(levels(as.factor(c(n1,n2))))
  
  
	x<-encode.ALINE(x1, mark=mark,m1,m2)
	y<-encode.ALINE(y1, mark=mark,m1,m2)
	
	names(x)<-NULL
	names(y)<-NULL
	`%nin%` <- function (x, table) match(x, table, nomatch = 0L) == 0L   
	
	if (length(x) == 0 | length (y) == 0 ) stop("Vectors of length zero not allowed.") 
	if (length(x) != length (y) ) stop("Vectors not the same length.")   
	if (NA %in% x | NA %in% y ) stop("NA not permitted.")   
	
	a<-apply(cbind(x,y),1,FUN=raw.alignment,...)    
	names(a)<-NULL
	num=length(a)

	if (sim==F) {
		b<-apply(cbind(x,x),1,FUN=raw.alignment,...)         
		c<-apply(cbind(y,y),1,FUN=raw.alignment,...) 
    scores<-1-((2*a[[1]][[2]])/(b[[1]][[2]]+c[[1]][[2]]))
    
    if(!num==1){
    for (p in 2:num){
		scores<-c(scores,1-((2*a[[p]][[2]])/(b[[p]][[2]]+c[[p]][[2]]))) # ALINE Distance as per Downey et al. 2008
    }}
    
	} else 
    {
      scores<-a[[1]][[2]]
      if(!num==1){
      for (p in 2:num){
      scores<-c(scores,a[[p]][[2]]) 
      }}
    }
  
  if(length(ind)>0){
    for(i in 1:length(ind)){
      scores<-append(scores,NA,ind[i]-1)
    }
  }

	if(alignment==T){
		
    for(i in 1:length(x)){
      a[[i]]$alignment1<-decode.ALINE(x1[i],a[[i]]$alignment1,m1,m2)
      a[[i]]$alignment2<-decode.ALINE(y1[i],a[[i]]$alignment2,m1,m2)
    }
    
    p=1
    a1<-vector()
    a2<-vector()
    for(i in 1:length(scores)){
      if(is.na(scores[i])){
        a1[i]<-NA
        a2[i]<-NA
      }
      else{
        a1[i]<-a[[p]]$alignment1
        a2[i]<-a[[p]]$alignment2
        p=p+1
      }
    }
    names<-vector()
    for(i in 1:length(scores)){
      names[i]<-paste("pair",i,sep="")
    }
    M=rbind(w1,w2,scores,a1,a2)
    ###
    M<-as.data.frame(M)
    colnames(M)<-names
    return(M)
	} else return(scores)
}
