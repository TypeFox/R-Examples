# Copyright (c) 2015 by Sean Downey
# Authors: Sean Downey (sean@codexdata.com) and Guowei Sun (gwsun@umd.edu)
# This software is distributed under GPL-3.

map<-function(m1,m2){     
  
  if(length(m1)!=length(m2)) stop("vector has different length")
  #ALINE.map<-NULL
  #data(ALINE.map,envir=environment())
  map<-ALINE.map
  #map<-show.map()
  Aline.v<-c('a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z')
  Aline.f<-c("D","V","X","P","S","N","A","H","F","C","Z")       
  
  f<-function(x){
    y<-unlist(strsplit(x,""))
    paste(as.character(sapply(y,utf8ToInt)),sep='', collapse=' ')
  }
  
  if(length(m2)!=0){
  for(i in 1:length(m2)){
    val<-TRUE
    uc<-strsplit(m2[i],split="")[[1]]
    if(sum(uc[1]==Aline.v)==0) {
      val<-FALSE
    }
    if(length(uc)>1){
    for(j in 2:length(uc)){
      if((sum(uc[j]==Aline.f))==0){
        val<-FALSE
      }
    }
    }
    
    if(val){
    m<-data.frame(
      #IPA=m1[i],
      Aline=m2[i],
      U.Val=utf8ToInt(m1[i]),
      A.Val=f(m2[i]),
      stringsAsFactors=FALSE, row.names=NULL)
      map<-rbind(m,map, deparse.level = 0)   # bug found by Vincent
    }
    else{
      message(paste("Invalide mapping:",m1[i],"to",m2[i]))
    } 
  }
  }
  else{
    map=ALINE.map
  }
  return(map)
}