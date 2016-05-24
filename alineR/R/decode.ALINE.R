# Copyright (c) 2015 by Sean Downey
# Authors: Sean Downey (sean@codexdata.com) and Guowei Sun (gwsun@umd.edu)
# This software is distributed under GPL-3.

decode.ALINE <-
function(x,y,m1=NULL,m2=NULL) {
  map<-map(m1,m2)
  splitx<-strsplit(x,"")
  splity<-strsplit(y,"")    
  t<-splitx[[1]]   
  for(j in 1:length(splitx[[1]])){
    catch<-FALSE
    for(i in 1:length(map$U.Val)){
      if(splitx[[1]][j]==intToUtf8(map$U.Val[i])){
        catch=TRUE
        break
      }
    }
    if(!catch){
      message(paste("Invalid character:",splitx[[1]][j],"dropped in alignment"))
      t<-splitx[[1]][-j]
    }
  }
  splitx[[1]]<-t
  num<-length(splity[[1]])
  J<-rep(NA,length(splitx[[1]]))  
  j=1    
  
  for( i in 1:num-1)
  {
    if (!splity[[1]][i]==" " && splity[[1]][i+1]==" ")
    {
      J[j]=i
      j=j+1    
    }
  }

 aligned_y=rep(" ",num)
  
  p=0
  
 for (i in 1:(length(J)))
 { 
   if(splity[[1]][J[i]]=="|"||splity[[1]][J[i]]=="<"||splity[[1]][J[i]]=="-")   
   {
     aligned_y[J[i]]<-splity[[1]][J[i]]
     p=p+1
   }
   else{
   aligned_y[J[i]]=splitx[[1]][(i-p)]   
   }
 }
  z<-paste(aligned_y,sep=' ',collapse='')
  return(z)
}
