"build.plot" <-
function(parent){
	    nmut<- length(parent$child)
	    level<- rep(0,nmut)
	    for(i in 1:nmut) {
	    if (is.na(parent$parent[i]))
	        parent$parent[i]<-""
	  #  if (is.na(parent[3,i]))
	  #      parent[3,i]<-""    
	    if (level[i]==0 & (parent$parent[i] == "Root")) 
	        level[i]<-2             
	        else if (parent$child[i] == "Root")
	                level[i]<-1
	        else if (parent$parent[i]=="")
	                level[i]<- 999
	    }
	    
	    while (min(level) ==0) {   
	        for(i in 1:nmut){
	        if (level[i]==0) {
	        x<- match(parent$parent[i], parent$child)
	        if (level[x] > 0) level[i] <- level[x] +1
	        }
	    }    
	 }
	    
	    leaves<- rep(0,nmut)
	    for (i in 1:nmut) {
	        x<-match(parent$child[i], parent$parent)
	        leaves[i] <- is.na(x) 
	    }
	
	
	levelnodes<-rep(0,max(level))
	for (i in 1:max(level)){
	    c<- 0
	    if (i==1) levelnodes[i]<-1
	    else {
	        for (j in 1:nmut) {
	            if (level[j]== i) {
	            c<- c+1
	            levelnodes[i] <-c }
	            }
	        }
	    }
	
	levelgrp<-array(0, dim=c(max(level),sum(leaves)))
	l<-1
	for (j in 1:nmut) {
	    if (level[j]== 2) {
	        levelgrp[2,l]<-parent$child[j]
	        l<-l+1
	        } 
	    }
	for (lev in 2:max(level)){
	    levelgrp[1,1]<-"Root"
	    nodeadd<-1
	    for (l in 1:levelnodes[lev]) {
	        for (j in 1:nmut) {
	            if (parent$parent[j]==levelgrp[lev,l]) {
	                levelgrp[(lev+1),nodeadd]<-parent$child[j]
	                nodeadd<- nodeadd+1
	                }
	            }
	        }
	    }
  plotinfo<-list(levelgrp=levelgrp, nmut=nmut, level=level)
  return(plotinfo)
}

