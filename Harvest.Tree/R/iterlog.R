iterlog <- function(modif.b, simutry, newsim, varn){ 
  flag=0
  aa <- rep(0, length(modif.b))
  for (i in 1:length(modif.b)){     
    aa[i] <- sum(modif.b[[i]]) 
  }
  if (max(aa) != -Inf){    
    delta <- max(aa)
    nodemax <- order(aa, decreasing=T)[1]    
    noden <- simutry[[nodemax]]
    label <- noden$label  
    d.name <- names(modif.b[[nodemax]])    
    r.name <- varn[!(varn%in%d.name)]    
    har.rule <- data.frame(rep(0,max(length(r.name),1)),rep(0,max(length(r.name),1)),rep(0,max(length(r.name),1)))
    
   #debug #2
   if(length(r.name)>0){
     rownames(har.rule) <- r.name
     colnames(har.rule) <- c("lower bound (numeric)","upper bound (numeric)","categrical")
     #debug  end
    for (t in 1:length(r.name)){            
      for (j in 1:length(varn)){         
        if (varn[j] == r.name[t]) har.rule[t,] <- noden$bounds[j,]  	
      }
    }
   }
   else{
     rownames(har.rule)= 'no constraint'
     colnames(har.rule) <- c("lower bound (numeric)","upper bound (numeric)","categrical")
     har.rule[1:length(varn)]=rep(NA,length(varn))
     flag=1
   }
  }
  else{     
    har.rule <- 0     
    delta <- -Inf
    label <- 0
  }             
  return(list(har.rule=har.rule, delta=delta,label=label,flag=flag))
}	