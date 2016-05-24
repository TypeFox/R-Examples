`posterior.comparison.freqMAP` <-
function(group1,group2)
{
  if(!inherits(group1,"freqMAP")){
    stop("group1 must be of class freqMAP")
  }
  if(!inherits(group2,"freqMAP")){
    stop("group2 must be of class freqMAP")
  }
  if(group1$x.label!=group2$x.label){
    stop("x.label must be identical in the two groups")
  }
  if(sum(group1$cat.ma[,group1$x.label]!=group2$cat.ma[,group2$x.label])>0){
    stop("X columns must be identical in the two groups.")
  }

  if(any(dim(group1$post.samples) != dim(group2$post.samples))){
    stop("Posterior samples arrays in group1 and group2 are not of equal dimension")
  }

  if(any(group1$cat.names != group2$cat.names)){
    stop("cat.names must be the same in both groups.")
  }
  cat.names <- group1$cat.names
  
  post.comparison <- data.frame(x=group1$cat.ma[,group1$x.label],
                                n1=group1$cat.ma$n,
                                n2=group2$cat.ma$n)
  names(post.comparison)[1] <- group1$x.label
  
  #Calculate posterior probability that frequency in group 1 > freq in group 2,
  #for each category.
  for(a in cat.names){
    post.comparison[,paste(a,".gr1.gt.gr2",sep="")] <- NA
  }

  #The number of posterior samples
  G <- dim(group1$post.samples)[1]
  
  for(i in 1:nrow(post.comparison)){
    for(j in 1:length(cat.names)){
      post.comparison[i,paste(cat.names[j],".gr1.gt.gr2",sep="")] <-
        sum(group1$post.samples[,j,i]>
            group2$post.samples[,j,i])/G
    }
  }

  #Calculate posterior intervals and probabilities on log odds ratios
  for(a in 1:nrow(post.comparison)){
    
    for(i in 1:(length(cat.names)-1)){
      for(j in (i+1):length(cat.names)){
        
        #column names for stats on lor: mean, lower and upper
        #post. bounds, and Prob(lor>0)
        cn <- paste(cat.names[j],".",cat.names[i],".lor",
                    c(".mean",".lpi",".upi",".p.gt.0"),sep="")
        
        for(cni in cn) post.comparison[a,cni] <- NA
        
        if(post.comparison$n1[a] > 1 && post.comparison$n2[a] > 1){
          
          lor <- log(group1$post.samples[,j,a]*group2$post.samples[,i,a]/
                     group2$post.samples[,j,a]/group1$post.samples[,i,a])
          post.comparison[a,cn] <-
            c(
              mean(lor),
              quantile(lor,probs=c(.025,.975),na.rm=TRUE),
              sum(lor>0)/G
              )
        }      
        
      } #for j
    } #for i
    
  } #for a
  
  return(post.comparison)
}

